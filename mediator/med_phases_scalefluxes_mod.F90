module med_phases_scalefluxes_mod

  use med_kind_mod          , only : CX=>SHR_KIND_CX, CS=>SHR_KIND_CS, CL=>SHR_KIND_CL, R8=>SHR_KIND_R8
  use med_constants_mod     , only : dbug_flag     => med_constants_dbug_flag
  use med_utils_mod         , only : chkerr        => med_utils_ChkErr
  use med_methods_mod       , only : fldchk     => med_methods_FB_FldChk
  use med_methods_mod       , only : FB_GetFldPtr  => med_methods_FB_GetFldPtr
  use med_methods_mod       , only : fldbun_getdata1d => med_methods_FB_getdata1d
  use med_internalstate_mod , only : InternalState, maintask, logunit
  use med_internalstate_mod , only : compocn, compice, compatm, comprof
  use perf_mod              , only : t_startf, t_stopf
  use shr_log_mod           , only : shr_log_error
  use ESMF                  , only : ESMF_SUCCESS

  implicit none
  private

  public :: med_phases_scalefreshwater_run  ! called from run sequence

  private :: scalefreshwater_get_precip
  private :: scalefreshwater_scale_precip

  character(*), parameter :: u_FILE_u  = &
       __FILE__

!-----------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------

  subroutine med_phases_scalefreshwater_run(gcomp, rc)

  !---------------------------------------
  ! balance freshwater fluxes between atmosphere & ocean (+ sea ice) to zero, 
  ! by scaling precip, such that the sum of precip, runoff and evap is zero
  ! this adjusts FBImp fields - so needs to be run :
  ! - before components merge freshwater fluxes in med_phases_prep_ocn_accum and med_phases_prep_ice
  ! - after evap is calculated in aoflux_run
  ! the sea ice fraction and evaporation haven't been calculated yet for the current
  ! timestep, so the sea ice freshwater balance is not an exact match for that added 
  ! during the ICE run stage
  !---------------------------------------

  use ESMF , only : ESMF_GridComp, ESMF_VMGet, ESMF_VMAllreduce, ESMF_REDUCE_SUM
  use ESMF , only : ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_SUCCESS
  use ESMF , only : ESMF_GridCompGet, ESMF_VM
  use med_constants_mod     , only : shr_const_pi
  use shr_reprosum_mod      , only : shr_reprosum_calc

  ! input/output variables
  type(ESMF_GridComp)    :: gcomp
  integer, intent(out) :: rc

  ! local variables
  type(InternalState) :: is_local
  type(ESMF_VM) :: vm
  integer             :: i, comm
  real(r8), pointer   :: evap(:), evap_si(:), rofl(:), rofi(:)
  real(r8), allocatable :: ocn_precip_sum(:), ocn_sum_weighted(:,:) ! local ocean sums
  real(r8), allocatable :: ice_precip_sum(:), ice_sum_weighted(:,:) ! local ice sums
  real(r8)            :: ocn_global_sum(2), ice_global_sum(2)       ! global ocean, ice sums
  real(r8)            :: local_sum(1), global_sum(2)            ! global ocean+ice sums
  real(r8)            :: precip_fact
  real(r8), pointer   :: ocn_areas(:), ice_areas(:)
  real(r8), pointer   :: ifrac(:)  ! ice fraction in ocean grid cell
  real(r8), pointer   :: ofrac(:)  ! non-ice fraction in ocean grid cell
  logical             :: first_call = .true. , sum_precip = .true.
  integer, parameter  :: ip=1, ifw=2 ! index for precip, freshwater
  integer, parameter  :: dbug_threshold = 20 ! threshold for writing debug information in this subroutine
  real(r8), parameter :: eps = 10.0_r8 * tiny(0.0_r8) ! threshold for zero
  character(len=*), parameter    :: subname='(med_phases_scalefreshwater_run)'
  !---------------------------------------

  call t_startf('MED:'//subname)
  if (dbug_flag > dbug_threshold) then
      call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)
  end if
  rc = ESMF_SUCCESS

  call ESMF_GridCompGet(gcomp, vm=vm, rc=rc)
  if (ChkErr(rc,__LINE__,u_FILE_u)) return

  ! Get the MPI communicator from the VM
  call ESMF_VMGet(vm, mpiCommunicator=comm, rc=rc)
  if (ChkErr(rc,__LINE__,u_FILE_u)) return

  ! Get the internal state
  nullify(is_local%wrap)
  call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
  if (ChkErr(rc,__LINE__,u_FILE_u)) return

  ocn_areas => is_local%wrap%mesh_info(compocn)%areas

  if ( first_call ) then
    ! test required fields exist
    if (fldchk(is_local%wrap%FBImp(compatm,compocn), 'Faxa_rainl', rc=rc) .and. &
        fldchk(is_local%wrap%FBImp(compatm,compocn), 'Faxa_rainc', rc=rc) .and. &
        fldchk(is_local%wrap%FBImp(compatm,compocn), 'Faxa_snowl', rc=rc) .and. &
        fldchk(is_local%wrap%FBImp(compatm,compocn), 'Faxa_snowc', rc=rc)) then
      sum_precip = .true.
    else if (fldchk(is_local%wrap%FBImp(compatm,compocn), 'Faxa_rain', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compatm,compocn), 'Faxa_snow', rc=rc)) then
      sum_precip = .false.
    else
      call shr_log_error(trim(subname)//": ERROR imported rain ocean fields for med_phases_scalefreshwater_run are missing ", &
        line=__LINE__, file=u_FILE_u, rc=rc)
      return
    end if

    if ( .not. (fldchk(is_local%wrap%FBMed_aoflux_o, 'Faox_evap', rc=rc) .and. &
        fldchk(is_local%wrap%FBImp(comprof,compocn), 'Forr_rofl', rc=rc) .and. &
        fldchk(is_local%wrap%FBImp(comprof,compocn), 'Forr_rofi', rc=rc))) then
      call shr_log_error(trim(subname)//": ERROR some ocean fields for med_phases_scalefreshwater_run are missing ", &
        line=__LINE__, file=u_FILE_u, rc=rc)
      return
    endif

    if (is_local%wrap%comp_present(compice) .and. &
        .not. ( fldchk(is_local%wrap%FBImp(compice,compice), 'Faii_evap', rc=rc))) then
      call shr_log_error(trim(subname)//": ERROR some ice fields for med_phases_scalefreshwater_run are missing ", &
        line=__LINE__, file=u_FILE_u, rc=rc)
      return
    endif
    first_call = .false.
  endif

  call fldbun_getdata1d(is_local%wrap%FBfrac(compocn), 'ofrac', ofrac, rc=rc)
  if (ChkErr(rc,__LINE__,u_FILE_u)) return

  allocate(ocn_precip_sum(size(ofrac)))
  allocate(ocn_sum_weighted(size(ofrac),2))

  ! First, get the precip fields
  call scalefreshwater_get_precip(sum_precip, is_local%wrap%FBImp(compatm,compocn), ocn_precip_sum, rc=rc)
  if (ChkErr(rc,__LINE__,u_FILE_u)) return

  ocn_sum_weighted(:,ip) = ocn_areas*ofrac*ocn_precip_sum   ! convert from a flux (km/m2/s) to a mass rate (kg/s)

  ! If cice IS PRESENT
  if (is_local%wrap%comp_present(compice)) then

    ice_areas => is_local%wrap%mesh_info(compice)%areas

    call fldbun_getdata1d(is_local%wrap%FBfrac(compice), 'ifrac', ifrac, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    allocate(ice_precip_sum(size(ifrac)))
    allocate(ice_sum_weighted(size(ifrac),2))

    call scalefreshwater_get_precip(sum_precip, is_local%wrap%FBImp(compatm,compice), ice_precip_sum, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ice_sum_weighted(:,ip) = ice_areas*ifrac*ice_precip_sum

  endif

  ! Second, add the runoff and evaporation to get a total freshwater flux

  call FB_GetFldPtr(is_local%wrap%FBMed_aoflux_o, 'Faox_evap' , evap, rc=rc)
  if (ChkErr(rc,__LINE__,u_FILE_u)) return

  call FB_GetFldPtr(is_local%wrap%FBImp(comprof,compocn), 'Forr_rofl' , rofl, rc=rc)
  if (ChkErr(rc,__LINE__,u_FILE_u)) return

  call FB_GetFldPtr(is_local%wrap%FBImp(comprof,compocn), 'Forr_rofi' , rofi, rc=rc)
  if (ChkErr(rc,__LINE__,u_FILE_u)) return

  ocn_sum_weighted(:,ifw) = ocn_sum_weighted(:,ip)+ocn_areas*(ofrac*evap + rofl + rofi)

  ! Sum runoff and total freshwater flux globally
  call shr_reprosum_calc(ocn_sum_weighted, ocn_global_sum, size(ofrac), size(ofrac), 2, &
                               commid=comm)

  global_sum = ocn_global_sum

  if (is_local%wrap%comp_present(compice)) then
    call FB_GetFldPtr(is_local%wrap%FBImp(compice,compice), 'Faii_evap' , evap_si, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ice_sum_weighted(:,ifw) = ice_sum_weighted(:,ip) + ice_areas*ifrac*evap_si

    call shr_reprosum_calc(ice_sum_weighted, ice_global_sum, size(ifrac), size(ifrac), 2, &
                              commid=comm)

    global_sum = global_sum + ice_global_sum
  endif

  if (maintask .and. (dbug_flag > dbug_threshold)) then
    write(logunit,'(a,ES26.18)') &
      trim(subname)//': global_precip_sum ', global_sum(ip)/(4.0_r8*shr_const_pi)
    write(logunit,'(a,ES26.18)') &
      trim(subname)//': global_fw_sum ', global_sum(ifw)/(4.0_r8*shr_const_pi)
  endif

  if (abs(global_sum(ip)) > eps) then
    ! Scale total freshwater to zero
    precip_fact = 1.0_r8 - (global_sum(ifw)/global_sum(ip))
  else
    if (maintask) write(logunit,'(a)') trim(subname)//': WARNING: global precip is zero, skipping scaling'
    precip_fact = 1.0_r8
  end if

  if (precip_fact < eps) then
    call shr_log_error(trim(subname)//": ERROR global freshwater flux is greater than global precip flux", &
      line=__LINE__, file=u_FILE_u, rc=rc)
    return
  end if

  if (maintask .and. (dbug_flag > dbug_threshold)) then
    write(logunit,'(a,ES26.18)') &
      trim(subname)//': Scaling rain & snow by non-unity precip_fact ', precip_fact
  endif

  call scalefreshwater_scale_precip(sum_precip, is_local%wrap%FBImp(compatm,compocn), precip_fact, rc=rc)
  if (ChkErr(rc,__LINE__,u_FILE_u)) return
  if (is_local%wrap%comp_present(compice)) then
    call scalefreshwater_scale_precip(sum_precip, is_local%wrap%FBImp(compatm,compice), precip_fact, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
  endif

  ! Check total freshwater flux globally is zero
  if (dbug_flag > dbug_threshold) then
    !check new global_fw_sum
    local_sum = 0
    call scalefreshwater_get_precip(sum_precip, is_local%wrap%FBImp(compatm,compocn), ocn_precip_sum, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    do i = 1, size(ofrac)
      local_sum(1) = local_sum(1) + ocn_areas(i)*(ofrac(i)*(ocn_precip_sum(i) + evap(i)) + rofl(i) + rofi(i))
    end do

    if (is_local%wrap%comp_present(compice)) then
      call scalefreshwater_get_precip(sum_precip, is_local%wrap%FBImp(compatm,compice), ice_precip_sum, rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return
      do i = 1, size(ifrac)
        local_sum(1) = local_sum(1) + ice_areas(i)*ifrac(i)*(ice_precip_sum(i) + evap_si(i))
      end do
    endif

    call ESMF_VMAllreduce(is_local%wrap%vm, senddata=local_sum, recvdata=global_sum, count=1, &
      reduceflag=ESMF_REDUCE_SUM, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (maintask) then
      write(logunit,'(a,ES26.18)') &
        trim(subname)//': global_fw_sum ', global_sum(1)/(4.0_r8*shr_const_pi)
    endif
  endif

  if (dbug_flag > dbug_threshold) then
    call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)
  end if
  call t_stopf('MED:'//subname)

  end subroutine med_phases_scalefreshwater_run

  subroutine scalefreshwater_get_precip(sum_precip, FB_comp, precip_sum, rc)

    ! For a given field bundle, return a combined precipitation array
    use ESMF                  , only : ESMF_FieldBundle

    logical, intent(in)                :: sum_precip
    type(ESMF_FieldBundle), intent(in) :: FB_comp
    real(r8), intent(out)              :: precip_sum(:)
    integer, intent(out)               :: rc

    ! local
    real(r8), pointer   :: rain(:),snow(:),rainl(:),rainc(:),snowl(:),snowc(:)

    rc = ESMF_SUCCESS

    if (sum_precip) then
      ! Sum rainc and rainl, following 
      ! https://github.com/access-nri/cmeps/blob/f84dd460bc0a14c55436e47bf6919fdba6633637/mediator/med_merge_mod.F90#L414-L419

      call FB_GetFldPtr(FB_comp, 'Faxa_rainl' , rainl, rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return

      call FB_GetFldPtr(FB_comp, 'Faxa_rainc' , rainc, rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return

      call FB_GetFldPtr(FB_comp, 'Faxa_snowl' , snowl, rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return

      call FB_GetFldPtr(FB_comp, 'Faxa_snowc' , snowc, rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return

      precip_sum = rainl + rainc + snowl + snowc

    else

      call FB_GetFldPtr(FB_comp, 'Faxa_rain' , rain, rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return

      call FB_GetFldPtr(FB_comp, 'Faxa_snow' , snow, rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return

      precip_sum = rain + snow

    endif

  end subroutine

  subroutine scalefreshwater_scale_precip(sum_precip, FB_comp, precip_fact,rc)

    ! For a given field bundle, scale precipitation by precip_fact
    use ESMF                  , only : ESMF_FieldBundle

    logical, intent(in)                :: sum_precip
    type(ESMF_FieldBundle), intent(in) :: FB_comp
    real(r8), intent(in)               :: precip_fact
    integer, intent(out)               :: rc

    ! local
    real(r8), pointer   :: rain(:),snow(:),rainl(:),rainc(:),snowl(:),snowc(:)

    rc = ESMF_SUCCESS

    if (sum_precip) then
      ! Sum rainc and rainl, following 
      ! https://github.com/access-nri/cmeps/blob/f84dd460bc0a14c55436e47bf6919fdba6633637/mediator/med_merge_mod.F90#L414-L419

      call FB_GetFldPtr(FB_comp, 'Faxa_rainl' , rainl, rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return

      call FB_GetFldPtr(FB_comp, 'Faxa_rainc' , rainc, rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return

      call FB_GetFldPtr(FB_comp, 'Faxa_snowl' , snowl, rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return

      call FB_GetFldPtr(FB_comp, 'Faxa_snowc' , snowc, rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return

      rainl = precip_fact * rainl
      rainc = precip_fact * rainc
      snowl = precip_fact * snowl
      snowc = precip_fact * snowc

    else

      call FB_GetFldPtr(FB_comp, 'Faxa_rain' , rain, rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return

      call FB_GetFldPtr(FB_comp, 'Faxa_snow' , snow, rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return

      rain = precip_fact * rain
      snow = precip_fact * snow

    endif

  end subroutine

end module med_phases_scalefluxes_mod
