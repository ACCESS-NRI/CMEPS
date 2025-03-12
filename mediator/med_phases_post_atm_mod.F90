module med_phases_post_atm_mod

  !-----------------------------------------------------------------------------
  ! Mediator phase for post atm calculations, maps atm->ice, atm->lnd, atm->ocn
  ! and atm->wav
  !-----------------------------------------------------------------------------
 
  implicit none
  private

  public :: med_phases_post_atm

  character(*), parameter :: u_FILE_u  = &
       __FILE__

!-----------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------

  subroutine med_phases_post_atm(gcomp, rc)

    !---------------------------------------
    ! map atm to ocn and atm to ice and atm to land
    !---------------------------------------

    use NUOPC_Mediator        , only : NUOPC_MediatorGet
    use ESMF                  , only : ESMF_Clock, ESMF_ClockIsCreated
    use ESMF                  , only : ESMF_GridComp, ESMF_GridCompGet, ESMF_FieldBundleGet
    use ESMF                  , only : ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_SUCCESS
    use med_kind_mod          , only : CX=>SHR_KIND_CX, CS=>SHR_KIND_CS, CL=>SHR_KIND_CL, R8=>SHR_KIND_R8
    use med_internalstate_mod , only : InternalState
    use med_phases_history_mod, only : med_phases_history_write_comp
    use med_map_mod           , only : med_map_field_packed
    use med_constants_mod     , only : dbug_flag => med_constants_dbug_flag
    use med_utils_mod         , only : chkerr    => med_utils_ChkErr
    use med_internalstate_mod , only : compocn, compatm, compice, complnd, compwav, coupling_mode
    use perf_mod              , only : t_startf, t_stopf

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(InternalState) :: is_local
    type(ESMF_Clock)    :: dClock
    character(len=*), parameter :: subname='(med_phases_post_atm)'
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    call t_startf('MED:'//subname)
    if (dbug_flag > 20) then
       call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)
    end if

    ! Get the internal state
    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (trim(coupling_mode) == 'access') then
      call med_phases_post_atm_custom_access(gcomp, rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    ! map atm to ocn
    if (is_local%wrap%med_coupling_active(compatm,compocn)) then
       call t_startf('MED:'//trim(subname)//' map_atm2ocn')
       call med_map_field_packed( &
            FBSrc=is_local%wrap%FBImp(compatm,compatm), &
            FBDst=is_local%wrap%FBImp(compatm,compocn), &
            FBFracSrc=is_local%wrap%FBFrac(compatm), &
            field_normOne=is_local%wrap%field_normOne(compatm,compocn,:), &
            packed_data=is_local%wrap%packed_data(compatm,compocn,:), &
            routehandles=is_local%wrap%RH(compatm,compocn,:), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call t_stopf('MED:'//trim(subname)//' map_atm2ocn')
    end if
    ! map atm->ice
    if (is_local%wrap%med_coupling_active(compatm,compice)) then
       call t_startf('MED:'//trim(subname)//' map_atm2ice')
       call med_map_field_packed( &
            FBSrc=is_local%wrap%FBImp(compatm,compatm), &
            FBDst=is_local%wrap%FBImp(compatm,compice), &
            FBFracSrc=is_local%wrap%FBFrac(compatm), &
            field_normOne=is_local%wrap%field_normOne(compatm,compice,:), &
            packed_data=is_local%wrap%packed_data(compatm,compice,:), &
            routehandles=is_local%wrap%RH(compatm,compice,:), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call t_stopf('MED:'//trim(subname)//' map_atm2ice')
    end if
    ! map atm->lnd
    if (is_local%wrap%med_coupling_active(compatm,complnd)) then
       call t_startf('MED:'//trim(subname)//' map_atm2lnd')
       call med_map_field_packed( &
            FBSrc=is_local%wrap%FBImp(compatm,compatm), &
            FBDst=is_local%wrap%FBImp(compatm,complnd), &
            FBFracSrc=is_local%wrap%FBFrac(compatm), &
            field_normOne=is_local%wrap%field_normOne(compatm,complnd,:), &
            packed_data=is_local%wrap%packed_data(compatm,complnd,:), &
            routehandles=is_local%wrap%RH(compatm,complnd,:), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call t_stopf('MED:'//trim(subname)//' map_atm2lnd')
    end if
    ! map atm->wav
    if (is_local%wrap%med_coupling_active(compatm,compwav)) then
       call t_startf('MED:'//trim(subname)//' map_atm2wav')
       call med_map_field_packed( &
            FBSrc=is_local%wrap%FBImp(compatm,compatm), &
            FBDst=is_local%wrap%FBImp(compatm,compwav), &
            FBFracSrc=is_local%wrap%FBFrac(compatm), &
            field_normOne=is_local%wrap%field_normOne(compatm,compwav,:), &
            packed_data=is_local%wrap%packed_data(compatm,compwav,:), &
            routehandles=is_local%wrap%RH(compatm,compwav,:), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call t_stopf('MED:'//trim(subname)//' map_atm2wav')
    end if

    ! Write atm inst, avg or aux if requested in mediator attributes
    call NUOPC_MediatorGet(gcomp, driverClock=dClock, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (ESMF_ClockIsCreated(dclock)) then
       call med_phases_history_write_comp(gcomp, compatm, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    if (dbug_flag > 20) then
       call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)
    end if
    call t_stopf('MED:'//subname)

  end subroutine med_phases_post_atm

  subroutine med_phases_post_atm_custom_access(gcomp, rc)
   use med_kind_mod          , only : CX=>SHR_KIND_CX, CS=>SHR_KIND_CS, CL=>SHR_KIND_CL, R8=>SHR_KIND_R8
   use med_internalstate_mod , only : compocn, compatm, compice, coupling_mode
   use med_internalstate_mod , only : InternalState, maintask, logunit
   use ESMF , only : ESMF_GridComp, ESMF_FieldBundleGet, ESMF_FieldCreate
   use ESMF , only : ESMF_FieldGet, ESMF_Field, ESMF_Mesh, ESMF_MESHLOC_ELEMENT, ESMF_TYPEKIND_R8
   use ESMF , only : ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_SUCCESS
   use ESMF , only : ESMF_FAILURE,  ESMF_LOGMSG_ERROR
   use med_map_mod           , only : med_map_field
   use med_internalstate_mod , only : mapconsf
   use med_constants_mod     , only : dbug_flag => med_constants_dbug_flag
   use med_methods_mod       , only : FB_GetFldPtr  => med_methods_FB_GetFldPtr

   ! input/output variables
   type(ESMF_GridComp)  :: gcomp
   integer, intent(out) :: rc

   ! local variables
   type(InternalState) :: is_local
   real(R8), pointer   :: ice_frac_cat_ptr(:, :), ice_flux_cat_ptr(:, :)
   type(ESMF_Field) :: ice_frac_cat, ice_flux_cat
   integer             :: lsize1, lsize2, i, j, n
   character(len=*), parameter    :: subname='(med_phases_post_atm_custom_access)'
   character(len=CS) :: fld_names(4)
   !---------------------------------------

   rc = ESMF_SUCCESS

   ! Get the internal state
   nullify(is_local%wrap)
   call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
   
   call ESMF_FieldBundleGet(is_local%wrap%FBImp(compice, compatm), fieldName='ia_aicen', field=ice_frac_cat, rc=rc)
   call ESMF_FieldGet(ice_frac_cat, farrayptr=ice_frac_cat_ptr)

   lsize1 = size(ice_frac_cat_ptr, dim=1)
   lsize2 = size(ice_frac_cat_ptr, dim=2)

   fld_names = [character(len=CS) :: &
   'topmelt', &
   'botmelt', &
   'sublim', &
   'pen_rad']
   
   do n = 1,size(fld_names)
      
      call ESMF_FieldBundleGet(is_local%wrap%FBImp(compatm, compatm), fieldName=trim(fld_names(n)), field=ice_flux_cat, rc=rc)
      call ESMF_FieldGet(ice_flux_cat, farrayptr=ice_flux_cat_ptr)

      do j = 1,lsize2
         do i = 1,lsize1
            if (ice_frac_cat_ptr(i, j) > 0.0) then
               ice_flux_cat_ptr(i, j) = ice_flux_cat_ptr(i, j) / ice_frac_cat_ptr(i, j)
            end if
         end do
      end do
      
   end do

  end subroutine med_phases_post_atm_custom_access

end module med_phases_post_atm_mod
