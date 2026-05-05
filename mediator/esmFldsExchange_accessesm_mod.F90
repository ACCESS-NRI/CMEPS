module esmFldsExchange_accessesm_mod

    use ESMF
    use NUOPC
    use med_utils_mod         , only : chkerr => med_utils_chkerr
    use med_kind_mod          , only : CX=>SHR_KIND_CX, CS=>SHR_KIND_CS, CL=>SHR_KIND_CL, R8=>SHR_KIND_R8
    use med_internalstate_mod , only : compmed, compatm, compocn, compice, ncomps
    use med_internalstate_mod , only : coupling_mode
    use esmFlds               , only : med_fldList_type
    use esmFlds               , only : addfld_to => med_fldList_addfld_to
    use esmFlds               , only : addmrg_to => med_fldList_addmrg_to
    use esmFlds               , only : addfld_from => med_fldList_addfld_from
    use esmFlds               , only : addmap_from => med_fldList_addmap_from

    !---------------------------------------------------------------------
    ! This is a mediator specific routine that determines ALL possible
    ! fields exchanged between components and their associated routing,
    ! mapping and merging

    ! In general, first order conservative remapping is used for flux fields
    ! while bilinear remapping is applied to state fields.
    ! Velocity and stress fields are remapped using a higher order patch method.
    ! Some sea-ice related fields are conservatively remapped with weighting by ice fraction.
    !---------------------------------------------------------------------

    implicit none
    public

    public :: esmFldsExchange_accessesm

    character(*), parameter :: u_FILE_u = &
         __FILE__

  !===============================================================================
  contains
  !===============================================================================

    subroutine esmFldsExchange_accessesm(gcomp, phase, rc)

      ! input/output parameters:
      type(ESMF_GridComp)              :: gcomp
      character(len=*) , intent(in)    :: phase
      integer          , intent(inout) :: rc

      ! local variables:
      character(len=*) , parameter   :: subname='(esmFldsExchange_accessesm)'
      !--------------------------------------

      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
      rc = ESMF_SUCCESS

      if (phase == 'advertise') then
        call esmFldsExchange_accessesm_advt(gcomp, phase, rc)
        if (chkerr(rc,__LINE__,u_FILE_u)) return
      elseif (phase == 'initialize') then
        call esmFldsExchange_accessesm_init(gcomp, phase, rc)
        if (chkerr(rc,__LINE__,u_FILE_u)) return
      else
        call ESMF_LogSetError(ESMF_FAILURE, &
           msg=trim(subname)//": Phase is set to "//trim(phase), &
           line=__LINE__, file=__FILE__, rcToReturn=rc)
        return  ! bail out
      endif

      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)

    end subroutine esmFldsExchange_accessesm

    !-----------------------------------------------------------------------------

    subroutine esmFldsExchange_accessesm_advt(gcomp, phase, rc)

      ! input/output parameters:
      type(ESMF_GridComp)              :: gcomp
      character(len=*) , intent(in)    :: phase
      integer          , intent(inout) :: rc

      ! local variables:
      integer             :: num, i, n
      logical             :: isPresent
      character(len=CL)   :: cvalue
      character(len=CS)   :: name, fldname
      character(len=CS)   :: fldname1, fldname2
      character(len=CS), allocatable :: flds(:)
      character(len=CS), allocatable :: S_flds(:)
      character(len=CS), allocatable :: F_flds(:,:)
      character(len=CS), allocatable :: suffix(:)
      character(len=*) , parameter   :: subname='(esmFldsExchange_accessesm_advt)'
      !--------------------------------------

      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
      rc = ESMF_SUCCESS

      !=====================================================================
      ! scalar information
      !=====================================================================

      call NUOPC_CompAttributeGet(gcomp, name='ScalarFieldName', &
         isPresent=isPresent, rc=rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return
      if (isPresent) then
         call NUOPC_CompAttributeGet(gcomp, name="ScalarFieldName", &
            value=cvalue, rc=rc)
         if (ChkErr(rc,__LINE__,u_FILE_u)) return
         do n = 1,ncomps
            call addfld_from(n, trim(cvalue))
            call addfld_to(n, trim(cvalue))
         end do
      end if


      !=====================================================================
      ! FIELDS TO MEDIATOR component (for fractions and atm/ocn flux calculation)
      !=====================================================================

      !----------------------------------------------------------
      ! to med: masks from components
      !----------------------------------------------------------
      call addfld_from(compocn, 'So_omask')
      call addfld_from(compice, 'Si_imask')

      !=====================================================================
      ! FIELDS TO ATMOSPHERE
      !=====================================================================

      call addfld_to(compatm, 'So_ofrac')
      call addfld_to(compatm, 'Si_ifrac')

      ! ---------------------------------------------------------------------
      ! to atm: from ocn
      ! ---------------------------------------------------------------------
      allocate(S_flds(3))
      S_flds = [character(len=CS) :: 'So_t', 'So_u', 'So_v']
      do n = 1,size(S_flds)
        fldname = trim(S_flds(n))
        call addfld_from(compocn, trim(fldname))
        call addfld_to(compatm, trim(fldname))
      end do
      deallocate(S_flds)

      ! ---------------------------------------------------------------------
      ! to atm: from ice
      ! ---------------------------------------------------------------------
      allocate(S_flds(9))
      S_flds = [character(len=CS) :: 'Si_t', &
               'Si_ifrac_n', &
               'Si_vsno_n', &
               'Si_vice_n', &
               'Si_topt', &
               'Si_topk', &
               'Si_pndf_n', &
               'Si_pndt_n', &
               'Si_Tf' &
               ]
      do n = 1,size(S_flds)
        fldname = trim(S_flds(n))
        call addfld_from(compice, trim(fldname))
        call addfld_to(compatm, trim(fldname))
      end do
      deallocate(S_flds)

      !=====================================================================
      ! FIELDS TO OCEAN (compocn)
      !=====================================================================

      ! ---------------------------------------------------------------------
      ! to ocn: state fields from atm
      ! ---------------------------------------------------------------------
      allocate(S_flds(2))
      S_flds = [character(len=CS) :: 'Sa_pslv', &
               'So_duu10n' ]
      do n = 1,size(S_flds)
         fldname = trim(S_flds(n))
         call addfld_from(compatm, trim(fldname))
         call addfld_to(compocn, trim(fldname))
      end do
      deallocate(S_flds)

      ! ---------------------------------------------------------------------
      ! to ocn: flux fields from atm
      ! ---------------------------------------------------------------------

      allocate(F_flds(11, 2))
      F_flds(1,:) = [character(len=CS) :: 'Faxa_taux ', 'Foxx_taux']
      F_flds(2,:) = [character(len=CS) :: 'Faxa_tauy ', 'Foxx_tauy']
      F_flds(3,:) = [character(len=CS) :: 'Faoa_sen', 'Foxx_sen']
      F_flds(4,:) = [character(len=CS) :: 'Faoa_evap', 'Foxx_evap']
      F_flds(5,:) = [character(len=CS) :: 'Faoa_lwnet', 'Foxx_lwnet']
      F_flds(6,:) = [character(len=CS) :: 'Faoa_swnet_vdr', 'Foxx_swnet_vdr']
      F_flds(7,:) = [character(len=CS) :: 'Faoa_swnet_vdf', 'Foxx_swnet_vdf']
      F_flds(8,:) = [character(len=CS) :: 'Faoa_swnet_idr', 'Foxx_swnet_idr']
      F_flds(9,:) = [character(len=CS) :: 'Faoa_swnet_idf', 'Foxx_swnet_idf']
      F_flds(10,:) = [character(len=CS) :: 'Faoa_rofl', 'Foxx_rofl']  ! mean runoff rate (liquid)
      F_flds(11,:) = [character(len=CS) :: 'Faoa_rofi', 'Foxx_rofi']  ! mean runnof rate (frozen)

      do n = 1,size(F_flds,1)
         fldname1 = trim(F_flds(n,1))
         fldname2 = trim(F_flds(n,2))
         call addfld_from(compatm, trim(fldname1))
         call addfld_to(compocn, trim(fldname2))
      end do
      deallocate(F_flds)


      ! ---------------------------------------------------------------------
      ! to ocn: fields from ice
      ! ---------------------------------------------------------------------

      allocate(F_flds(6, 2))
      F_flds(1,:) = [character(len=CS) :: 'Fioi_salt', 'Fioi_salt'] ! salt flux sea-ice to ocean
      F_flds(2,:) = [character(len=CS) :: 'Si_ifrac', 'Si_ifrac'] ! ice_fraction
      F_flds(3,:) = [character(len=CS) :: 'Fioi_meltw', 'Fioi_meltw'] ! freshwater flux sea-ice to ocean
      F_flds(4,:) = [character(len=CS) :: 'Fioi_melth', 'Fioi_melth'] ! heat flux sea-ice to ocean
      F_flds(5,:) = [character(len=CS) :: 'Fioi_taux', 'Foxx_taux']
      F_flds(6,:) = [character(len=CS) :: 'Fioi_tauy', 'Foxx_tauy'] ! surface stress sea-ice to ocean
      do n = 1,size(F_flds,1)
         fldname1 = trim(F_flds(n,1))
         fldname2 = trim(F_flds(n,2))
         call addfld_from(compice, trim(fldname1))
         call addfld_to(compocn, trim(fldname2))
      end do
      deallocate(F_flds)

      !=====================================================================
      ! FIELDS TO ICE (compice)
      !=====================================================================

      ! ---------------------------------------------------------------------
      ! to ice: fields from ocn
      ! ---------------------------------------------------------------------

      allocate(S_flds(7))
      S_flds = [character(len=CS) :: 'So_dhdx', &
                 'So_dhdy', &
                 'So_t', &
                 'So_s', &
                 'So_u', &
                 'So_v', &
                 'Fioo_q' ]
      do n = 1,size(S_flds)
         fldname = trim(S_flds(n))
         call addfld_from(compocn, trim(fldname))
         call addfld_to(compice, trim(fldname))
      end do
      deallocate(S_flds)

      ! ---------------------------------------------------------------------
      ! to ice: fields from atm
      ! ---------------------------------------------------------------------

      allocate(F_flds(5, 2))
      F_flds(1,:) = [character(len=CS) :: 'Faxa_swpen_n', 'Faxa_swpen_n']
      F_flds(2,:) = [character(len=CS) :: 'Faxa_melthtop_n', 'Faxa_melthtop_n']
      F_flds(3,:) = [character(len=CS) :: 'Faxa_condtop_n', 'Faxa_condtop_n']
      F_flds(4,:) = [character(len=CS) :: 'Sa_tskn_n', 'Sa_tskn_n']
      F_flds(5,:) = [character(len=CS) :: 'Faxa_sublim_n', 'Faxa_sublim_n']
      do n = 1,size(F_flds,1)
         fldname1 = trim(F_flds(n,1))
         fldname2 = trim(F_flds(n,2))
         call addfld_from(compatm, trim(fldname1))
         call addfld_to(compice, trim(fldname2))
      end do
      deallocate(F_flds)

      ! ---------------------------------------------------------------------
      ! precipitation
      ! ---------------------------------------------------------------------

      call addfld_from(compatm, 'Faxa_rainc')
      call addfld_from(compatm, 'Faxa_snowc')
      call addfld_from(compatm, 'Faxa_rainl')
      call addfld_from(compatm, 'Faxa_snowl')

      call addfld_to(compocn, 'Faxa_rain')
      call addfld_to(compocn, 'Faxa_snow')

      call addfld_to(compice, 'Faxa_rain')
      call addfld_to(compice, 'Faxa_snow')

      ! ---------------------------------------------------------------------
      ! atm/ice wind stress
      ! ---------------------------------------------------------------------

      call addfld_to(compice, 'Faia_taux')
      call addfld_to(compice, 'Faia_tauy')

      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)

    end subroutine esmFldsExchange_accessesm_advt

    !-----------------------------------------------------------------------------

    subroutine esmFldsExchange_accessesm_init(gcomp, phase, rc)

      use med_methods_mod       , only : fldchk => med_methods_FB_FldChk
      use med_internalstate_mod , only : InternalState
      use med_internalstate_mod , only : mapbilnr, mapconsf, mapconsd, mappatch
      use med_internalstate_mod , only : mapfcopy, mapnstod, mapnstod_consd
      use med_internalstate_mod , only : mapfillv_bilnr
      use med_internalstate_mod , only : mapnstod_consf

      ! input/output parameters:
      type(ESMF_GridComp)              :: gcomp
      character(len=*) , intent(in)    :: phase
      integer          , intent(inout) :: rc

      ! local variables:
      type(InternalState) :: is_local
      integer             :: num, i, n
      integer             :: n1, n2, n3, n4
      character(len=CL)   :: cvalue
      character(len=CS)   :: name, fldname
      character(len=CS)   :: fldname1, fldname2
      character(len=CS), allocatable :: flds(:)
      character(len=CS), allocatable :: S_flds(:)
      character(len=CS), allocatable :: F_flds(:,:)
      character(len=CS), allocatable :: suffix(:)
      character(len=*) , parameter   :: subname='(esmFldsExchange_accessesm_init)'
      !--------------------------------------

      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
      rc = ESMF_SUCCESS

      !---------------------------------------
      ! Get the internal state
      !---------------------------------------
      nullify(is_local%wrap)
      call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return

      !=====================================================================
      ! FIELDS TO ATMOSPHERE
      !=====================================================================

      fldname = trim('So_t')
      if (fldchk(is_local%wrap%FBExp(compatm), trim(fldname), rc=rc) .and. &
         fldchk(is_local%wrap%FBImp(compocn, compocn), trim(fldname), rc=rc) &
        ) then
        call addmap_from(compocn, trim(fldname), compatm, mapbilnr, 'one', 'unset')
        call addmrg_to(compatm, trim(fldname), mrg_from=compocn, mrg_fld=trim(fldname), mrg_type='copy')
     end if

      allocate(S_flds(2))
      S_flds = [character(len=CS) :: 'So_u', 'So_v']
      do n = 1,size(S_flds)
         fldname = trim(S_flds(n))
         if (fldchk(is_local%wrap%FBExp(compatm), trim(fldname), rc=rc) .and. &
             fldchk(is_local%wrap%FBImp(compocn, compocn), trim(fldname), rc=rc) &
            ) then
            call addmap_from(compocn, trim(fldname), compatm, mappatch, 'one', 'unset')
            call addmrg_to(compatm, trim(fldname), mrg_from=compocn, mrg_fld=trim(fldname), mrg_type='copy')
         end if
      end do
      deallocate(S_flds)

      allocate(S_flds(9))
      S_flds = [character(len=CS) :: 'Si_t', &
               'Si_ifrac_n', &
               'Si_vsno_n', &
               'Si_vice_n', &
               'Si_topt', &
               'Si_topk', &
               'Si_pndf_n', &
               'Si_pndt_n', &
               'Si_Tf' &
               ]
      do n = 1,size(S_flds)
        fldname = trim(S_flds(n))
        if (fldchk(is_local%wrap%FBExp(compatm), trim(fldname), rc=rc) .and. &
             fldchk(is_local%wrap%FBImp(compice, compice), trim(fldname), rc=rc) &
           ) then
           call addmap_from(compice, trim(fldname), compatm, mapconsf, 'none', 'unset')
           call addmrg_to(compatm, trim(fldname), mrg_from=compice, mrg_fld=trim(fldname), mrg_type='copy')
        end if
      end do
      deallocate(S_flds)

      !=====================================================================
      ! FIELDS TO OCEAN (compocn)
      !=====================================================================

      ! ---------------------------------------------------------------------
      ! to ocn: state fields
      ! ---------------------------------------------------------------------
      allocate(S_flds(2))
      S_flds = [character(len=CS) :: 'Sa_pslv', &
               'So_duu10n' ]
      do n = 1,size(S_flds)
         fldname = trim(S_flds(n))
         if (fldchk(is_local%wrap%FBExp(compocn), trim(fldname), rc=rc) .and. &
             fldchk(is_local%wrap%FBImp(compatm, compatm), trim(fldname), rc=rc) &
            ) then

            call addmap_from(compatm, trim(fldname), compocn, mapbilnr, 'one', 'unset')
            call addmrg_to(compocn, trim(fldname), mrg_from=compatm, mrg_fld=trim(fldname), mrg_type='copy')

         end if
      end do
      deallocate(S_flds)

      ! ---------------------------------------------------------------------
      ! to ocn: flux fields from atm
      ! ---------------------------------------------------------------------

      allocate(F_flds(9, 2))
      F_flds(1,:) = [character(len=CS) :: 'Faoa_sen', 'Foxx_sen']
      F_flds(2,:) = [character(len=CS) :: 'Faoa_evap', 'Foxx_evap']
      F_flds(3,:) = [character(len=CS) :: 'Faoa_lwnet', 'Foxx_lwnet']
      F_flds(4,:) = [character(len=CS) :: 'Faoa_swnet_vdr', 'Foxx_swnet_vdr']
      F_flds(5,:) = [character(len=CS) :: 'Faoa_swnet_vdf', 'Foxx_swnet_vdf']
      F_flds(6,:) = [character(len=CS) :: 'Faoa_swnet_idr', 'Foxx_swnet_idr']
      F_flds(7,:) = [character(len=CS) :: 'Faoa_swnet_idf', 'Foxx_swnet_idf']
      F_flds(8,:) = [character(len=CS) :: 'Faoa_rofl', 'Foxx_rofl']  ! mean runoff rate (liquid)
      F_flds(9,:) = [character(len=CS) :: 'Faoa_rofi', 'Foxx_rofi']  ! mean runnof rate (frozen)

      do n = 1,size(F_flds,1)
         fldname1 = trim(F_flds(n,1))
         fldname2 = trim(F_flds(n,2))
         if (fldchk(is_local%wrap%FBExp(compocn), trim(fldname2), rc=rc) .and. &
             fldchk(is_local%wrap%FBImp(compatm, compatm), trim(fldname1), rc=rc) &
            ) then
            call addmap_from(compatm, trim(fldname1), compocn, mapconsf, 'one', 'unset')
            call addmrg_to(compocn, trim(fldname2), mrg_from=compatm, mrg_fld=trim(fldname1), mrg_type='copy')
         end if
      end do
      deallocate(F_flds)

      ! precip
      if (fldchk(is_local%wrap%FBExp(compocn), trim('Faxa_rain'), rc=rc) .and. &
          fldchk(is_local%wrap%FBImp(compatm, compatm), trim('Faxa_rainc'),rc=rc) .and. &
          fldchk(is_local%wrap%FBImp(compatm, compatm), trim('Faxa_rainl'),rc=rc) &
         ) then
         call addmap_from(compatm, 'Faxa_rainc', compocn, mapconsf, 'one', 'unset')
         call addmap_from(compatm, 'Faxa_rainl', compocn, mapconsf, 'one', 'unset')
         call addmrg_to(compocn, 'Faxa_rain' , mrg_from=compatm, mrg_fld='Faxa_rainc:Faxa_rainl', &
               mrg_type='sum_with_weights', mrg_fracname='ofrac')
      end if

      if (fldchk(is_local%wrap%FBExp(compocn), trim('Faxa_snow'), rc=rc) .and. &
          fldchk(is_local%wrap%FBImp(compatm, compatm), trim('Faxa_snowc'),rc=rc) .and. &
          fldchk(is_local%wrap%FBImp(compatm, compatm), trim('Faxa_snowl'),rc=rc) &
         ) then
         call addmap_from(compatm, 'Faxa_snowc', compocn, mapconsf, 'one', 'unset')
         call addmap_from(compatm, 'Faxa_snowl', compocn, mapconsf, 'one', 'unset')
         call addmrg_to(compocn, 'Faxa_snow' , mrg_from=compatm, mrg_fld='Faxa_snowc:Faxa_snowl', &
               mrg_type='sum_with_weights', mrg_fracname='ofrac')
      end if

      ! ---------------------------------------------------------------------
      ! to ocn: fields from ice
      ! ---------------------------------------------------------------------

      allocate(F_flds(4, 2))
      F_flds(1,:) = [character(len=CS) :: 'Fioi_salt', 'Fioi_salt']
      F_flds(2,:) = [character(len=CS) :: 'Fioi_meltw', 'Fioi_meltw']
      F_flds(3,:) = [character(len=CS) :: 'Fioi_melth', 'Fioi_melth']
      F_flds(4,:) = [character(len=CS) :: 'Si_ifrac', 'Si_ifrac']
      do n = 1,size(F_flds,1)
         fldname1 = trim(F_flds(n,1))
         fldname2 = trim(F_flds(n,2))
         if (fldchk(is_local%wrap%FBExp(compocn), trim(fldname2), rc=rc) .and. &
             fldchk(is_local%wrap%FBImp(compice, compice), trim(fldname1),rc=rc) &
            ) then
            call addmap_from(compice, trim(fldname1), compocn, mapfcopy, 'unset', 'unset')
            call addmrg_to(compocn, trim(fldname2), mrg_from=compice, mrg_fld=trim(fldname1), mrg_type='copy')
         end if
      end do
      deallocate(F_flds)

      ! momentum transfer
      if (fldchk(is_local%wrap%FBExp(compocn), trim('Foxx_taux'), rc=rc) .and. &
          fldchk(is_local%wrap%FBImp(compice, compice), trim('Fioi_taux'),rc=rc) &
         ) then
         call addmap_from(compice, trim('Fioi_taux'), compocn, mapfcopy, 'unset', 'unset')
         call addmrg_to(compocn, trim('Foxx_taux'), mrg_from=compice, mrg_fld=trim('Fioi_taux'), mrg_type='merge', mrg_fracname='ifrac')
      end if

      if (fldchk(is_local%wrap%FBExp(compocn), trim('Foxx_tauy'), rc=rc) .and. &
          fldchk(is_local%wrap%FBImp(compice, compice), trim('Fioi_tauy'),rc=rc) &
         ) then
         call addmap_from(compice, trim('Fioi_tauy'), compocn, mapfcopy, 'unset', 'unset')
         call addmrg_to(compocn, trim('Foxx_tauy'), mrg_from=compice, mrg_fld=trim('Fioi_tauy'), mrg_type='merge', mrg_fracname='ifrac')
      end if

      if (fldchk(is_local%wrap%FBExp(compocn), trim('Foxx_taux'), rc=rc) .and. &
          fldchk(is_local%wrap%FBImp(compatm, compatm), trim('Faxa_taux'),rc=rc) &
         ) then
         call addmap_from(compatm, trim('Faxa_taux'), compocn, mappatch, 'one', 'unset')
         call addmrg_to(compocn, trim('Foxx_taux'), mrg_from=compatm, mrg_fld=trim('Faxa_taux'), mrg_type='merge', mrg_fracname='ofrac')
      end if

      if (fldchk(is_local%wrap%FBExp(compocn), trim('Foxx_tauy'), rc=rc) .and. &
          fldchk(is_local%wrap%FBImp(compatm, compatm), trim('Faxa_tauy'),rc=rc) &
         ) then
         call addmap_from(compatm, trim('Faxa_tauy'), compocn, mappatch, 'one', 'unset')
         call addmrg_to(compocn, trim('Foxx_tauy'), mrg_from=compatm, mrg_fld=trim('Faxa_tauy'), mrg_type='merge', mrg_fracname='ofrac')
      end if

      !=====================================================================
      ! FIELDS TO ICE (compice)
      !=====================================================================

      ! ---------------------------------------------------------------------
      ! to ice: fields from ocn
      ! ---------------------------------------------------------------------

      allocate(S_flds(7))
      S_flds = [character(len=CS) :: 'So_dhdx', & ! sea_surface_slope_zonal
                 'So_dhdy', & ! sea_surface_slope_merid
                 'So_t ', & ! sea_surface_temperature
                 'So_s ', & ! sea surface salinity
                 'So_u', & ! ocean surface zonal current
                 'So_v', & ! ocean surface meridional current
                 'Fioo_q' ] ! Freezing/melting potential
      do n = 1,size(S_flds)
         fldname = trim(S_flds(n))
         if (fldchk(is_local%wrap%FBExp(compice),trim(fldname),rc=rc) .and. &
             fldchk(is_local%wrap%FBImp(compocn, compocn), trim(fldname),rc=rc) &
            ) then

            call addmap_from(compocn, trim(fldname), compice, mapfcopy, 'unset', 'unset')
            call addmrg_to(compice, trim(fldname), mrg_from=compocn, mrg_fld=trim(fldname), mrg_type='copy')

         end if
      end do
      deallocate(S_flds)

      ! ---------------------------------------------------------------------
      ! to ice: fields from atm
      ! ---------------------------------------------------------------------

      allocate(F_flds(5, 2))
      F_flds(1,:) = [character(len=CS) :: 'Faxa_swpen_n', 'Faxa_swpen_n']
      F_flds(2,:) = [character(len=CS) :: 'Faxa_melthtop_n', 'Faxa_melthtop_n']
      F_flds(3,:) = [character(len=CS) :: 'Faxa_condtop_n', 'Faxa_condtop_n']
      F_flds(4,:) = [character(len=CS) :: 'Sa_tskn_n', 'Sa_tskn_n']
      F_flds(5,:) = [character(len=CS) :: 'Faxa_sublim_n', 'Faxa_sublim_n']

      do n = 1,size(F_flds,1)
         fldname1 = trim(F_flds(n,1))
         fldname2 = trim(F_flds(n,2))
         if (fldchk(is_local%wrap%FBExp(compice), trim(fldname2), rc=rc) .and. &
             fldchk(is_local%wrap%FBImp(compatm, compatm), trim(fldname1), rc=rc) &
            ) then

            call addmap_from(compatm, trim(fldname1), compice, mapconsf, 'one', 'unset')
            call addmrg_to(compice, trim(fldname2), mrg_from=compatm, mrg_fld=trim(fldname1), mrg_type='copy')

         end if
      end do
      deallocate(F_flds)

      ! wind stress
      if (fldchk(is_local%wrap%FBExp(compice), trim('Faia_taux'), rc=rc) .and. &
          fldchk(is_local%wrap%FBImp(compatm, compatm), trim('Faxa_taux'),rc=rc) &
         ) then
         call addmap_from(compatm, trim('Faxa_taux'), compice, mappatch, 'one', 'unset')
         call addmrg_to(compice, trim('Faia_taux'), mrg_from=compatm, mrg_fld=trim('Faxa_taux'), mrg_type='copy')
      end if

      if (fldchk(is_local%wrap%FBExp(compice), trim('Faia_tauy'), rc=rc) .and. &
          fldchk(is_local%wrap%FBImp(compatm, compatm), trim('Faxa_tauy'),rc=rc) &
         ) then
         call addmap_from(compatm, trim('Faxa_tauy'), compice, mappatch, 'one', 'unset')
         call addmrg_to(compice, trim('Faia_tauy'), mrg_from=compatm, mrg_fld=trim('Faxa_tauy'), mrg_type='copy')
      end if

      ! precip
      if (fldchk(is_local%wrap%FBExp(compice), trim('Faxa_rain'), rc=rc) .and. &
          fldchk(is_local%wrap%FBImp(compatm, compatm), trim('Faxa_rainc'),rc=rc) .and. &
          fldchk(is_local%wrap%FBImp(compatm, compatm), trim('Faxa_rainl'),rc=rc) &
         ) then
         call addmap_from(compatm, 'Faxa_rainc', compice, mapconsf, 'one', 'unset')
         call addmap_from(compatm, 'Faxa_rainl', compice, mapconsf, 'one', 'unset')
         call addmrg_to(compice, 'Faxa_rain' , mrg_from=compatm, mrg_fld='Faxa_rainc:Faxa_rainl', &
               mrg_type='sum')
      end if

      if (fldchk(is_local%wrap%FBExp(compice), trim('Faxa_snow'), rc=rc) .and. &
          fldchk(is_local%wrap%FBImp(compatm, compatm), trim('Faxa_snowc'),rc=rc) .and. &
          fldchk(is_local%wrap%FBImp(compatm, compatm), trim('Faxa_snowl'),rc=rc) &
         ) then
         call addmap_from(compatm, 'Faxa_snowc', compice, mapconsf, 'one', 'unset')
         call addmap_from(compatm, 'Faxa_snowl', compice, mapconsf, 'one', 'unset')
         call addmrg_to(compice, 'Faxa_snow' , mrg_from=compatm, mrg_fld='Faxa_snowc:Faxa_snowl', &
               mrg_type='sum')
      end if

      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)

    end subroutine esmFldsExchange_accessesm_init

    !-----------------------------------------------------------------------------

  end module esmFldsExchange_accessesm_mod
