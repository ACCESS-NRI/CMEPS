module esmFldsExchange_access_mod

    use ESMF
    use NUOPC
    use med_utils_mod         , only : chkerr => med_utils_chkerr
    use med_kind_mod          , only : CX=>SHR_KIND_CX
    use med_kind_mod          , only : CS=>SHR_KIND_CS
    use med_kind_mod          , only : CL=>SHR_KIND_CL
    use med_kind_mod          , only : R8=>SHR_KIND_R8
    use med_internalstate_mod , only : compmed, compatm, compocn, compwav, compice
    use med_internalstate_mod , only : ncomps
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
    !---------------------------------------------------------------------

    implicit none
    public

    public :: esmFldsExchange_access

    character(*), parameter :: u_FILE_u = &
         __FILE__

  !===============================================================================
  contains
  !===============================================================================

    subroutine esmFldsExchange_access(gcomp, phase, rc)

      ! input/output parameters:
      type(ESMF_GridComp)              :: gcomp
      character(len=*) , intent(in)    :: phase
      integer          , intent(inout) :: rc

      ! local variables:
      character(len=*) , parameter   :: subname='(esmFldsExchange_access)'
      !--------------------------------------

      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
      rc = ESMF_SUCCESS

      if (phase == 'advertise') then
        call esmFldsExchange_access_advt(gcomp, phase, rc)
        if (chkerr(rc,__LINE__,u_FILE_u)) return
      elseif (phase == 'fieldcheck') then
        call esmFldsExchange_access_fchk(gcomp, phase, rc)
        if (chkerr(rc,__LINE__,u_FILE_u)) return
      elseif (phase == 'initialize') then
        call esmFldsExchange_access_init(gcomp, phase, rc)
        if (chkerr(rc,__LINE__,u_FILE_u)) return
      else
        call ESMF_LogSetError(ESMF_FAILURE, &
           msg=trim(subname)//": Phase is set to "//trim(phase), &
           line=__LINE__, file=__FILE__, rcToReturn=rc)
        return  ! bail out
      endif

      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)

    end subroutine esmFldsExchange_access

    !-----------------------------------------------------------------------------

    subroutine esmFldsExchange_access_advt(gcomp, phase, rc)

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
      character(len=*) , parameter   :: subname='(esmFldsExchange_access_advt)'
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
      S_flds = (/'So_t', 'So_u', 'So_v'/) ! sea_surface_temperature
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
      S_flds = (/'Si_t', &
                  'ia_aicen', &
                  'ia_snown', &
                  'ia_thikn', &
                  'ia_itopt', &
                  'ia_itopk', &
                  'ia_pndfn', &
                  'ia_pndtn', &
                  'sstfrz' &
               /) ! sea_surface_temperature
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
      ! to ocn: state fields
      ! ---------------------------------------------------------------------
      allocate(S_flds(2))
      S_flds = (/'Sa_pslv', & ! inst_zonal_wind_height10m
                  'So_duu10n' /) ! inst_temp_height_surface
      do n = 1,size(S_flds)
         fldname = trim(S_flds(n))
         call addfld_from(compatm, trim(fldname))
         call addfld_to(compocn, trim(fldname))
      end do
      deallocate(S_flds)

      ! ---------------------------------------------------------------------
      ! to ocn: flux fields
      ! ---------------------------------------------------------------------

      ! from atm
      allocate(F_flds(13, 2))
      F_flds(1,:) = (/'Faxa_taux ', 'Foxx_taux'/)
      F_flds(2,:) = (/'Faxa_tauy ', 'Foxx_tauy'/)
      F_flds(3,:) = (/'Foxx_sen', 'Foxx_sen'/)
      F_flds(4,:) = (/'Foxx_evap', 'Foxx_evap'/)
      F_flds(5,:) = (/'Foxx_lwnet', 'Foxx_lwnet'/)
      F_flds(6,:) = (/'Foxx_swnet_vdr', 'Foxx_swnet_vdr'/)
      F_flds(7,:) = (/'Foxx_swnet_vdf', 'Foxx_swnet_vdf'/)
      F_flds(8,:) = (/'Foxx_swnet_idr', 'Foxx_swnet_idr'/)
      F_flds(9,:) = (/'Foxx_swnet_idf', 'Foxx_swnet_idf'/)
      F_flds(10,:) = (/'Faxa_rainc', 'Faxa_rain'/)
      F_flds(11,:) = (/'Faxa_snowc', 'Faxa_snow'/)
      F_flds(12,:) = (/'Foxx_rofl', 'Foxx_rofl'/)  ! mean runoff rate (liquid)
      F_flds(13,:) = (/'Foxx_rofi', 'Foxx_rofi'/)  ! mean runnof rate (frozen)

      do n = 1,size(F_flds,1)
         fldname1 = trim(F_flds(n,1))
         fldname2 = trim(F_flds(n,2))
         call addfld_from(compatm, trim(fldname1))
         call addfld_to(compocn, trim(fldname2))
      end do
      deallocate(F_flds)

      call addfld_from(compatm, 'Faxa_rainc')
      call addfld_from(compatm, 'Faxa_snowc')

      ! from ice
      allocate(F_flds(6, 2))
      F_flds(1,:) = (/'Fioi_salt', 'Fioi_salt'/)
      F_flds(2,:) = (/'Si_ifrac', 'Si_ifrac'/) ! ice_fraction
      F_flds(3,:) = (/'Fioi_meltw', 'Fioi_meltw'/)
      F_flds(4,:) = (/'Fioi_melth', 'Fioi_melth'/) ! heat flux sea-ice to ocean
      F_flds(5,:) = (/'Fioi_taux', 'Foxx_taux'/)
      F_flds(6,:) = (/'Fioi_tauy', 'Foxx_tauy'/) ! heat flux sea-ice to ocean
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
      ! to ice: state fields
      ! ---------------------------------------------------------------------

      ! from atm
      allocate(S_flds(10))
      S_flds = (/'Sa_z', &
                  'Sa_u', &
                  'Sa_v', &
                  'Sa_shum', &
                  'Sa_tbot', &
                  'Sa_pbot', &
                  'Sa_dens', &
                  'Sa_ptem', &
                  'um_icesth', &
                  'um_icenth' /)
      do n = 1,size(S_flds)
         fldname = trim(S_flds(n))
         call addfld_from(compatm, trim(fldname))
         call addfld_to(compice, trim(fldname))
      end do
      deallocate(S_flds)

      ! from ocn
      allocate(S_flds(7))
      S_flds = (/'So_dhdx', & 
                 'So_dhdy', &
                 'So_t', & 
                 'So_s', & 
                 'So_u', & 
                 'So_v', & 
                 'Fioo_q' /) 
      do n = 1,size(S_flds)
         fldname = trim(S_flds(n))
         call addfld_from(compocn, trim(fldname))
         call addfld_to(compice, trim(fldname))
      end do
      deallocate(S_flds)

      ! ---------------------------------------------------------------------
      ! to ice: flux fields
      ! ---------------------------------------------------------------------

      allocate(F_flds(12, 2))
      F_flds(1,:) = (/'Faxa_swvdr ', 'Faxa_swvdr '/)
      F_flds(2,:) = (/'Faxa_swndr ', 'Faxa_swndr '/)
      F_flds(3,:) = (/'Faxa_swvdf', 'Faxa_swvdf'/)
      F_flds(4,:) = (/'Faxa_swndf', 'Faxa_swndf'/)
      F_flds(5,:) = (/'Faxa_lwdn', 'Faxa_lwdn'/)
      F_flds(6,:) = (/'pen_rad', 'pen_rad'/)
      F_flds(7,:) = (/'topmelt', 'topmelt'/)
      F_flds(8,:) = (/'botmelt', 'botmelt'/)
      F_flds(9,:) = (/'tstar_sice', 'tstar_sice'/)
      F_flds(10,:) = (/'sublim', 'sublim'/)
      F_flds(11,:) = (/'Foxx_sen', 'Foxx_sen'/)
      F_flds(12,:) = (/'Faxa_swdn', 'Faxa_swdn'/)
      do n = 1,size(F_flds,1)
         fldname1 = trim(F_flds(n,1))
         fldname2 = trim(F_flds(n,2))
         call addfld_from(compatm, trim(fldname1))
         call addfld_to(compice, trim(fldname2))
      end do
      deallocate(F_flds)

      call addfld_from(compatm, 'Faxa_rainc')
      call addfld_from(compatm, 'Faxa_snowc')
      call addfld_from(compatm, 'Faxa_rainl')
      call addfld_from(compatm, 'Faxa_snowl')

      call addfld_to(compice, 'Faxa_rain')
      call addfld_to(compice, 'Faxa_snow')
      
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)

    end subroutine esmFldsExchange_access_advt

    !-----------------------------------------------------------------------------

    subroutine esmFldsExchange_access_fchk(gcomp, phase, rc)

      use med_methods_mod       , only : fldchk => med_methods_FB_FldChk
      use med_internalstate_mod , only : InternalState

      ! input/output parameters:
      type(ESMF_GridComp)              :: gcomp
      character(len=*) , intent(in)    :: phase
      integer          , intent(inout) :: rc

      ! local variables:
      type(InternalState) :: is_local
      character(len=*) , parameter   :: subname='(esmFldsExchange_access_fchk)'
      !--------------------------------------

      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
      rc = ESMF_SUCCESS

      !---------------------------------------
      ! Get the internal state
      !---------------------------------------
      nullify(is_local%wrap)
      call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return

      if (fldchk(is_local%wrap%FBImp(compocn,compocn),'So_omask',rc=rc)) then
         call ESMF_LogWrite(trim(subname)//": Field connected "//"So_omask", &
            ESMF_LOGMSG_INFO)
      else
         call ESMF_LogSetError(ESMF_FAILURE, &
            msg=trim(subname)//": Field is not connected "//"So_omask", &
            line=__LINE__, file=__FILE__, rcToReturn=rc)
         return  ! bail out
      endif

      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)

    end subroutine esmFldsExchange_access_fchk

    !-----------------------------------------------------------------------------

    subroutine esmFldsExchange_access_init(gcomp, phase, rc)

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
      character(len=*) , parameter   :: subname='(esmFldsExchange_access_init)'
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

      ! ---------------------------------------------------------------------
      ! to atm: sea surface temperature
      ! ---------------------------------------------------------------------
      call addmap_from(compocn, 'So_t', compatm, mapconsf, 'ofrac', 'unset')
      call addmrg_to(compatm, 'So_t', mrg_from=compocn, mrg_fld='So_t', mrg_type='copy')
      call addmap_from(compocn, 'So_u', compatm, mapconsf, 'ofrac', 'unset')
      call addmrg_to(compatm, 'So_u', mrg_from=compocn, mrg_fld='So_u', mrg_type='copy')
      call addmap_from(compocn, 'So_v', compatm, mapconsf, 'ofrac', 'unset')
      call addmrg_to(compatm, 'So_v', mrg_from=compocn, mrg_fld='So_v', mrg_type='copy')

      call addmap_from(compice, 'Si_t', compatm, mapconsd, 'ifrac', 'unset')
      call addmrg_to(compatm, 'Si_t', mrg_from=compice, mrg_fld='Si_t', mrg_type='copy')

      call addmap_from(compice, 'sstfrz', compatm, mapconsf, 'none', 'unset')
      call addmrg_to(compatm, 'sstfrz', mrg_from=compice, mrg_fld='sstfrz', mrg_type='copy')
      
      allocate(S_flds(7))
      S_flds = (/'ia_aicen', &
                  'ia_snown', &
                  'ia_thikn', &
                  'ia_itopt', &
                  'ia_itopk', &
                  'ia_pndfn', &
                  'ia_pndtn'/) ! sea_surface_temperature
      do n = 1,size(S_flds)
        fldname = trim(S_flds(n))
        call addmap_from(compice, trim(fldname), compatm, mapconsf, 'none', 'unset')
        call addmrg_to(compatm, trim(fldname), mrg_from=compice, mrg_fld=trim(fldname), mrg_type='copy')
      end do
      deallocate(S_flds)
      ! call addmap(fldListFr(compice)%flds, 'Si_t', compatm, mapconsf, 'ifrac', 'unset')
      ! call addmrg(fldListTo(compatm)%flds, 'Si_t', mrg_from=compice, mrg_fld='Si_t', mrg_type='copy')

      !=====================================================================
      ! FIELDS TO OCEAN (compocn)
      !=====================================================================

      ! ---------------------------------------------------------------------
      ! to ocn: state fields
      ! ---------------------------------------------------------------------
      allocate(S_flds(2))
      S_flds = (/'Sa_pslv', & ! inst_zonal_wind_height10m
                  'So_duu10n' /) ! inst_temp_height_surface
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
      ! to ocn: flux fields
      ! ---------------------------------------------------------------------

      ! from atm
      allocate(F_flds(9, 2))
      F_flds(1,:) = (/'Foxx_sen', 'Foxx_sen'/)
      F_flds(2,:) = (/'Foxx_evap', 'Foxx_evap'/)
      F_flds(3,:) = (/'Foxx_lwnet', 'Foxx_lwnet'/)
      F_flds(4,:) = (/'Foxx_swnet_vdr', 'Foxx_swnet_vdr'/)
      F_flds(5,:) = (/'Foxx_swnet_vdf', 'Foxx_swnet_vdf'/)
      F_flds(6,:) = (/'Foxx_swnet_idr', 'Foxx_swnet_idr'/)
      F_flds(7,:) = (/'Foxx_swnet_idf', 'Foxx_swnet_idf'/)
      F_flds(8,:) = (/'Foxx_rofl', 'Foxx_rofl'/)  ! mean runoff rate (liquid)
      F_flds(9,:) = (/'Foxx_rofi', 'Foxx_rofi'/)  ! mean runnof rate (frozen)

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
      call addmap_from(compatm, 'Faxa_rainc', compocn, mapconsf, 'one', 'unset') ! TODO: weight by ocean fraction
      call addmap_from(compatm, 'Faxa_rainl', compocn, mapconsf, 'one', 'unset')
      call addmrg_to(compocn, 'Faxa_rain' , mrg_from=compatm, mrg_fld='Faxa_rainc:Faxa_rainl', &
               mrg_type='sum_with_weights', mrg_fracname='ofrac')
      
      call addmap_from(compatm, 'Faxa_snowc', compocn, mapconsf, 'one', 'unset')
      call addmap_from(compatm, 'Faxa_snowl', compocn, mapconsf, 'one', 'unset')
      call addmrg_to(compocn, 'Faxa_snow' , mrg_from=compatm, mrg_fld='Faxa_snowc:Faxa_snowl', &
               mrg_type='sum_with_weights', mrg_fracname='ofrac')
      
      ! from ice
      call addmap_from(compice, 'Si_ifrac', compocn, mapfcopy, 'unset', 'unset')
      call addmrg_to(compocn, 'Si_ifrac', mrg_from=compice, mrg_fld='Si_ifrac', mrg_type='copy')

      allocate(F_flds(3, 2))
      F_flds(1,:) = (/'Fioi_salt', 'Fioi_salt'/)
      F_flds(2,:) = (/'Fioi_meltw', 'Fioi_meltw'/)
      F_flds(3,:) = (/'Fioi_melth', 'Fioi_melth'/) ! heat flux sea-ice to ocean
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
      call addmap_from(compice, 'Fioi_taux', compocn, mapfcopy, 'unset', 'unset')
      call addmrg_to(compocn, 'Foxx_taux', mrg_from=compice, mrg_fld='Fioi_taux', mrg_type='merge', mrg_fracname='ifrac')
      call addmap_from(compatm, 'Faxa_taux', compocn, mapconsf, 'one', 'unset')
      call addmrg_to(compocn, 'Foxx_taux', mrg_from=compatm, mrg_fld='Faxa_taux', mrg_type='merge', mrg_fracname='ofrac')

      call addmap_from(compice, 'Fioi_tauy', compocn, mapfcopy, 'unset', 'unset')
      call addmrg_to(compocn, 'Foxx_tauy', mrg_from=compice, mrg_fld='Fioi_tauy', mrg_type='merge', mrg_fracname='ifrac')
      call addmap_from(compatm, 'Faxa_tauy', compocn, mapconsf, 'one', 'unset')
      call addmrg_to(compocn, 'Foxx_tauy', mrg_from=compatm, mrg_fld='Faxa_tauy', mrg_type='merge', mrg_fracname='ofrac')

      !=====================================================================
      ! FIELDS TO ICE (compice)
      !=====================================================================

      ! ---------------------------------------------------------------------
      ! to ice: state fields
      ! ---------------------------------------------------------------------

      ! from atm
      allocate(S_flds(10))
      S_flds = (/'Sa_z', &
                  'Sa_u', &
                  'Sa_v', &
                  'Sa_shum', &
                  'Sa_tbot', &
                  'Sa_pbot', &
                  'Sa_dens', &
                  'Sa_ptem', &
                  'um_icesth', &
                  'um_icenth' /)

      do n = 1,size(S_flds)
         fldname = trim(S_flds(n))
         if (fldchk(is_local%wrap%FBExp(compice), trim(fldname),rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(compatm, compatm), trim(fldname),rc=rc) &
            ) then

            call addmap_from(compatm, trim(fldname), compice, mapbilnr, 'one', 'unset')
            call addmrg_to(compice, trim(fldname), mrg_from=compatm, mrg_fld=trim(fldname), mrg_type='copy')

         end if
      end do
      deallocate(S_flds)

      ! from ocn
      allocate(S_flds(7))
      S_flds = (/'So_dhdx', & ! inst_zonal_wind_height10m
                 'So_dhdy', & ! inst_merid_wind_height10m
                 'So_t ', & ! inst_temp_height2m
                 'So_s ', & ! inst_spec_humid_height2m
                 'So_u', & ! Sa_pslv
                 'So_v', & ! Sa_pslv
                 'Fioo_q' /) ! inst_temp_height_surface
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
      ! to ice: flux fields
      ! ---------------------------------------------------------------------

      ! from atm
      allocate(F_flds(12, 2))
      F_flds(1,:) = (/'pen_rad', 'pen_rad'/)
      F_flds(2,:) = (/'topmelt', 'topmelt'/)
      F_flds(3,:) = (/'botmelt', 'botmelt'/)
      F_flds(4,:) = (/'tstar_sice', 'tstar_sice'/)
      F_flds(5,:) = (/'sublim', 'sublim'/)
      F_flds(6,:) = (/'Faxa_swvdr ', 'Faxa_swvdr '/)
      F_flds(7,:) = (/'Faxa_swndr ', 'Faxa_swndr '/)
      F_flds(8,:) = (/'Faxa_swvdf', 'Faxa_swvdf'/)
      F_flds(9,:) = (/'Faxa_swndf', 'Faxa_swndf'/)
      F_flds(10,:) = (/'Faxa_lwdn', 'Faxa_lwdn'/)
      F_flds(11,:) = (/'Foxx_sen', 'Foxx_sen'/)
      F_flds(12,:) = (/'Faxa_swdn', 'Faxa_swdn'/)

      do n = 1,size(F_flds,1)
         fldname1 = trim(F_flds(n,1))
         fldname2 = trim(F_flds(n,2))
         if (fldchk(is_local%wrap%FBExp(compice), trim(fldname2), rc=rc) .and. &
             fldchk(is_local%wrap%FBImp(compatm, compatm), trim(fldname1), rc=rc) &
            ) then

            call addmap_from(compatm, trim(fldname1), compice, mapconsf, 'one', 'unset') ! mapping with total ifrac, should use category fractions
            call addmrg_to(compice, trim(fldname2), mrg_from=compatm, mrg_fld=trim(fldname1), mrg_type='copy')

         end if
      end do
      deallocate(F_flds)

      ! precip
      call addmap_from(compatm, 'Faxa_rainc', compice, mapconsf, 'one', 'unset')
      call addmap_from(compatm, 'Faxa_rainl', compice, mapconsf, 'one', 'unset')
      call addmrg_to(compice, 'Faxa_rain' , mrg_from=compatm, mrg_fld='Faxa_rainc:Faxa_rainl', &
               mrg_type='sum')
      
      call addmap_from(compatm, 'Faxa_snowc', compice, mapconsf, 'one', 'unset')
      call addmap_from(compatm, 'Faxa_snowl', compice, mapconsf, 'one', 'unset')
      call addmrg_to(compice, 'Faxa_snow' , mrg_from=compatm, mrg_fld='Faxa_snowc:Faxa_snowl', &
               mrg_type='sum')

      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)

    end subroutine esmFldsExchange_access_init

    !-----------------------------------------------------------------------------

  end module esmFldsExchange_access_mod
