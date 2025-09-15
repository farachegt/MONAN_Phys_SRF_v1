module module_bl_shocmf
 use mpas_atmphys_constants
 use shoc, only: shoc_init, shoc_main, length_fac, mintke, linear_interp
 use mpas_pool_routines


 implicit none
 private
 public:: bl_shocmf_init, bl_shocmf_driver

 real(kind=RKIND) :: max_shoc_timestep = 150.0  ! Default SHOC timestep set in namelist

 contains

 subroutine calculate_surface_momentum_fluxes(ims, ime, jms, jme, ust_p, u_p, v_p, wspd_p, rho_p, uw_sfc_p, vw_sfc_p)
    ! Calculate the surface momentum flux
    ! Reference: Stull's Book - An Introduction to Boundary Layer Meteorology, chapter 7.4.1.
    integer, intent(in) :: ims, ime, jms, jme
    real(kind=RKIND), dimension(:,:), intent(in) :: ust_p, u_p, v_p, wspd_p, rho_p
    real(kind=RKIND), dimension(:,:), intent(out) :: uw_sfc_p, vw_sfc_p
    
    real(kind=RKIND), dimension(ime, jme) :: fmom

    integer :: i, j

    do i=ims,ime
    do j=jms,jme
        fmom(i,j) =  rho_p(i,j)*((ust_p(i,j)*ust_p(i,j))/max(wspd_p(i,j),0.001))
        uw_sfc_p(i,j) = -fmom(i,j)*u_p(i,j)/rho_p(i,j)
        vw_sfc_p(i,j) = -fmom(i,j)*v_p(i,j)/rho_p(i,j)
    enddo
    enddo
 end subroutine calculate_surface_momentum_fluxes


 subroutine calculate_dry_static_energy(ims, ime, kms, kme, jms, jme, t, z, gravity_in, phis, cp_in, dse)
    ! Calculate the dry static energy considering the terrain height.
    ! This is the calculation waited by the SHOCMF module. 
    ! Reference: https://glossary.ametsoc.org/wiki/Dry_static_energy
    integer, intent(in) :: ims, ime, kms, kme, jms, jme
    real(kind=RKIND), dimension(:,:,:), intent(in) :: t, z
    real(kind=RKIND), dimension(:,:), intent(in) :: phis
    real(kind=RKIND), intent(in) :: gravity_in, cp_in
    real(kind=RKIND), dimension(:,:,:), intent(out) :: dse

    integer :: i, j, k
    
    do i=ims,ime
    do k=kms,kme-1
    do j=jms,jme
        dse(i,k,j) = cp_in*t(i,k,j) + gravity_in*z(i,k,j) + phis(i,j)
    enddo
    enddo
    enddo
 end subroutine calculate_dry_static_energy


 subroutine calculate_presdiff(ims, ime, kms, kme, jms, jme, pres_levels, presdiff_levels)
    ! Calculate the pressure difference between the interfaces
    ! pres_levels is the pressure at the levels
    ! presdiff_levels is the pressure difference between the interfaces
    integer, intent(in) :: ims,ime,kms,kme,jms,jme
    real(kind=RKIND), dimension(:,:,:), intent(in) :: pres_levels
    real(kind=RKIND), dimension(:,:,:), intent(inout) :: presdiff_levels

    integer :: i, j, k

    ! As MPAS vertical coordinates has its origin at the surface,
    ! the pressure difference is calculated between the levels
    ! and the level above it.
    do j = jms,jme
    do k = kms,kme-1
    do i = ims,ime
        presdiff_levels(i,k,j) = pres_levels(i,k,j) - pres_levels(i,k+1,j)
    enddo
    enddo
    enddo
 end subroutine calculate_presdiff


 subroutine calculate_thetal(ims, ime, jms, jme, kms, kme, cp_in, xlv_in, theta, qc, inv_exner, thetal)
    ! Calculate liquid water potential temperature
    ! Reference: shoc_intr.F90 in E3SM code and
    ! https://glossary.ametsoc.org/wiki/Liquid_water_potential_temperature
    integer, intent(in) :: ims, ime, jms, jme, kms, kme
    real(kind=RKIND), intent(in) :: cp_in, xlv_in
    real(kind=RKIND), dimension(:,:,:), intent(in) :: theta, qc, inv_exner
    real(kind=RKIND), dimension(:,:,:), intent(inout) :: thetal
    
    integer :: i, j, k
    do j = jms,jme
    do k = kms,kme-1
    do i = ims,ime
        thetal(i,k,j) = theta(i,k,j) - inv_exner(i,k,j)*(xlv_in/cp_in)*qc(i,k,j)
        ! thetal(i,k,j) = theta(i,k,j)*(1. + xlv_in * qc(i,k,j)/(cp_in * max(t(i,k,j),253.)))**(-1.0) ! BRAMS
    enddo
    enddo
    enddo
 end subroutine calculate_thetal


 subroutine calculate_theta(ims, ime, jms, jme, kms, kme, cp_in, xlv_in, thetal, qc, inv_exner, theta)
    ! Calculate potential temperature from liquid water potential temperature.
    ! Inverse function of calculate_thetal.
    integer, intent(in) :: ims, ime, jms, jme, kms, kme
    real(kind=RKIND), intent(in) :: cp_in, xlv_in
    real(kind=RKIND), dimension(:,:,:), intent(in) :: thetal, qc, inv_exner
    real(kind=RKIND), dimension(:,:,:), intent(out) :: theta

    integer :: i, j, k
    do j = jms,jme
    do k = kms,kme-1
    do i = ims,ime
        theta(i,k,j) = thetal(i,k,j) + inv_exner(i,k,j)*(xlv_in/cp_in)*qc(i,k,j)
        ! theta(i,k,j) = thetal(i,k,j)*(1. + xlv_in * qc(i,k,j)/(cp_in * max(t(i,k,j),253.))) ! BRAMS, usar qc e t atualizado pelo shoc
    enddo
    enddo
    enddo
 end subroutine calculate_theta


 subroutine calculate_thetav(ims, ime, jms, jme, kms, kme, theta, qv, qc, theta_v)
    ! Calculate virtual potential temperature as SHOCMF expects it
    ! This is called before the SHOCMF call
    ! Reference: shoc_intr.F90 in E3SM code
    integer, intent(in) :: ims, ime, jms, jme, kms, kme
    real(kind=RKIND), dimension(:,:,:), intent(in) :: theta, qv, qc
    real(kind=RKIND), dimension(:,:,:), intent(out) :: theta_v

    integer :: i, j, k
    do j = jms,jme
    do k = kms,kme-1
    do i = ims,ime
        theta_v(i,k,j) = theta(i,k,j)*(1.0_RKIND + 0.608_RKIND*qv(i,k,j) - qc(i,k,j))
    enddo
    enddo
    enddo
 end subroutine calculate_thetav


 subroutine calculate_qv(ims, ime, jms, jme, kms, kme, qw, qc, qv)
    ! Calculate the water vapor mixing ratio
    integer, intent(in) :: ims, ime, jms, jme, kms, kme
    real(kind=RKIND), dimension(:,:,:), intent(in) :: qw, qc
    real(kind=RKIND), dimension(:,:,:), intent(out) :: qv
    
    integer :: i, j, k
    do j = jms,jme
    do k = kms,kme-1
    do i = ims,ime
        qv(i,k,j) = qw(i,k,j) - qc(i,k,j)
    enddo
    enddo
    enddo
 end subroutine calculate_qv


 subroutine calculate_qw(ims, ime, jms, jme, kms, kme, qc, qv, qw)
    ! Calculate total water mixing ratio
    integer, intent(in) :: ims, ime, jms, jme, kms, kme
    real(kind=RKIND), dimension(:,:,:), intent(in) :: qc, qv
    real(kind=RKIND), dimension(:,:,:), intent(out) :: qw

    integer :: i, j, k
    do j = jms,jme
    do k = kms,kme-1
    do i = ims,ime
        qw(i,k,j) = qc(i,k,j) + qv(i,k,j)
    enddo
    enddo
    enddo
 end subroutine calculate_qw


 subroutine flip_vertical_coordinates(&
    kms, kme, rino_p, shoc_mix_p, isotropy_p, tke_dissip_p, tke_prod_b_p, &
    tke_prod_sh_p, w_sec_p, thl_sec_p, qw_sec_p, qwthl_sec_p, &
    wthl_sec_p, wqw_sec_p, wtke_sec_p, uw_sec_p, vw_sec_p, w3_p, &
    wqls_sec_p, brunt_p, shoc_ql2_p, mf_dry_a_p, mf_moist_a_p, &
    mf_dry_w_p, mf_moist_w_p, mf_dry_qt_p, mf_moist_qt_p, mf_dry_thl_p, &
    mf_moist_thl_p, mf_dry_u_p, mf_moist_u_p, mf_dry_v_p, mf_moist_v_p, &
    mf_moist_qc_p, mf_thlflx_p, mf_qtflx_p, mf_ae_p, mf_aw_p, &
    mf_awthv_p, mf_awthl_p, mf_awqt_p, mf_awql_p, mf_awqi_p, mf_awu_p, &
    mf_awv_p, cfl_mf_p, wthl_sec_ed_p, wthl_sec_mf_p, wqw_sec_ed_p, &
    wqw_sec_mf_p, dse_p, qw_p, u_p, v_p, wthv_sec_p, tke_shmf_p, kzh_p, &
    kzm_p, qc_p, cldfrabl_p, thv_p, height_interface_p, &
    height_mid_p, pres_p, pres2_p, w_p, inv_exner_p, &
    thetal_p, presdiff_p &
  )
    integer, intent(in) :: kms, kme
    real(kind=RKIND),dimension(:,:,:), intent(inout) :: rino_p, shoc_mix_p, isotropy_p, &
                                                        tke_dissip_p, tke_prod_b_p, &
                                                        tke_prod_sh_p, w_sec_p, thl_sec_p, qw_sec_p, &
                                                        qwthl_sec_p, wthl_sec_p, wqw_sec_p, wtke_sec_p, &
                                                        uw_sec_p, vw_sec_p, w3_p, wqls_sec_p, brunt_p, shoc_ql2_p, &
                                                        mf_dry_a_p, mf_moist_a_p, mf_dry_w_p, mf_moist_w_p, &
                                                        mf_dry_qt_p, mf_moist_qt_p, mf_dry_thl_p, mf_moist_thl_p, &
                                                        mf_dry_u_p, mf_moist_u_p, mf_dry_v_p, mf_moist_v_p, &
                                                        mf_moist_qc_p, mf_thlflx_p, mf_qtflx_p, mf_ae_p, &
                                                        mf_aw_p, mf_awthv_p, mf_awthl_p, mf_awqt_p, mf_awql_p, &
                                                        mf_awqi_p, mf_awu_p, mf_awv_p, cfl_mf_p, wthl_sec_ed_p, &
                                                        wthl_sec_mf_p, wqw_sec_ed_p, wqw_sec_mf_p, dse_p, qw_p, u_p, &
                                                        v_p, wthv_sec_p, tke_shmf_p, kzh_p, kzm_p, qc_p, cldfrabl_p, &
                                                        thv_p, height_interface_p, height_mid_p, &
                                                        pres_p, pres2_p, w_p, inv_exner_p, thetal_p, presdiff_p

    ! Flip the vertical coordinates of the SHOC-MF variables
    ! This is necessary because SHOC-MF expects the vertical
    ! coordinate starting from the top of the atmosphere.
    ! This subroutine is called before and after the SHOC-MF call
    thl_sec_p          = thl_sec_p(:, kme:kms:-1, :)
    qw_sec_p           = qw_sec_p(:, kme:kms:-1, :)
    qwthl_sec_p        = qwthl_sec_p(:, kme:kms:-1, :)
    wthl_sec_p         = wthl_sec_p(:, kme:kms:-1, :)
    wqw_sec_p          = wqw_sec_p(:, kme:kms:-1, :)
    wtke_sec_p         = wtke_sec_p(:, kme:kms:-1, :)
    uw_sec_p           = uw_sec_p(:, kme:kms:-1, :)
    vw_sec_p           = vw_sec_p(:, kme:kms:-1, :)
    w3_p               = w3_p(:, kme:kms:-1, :)
    mf_dry_a_p         = mf_dry_a_p(:, kme:kms:-1, :)
    mf_moist_a_p       = mf_moist_a_p(:, kme:kms:-1, :)
    mf_dry_w_p         = mf_dry_w_p(:, kme:kms:-1, :)
    mf_moist_w_p       = mf_moist_w_p(:, kme:kms:-1, :)
    mf_dry_qt_p        = mf_dry_qt_p(:, kme:kms:-1, :)
    mf_moist_qt_p      = mf_moist_qt_p(:, kme:kms:-1, :)
    mf_dry_thl_p       = mf_dry_thl_p(:, kme:kms:-1, :)
    mf_moist_thl_p     = mf_moist_thl_p(:, kme:kms:-1, :)
    mf_dry_u_p         = mf_dry_u_p(:, kme:kms:-1, :)
    mf_moist_u_p       = mf_moist_u_p(:, kme:kms:-1, :)
    mf_dry_v_p         = mf_dry_v_p(:, kme:kms:-1, :)
    mf_moist_v_p       = mf_moist_v_p(:, kme:kms:-1, :)
    mf_moist_qc_p      = mf_moist_qc_p(:, kme:kms:-1, :)
    mf_thlflx_p        = mf_thlflx_p(:, kme:kms:-1, :)
    mf_qtflx_p         = mf_qtflx_p(:, kme:kms:-1, :)
    mf_ae_p            = mf_ae_p(:, kme:kms:-1, :)
    mf_aw_p            = mf_aw_p(:, kme:kms:-1, :)
    mf_awthv_p         = mf_awthv_p(:, kme:kms:-1, :)
    mf_awthl_p         = mf_awthl_p(:, kme:kms:-1, :)
    mf_awqt_p          = mf_awqt_p(:, kme:kms:-1, :)
    mf_awql_p          = mf_awql_p(:, kme:kms:-1, :)
    mf_awqi_p          = mf_awqi_p(:, kme:kms:-1, :)
    mf_awu_p           = mf_awu_p(:, kme:kms:-1, :)
    mf_awv_p           = mf_awv_p(:, kme:kms:-1, :)
    cfl_mf_p           = cfl_mf_p(:, kme:kms:-1, :)
    wthl_sec_ed_p      = wthl_sec_ed_p(:, kme:kms:-1, :)
    wthl_sec_mf_p      = wthl_sec_mf_p(:, kme:kms:-1, :)
    wqw_sec_ed_p       = wqw_sec_ed_p(:, kme:kms:-1, :)
    wqw_sec_mf_p       = wqw_sec_mf_p(:, kme:kms:-1, :)
    height_interface_p = height_interface_p(:, kme:kms:-1, :)
    pres2_p            = pres2_p(:, kme:kms:-1, :)

    ! All variables were allocated with kme levels, but some
    ! of them are only kme-1 levels (thermodynamic ones).
    ! So, in these cases we need to flip only the vertical
    ! coordinates of the kme-1 levels
    rino_p(:,kms:kme-1,:)        = rino_p(:,kme-1:kms:-1,:)
    thv_p(:,kms:kme-1,:)         = thv_p(:,kme-1:kms:-1,:)
    shoc_mix_p(:,kms:kme-1,:)    = shoc_mix_p(:,kme-1:kms:-1,:)
    isotropy_p(:,kms:kme-1,:)    = isotropy_p(:,kme-1:kms:-1,:)
    tke_dissip_p(:,kms:kme-1,:)  = tke_dissip_p(:,kme-1:kms:-1,:)
    tke_prod_b_p(:,kms:kme-1,:)  = tke_prod_b_p(:,kme-1:kms:-1,:)
    tke_prod_sh_p(:,kms:kme-1,:) = tke_prod_sh_p(:,kme-1:kms:-1,:)
    w_sec_p(:,kms:kme-1,:)       = w_sec_p(:,kme-1:kms:-1,:)
    wqls_sec_p(:,kms:kme-1,:)    = wqls_sec_p(:,kme-1:kms:-1,:)
    brunt_p(:,kms:kme-1,:)       = brunt_p(:,kme-1:kms:-1,:)
    shoc_ql2_p(:,kms:kme-1,:)    = shoc_ql2_p(:,kme-1:kms:-1,:)
    dse_p(:,kms:kme-1,:)         = dse_p(:,kme-1:kms:-1,:)
    u_p(:,kms:kme-1,:)           = u_p(:,kme-1:kms:-1,:)
    v_p(:,kms:kme-1,:)           = v_p(:,kme-1:kms:-1,:)
    wthv_sec_p(:,kms:kme-1,:)    = wthv_sec_p(:,kme-1:kms:-1,:)
    tke_shmf_p(:,kms:kme-1,:)    = tke_shmf_p(:,kme-1:kms:-1,:)
    kzh_p(:,kms:kme-1,:)         = kzh_p(:,kme-1:kms:-1,:)
    kzm_p(:,kms:kme-1,:)         = kzm_p(:,kme-1:kms:-1,:)
    qc_p(:,kms:kme-1,:)          = qc_p(:,kme-1:kms:-1,:)
    cldfrabl_p(:,kms:kme-1,:)    = cldfrabl_p(:,kme-1:kms:-1,:)
    inv_exner_p(:,kms:kme-1,:)   = inv_exner_p(:,kme-1:kms:-1,:)
    qw_p(:,kms:kme-1,:)          = qw_p(:,kme-1:kms:-1,:)
    pres_p(:,kms:kme-1,:)        = pres_p(:,kme-1:kms:-1,:)
    height_mid_p(:,kms:kme-1,:)  = height_mid_p(:,kme-1:kms:-1,:)
    thetal_p(:,kms:kme-1,:)      = thetal_p(:,kme-1:kms:-1,:)
    presdiff_p(:,kms:kme-1,:)    = presdiff_p(:,kme-1:kms:-1,:)
    w_p(:,kms:kme-1,:)           = w_p(:,kme-1:kms:-1,:)
 end subroutine flip_vertical_coordinates

  
 subroutine convert_heat_fluxes(ims, ime, jms, jme, hfx_p, lh_p, rho_p, cp_in, xlv_in, hfx_converted_p, lh_converted_p)
    ! Convert the surface heat fluxes from W/m^2 to both K*m/s and kg/kg m/s
    ! hfx_p is the sensible heat flux
    ! lh_p is the latent heat flux
    integer, intent(in) :: ims, ime, jms, jme
    real(kind=RKIND), intent(in) :: cp_in, xlv_in
    real(kind=RKIND), dimension(:,:), intent(in) :: hfx_p, lh_p
    real(kind=RKIND), dimension(:,:,:), intent(in) :: rho_p
    real(kind=RKIND), dimension(:,:), intent(inout) :: hfx_converted_p, lh_converted_p

    integer :: i, j

    do j = jms,jme
    do i = ims,ime
        ! Converting hfx_p from W/m^2 to K*m/s
        hfx_converted_p(i,j) = hfx_p(i,j)/(cp_in*rho_p(i,1,j)) ! [K m/s] = [W/m^2]/([J kg-1 K-1]*[kg m^{-3}])
        ! hfx_converted_p(i,j) = hfx_converted_p(i,j)*inv_exner_p(i,1,j) ! Got from E3SM shoc_intr.F90 code
        ! Converting lh_p from W/m^2 to kg/kg m/s
        lh_converted_p(i,j) = lh_p(i,j)/(xlv_in*rho_p(i,1,j)) ! [kg/kg m/s] = [W/m^2]/([J kg-1]*[kg m^{-3}])
        ! xls instead of xlv - shoc_energy_total_fixer method in shoc
        ! lh_p/xlv = qfx_p
    enddo
    enddo
 end subroutine convert_heat_fluxes


 subroutine calculate_presref_mid(ims,ime,kms,kme,jms,jme,its,ite,kts,kte,jts,jte,diag,presref_mid)
    ! Calculate the reference pressure at the mid-levels
    ! For now it is calculated as the average of the pressure at the levels
    integer, intent(in) :: ims,ime,kms,kme,jms,jme,its,ite,kts,kte,jts,jte
    type(mpas_pool_type), intent(in) :: diag
    real(RKIND), dimension(:), intent(inout) :: presref_mid

    integer :: i, j, k
    real(kind=RKIND), dimension(:,:,:), allocatable :: pres_p
    real(kind=RKIND),dimension(:,:), pointer :: pressure_b,pressure_p

    call mpas_pool_get_array(diag,'pressure_base',pressure_b)
    call mpas_pool_get_array(diag,'pressure_p'   ,pressure_p)
    if(.not.allocated(pres_p) ) allocate(pres_p(ims:ime,kms:kme,jms:jme))

    do k=kts,kte
    do j=jts,jte
    do i=its,ite
        pres_p(i,k,j) = pressure_p(k,i) + pressure_b(k,i)
    enddo
    enddo
    enddo

    presref_mid = 0.0_RKIND

    do k=kts,kte
    do j=jts,jte
    do i=its,ite
        presref_mid(k) = presref_mid(k) + pres_p(i,k,j)
    enddo
    enddo
    enddo
    presref_mid = presref_mid/ite

    deallocate(pres_p)
 end subroutine calculate_presref_mid

 subroutine calculate_heights(ims, ime, kms, kme, jms, jme, z_p, zmid_p, height_interface_p, height_mid_p)
    ! From z_p and zmid_p which are the heights in reference to mean sea level,
    ! height_i_p and height_mid_p are calculated in reference to surface
    integer, intent(in) :: ims, ime, kms, kme, jms, jme
    real(RKIND), dimension(:,:,:), intent(in) :: z_p, zmid_p
    real(RKIND), dimension(:,:,:), intent(out) :: height_interface_p, height_mid_p

    integer :: i, j, k
    
    do i=ims,ime
    do j=jms,jme
    do k=kms,kme
        height_interface_p(i,k,j) = z_p(i,k,j)    - z_p(i,1,j)

        if (k .ne. kme) then
            height_mid_p(i,k,j)   = zmid_p(i,k,j) - z_p(i,1,j)
        end if
    enddo
    enddo
    enddo
 end subroutine calculate_heights

 subroutine calculate_kpbl(ims, ime, jms, jme, kms, kme, hpbl_p, height_mid_p, kpbl_p)
    ! Calculate the PBL height kpbl_p from hpbl_p and height_mid_p
    ! We get the last vertical level that is below the PBL height
    ! hpbl_p is the PBL height in reference to surface
    ! height_mid_p is the height thermodynamic grid in reference to surface
    integer, intent(in) :: ims, ime, jms, jme, kms, kme
    real(RKIND), dimension(:,:), intent(in) :: hpbl_p
    real(RKIND), dimension(:,:,:), intent(in) :: height_mid_p
    integer, dimension(:,:), intent(inout) :: kpbl_p

    integer :: i, j, k
    
    do i=ims,ime
    do j=jms,jme
    do k=kms,kme-1
        if (height_mid_p(i,k,j) > hpbl_p(i,j)) then
            kpbl_p(i,j) = k-1
            exit
        endif
    enddo
    enddo
    enddo
 end subroutine calculate_kpbl

 subroutine bl_shocmf_init( &
    ims,ime,kms,kme,jms,jme,its,ite,kts,kte,jts,jte, gravity_in, Rd_in, Rv_in, cp_in, ep_1_in, xlv_in, &
    xlf_in, karman_in, diag &
 )
    ! Initialize the SHOC-MF module
    ! This subroutine is called from the main driver
    integer, intent(in) :: ims,ime,kms,kme,jms,jme,its,ite,kts,kte,jts,jte
    type(mpas_pool_type), intent(in) :: diag
    real(RKIND), intent(in) :: gravity_in, Rd_in, Rv_in, cp_in, ep_1_in, xlv_in, xlf_in, &
                               karman_in

    integer:: ntop_shoc, nbot_shoc
    real(RKIND), dimension(:), allocatable :: presref_mid
    allocate(presref_mid(kms:kme))
    
    call calculate_presref_mid(ims,ime,kms,kme,jms,jme,its,ite,kts,kte,jts,jte,diag,presref_mid)

    ntop_shoc = kte
    nbot_shoc = 1

    call shoc_init(kte, gravity_in, Rd_in, Rv_in, cp_in, ep_1_in, xlv_in, xlf_in, karman_in, &
                   presref_mid, nbot_shoc, ntop_shoc, length_fac)

    deallocate(presref_mid)
 end subroutine bl_shocmf_init

 subroutine bl_shocmf_driver( &
    itimestep, dt_pbl, dx_p, lh_p, hfx_p, ht_p, theta_p, &
    z_p, zmid_p, pres_p, pres2_p, uw_sfc_p, vw_sfc_p, w_p, pi_p, dse_p, tke_shmf_p, qw_p, &
    u_p, v_p, wthv_sec_p, kzh_p, kzm_p, qc_p, cldfrabl_p, hpbl_p, pblh_wthl_p, pblh_tke_p, &
    rino_p, shoc_mix_p, isotropy_p, tke_dissip_p, tke_prod_b_p, tke_prod_sh_p, w_sec_p, &
    thl_sec_p, qw_sec_p, qwthl_sec_p, wthl_sec_p, wqw_sec_p, wtke_sec_p, uw_sec_p, &
    vw_sec_p, w3_p, wqls_sec_p, brunt_p, shoc_ql2_p, mf_dry_a_p, mf_moist_a_p, mf_dry_w_p, &
    mf_moist_w_p, mf_dry_qt_p, mf_moist_qt_p, mf_dry_thl_p, mf_moist_thl_p, mf_dry_u_p, &
    mf_moist_u_p, mf_dry_v_p, mf_moist_v_p, mf_moist_qc_p, mf_thlflx_p, mf_qtflx_p, &
    mf_ae_p, mf_aw_p, mf_awthv_p, mf_awthl_p, mf_awqt_p, mf_awql_p, mf_awqi_p, mf_awu_p, &
    mf_awv_p, cfl_mf_p, mf_dry_freq_p, mf_moist_freq_p, plumeheight_p, plume_dry_height_p, &
    ustar_p, wspd_p, ztop_b_p, wthl_sec_ed_p, wthl_sec_mf_p, wqw_sec_ed_p, wqw_sec_mf_p, &
    kpbl_p, rho_p, t_p, qv_p, rublten_p, rvblten_p, rthblten_p, &
    rqvblten_p, rqcblten_p, rqiblten_p, rtkeblten_p, &
    ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, dminfo &
 )
    ! Driver for the SHOC-MF module
    ! This subroutine is called from the main driver
    ! It calls the SHOC-MF module to calculate the turbulence
    ! and the EDMF module to calculate the subgrid-scale clouds
    integer, intent(in) :: itimestep, &
                           ims, ime, jms, jme, kms, kme, its, ite, &
                           jts, jte, kts, kte
    type(dm_info), intent(in) :: dminfo
    real(kind=RKIND), intent(in) :: dt_pbl
    real(kind=RKIND),dimension(:,:), intent(in) :: dx_p, lh_p, hfx_p, ht_p

    real(kind=RKIND),dimension(:,:), intent(inout) :: uw_sfc_p, vw_sfc_p, hpbl_p, pblh_wthl_p, pblh_tke_p, &
                                                      mf_dry_freq_p, mf_moist_freq_p, plumeheight_p, & 
                                                      plume_dry_height_p, ztop_b_p, ustar_p, wspd_p
    integer, dimension(:,:), intent(inout) :: kpbl_p
    real(kind=RKIND),dimension(:,:,:), intent(inout) :: rino_p, shoc_mix_p, isotropy_p, &
        tke_dissip_p, tke_prod_b_p, tke_prod_sh_p, w_sec_p, thl_sec_p, qw_sec_p, &
        qwthl_sec_p, wthl_sec_p, wqw_sec_p, wtke_sec_p, uw_sec_p, vw_sec_p, w3_p, &
        wqls_sec_p, brunt_p, shoc_ql2_p, mf_dry_a_p, mf_moist_a_p, mf_dry_w_p, &
        mf_moist_w_p, mf_dry_qt_p, mf_moist_qt_p, mf_dry_thl_p, mf_moist_thl_p, &
        mf_dry_u_p, mf_moist_u_p, mf_dry_v_p, mf_moist_v_p, mf_moist_qc_p, mf_thlflx_p, &
        mf_qtflx_p, mf_ae_p, mf_aw_p, mf_awthv_p, mf_awthl_p, mf_awqt_p, mf_awql_p, &
        mf_awqi_p, mf_awu_p, mf_awv_p, cfl_mf_p, wthl_sec_ed_p, wthl_sec_mf_p, wqw_sec_ed_p, &
        wqw_sec_mf_p, dse_p, qw_p, u_p, v_p, wthv_sec_p, tke_shmf_p, kzh_p, kzm_p, qc_p, &
        cldfrabl_p, theta_p, t_p, z_p, zmid_p, pres_p, pres2_p, w_p, pi_p, qv_p, &
        rho_p, rublten_p, rvblten_p, rthblten_p, rqvblten_p, rtkeblten_p, &
        rqcblten_p, rqiblten_p

    real(kind=RKIND), dimension(:,:,:), allocatable :: presdiff_p, thetal_p, qtracers, height_interface_p, height_mid_p
    real(kind=RKIND), dimension(:,:,:), allocatable :: um_p, vm_p, thetam_p, qvm_p, qcm_p, tkem_p, w_zt
    real(kind=RKIND), dimension(:,:,:), allocatable :: inv_exner_p, rho_interface_p, thv_p
    real(kind=RKIND), dimension(:,:)  , allocatable :: wtracer_sfc, dy_p, hfx_converted_p, lh_converted_p, phis_p

    integer :: num_qtracers, nadv, i, j, k
    real(kind=RKIND) :: shoc_timestep

    allocate(                                  &
        um_p(ims:ime, kms:kme, jms:jme)      , &
        vm_p(ims:ime, kms:kme, jms:jme)      , &
        thetam_p(ims:ime, kms:kme, jms:jme), &
        qvm_p(ims:ime, kms:kme, jms:jme)     , &
        qcm_p(ims:ime, kms:kme, jms:jme)     , &
        tkem_p(ims:ime, kms:kme, jms:jme)      &
    )
    ! Deep copying the input variables
    um_p = u_p
    vm_p = v_p
    qvm_p = qv_p
    qcm_p = qc_p
    tkem_p = tke_shmf_p
    thetam_p = theta_p

    allocate(inv_exner_p(ims:ime, kms:kme, jms:jme))
    do i=ims,ime
    do k=kms,kme-1
    do j=jms,jme
        inv_exner_p(i,k,j) = 1.0_RKIND/pi_p(i,k,j) ! inv_exner_p = theta_p/t_p
        tke_shmf_p(i,k,j) = max(mintke, tke_shmf_p(i,k,j))
    enddo
    enddo
    enddo

    allocate(dy_p(ims:ime, jms:jme))
    dy_p = dx_p ! Assuming dx_p is the same for both x and y directions

    allocate(thv_p(ims:ime, kms:kme, jms:jme))
    call calculate_thetav(ims, ime, jms, jme, kms, kme, theta_p, qv_p, qc_p, thv_p)

    allocate(thetal_p(ims:ime, kms:kme, jms:jme))
    call calculate_thetal(ims, ime, jms, jme, kms, kme, cp, xlv, theta_p, qc_p, inv_exner_p, thetal_p)

    allocate(presdiff_p(ims:ime, kms:kme, jms:jme))
    call calculate_presdiff(ims, ime, kms, kme, jms, jme, pres2_p, presdiff_p)

    call calculate_qw(ims, ime, jms, jme, kms, kme, qc_p, qv_p, qw_p)

    allocate(height_interface_p(ims:ime, kms:kme, jms:jme), height_mid_p(ims:ime, kms:kme, jms:jme))
    call calculate_heights(ims, ime, kms, kme, jms, jme, z_p, zmid_p, height_interface_p, height_mid_p)

    allocate(phis_p(ims:ime, jms:jme))
    do i=ims,ime
    do j=jms,jme
        phis_p(i,j) = gravity*ht_p(i,j)
    enddo
    enddo
    call calculate_dry_static_energy(ims, ime, kms, kme, jms, jme, t_p, height_mid_p, gravity, phis_p, cp, dse_p)

    call calculate_surface_momentum_fluxes(ims, ime, jms, jme, ustar_p, u_p(:,1,:), v_p(:,1,:), wspd_p, rho_p(:,1,:), uw_sfc_p, vw_sfc_p)

    allocate(w_zt(ims:ime, kms:kme, jms:jme))
    call linear_interp(height_interface_p,height_mid_p,w_p,w_zt,kme,kme-1,ime,-200.0_RKIND)

    allocate(hfx_converted_p(ims:ime, jms:jme), lh_converted_p(ims:ime, jms:jme), rho_interface_p(ims:ime, kms:kme, jms:jme))
    ! Get density on the interface grid
    call linear_interp(height_mid_p,height_interface_p,rho_p,rho_interface_p,kme-1,kme,ime,0.0_RKIND)
    call convert_heat_fluxes(ims, ime, jms, jme, hfx_p, lh_p, rho_interface_p, cp, xlv, hfx_converted_p, lh_converted_p)

    nadv = ceiling(dt_pbl/max_shoc_timestep)
    if (nadv .gt. 1) then
        shoc_timestep = dt_pbl/nadv
    else
        shoc_timestep = dt_pbl
    endif

    ! If mpi rank is 0, print some debugging information
    if (dminfo % my_proc_id == 5) then
        print *, "1. Before flipping vertical coordinates:"
        print *, ""
        print *, "dx_p(1,1)", dx_p(1,1)
        print *, ""
        print *, "height_interface_p(1, :, 1): ", height_interface_p(1,:,1)
        print *, ""
        print *, "height_mid_p(1, :, 1): ", height_mid_p(1,:,1)
        print *, ""
        print *, "w_p(1, :, 1): ", w_p(1,:,1)
        print *, ""
        print *, "w_zt(1, :, 1): ", w_zt(1,:,1)
        print *, ""
        print *, "pres_p(1, :, 1): ", pres_p(1,:,1)
        print *, ""
        print *, "pres2_p(1, :, 1): ", pres2_p(1,:,1)
        print *, ""
        print *, "presdiff_p(1, :, 1): ", presdiff_p(1,:,1)
        print *, ""
        print *, "rho_p(1, :, 1): ", rho_p(1,:,1)
        print *, ""
        print *, "rho_interface_p(1, :, 1): ", rho_interface_p(1,:,1)
        print *, ""
        print *, "theta_p(1, :, 1): ", theta_p(1,:,1)
        print *, ""
        print *, "thetal_p(1, :, 1): ", thetal_p(1,:,1)
        print *, ""
        print *, "thetav(1, :, 1): ", thv_p(1,:,1)
        print *, ""
        print *, "qv_p(1, :, 1): ", qv_p(1,:,1)
        print *, ""
        print *, "qc_p(1, :, 1): ", qc_p(1,:,1)
        print *, ""
        print *, "qw_p(1, :, 1): ", qw_p(1,:,1)
        print *, ""
        print *, "hfx_p(1, 1): ", hfx_p(1,1)
        print *, ""
        print *, "hfx_converted_p(1, 1): ", hfx_converted_p(1,1)
        print *, ""
        print *, "lh_p(1, 1): ", lh_p(1,1)
        print *, ""
        print *, "lh_converted_p(1, 1): ", lh_converted_p(1,1)
        print *, ""
        print *, "uw_sfc_p(1, 1): ", uw_sfc_p(1,1)
        print *, ""
        print *, "vw_sfc_p(1, 1): ", vw_sfc_p(1,1)
    endif

    ! If resulting host model time step and SHOC time step do not divide evenly
    ! into each other, have model throw a fit.
    if (mod(dt_pbl,shoc_timestep) .ne. 0) then
        stop "SHOC time step and HOST time step NOT compatible"
    endif

    num_qtracers = 0
    allocate(wtracer_sfc(ims:ime, max(0, num_qtracers)))
    wtracer_sfc = 0.0_RKIND
    allocate(qtracers(ims:ime, kms:kme, max(0,num_qtracers)))
    qtracers = 0.0_RKIND

    call flip_vertical_coordinates(&
        kms, kme, rino_p, shoc_mix_p, isotropy_p, tke_dissip_p, tke_prod_b_p, &
        tke_prod_sh_p, w_sec_p, thl_sec_p, qw_sec_p, qwthl_sec_p, &
        wthl_sec_p, wqw_sec_p, wtke_sec_p, uw_sec_p, vw_sec_p, w3_p, &
        wqls_sec_p, brunt_p, shoc_ql2_p, mf_dry_a_p, mf_moist_a_p, &
        mf_dry_w_p, mf_moist_w_p, mf_dry_qt_p, mf_moist_qt_p, mf_dry_thl_p, &
        mf_moist_thl_p, mf_dry_u_p, mf_moist_u_p, mf_dry_v_p, mf_moist_v_p, &
        mf_moist_qc_p, mf_thlflx_p, mf_qtflx_p, mf_ae_p, mf_aw_p, &
        mf_awthv_p, mf_awthl_p, mf_awqt_p, mf_awql_p, mf_awqi_p, mf_awu_p, &
        mf_awv_p, cfl_mf_p, wthl_sec_ed_p, wthl_sec_mf_p, wqw_sec_ed_p, &
        wqw_sec_mf_p, dse_p, qw_p, u_p, v_p, wthv_sec_p, tke_shmf_p, kzh_p, &
        kzm_p, qc_p, cldfrabl_p, thv_p, height_interface_p, &
        height_mid_p, pres_p, pres2_p, w_zt, inv_exner_p, &
        thetal_p, presdiff_p &
    )

    if (dminfo % my_proc_id == 5) then
        print *, ""
        print *, "2. After flipping vertical coordinates:"
        print *, "height_interface_p(1, :, 1): ", height_interface_p(1,:,1)
        print *, ""
        print *, "height_mid_p(1, :, 1): ", height_mid_p(1,:,1)
        print *, ""
        print *, "w_p(1, :, 1): ", w_p(1,:,1)
        print *, ""
        print *, "w_zt(1, :, 1): ", w_zt(1,:,1)
        print *, ""
        print *, "pres_p(1, :, 1): ", pres_p(1,:,1)
        print *, ""
        print *, "pres2_p(1, :, 1): ", pres2_p(1,:,1)
        print *, ""
        print *, "presdiff_p(1, :, 1): ", presdiff_p(1,:,1)
        print *, ""
        print *, "rho_p(1, :, 1): ", rho_p(1,:,1)
        print *, ""
        print *, "rho_interface_p(1, :, 1): ", rho_interface_p(1,:,1)
        print *, ""
        print *, "theta_p(1, :, 1): ", theta_p(1,:,1)
        print *, ""
        print *, "thetal_p(1, :, 1): ", thetal_p(1,:,1)
        print *, ""
        print *, "thetav(1, :, 1): ", thv_p(1,:,1)
        print *, ""
        print *, "qv_p(1, :, 1): ", qv_p(1,:,1)
        print *, ""
        print *, "qc_p(1, :, 1): ", qc_p(1,:,1)
        print *, ""
        print *, "qw_p(1, :, 1): ", qw_p(1,:,1)
    endif
    ! As last position of the array is a extra position, to not letting a
    ! zero value, we are going to extrapolate it considering that
    ! the kme-1 value is the average among kme-2 and kme
    ! thetal_p(:,kme-1,:) = (thetal_p(:,kme-2,:) + thetal_p(:,kme,:))/2
    thetal_p(:,kme,:) = 2*thetal_p(:,kme-1,:) - thetal_p(:,kme-2,:)

    call shoc_main( &
        itimestep, ime, kme-1, kme, shoc_timestep, nadv, &                           ! Input
        dx_p(:ime,1), dy_p(:ime,1), thv_p(:ime,:kme-1,1),&                           ! Input
        height_mid_p(:ime,:kme-1,1), height_interface_p(:ime,:kme,1), &              ! Input
        pres_p(:ime,:kme-1,1), pres2_p(:ime,:kme,1), &                               ! Input
        presdiff_p(:ime,:kme-1,1), hfx_converted_p(:ime,1), &                        ! Input
        lh_converted_p(:ime,1), uw_sfc_p(:ime,1), vw_sfc_p(:ime,1), &                ! Input
        wtracer_sfc(:ime,:num_qtracers), num_qtracers, w_zt(:ime,:kme-1,1), &        ! Input
        inv_exner_p(:ime,:kme-1,1),phis_p(:ime,1), &                                 ! Input
        dse_p(:ime,:kme-1,1), tke_shmf_p(:ime,:kme-1,1), &                           ! Input/Output
        thetal_p(:ime,:kme-1,1), qw_p(:ime,:kme-1,1), &                              ! Input/Output
        u_p(:ime,:kme-1,1), v_p(:ime,:kme-1,1), qtracers(:,:kme-1,:num_qtracers), &  ! Input/Output
        wthv_sec_p(:ime,:kme-1,1), kzh_p(:ime,:kme-1,1),kzm_p(:ime,:kme-1,1), &      ! Input/Output
        qc_p(:ime,:kme-1,1),cldfrabl_p(:ime,:kme-1,1), &                             ! Input/Output
        hpbl_p(:ime,1),pblh_wthl_p(:ime,1),pblh_tke_p(:ime,1), &                     ! Output
        rino_p(:ime,:kme-1,1), &                                                     ! Output
        shoc_mix_p(:ime,:kme-1,1), isotropy_p(:ime,:kme-1,1), &                      ! Output (diagnostic)
        tke_dissip_p(:ime,:kme-1,1), tke_prod_b_p(:ime,:kme-1,1), &                  ! Output (diagnostic)
        tke_prod_sh_p(:ime,:kme-1,1), &                                              ! Output (diagnostic)
        w_sec_p(:ime,:kme-1,1), thl_sec_p(:ime,:kme,1), qw_sec_p(:ime,:kme,1), &     ! Output (diagnostic)
        qwthl_sec_p(:ime,:kme,1), &                                                  ! Output (diagnostic)
        wthl_sec_p(:ime,:kme,1), wqw_sec_p(:ime,:kme,1), wtke_sec_p(:ime,:kme,1), &  ! Output (diagnostic)
        uw_sec_p(:ime,:kme,1), vw_sec_p(:ime,:kme,1), w3_p(:ime,:kme,1), &           ! Output (diagnostic)
        wqls_sec_p(:ime,:kme-1,1),brunt_p(:ime,:kme-1,1),shoc_ql2_p(:ime,:kme-1,1), &! Output (diagnostic)
        mf_dry_a_p(:ime,:kme,1), mf_moist_a_p(:ime,:kme,1), &                        ! Output (EDMF diagnostic)
        mf_dry_w_p(:ime,:kme,1), mf_moist_w_p(:ime,:kme,1), &                        ! Output (EDMF diagnostic)
        mf_dry_qt_p(:ime,:kme,1), mf_moist_qt_p(:ime,:kme,1), &                      ! Output (EDMF diagnostic)
        mf_dry_thl_p(:ime,:kme,1), mf_moist_thl_p(:ime,:kme,1), &                    ! Output (EDMF diagnostic)
        mf_dry_u_p(:ime,:kme,1),   mf_moist_u_p(:ime,:kme,1), &                      ! Output (EDMF diagnostic)
        mf_dry_v_p(:ime,:kme,1),   mf_moist_v_p(:ime,:kme,1), &                      ! Output (EDMF diagnostic)
        mf_moist_qc_p(:ime,:kme,1), &                                                ! Output (EDMF diagnostic)
        mf_thlflx_p(:ime,:kme,1), mf_qtflx_p(:ime,:kme,1), &                         ! Output (EDMF diagnostic)
        mf_ae_p(:ime,:kme,1), mf_aw_p(:ime,:kme,1), &                                ! Output (EDMF diagnostic)
        mf_awthv_p(:ime,:kme,1), &                                                   ! Output (EDMF diagnostic)
        mf_awthl_p(:ime,:kme,1), mf_awqt_p(:ime,:kme,1), &                           ! Output (EDMF diagnostic)
        mf_awql_p(:ime,:kme,1), mf_awqi_p(:ime,:kme,1), &                            ! Output (EDMF diagnostic)
        mf_awu_p(:ime,:kme,1), mf_awv_p(:ime,:kme,1), cfl_mf_p(:ime,:kme,1), &       ! Output (EDMF diagnostic)
        mf_dry_freq_p(:,1), mf_moist_freq_p(:,1), plumeheight_p(:,1), &              ! Output (EDMF diagnostic)
        plume_dry_height_p(:,1), ztop_b_p(:,1), &                                    ! Output (EDMF diagnostic)
        wthl_sec_ed_p(:ime,:kme,1), wthl_sec_mf_p(:ime,:kme,1), &                    ! Output (EDMF diagnostic)
        wqw_sec_ed_p(:ime,:kme,1), wqw_sec_mf_p(:ime,:kme,1) &                       ! Output (EDMF diagnostic)  
    )
    call flip_vertical_coordinates(&
        kms, kme, rino_p, shoc_mix_p, isotropy_p, tke_dissip_p, tke_prod_b_p, &
        tke_prod_sh_p, w_sec_p, thl_sec_p, qw_sec_p, qwthl_sec_p, &
        wthl_sec_p, wqw_sec_p, wtke_sec_p, uw_sec_p, vw_sec_p, w3_p, &
        wqls_sec_p, brunt_p, shoc_ql2_p, mf_dry_a_p, mf_moist_a_p, &
        mf_dry_w_p, mf_moist_w_p, mf_dry_qt_p, mf_moist_qt_p, mf_dry_thl_p, &
        mf_moist_thl_p, mf_dry_u_p, mf_moist_u_p, mf_dry_v_p, mf_moist_v_p, &
        mf_moist_qc_p, mf_thlflx_p, mf_qtflx_p, mf_ae_p, mf_aw_p, &
        mf_awthv_p, mf_awthl_p, mf_awqt_p, mf_awql_p, mf_awqi_p, mf_awu_p, &
        mf_awv_p, cfl_mf_p, wthl_sec_ed_p, wthl_sec_mf_p, wqw_sec_ed_p, &
        wqw_sec_mf_p, dse_p, qw_p, u_p, v_p, wthv_sec_p, tke_shmf_p, kzh_p, &
        kzm_p, qc_p, cldfrabl_p, thv_p, height_interface_p, &
        height_mid_p, pres_p, pres2_p, w_zt, inv_exner_p, &
        thetal_p, presdiff_p &
    )
    ! As last position of the array is a extra position,
    ! we are going to extrapolate it considering that
    ! the kme-1 value is the average among kme-2 and kme
    ! thetal_p(:,kme-1,:) = (thetal_p(:,kme-2,:) + thetal_p(:,kme,:))/2
    thetal_p(:,kme,:) = 2*thetal_p(:,kme-1,:) - thetal_p(:,kme-2,:)

    call calculate_theta(ims, ime, jms, jme, kms, kme, cp, xlv, thetal_p, qc_p, inv_exner_p, theta_p)
    call calculate_qv(ims, ime, jms, jme, kms, kme, qw_p, qc_p, qv_p)
    call calculate_kpbl(ims, ime, jms, jme, kms, kme, hpbl_p, height_mid_p, kpbl_p)

    if (dminfo % my_proc_id == 5) then
        print *, ""
        print *, "3. After reflipping vertical coordinates:"
        print *, "height_interface_p(1, :, 1): ", height_interface_p(1,:,1)
        print *, ""
        print *, "height_mid_p(1, :, 1): ", height_mid_p(1,:,1)
        print *, ""
        print *, "w_p(1, :, 1): ", w_p(1,:,1)
        print *, ""
        print *, "w_zt(1, :, 1): ", w_zt(1,:,1)
        print *, ""
        print *, "pres_p(1, :, 1): ", pres_p(1,:,1)
        print *, ""
        print *, "pres2_p(1, :, 1): ", pres2_p(1,:,1)
        print *, ""
        print *, "presdiff_p(1, :, 1): ", presdiff_p(1,:,1)
        print *, ""
        print *, "rho_p(1, :, 1): ", rho_p(1,:,1)
        print *, ""
        print *, "rho_interface_p(1, :, 1): ", rho_interface_p(1,:,1)
        print *, ""
        print *, "theta_p(1, :, 1): ", theta_p(1,:,1)
        print *, ""
        print *, "thetal_p(1, :, 1): ", thetal_p(1,:,1)
        print *, ""
        print *, "thetav(1, :, 1): ", thv_p(1,:,1)
        print *, ""
        print *, "qv_p(1, :, 1): ", qv_p(1,:,1)
        print *, ""
        print *, "qc_p(1, :, 1): ", qc_p(1,:,1)
        print *, ""
        print *, "qw_p(1, :, 1): ", qw_p(1,:,1)
        print *, ""
        print *, "hpbl_p(1, 1): ", hpbl_p(1,1)
        print *, ""
        print *, "kpbl_p(1, 1): ", kpbl_p(1,1)
    endif

    ! Calculating tendencies
    do i=ims,ime
    do k=kms,kme-1
    do j=jms,jme
        rublten_p(i,k,j)   = (u_p(i,k,j) - um_p(i,k,j))/dt_pbl
        rvblten_p(i,k,j)   = (v_p(i,k,j) - vm_p(i,k,j))/dt_pbl
        rthblten_p(i,k,j)  = (theta_p(i,k,j) - thetam_p(i,k,j))/dt_pbl
        rqvblten_p(i,k,j)  = (qv_p(i,k,j) - qvm_p(i,k,j))/dt_pbl
        rqcblten_p(i,k,j)  = (qc_p(i,k,j) - qcm_p(i,k,j))/dt_pbl
        rtkeblten_p(i,k,j) = (tke_shmf_p(i,k,j) - tkem_p(i,k,j))/dt_pbl
    enddo
    enddo
    enddo
    
    rqiblten_p = 0.0_RKIND  ! TODO (GTF): qi_p and qs_p are not calculated by shocmf. Only qc_p and qw_p

    deallocate(presdiff_p, thetal_p, qtracers, wtracer_sfc, w_zt)
    deallocate(um_p, vm_p, thetam_p, qvm_p, qcm_p, tkem_p, phis_p)
    deallocate(hfx_converted_p, lh_converted_p)
    deallocate(dy_p)
    deallocate(height_interface_p, height_mid_p)
    deallocate(inv_exner_p)
  
 end subroutine bl_shocmf_driver

end module module_bl_shocmf
