module module_bl_shocmf
 use mpas_atmphys_constants
 use shoc, only: shoc_init, shoc_main, length_fac, mintke
 use mpas_pool_routines


 implicit none
 private
 public:: bl_shocmf_init, bl_shocmf_driver

 real(kind=RKIND) :: max_shoc_timestep = 300.0  ! Default SHOC timestep set in namelist

 contains

 subroutine calculate_surface_momentum_fluxes(ime, kme, cda_p, u_p, v_p, uw_sfc_p, vw_sfc_p)
    ! Calculate the surface momentum flux
    ! Reference: Stull's Book - An Introduction to Boundary Layer Meteorology, chapter 7.4.1.
    integer, intent(in) :: ime, kme
    real(kind=RKIND), dimension(:,:), intent(in) :: cda_p, u_p, v_p
    real(kind=RKIND), dimension(:,:), intent(out) :: uw_sfc_p, vw_sfc_p
    
    real(kind=RKIND), dimension(:,:), allocatable :: mag_velocity
    allocate(mag_velocity(ime, kme))

    mag_velocity = sqrt(u_p**2 + v_p**2)

    uw_sfc_p = -cda_p*mag_velocity*u_p
    vw_sfc_p = -cda_p*mag_velocity*v_p
 end subroutine calculate_surface_momentum_fluxes


 subroutine calculate_dry_static_energy(ims, ime, kms, kme, jts, jte, t, z, gravity_in, ht_p, cp_in, dse)
    ! Calculate the dry static energy considering the terrain height.
    ! This is the calculation waited by the SHOCMF module. 
    ! Reference: https://glossary.ametsoc.org/wiki/Dry_static_energy
    integer, intent(in) :: ims, ime, kms, kme, jts, jte
    real(kind=RKIND), dimension(:,:,:), intent(in) :: t, z
    real(kind=RKIND), dimension(:,:), intent(in) :: ht_p
    real(kind=RKIND), intent(in) :: gravity_in, cp_in
    real(kind=RKIND), dimension(:,:,:), intent(out) :: dse

    integer :: i, j, k
    
    do i=ims,ime
    do k=kms,kme
    do j=jts,jte
        dse(i,k,j) = cp_in*t(i,k,j) + gravity_in*z(i,k,j) + ht_p(i,j)
    enddo
    enddo
    enddo
 end subroutine calculate_dry_static_energy


 subroutine calculate_presdiff(pres_levels, ims, ime, kms, kme, jms, jme, presdiff_levels)
    ! Calculate the pressure difference between the interfaces
    ! pres_levels is the pressure at the levels
    ! presdiff_levels is the pressure difference between the interfaces
    real(kind=RKIND), dimension(:,:,:), intent(in) :: pres_levels
    integer, intent(in) :: ims,ime,kms,kme,jms,jme
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


 subroutine calculate_thetal(ime, jme, kme, theta, qc, inv_exner, thetal)
    ! Calculate liquid water potential temperature
    ! Reference: shoc_intr.F90 in E3SM code and
    ! https://glossary.ametsoc.org/wiki/Liquid_water_potential_temperature
    integer, intent(in) :: ime, jme, kme
    real(kind=RKIND), dimension(:,:,:), intent(in) :: theta, qc, inv_exner
    real(kind=RKIND), dimension(:,:,:), intent(inout) :: thetal
    
    thetal = theta - inv_exner*(xlv/cp)*qc
 end subroutine calculate_thetal


 subroutine calculate_theta(thetal, qc, inv_exner, theta)
    ! Calculate potential temperature from liquid water potential temperature.
    ! Inverse function of calculate_thetal.
    real(kind=RKIND), dimension(:,:,:), intent(in) :: thetal, qc, inv_exner
    real(kind=RKIND), dimension(:,:,:), intent(out) :: theta

    theta = thetal + inv_exner*(xlv/cp)*qc
 end subroutine calculate_theta


 subroutine calculate_thetam_to_shocmf(theta, qv, qc, R_v_in, R_d_in, theta_m)
    ! Update the moist potential temperature as SHOCMF expects it
    ! This is called before the SHOCMF call
    ! Reference: shoc_intr.F90 in E3SM code
    real(kind=RKIND), dimension(:,:,:), intent(in) :: theta, qv, qc
    real(kind=RKIND), intent(in) :: R_v_in, R_d_in
    real(kind=RKIND), dimension(:,:,:), intent(out) :: theta_m

    theta_m = theta*(1.0_RKIND + (R_v_in/R_d_in)*qv - qc)
 end subroutine calculate_thetam_to_shocmf


 subroutine calculate_thetam_to_mpas(theta, qv, R_v_in, R_d_in, theta_m)
   ! Calculate moist potential temperature from potential temperature.
   ! Equation got from Skamarock et al (2012)
   ! DOI: https://doi.org/10.1175/MWR-D-11-00215.1
   real(kind=RKIND), dimension(:,:,:), intent(in) :: theta, qv
   real(kind=RKIND), intent(in) :: R_v_in, R_d_in
   real(kind=RKIND), dimension(:,:,:), intent(out) :: theta_m

   theta_m = theta*(1.0_RKIND + (R_v_in/R_d_in)*qv)
 end subroutine calculate_thetam_to_mpas


 subroutine calculate_qv(qw, qc, qs, qi, qg, qr, qv)
    ! Calculate the water vapor mixing ratio
    real(kind=RKIND), dimension(:,:,:), intent(in) :: qw, qc, qs, qi, qg, qr
    real(kind=RKIND), dimension(:,:,:), intent(out) :: qv

    qv = qw - qc !- qs - qi - qg - qr
 end subroutine calculate_qv


 subroutine calculate_qw(qc, qv, qs, qi, qg, qr, qw)
    ! Calculate total water mixing ratio
    real(kind=RKIND), dimension(:,:,:), intent(in) :: qc, qv, qs, qi, qg, qr
    real(kind=RKIND), dimension(:,:,:), intent(out) :: qw

    qw = qc + qv !+ qs + qi + qg + qr
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
    kzm_p, qc_p, cldfrabl_p, theta_m_p, height_interface_p, &
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
                                                        theta_m_p, height_interface_p, height_mid_p, &
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
    ! of them are only kme-1 levels. So, in these cases we
    ! need to flip only the vertical coordinates of the kme-1 levels
    rino_p(:,kms:kme-1,:)        = rino_p(:,kme-1:kms:-1,:)
    theta_m_p(:,kms:kme-1,:)     = theta_m_p(:,kme-1:kms:-1,:)
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

 
 subroutine define_horizontal_discretization(ime, jme, area_p, dx_p, dy_p)
    ! Define the horizontal discretization
    ! dx_p is the horizontal grid spacing in the x direction
    ! dy_p is the horizontal grid spacing in the y direction
    integer, intent(in) :: ime, jme
    real(kind=RKIND), dimension(:,:), intent(in) :: area_p
    real(kind=RKIND), dimension(:,:), intent(inout) :: dx_p, dy_p

    ! dx_p and dy_p are calculated as sqrt(area_p) in the main driver
    dx_p = sqrt(area_p)
    dy_p = sqrt(area_p)
 end subroutine define_horizontal_discretization

 
 subroutine convert_heat_fluxes(ims, ime, jms, jme, hfx_p, lh_p, rho_p, cp_in, xlv_in, hfx_converted_p, lh_converted_p)
    ! Convert the surface heat fluxes from W/m^2 to both K*m/s and kg/kg m/s
    ! hfx_p is the sensible heat flux
    ! lh_p is the latent heat flux
    integer, intent(in) :: ims, ime, jms, jme
    real(kind=RKIND), intent(in) :: cp_in, xlv_in
    real(kind=RKIND), dimension(:,:), intent(in) :: hfx_p, lh_p
    real(kind=RKIND), dimension(:,:,:), intent(in) :: rho_p
    real(kind=RKIND), dimension(:,:), allocatable, intent(out) :: hfx_converted_p, lh_converted_p

    integer :: i, j

    allocate(hfx_converted_p(ime, jme), lh_converted_p(ime, jme))
    do j = jms,jme
    do i = ims,ime
        ! Converting hfx_p from W/m^2 to K*m/s
        hfx_converted_p(i,j) = hfx_p(i,j)/(cp_in*rho_p(i,1,j)) ! [K m/s] = [W/m^2]/([J kg-1 K-1]*[kg m^{-3}])
      
        ! Converting lh_p from W/m^2 to kg/kg m/s
        lh_converted_p(i,j) = lh_p(i,j)/(xlv_in*rho_p(i,1,j)) ! [kg/kg m/s] = [W/m^2]/([J kg-1]*[kg m^{-3}])
    enddo
    enddo
 end subroutine convert_heat_fluxes


 subroutine calculate_presref_mid(ims,ime,kms,kme,jms,jme,its,ite,kts,kte,jts,jte,diag,presref_mid)
    ! Calculate the reference pressure at the mid-levels
    ! For now it is calculated as the average of the pressure at the levels
    integer, intent(in) :: ims,ime,kms,kme,jms,jme,its,ite,kts,kte,jts,jte
    type(mpas_pool_type), intent(in) :: diag
    real(RKIND), dimension(:), allocatable, intent(inout) :: presref_mid

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

    allocate(presref_mid(kme))
    presref_mid = 0.0_RKIND

    do k=kts,kte
    do j=jts,jte
    do i=its,ite
        presref_mid(k) = presref_mid(k) + pres_p(i,k,j)
    enddo
    enddo
    enddo
    presref_mid = presref_mid/ime

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
        height_interface_p(i,k,j) = z_p(i,k,j) - z_p(i,1,j)
        height_mid_p(i,k,j) = zmid_p(i,k,j) - z_p(i,1,j)
    enddo
    enddo
    enddo
 end subroutine calculate_heights

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

    call calculate_presref_mid(ims,ime,kms,kme,jms,jme,its,ite,kts,kte,jts,jte,diag,presref_mid)

    ntop_shoc = 1    ! if >1, must be <= nbot_molec
    nbot_shoc = kte ! currently always pver 

    call mpas_log_write('Calling shoc_init from SHOCMF module')
    call shoc_init(kte, gravity_in, Rd_in, Rv_in, cp_in, ep_1_in, xlv_in, xlf_in, karman_in, &
                   presref_mid, nbot_shoc, ntop_shoc, length_fac)
    call mpas_log_write('!! returned from shoc_init !!')
 end subroutine bl_shocmf_init

 subroutine bl_shocmf_driver( &
    itimestep, dt_pbl, area_p, cda_p, lh_p, hfx_p, ht_p, theta_m_p, theta_p, &
    z_p, zmid_p, pres_p, pres2_p, uw_sfc_p, vw_sfc_p, w_p, pi_p, dse_p, tke_shmf_p, qw_p, &
    u_p, v_p, wthv_sec_p, kzh_p, kzm_p, qc_p, cldfrabl_p, hpbl_p, pblh_wthl_p, pblh_tke_p, &
    rino_p, shoc_mix_p, isotropy_p, tke_dissip_p, tke_prod_b_p, tke_prod_sh_p, w_sec_p, &
    thl_sec_p, qw_sec_p, qwthl_sec_p, wthl_sec_p, wqw_sec_p, wtke_sec_p, uw_sec_p, &
    vw_sec_p, w3_p, wqls_sec_p, brunt_p, shoc_ql2_p, mf_dry_a_p, mf_moist_a_p, mf_dry_w_p, &
    mf_moist_w_p, mf_dry_qt_p, mf_moist_qt_p, mf_dry_thl_p, mf_moist_thl_p, mf_dry_u_p, &
    mf_moist_u_p, mf_dry_v_p, mf_moist_v_p, mf_moist_qc_p, mf_thlflx_p, mf_qtflx_p, &
    mf_ae_p, mf_aw_p, mf_awthv_p, mf_awthl_p, mf_awqt_p, mf_awql_p, mf_awqi_p, mf_awu_p, &
    mf_awv_p, cfl_mf_p, mf_dry_freq_p, mf_moist_freq_p, plumeheight_p, plume_dry_height_p, &
    ztop_b_p, wthl_sec_ed_p, wthl_sec_mf_p, wqw_sec_ed_p, wqw_sec_mf_p, rho_p, t_p, qv_p, &
    qs_p, qi_p, qg_p, qr_p, rublten_p, rvblten_p, rthblten_p, rqvblten_p, rqcblten_p, &
    rqiblten_p, rtkeblten_p, &
    ids, ide, jds, jde, kds, kde, &
    ims, ime, jms, jme, kms, kme, &
    its, ite, jts, jte, kts, kte  &
 )
    ! Driver for the SHOC-MF module
    ! This subroutine is called from the main driver
    ! It calls the SHOC-MF module to calculate the turbulence
    ! and the EDMF module to calculate the subgrid-scale clouds
    integer, intent(in) :: itimestep, ids, ide, jds, jde, kds, kde, &
                           ims, ime, jms, jme, kms, kme, its, ite, &
                           jts, jte, kts, kte
    real(kind=RKIND), intent(in) :: dt_pbl
    real(kind=RKIND),dimension(:,:), intent(in) :: area_p, cda_p, lh_p, hfx_p, ht_p

    real(kind=RKIND),dimension(:,:), intent(inout) :: uw_sfc_p, vw_sfc_p, hpbl_p, pblh_wthl_p, pblh_tke_p, &
                                                      mf_dry_freq_p, mf_moist_freq_p, plumeheight_p, & 
                                                      plume_dry_height_p, ztop_b_p
    real(kind=RKIND),dimension(:,:,:), intent(inout) :: rino_p, shoc_mix_p, isotropy_p, &
        tke_dissip_p, tke_prod_b_p, tke_prod_sh_p, w_sec_p, thl_sec_p, qw_sec_p, &
        qwthl_sec_p, wthl_sec_p, wqw_sec_p, wtke_sec_p, uw_sec_p, vw_sec_p, w3_p, &
        wqls_sec_p, brunt_p, shoc_ql2_p, mf_dry_a_p, mf_moist_a_p, mf_dry_w_p, &
        mf_moist_w_p, mf_dry_qt_p, mf_moist_qt_p, mf_dry_thl_p, mf_moist_thl_p, &
        mf_dry_u_p, mf_moist_u_p, mf_dry_v_p, mf_moist_v_p, mf_moist_qc_p, mf_thlflx_p, &
        mf_qtflx_p, mf_ae_p, mf_aw_p, mf_awthv_p, mf_awthl_p, mf_awqt_p, mf_awql_p, &
        mf_awqi_p, mf_awu_p, mf_awv_p, cfl_mf_p, wthl_sec_ed_p, wthl_sec_mf_p, wqw_sec_ed_p, &
        wqw_sec_mf_p, dse_p, qw_p, u_p, v_p, wthv_sec_p, tke_shmf_p, kzh_p, kzm_p, qc_p, &
        cldfrabl_p, theta_p, t_p, theta_m_p, z_p, zmid_p, pres_p, pres2_p, w_p, pi_p, qv_p, &
        qs_p, qi_p, qg_p, qr_p, rho_p, rublten_p, rvblten_p, rthblten_p, rqvblten_p, rtkeblten_p, &
        rqcblten_p, rqiblten_p

    real(kind=RKIND), dimension(:,:,:), allocatable :: presdiff_p, thetal_p, qtracers, height_interface_p, height_mid_p
    real(kind=RKIND), dimension(:,:,:), allocatable :: um_p, vm_p, theta_mm_p, qvm_p, qcm_p, qim_p, qsm_p, tkem_p, w_zt
    real(kind=RKIND), dimension(:,:,:), allocatable :: inv_exner_p
    real(kind=RKIND), dimension(:,:)  , allocatable :: wtracer_sfc, dx_p, dy_p, hfx_converted_p, lh_converted_p

    integer :: num_qtracers, nadv, i, j, k
    real(kind=RKIND) :: shoc_timestep

    call mpas_log_write('!!!!!!!!!!! Entered in bl_shocmf_driver !!!!!!!!!!!!')

    tke_shmf_p = max(mintke, tke_shmf_p)

    allocate(                                  &
        um_p(ims:ime, kms:kme, jms:jme)      , &
        vm_p(ims:ime, kms:kme, jms:jme)      , &
        theta_mm_p(ims:ime, kms:kme, jms:jme), &
        qvm_p(ims:ime, kms:kme, jms:jme)     , &
        qcm_p(ims:ime, kms:kme, jms:jme)     , &
        tkem_p(ims:ime, kms:kme, jms:jme)      &
    )
    ! Deep copying the input variables
    um_p = u_p
    vm_p = v_p
    theta_mm_p = theta_m_p
    qvm_p = qv_p
    qcm_p = qc_p
    tkem_p = tke_shmf_p

    allocate(inv_exner_p(ims:ime, kms:kme, jms:jme))
    do i=ims,ime
    do k=kms,kme-1
    do j=jms,jme
        inv_exner_p(i,k,j) = 1.0_RKIND/pi_p(i,k,j) ! inv_exner_p = theta_p/t_p
    enddo
    enddo
    enddo

    call calculate_thetam_to_shocmf(theta_p, qv_p, qc_p, R_v, R_d, theta_m_p)

    allocate(thetal_p(ims:ime, kms:kme, jms:jme))
    call calculate_thetal(ime, jme, kme, theta_p, qc_p, inv_exner_p, thetal_p)

    allocate(presdiff_p(ims:ime, kms:kme, jms:jme))
    call calculate_presdiff(pres2_p, ims, ime, kms, kme, jms, jme, presdiff_p)

    call calculate_qw(qc_p, qv_p, qs_p, qi_p, qg_p, qr_p, qw_p)

    allocate(height_interface_p(ims:ime,kms:kme,jms:jme), height_mid_p(ims:ime,kms:kme,jms:jme))
    call calculate_heights(ims, ime, kms, kme, jms, jme, z_p, zmid_p, height_interface_p, height_mid_p)
    call calculate_dry_static_energy(ims, ime, kms, kme, jms, jme, t_p, height_interface_p, gravity, ht_p, cp, dse_p)
    
    call calculate_surface_momentum_fluxes(ime, kme, cda_p, u_p(:,1,:), v_p(:,1,:), uw_sfc_p, vw_sfc_p)
    ! TODO (GTF): w should be in thermodynamic grid 
    !call calculate_w_zt(ims, ime, kms, kme, jms, jme, w_p, gravity, rho_p, w_zt)

    allocate(dx_p(ims:ime, jms:jme), dy_p(ims:ime, jms:jme))
    call define_horizontal_discretization(ime, jme, area_p, dx_p, dy_p)

    call convert_heat_fluxes(ims, ime, jms, jme, hfx_p, lh_p, rho_p, cp, xlv, hfx_converted_p, lh_converted_p)

    nadv = ceiling(dt_pbl/max_shoc_timestep)
    if (nadv .gt. 1) then
        shoc_timestep = dt_pbl/nadv
    else
        shoc_timestep = dt_pbl
    endif

    print *, 'max timestep', max_shoc_timestep, 'dt_pbl', dt_pbl, 'shoc timestep', shoc_timestep, 'nadv', nadv

    ! If resulting host model time step and SHOC time step do not divide evenly
    ! into each other, have model throw a fit.
    if (mod(dt_pbl,shoc_timestep) .ne. 0) then
        stop "SHOC time step and HOST time step NOT compatible"
    endif

    num_qtracers = 0
    allocate(wtracer_sfc(ime, jme))
    wtracer_sfc = 0.0_RKIND
    allocate(qtracers(ime, kme, jme))
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
        kzm_p, qc_p, cldfrabl_p, theta_m_p, height_interface_p, &
        height_mid_p, pres_p, pres2_p, w_p, inv_exner_p, &
        thetal_p, presdiff_p &
    )
    ! As last position of the array is a extra position, to not letting a
    ! zero value, we are going to extrapolate it considering that
    ! the kme-1 value is the average among kme-2 and kme
    ! thetal_p(:,kme-1,:) = (thetal_p(:,kme-2,:) + thetal_p(:,kme,:))/2
    thetal_p(:,kme,:) = 2*thetal_p(:,kme-1,:) - thetal_p(:,kme-2,:)

    call mpas_log_write('!!!!!!!!!!! Calling SHOC-MF main module !!!!!!!!!!!!')
    call shoc_main( &
        itimestep, ime, kme-1, kme, shoc_timestep, nadv, &                                    ! Input
        dx_p(:,1), dy_p(:,1), theta_m_p(:,:,1),&                                              ! Input
        height_mid_p(:,:,1), height_interface_p(:,:,1), pres_p(:,:,1), pres2_p(:,:,1), &      ! Input
        presdiff_p(:,:,1), hfx_converted_p(:,1), lh_converted_p(:,1), &                       ! Input
        uw_sfc_p(:,1), vw_sfc_p(:,1), &                                                       ! Input
        wtracer_sfc(:,1), num_qtracers, w_p(:,:,1), &                                         ! Input
        inv_exner_p(:,:,1),ht_p(:,1), &                                                              ! Input
        dse_p(:,:,1), tke_shmf_p(:,:,1), thetal_p(:,:,1), qw_p(:,:,1), &                      ! Input/Ouput
        u_p(:,:,1), v_p(:,:,1), qtracers(:,:,:), &                                            ! Input/Output
        wthv_sec_p(:,:,1), kzh_p(:,:,1),kzm_p(:,:,1), &                                       ! Input/Output
        qc_p(:,:,1),cldfrabl_p(:,:,1), &                                                      ! Input/Output
        hpbl_p(:,1),pblh_wthl_p(:,1),pblh_tke_p(:,1), &                                       ! Output
        rino_p(:,:,1), &                                                                      ! Output
        shoc_mix_p(:,:,1), isotropy_p(:,:,1), &                                               ! Output (diagnostic)
        tke_dissip_p(:,:,1), tke_prod_b_p(:,:,1), tke_prod_sh_p(:,:,1), &                     ! Output (diagnostic)
        w_sec_p(:,:,1), thl_sec_p(:,:,1), qw_sec_p(:,:,1), qwthl_sec_p(:,:,1), &              ! Output (diagnostic)
        wthl_sec_p(:,:,1), wqw_sec_p(:,:,1), wtke_sec_p(:,:,1), &                             ! Output (diagnostic)
        uw_sec_p(:,:,1), vw_sec_p(:,:,1), w3_p(:,:,1), &                                      ! Output (diagnostic)
        wqls_sec_p(:,:,1),brunt_p(:,:,1),shoc_ql2_p(:,:,1), &                                 ! Output (diagnostic)
        mf_dry_a_p(:,:,1), mf_moist_a_p(:,:,1), &                                             ! Output (EDMF diagnostic)
        mf_dry_w_p(:,:,1), mf_moist_w_p(:,:,1), &                                             ! Output (EDMF diagnostic)
        mf_dry_qt_p(:,:,1), mf_moist_qt_p(:,:,1), &                                           ! Output (EDMF diagnostic)
        mf_dry_thl_p(:,:,1), mf_moist_thl_p(:,:,1), &                                         ! Output (EDMF diagnostic)
        mf_dry_u_p(:,:,1),   mf_moist_u_p(:,:,1), &                                           ! Output (EDMF diagnostic)
        mf_dry_v_p(:,:,1),   mf_moist_v_p(:,:,1), &                                           ! Output (EDMF diagnostic)
        mf_moist_qc_p(:,:,1), &                                                               ! Output (EDMF diagnostic)
        mf_thlflx_p(:,:,1), mf_qtflx_p(:,:,1), &                                              ! Output (EDMF diagnostic)
        mf_ae_p(:,:,1), mf_aw_p(:,:,1), &                                                     ! Output (EDMF diagnostic)
        mf_awthv_p(:,:,1), &                                                                  ! Output (EDMF diagnostic)
        mf_awthl_p(:,:,1), mf_awqt_p(:,:,1), &                                                ! Output (EDMF diagnostic)
        mf_awql_p(:,:,1), mf_awqi_p(:,:,1), &                                                 ! Output (EDMF diagnostic)
        mf_awu_p(:,:,1), mf_awv_p(:,:,1), cfl_mf_p(:,:,1), &                                  ! Output (EDMF diagnostic)
        mf_dry_freq_p(:,1), mf_moist_freq_p(:,1), plumeheight_p(:,1), &                       ! Output (EDMF diagnostic)
        plume_dry_height_p(:,1), ztop_b_p(:,1), &                                             ! Output (EDMF diagnostic)
        wthl_sec_ed_p(:,:,1), wthl_sec_mf_p(:,:,1), &                                         ! Output (EDMF diagnostic)
        wqw_sec_ed_p(:,:,1), wqw_sec_mf_p(:,:,1))                                             ! Output (EDMF diagnostic)

    call mpas_log_write('!!!!!!!!!!! Returned from SHOC-MF main module !!!!!!!!!!!!')
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
        kzm_p, qc_p, cldfrabl_p, theta_m_p, height_interface_p, &
        height_mid_p, pres_p, pres2_p, w_p, inv_exner_p, &
        thetal_p, presdiff_p &
    )
    ! As last position of the array is a extra position,
    ! we are going to extrapolate it considering that
    ! the kme-1 value is the average among kme-2 and kme
    ! thetal_p(:,kme-1,:) = (thetal_p(:,kme-2,:) + thetal_p(:,kme,:))/2
    thetal_p(:,kme,:) = 2*thetal_p(:,kme-1,:) - thetal_p(:,kme-2,:)

    call calculate_theta(thetal_p, qc_p, inv_exner_p, theta_p)
    call calculate_qv(qw_p, qc_p, qs_p, qi_p, qg_p, qr_p, qv_p)
    call calculate_thetam_to_mpas(theta_p, qv_p, R_v, R_d, theta_m_p)

    ! Calculating tendencies
    rublten_p = (u_p - um_p)/dt_pbl
    rvblten_p = (v_p - vm_p)/dt_pbl
    rthblten_p = (theta_m_p - theta_mm_p)/dt_pbl
    rqvblten_p = (qv_p - qvm_p)/dt_pbl
    rqcblten_p = (qc_p - qcm_p)/dt_pbl
    rtkeblten_p = (tke_shmf_p - tkem_p)/dt_pbl

    rqiblten_p = 0.0_RKIND  ! TODO (GTF): qi_p and qs_p are not calculated by shocmf. Only qc_p and qw_p
    !rqsblten_p = 0.0_RKIND
    !rniblten_p = 0.0_RKIND ! TODO (GTF): Is it possible to calculate it?

    deallocate(presdiff_p, thetal_p, qtracers, wtracer_sfc)
    deallocate(um_p, vm_p, theta_mm_p, qvm_p, qcm_p, tkem_p)
    deallocate(hfx_converted_p, lh_converted_p)
    deallocate(dx_p, dy_p)
    deallocate(height_interface_p, height_mid_p)
    deallocate(inv_exner_p)
    call mpas_log_write('!!!!!!!!!!! Exiting bl_shocmf_driver !!!!!!!!!!!!')
  
 end subroutine bl_shocmf_driver

end module module_bl_shocmf
