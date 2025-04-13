module module_bl_shocmf
 use mpas_atmphys_constants
 use shoc, only: shoc_init, shoc_main, length_fac, mintke
 use mpas_pool_routines


 implicit none
 private
 public:: bl_shocmf_init, bl_shocmf_driver

 real(kind=RKIND) :: max_shoc_timestep = 300.0  ! Default SHOC timestep set in namelist

 contains

 subroutine calculate_surface_momentum_flux(cda_p, velocity_p, surface_momentum_flux)
    ! Calculate the surface momentum flux
    ! Reference: Kubota's class notes
    ! cda_p is the drag coefficient
    ! velocity_p is the velocity (zonal or meridional) at the surface
    ! surface_momentum_flux is the surface momentum flux
    real(kind=RKIND), dimension(:,:), intent(in) :: cda_p, velocity_p
    real(kind=RKIND), dimension(:,:), intent(out) :: surface_momentum_flux

    surface_momentum_flux(:,:) = -cda_p(:,:)*velocity_p(:,:)*velocity_p(:,:)
 end subroutine calculate_surface_momentum_flux


 subroutine calculate_dry_static_energy(t, z, gravity_in, dse)
    ! Calculate the dry static energy considering the terrain height.
    ! This is the calculation waited by the SHOCMF module. 
    ! Reference: https://glossary.ametsoc.org/wiki/Dry_static_energy
    real(kind=RKIND), dimension(:,:,:), intent(in) :: t, z
    real(kind=RKIND), intent(in) :: gravity_in
    real(kind=RKIND), dimension(:,:,:), intent(out) :: dse
  
    dse(:,:,:) = cp*t(:,:,:) + gravity_in*z(:,:,:)
 end subroutine calculate_dry_static_energy


 subroutine calculate_presdiff(pres_levels, its, ite, kts, kte, jts, jte, presdiff_levels)
    ! Calculate the pressure difference between the interfaces
    ! pres_levels is the pressure at the levels
    ! presdiff_levels is the pressure difference between the interfaces
    real(kind=RKIND), dimension(:,:,:), intent(in) :: pres_levels
    integer, intent(in) :: its, ite, kts, kte, jts, jte
    real(kind=RKIND), dimension(:,:,:), intent(out) :: presdiff_levels

    integer :: i, j, k
   
    do j = jts,jte
    do k = kts,kte-1
    do i = its,ite
        presdiff_levels(i,k,j) = pres_levels(i,k+1,j) - pres_levels(i,k,j)
    enddo
    enddo
    enddo
 end subroutine calculate_presdiff


 subroutine calculate_thetal(theta, temperature, ql, thetal)
    ! Calculate liquid water potential temperature as SHOC calculates it
    real(kind=RKIND), dimension(:,:,:), intent(in) :: theta, temperature, ql
    real(kind=RKIND), dimension(:,:,:), intent(out) :: thetal

    thetal(:,:,:) = theta(:,:,:) - (xlv/cp)*ql(:,:,:)
 end subroutine calculate_thetal


 subroutine calculate_theta(thetal, temperature, ql, theta)
    ! Calculate potential temperature from liquid water potential temperature.
    ! Inverse function of calculate_thetal.
    real(kind=RKIND), dimension(:,:,:), intent(in) :: thetal, temperature, ql
    real(kind=RKIND), dimension(:,:,:), intent(out) :: theta

    theta(:,:,:) = thetal(:,:,:) + (xlv/cp)*ql(:,:,:)
end subroutine calculate_theta


 subroutine calculate_qw(qc, qv, qs, qi, qg, qr, qw)
    ! Calculate total water mixing ratio
    real(kind=RKIND), dimension(:,:,:), intent(in) :: qc, qv, qs, qi, qg, qr
    real(kind=RKIND), dimension(:,:,:), intent(out) :: qw

    qw(:,:,:) = qc(:,:,:) + qv(:,:,:) + qs(:,:,:) + qi(:,:,:) + qg(:,:,:) + qr(:,:,:)
 end subroutine calculate_qw 


 subroutine flip_vertical_coordinates(&
    kte, rino_p, shoc_mix_p, isotropy_p, tke_dissip_p, tke_prod_b_p, &
    tke_prod_sh_p, w_sec_p, thl_sec_p, qw_sec_p, qwthl_sec_p, &
    wthl_sec_p, wqw_sec_p, wtke_sec_p, uw_sec_p, vw_sec_p, w3_p, &
    wqls_sec_p, brunt_p, shoc_ql2_p, mf_dry_a_p, mf_moist_a_p, &
    mf_dry_w_p, mf_moist_w_p, mf_dry_qt_p, mf_moist_qt_p, mf_dry_thl_p, &
    mf_moist_thl_p, mf_dry_u_p, mf_moist_u_p, mf_dry_v_p, mf_moist_v_p, &
    mf_moist_qc_p, mf_thlflx_p, mf_qtflx_p, mf_ae_p, mf_aw_p, &
    mf_awthv_p, mf_awthl_p, mf_awqt_p, mf_awql_p, mf_awqi_p, mf_awu_p, &
    mf_awv_p, cfl_mf_p, wthl_sec_ed_p, wthl_sec_mf_p, wqw_sec_ed_p, &
    wqw_sec_mf_p, dse_p, qw_p, u_p, v_p, wthv_sec_p, tke_shmf_p, kzh_p, &
    kzm_p, qc_p, cldfrabl_p, theta_p, t_p, theta_m_p, z_p, zgrid_p, pres_p, &
    pres2_p, w_p, pi_p, qv_p, qs_p, qi_p, qg_p, qr_p, rho_p &
  )
    integer, intent(in) :: kte
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
                                                        theta_p, t_p, theta_m_p, z_p, zgrid_p, pres_p, &
                                                        pres2_p, w_p, pi_p, qv_p, qs_p, &
                                                        qi_p, qg_p, qr_p, rho_p

    ! Flip the vertical coordinates of the SHOC-MF variables
    ! This is necessary because SHOC-MF expects the vertical
    ! coordinate starting from the top of the atmosphere.
    ! This subroutine is called before and after the SHOC-MF call
    rino_p = rino_p(:, kte:1:-1, :)
    shoc_mix_p = shoc_mix_p(:, kte:1:-1, :)
    isotropy_p = isotropy_p(:, kte:1:-1, :)
    tke_dissip_p = tke_dissip_p(:, kte:1:-1, :)
    tke_prod_b_p = tke_prod_b_p(:, kte:1:-1, :)
    tke_prod_sh_p = tke_prod_sh_p(:, kte:1:-1, :)
    w_sec_p = w_sec_p(:, kte:1:-1, :)
    thl_sec_p = thl_sec_p(:, kte:1:-1, :)
    qw_sec_p = qw_sec_p(:, kte:1:-1, :)
    qwthl_sec_p = qwthl_sec_p(:, kte:1:-1, :)
    wthl_sec_p = wthl_sec_p(:, kte:1:-1, :)
    wqw_sec_p = wqw_sec_p(:, kte:1:-1, :)
    wtke_sec_p = wtke_sec_p(:, kte:1:-1, :)
    uw_sec_p = uw_sec_p(:, kte:1:-1, :)
    vw_sec_p = vw_sec_p(:, kte:1:-1, :)
    w3_p = w3_p(:, kte:1:-1, :)
    wqls_sec_p = wqls_sec_p(:, kte:1:-1, :)
    brunt_p = brunt_p(:, kte:1:-1, :)
    shoc_ql2_p = shoc_ql2_p(:, kte:1:-1, :)
    mf_dry_a_p = mf_dry_a_p(:, kte:1:-1, :)
    mf_moist_a_p = mf_moist_a_p(:, kte:1:-1, :)
    mf_dry_w_p = mf_dry_w_p(:, kte:1:-1, :)
    mf_moist_w_p = mf_moist_w_p(:, kte:1:-1, :)
    mf_dry_qt_p = mf_dry_qt_p(:, kte:1:-1, :)
    mf_moist_qt_p = mf_moist_qt_p(:, kte:1:-1, :)
    mf_dry_thl_p = mf_dry_thl_p(:, kte:1:-1, :)
    mf_moist_thl_p = mf_moist_thl_p(:, kte:1:-1, :)
    mf_dry_u_p = mf_dry_u_p(:, kte:1:-1, :)
    mf_moist_u_p = mf_moist_u_p(:, kte:1:-1, :)
    mf_dry_v_p = mf_dry_v_p(:, kte:1:-1, :)
    mf_moist_v_p = mf_moist_v_p(:, kte:1:-1, :)
    mf_moist_qc_p = mf_moist_qc_p(:, kte:1:-1, :)
    mf_thlflx_p = mf_thlflx_p(:, kte:1:-1, :)
    mf_qtflx_p = mf_qtflx_p(:, kte:1:-1, :)
    mf_ae_p = mf_ae_p(:, kte:1:-1, :)
    mf_aw_p = mf_aw_p(:, kte:1:-1, :)
    mf_awthv_p = mf_awthv_p(:, kte:1:-1, :)
    mf_awthl_p = mf_awthl_p(:, kte:1:-1, :)
    mf_awqt_p = mf_awqt_p(:, kte:1:-1, :)
    mf_awql_p = mf_awql_p(:, kte:1:-1, :)
    mf_awqi_p = mf_awqi_p(:, kte:1:-1, :)
    mf_awu_p = mf_awu_p(:, kte:1:-1, :)
    mf_awv_p = mf_awv_p(:, kte:1:-1, :)
    cfl_mf_p = cfl_mf_p(:, kte:1:-1, :)
    wthl_sec_ed_p = wthl_sec_ed_p(:, kte:1:-1, :)
    wthl_sec_mf_p = wthl_sec_mf_p(:, kte:1:-1, :)
    wqw_sec_ed_p = wqw_sec_ed_p(:, kte:1:-1, :)
    wqw_sec_mf_p = wqw_sec_mf_p(:, kte:1:-1, :)
    dse_p = dse_p(:, kte:1:-1, :)
    qw_p = qw_p(:, kte:1:-1, :)
    u_p = u_p(:, kte:1:-1, :)
    v_p = v_p(:, kte:1:-1, :)
    wthv_sec_p = wthv_sec_p(:, kte:1:-1, :)
    tke_shmf_p = tke_shmf_p(:, kte:1:-1, :)
    kzh_p = kzh_p(:, kte:1:-1, :)
    kzm_p = kzm_p(:, kte:1:-1, :)
    qc_p = qc_p(:, kte:1:-1, :)
    cldfrabl_p = cldfrabl_p(:, kte:1:-1, :)
    theta_p = theta_p(:, kte:1:-1, :)
    t_p = t_p(:, kte:1:-1, :)
    theta_m_p = theta_m_p(:, kte:1:-1, :)
    z_p = z_p(:, kte:1:-1, :)
    zgrid_p = zgrid_p(:, kte:1:-1, :)
    pres_p = pres_p(:, kte:1:-1, :)
    pres2_p = pres2_p(:, kte:1:-1, :)
    w_p = w_p(:, kte:1:-1, :)
    pi_p = pi_p(:, kte:1:-1, :)
    qv_p = qv_p(:, kte:1:-1, :)
    qs_p = qs_p(:, kte:1:-1, :)
    qi_p = qi_p(:, kte:1:-1, :)
    qg_p = qg_p(:, kte:1:-1, :)
    qr_p = qr_p(:, kte:1:-1, :)
    rho_p = rho_p(:, kte:1:-1, :)
 end subroutine flip_vertical_coordinates

 
 subroutine define_horizontal_discretization(ite, jte, area_p, dx_p, dy_p)
    ! Define the horizontal discretization
    ! dx_p is the horizontal grid spacing in the x direction
    ! dy_p is the horizontal grid spacing in the y direction
    integer, intent(in) :: ite, jte
    real(kind=RKIND), dimension(:,:), intent(in) :: area_p
    real(kind=RKIND), dimension(:,:), allocatable, intent(out) :: dx_p, dy_p

    allocate(dx_p(ite, jte), dy_p(ite, jte))
    ! dx_p and dy_p are calculated as sqrt(area_p) in the main driver
    dx_p = sqrt(area_p)
    dy_p = sqrt(area_p)
 end subroutine define_horizontal_discretization

 
 subroutine convert_heat_fluxes(its, ite, jts, jte, hfx_p, lh_p, rho_p, hfx_converted_p, lh_converted_p)
    ! Convert the surface heat fluxes from W/m^2 to m/s
    ! hfx_p is the sensible heat flux
    ! lh_p is the latent heat flux
    integer, intent(in) :: its, ite, jts, jte
    real(kind=RKIND), dimension(:,:), intent(in) :: hfx_p, lh_p
    real(kind=RKIND), dimension(:,:,:), intent(in) :: rho_p
    real(kind=RKIND), dimension(:,:), allocatable, intent(out) :: hfx_converted_p, lh_converted_p

    integer :: i, j

    allocate(hfx_converted_p(ite, jte), lh_converted_p(ite, jte))
    do j = jts,jte
    do i = its,ite
        ! Converting hfx_p from W/m^2 to K*m/s
        hfx_converted_p(i,j) = hfx_p(i,j)/(cp*rho_p(i,j,1))
      
        ! Converting lh_p from W/m^2 to kg/kg m/s
        lh_converted_p(i,j) = lh_p(i,j)/(xlv*rho_p(i,j,1))
    enddo
    enddo
 end subroutine convert_heat_fluxes

 subroutine calculate_thetam_from_thetal(thetal, ql, qv, ep_1_in, xlv_in, cp_in, theta_m)
    ! Calculate moist potential temperature from liquid water potential temperature.
    ! Inverse function of calculate_thetal.
    real(kind=RKIND), dimension(:,:,:), intent(in) :: thetal, ql, qv
   real(kind=RKIND), intent(in) :: ep_1_in, xlv_in, cp_in
    real(kind=RKIND), dimension(:,:,:), intent(out) :: theta_m

    theta_m(:,:,:) = (thetal(:,:,:) + (xlv_in/cp_in)*ql(:,:,:))*(1.0_RKIND + ep_1_in*qv(:,:,:) - ql(:,:,:))
 end subroutine calculate_thetam_from_thetal

 subroutine bl_shocmf_init( &
    kte, ims, ime, kms, kme, jms, jme, gravity_in, Rd_in, Rv_in, cp_in, ep_1_in, xlv_in, &
    xlf_in, karman_in, pres_p, tke_shmf_p, wthv_sec_p, &
    kzm_p, kzh_p &
 )
    ! Initialize the SHOC-MF module
    ! This subroutine is called from the main driver
    integer, intent(in) :: kte, ims, ime, kms, kme, jms, jme
    real(RKIND), intent(in) :: gravity_in, Rd_in, Rv_in, cp_in, ep_1_in, xlv_in, xlf_in, &
                               karman_in
    real(RKIND), dimension(:,:,:), intent(in) :: pres_p
    real(RKIND), dimension(:,:,:), allocatable, intent(inout) :: tke_shmf_p, wthv_sec_p, &
                                                    kzm_p, kzh_p

    integer:: ntop_shoc, nbot_shoc

    ntop_shoc = 1    ! if >1, must be <= nbot_molec
    nbot_shoc = kte ! currently always pver 
    
    call mpas_log_write('Before allocation of SHOC-MF variables')
    if(.not.allocated(tke_shmf_p)) allocate(tke_shmf_p(ims:ime,kms:kme,jms:jme))
    if(.not.allocated(wthv_sec_p)) allocate(wthv_sec_p(ims:ime,kms:kme,jms:jme))
    if(.not.allocated(kzh_p)) allocate(kzh_p(ims:ime,kms:kme,jms:jme))
    if(.not.allocated(kzm_p)) allocate(kzm_p(ims:ime,kms:kme,jms:jme))
    call mpas_log_write('After allocation of SHOC-MF variables')

    tke_shmf_p(:,:,:) = mintke
    wthv_sec_p(:,:,:) = 0._RKIND
    kzm_p(:,:,:) = 0._RKIND
    kzh_p(:,:,:) = 0._RKIND

    call mpas_log_write('Calling shoc_init from SHOCMF module')
    ! TODO (GTF): Review pres_p and its usage in the shoc_init method
    call shoc_init(kte, gravity_in, Rd_in, Rv_in, cp_in, ep_1_in, xlv_in, xlf_in, karman_in, &
                  pres_p, nbot_shoc, ntop_shoc, length_fac)
    call mpas_log_write('!! returned from shoc_init !!')
 end subroutine bl_shocmf_init


 subroutine bl_shocmf_driver( &
    itimestep, dt_pbl, area_p, cda_p, lh_p, hfx_p, ht_p, theta_m_p, theta_p, &
    z_p, zgrid_p, pres_p, pres2_p, uw_sfc_p, vw_sfc_p, w_p, pi_p, dse_p, tke_shmf_p, qw_p, &
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
    rqiblten_p, rqsblten_p, rniblten_p, rtkeblten_p, its, ite, jts, jte, kts, kte &
 )
    ! Driver for the SHOC-MF module
    ! This subroutine is called from the main driver
    ! It calls the SHOC-MF module to calculate the turbulence
    ! and the EDMF module to calculate the subgrid-scale clouds
    integer, intent(in) :: itimestep, its, ite, jts, jte, kts, kte
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
        cldfrabl_p, theta_p, t_p, theta_m_p, z_p, zgrid_p, pres_p, pres2_p, w_p, pi_p, qv_p, &
        qs_p, qi_p, qg_p, qr_p, rho_p, rublten_p, rvblten_p, rthblten_p, rqvblten_p, rtkeblten_p, &
        rqcblten_p, rqiblten_p, rqsblten_p, rniblten_p

    real(kind=RKIND), dimension(:,:,:), allocatable :: presdiff_p, thetal_p, qtracers
    real(kind=RKIND), dimension(:,:,:), allocatable :: um_p, vm_p, theta_mm_p, qvm_p, qcm_p, qim_p, qsm_p, nim_p, tkem_p 
    real(kind=RKIND), dimension(:,:), allocatable :: wtracer_sfc, dx_p, dy_p, hfx_converted_p, lh_converted_p

    integer :: num_qtracers, nadv, i, j
    real(kind=RKIND) :: shoc_timestep

    call mpas_log_write('!!!!!!!!!!! Entered in bl_shocmf_driver !!!!!!!!!!!!')

    allocate(um_p(ite, kte, jte)    , &
             vm_p(ite, kte, jte)    , &
             theta_mm_p(ite, kte, jte), &
             qvm_p(ite, kte, jte)   , &
             qcm_p(ite, kte, jte)   , &
             qim_p(ite, kte, jte)   , &
             qsm_p(ite, kte, jte)   , &
             nim_p(ite, kte, jte)   , &
             tkem_p(ite, kte, jte))
    ! Deep copying the input variables
    um_p = u_p
    vm_p = v_p
    theta_mm_p = theta_m_p
    qvm_p = qv_p
    qcm_p = qc_p
    qim_p = qi_p
    qsm_p = qs_p
    nim_p = qi_p
    tkem_p = tke_shmf_p

    call flip_vertical_coordinates(&
        kte, rino_p, shoc_mix_p, isotropy_p, tke_dissip_p, tke_prod_b_p, &
        tke_prod_sh_p, w_sec_p, thl_sec_p, qw_sec_p, qwthl_sec_p, &
        wthl_sec_p, wqw_sec_p, wtke_sec_p, uw_sec_p, vw_sec_p, w3_p, &
        wqls_sec_p, brunt_p, shoc_ql2_p, mf_dry_a_p, mf_moist_a_p, &
        mf_dry_w_p, mf_moist_w_p, mf_dry_qt_p, mf_moist_qt_p, mf_dry_thl_p, &
        mf_moist_thl_p, mf_dry_u_p, mf_moist_u_p, mf_dry_v_p, mf_moist_v_p, &
        mf_moist_qc_p, mf_thlflx_p, mf_qtflx_p, mf_ae_p, mf_aw_p, &
        mf_awthv_p, mf_awthl_p, mf_awqt_p, mf_awql_p, mf_awqi_p, mf_awu_p, &
        mf_awv_p, cfl_mf_p, wthl_sec_ed_p, wthl_sec_mf_p, wqw_sec_ed_p, &
        wqw_sec_mf_p, dse_p, qw_p, u_p, v_p, wthv_sec_p, tke_shmf_p, kzh_p, &
        kzm_p, qc_p, cldfrabl_p, theta_p, t_p, theta_m_p, z_p, zgrid_p, pres_p, &
        pres2_p, w_p, pi_p, qv_p, qs_p, qi_p, qg_p, qr_p, rho_p &
    )
    call calculate_thetal(theta_p, t_p, qc_p, thetal_p)
    call calculate_presdiff(pres_p, its, ite, kts, kte, jts, jte, presdiff_p)
    call calculate_qw(qc_p, qv_p, qs_p, qi_p, qg_p, qr_p, qw_p)
    call calculate_dry_static_energy(t_p, z_p, gravity, dse_p)
    call calculate_surface_momentum_flux(cda_p, u_p(:, 1, :), uw_sfc_p)
    call calculate_surface_momentum_flux(cda_p, v_p(:, 1, :), vw_sfc_p)
    call define_horizontal_discretization(ite, jte, area_p, dx_p, dy_p)

    call convert_heat_fluxes(its, ite, jts, jte, hfx_p, lh_p, rho_p, hfx_converted_p, lh_converted_p)

    nadv = ceiling(dt_pbl/max_shoc_timestep)
    if (nadv .gt. 1) then
        shoc_timestep = dt_pbl/nadv
    else
        shoc_timestep = dt_pbl
    endif

    ! If resulting host model time step and SHOC time step do not divide evenly
    ! into each other, have model throw a fit.
    if (mod(dt_pbl,shoc_timestep) .ne. 0) then
        stop "SHOC time step and HOST time step NOT compatible"
    endif

    num_qtracers = 0
    allocate(wtracer_sfc(ite, jte))
    wtracer_sfc = 0.0_RKIND
    allocate(qtracers(ite, kte, jte))
    qtracers = 0.0_RKIND

    call mpas_log_write('!!!!!!!!!!! Calling SHOC-MF main module !!!!!!!!!!!!')
    call shoc_main( &
        itimestep, ite, kte, kte+1, shoc_timestep, nadv, &                                    ! Input
        dx_p(:,1), dy_p(:,1), theta_m_p(:,:,1),&                                              ! Input
        z_p(:,:,1), zgrid_p(:,:,1), pres_p(:,:,1), pres2_p(:,:,1), &                          ! Input
        presdiff_p(:,:,1), hfx_converted_p(:,1), lh_converted_p(:,1), &                       ! Input
        uw_sfc_p(:,1), vw_sfc_p(:,1), &                                                       ! Input
        wtracer_sfc(:,1), num_qtracers, w_p(:,:,1), &                                         ! Input
        pi_p(:,:,1),ht_p(:,1), &                                                              ! Input
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
    call calculate_theta(thetal_p, t_p, qc_p, theta_p)
    call calculate_thetam_from_thetal(thetal_p, qc_p, qv_p, ep_1, xlv, cp, theta_m_p)

    call flip_vertical_coordinates(&
        kte, rino_p, shoc_mix_p, isotropy_p, tke_dissip_p, tke_prod_b_p, &
        tke_prod_sh_p, w_sec_p, thl_sec_p, qw_sec_p, qwthl_sec_p, &
        wthl_sec_p, wqw_sec_p, wtke_sec_p, uw_sec_p, vw_sec_p, w3_p, &
        wqls_sec_p, brunt_p, shoc_ql2_p, mf_dry_a_p, mf_moist_a_p, &
        mf_dry_w_p, mf_moist_w_p, mf_dry_qt_p, mf_moist_qt_p, mf_dry_thl_p, &
        mf_moist_thl_p, mf_dry_u_p, mf_moist_u_p, mf_dry_v_p, mf_moist_v_p, &
        mf_moist_qc_p, mf_thlflx_p, mf_qtflx_p, mf_ae_p, mf_aw_p, &
        mf_awthv_p, mf_awthl_p, mf_awqt_p, mf_awql_p, mf_awqi_p, mf_awu_p, &
        mf_awv_p, cfl_mf_p, wthl_sec_ed_p, wthl_sec_mf_p, wqw_sec_ed_p, &
        wqw_sec_mf_p, dse_p, qw_p, u_p, v_p, wthv_sec_p, tke_shmf_p, kzh_p, &
        kzm_p, qc_p, cldfrabl_p, theta_p, t_p, theta_m_p, z_p, zgrid_p, pres_p, &
        pres2_p, w_p, pi_p, qv_p, qs_p, qi_p, qg_p, qr_p, rho_p &
    )

    ! Calculating tendencies
    rublten_p = (u_p - um_p)/dt_pbl
    rvblten_p = (v_p - vm_p)/dt_pbl
    rthblten_p = (theta_m_p - theta_mm_p)/dt_pbl
    rqvblten_p = (qv_p - qvm_p)/dt_pbl
    rqcblten_p = (qc_p - qcm_p)/dt_pbl
    !rqiblten_p = (qi_p - qim_p)/dt_pbl ! TODO: qi_p and qs_p are not calculated by shocmf. Only qc_p and qw_p
    !rqsblten_p = (qs_p - qsm_p)/dt_pbl
    rtkeblten_p = (tke_shmf_p - tkem_p)/dt_pbl
    !rniblten_p = 0.0_RKIND ! TODO: Is it possible to calculate it?
  
 end subroutine bl_shocmf_driver

end module module_bl_shocmf
