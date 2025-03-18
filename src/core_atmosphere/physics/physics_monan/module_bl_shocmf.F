module module_bl_shocmf
 use mpas_atmphys_constants

 use shoc, only: shoc_init, shoc_main, length_fac, mintke


 implicit none
 private
 public:: bl_shocmf_init, bl_shocmf_driver

 real(kind=RKIND) :: shoc_timestep = 150.0  ! Default SHOC timestep set in namelist

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


 subroutine calculate_presdiff(pres_levels, jts, jte, kts, kte, its, ite, presdiff_levels)
   ! Calculate the pressure difference between the interfaces
   ! pres_levels is the pressure at the levels
   ! presdiff_levels is the pressure difference between the interfaces
   real(kind=RKIND), dimension(:,:,:), intent(in) :: pres_levels
   integer, intent(in) :: jts, jte, kts, kte, its, ite
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
   ! Calculate liquid water potential temperature
   ! Reference: https://glossary.ametsoc.org/wiki/Liquid_water_potential_temperature
   real(kind=RKIND), dimension(:,:,:), intent(in) :: theta, temperature, ql
   real(kind=RKIND), dimension(:,:,:), intent(out) :: thetal

   thetal(:,:,:) = theta(:,:,:)*(1.0_RKIND - (xlv/cp)*ql(:,:,:)/temperature(:,:,:))
 end subroutine calculate_thetal


 subroutine calculate_theta(thetal, temperature, ql, theta)
   ! Calculate potential temperature. Inverse function of calculate_thetal.
   real(kind=RKIND), dimension(:,:,:), intent(in) :: thetal, temperature, ql
   real(kind=RKIND), dimension(:,:,:), intent(out) :: theta

   theta(:,:,:) = thetal(:,:,:)/(1.0_RKIND - (xlv/cp)*ql(:,:,:)/temperature(:,:,:))
end subroutine calculate_theta


 subroutine calculate_qw(qc, qv, qs, qi, qg, qr, qw)
   ! Calculate total water mixing ratio
   real(kind=RKIND), dimension(:,:,:), intent(in) :: qc, qv, qs, qi, qg, qr
   real(kind=RKIND), dimension(:,:,:), intent(out) :: qw

   qw(:,:,:) = qc(:,:,:) + qv(:,:,:) + qs(:,:,:) + qi(:,:,:) + qg(:,:,:) + qr(:,:,:)
 end subroutine calculate_qw 


 subroutine bl_shocmf_init( &
      kte, gravity_in, Rd_in, Rv_in, cp_in, ep_1_in, xlv_in, &
      xlf_in, karman_in, pres_p, tke_shmf_p, wthv_sec_p, &
      tk_p, tkh_p &
  )
   ! Initialize the SHOC-MF module
   ! This subroutine is called from the main driver
   integer, intent(in) :: kte
   real(RKIND), intent(in) :: gravity_in, Rd_in, Rv_in, cp_in, ep_1_in, xlv_in, xlf_in, &
                               karman_in
   real(RKIND), dimension(:,:,:), intent(in) :: pres_p
   real(RKIND), dimension(:,:,:), intent(inout) :: tke_shmf_p, wthv_sec_p, &
                                                    tk_p, tkh_p

   integer:: ntop_shoc, nbot_shoc

   ntop_shoc = 1    ! if >1, must be <= nbot_molec
   nbot_shoc = kte ! currently always pver

   tke_shmf_p(:,:,:) = mintke
   wthv_sec_p(:,:,:) = 0._RKIND
   tk_p(:,:,:) = 0._RKIND
   tkh_p(:,:,:) = 0._RKIND
  
   ! Revise pres
   call shoc_init(kte, gravity_in, Rd_in, Rv_in, cp_in, ep_1_in, xlv_in, xlf_in, karman_in, &
                  pres_p, nbot_shoc, ntop_shoc, length_fac)
 end subroutine bl_shocmf_init


 subroutine bl_shocmf_driver( &
      itimestep, dt_dyn, area_p, theta_m, theta_p, z_p, zgrid_p, pres_p, pres2_p, &
      hfx_p, lh_p, uw_sfc_p, vw_sfc_p, w_p, pi_p, ht_p, dse_p, tke_shmf_p, qw_p, &
      u_p, v_p, wthv_sec_p, kzh_p, kzm_p, qc_p, cldfrabl_p, hpbl_p, pblh_wthl_p, pblh_tke_p, &
      rino_p, shoc_mix_p, isotropy_p, tke_dissip_p, tke_prod_b_p, tke_prod_sh_p, w_sec_p, &
      thl_sec_p, qw_sec_p, qwthl_sec_p, wthl_sec_p, wqw_sec_p, wtke_sec_p, uw_sec_p, &
      vw_sec_p, w3_p, wqls_p, brunt_p, shoc_ql2_p, mf_dry_a_p, mf_moist_a_p, mf_dry_w_p, &
      mf_moist_w_p, mf_dry_qt_p, mf_moist_qt_p, mf_dry_thl_p, mf_moist_thl_p, mf_dry_u_p, &
      mf_moist_u_p, mf_dry_v_p, mf_moist_v_p, mf_moist_qc_p, mf_thlflx_p, mf_qtflx_p, &
      mf_ae_p, mf_aw_p, mf_awthv_p, mf_awthl_p, mf_awqt_p, mf_awql_p, mf_awqi_p, mf_awu_p, &
      mf_awv_p, cfl_mf_p, mf_dry_freq_p, mf_moist_freq_p, plumeheight_p, plume_dry_height_p, &
      ztop_b_p, wthl_sec_ed_p, wthl_sec_mf_p, wqw_sec_ed_p, wqw_sec_mf_p, rho_p, cda_p, t_p, &
      qv_p, qs_p, qi_p, qg_p, qr_p, its, ite, jts, jte, kts, kte &
   )
   ! Driver for the SHOC-MF module
   ! This subroutine is called from the main driver
   ! It calls the SHOC-MF module to calculate the turbulence
   ! and the EDMF module to calculate the subgrid-scale clouds
   integer, intent(in) :: itimestep, its, ite, jts, jte, kts, kte
   real(kind=RKIND), intent(in) :: dt_dyn
   real(kind=RKIND),dimension(:,:), intent(in) :: area_p, mf_dry_freq_p, mf_moist_freq_p, plumeheight_p, &
                                                  plume_dry_height_p, ztop_b_p, pblh_wthl_p, pblh_tke_p, &
                                                  cda_p, lh_p, hfx_p, ht_p, hpbl_p
   real(kind=RKIND),dimension(:,:,:), intent(in) :: t_p, theta_m, z_p, zgrid_p, pres_p, &
                                                    pres2_p, w_p, pi_p, &
                                                    tke_shmf_p, &
                                                    rino_p, &
                                                    shoc_mix_p, isotropy_p, tke_dissip_p, tke_prod_b_p, &
                                                    tke_prod_sh_p, w_sec_p, thl_sec_p, qw_sec_p, &
                                                    qwthl_sec_p, wthl_sec_p, wqw_sec_p, wtke_sec_p, &
                                                    uw_sec_p, vw_sec_p, w3_p, wqls_p, brunt_p, shoc_ql2_p, &
                                                    mf_dry_a_p, mf_moist_a_p, mf_dry_w_p, mf_moist_w_p, &
                                                    mf_dry_qt_p, mf_moist_qt_p, mf_dry_thl_p, mf_moist_thl_p, &
                                                    mf_dry_u_p, mf_moist_u_p, mf_dry_v_p, mf_moist_v_p, &
                                                    mf_moist_qc_p, mf_thlflx_p, mf_qtflx_p, mf_ae_p, &
                                                    mf_aw_p, mf_awthv_p, mf_awthl_p, mf_awqt_p, mf_awql_p, &
                                                    mf_awqi_p, mf_awu_p, mf_awv_p, cfl_mf_p, wthl_sec_ed_p, &
                                                    wthl_sec_mf_p, wqw_sec_ed_p, wqw_sec_mf_p, rho_p, qv_p, qs_p, &
                                                    qi_p, qg_p, qr_p

   real(kind=RKIND),dimension(:,:), intent(inout) :: uw_sfc_p, vw_sfc_p
   real(kind=RKIND),dimension(:,:,:), intent(inout) :: dse_p, qw_p, u_p, v_p, wthv_sec_p, &
                                                       kzh_p, kzm_p, qc_p, cldfrabl_p, theta_p
   real(kind=RKIND), dimension(:,:,:), allocatable:: presdiff_p, thetal_p, qtracers
   real(kind=RKIND), dimension(:,:), allocatable:: wtracer_sfc, dx_p, dy_p,  hfx_converted_p, lh_converted_p

   integer :: num_qtracers, nadv, i, j

   call calculate_thetal(theta_p, t_p, qc_p, thetal_p)
   call calculate_presdiff(pres_p, kte, presdiff_p)
   call calculate_qw(qc_p, qv_p, qs_p, qi_p, qg_p, qr_p, qw_p)
   call calculate_dry_static_energy(t_p, z_p, gravity, dse_p)
   call calculate_surface_momentum_flux(cda_p, u_p(:, 1, :), uw_sfc_p)
   call calculate_surface_momentum_flux(cda_p, v_p(:, 1, :), vw_sfc_p)

   nadv = floor(dt_dyn/shoc_timestep) ! TODO (GTF): Floor or ceiling? It should be integer

   num_qtracers = 0
   allocate(wtracer_sfc(ite, jte))
   wtracer_sfc = 0.0_RKIND
   allocate(qtracers(ite, jte, kte))
   qtracers = 0.0_RKIND

   allocate(dx_p(ite, jte), dy_p(ite, jte))
   ! dx_p and dy_p are calculated as sqrt(area_p) in the main driver
   dx_p = sqrt(area_p)
   dy_p = sqrt(area_p)

   allocate(hfx_converted_p(ite, jte), lh_converted_p(ite, jte))
   do j = jts,jte
   do i = its,ite
      ! Converting hfx_p from W/m^2 to K*m/s
      hfx_converted_p(i,j) = hfx_p(i,j)/(cp*rho_p(i,j,1))
      
      ! Converting lh_p from W/m^2 to kg/kg m/s
      lh_converted_p(i,j) = lh_p(i,j)/(xlv*rho_p(i,j,1))
   enddo
   enddo

   call shoc_main( &
      itimestep, ite, kte, kte+1, shoc_timestep, nadv, &                                    ! Input
	   dx_p(:,1), dy_p(:,1), theta_m(:,:,1),&                                                ! Input
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
      wqls_p(:,:,1),brunt_p(:,:,1),shoc_ql2_p(:,:,1), &                                     ! Output (diagnostic)
      mf_dry_a_p(:,:,1), mf_moist_a_p(:,:,1), &                                             ! Output (EDMF diagnostic)
      mf_dry_w_p(:,:,1), mf_moist_w_p(:,:,1), &                                             ! Output (EDMF diagnostic)
      mf_dry_qt_p(:,:,1), mf_moist_qt_p(:,:,1), &                                           ! Output (EDMF diagnostic)
      mf_dry_thl_p(:,:,1), mf_moist_thl_p(:,:,1), &                                         ! Output (EDMF diagnostic)
      mf_dry_u_p(:,:,1),   mf_moist_u_p(:,:,1), &                                           ! Output (EDMF diagnostic)
      mf_dry_v_p(:,:,1),   mf_moist_v_p(:,:,1), &                                           ! Output (EDMF diagnostic)
                           mf_moist_qc_p(:,:,1), &                                          ! Output (EDMF diagnostic)
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

   call calculate_theta(thetal_p, t_p, qc_p, theta_p)
 end subroutine bl_shocmf_driver

end module module_bl_shocmf
