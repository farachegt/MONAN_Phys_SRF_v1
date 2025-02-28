module module_bl_shocmf
 use mpas_atmphys_constants, only: gravity, R_d, R_v, cp, ep_1, xlv, xlf, karman
 use mpas_atmphys_vars

 use shoc, only: shoc_init, shoc_main, length_fac

 implicit none
 private
 public:: bl_shocmf_init, &
          bl_shocmf_driver

 real(r8) :: shoc_timestep = 150.0  ! Default SHOC timestep set in namelist

 contains

 ! SHOCMF -> MPAS
 ! thetav -> theta_m
 ! ql     -> qc
 ! zt     -> zmid_p
 ! qv     -> qv
 ! rho_i  -> rho_interfaces (criar)
 ! rho    -> rho
 ! zm     -> zgrid
 ! pbl_h  -> hpbl


 subroutine calculate_rho_i(pres_i, t_i, rho_i)
   ! Calculate the density at the interfaces
   ! pres_i and t_i are the pressure and temperature at the interfaces
   ! which in MONAN vars module are called pres2_p and t2_p
   real(kind=RKIND), intent(in) :: pres_i(:), t_i(:)
   real(kind=RKIND), intent(out) :: rho_i(:)
   integer :: k

   ! p = rho * R_d * T
   ! rho = p / (R_d * T)
   do k = 1, kte
      rho_i(k) = pres_i(k) / (R_d * t_i(k))
   enddo
   
 end subroutine calculate_rho_i

 subroutine calculate_theta_l(theta, temperature, ql, theta_l)
   ! Calculate liquid water potential temperature
   ! Reference: https://glossary.ametsoc.org/wiki/Liquid_water_potential_temperature
   real(kind=RKIND), intent(in) :: theta(:), temperature(:), ql(:)
   real(kind=RKIND), intent(out) :: theta_l(:)
   integer :: k

   do k = 1, kte
      theta_l(k) = theta(k) + (xlv/cp)*ql(k)/temperature(k)
   enddo
 end subroutine calculate_theta_l

 subroutine calculate_qw(ql, qv, qs, qi, qg, qr, qw)
   ! Calculate total water mixing ratio
   real(kind=RKIND), intent(in) :: ql(:), qv(:), qs(:), qi(:), qg(:), qr(:)
   real(kind=RKIND), intent(out) :: qw(:)
   integer :: k

   do k = 1, kte
      qw(k) = ql(k) + qv(k) + qs(k) + qi(k) + qg(k) + qr(k)
   enddo
 end subroutine calculate_qw


 subroutine bl_shocmf_init()
    integer:: ntop_shoc, nbot_shoc
    ntop_shoc = 1    ! if >1, must be <= nbot_molec
    nbot_shoc = kte ! currently always pver

    ! TODO: pref_mid == pres_p (?)
    call shoc_init(kte, gravity, R_d, R_v, cp, ep_1, xlv, xlf, karman, pres_p, nbot_shoc, &
                    ntop_shoc, length_fac)
 end subroutine bl_shocmf_init

 subroutine bl_shocmf_driver

   call shoc_main( &
      ite, kte, kte+1, shoc_timestep, nadv, &                                                         ! Input
	   dx_p(:ncol), dy_p(:ncol), thv(:ncol,:),&                                      ! Input
      zt_g(:ncol,:), zi_g(:ncol,:), state%pmid(:ncol,:pver), state%pint(:ncol,:pverp), &        ! Input
      state1%pdel(:ncol,:pver), &                                                               ! Input
	   wpthlp_sfc(:ncol), wprtp_sfc(:ncol), upwp_sfc(:ncol), vpwp_sfc(:ncol), &                  ! Input
	   wtracer_sfc(:ncol,:), edsclr_dim, wm_zt(:ncol,:), &                                       ! Input
	   exner(:ncol,:),state1%phis(:ncol), &                                                      ! Input
	   shoc_s(:ncol,:), tke_zt(:ncol,:), thlm(:ncol,:), rtm(:ncol,:), &                          ! Input/Ouput
	   um(:ncol,:), vm(:ncol,:), edsclr_in(:ncol,:,:), &                                         ! Input/Output
	   wthv(:ncol,:),tkh(:ncol,:),tk(:ncol,:), &                                                 ! Input/Output
	   rcm(:ncol,:),cloud_frac(:ncol,:), &                                                       ! Input/Output
      pblh(:ncol),pblh_wthl(:ncol),pblh_tke(:ncol),&                                            ! Output
      rino(:ncol,:), &                                                                          ! Output
      shoc_mix_out(:ncol,:), isotropy_out(:ncol,:), &                                           ! Output (diagnostic)
      a_diss_out(:ncol,:), a_prod_bu_out(:ncol,:), a_prod_sh_out(:ncol,:), &                    ! Output (diagnostic)
      w_sec_out(:ncol,:), thl_sec_out(:ncol,:), qw_sec_out(:ncol,:), qwthl_sec_out(:ncol,:), &  ! Output (diagnostic)
      wthl_sec_out(:ncol,:), wqw_sec_out(:ncol,:), wtke_sec_out(:ncol,:), &                     ! Output (diagnostic)
      uw_sec_out(:ncol,:), vw_sec_out(:ncol,:), w3_out(:ncol,:), &                              ! Output (diagnostic)
      wqls_out(:ncol,:),brunt_out(:ncol,:),rcm2(:ncol,:), &                                     ! Output (diagnostic)
      mf_dry_a(:ncol,:), mf_moist_a(:ncol,:), &                                                 ! Output (EDMF diagnostic)
      mf_dry_w(:ncol,:), mf_moist_w(:ncol,:), &                                                 ! Output (EDMF diagnostic)
      mf_dry_qt(:ncol,:), mf_moist_qt(:ncol,:), &                                               ! Output (EDMF diagnostic)
      mf_dry_thl(:ncol,:), mf_moist_thl(:ncol,:), &                                             ! Output (EDMF diagnostic)
      mf_dry_u(:ncol,:),   mf_moist_u(:ncol,:), &                                               ! Output (EDMF diagnostic)
      mf_dry_v(:ncol,:),   mf_moist_v(:ncol,:), &                                               ! Output (EDMF diagnostic)
                           mf_moist_qc(:ncol,:), &                                              ! Output (EDMF diagnostic)
      mf_thlflx(:ncol,:), mf_qtflx(:ncol,:), &                                                  ! Output (EDMF diagnostic)
      mf_ae(:ncol,:), mf_aw(:ncol,:), &                                                         ! Output (EDMF diagnostic)
      mf_awthv(:ncol,:), &                                                                      ! Output (EDMF diagnostic)
      mf_awthl(:ncol,:), mf_awqt(:ncol,:), &                                                    ! Output (EDMF diagnostic)
      mf_awql(:ncol,:), mf_awqi(:ncol,:), &                                                     ! Output (EDMF diagnostic)
      mf_awu(:ncol,:), mf_awv(:ncol,:), cfl_mf(:ncol,:), &                                      ! Output (EDMF diagnostic)
      mf_dry_freq(:), mf_moist_freq(:), plumeheight(:), plume_dry_height(:), ztop_B(:), &       ! Output (EDMF diagnostic)
      wthl_sec_ed_out(:ncol,:), wthl_sec_mf_out(:ncol,:), &                                     ! Output (EDMF diagnostic)
      wqw_sec_ed_out(:ncol,:), wqw_sec_mf_out(:ncol,:))                   	

 end subroutine bl_shocmf_driver

end module module_bl_shocmf
