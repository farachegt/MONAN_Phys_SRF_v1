module module_bl_shocmf
 use mpas_atmphys_constants

 use shoc, only: shoc_init, shoc_main, length_fac

 implicit none
 private
 public:: bl_shocmf_init, &
          bl_shocmf_driver

 real(RKIND) :: shoc_timestep = 150.0  ! Default SHOC timestep set in namelist

 contains

 subroutine calculate_dry_static_energy()
! TODO (G): Write the code
 end subroutine calculate_dry_static_energy

 subroutine calculate_presdiff(pres_levels, nlevels, presdiff_levels)
   ! Calculate the pressure difference between the interfaces
   ! pres_levels is the pressure at the levels
   ! presdiff_levels is the pressure difference between the interfaces
   real(kind=RKIND), dimension(:,:,:), intent(in) :: pres_levels
   integer, intent(in) :: nlevels
   real(kind=RKIND), dimension(:,:,:), intent(out) :: presdiff_levels
   integer :: k

   do k=1, nlevels-1
      presdiff_levels(:,k,:) = pres_levels(:,k+1,:) - pres_levels(:,k,:)
   enddo

 end subroutine calculate_presdiff

 subroutine calculate_thetal(theta, temperature, ql, thetal)
   ! Calculate liquid water potential temperature
   ! Reference: https://glossary.ametsoc.org/wiki/Liquid_water_potential_temperature
   real(kind=RKIND), dimension(:,:,:), intent(in) :: theta, temperature, ql
   real(kind=RKIND), dimension(:,:,:), intent(out) :: thetal

   thetal(:,:,:) = theta(:,:,:)*(1.0_RKIND - (xlv/cp)*ql(:,:,:)/temperature(:,:,:))
 end subroutine calculate_thetal

 subroutine calculate_qw(ql, qv, qs, qi, qg, qr, qw)
   ! Calculate total water mixing ratio
   real(kind=RKIND), dimension(:,:,:), intent(in) :: ql, qv, qs, qi, qg, qr
   real(kind=RKIND), dimension(:,:,:), intent(out) :: qw

   qw(:,:,:) = ql(:,:,:) + qv(:,:,:) + qs(:,:,:) + qi(:,:,:) + qg(:,:,:) + qr(:,:,:)
 end subroutine calculate_qw


 subroutine bl_shocmf_init(kte, gravity, R_d, R_v, cp, ep_1, xlv, &
                           xlf, karman, pres_p)
    integer, intent(in) :: kte
    real(RKIND), intent(in) :: gravity, R_d, R_v, cp, ep_1, xlv, xlf, &
                               karman
    real(RKIND), dimension(:,:), intent(in) :: pres_p
  
    integer:: ntop_shoc, nbot_shoc
    ntop_shoc = 1    ! if >1, must be <= nbot_molec
    nbot_shoc = kte ! currently always pver

    call shoc_init(kte, gravity, R_d, R_v, cp, ep_1, xlv, xlf, karman, pres_p, nbot_shoc, &
                    ntop_shoc, length_fac)
 end subroutine bl_shocmf_init

 subroutine bl_shocmf_driver(nlevels)
   ! Driver for the SHOC-MF module
   ! This subroutine is called from the main driver
   ! It calls the SHOC-MF module to calculate the turbulence
   ! and the EDMF module to calculate the subgrid-scale clouds

   real(kind=RKIND), dimension(:,:,:), pointer:: theta_p, t_p, qw_p, qc_p, &
                                                  qv_p, qs_p, qi_p, qg_p, qr_p, &
                                                  thetal_p, pres_p, presdiff_p
   integer, intent(in) :: nlevels
   call calculate_thetal(theta_p, t_p, qc_p, thetal_p)
   call calculate_presdiff(pres_p, nlevels, presdiff_p)
   call calculate_qw(qc_p, qv_p, qs_p, qi_p, qg_p, qr_p, qw_p)

!    call shoc_main( &
!       ite, kte, kte+1, shoc_timestep, nadv, &                                                         ! Input
! 	   dx_p(:ncol), dy_p(:ncol), thv(:ncol,:),&                                      ! Input
!       zt_g(:ncol,:), zi_g(:ncol,:), state%pmid(:ncol,:pver), state%pint(:ncol,:pverp), &        ! Input
!       state1%pdel(:ncol,:pver), &                                                               ! Input
! 	   wpthlp_sfc(:ncol), wprtp_sfc(:ncol), upwp_sfc(:ncol), vpwp_sfc(:ncol), &                  ! Input
! 	   wtracer_sfc(:ncol,:), edsclr_dim, wm_zt(:ncol,:), &                                       ! Input
! 	   exner(:ncol,:),state1%phis(:ncol), &                                                      ! Input
! 	   shoc_s(:ncol,:), tke_zt(:ncol,:), thlm(:ncol,:), rtm(:ncol,:), &                          ! Input/Ouput
! 	   um(:ncol,:), vm(:ncol,:), edsclr_in(:ncol,:,:), &                                         ! Input/Output
! 	   wthv(:ncol,:),tkh(:ncol,:),tk(:ncol,:), &                                                 ! Input/Output
! 	   rcm(:ncol,:),cloud_frac(:ncol,:), &                                                       ! Input/Output
!       pblh(:ncol),pblh_wthl(:ncol),pblh_tke(:ncol),&                                            ! Output
!       rino(:ncol,:), &                                                                          ! Output
!       shoc_mix_out(:ncol,:), isotropy_out(:ncol,:), &                                           ! Output (diagnostic)
!       a_diss_out(:ncol,:), a_prod_bu_out(:ncol,:), a_prod_sh_out(:ncol,:), &                    ! Output (diagnostic)
!       w_sec_out(:ncol,:), thl_sec_out(:ncol,:), qw_sec_out(:ncol,:), qwthl_sec_out(:ncol,:), &  ! Output (diagnostic)
!       wthl_sec_out(:ncol,:), wqw_sec_out(:ncol,:), wtke_sec_out(:ncol,:), &                     ! Output (diagnostic)
!       uw_sec_out(:ncol,:), vw_sec_out(:ncol,:), w3_out(:ncol,:), &                              ! Output (diagnostic)
!       wqls_out(:ncol,:),brunt_out(:ncol,:),rcm2(:ncol,:), &                                     ! Output (diagnostic)
!       mf_dry_a(:ncol,:), mf_moist_a(:ncol,:), &                                                 ! Output (EDMF diagnostic)
!       mf_dry_w(:ncol,:), mf_moist_w(:ncol,:), &                                                 ! Output (EDMF diagnostic)
!       mf_dry_qt(:ncol,:), mf_moist_qt(:ncol,:), &                                               ! Output (EDMF diagnostic)
!       mf_dry_thl(:ncol,:), mf_moist_thl(:ncol,:), &                                             ! Output (EDMF diagnostic)
!       mf_dry_u(:ncol,:),   mf_moist_u(:ncol,:), &                                               ! Output (EDMF diagnostic)
!       mf_dry_v(:ncol,:),   mf_moist_v(:ncol,:), &                                               ! Output (EDMF diagnostic)
!                            mf_moist_qc(:ncol,:), &                                              ! Output (EDMF diagnostic)
!       mf_thlflx(:ncol,:), mf_qtflx(:ncol,:), &                                                  ! Output (EDMF diagnostic)
!       mf_ae(:ncol,:), mf_aw(:ncol,:), &                                                         ! Output (EDMF diagnostic)
!       mf_awthv(:ncol,:), &                                                                      ! Output (EDMF diagnostic)
!       mf_awthl(:ncol,:), mf_awqt(:ncol,:), &                                                    ! Output (EDMF diagnostic)
!       mf_awql(:ncol,:), mf_awqi(:ncol,:), &                                                     ! Output (EDMF diagnostic)
!       mf_awu(:ncol,:), mf_awv(:ncol,:), cfl_mf(:ncol,:), &                                      ! Output (EDMF diagnostic)
!       mf_dry_freq(:), mf_moist_freq(:), plumeheight(:), plume_dry_height(:), ztop_B(:), &       ! Output (EDMF diagnostic)
!       wthl_sec_ed_out(:ncol,:), wthl_sec_mf_out(:ncol,:), &                                     ! Output (EDMF diagnostic)
!       wqw_sec_ed_out(:ncol,:), wqw_sec_mf_out(:ncol,:))                   	

 end subroutine bl_shocmf_driver

end module module_bl_shocmf
