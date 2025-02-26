module module_bl_shocmf
 use mpas_atmphys_constants, only: gravity, R_d, R_v, cp, ep_1, xlv, xlf, karman
 use mpas_atmphys_vars

 use shoc, only: shoc_init, shoc_main, length_fac

 implicit none
 private
 public:: bl_shocmf_init, &
          bl_shocmf_driver

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

 end subroutine bl_shocmf_driver

end module module_bl_shocmf
