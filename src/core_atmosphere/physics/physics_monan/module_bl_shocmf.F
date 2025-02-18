module module_bl_shocmf
 use mpas_atmphys_constants, only: gravity, R_d, R_v, cp, ep_1, xlv, xlf, karman
 use mpas_atmphys_vars

 use shoc, only: shoc_init, shoc_main, length_fac
 implicit none
 private
 public:: shocmf_bl_init, &
          shocmf_bl_driver

 contains

 subroutine shocmf_bl_init()
    integer:: ntop_shoc, nbot_shoc
    ntop_shoc = 1    ! if >1, must be <= nbot_molec
    nbot_shoc = kte ! currently always pver

    call shoc_init(kte, gravity, R_d, R_v, cp, ep_1, xlv, xlf, karman, pres_p, nbot_shoc, &
                    ntop_shoc, length_fac)
 end subroutine shocmf_bl_init

 subroutine shocmf_bl_driver

 end subroutine shocmf_bl_driver

end module module_bl_shocmf
