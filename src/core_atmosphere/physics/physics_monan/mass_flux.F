module edmf

  use mpas_kind_types, only: rtype=>RKIND, rtype8=>R8KIND

  use mpas_atmphys_constants, only: lcond => xlv, & ! Latent heat of vaporization (J/kg) 
                                    lice => xlf, & ! Latent heat of fusion (J/kg) 
                                    eps => ep_1 ! (rh2o/rair) - 1 
                                    
  use mpas_constants, only: ggr => gravity, & ! gravitational acceleration (m/s**2)
                            rgas, &           ! Dry air gas constant ~ J/K/kg
                            cp                ! specific heat of dry air (J/K/kg)
!  use spmd_utils,    only: masterproc ! TODO (GTF): discover where to get master process in MONAN project

  implicit none
  
  private
  
  public :: integrate_mf, calc_mf_vertflux, compute_tmpi3
  public :: mf_nup, mf_L0, mf_ent0, do_edmf, do_condensation, do_mf_diag, do_wthv_mf
!  public :: mf_readnl, mf_nup, mf_L0, mf_ent0, do_edmf, do_condensation, do_mf_diag, do_wthv_mf
  
  real(rtype) :: mf_L0   = 50.0_rtype    ! Default in namelist_defaults_cam.xml: 50 m
  real(rtype) :: mf_ent0 = 0.22_rtype    ! Default in namelist_defaults_cam.xml: 0.2
  integer     :: mf_nup  = 10           ! Default in namelist_defaults_cam.xml: 10 plumes
  real(rtype) :: mf_a = 1._rtype       ! Default in namelist_defaults_cam.xml: 1
  real(rtype) :: mf_b = 0.5_rtype       ! Default in namelist_defaults_cam.xml: 0.5
  real(rtype) :: mf_c = 0.5_rtype       ! Default in namelist_defaults_cam.xml: 0.5

  logical, protected :: do_edmf         = .false.
  logical, protected :: do_condensation = .true.
  logical, protected :: do_mf_diag      = .false.
  logical, protected :: do_wthv_mf      = .true.


contains

  ! TODO: Consider an "edmf_init" routine to initialize physical constants

  ! =============================================================================== !
  !  Eddy-diffusivity mass-flux routine                                                                               !
  ! =============================================================================== !
!   subroutine mf_readnl(nlfile) ! TODO: Reimplement this method to read a MONAN namelist
!   ! =============================================================================== !
!   ! MF namelists                                                                    !
!   ! =============================================================================== !
!     ! use namelist_utils,  only: find_group_name ! TODO (GTF): Implement find_group_name in MONAN project if necessary to read namelists.
!     ! use cam_abortutils,  only: endrun
!     ! TODO (GTF): Get similar to mpiint, mpii8, ... in MONAN project
!     ! use mpishorthand  ! in spmd_utils.F90: mpiint, mpii8, mpichar, mpilog, mpipk, mpic16, mpir8, mpir4, mpicom, mpimax 
    
!     character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

!     integer :: iunit, read_status

!     namelist /shoc_mf_nl/ mf_L0, mf_ent0, mf_nup, mf_a, mf_b, mf_c, do_edmf, do_condensation, do_mf_diag, do_wthv_mf

!     !  Read namelist to determine if SHOC history should be called
!     if (masterproc) then ! TODO (GTF): voltar
!       iunit = getunit()
!       open( newunit=iunit, file=trim(nlfile), status='old' )

!       call find_group_name(iunit, 'shoc_mf_nl', status=read_status)
!       if (read_status == 0) then
!          read(iunit, nml=shoc_mf_nl, iostat=read_status)
!          if (read_status /= 0) then
!           !  call endrun('mf_readnl:  error reading namelist')
!             stop "mf_readnl:  error reading namelist"
!          end if
!       end if
!       close(iunit) 
!       !close(unit=iunit)
!       !call freeunit(iunit)
!     end if

! #ifdef SPMD
! ! Broadcast namelist variables
!       call mpibcast(mf_L0,            1,    mpir8,   0, mpicom)
!       call mpibcast(mf_ent0,          1,    mpir8,   0, mpicom)
!       call mpibcast(mf_nup,           1,   mpiint,   0, mpicom)
!       call mpibcast(mf_a,             1,    mpir8,   0, mpicom)
!       call mpibcast(mf_b,             1,    mpir8,   0, mpicom)
!       call mpibcast(mf_c,             1,    mpir8,   0, mpicom)
!       call mpibcast(do_edmf,          1,   mpilog,   0, mpicom)
!       call mpibcast(do_condensation,  1,   mpilog,   0, mpicom)  
!       call mpibcast(do_mf_diag,       1,   mpilog,   0, mpicom)  
!       call mpibcast(do_wthv_mf,       1,   mpilog,   0, mpicom)  
! #endif ! TODO (GTF): voltar

!     if ((.not. do_edmf) .and. do_mf_diag ) then
!       call endrun('shoc_mf_nl: Error - cannot turn on do_mf_diag without also turning on do_edmf')
!     end if
    
!   end subroutine mf_readnl

  subroutine integrate_mf(itimestep, shcol, nz, nzi, dt,           & ! input
                 rho_zt_in, rho_zi_in,                             & ! input
                 zt_in, zi_in, dz_zt_in, p_in, thv_zi_in,          & ! input
                 u_in,   v_in,   thl_in,   thv_in, qt_in,          & ! input
                 ust,    wthl,   wqt,   qc_in,                     & ! input
                 pblh, pblh_wthl, wthv_sec_in,                     & ! input
                 dry_a_out,   moist_a_out,                         & ! output: updraft properties for diagnostics
                 dry_w_out,   moist_w_out,                         & ! output: updraft properties for diagnostics
                 dry_qt_out,  moist_qt_out,                        & ! output: updraft properties for diagnostics
                 dry_thl_out, moist_thl_out,                       & ! output: updraft properties for diagnostics
                 dry_u_out,   moist_u_out,                         & ! output: updraft properties for diagnostics
                 dry_v_out,   moist_v_out,                         & ! output: updraft properties for diagnostics
                              moist_qc_out,                        & ! output: updraft properties for diagnostics
                 ae_out, aw_out,                                   & ! output: variables needed for  diffusion solver
                 awthv_out,                                        & ! output: variable needed for total wthv
                 awthl_out, awqt_out,                              & ! output: variables needed for  diffusion solver
                 awql_out, awqi_out,                               & ! output: variables needed for  diffusion solver
                 awu_out, awv_out,                                 & ! output: variables needed for  diffusion solver
                 freq_dry, freq_moist, plumeheight,                & ! output
                 plume_dry_height, ztop_B, cfl)                      ! output: frequency of plume activation (2D)
  
  ! ================================================================================= !
  ! Original author: Marcin Kurowski, JPL
  ! Modified heavily by Maria Chinita and Mikael Witte, UCLA/JPL/NPS for implementation in E3SM/SCREAM
  ! email: maria.j.chinita.candeias@jpl.nasa.gov
  !
  ! Variables needed for solver:
  ! ae = sum_i (1-a_i)
  ! aw = sum (a_i w_i)
  ! awthl = sum(a_i w_i*thl_i)
  ! awqt  = sum(a_i w_i*qt_i)
  ! awql,awqi,awu,awv similar to above except for different variables - not currently coupled to SHOC diffusion solver
  !
  !
  ! - mass flux variables are computed on edges (i.e. momentum grid):
  !  upa,upw,upqt,... 1:nzi
  !  dry_a,moist_a,dry_w,moist_w, ... 1:nzi
  ! ================================================================================= !
  
     ! ============================================================================== ! 
     ! INPUTS   
     ! physics controls
     integer, intent(in) :: shcol,nz,nzi
     real(rtype), dimension(shcol,nz),  intent(in) :: zt_in,   dz_zt_in, rho_zt_in

     real(rtype), dimension(shcol,nzi), intent(in) :: zi_in, p_in, thv_zi_in, rho_zi_in
     real(rtype), dimension(shcol,nz),  intent(in) :: u_in,v_in,thl_in,qt_in,qc_in,thv_in, wthv_sec_in  ! all on thermodynamic/midpoint levels

     real(rtype), dimension(shcol), intent(in) :: ust,   wthl,   wqt
     real(rtype), dimension(shcol), intent(in) :: pblh, pblh_wthl
     !time step [s]   
     real(rtype), intent(in) :: dt
     integer, intent(in) :: itimestep
     ! ============================================================================== !
     ! OUTPUTS
     ! updraft properties
     real(rtype),dimension(shcol,nzi), intent(out) :: dry_a_out,    moist_a_out,     &
                                                      dry_w_out,    moist_w_out,     &
                                                      dry_qt_out,   moist_qt_out,    &
                                                      dry_thl_out,  moist_thl_out,   &
                                                      dry_u_out,    moist_u_out,     &
                                                      dry_v_out,    moist_v_out,     &
                                                      moist_qc_out
     ! variables needed for diffusion solver
     real(rtype),dimension(shcol,nzi), intent(out) :: ae_out,       aw_out,          &
                                                      awthv_out,    awthl_out,       &
                                                      awqt_out,     awql_out,        &
                                                      awqi_out,     awu_out,         &
                                                      awv_out,      cfl

     ! plume activation frequency
     real(rtype),dimension(shcol),     intent(out) :: freq_dry,     freq_moist
     
     ! plume height from one plume test
     real(rtype), dimension(shcol), intent(out) :: plumeheight, plume_dry_height, ztop_B

     
     ! ============================================================================== !
     ! INTERNAL VARIABLES

     ! flipped variables (i.e. here index 1 is at surface)
     real(rtype), dimension(shcol,nz)  :: zt, dz_zt, rho_zt
     real(rtype), dimension(shcol,nzi) :: zi, p, thv_zi, rho_zi
     real(rtype), dimension(shcol,nz)  :: u, v, thl, qt, qc, thv, wthv_sec
     
     ! flipped updraft properties (i.e. index 1 is at surface)
     real(rtype), dimension(shcol,nzi) :: dry_a,     moist_a,      &
                                          dry_w,     moist_w,      &
                                          dry_qt,    moist_qt,     &
                                          dry_thl,   moist_thl,    &
                                          dry_u,     moist_u,      &
                                          dry_v,     moist_v,      &
                                                     moist_qc,     &
                                          dry_th,    moist_th,     &           
                                          ae,        aw,           &
                                          awu,       awv,          &
                                          awthv,     awthl,        &
                                          awqt,      awql,         &
                                          awqv,      awth,         & 
                                          awqi,      awqc                    
                                                          

     ! updraft properties
     real(rtype), dimension(nzi,mf_nup) :: upw,      upa,      &
                                           upthl,    upthv,    &
                                                     upth,     &
                                           upqt,     upqc,     &
                                           upql,     upqv,     &
                                           upqi,     ups,      &
                                           upu,      upv
                                           

     ! entrainment variables
     real(rtype), dimension(nz,mf_nup) :: entf, ent 
     integer,     dimension(nz,mf_nup) :: enti
     real(rtype), dimension(nz) :: ent_oneplume
     
     
     ! other variables
     integer     :: k,j,i,index_top
     real        :: cfl_zt
     
     
     real(rtype) :: wthv,       wstar,     qstar,   thstar,    &
                    sigmaw,   sigmaqt,   sigmath,       z0,    &
                    wmin,        wmax,       wlv,      wtv,    &
                    wp, wstar_aux
                    
                    
     real(rtype) :: pbj,            b,       qtn,     thln,    &
                    thvn,         thn,       qcn,      qln,    &
                    qin,           un,        vn,      wn2,    &
                    entexp,   entexpu,      entw,      iexh

     
     real(rtype) :: plumeheight_aux 
     real(rtype), dimension(shcol) :: plume_top_height
     ! internal surface cont
     real(rtype) :: dzt(nz) 
     ! Dynamic L0 and ztop
     real(rtype) :: dynamic_L0, ztop

     ! w parameters
     ! virtual mass coefficients for w-eqn after Suselj etal 2019
     real(rtype),parameter :: wa = 1._rtype,     &
                              wb = 1.5_rtype


     ! parameters defining initial conditions for updrafts
     real(rtype),parameter :: pwmin = 1.5_rtype, &
                              pwmax = 3._rtype

     ! min values to avoid singularities
     real(rtype),parameter :: wstarmin = 1.e-3_rtype,  &
                              pblhmin  = 100._rtype
                              
    ! fixed entrainment rate (to debug only)
    ! TODO (GTF): Check if this parameter is necessary. It says "to debug only" but variable
    ! TODO (GTF): ent_oneplume is assigned with this value.
     real(rtype),parameter :: fixent = 1.e-3_rtype
     
     logical ::check
     
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !!!!!!!!!!!!!!!!!!!!!! BEGIN CODE !!!!!!!!!!!!!!!!!!!!!!!
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     ! Flip vertical coordinates and all input variables
     do k=1,nz
       ! thermodynamic grid variables
       zt(:,k)    =  zt_in(:,nz-k+1)
       dz_zt(:,k) =  dz_zt_in(:,nz-k+1)

       u(:,k)     =  u_in(:,nz-k+1)
       v(:,k)     =  v_in(:,nz-k+1)
       
       thl(:,k)   =  thl_in(:,nz-k+1)
       thv(:,k)   =  thv_in(:,nz-k+1)
       thv_zi(:,k)   =  thv_zi_in(:,nz-k+1)
       
       qt(:,k)    =  qt_in(:,nz-k+1)
       qc(:,k)    =  qc_in(:,nz-k+1)

       wthv_sec(:,k)  =  wthv_sec_in(:,nz-k+1)

       rho_zt(:,k) =  rho_zt_in(:,nz-k+1)
       rho_zi(:,k) =  rho_zi_in(:,nz-k+1)

       ! momentum altitude grid
       zi(:,k)    = zi_in(:,nzi-k+1)
       p(:,k)     =  p_in(:,nzi-k+1)
     enddo
     zi(:,nzi) = zi_in(:,1)
     p(:,nzi)  =  p_in(:,1)
     
     rho_zi(:,nzi) = rho_zi_in(:,1)

     ! INITIALIZE OUTPUT VARIABLES
     ! set updraft properties to zero
     dry_a     = 0._rtype
     moist_a   = 0._rtype
     dry_w     = 0._rtype
     moist_w   = 0._rtype
     dry_qt    = 0._rtype
     moist_qt  = 0._rtype
     dry_thl   = 0._rtype
     moist_thl = 0._rtype
     dry_u     = 0._rtype
     moist_u   = 0._rtype
     dry_v     = 0._rtype
     moist_v   = 0._rtype
     moist_qc  = 0._rtype
     moist_th  = 0._rtype
     dry_th    = 0._rtype
     ! outputs - variables needed for solver
     aw        = 0._rtype
     awthv     = 0._rtype
     awthl     = 0._rtype
     awqt      = 0._rtype
     awqc      = 0._rtype
     awqv      = 0._rtype
     awql      = 0._rtype
     awqi      = 0._rtype
     awu       = 0._rtype
     awv       = 0._rtype
     awth      = 0._rtype

     ! this is the environmental area - by default 1.
     ae = 1._rtype
     ! CFL number
     cfl = 0._rtype

     ! START MAIN COMPUTATION
     ! NOTE: SHOC does not invert the vertical coordinate, which by default is ordered from lowest to highest pressure
     ! (i.e. top of atmosphere to bottom) so surface-based do loops are performed in reverse (i.e. from nz to 1)
     do j=1,shcol
       ! zero out plume properties
       upw   = 0._rtype
       upthl = 0._rtype
       upthv = 0._rtype
       upqt  = 0._rtype
       upa   = 0._rtype
       upu   = 0._rtype
       upv   = 0._rtype
       upqc  = 0._rtype
       ent   = 0._rtype
       upth  = 0._rtype
       upql  = 0._rtype
       upqi  = 0._rtype
       upqv  = 0._rtype

       ! pblh variable is very off from the correct height according to eye estimation of the profiles
       ! let's calculate pblh here as
       !the level of the minimum of qt gradient. Loop starts at 3 so we don't pick near
       !surface values
              
       ! Using our PBL height:
       pbj = max(pblh(j),pblhmin)
       wthv = wthl(j)+eps*thv(j,1)*wqt(j)
       
       ! if surface buoyancy is positive then do mass-flux, otherwise not
       if (wthv>0.0) then
         dzt = dz_zt(j,:)
                                 
         ! calculate mixed layer height using wthv_sec profile (its minimum/inflection point)
         ztop_B(j) = 0._rtype   
         
         ! TODO (GTF): Changed by me, it was "do k=1,nz-1" because in case k=1 it violates memory access
         do k=2,nz-1
            if (zt(j,k) < 100._rtype) cycle
            if (wthv_sec(j,k) < 0 .and. wthv_sec(j,k-1)>wthv_sec(j,k) .and. wthv_sec(j,k+1)>wthv_sec(j,k) ) then
               ztop_B(j) = zt(j,k)
               exit
            endif
         enddo         
         
         if (plume_dry_height(j) > 0._rtype .and. plume_dry_height(j) > pblhmin) then
            pbj = plume_dry_height(j)
         else
            ! For the first time step we don't have a plume_dry_height value
            pbj = max(ztop_B(j),pblhmin)
         endif
         
         wstar  = max( wstarmin, (ggr/thv(j,1)*wthv*pbj)**(1._rtype/3._rtype) )
         qstar  = wqt(j) / wstar
         thstar = wthl(j) / wstar

         sigmaw  = 0.572_rtype * wstar    
         sigmaqt = 2.89_rtype * abs(qstar) 
         sigmath = 2.89_rtype * abs(thstar)
         
         wmin = sigmaw * pwmin
         wmax = sigmaw * pwmax
              
         ent_oneplume = fixent    
         plumeheight_aux = 0._rtype        
         call oneplume( nz, nzi, zi(j,:), dzt, ent_oneplume, p(j,:), qt(j,:), thl(j,:), thv(j,:),  &
                       thv_zi(j,:), wmax, wmin, sigmaw, sigmaqt, sigmath, wa, wb, &
                       do_condensation, plumeheight_aux)
         
         plumeheight(j) = plumeheight_aux      
         
         ! compute entrainment coefficient
         ! get dz/L0
         ztop = max(pblh_wthl(j),ztop_B(j),pblhmin)
         dynamic_L0 = mf_a*(ztop**mf_b)
         
         do i=1,mf_nup
           do k=1,nz
             !entf(k,i) = dzt(k) / mf_L0
             entf(k,i) = dzt(k) / dynamic_L0
           enddo
         enddo
         
         call Poisson( nz, mf_nup, entf, enti, 69._rtype, itimestep)
          
         ! entrainment: Ent=Ent0/dz*P(dz/L0)
         do i=1,mf_nup
           do k=1,nz
             ent(k,i) = real( enti(k,i))*mf_ent0/dzt(k)
           enddo
         enddo
    
                                   
         do i=1,mf_nup
           ! wlv = w_min_i
           wlv = wmin + (wmax-wmin) / (real(mf_nup)) * (real(i)-1._rtype)
           ! wtv = w_max_i
           wtv = wmin + (wmax-wmin) / (real(mf_nup)) * real(i)

           ! Surface vertical velocity of updraft i: w_i
           upw(1,i) = 0.5_rtype * (wlv+wtv)
           ! Surface area of updraft i: a_i
           upa(1,i) = 0.5_rtype * erf( wtv/(sqrt(2._rtype)*sigmaw) ) &
                      - 0.5_rtype * erf( wlv/(sqrt(2._rtype)*sigmaw) )

           upu(1, i) = u(j,1)
           upv(1, i) = v(j,1)

           upqt(1,i)  = qt(j,1)  + 0.32_rtype * upw(1,i) * sigmaqt/sigmaw
           upthv(1,i) = thv(j,1) + 0.58_rtype * upw(1,i) * sigmath/sigmaw
           upthl(1,i) = upthv(1,i) / (1._rtype+eps*upqt(1,i))
           
           upqv(1,i)  = upqt(1,i)
           
           if (do_condensation) then
                iexh = (1.e5_rtype / p(j,1))**(rgas/cp)
                call condensation_mf(upqt(1,i), upthl(1,i), p(j,1), iexh, &
                                     thvn, qcn, thn, qln, qin)
                upthv(1,i) = thvn
                upqc(1,i) = qcn
                upql(1,i) = qln
                upqi(1,i) = qin
                upth(1,i) = thn
           else
                upqc(1,i) = 0._rtype
                upql(1,i) = 0._rtype
                upqi(1,i) = 0._rtype
                !upthv(1,i) = upthl(1,i)*(1._rtype+eps*upqt(1,i))
                upth(1,i)  = upthl(1,i)
           end if


         enddo


         ! Integrate updrafts
         do i=1,mf_nup
           do k=2,nzi

             entexp  = exp(-ent(k-1,i)*dzt(k-1))
             entexpu = exp(-ent(k-1,i)*dzt(k-1)/3._rtype)
        
             qtn  = qt(j,k-1) *(1._rtype-entexp ) + upqt (k-1,i)*entexp
             thln = thl(j,k-1)*(1._rtype-entexp ) + upthl(k-1,i)*entexp
             un   = u(j,k-1)  *(1._rtype-entexpu) + upu  (k-1,i)*entexpu
             vn   = v(j,k-1)  *(1._rtype-entexpu) + upv  (k-1,i)*entexpu
             iexh = (1.e5_rtype / p(j,k))**(rgas/cp) 

             ! Condensation within updrafts, input/output at full levels:
             if (do_condensation) then
                call condensation_mf(qtn, thln, p(j,k), iexh, &
                                    thvn, qcn, thn, qln, qin)
             else
                thvn = thln*(1._rtype+eps*qtn)
                thn = thln                     
                qcn = 0._rtype
                qin = 0._rtype
                qln = 0._rtype
                
             end if

             ! To avoid singularities w equation has to be computed diferently if wp==0
             b=ggr*(0.5_rtype*(thvn+upthv(k-1,i))/thv(j,k-1)-1._rtype)
             wp = wb*ent(k-1,i)
             if (wp==0._rtype) then
               wn2 = upw(k-1,i)**2._rtype+2._rtype*wa*b*dzt(k-1)
             else
               entw = exp(-2._rtype*wp*dzt(k-1))
               wn2 = entw*upw(k-1,i)**2._rtype+wa*b/(wb*ent(k-1,i))*(1._rtype-entw)
             end if

             if (wn2>0._rtype) then
               upw(k,i)   = sqrt(wn2)
               upthv(k,i) = thvn
               upthl(k,i) = thln
               upqt(k,i)  = qtn
               upqc(k,i)  = qcn
               upu(k,i)   = un
               upv(k,i)   = vn
               upa(k,i)   = upa(k-1,i)
               upth(k,i)  = thn
               upql(k,i)  = qln
               upqi(k,i)  = qin
               upqv(k,i)  = qtn - qcn
             else
               exit
             end if
           enddo
         enddo

         ! writing updraft properties for output
         ! all variables, except areas (moist_a and dry_a) are now multipled by the area
         do k=1,nzi

           ! first sum over all i-updrafts
           do i=1,mf_nup
             if (upqc(k,i)>0._rtype) then
               moist_a(j,k)   = moist_a(j,k)   + upa(k,i)
               moist_w(j,k)   = moist_w(j,k)   + upa(k,i)*upw(k,i)
               moist_qt(j,k)  = moist_qt(j,k)  + upa(k,i)*upqt(k,i)
               moist_thl(j,k) = moist_thl(j,k) + upa(k,i)*upthl(k,i)
               moist_th(j,k)  = moist_th(j,k)  + upa(k,i)*upth(k,i)
               moist_u(j,k)   = moist_u(j,k)   + upa(k,i)*upu(k,i)
               moist_v(j,k)   = moist_v(j,k)   + upa(k,i)*upv(k,i)
               moist_qc(j,k)  = moist_qc(j,k)  + upa(k,i)*upqc(k,i)
             else
               dry_a(j,k)     = dry_a(j,k)     + upa(k,i)
               dry_w(j,k)     = dry_w(j,k)     + upa(k,i)*upw(k,i)
               dry_qt(j,k)    = dry_qt(j,k)    + upa(k,i)*upqt(k,i)
               dry_thl(j,k)   = dry_thl(j,k)   + upa(k,i)*upthl(k,i)
               dry_th(j,k)    = dry_th(j,k)    + upa(k,i)*upth(k,i)
               dry_u(j,k)     = dry_u(j,k)     + upa(k,i)*upu(k,i)
               dry_v(j,k)     = dry_v(j,k)     + upa(k,i)*upv(k,i)
             endif
           enddo

           if ( dry_a(j,k) > 0._rtype ) then
             dry_w(j,k)   = dry_w(j,k)   / dry_a(j,k)
             dry_qt(j,k)  = dry_qt(j,k)  / dry_a(j,k)
             dry_thl(j,k) = dry_thl(j,k) / dry_a(j,k)
             dry_th(j,k)  = dry_th(j,k)  / dry_a(j,k)
             dry_u(j,k)   = dry_u(j,k)   / dry_a(j,k)
             dry_v(j,k)   = dry_v(j,k)   / dry_a(j,k)
           else
             dry_w(j,k)   = 0._rtype
             dry_qt(j,k)  = 0._rtype
             dry_thl(j,k) = 0._rtype
             dry_th(j,k)  = 0._rtype
             dry_u(j,k)   = 0._rtype
             dry_v(j,k)   = 0._rtype
           endif

           if ( moist_a(j,k) > 0._rtype ) then
             moist_w(j,k)   = moist_w(j,k)   / moist_a(j,k)
             moist_qt(j,k)  = moist_qt(j,k)  / moist_a(j,k)
             moist_thl(j,k) = moist_thl(j,k) / moist_a(j,k)
             moist_th(j,k)  = moist_th(j,k)  / moist_a(j,k)
             moist_u(j,k)   = moist_u(j,k)   / moist_a(j,k)
             moist_v(j,k)   = moist_v(j,k)   / moist_a(j,k)
             moist_qc(j,k)  = moist_qc(j,k)  / moist_a(j,k)
           else
             moist_w(j,k)   = 0._rtype
             moist_qt(j,k)  = 0._rtype
             moist_thl(j,k) = 0._rtype
             moist_th(j,k)  = 0._rtype
             moist_u(j,k)   = 0._rtype
             moist_v(j,k)   = 0._rtype
             moist_qc(j,k)  = 0._rtype
           endif
             
         enddo

         do k=1,nzi
           do i=1,mf_nup
             ae  (j,k) = ae  (j,k) - upa(k,i)
             aw  (j,k) = aw  (j,k) + upa(k,i)*upw(k,i)
             awu (j,k) = awu (j,k) + upa(k,i)*upw(k,i)*upu(k,i)
             awv (j,k) = awv (j,k) + upa(k,i)*upw(k,i)*upv(k,i)
             awthv(j,k)= awthv(j,k)+ upa(k,i)*upw(k,i)*upthv(k,i)
             awthl(j,k)= awthl(j,k)+ upa(k,i)*upw(k,i)*upthl(k,i) !*cpair/iexh
             awth(j,k) = awth(j,k) + upa(k,i)*upw(k,i)*upth(k,i) !*cpair/iexh
             awqt(j,k) = awqt(j,k) + upa(k,i)*upw(k,i)*upqt(k,i)
             awqc(j,k) = awqc(j,k) + upa(k,i)*upw(k,i)*upqc(k,i)             
             awqv(j,k) = awqv(j,k) + upa(k,i)*upw(k,i)*upqv(k,i)
             awql(j,k) = awql(j,k) + upa(k,i)*upw(k,i)*upql(k,i)
             awqi(j,k) = awqi(j,k) + upa(k,i)*upw(k,i)*upqi(k,i)
           enddo
         enddo

         ! Find highest vertical level where the plume ensemble is dry (i.e., moist_qc = 0)
         check  = .true. 
         plume_dry_height(j) = 0._rtype
         index_top = 1
         do k=1,nzi
            if (check .and. aw(j,k) > 0._rtype .and. moist_qc(j,k) .EQ. 0._rtype) then
               plume_dry_height(j) = zi(j,k)
               index_top = k
            else
               check = .false.  
            endif
         enddo
         
         !! Highest vertical level reached by the moist plumes 
         !(analyzed from the top of the dry CBL given by plume_dry_height)
         check  = .true. 
         plume_top_height(j) = 0._rtype
         do k=index_top+1,nzi
            if (check .and. aw(j,k) > 0._rtype .and. moist_qc(j,k) > 0._rtype) then
               plume_top_height(j) = zi(j,k)
            else
               check = .false.  
            endif
         enddo

         ! Check CLF condition on mass-flux (aw)
         do k=1,nz
           cfl_zt = (2._rtype/dt)*rho_zt(j,k)*dz_zt(j,k)
           if (zi(i,k) < ztop*1.5_rtype) then
           
              if (aw(j,k) > (cfl_zt/rho_zi(j,k)) ) then
                 print*,'WARNING: aw > CFL'
                 print*,'aw(j,k) = ',aw(j,k)
                 print*,'CFL = ',cfl_zt/rho_zi(j,k)
                 print*,'k index = ',k
              endif
             
           endif
           cfl(i,k) = cfl_zt/rho_zi(j,k)   
         enddo
         
       end if  ! ( wthv > 0.0 )

       if (ANY(dry_a  (j,:)>0._rtype)) freq_dry(j)   = 1._rtype
       if (ANY(moist_a(j,:)>0._rtype)) freq_moist(j) = 1._rtype
     end do ! j=1,shcol

     ! flip output variables so index 1 = model top (i.e. lowest pressure)
     do k=1,nzi
       dry_a_out(:,nzi-k+1) = dry_a(:,k)
       dry_w_out(:,nzi-k+1) = dry_w(:,k)
       dry_qt_out(:,nzi-k+1) = dry_qt(:,k)
       dry_thl_out(:,nzi-k+1) = dry_thl(:,k)
       dry_u_out(:,nzi-k+1) = dry_u(:,k)
       dry_v_out(:,nzi-k+1) = dry_v(:,k)

       moist_a_out(:,nzi-k+1) = moist_a(:,k)
       moist_w_out(:,nzi-k+1) = moist_w(:,k)
       moist_qt_out(:,nzi-k+1) = moist_qt(:,k)
       moist_thl_out(:,nzi-k+1) = moist_thl(:,k)
       moist_u_out(:,nzi-k+1) = moist_u(:,k)
       moist_v_out(:,nzi-k+1) = moist_v(:,k)
       moist_qc_out(:,nzi-k+1) = moist_qc(:,k)

       ae_out(:,nzi-k+1) = ae(:,k)
       aw_out(:,nzi-k+1) = aw(:,k)
       awthv_out(:,nzi-k+1) = awthv(:,k)
       awthl_out(:,nzi-k+1) = awthl(:,k)
       awqt_out(:,nzi-k+1) = awqt(:,k)
       awql_out(:,nzi-k+1) = awql(:,k)
       awqi_out(:,nzi-k+1) = awqi(:,k)
       awu_out(:,nzi-k+1) = awu(:,k)
       awv_out(:,nzi-k+1) = awv(:,k)


     end do


  end subroutine integrate_mf

                       
  subroutine oneplume( nz, nzi, zi, dzt, ent, p, qt, thl, thv,   &
                       thv_zi, wmax, wmin, sigmaw, sigmaqt, sigmath, wa, wb, &
                       do_condensation, plumeheight )
  !**********************************************************************
  ! Calculate a single plume with zero entrainment
  ! to be used for a dynamic mixing length calculation
  ! By Rachel Storer
  !**********************************************************************

    integer,  intent(in)                    :: nz, nzi

    real(rtype), intent(in)                 :: wmax, wmin, sigmaw, sigmaqt, sigmath, wa, wb

    real(rtype), dimension(nz),  intent(in) ::  dzt, qt, thl, thv, ent
    real(rtype), dimension(nzi), intent(in) ::  zi, p, thv_zi

                                                      
    logical, intent(in)                       :: do_condensation

    real(rtype), intent(inout) :: plumeheight
     
    !local variables
    integer                        :: k
    real(rtype)                    :: thvn, qtn, thln, qcn, thn, qln, qin, wn2
    real(rtype)                    :: iexh, entexp, entexpu, wp, entw
    real(rtype)                    :: ztop
    real(rtype), dimension(nzi)    :: upw, upa, upqt, upthv, upthl, upth, &
                                      upqc, upql, upqi, b, thvflx
                                      
        

    thvflx  = 0._rtype
    b     = 0._rtype
    upw   = 0._rtype
    upthl = 0._rtype
    upthv = 0._rtype
    upqt  = 0._rtype
    upa   = 0._rtype
    upqc  = 0._rtype
    upth  = 0._rtype
    upql  = 0._rtype
    upqi  = 0._rtype
       
    upw(1) = 0.5_rtype * (wmax+wmin)
    upa(1) = 0.5_rtype * erf( wmax/(sqrt(2.5_rtype)*sigmaw) ) &
                      - 0.5_rtype * erf( wmin/(sqrt(2._rtype)*sigmaw) )
  
    upqt(1)  = qt(1)  + 0.32_rtype * upw(1) * sigmaqt/sigmaw
    upthv(1) = thv(1) + 0.58_rtype * upw(1) * sigmath/sigmaw
           
    upthl(1) = upthv(1) / (1._rtype+eps*upqt(1))
    upth(1)  = upthl(1)
  
    ! get cloud, lowest momentum level 
    if (do_condensation) then
      iexh = (1.e5_rtype / p(1))**(rgas/cp)
      call condensation_mf(upqt(1), upthl(1), p(1), iexh, &
                           thvn, qcn, thn, qln, qin)
      upthv(1) = thvn
      upqc(1)  = qcn
      upql(1)  = qln
      upqi(1)  = qin
      upth(1)  = thn
    else
      ! assume no cldliq
      upthv(1) = upthl(1)*(1._rtype+eps*upqt(1))
      upth(1)  = upthl(1)

    end if
  
    do k=2,nzi
   
      entexp  = exp(-ent(k-1)*dzt(k-1))
      entexpu = exp(-ent(k-1)*dzt(k-1)/3._rtype)
              
      ! integrate updraft
      qtn  = qt(k-1) *(1._rtype-entexp ) + upqt (k-1)*entexp
      thln = thl(k-1)*(1._rtype-entexp ) + upthl(k-1)*entexp
                      
      ! get cloud, momentum levels
      if (do_condensation) then
        iexh = (1.e5_rtype / p(k))**(rgas/cp)
        call condensation_mf(qtn, thln, p(k), iexh, &
                             thvn, qcn, thn, qln, qin)
      else
        thvn = thln*(1._rtype+eps*qtn)
        ! TODO (GTF): Check with Chinita this thn variable?
        thn = thln                       ! THIS NEEEDS TO BE FIXED!! CALCULATED CORRECTLY
        qcn = 0._rtype
        qin = 0._rtype
        qln = 0._rtype
      end if
      ! get buoyancy
      b(k)=ggr*(0.5_rtype*(thvn+upthv(k-1))/thv(k-1)-1._rtype)
      wp = wb*ent(k-1)
      if (wp==0._rtype) then
         wn2 = upw(k-1)**2._rtype+2._rtype*wa*b(k)*dzt(k-1)
      else
         entw = exp(-2._rtype*wp*dzt(k-1))
         wn2 = entw*upw(k-1)**2._rtype+wa*b(k)/(wb*ent(k-1))*(1._rtype-entw)
      endif   
  
      if (wn2>0._rtype) then
        upw(k)   = sqrt(wn2)
        upthv(k) = thvn
        upthl(k) = thln
        upqt(k)  = qtn
        upqc(k)  = qcn
        upa(k)   = upa(k-1)
        upql(k)  = qln
        upqi(k)  = qin
        upth(k)  = thn
        plumeheight = zi(k)
      else
        exit
      end if
      
    enddo
     

  end subroutine oneplume
  


  subroutine condensation_mf( qt, thl, p, iex, thv, qc, th, ql, qi)
  !
  ! zero or one condensation for edmf: calculates thv and qc
  !
     !!!!!  use wv_saturation,      only : qsat

       real(rtype),intent(in) :: qt,thl,p,iex
       real(rtype),intent(out):: thv,qc,th,ql,qi

       !local variables
       integer :: niter,i
       real(rtype) :: diff,t,qs,qcold,es,wf

       ! max number of iterations
       niter=50
       ! minimum difference
       diff=2.e-5_rtype

       qc=0._rtype
       t=thl/iex

  !by definition:
  ! T   = Th*Exner, Exner=(p/p0)^(R/cp)   (1)
  ! Thl = Th - L/cp*ql/Exner              (2)
  !so:
  ! Th  = Thl + L/cp*ql/Exner             (3)
  ! T   = Th*Exner=(Thl+L/cp*ql/Exner)*Exner    (4)
  !     = Thl*Exner + L/cp*ql
       do i=1,niter
         wf = get_watf(t)
         t = thl/iex+get_alhl(wf)/cp*qc   !as in (4)

         ! qsat, p is in pascal (check!)
         !--- call qsat(t,p,es,qs)
         qs=satur_spec_hum(t,p)
         
         qcold = qc
         qc = max(0.5_rtype*qc+0.5_rtype*(qt-qs),0._rtype)
         if (abs(qc-qcold)<diff) exit
       enddo

       wf = get_watf(t)
       t = thl/iex+get_alhl(wf)/cp*qc


      !--- call qsat(t,p,es,qs)
       qs=satur_spec_hum(t,p)
      !----
       qc = max(qt-qs,0._rtype)
       thv = (thl+get_alhl(wf)/cp*iex*qc)*(1.+eps*(qt-qc)-qc)
       th = t*iex
       qi = qc*(1.-wf)
       ql = qc*wf

       contains

       function get_watf(t)
         real(rtype) :: t,get_watf,tc
         real(rtype), parameter :: &
         tmax=-10._rtype, &
         tmin=-40._rtype

         tc=t-273.16_rtype

         if (tc>tmax) then
           get_watf=1._rtype
         else if (tc<tmin) then
           get_watf=0._rtype
         else
           get_watf=(tc-tmin)/(tmax-tmin);
         end if

       end function get_watf


       function get_alhl(wf)
       !latent heat of the mixture based on water fraction
         real(rtype) :: get_alhl,wf

         get_alhl = wf*lcond+(1._rtype-wf)*(lcond+lice)

       end function get_alhl

  end subroutine condensation_mf

  subroutine calc_mf_vertflux(shcol,nlev,nlevi,aw,awvar,var,var_zi,varflx)

    implicit none

    ! INPUT VARIABLES
    ! number of SHOC columns
    integer, intent(in) :: shcol
    ! number of midpoint levels
    integer, intent(in) :: nlev
    ! number of interface levels
    integer, intent(in) :: nlevi
    ! Sum plume (a_i*w_i) [m/s]
    real(rtype), intent(in) :: aw(shcol,nlevi)
    ! Sum plume vertical flux of generic variable var (a_i*w_i*var_i) [units vary]
    real(rtype), intent(in) :: awvar(shcol,nlevi)
    ! Input variable on thermo/full grid [units vary]
    real(rtype), intent(in) :: var_zi(shcol,nlevi) ! NOTE: var is interpolated to zi, so has dim nzi
    real(rtype), intent(in) :: var(shcol,nlev)

    ! OUTPUT VARIABLE
    real(rtype), intent(out) :: varflx(shcol,nlevi)

    ! INTERNAL VARIABLES
    integer :: i,k

    ! diagnose MF fluxes
    varflx(:shcol,1) = 0._rtype
    do k=2,nlev
      do i=1,shcol
        varflx(i,k)= awvar(i,k) - aw(i,k)*0.5*(var(i,k-1)+var(i,k)) ! centered differences
        !varflx(i,k)= awvar(i,k) - aw(i,k)*var(i,k)   ! upwind scheme (in reference to the surface)
        !varflx(i,k)= awvar(i,k) - aw(i,k)*var(i,k-1) ! downwind scheme (in reference to the surface)
        !varflx(i,k)= awvar(i,k) - aw(i,k)*var_zi(i,k)

      end do
    end do
    varflx(:shcol,nlevi) = 0._rtype
    
  end subroutine calc_mf_vertflux

  subroutine compute_tmpi3(nlevi, shcol, dtime, rho_zi, tmpi3)

    !intent-ins
    integer,     intent(in) :: nlevi, shcol
    !time step [s]
    real(rtype), intent(in) :: dtime
    !air density at interfaces [kg/m3]
    real(rtype), intent(in) :: rho_zi(shcol,nlevi)

    !intent-out
    real(rtype), intent(out) :: tmpi3(shcol,nlevi)

    !local vars
    integer :: i, k

    tmpi3(:,1) = 0._rtype
    ! eqn: tmpi3 = dt*g*rho
    do k = 2, nlevi
      do i = 1, shcol
         tmpi3(i,k) = dtime *  ggr*rho_zi(i,k)
      enddo
    enddo

  end subroutine compute_tmpi3

    subroutine Poisson(nz,nup,mu,POI,seed,itimestep)
      implicit none
      integer, intent(in) :: nz,nup, itimestep
      real(rtype),dimension(nz,nup),intent(in) :: MU
      integer, dimension(nz,nup), intent(out) :: POI
      real(rtype), intent(in) :: seed
      integer :: seed_len,i,j
      integer,allocatable:: the_seed(:)
      !!
      integer, dimension(nz,nup) :: poi_aux

      if (is_first_step(itimestep)) then
        call random_seed(SIZE=seed_len)
        allocate(the_seed(seed_len))
        the_seed(1) = int((seed - int(seed)) * 1000000000._rtype)
        if (seed_len > 1) the_seed(2:) = the_seed(1)

        call random_seed(put=the_seed)
      end if


      do i=1,nz
        do j=1,nup
           poi(i,j) = int(poidev(mu(i,j)))   
        enddo
      enddo
      
    end subroutine Poisson

    FUNCTION poidev(xm)
      IMPLICIT NONE
      REAL(rtype), INTENT(IN) :: xm
      REAL(rtype) :: poidev
      REAL(rtype), PARAMETER :: PI=3.141592653589793238462643383279502884197_rtype
      !Returns as a floating-point number an integer value that is a random deviate drawn from a
      !Poisson distribution of mean xm, using ran1 as a source of uniform random deviates.
      REAL(rtype) :: em,harvest,t,y
      REAL(rtype), SAVE :: alxm,g,oldm=-1.0_rtype,sq

      !Begin code
      !oldm is a flag for whether xm has changed since last call.
      if (xm < 12.0) then !Use direct method.
        if (xm /= oldm) then
          oldm=xm
          g=exp(-xm) !If xm is new, compute the exponential.
        end if
        em=-1
        t=1.0_rtype
        do
          em=em+1.0_rtype     !Instead of adding exponential deviates it is
                           !equivalent to multiply uniform deviates.
                           !We never actually have to take the log;
                           !merely compare to the pre-computed exponential.
          call random_number(harvest)
          t=t*harvest
          if (t <= g) exit
        end do
      else      !    Use rejection method.
        if (xm /= oldm) then  !If xm has changed since the last call, then precompute
                               !some functions that occur below.
          oldm=xm
          sq=sqrt(2.0_rtype*xm)
          alxm=log(xm)
          g=xm*alxm-gammln_s(xm+1.0_rtype) ! The function gammln is the natural log of the
                                        ! gamma function, as given in §6.1.
        end if
        do
          do
            call random_number(harvest)  !y is a deviate from a Lorentzian comparison
            y=tan(PI*harvest)   !function.
            em=sq*y+xm          !em is y, shifted and scaled.
            if (em >= 0.0) exit !Reject if in regime of zero probability.
          end do
          em=int(em)          ! The trick for integer-valued distributions.
          t=0.9_rtype*(1.0_rtype+y**2)*exp(em*alxm-gammln_s(em+1.0_rtype)-g)
          !The ratio of the desired distribution to the comparison function; we accept or reject
          !by comparing it to another uniform deviate. The factor 0.9 is chosen so that t never
          !exceeds 1.
          call random_number(harvest)
          if (harvest <= t) exit
        end do
      end if
      poidev=em
    END FUNCTION poidev

    FUNCTION arth_d(first,increment,n)
      implicit none
      REAL(rtype8), INTENT(IN) :: first,increment
      INTEGER, PARAMETER :: NPAR_ARTH=16,NPAR2_ARTH=8
      INTEGER, INTENT(IN) :: n
      REAL(rtype8), DIMENSION(n) :: arth_d
      INTEGER :: k,k2
      REAL(rtype8) :: temp

      ! Begin code
      if (n > 0) arth_d(1)=first
      if (n <= NPAR_ARTH) then
        do k=2,n
          arth_d(k)=arth_d(k-1)+increment
        end do
      else
        do k=2,NPAR2_ARTH
          arth_d(k)=arth_d(k-1)+increment
        end do
        temp=increment*NPAR2_ARTH
        k=NPAR2_ARTH
        do
          if (k >= n) exit
          k2=k+k
          arth_d(k+1:min(k2,n))=temp+arth_d(1:min(k,n-k))
          temp=temp+temp
          k=k2
        end do
      end if
    END FUNCTION arth_d

    FUNCTION gammln_s(xx)
      IMPLICIT NONE
      REAL(rtype), INTENT(IN) :: xx
      REAL(rtype) :: gammln_s
      !Returns the value ln[Γ(xx)] for xx > 0.
      REAL(rtype8) :: tmp,x
      !Internal arithmetic will be done in double precision, a nicety that you can omit if five-figure
      !accuracy is good enough.
      REAL(rtype8) :: stp = 2.5066282746310005_rtype8
      REAL(rtype8), DIMENSION(6) :: coef = (/76.18009172947146_rtype8,&
        -86.50532032941677_rtype8,24.01409824083091_rtype8,&
        -1.231739572450155_rtype8,0.1208650973866179e-2_rtype8,&
        -0.5395239384953e-5_rtype8/)

      !Begin code
      !call assert(xx > 0.0, ’gammln_s arg’)
      !if (xx .le. 0.) print *,'gammaln fails'
      x=xx
      tmp=x+5.5_rtype8
      tmp=(x+0.5_rtype8)*log(tmp)-tmp
      gammln_s=tmp+log(stp*(1.000000000190015_rtype8+&
      sum(coef(:)/arth_d(x+1.0_rtype8,1.0_rtype8,size(coef))))/x)
    END FUNCTION gammln_s
    
    !---------------------------------------------------------------------------------------------------
    real function satur_spec_hum(pt,press) result(qs)
        implicit none
        real, parameter :: p_max_qsat = 0.5
        real, intent(in) :: pt,press ! Kelvin, hPa
        !---locals
        real :: zew,zqs,zcor,foealfcu,foeewmcu,press_hpa
        real, parameter ::            &
            c_rd=287.06               &
            ,c_rv=461.52              &
            ,c_rtt=273.16             &
            ,c_retv=c_rv/c_rd-1.0     &
            ,c_r2es=611.21*c_rd/c_rv  &
            ,c_r3les=17.502           &
            ,c_r3ies=22.587           &
            ,c_r4les=32.19            &
            ,c_r4ies=-0.7             &
            ,c_rtwat=c_rtt            &
            ,c_rtice=c_rtt-23.        &
            ,c_rticecu=c_rtt-23.      &
            ,c_rtwat_rtice_r=1./(c_rtwat-c_rtice)   &
            ,c_rtwat_rticecu_r=1./(c_rtwat-c_rticecu)
        
        ! press is in Pa, convert to hPa
        press_hpa = press / 100.0

        foealfcu = min(1.0,((max(c_rticecu,min(c_rtwat,pt))-c_rticecu)*c_rtwat_rticecu_r)**2)
        foeewmcu = c_r2es *(    foealfcu *exp(c_r3les*(pt-c_rtt)/(pt-c_r4les))+ &
                           (1.0-foealfcu)*exp(c_r3ies*(pt-c_rtt)/(pt-c_r4ies)))

        zew  = foeewmcu
        zqs  = zew/(100.*press_hpa)
        if(1.0-c_retv*zqs > 0. )then
            zcor = 1.0/(1.0-c_retv*zqs)  ! divide by zero
            qs= zqs*zcor
        else
            qs= p_max_qsat
        endif
    end function satur_spec_hum

      logical function is_first_step(itimestep)
        ! Return true on first step of initial run only.
        implicit none
        integer, intent(in) :: itimestep

        is_first_step = (itimestep == 1)
    end function is_first_step

end module edmf

! End of EDMF module. Thanks for visiting.
