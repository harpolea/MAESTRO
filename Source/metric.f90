module metric_module
    use bl_types           , only: dp_t
    use bl_constants_module, only: ZERO, HALF, TWO
    use multifab_module !   , only: multifab, multifab_build, multifab_build_edge, &
                        !           destroy, setval, nghost, &
                        !           extent, multifab_volume, nboxes, &
                        !           multifab_copy, multifab_copy_c, &
                        !           multifab_sub_sub, multifab_div_div_s, multifab_plus_plus
    use ml_layout_module   , only: ml_layout
    use define_bc_module   , only: bc_tower
    use parallel           , only: parallel_IOProcessor, parallel_IOProcessorNode, &
                                   parallel_wtime, parallel_reduce, parallel_barrier, &
                                   MPI_MAX
    use particle_module

    implicit none

    private :: root_find_on_me, Newton_Raphson

    public :: make_weak_field, inverse_metric, calcW, calcu0, &
              g, christoffels, cons_to_prim

contains

    subroutine make_weak_field(nr, rs, alpha, beta, gam, mla)
        ! NOTE: gam is a vector as will assume it is a diagonal matrix
        ! Returns alpha, beta and gamma of weak field metric
        ! FIXME: I have no idea how to implement this in a useful way.

        use probin_module, only : g, Rr, c

        type(multifab) , intent(inout) :: alpha(:)
        type(multifab) , intent(inout) :: beta(:)
        type(multifab) , intent(inout) :: gam(:)
        integer, intent(in)    :: nr(:)
        real(kind=dp_t)   , intent(in) :: rs(:,:)
        type(ml_layout), intent(inout) :: mla

        real(kind=dp_t), pointer :: alphap(:,:,:,:)
        real(kind=dp_t), pointer :: betap(:,:,:,:)
        real(kind=dp_t), pointer :: gamp(:,:,:,:)

        integer :: n, i, j, dm, nlevs

        nlevs = mla%nlevel

        do n = 1, nlevs
            do j = 1, nfabs(alpha(n))

                alphap => dataptr(alpha(n),j)
                betap => dataptr(beta(n),j)
                gamp => dataptr(gam(n),j)

                betap = 0.d0

                do i = 0, nr(n)
                    alphap(:,:,i,1) = sqrt(1.d0 - 2.d0 * (1.d0 - rs(n,i) / Rr) * g / c**2)
                    gamp(:,:,i,:) = 1.d0 + 2.d0 * (1.d0 - rs(n,i) / Rr) * g / c**2
                enddo

            enddo
        enddo


    end subroutine make_weak_field

    subroutine inverse_metric(alpha, beta, gam, g_up)
        ! Return components of the inverse metric. Works for any 3+1 metric.

        real(kind=dp_t)   , intent(inout) :: alpha(0:,0:,0:)
        real(kind=dp_t)   , intent(inout) :: beta(:,0:,0:,0:)
        real(kind=dp_t)   , intent(inout) :: gam(:,0:,0:,0:)

        real(kind=dp_t)   , intent(out) :: g_up(0:,0:,0:,0:,0:)

        integer             :: i, j
        real(kind=dp_t)     :: eye(1:3,1:3)

        ! make identity matrix
        eye(:,:) = 0.d0
        do i = 1,3
            eye(i,i) = 1.d0
        enddo

        g_up(0,0,:,:,:) = -1 / alpha(:,:,:)**2

        do i = 1, 3
            g_up(0,i,:,:,:) = beta(i,:,:,:) / alpha(:,:,:)**2
            g_up(i,0,:,:,:) = beta(i,:,:,:) / alpha(:,:,:)**2

            do j = 1, 3
                g_up(i,j,:,:,:) = eye(i,j) * gam(i,:,:,:) - beta(i,:,:,:) * beta(j,:,:,:) / alpha(:,:,:)**2
            enddo
        enddo

    end subroutine inverse_metric

    subroutine calcW(alpha, beta, gam, u, mla, W_lor)
        ! Calculates the Lorentz factor
        use probin_module, only : c

        type(multifab)    , intent(in)    :: alpha(:)
        type(multifab)    , intent(in)    :: beta(:)
        type(multifab)    , intent(in)    :: gam(:)
        type(multifab)    , intent(in)    :: u(:)
        type(ml_layout)   , intent(inout) :: mla

        type(multifab)    , intent(in)    :: W_lor(:)

        real(kind=dp_t), pointer:: alphap(:,:,:,:)
        real(kind=dp_t), pointer:: betap(:,:,:,:)
        real(kind=dp_t), pointer:: gamp(:,:,:,:)
        real(kind=dp_t), pointer:: up(:,:,:,:)
        real(kind=dp_t), pointer:: Wp(:,:,:,:)
        integer     :: dm, nlevs, n, i, j, k
        real(kind=dp_t)     :: eye(1:3,1:3)

        ! make identity matrix
        eye(:,:) = 0.d0
        do i = 1,3
            eye(i,i) = 1.d0
        enddo

        do n = 1, nlevs
            do k = 1, nfabs(alpha(n))
                alphap => dataptr(alpha(n),k)
                betap => dataptr(beta(n),k)
                gamp => dataptr(gam(n),k)
                up => dataptr(u(n),k)
                Wp => dataptr(W_lor(n),k)

                Wp(:,:,:,:) = 0.d0

                do i = 1, 3
                    do j = 1, 3
                        Wp(:,:,:,1) = Wp(:,:,:,1) + eye(i,j) * &
                        gamp(i,:,:,:) * (up(:,:,:,i) + betap(:,:,:,i)) * &
                        (up(:,:,:,j) + betap(:,:,:,j)) / alphap(:,:,:,1)**2
                    enddo
                enddo
            enddo
        enddo

        Wp(:,:,:,:) = 1 / sqrt(1 - Wp(:,:,:,:)/c**2)

    end subroutine calcW

    subroutine calcu0(alpha, beta, gam, W_lor, u0, mla)
        ! Calculates timelike coordinate of 3+1 velocity
        ! using Lorentz factor and alpha:
        ! W = alpha * u0

        type(multifab)    , intent(in)    :: alpha(:)
        type(multifab)    , intent(in)    :: beta(:)
        type(multifab)    , intent(in)    :: gam(:)
        type(multifab)    , intent(in)    :: W_lor(:)
        type(multifab)   , intent(inout) :: u0(:)
        type(ml_layout)   , intent(inout) :: mla

        real(kind=dp_t), pointer :: alphap(:,:,:,:)
        real(kind=dp_t), pointer :: betap(:,:,:,:)
        real(kind=dp_t), pointer :: gamp(:,:,:,:)
        real(kind=dp_t), pointer :: Wp(:,:,:,:)
        real(kind=dp_t), pointer :: u0p(:,:,:,:)
        integer     :: nlevs, n, i

        nlevs = mla%nlevel

        do n = 1, nlevs
            do i = 1, nfabs(alpha(n))
                alphap => dataptr(alpha(n),i)
                betap => dataptr(beta(n),i)
                gamp => dataptr(gam(n),i)
                u0p => dataptr(u0(n),i)
                Wp => dataptr(W_lor(n),i)

                u0p(:,:,:,1) = Wp(:,:,:,1) / alphap(:,:,:,1)
            enddo
        enddo

    end subroutine calcu0

    subroutine g(alpha, beta, gam, x, met)
        ! Returns metric components at coordinate x

        real(kind=dp_t)   , intent(inout) :: alpha(0:,0:,0:)
        real(kind=dp_t)   , intent(inout) :: beta(:,0:,0:,0:)
        real(kind=dp_t)   , intent(inout) :: gam(:,0:,0:,0:)
        integer           , intent(inout) :: x(:)

        real(kind=dp_t)   , intent(out) :: met(0:,0:)

        integer             :: i, j, k, m, n
        real(kind=dp_t)     :: eye(1:3,1:3)

        ! make identity matrix
        eye(:,:) = 0.d0
        do i = 1,3
            eye(i,i) = 1.d0
        enddo

        i = x(1)
        j = x(2)
        k = x(3)

        ! initialise
        met(:,:) = 0

        met(0,0) = -alpha(i,j,k)**2

        do m = 1, 3
            do n = 1, 3
                met(0,0) = met(0,0) + beta(m,i,j,k) * beta(n,i,j,k)
                met(m,n) = eye(m,n) * gam(m,i,j,k)
            enddo
        enddo

        met(0,1:) = beta(:,i,j,k)
        met(1:,0) = beta(:,i,j,k)

    end subroutine g

    subroutine christoffels(alpha, beta, gam, chrls, mla)
        ! Calculates christoffel symbols given metric functions.
        ! These are metric dependent, so have only implemented cartesian
        ! weak field here.
        ! TODO: Need instead to provide these in initial data, or at least
        ! have the initial data provide an analytic form of the derivatives
        ! of alpha and beta, K and 3-metric christoffels.
        use probin_module, only : g, Rr, c

        type(multifab)    , intent(in)    :: alpha(:)
        type(multifab)    , intent(in)    :: beta(:)
        type(multifab)    , intent(in)    :: gam(:)
        type(ml_layout)   , intent(inout) :: mla

        real(kind=dp_t)   , intent(inout) :: chrls(:,0:,0:,0:,0:,0:,0:)

        real(kind=dp_t), pointer:: alphap(:,:,:,:)
        real(kind=dp_t), pointer:: betap(:,:,:,:)
        real(kind=dp_t), pointer:: gamp(:,:,:,:)
        integer     :: nlevs, n, i

        nlevs = mla%nlevel

        do n = 1, nlevs
            do i = 1, nfabs(alpha(n))
                alphap => dataptr(alpha(n),i)
                betap => dataptr(beta(n),i)
                gamp => dataptr(gam(n),i)

                ! cartesian weak field
                ! t_tr
                chrls(n,0,0,3,:,:,:) = g / (Rr * alphap(:,:,:,1)**2)
                ! t_rt
                chrls(n,0,3,0,:,:,:) = chrls(n,0,0,3,:,:,:)
                ! r_tt
                chrls(n,3,0,0,:,:,:) = g * alphap(:,:,:,1)**2 / (c**2 * Rr)
                ! r_xx
                chrls(n,3,1,1,:,:,:) = chrls(n,0,0,3,:,:,:)
                ! r_yy
                chrls(n,3,2,2,:,:,:) = chrls(n,0,0,3,:,:,:)
                !r_rr
                chrls(n,3,3,3,:,:,:) = -chrls(n,0,0,3,:,:,:)
                ! r_xr
                chrls(n,3,1,3,:,:,:) = -chrls(n,0,0,3,:,:,:)
                ! r_rx
                chrls(n,3,3,1,:,:,:) = chrls(n,3,1,3,:,:,:)
                ! x_rx
                chrls(n,1,3,1,:,:,:) = -chrls(n,0,0,3,:,:,:)
                ! x_xr
                chrls(n,1,1,3,:,:,:) = chrls(n,1,3,1,:,:,:)
                ! y_yr
                chrls(n,2,2,3,:,:,:) = -chrls(n,0,0,3,:,:,:)
                ! y_ry
                chrls(n,2,3,2,:,:,:) = chrls(n,2,2,3,:,:,:)
            enddo
        enddo

    end subroutine christoffels

    subroutine cons_to_prim(s, u, alpha, s_prim, u_prim, mla)
        use variables, only: rho_comp, rhoh_comp, spec_comp
        ! FIXME: need to import thermal gamma
        use probin_module, only : gamm_therm

        type(multifab), intent(in) :: s(:)
        type(multifab), intent(in) :: u(:)
        type(multifab), intent(in) :: alpha(:)
        type(multifab), intent(in) :: s_prim(:)
        type(multifab), intent(in) :: u_prim(:)
        type(ml_layout)   , intent(inout) :: mla

        real(kind=dp_t), pointer ::  sp(:,:,:,:)
        real(kind=dp_t), pointer ::  up(:,:,:,:)
        real(kind=dp_t), pointer ::  alphap(:,:,:,:)
        real(kind=dp_t), pointer ::  s_primp(:,:,:,:)
        real(kind=dp_t), pointer ::  u_primp(:,:,:,:)

        real(kind=dp_t) :: pmin, pmax, pbar, u0

        integer :: i,j,k,m,n,nlevs,lo(mla%dim),hi(mla%dim)
        logical :: fail

        nlevs = mla%nlevel

        do n=1,nlevs
           do m=1,nfabs(s(n))
               sp => dataptr(s(n), m)
               up => dataptr(u(n), m)
               alphap => dataptr(alpha(n), m)
               s_primp => dataptr(s_prim(n), m)
               u_primp => dataptr(u_prim(n), m)

               do k = lo(3), hi(3)
                   do j = lo(2), hi(2)
                       do i = lo(1), hi(1)
                           pmin = 0.d0
                           pmax = 0.d0
                           pbar = 0.d0

                           call Newton_Raphson(pmin, pmax, pbar, root_find_on_me, &
                           sp(i,j,k,rho_comp), up(i,j,k,:), sp(i,j,k,rhoh_comp), alphap(i,j,k,1), fail)

                           u0 = (sp(i,j,k,rhoh_comp) - sp(i,j,k,rho_comp)) * (gamm_therm - 1.d0) / (gamm_therm * pbar)

                           s_primp(i,j,k,rho_comp) = sp(i,j,k,rho_comp) / u0
                           s_primp(i,j,k,rhoh_comp) = sp(i,j,k,rhoh_comp) / u0
                           u_primp(i,j,k,:) = up(i,j,k,:) * u0
                       enddo
                   enddo
               enddo
           enddo
       enddo


    end subroutine cons_to_prim

    subroutine root_find_on_me(pnew, pbar, dp, D, U, Dh, alpha)
        ! FIXME: need to import thermal gamma
        use probin_module, only : gamm_therm

        real(kind=dp_t), intent(out) :: pnew
        real(kind=dp_t), intent(in) :: pbar
        real(kind=dp_t), intent(inout) :: dp
        real(kind=dp_t), intent(in) :: D
        real(kind=dp_t), intent(in) :: U(3)
        real(kind=dp_t), intent(in) :: Dh
        real(kind=dp_t), intent(in) :: alpha

        integer :: i
        real(kind=dp_t) :: usq

        usq = 0.d0

        do i = 1,3
            usq = usq + U(i)*U(i)
        enddo

        pnew = (gamm_therm - 1.d0) * (Dh/D - 1.d0) * D * sqrt(alpha**2 - usq) / gamm_therm - pbar

        dp = usq * gamm_therm * pnew / Dh - 1.d0


    end subroutine root_find_on_me

    subroutine Newton_Raphson(xl,xh,xm,NonLinearEquation, D, U, Dh, alpha,fail)
        ! Root finder taken from multimodel code
      real(kind=dp_t), intent(inout) :: xl     ! root bracket 1
      real(kind=dp_t), intent(inout) :: xh     ! root bracket 2
      real(kind=dp_t), intent(out)   :: xm     ! root
      external NonLinearEquation        ! subroutine for non-linear equation
      real(kind=dp_t), intent(in) :: D
      real(kind=dp_t), intent(in) :: U(3)
      real(kind=dp_t), intent(in) :: Dh
      real(kind=dp_t), intent(in) :: alpha
      logical, intent(out) :: fail

      integer  :: k                      ! iteration count
      real(kind=dp_t) :: xn, f, df, error
      real(kind=dp_t), parameter :: zero=0.0d0, half=0.5d0
      integer, parameter :: ITERMAX=3000 ! max number of iteration
      real(kind=dp_t), parameter :: TOLERANCE=1.d-14 ! tolerance for error

      call NonLinearEquation(xm,xl,df,D, U, Dh, alpha)  ! Set xm=f(xl)
      call NonLinearEquation(xn,xh,df,D, U, Dh, alpha)  ! Set xn=f(xh)
      if (xm*xn.gt.zero) then           ! check if function changes sign
         xl = 0.0d0
         xh = 10.0d0**(14)
      endif
      if (xm > zero) then               ! Rearrange so that f(xl)< 0.d0 < f(xh)
         xm = xl
         xl = xh
         xh = xm
      endif

      error = abs(xh-xl)                ! error is width of bracketing interval
      xm = half*(xl+xh)                 ! Initialize guess for root
      k = 0                             ! initialize iteration count
      do while (error > TOLERANCE .and. k < ITERMAX) ! iterate
         k = k+1                         ! increment iteration count
         call NonLinearEquation(f,xm,df,D, U, Dh, alpha) ! calculate f(xm), df(xm)
         if (f > zero) then              ! Update root bracketing
            xh = xm                       ! update high
         else
            xl = xm                       ! update low
         endif
         xn = xm - f/df                 ! Tentative newton-Raphson step
         if ( (xn-xl)*(xn-xh) > zero ) then ! check if new root falls within bracket
            xm = half* (xh+xl)            ! if no use a Bisection step
            error = abs(xh-xl)            ! error is width of interval
         else
            error = abs(xn-xm)            ! if within bracket: error is change in root
            xm = xn                       ! update successful Newton-Raphson step
         endif
      enddo
      fail = .true.
      if (error > TOLERANCE) then       ! Check if solution converged
         write(*,*) 'solution did not converge, xn, NonLinearEquation(xn),D(xn)'
         write(*,*) xn, f, error
         fail = .true.
      endif

    end subroutine Newton_Raphson

end module metric_module
