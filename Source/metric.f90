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

    use define_bc_module
    use ml_restrict_fill_module

    implicit none

    private :: root_find_on_me, Newton_Raphson

    public :: make_weak_field, inverse_metric, calcW, calcu0, calcu0_1d, &
              g, christoffels, cons_to_prim, prim_to_cons, prim_to_cons_1d, &
              cons_to_prim_a, prim_to_cons_a, cons_to_prim_edge, &
              prim_to_cons_edge

contains

    subroutine make_weak_field(alpha, beta, gam, mla,dx,the_bc_level)
        ! NOTE: gam is a vector as will assume it is a diagonal matrix
        ! Returns alpha, beta and gamma of weak field metric
        ! FIXME: I have no idea how to implement this in a useful way.

        use probin_module, only : g, Rr, c, prob_lo

        type(multifab) , intent(inout) :: alpha(:)
        type(multifab) , intent(inout) :: beta(:)
        type(multifab) , intent(inout) :: gam(:)
        type(ml_layout), intent(inout) :: mla
        real(kind=dp_t), intent(in   ) :: dx(:,:)
        type(bc_level)    , intent(in   ) :: the_bc_level(:)

        real(kind=dp_t), pointer :: alphap(:,:,:,:)
        real(kind=dp_t), pointer :: betap(:,:,:,:)
        real(kind=dp_t), pointer :: gamp(:,:,:,:)

        real(kind=dp_t) :: rad
        integer :: n, i, j, dm, nlevs
        integer :: lo(mla%dim), hi(mla%dim)

        nlevs = mla%nlevel
        dm = mla%dim

        do n = 1, nlevs
            do j = 1, nfabs(alpha(n))

                alphap => dataptr(alpha(n),j)
                betap => dataptr(beta(n),j)
                gamp => dataptr(gam(n),j)
                lo = lwb(get_box(alpha(n), j))
                hi = upb(get_box(alpha(n), j))

                betap = 0.d0

                do i = lo(3), hi(3)
                    rad = prob_lo(3) + (dble(i) + HALF) * dx(n,3)
                    !rad = (dble(i) + HALF) * dx(n,3)
                    !print *, 'radial coord: ', rad
                    !print *,  2.d0 * (1.d0 - rad / Rr) * g / c**2
                    alphap(:,:,i,1) = c * sqrt(1.d0 - 2.d0 * (1.d0 - rad / Rr) * g / c**2)
                    gamp(:,:,i,:) = 1.d0 + 2.d0 * (1.d0 - rad / Rr) * g / c**2
                enddo

            enddo
        enddo

        ! fill ghosts
        call ml_restrict_and_fill(nlevs,alpha,mla%mba%rr,the_bc_level, &
                                  icomp=1, &
                                  bcomp=1, &
                                  nc=1, &
                                  ng=alpha(1)%ng)
        call ml_restrict_and_fill(nlevs,beta,mla%mba%rr,the_bc_level, &
                                icomp=1, &
                                bcomp=1, &
                                nc=dm, &
                                ng=beta(1)%ng)
        call ml_restrict_and_fill(nlevs,gam,mla%mba%rr,the_bc_level, &
                                  icomp=1, &
                                  bcomp=1, &
                                  nc=dm, &
                                  ng=gam(1)%ng)

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

    subroutine calcW(alpha, beta, gam, u, mla, W_lor,the_bc_level)
        ! Calculates the Lorentz factor
        use probin_module, only : c

        type(multifab)    , intent(in)    :: alpha(:)
        type(multifab)    , intent(in)    :: beta(:)
        type(multifab)    , intent(in)    :: gam(:)
        type(multifab)    , intent(in)    :: u(:)
        type(ml_layout)   , intent(in) :: mla
        type(bc_level)    , intent(in   ) :: the_bc_level(:)

        type(multifab)    , intent(inout)    :: W_lor(:)

        real(kind=dp_t), pointer:: alphap(:,:,:,:)
        real(kind=dp_t), pointer:: betap(:,:,:,:)
        real(kind=dp_t), pointer:: gamp(:,:,:,:)
        real(kind=dp_t), pointer:: up(:,:,:,:)
        real(kind=dp_t), pointer:: Wp(:,:,:,:)
        integer     :: dm, nlevs, n, i, j, k, lo(mla%dim), hi(mla%dim)
        real(kind=dp_t)     :: eye(1:3,1:3)

        ! make identity matrix
        eye(:,:) = 0.d0
        do i = 1,3
            eye(i,i) = 1.d0
        enddo

        nlevs = mla%nlevel
        dm = mla%dim

        do n = 1, nlevs
            do k = 1, nfabs(alpha(n))
                alphap => dataptr(alpha(n),k)
                betap => dataptr(beta(n),k)
                gamp => dataptr(gam(n),k)
                up => dataptr(u(n),k)
                Wp => dataptr(W_lor(n),k)
                lo = lwb(get_box(alpha(n), k))
                hi = upb(get_box(alpha(n), k))

                Wp(:,:,:,1) = ZERO

                do i = 1, 3
                    do j = 1, 3
                        Wp(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1) = Wp(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1) + &
                        eye(i,j) * &
                        gamp(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),i) * &
                        (up(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),i) + &
                        betap(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),i)) * &
                        (up(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),j) + &
                        betap(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),j)) / &
                        alphap(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1)**2
                    enddo
                enddo
                Wp(:,:,:,1) = ONE / sqrt(1.d0 - Wp(:,:,:,1)/c**2)
            enddo
        enddo

        ! fill ghosts
        call ml_restrict_and_fill(nlevs,W_lor,mla%mba%rr,the_bc_level, &
                                  icomp=1, &
                                  bcomp=1, &
                                  nc=1, &
                                  ng=W_lor(1)%ng)

    end subroutine calcW

    subroutine calcW_mac(alpha, beta, gam, u, mla, W_lor,the_bc_level)
        ! Calculates the Lorentz factor
        use probin_module, only : c

        type(multifab)    , intent(in)    :: alpha(:)
        type(multifab)    , intent(in)    :: beta(:)
        type(multifab)    , intent(in)    :: gam(:)
        type(multifab)    , intent(in)    :: u(:,:)
        type(ml_layout)   , intent(in) :: mla
        type(bc_level)    , intent(in   ) :: the_bc_level(:)

        type(multifab)    , intent(inout)    :: W_lor(:)

        real(kind=dp_t), pointer:: alphap(:,:,:,:)
        real(kind=dp_t), pointer:: betap(:,:,:,:)
        real(kind=dp_t), pointer:: gamp(:,:,:,:)
        real(kind=dp_t), pointer:: up(:,:,:,:)
        real(kind=dp_t), pointer:: vp(:,:,:,:)
        real(kind=dp_t), pointer:: wp(:,:,:,:)
        real(kind=dp_t), pointer:: W_lorp(:,:,:,:)
        integer     :: dm, nlevs, n, i, j, k, lo(mla%dim), hi(mla%dim)
        real(kind=dp_t)     :: eye(1:3,1:3)
        real(kind=dp_t), allocatable   :: v(:,:,:,:)

        ! make identity matrix
        eye(:,:) = 0.d0
        do i = 1,3
            eye(i,i) = 1.d0
        enddo

        nlevs = mla%nlevel
        dm = mla%dim

        do n = 1, nlevs
            do k = 1, nfabs(alpha(n))
                alphap => dataptr(alpha(n),k)
                betap => dataptr(beta(n),k)
                gamp => dataptr(gam(n),k)
                up => dataptr(u(n,1),k)
                vp => dataptr(u(n,2),k)
                wp => dataptr(u(n,3),k)
                W_lorp => dataptr(W_lor(n),k)
                lo = lwb(get_box(alpha(n), k))
                hi = upb(get_box(alpha(n), k))

                allocate(v(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1:3))

                W_lorp(:,:,:,1) = ZERO

                v(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1) = up(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1)

                v(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),2) = vp(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1)

                v(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),3) = wp(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1)

                do i = 1, 3
                    do j = 1, 3
                        W_lorp(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1) = W_lorp(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1) + &
                        eye(i,j) * &
                        gamp(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),i) * &
                        (v(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),i) + &
                        betap(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),i)) * &
                        (v(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),j) + &
                        betap(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),j)) / &
                        alphap(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1)**2
                    enddo
                enddo
                W_lorp(:,:,:,1) = ONE / sqrt(1.d0 - W_lorp(:,:,:,1)/c**2)

                deallocate(v)
            enddo
        enddo

        ! fill ghosts
        call ml_restrict_and_fill(nlevs,W_lor,mla%mba%rr,the_bc_level, &
                                  icomp=1, &
                                  bcomp=1, &
                                  nc=1, &
                                  ng=W_lor(1)%ng)

    end subroutine calcW_mac

    subroutine calcu0(alpha, beta, gam, W_lor, u0, mla,the_bc_level)
        ! Calculates timelike coordinate of 3+1 velocity
        ! using Lorentz factor and alpha:
        ! W = alpha * u0

        type(multifab)    , intent(in)    :: alpha(:)
        type(multifab)    , intent(in)    :: beta(:)
        type(multifab)    , intent(in)    :: gam(:)
        type(multifab)    , intent(in)    :: W_lor(:)
        type(multifab)   , intent(inout) :: u0(:)
        type(ml_layout)   , intent(inout) :: mla
        type(bc_level)    , intent(in   ) :: the_bc_level(:)

        real(kind=dp_t), pointer :: alphap(:,:,:,:)
        real(kind=dp_t), pointer :: betap(:,:,:,:)
        real(kind=dp_t), pointer :: gamp(:,:,:,:)
        real(kind=dp_t), pointer :: Wp(:,:,:,:)
        real(kind=dp_t), pointer :: u0p(:,:,:,:)
        integer     :: nlevs, n, i, dm

        nlevs = mla%nlevel
        dm = mla%dim

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

        ! fill ghosts
        call ml_restrict_and_fill(nlevs,u0,mla%mba%rr,the_bc_level, &
                                  icomp=1, &
                                  bcomp=1, &
                                  nc=1, &
                                  ng=u0(1)%ng)

    end subroutine calcu0

    subroutine calcu0_1d(alpha, beta, gam, w0, u0_1d, mla,the_bc_level, dx)
        ! Calculates timelike coordinate of 3+1 velocity
        ! using Lorentz factor and alpha:
        ! W = alpha * u0
        use geometry   , only : nlevs_radial, nr_fine
        use probin_module, only : c
        use average_module   , only : average

        type(multifab)    , intent(in)    :: alpha(:)
        type(multifab)    , intent(in)    :: beta(:)
        type(multifab)    , intent(in)    :: gam(:)
        real(dp_t)        ,  intent(in)   :: w0(:,0:)
        real(dp_t)        ,  intent(inout):: u0_1d(:,0:)
        type(ml_layout)   , intent(inout) :: mla
        type(bc_level)    , intent(in   ) :: the_bc_level(:)
        real(dp_t)    ,  intent(in   ) :: dx(:,:)

        real(kind=dp_t), allocatable :: W_lor(:,:)
        real(kind=dp_t), allocatable :: alpha_1d(:,:)
        real(kind=dp_t), allocatable :: beta_1d(:,:,:)
        real(kind=dp_t), allocatable :: gam_1d(:,:,:)
        integer     :: i, dm

        dm = mla%dim

        allocate(              W_lor(nlevs_radial,0:nr_fine-1))
        allocate(              alpha_1d(nlevs_radial,0:nr_fine-1))
        allocate(              beta_1d(nlevs_radial,0:nr_fine-1,dm))
        allocate(              gam_1d(nlevs_radial,0:nr_fine-1,dm))

        call average(mla, alpha, alpha_1d, dx, 1)

        !print *, 'alpha_1d', alpha_1d

        do i = 1, dm
            call average(mla, beta, beta_1d(:,:,i), dx, i)
            call average(mla, gam, gam_1d(:,:,i), dx, i)
        enddo

        W_lor(:,:) = gam_1d(:,:,3) * (w0(:,:nr_fine-1) + beta_1d(:,:,3))**2 / &
                     alpha_1d(:,:)**2

        W_lor(:,:) = ONE / sqrt(ONE - W_lor(:,:)/c**2)

        u0_1d(:,:) = W_lor(:,:) / alpha_1d(:,:)

    end subroutine calcu0_1d

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
        integer     :: nlevs, n, m, i, j, k, lo(mla%dim), hi(mla%dim)

        nlevs = mla%nlevel

        chrls(:,:,:,:,:,:,:) = ZERO

        do n = 1, nlevs
            do m = 1, nfabs(alpha(n))
                alphap => dataptr(alpha(n),m)
                betap => dataptr(beta(n),m)
                gamp => dataptr(gam(n),m)
                lo =  lwb(get_box(alpha(n), m))
                hi =  upb(get_box(alpha(n), m))

                do k = lo(3), hi(3)
                    do j = lo(2), hi(2)
                        do i = lo(1), hi(1)
                            ! t_tr
                            chrls(n,0,0,3,i,j,k) = g / (Rr * alphap(i,j,k,1)**2)
                            ! r_tt
                            chrls(n,3,0,0,i,j,k) = g * alphap(i,j,k,1)**2 / (c**2 * Rr)
                        enddo
                    enddo
                enddo

                ! cartesian weak field
                ! t_rt
                chrls(n,0,3,0,:,:,:) = chrls(n,0,0,3,:,:,:)
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

    subroutine cons_to_prim(s, u, alpha, beta, gam, s_prim, u_prim, mla,the_bc_level)
        use variables, only: rho_comp, rhoh_comp, spec_comp, nspec, nscal
        ! FIXME: need to import thermal gamma
        use probin_module, only : gamm_therm

        type(multifab), intent(in) :: s(:)
        type(multifab), intent(in) :: u(:)
        type(multifab), intent(in) :: alpha(:)
        type(multifab), intent(in) :: beta(:)
        type(multifab), intent(in) :: gam(:)
        type(multifab), intent(inout) :: s_prim(:)
        type(multifab), intent(inout) :: u_prim(:)
        type(ml_layout)   , intent(in) :: mla
        type(bc_level)    , intent(in   ) :: the_bc_level(:)

        real(kind=dp_t), pointer ::  sp(:,:,:,:)
        real(kind=dp_t), pointer ::  up(:,:,:,:)
        real(kind=dp_t), pointer ::  alphap(:,:,:,:)
        real(kind=dp_t), pointer ::  betap(:,:,:,:)
        real(kind=dp_t), pointer ::  gamp(:,:,:,:)
        real(kind=dp_t), pointer ::  s_primp(:,:,:,:)
        real(kind=dp_t), pointer ::  u_primp(:,:,:,:)

        real(kind=dp_t) :: pmin, pmax, pbar, u0

        integer :: i,j,k,l,m,n,nlevs,lo(mla%dim),hi(mla%dim),dm
        logical :: fail

        nlevs = mla%nlevel
        dm = mla%dim

        do n=1,nlevs
           do m=1,nfabs(s(n))

               sp => dataptr(s(n), m)
               up => dataptr(u(n), m)
               alphap => dataptr(alpha(n), m)
               betap => dataptr(beta(n), m)
               gamp => dataptr(gam(n), m)
               s_primp => dataptr(s_prim(n), m)
               u_primp => dataptr(u_prim(n), m)
               lo =  lwb(get_box(s(n), m))
               hi =  upb(get_box(s(n), m))

               !print *, 's', sp(:,:,:,rhoh_comp)

               do k = lo(3), hi(3)
                   do j = lo(2), hi(2)
                       do i = lo(1), hi(1)

                           s_primp(i,j,k,:) = sp(i,j,k,:)
                           !pmin = (Dh - D) * (gamma - 1.) / gamma
                           !pmax = (gamma - 1.) * Dh
                           pmin = (sp(i,j,k,rhoh_comp) - sp(i,j,k,rho_comp)) * (gamm_therm - 1.d0) / gamm_therm
                           if (pmin < 0.d0) then
                               pmin = 0.d0
                           endif
                           pmax = (gamm_therm - 1.d0) * sp(i,j,k,rhoh_comp)
                           pbar = 0.5d0 * (pmin + pmax)

                           !print *, 'pmax', pmax

                           call Newton_Raphson(pmin, pmax, pbar, root_find_on_me, &
                           sp(i,j,k,rho_comp), up(i,j,k,:), &
                           sp(i,j,k,rhoh_comp), alphap(i,j,k,1), &
                           betap(i,j,k,:), gamp(i,j,k,:), fail)

                           u0 = (sp(i,j,k,rhoh_comp) - sp(i,j,k,rho_comp)) * (gamm_therm - 1.d0) / (gamm_therm * pbar)

                           !print *, 'pmax', pmax

                           s_primp(i,j,k,rho_comp) = sp(i,j,k,rho_comp) / u0
                           s_primp(i,j,k,rhoh_comp) = sp(i,j,k,rhoh_comp) / u0
                           do l = 0, nspec-1
                               s_primp(i,j,k,spec_comp+l) = sp(i,j,k,spec_comp+l) / u0
                           enddo
                           u_primp(i,j,k,:) = up(i,j,k,:) * u0
                       enddo
                   enddo
               enddo
           enddo
       enddo

       !print *, 'sprim', s_primp(:,:,:,1)
       !print *, 'gamm_therm', gamm_therm

       ! fill ghosts
       call ml_restrict_and_fill(nlevs,s_prim,mla%mba%rr,the_bc_level, &
                                 icomp=rho_comp, &
                                 bcomp=dm+rho_comp, &
                                 nc=nscal, &
                                 ng=s_prim(1)%ng)

       call ml_restrict_and_fill(nlevs,u_prim,mla%mba%rr,the_bc_level, &
                                 icomp=1, &
                                 bcomp=1, &
                                 nc=dm, &
                                 ng=u_prim(1)%ng)

        !print *, 'sprim', s_primp(:,:,:,1)

    end subroutine cons_to_prim

    subroutine cons_to_prim_edge(sedge, umac, u, alpha, beta, gam, s_prim_edge, u_prim_edge, mla,the_bc_level)
        use variables, only: rho_comp, rhoh_comp, spec_comp, nspec, nscal
        ! FIXME: need to import thermal gamma
        use probin_module, only : gamm_therm

        type(multifab), intent(in) :: sedge(:,:)
        type(multifab), intent(in) :: umac(:,:)
        type(multifab), intent(in) :: u(:)
        type(multifab), intent(in) :: alpha(:)
        type(multifab), intent(in) :: beta(:)
        type(multifab), intent(in) :: gam(:)
        type(multifab), intent(inout) :: s_prim_edge(:,:)
        type(multifab), intent(inout) :: u_prim_edge(:,:)
        type(ml_layout)   , intent(in) :: mla
        type(bc_level)    , intent(in   ) :: the_bc_level(:)

        real(kind=dp_t), pointer ::  sp(:,:,:,:)
        real(kind=dp_t), pointer ::  umacp(:,:,:,:)
        real(kind=dp_t), pointer ::  up(:,:,:,:)
        real(kind=dp_t), pointer ::  alphap(:,:,:,:)
        real(kind=dp_t), pointer ::  betap(:,:,:,:)
        real(kind=dp_t), pointer ::  gamp(:,:,:,:)
        real(kind=dp_t), pointer ::  s_primp(:,:,:,:)
        real(kind=dp_t), pointer ::  u_primp(:,:,:,:)

        real(kind=dp_t) :: pmin, pmax, pbar, u0
        real(kind=dp_t), allocatable :: uedge(:,:,:,:)

        integer :: i,j,k,l,m,n,nlevs,lo(mla%dim),hi(mla%dim),dm,comp
        logical :: fail

        nlevs = mla%nlevel
        dm = mla%dim

        do n=1,nlevs
           do comp=1,dm
               do m=1,nfabs(sedge(n,comp))
                   sp => dataptr(sedge(n,comp), m)
                   umacp => dataptr(umac(n,comp), m)
                   up => dataptr(u(n), m)
                   alphap => dataptr(alpha(n), m)
                   betap => dataptr(beta(n), m)
                   gamp => dataptr(gam(n), m)
                   s_primp => dataptr(s_prim_edge(n,comp), m)
                   u_primp => dataptr(u_prim_edge(n,comp), m)
                   lo =  lwb(get_box(sedge(n,comp), m))
                   hi =  upb(get_box(sedge(n,comp), m))

                   allocate(uedge(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3), 1:dm))

                   !print *, 's', sp(:,:,:,rhoh_comp)

                   do k = lo(3), hi(3)
                       do j = lo(2), hi(2)
                           do i = lo(1), hi(1)

                               uedge(i,j,k,:) = up(i,j,k,:)
                               uedge(i,j,k,comp) = umacp(i,j,k,1)

                               s_primp(i,j,k,:) = sp(i,j,k,:)
                               !u_primp(i,j,k,:) = up(i,j,k,:)
                               !pmin = (Dh - D) * (gamma - 1.) / gamma
                               !pmax = (gamma - 1.) * Dh
                               pmin = (sp(i,j,k,rhoh_comp) - sp(i,j,k,rho_comp)) * (gamm_therm - 1.d0) / gamm_therm
                               if (pmin < 0.d0) then
                                   pmin = 0.d0
                               endif
                               pmax = (gamm_therm - 1.d0) * sp(i,j,k,rhoh_comp)
                               pbar = 0.5d0 * (pmin + pmax)

                               !print *, 'pmax', pmax

                               call Newton_Raphson(pmin, pmax, pbar, root_find_on_me, &
                               sp(i,j,k,rho_comp), uedge(i,j,k,:), &
                               sp(i,j,k,rhoh_comp), alphap(i,j,k,1), &
                               betap(i,j,k,:), gamp(i,j,k,:), fail)
                               !sp(i,j,k,rho_comp), umacp(i,j,k,:), &
                               !sp(i,j,k,rhoh_comp), alphap(i,j,k,1), &
                               !betap(i,j,k,:), gamp(i,j,k,:), fail)

                               u0 = (sp(i,j,k,rhoh_comp) - sp(i,j,k,rho_comp)) * (gamm_therm - 1.d0) / (gamm_therm * pbar)

                               !print *, 'pmax', pmax

                               s_primp(i,j,k,rho_comp) = sp(i,j,k,rho_comp) / u0
                               s_primp(i,j,k,rhoh_comp) = sp(i,j,k,rhoh_comp) / u0
                               do l = 0, nspec-1
                                   s_primp(i,j,k,spec_comp+l) = sp(i,j,k,spec_comp+l) / u0
                               enddo
                               u_primp(i,j,k,1) = umacp(i,j,k,1) * u0
                           enddo
                       enddo
                   enddo

                   deallocate(uedge)
               enddo
           enddo
        enddo



       !print *, 'sprim', s_primp(:,:,:,1)
       !print *, 'gamm_therm', gamm_therm
       do comp=1,dm
           ! fill ghosts
           call ml_restrict_and_fill(nlevs,s_prim_edge(:,comp),mla%mba%rr,the_bc_level, &
                                     icomp=rho_comp, &
                                     bcomp=dm+rho_comp, &
                                     nc=nscal, &
                                     ng=s_prim_edge(1,comp)%ng)

           call ml_restrict_and_fill(nlevs,u_prim_edge(:,comp),mla%mba%rr,the_bc_level, &
                                     icomp=1, &
                                     bcomp=1, &
                                     nc=1, &
                                     ng=u_prim_edge(1,comp)%ng)

            !print *, 'sprim', s_primp(:,:,:,1)
        enddo

    end subroutine cons_to_prim_edge

    subroutine cons_to_prim_a(s, u, alpha, beta, gam, s_prim,lo,hi,ng_s,ng_u)
        ! Works on arrays
        ! Must fill primitive ghosts after calling this function
        ! Does not convert u as not needed in any of the functions that call
        ! this.
        use variables, only: rho_comp, rhoh_comp, spec_comp, nspec
        ! FIXME: need to import thermal gamma
        use probin_module, only : gamm_therm


        integer, intent(in) :: lo(:),hi(:),ng_s,ng_u
        real(kind=dp_t), intent(in) :: s(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:,:)
        real(kind=dp_t), intent(in) :: u(lo(1)-ng_u:,lo(2)-ng_u:,lo(3)-ng_u:,:)
        real(kind=dp_t), intent(in) :: alpha(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:)
        real(kind=dp_t), intent(in) :: beta(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:,:)
        real(kind=dp_t), intent(in) :: gam(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:,:)
        real(kind=dp_t), intent(inout) :: s_prim(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:,:)

        real(kind=dp_t) :: pmin, pmax, pbar, u0

        integer :: i,j,k,l
        logical :: fail


       !print *, 's', sp(:,:,:,rhoh_comp)

       ! initialise by simply copying - assume rest of the components are unchanged
       s_prim(:,:,:,:) = s(:,:,:,:)

       do k = lo(3), hi(3)
           do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                   !pmin = (Dh - D) * (gamma - 1.) / gamma
                   !pmax = (gamma - 1.) * Dh
                   pmin = (s(i,j,k,rhoh_comp) - s(i,j,k,rho_comp)) * (gamm_therm - 1.d0) / gamm_therm
                   if (pmin < 0.d0) then
                       pmin = 0.d0
                   endif
                   pmax = (gamm_therm - 1.d0) * s(i,j,k,rhoh_comp)
                   pbar = 0.5d0 * (pmin + pmax)

                   !print *, 'pmax', pmax

                   !print *, 'hi'

                   call Newton_Raphson(pmin, pmax, pbar, root_find_on_me, &
                   s(i,j,k,rho_comp), u(i,j,k,:), &
                   s(i,j,k,rhoh_comp), alpha(i,j,k), &
                   beta(i,j,k,:), gam(i,j,k,:), fail)

                   u0 = (s(i,j,k,rhoh_comp) - s(i,j,k,rho_comp)) * (gamm_therm - 1.d0) / (gamm_therm * pbar)

                   !print *, 'pmax', pmax

                   s_prim(i,j,k,rho_comp) = s(i,j,k,rho_comp) / u0
                   s_prim(i,j,k,rhoh_comp) = s(i,j,k,rhoh_comp) / u0
                   do l = 0, nspec-1
                       s_prim(i,j,k,spec_comp+l) = s(i,j,k,spec_comp+l) / u0
                   enddo
               enddo
           enddo
       enddo


    end subroutine cons_to_prim_a

    subroutine root_find_on_me(pnew, pbar, dp, D, U, Dh, alpha, gam, beta)
        ! FIXME: need to import thermal gamma
        use probin_module, only : gamm_therm, c

        real(kind=dp_t), intent(out) :: pnew
        real(kind=dp_t), intent(in) :: pbar
        real(kind=dp_t), intent(inout) :: dp
        real(kind=dp_t), intent(in) :: D
        real(kind=dp_t), intent(in) :: U(3)
        real(kind=dp_t), intent(in) :: Dh
        real(kind=dp_t), intent(in) :: alpha
        real(kind=dp_t), intent(in) :: gam(3)
        real(kind=dp_t), intent(in) :: beta(3)

        integer :: i,j
        real(kind=dp_t) :: usq, h, rho, W_lor, u0
        real(kind=dp_t)     :: eye(1:3,1:3)

        ! make identity matrix
        eye(:,:) = 0.d0
        do i = 1,3
            eye(i,i) = 1.d0
        enddo

        if (pbar > 0.d0) then

            usq = 0.d0

            usq = sum(U**2)

            if (usq < c**2) then
                W_lor = 0.d0
                do i=1,3
                    do j=1,3
                        W_lor = W_lor + eye(i,j) * gam(i) * (u(i) + beta(i)) &
                                * (u(j) + beta(j)) / alpha**2
                    enddo
                enddo
                W_lor = ONE / sqrt(1.d0 - W_lor/c**2)
                u0 = W_lor / alpha
            else
                u0 = 1.d0 / c
            endif

            !print *, alpha

            rho = D / u0
            h = Dh / D

            pnew = (gamm_therm - 1.d0) * (h - 1.d0) * rho / gamm_therm - pbar

            dp = usq * gamm_therm * pnew / Dh - 1.d0
        else
            pnew = c
            dp = c

        endif


    end subroutine root_find_on_me

    subroutine Newton_Raphson(xl,xh,xm,NonLinearEquation, D, U, Dh, alpha,beta,gam,fail)
        ! Root finder taken from multimodel code
      real(kind=dp_t), intent(inout) :: xl     ! root bracket 1
      real(kind=dp_t), intent(inout) :: xh     ! root bracket 2
      real(kind=dp_t), intent(out)   :: xm     ! root
      external NonLinearEquation        ! subroutine for non-linear equation
      real(kind=dp_t), intent(in) :: D
      real(kind=dp_t), intent(in) :: U(3)
      real(kind=dp_t), intent(in) :: Dh
      real(kind=dp_t), intent(in) :: alpha
      real(kind=dp_t), intent(in) :: gam(3)
      real(kind=dp_t), intent(in) :: beta(3)
      logical, intent(out) :: fail

      integer  :: k                      ! iteration count
      real(kind=dp_t) :: xn, f, df, error
      real(kind=dp_t), parameter :: zero=0.0d0, half=0.5d0
      integer, parameter :: ITERMAX=3000 ! max number of iteration
      real(kind=dp_t), parameter :: TOLERANCE=1.d-14 ! tolerance for error

      call NonLinearEquation(xm,xl,df,D, U, Dh, alpha, beta,gam)  ! Set xm=f(xl)
      call NonLinearEquation(xn,xh,df,D, U, Dh, alpha, beta, gam)  ! Set xn=f(xh)
      if (xm*xn.gt.zero) then           ! check if function changes sign
          !print *,'OOPS'
         xl = 0.0d0
         xh = 10.0d0**(30)
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
         call NonLinearEquation(f,xm,df,D, U, Dh, alpha, beta, gam) ! calculate f(xm), df(xm)
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
      fail = .false.
      !if (error > TOLERANCE) then       ! Check if solution converged
        ! write(*,*) 'solution did not converge, xn, NonLinearEquation(xn),D(xn)'
         !write(*,*) xn, f, error
         !fail = .true.
      !endif

    end subroutine Newton_Raphson

    subroutine prim_to_cons(s, u, u0, s_prim, u_prim, mla,the_bc_level)
        use variables, only: rho_comp, rhoh_comp, spec_comp, nspec,nscal

        type(multifab), intent(inout) :: s(:)
        type(multifab), intent(inout) :: u(:)
        type(multifab), intent(in) :: u0(:)
        type(multifab), intent(in) :: s_prim(:)
        type(multifab), intent(in) :: u_prim(:)
        type(ml_layout)   , intent(in) :: mla
        type(bc_level)    , intent(in   ) :: the_bc_level(:)

        real(kind=dp_t), pointer ::  sp(:,:,:,:)
        real(kind=dp_t), pointer ::  up(:,:,:,:)
        real(kind=dp_t), pointer ::  u0p(:,:,:,:)
        real(kind=dp_t), pointer ::  s_primp(:,:,:,:)
        real(kind=dp_t), pointer ::  u_primp(:,:,:,:)

        integer :: i,m,n,nlevs,dm

        nlevs = mla%nlevel
        dm = mla%dim

        do n=1,nlevs
           do m=1,nfabs(s(n))
               sp => dataptr(s(n), m)
               up => dataptr(u(n), m)
               u0p => dataptr(u0(n), m)
               s_primp => dataptr(s_prim(n), m)
               u_primp => dataptr(u_prim(n), m)

               ! intialise by simply copying - assume rest of the components are unchanged
               sp(:,:,:,:) = s_primp(:,:,:,:)

               sp(:,:,:,rho_comp) = s_primp(:,:,:,rho_comp) * u0p(:,:,:,1)
               sp(:,:,:,rhoh_comp) = s_primp(:,:,:,rhoh_comp) * u0p(:,:,:,1)
               do i = 0, nspec-1
                   sp(:,:,:,spec_comp+i) = s_primp(:,:,:,spec_comp+i) * u0p(:,:,:,1)
               enddo
               do i = 1, mla%dim
                   up(:,:,:,i) = u_primp(:,:,:,i) / u0p(:,:,:,1)
               enddo
           enddo
       enddo

       ! fill ghosts
       call ml_restrict_and_fill(nlevs,s,mla%mba%rr,the_bc_level, &
                                 icomp=rho_comp, &
                                 bcomp=dm+rho_comp, &
                                 nc=nscal, &
                                 ng=s(1)%ng)

       call ml_restrict_and_fill(nlevs,u,mla%mba%rr,the_bc_level, &
                                 icomp=1, &
                                 bcomp=1, &
                                 nc=dm, &
                                 ng=u(1)%ng)


   end subroutine prim_to_cons

   subroutine prim_to_cons_edge(sedge, umac, u0, s_prim_edge, u_prim_edge, mla,the_bc_level)
       use variables, only: rho_comp, rhoh_comp, spec_comp, nspec,nscal

       type(multifab), intent(inout) :: sedge(:,:)
       type(multifab), intent(inout) :: umac(:,:)
       type(multifab), intent(in) :: u0(:)
       type(multifab), intent(in) :: s_prim_edge(:,:)
       type(multifab), intent(in) :: u_prim_edge(:,:)
       type(ml_layout)   , intent(inout) :: mla
       type(bc_level)    , intent(in   ) :: the_bc_level(:)

       real(kind=dp_t), pointer ::  sp(:,:,:,:)
       real(kind=dp_t), pointer ::  up(:,:,:,:)
       real(kind=dp_t), pointer ::  u0p(:,:,:,:)
       real(kind=dp_t), pointer ::  s_primp(:,:,:,:)
       real(kind=dp_t), pointer ::  u_primp(:,:,:,:)

       integer :: i,m,n,nlevs,dm,comp

       nlevs = mla%nlevel
       dm = mla%dim

       do n=1,nlevs
           do comp=1,dm
              do m=1,nfabs(sedge(n,comp))
                  sp => dataptr(sedge(n,comp), m)
                  up => dataptr(umac(n,comp), m)
                  u0p => dataptr(u0(n), m)
                  s_primp => dataptr(s_prim_edge(n,comp), m)
                  u_primp => dataptr(u_prim_edge(n,comp), m)

                  ! intialise by simply copying - assume rest of the components are unchanged
                  sp(:,:,:,:) = s_primp(:,:,:,:)

                  sp(:,:,:,rho_comp) = s_primp(:,:,:,rho_comp) * u0p(:,:,:,1)
                  sp(:,:,:,rhoh_comp) = s_primp(:,:,:,rhoh_comp) * u0p(:,:,:,1)
                  do i = 0, nspec-1
                      sp(:,:,:,spec_comp+i) = s_primp(:,:,:,spec_comp+i) * u0p(:,:,:,1)
                  enddo
                  up(:,:,:,1) = u_primp(:,:,:,1) / u0p(:,:,:,1)
              enddo
          enddo
      enddo

      do comp=1,dm
          ! fill ghosts
          call ml_restrict_and_fill(nlevs,sedge(:,comp),mla%mba%rr,the_bc_level, &
                                    icomp=rho_comp, &
                                    bcomp=dm+rho_comp, &
                                    nc=nscal, &
                                    ng=sedge(1,comp)%ng)

          call ml_restrict_and_fill(nlevs,umac(:,comp),mla%mba%rr,the_bc_level, &
                                    icomp=1, &
                                    bcomp=1, &
                                    nc=dm, &
                                    ng=umac(1,comp)%ng)
     enddo


  end subroutine prim_to_cons_edge

   subroutine prim_to_cons_a(s, u0, s_prim)
       use variables, only: rho_comp, rhoh_comp, spec_comp, nspec,nscal

       real(kind=dp_t), intent(inout) :: s(:,:,:,:)
       real(kind=dp_t), intent(in) :: u0(:,:,:)
       real(kind=dp_t), intent(in) :: s_prim(:,:,:,:)

       integer :: i

      ! intialise by simply copying - assume rest of the components are unchanged
      s(:,:,:,:) = s_prim(:,:,:,:)

      s(:,:,:,rho_comp) = s_prim(:,:,:,rho_comp) * u0(:,:,:)
      s(:,:,:,rhoh_comp) = s_prim(:,:,:,rhoh_comp) * u0(:,:,:)
      do i = 0, nspec-1
          s(:,:,:,spec_comp+i) = s_prim(:,:,:,spec_comp+i) * u0(:,:,:)
      enddo
  end subroutine prim_to_cons_a

   subroutine prim_to_cons_1d(s, u0_1d, s_prim)
       use variables, only: rho_comp, rhoh_comp, spec_comp, nspec

       real(kind=dp_t), intent(inout) :: s(:,:,:)
       real(kind=dp_t), intent(in) :: u0_1d(:,:)
       real(kind=dp_t), intent(in) :: s_prim(:,:,:)

       integer :: i

       ! intialise by simply copying - assume rest of the components are unchanged
       s(:,:,:) = s_prim(:,:,:)

       s(:,:,rho_comp) = s_prim(:,:,rho_comp) * u0_1d(:,:)
       s(:,:,rhoh_comp) = s_prim(:,:,rhoh_comp) * u0_1d(:,:)

       do i = 0, nspec-1
           s(:,:,spec_comp+i) = s_prim(:,:,spec_comp+i) * u0_1d(:,:)
       enddo


  end subroutine prim_to_cons_1d

end module metric_module
