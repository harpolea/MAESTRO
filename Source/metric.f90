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

    private

    public :: make_weak_field, inverse_metric, calcW, calcu0, g, christoffels

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

end module metric_module
