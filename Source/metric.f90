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

    public

contains

    subroutine inverse_metric(alpha, beta, gam, g_up)

        real(kind=dp_t)   , intent(inout) :: alpha(0:,0:,0:)
        real(kind=dp_t)   , intent(inout) :: beta(:,0:,0:,0:)
        real(kind=dp_t)   , intent(inout) :: gam(:,:,0:,0:,0:)

        real(kind=dp_t)   , intent(out) :: g_up(0:,0:,0:,0:,0:)

        integer             :: i, j

        g_up(0,0,:,:,:) = -1 / alpha(:,:,:)**2
        g_up(0,1:,:,:,:) = beta(1:,:,:,:) / alpha(:,:,:)**2
        g_up(1:,0,:,:,:) = beta(1:,:,:,:) / alpha(:,:,:)**2
        do i = 1, 3
            do j = 1, 3
                g_up(i,j,:,:,:) = gam(i,j,:,:,:) - beta(i,:,:,:) * beta(j,:,:,:) / alpha(:,:,:)**2
            enddo
        enddo

    end subroutine inverse_metric

    subroutine calcW(c, alpha, beta, gam, u, v, w, W)

        real(kind=dp_t)   , intent(in) :: c
        real(kind=dp_t)   , intent(inout) :: alpha(0:,0:,0:)
        real(kind=dp_t)   , intent(inout) :: beta(:,0:,0:,0:)
        real(kind=dp_t)   , intent(inout) :: gam(:,:,0:,0:,0:)
        real(kind=dp_t)   , intent(inout) :: u(:,0:,0:,0:)
        real(kind=dp_t)   , intent(inout) :: v(:,0:,0:,0:)
        real(kind=dp_t)   , intent(inout) :: w(:,0:,0:,0:)

        real(kind=dp_t)   , intent(out) :: W(0:,0:,0:)

        real(kind=dp_t)     :: V(:,0:,0:,0:)
        integer             :: i, j

        V(1,:,:,:) = (u(:,:,:) + beta(1,:,:,:)) / alpha(:,:,:)
        V(2,:,:,:) = (v(:,:,:) + beta(2,:,:,:)) / alpha(:,:,:)
        V(3,:,:,:) = (w(:,:,:) + beta(3,:,:,:)) / alpha(:,:,:)

        ! initialise
        W(:,:,:) = 0

        do i = 1, 3
            do j = 1, 3
                W(:,:,:) = W(:,:,:) + gam(i,j,:,:,:) * V(i,:,:,:) * V(j,:,:,:)
            enddo
        enddo

        W(:,:,:) = 1 / sqrt(1 - W(:,:,:)/c**2)

    end subroutine calcW

    subroutine calcu0(c, alpha, beta, gam, u, v, w, u0)

        real(kind=dp_t)   , intent(in) :: c
        real(kind=dp_t)   , intent(inout) :: alpha(0:,0:,0:)
        real(kind=dp_t)   , intent(inout) :: beta(:,0:,0:,0:)
        real(kind=dp_t)   , intent(inout) :: gam(:,:,0:,0:,0:)
        real(kind=dp_t)   , intent(inout) :: u(:,0:,0:,0:)
        real(kind=dp_t)   , intent(inout) :: v(:,0:,0:,0:)
        real(kind=dp_t)   , intent(inout) :: w(:,0:,0:,0:)

        real(kind=dp_t)   , intent(out) :: u0(0:,0:,0:)

        real(kind=dp_t)   :: W(0:,0:,0:)

        calcW(c, alpha, beta, gam, u, v, w, W)

        u0(:,:,:) = W(:,:,:) / alpha(:,:,:)

    end subroutine calcu0

    subroutine g(alpha, beta, gam, x, met)

        real(kind=dp_t)   , intent(inout) :: alpha(0:,0:,0:)
        real(kind=dp_t)   , intent(inout) :: beta(:,0:,0:,0:)
        real(kind=dp_t)   , intent(inout) :: gam(:,:,0:,0:,0:)
        integer           , intent(inout) :: x(:)

        real(kind=dp_t)   , intent(out) :: met(0:,0:)

        integer             :: i, j, k, m, n

        i = x(1)
        j = x(2)
        k = x(3)

        ! initialise
        met(:,:) = 0

        met(0,0) = -alpha(i,j,k)**2

        do m = 1, 3
            do n = 1, 3
                met(0,0) = met(0,0) + beta(m,i,j,k) * beta(n,i,j,k)
            enddo
        enddo

        met(0,1:) = beta(:,i,j,k)
        met(1:,0) = beta(:,i,j,k)
        met(1:,1:) = gam(:,:,i,j,k)

    end subroutine g

    subroutine christoffels(alpha, beta, gam, g, R, c, chrls)
        ! these are metric dependent

        real(kind=dp_t)   , intent(in) :: g, R, c

        real(kind=dp_t)   , intent(inout) :: alpha(0:,0:,0:)
        real(kind=dp_t)   , intent(inout) :: beta(:,0:,0:,0:)
        real(kind=dp_t)   , intent(inout) :: gam(:,:,0:,0:,0:)

        real(kind=dp_t)   , intent(out) :: chrls(0:,0:,0:,0:,0:,0:)

        real(kind=dp_t)         :: g_ralph2(0:,0:,0:)

        g_ralph2(:,:,:) = g / (R * alpha(:,:,:)**2)

        ! cartesian weak field
        ! t_tr
        chrls(0,0,3,:,:,:) = g_ralph2(:,:,:)
        ! t_rt
        chrls(0,3,0,:,:,:) = chrls(0,0,3,:,:,:)
        ! r_tt
        chrls(3,0,0,:,:,:) = g * alpha(:,:,:)**2 / (c**2 * R)
        ! r_xx
        chrls(3,1,1,:,:,:) = g_ralph2(:,:,:)
        ! r_yy
        chrls(3,2,2,:,:,:) = g_ralph2(:,:,:)
        !r_rr
        chrls(3,3,3,:,:,:) = -g_ralph2(:,:,:)
        ! r_xr
        chrls(3,1,3,:,:,:) = -g_ralph2(:,:,:)
        ! r_rx
        chrls(3,3,1,:,:,:) = chrls(3,1,3,:,:,:)
        ! x_rx
        chrls(1,3,1,:,:,:) = -g_ralph2(:,:,:)
        ! x_xr
        chrls(1,1,3,:,:,:) = chrls(1,3,1,:,:,:)
        ! y_yr
        chrls(2,2,3,:,:,:) = -g_ralph2(:,:,:)
        ! y_ry
        chrls(2,3,2,:,:,:) = chrls(2,2,3,:,:,:)

    end subroutine christoffels

end module metric_module
