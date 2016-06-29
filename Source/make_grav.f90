module make_grav_module

  use bl_types

  implicit none

  private

  public :: make_dpdr_cell, make_dpdr_edge, make_dpdr_edge_uniform

contains



  subroutine make_dpdr_cell(dpdr_cell, Dh0, p0, u0_1d)

    use bl_constants_module
    use geometry, only: spherical, nr_fine, r_cc_loc, r_edge_loc, &
         nlevs_radial, nr, numdisjointchunks, &
         r_start_coord, r_end_coord
    use probin_module, only: base_cutoff_density, &
         ref_ratio, g, Rr, c
    use fundamental_constants_module, only: Gconst
    use restrict_base_module
    use multifab_module

    ! compute the base state gravitational acceleration at the cell
    ! centers.  The base state uses 0-based indexing, so grav_cell
    ! does too.

    real(kind=dp_t), intent(in   ) ::      Dh0(:,0:)

    real(kind=dp_t), intent(  out) :: dpdr_cell(:,0:)
    real(kind=dp_t), intent(in) ::   p0(:,0:)
    real(kind=dp_t),  intent(in) ::   u0_1d(:,0:)

    ! Local variables
    integer                      :: r, n, i,j
    real(kind=dp_t), allocatable :: m(:,:)
    real(kind=dp_t)              :: term1, term2

    if (spherical .eq. 0) then

          ! we are doing a plane-parallel geometry with a 1/r**2
          ! gravitational acceleration.  The mass is assumed to be
          ! at the origin.  The mass in the computational domain
          ! does not contribute to the gravitational acceleration.
          do n=1,nlevs_radial
             do r = 0, nr(n)-1
                dpdr_cell(n,r) = -(Dh0(n,r) + p0(n,r) * u0_1d(n,r)) * g / &
                    ((c**2 * (Rr + TWO * r_cc_loc(n,r)) - TWO * g) * u0_1d(n,r))

                !-Gconst*planar_invsq_mass / r_cc_loc(n,r)**2
             enddo
          enddo

          !print *, 'c**2 * Dh0(n,r)', c**2 * Dh0(:,:)
          !print *, 'p0(n,r) * u0_1d(n,r)', p0(:,:) * u0_1d(:,:)
          !print *, 'dpdr cell', dpdr_cell


    else  ! spherical = 1

       allocate(m(1,0:nr_fine-1))

       m(1,0) = FOUR3RD*M_PI*Dh0(1,0)*r_cc_loc(1,0)**3
       dpdr_cell(1,0) = -Gconst * m(1,0) / r_cc_loc(1,0)**2

       do r=1,nr_fine-1

          ! the mass is defined at the cell-centers, so to compute
          ! the mass at the current center, we need to add the
          ! contribution of the upper half of the zone below us and
          ! the lower half of the current zone.

          ! don't add any contributions from outside the star --
          ! i.e.  rho < base_cutoff_density
          if (Dh0(1,r-1) > base_cutoff_density) then
             term1 = FOUR3RD*M_PI*Dh0(1,r-1) * &
                  (r_edge_loc(1,r) - r_cc_loc(1,r-1)) * &
                  (r_edge_loc(1,r)**2 + &
                   r_edge_loc(1,r)*r_cc_loc(1,r-1) + &
                   r_cc_loc(1,r-1)**2)
          else
             term1 = ZERO
          endif

          if (Dh0(1,r) > base_cutoff_density) then
             term2 = FOUR3RD*M_PI*Dh0(1,r  )*&
                  (r_cc_loc(1,r) - r_edge_loc(1,r  )) * &
                  (r_cc_loc(1,r)**2 + &
                   r_cc_loc(1,r)*r_edge_loc(1,r  ) + &
                   r_edge_loc(1,r  )**2)
          else
             term2 = ZERO
          endif

          m(1,r) = m(1,r-1) + term1 + term2

          dpdr_cell(1,r) = -Gconst * m(1,r) / r_cc_loc(1,r)**2

       enddo

       deallocate(m)

    end if

  end subroutine make_dpdr_cell

  subroutine make_dpdr_edge(dpdr_edge, Dh0, p0, u0_1d)

    use bl_constants_module
    use geometry, only: spherical, r_edge_loc, nr_fine, nlevs_radial, nr, &
         numdisjointchunks, r_start_coord, r_end_coord
    use probin_module, only: base_cutoff_density, &
         ref_ratio, g, Rr, c
    use fundamental_constants_module, only: Gconst
    use restrict_base_module

    ! compute the base state gravity at the cell edges (grav_edge(1)
    ! is the gravitational acceleration at the left edge of zone 1).
    ! The base state uses 0-based indexing, so grav_edge does too.

    real(kind=dp_t), intent(  out) :: dpdr_edge(:,0:)
    real(kind=dp_t), intent(in) ::   p0(:,0:)
    real(kind=dp_t), intent(in) ::   Dh0(:,0:)
    real(kind=dp_t), intent(in) ::   u0_1d(:,0:)

    ! Local variables
    integer                      :: r, n, i
    real(kind=dp_t)              :: mencl
    real(kind=dp_t), allocatable :: m(:,:)
    real(kind=dp_t)              :: Dh0_edge(nlevs_radial,1:nr_fine-1)
    real(kind=dp_t)              :: p0_edge(nlevs_radial,1:nr_fine-1)
    real(kind=dp_t)              :: u0_edge(nlevs_radial,1:nr_fine-1)

    do n=1,nlevs_radial
        do r = 1, nr(n)-1
          Dh0_edge(n,r) = (Dh0(n,r) + Dh0(n,r-1))/2.0d0
          p0_edge(n,r) = (p0(n,r) + p0(n,r-1))/2.0d0
          u0_edge(n,r) = (u0_1d(n,r) + u0_1d(n,r-1))/2.0d0
        enddo
    enddo

    if (spherical .eq. 0) then

          ! we are doing a plane-parallel geometry with a 1/r**2
          ! gravitational acceleration.  The mass is assumed to be
          ! at the origin.  The mass in the computational domain
          ! does not contribute to the gravitational acceleration.
          do n=1,nlevs_radial
             do r = 1, nr(n)-1
                dpdr_edge(n,r) = -(Dh0_edge(n,r) + p0_edge(n,r) * &
                    u0_edge(n,r)) * g / ((c**2 * &
                    (Rr + 2.0d0 * r_edge_loc(n,r)) - 2.0d0 * g) * u0_edge(n,r))
             enddo
          enddo

    else

       dpdr_edge(1,0) = zero
       mencl = ZERO

       do r=1,nr_fine-1

          ! only add to the enclosed mass if the density is
          ! > base_cutoff_density
          if (Dh0(1,r-1) > base_cutoff_density) then
             mencl = mencl + FOUR3RD*M_PI * &
                  (r_edge_loc(1,r) - r_edge_loc(1,r-1)) * &
                  (r_edge_loc(1,r)**2 + &
                   r_edge_loc(1,r)*r_edge_loc(1,r-1) + &
                   r_edge_loc(1,r-1)**2) * Dh0(1,r-1)
          endif

          dpdr_edge(1,r) = -Gconst * mencl / r_edge_loc(1,r)**2

       end do

    end if

  end subroutine make_dpdr_edge

  subroutine make_dpdr_edge_uniform(dpdr_edge_fine,Dh0_fine,p0_fine,u0_1d)

    ! a special version of the make_grav_edge routine that takes a
    ! uniformly-gridded, single level density array at the finest base
    ! state resolution and returns the uniformly-gridded, single-level
    ! gravity at the same resolution

    use bl_constants_module
    use geometry, only: spherical, r_edge_loc, nr_fine, nlevs_radial, nr
    use probin_module, only: base_cutoff_density, &
         g, Rr, c
    use fundamental_constants_module, only: Gconst

    ! compute the base state gravity at the cell edges (grav_edge(1)
    ! is the gravitational acceleration at the left edge of zone 1).
    ! The base state uses 0-based indexing, so grav_edge does too.

    real(kind=dp_t), intent(  out) :: dpdr_edge_fine(0:)
    real(kind=dp_t), intent(in   ) ::      Dh0_fine(0:)
    real(kind=dp_t), intent(in   ) ::      p0_fine(0:)
    real(kind=dp_t), intent(in   ) ::      u0_1d(0:)

    ! Local variables
    integer                      :: r
    real(kind=dp_t)              :: mencl

    if (spherical .eq. 0) then

          ! we are doing a plane-parallel geometry with a 1/r**2
          ! gravitational acceleration.  The mass is assumed to be
          ! at the origin.  The mass in the computational domain
          ! does not contribute to the gravitational acceleration.
          do r = 0, nr(nlevs_radial)-1
             dpdr_edge_fine(r) = -(Dh0_fine(r) + p0_fine(r) * &
                 u0_1d(r)) * g / ((c**2 * &
                 (Rr + 2.0d0 * r_edge_loc(nlevs_radial,r)) - 2.0d0 * g) * u0_1d(r))
          enddo

    else

       dpdr_edge_fine(0) = ZERO
       mencl = ZERO

       do r=1,nr_fine-1

          ! only add to the enclosed mass if the density is
          ! > base_cutoff_density
          if (Dh0_fine(r-1) > base_cutoff_density) then
             mencl = mencl + FOUR3RD*M_PI * &
                  (r_edge_loc(nlevs_radial,r) - r_edge_loc(nlevs_radial,r-1)) * &
                  (r_edge_loc(nlevs_radial,r)**2 + &
                   r_edge_loc(nlevs_radial,r)*r_edge_loc(nlevs_radial,r-1) + &
                   r_edge_loc(nlevs_radial,r-1)**2) * Dh0_fine(r-1)
          endif

          dpdr_edge_fine(r) = -Gconst * mencl / r_edge_loc(nlevs_radial,r)**2

       end do

    end if

  end subroutine make_dpdr_edge_uniform

end module make_grav_module
