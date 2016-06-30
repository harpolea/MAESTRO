! compute beta_0, the coefficient in our constraint equation,
! div{beta_0 U} = beta_0 S

module make_div_coeff_module

  use bl_types

  implicit none

  private

  public :: make_div_coeff

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine make_div_coeff(div_coeff,D0,Dh0,u0,p0,gamma1bar,dpdr_centre)

    use bl_constants_module
    use geometry, only: nr_fine, dr, anelastic_cutoff_coord, r_start_coord, r_end_coord, &
         nr, numdisjointchunks, nlevs_radial, r_cc_loc
    use restrict_base_module
    use probin_module, only: beta_type, use_linear_grav_in_beta, c, g, Rr

    real(kind=dp_t), intent(  out) :: div_coeff(:,0:)
    real(kind=dp_t), intent(in   ) :: D0(:,0:), Dh0(:,0:), u0(:,0:)
    real(kind=dp_t), intent(in   ) :: p0(:,0:), gamma1bar(:,0:)
    real(kind=dp_t), intent(in   ) :: dpdr_centre(:,0:)

    ! local
    integer :: r, n, i, refrat, j
    real(kind=dp_t) :: integral
    real(kind=dp_t) :: beta0_edge(nlevs_radial,0:nr_fine)
    real(kind=dp_t) :: lambda, mu, nu, kappa
    real(kind=dp_t) :: denom, coeff1, coeff2, coeff3
    real(kind=dp_t) :: del,dpls,dmin,slim,sflag
    real(kind=dp_t) :: offset

    !print *, 'D0', D0
    !print *, 'u0', u0
    !print *, 'p0', p0
    !print *, 'dpdr', dpdr_centre

    div_coeff = 0.d0

    if (beta_type .eq. 1) then

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! Compute beta0 on the edges and average to the center
       !
       ! Multilevel Outline:
       !
       ! First, compute beta0 on edges and centers at level 1 only
       ! Obtain the starting value from D0 at the bottom of the domain.
       ! do n=2,nlevs_radial
       !   Compute beta0 on edges and centers at level n
       !   Obtain the starting value of beta0_edge_lo from the coarser grid
       !   if n>1, compare the difference between beta0 at the top of level n to the
       !           corresponding point on level n-1
       !   do i=n-1,1,-1
       !     Offset the centered beta on level i above this point so the total integral
       !      is consistent
       !     Redo the anelastic cutoff part
       !   end do
       ! end do
       ! call restrict_base and fill_ghost_base
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       do n=1,nlevs_radial

          do j=1,numdisjointchunks(n)

             ! Compute beta0 on edges and centers at level n

             if (n .eq. 1) then
                beta0_edge(1,0) = D0(1,0)/u0(1,0)
             else
                ! Obtain the starting value of beta0_edge_lo from the coarser grid
                beta0_edge(n,r_start_coord(n,j)) = beta0_edge(n-1,r_start_coord(n,j)/2)
             end if

             do r=r_start_coord(n,j),r_end_coord(n,j)

                if (r < anelastic_cutoff_coord(n)) then

                   if (r .eq. 0 .or. r .eq. nr(n)-1) then

                      lambda = ZERO
                      mu = ZERO
                      nu = ZERO

                   else

                      ! piecewise linear reconstruction of Dh0,
                      ! gamma1bar, and p0 -- see paper III, appendix C
                      del    = HALF* (Dh0(n,r+1) - Dh0(n,r-1))/dr(n)
                      dpls   = TWO * (Dh0(n,r+1) - Dh0(n,r  ))/dr(n)
                      dmin   = TWO * (Dh0(n,r  ) - Dh0(n,r-1))/dr(n)
                      slim   = min(abs(dpls), abs(dmin))
                      slim   = merge(slim, zero, dpls*dmin.gt.ZERO)
                      sflag  = sign(ONE,del)
                      lambda = sflag*min(slim,abs(del))

                      del   = HALF* (gamma1bar(n,r+1) - gamma1bar(n,r-1))/dr(n)
                      dpls  = TWO * (gamma1bar(n,r+1) - gamma1bar(n,r  ))/dr(n)
                      dmin  = TWO * (gamma1bar(n,r  ) - gamma1bar(n,r-1))/dr(n)
                      slim  = min(abs(dpls), abs(dmin))
                      slim  = merge(slim, zero, dpls*dmin.gt.ZERO)
                      sflag = sign(ONE,del)
                      mu    = sflag*min(slim,abs(del))

                      del   = HALF* (p0(n,r+1) - p0(n,r-1))/dr(n)
                      dpls  = TWO * (p0(n,r+1) - p0(n,r  ))/dr(n)
                      dmin  = TWO * (p0(n,r  ) - p0(n,r-1))/dr(n)
                      slim  = min(abs(dpls), abs(dmin))
                      slim  = merge(slim, zero, dpls*dmin.gt.ZERO)
                      sflag = sign(ONE,del)
                      nu    = sflag*min(slim,abs(del))

                   end if

                   if (nu .eq. ZERO .or. mu .eq. ZERO .or. &
                        (nu*gamma1bar(n,r) - mu*p0(n,r)) .eq. ZERO .or. &
                        ((gamma1bar(n,r) + HALF*mu*dr(n))/ &
                        (gamma1bar(n,r) - HALF*mu*dr(n))) .le. ZERO .or. &
                        ((p0(n,r) + HALF*nu*dr(n))/ &
                        (p0(n,r) - HALF*nu*dr(n))) .le. ZERO) then

                      ! just do piecewise constant integration
                      !integral = !abs(grav_center(n,r))*D0(n,r)*u0(n,r)*dr(n)/(p0(n,r)*gamma1bar(n,r))
                      integral = dpdr_centre(n,r)*dr(n)/(p0(n,r)*gamma1bar(n,r))

                      !print *, 'integral 1', integral

                   else

                      if ( use_linear_grav_in_beta ) then

                         ! also do piecewise linear reconstruction of
                         ! gravity -- not documented in publication yet.
                         del   = HALF* (dpdr_centre(n,r+1) - dpdr_centre(n,r-1))/dr(n)
                         dpls  = TWO * (dpdr_centre(n,r+1) - dpdr_centre(n,r  ))/dr(n)
                         dmin  = TWO * (dpdr_centre(n,r  ) - dpdr_centre(n,r-1))/dr(n)
                         slim  = min(abs(dpls), abs(dmin))
                         slim  = merge(slim, zero, dpls*dmin.gt.ZERO)
                         sflag = sign(ONE,del)
                         kappa = sflag*min(slim,abs(del))

                         denom = nu*gamma1bar(n,r) - mu*p0(n,r)
                         coeff1 = (lambda*gamma1bar(n,r) - mu*Dh0(n,r)) * &
                              (kappa *gamma1bar(n,r) + mu*abs(dpdr_centre(n,r))) / &
                              (mu*mu*denom)
                         coeff2 = (lambda*p0(n,r) - nu*Dh0(n,r))* &
                              (-kappa*p0(n,r) - nu*abs(dpdr_centre(n,r))) / &
                              (nu*nu*denom)
                         coeff3 = kappa*lambda / (mu*nu)

                         integral =  &
                              coeff1*log( (gamma1bar(n,r) + HALF*mu*dr(n))/ &
                                          (gamma1bar(n,r) - HALF*mu*dr(n)) ) + &
                              coeff2*log( (p0(n,r) + HALF*nu*dr(n))/ &
                                          (p0(n,r) - HALF*nu*dr(n)) ) - &
                              coeff3*dr(n)

                        !print *, 'integral 2', integral

                      else

                         ! paper III, equation C2
                         denom = nu*gamma1bar(n,r) - mu*p0(n,r)
                         del   = HALF* (dpdr_centre(n,r+1) - dpdr_centre(n,r-1))/dr(n)
                         dpls  = TWO * (dpdr_centre(n,r+1) - dpdr_centre(n,r  ))/dr(n)
                         dmin  = TWO * (dpdr_centre(n,r  ) - dpdr_centre(n,r-1))/dr(n)
                         slim  = min(abs(dpls), abs(dmin))
                         slim  = merge(slim, zero, dpls*dmin.gt.ZERO)
                         sflag = sign(ONE,del)
                         lambda = sflag*min(slim,abs(del))

                         coeff1 = lambda*gamma1bar(n,r)/mu - dpdr_centre(n,r)
                         coeff2 = lambda*p0(n,r)/nu - dpdr_centre(n,r)

                         !coeff1 = -g*((lambda + nu) *gamma1bar(n,r) - &
                        !    mu*(Dh0(n,r)+p0(n,r)*u0(n,r))) / &
                        !    (TWO*gamma1bar(n,r)*c**2.d0 + mu*(TWO*g - &
                        !    c**2.d0*(TWO*r_cc_loc(n,r)+Rr)))
                         !coeff2 = g*(lambda*p0(n,r)*u0(n,r) - Dh0(n,r)*nu) / &
                        !    (c**2*(TWO*p0(n,r)*u0(n,r) - nu*(TWO* &
                        !    r_cc_loc(n,r)+Rr)) + TWO*g*nu)
                         !coeff3 = g*(c**2.d0*(-TWO*Dh0(n,r)+(lambda+nu) * &
                        !    (TWO*r_cc_loc(n,r)+Rr) - TWO*p0(n,r)*u0(n,r)) - &
                        !    TWO*g*(lambda+nu)) / &
                        !    (TWO*gamma1bar(n,r)*c**2.d0 + mu*(TWO*g - &
                        !    c**2*(TWO*r_cc_loc(n,r)+Rr)) * &
                        !    (c**2*(nu*(TWO*r_cc_loc(n,r)+Rr)- &
                        !    TWO*p0(n,r)*u0(n,r)) - TWO*g*nu))
                        !coeff3 = coeff3 * log((c**2.d0*(Rr+TWO*r_cc_loc(n,r)+dr(n))&
                        !            -TWO*g) / &
                        !           (c**2.d0*(Rr+TWO*r_cc_loc(n,r)-dr(n))&
                        !           -TWO*g))

                         integral = (1.d0/denom)* &
                              (coeff1*log((gamma1bar(n,r) + HALF*mu*dr(n))/ &
                                          (gamma1bar(n,r) - HALF*mu*dr(n))) - &
                               coeff2*log((p0(n,r) + HALF*nu*dr(n))/ &
                                          (p0(n,r) - HALF*nu*dr(n)))) !+ &
                              !coeff3*log((c**2.d0*(Rr+TWO*r_cc_loc(n,r)+dr(n))&
                                !          -TWO*g) / &
                                !         (c**2.d0*(Rr+TWO*r_cc_loc(n,r)-dr(n))&
                                !          -TWO*g))

                        !print *, 'integral 3', integral

                      end if
                   endif
                   ! FIXME: hack to get this to work:
                   if (abs(integral) .gt. 1.d0) then
                       !integral = integral / c
                       integral = 0.d0
                   endif

                   beta0_edge(n,r+1) = beta0_edge(n,r) * exp(integral)
                   div_coeff(n,r) = HALF*(beta0_edge(n,r) + beta0_edge(n,r+1))

                else ! r >= anelastic_cutoff

                   div_coeff(n,r) = div_coeff(n,r-1) * ((D0(n,r)/u0(n,r))/(D0(n,r-1)/u0(n,r-1)))
                   beta0_edge(n,r+1) = 2.d0*div_coeff(n,r) - beta0_edge(n,r)

                endif

             end do

             !! FIXME: previous loop is where div_coeff goes wrong.

             if (n .gt. 1) then

                ! Compare the difference between beta0 at the top of level n to the
                ! corresponding point on level n-1
                offset = beta0_edge(n,r_end_coord(n,j)+1) &
                     - beta0_edge(n-1,(r_end_coord(n,j)+1)/2)

                do i=n-1,1,-1

                   refrat = 2**(n-i)

                   ! Offset the centered beta on level i above this point so the total
                   ! integral is consistent
                   do r=r_end_coord(n,j)/refrat+1,nr(i)
                      div_coeff(i,r) = div_coeff(i,r) + offset
                   end do

                   ! Redo the anelastic cutoff part
                   do r=anelastic_cutoff_coord(i),nr(i)
                      if (D0(i,r-1)/u0(i,r-1) /= ZERO) then
                         div_coeff(i,r) = div_coeff(i,r-1) * (D0(i,r)/u0(i,r))/(D0(i,r-1)/u0(i,r-1))
                      endif
                   end do

                   ! This next piece of coded is needed for the case when the anelastic
                   ! cutoff coordinate lives on level n.  We first average div_coeff from
                   ! level i+1 to level i in the region between the anelastic cutoff and
                   ! the top of grid n.  Then recompute div_coeff at level i above the top
                   ! of grid n.
                   if (r_end_coord(n,j) .ge. anelastic_cutoff_coord(n)) then

                      do r=anelastic_cutoff_coord(i),(r_end_coord(n,j)+1)/refrat-1
                         div_coeff(i,r) = HALF*(div_coeff(i+1,2*r)+div_coeff(i+1,2*r+1))
                      end do

                      do r=(r_end_coord(n,j)+1)/refrat,nr(i)
                         if (D0(i,r-1)/u0(i,r-1) /= ZERO) then
                            div_coeff(i,r) = div_coeff(i,r-1) * (D0(i,r)/u0(i,r))/(D0(i,r-1)/u0(i,r-1))
                         endif
                      end do

                   end if

                end do ! end loop over i=n-1,1,-1

             end if ! end if (n .gt. 1)

          end do ! end loop over disjoint chunks

       end do ! end loop over levels

       ! zero the div_coeff where there is no corresponding full state array
       do n=2,nlevs_radial
          do j=1,numdisjointchunks(n)
             if (j .eq. numdisjointchunks(n)) then
                do r=r_end_coord(n,j)+1,nr(n)-1
                   div_coeff(n,r) = ZERO
                end do
             else
                do r=r_end_coord(n,j)+1,r_start_coord(n,j+1)-1
                   div_coeff(n,r) = ZERO
                end do
             end if
          end do
       end do

    else if (beta_type .eq. 2) then

       ! beta_0 = rho_0
       do n=1,nlevs_radial
          do j=1,numdisjointchunks(n)
             do r=r_start_coord(n,j),r_end_coord(n,j)
                div_coeff(n,r) = D0(n,r)/u0(n,r)
             end do
          end do
       end do

    else if (beta_type .eq. 3) then

       ! beta_0 = 1.d0
       do n=1,nlevs_radial
          do j=1,numdisjointchunks(n)
             do r=r_start_coord(n,j),r_end_coord(n,j)
                div_coeff = 1.d0
             end do
          end do
       end do

    end if

    print *, 'div_coeff', div_coeff

    ! FIXME: hack
    div_coeff(:,:) = div_coeff(:,:) / c

    call restrict_base(div_coeff,.true.)
    call fill_ghost_base(div_coeff,.true.)

  end subroutine make_div_coeff

end module make_div_coeff_module
