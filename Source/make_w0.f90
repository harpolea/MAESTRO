! compute w0 -- the base state velocity.  This is based on the average 
! heating in a layer (Sbar) and the mixing (the eta quantities).  The
! computation of w0 for plane-parallel atmospheres was first described
! in paper II, with modifications due to mixing in paper III.  For 
! spherical geometry, it was first described in paper III.

module make_w0_module

  use bl_types

  implicit none

  private

  public :: make_w0

contains

  subroutine make_w0(nlevs,vel,vel_old,f,Sbar_in,p0,rho0,gam1,eta,dt,dtold)

    use parallel
    use bl_prof_module
    use geometry, only: spherical, nr, dr
    use bl_constants_module
    use probin_module, only: verbose

    integer        , intent(in   ) :: nlevs
    real(kind=dp_t), intent(  out) :: vel(:,0:)
    real(kind=dp_t), intent(in   ) :: vel_old(:,0:)
    real(kind=dp_t), intent(in   ) :: eta(:,0:,:)
    real(kind=dp_t), intent(inout) :: f(:,0:)
    real(kind=dp_t), intent(in   ) :: p0(:,0:),rho0(:,0:),gam1(:,0:)
    real(kind=dp_t), intent(in   ) :: Sbar_in(:,0:)
    real(kind=dp_t), intent(in   ) :: dt,dtold

    integer         :: r,n
    real(kind=dp_t) :: max_vel

    type(bl_prof_timer), save :: bpt

    call build(bpt, "make_w0")

    f(:,:) = ZERO

    do n=1,nlevs
       if (spherical .eq. 0) then
          call make_w0_planar(n,vel(n,:),vel_old(n,:),Sbar_in(n,:),p0(n,:),rho0(n,:), &
                              gam1(n,:),eta(n,:,:),f(n,:),dt,dtold)
       else
          call make_w0_spherical(n,vel(n,:),Sbar_in(n,:),p0(n,:),rho0(n,:),gam1(n,:))
       endif

       max_vel = zero
       do r = 0,nr(n)
          max_vel = max(max_vel, abs(vel(n,r)))
       end do

       if (parallel_IOProcessor() .and. verbose .ge. 1) &
            write(6,*) '... max CFL of w0: ',max_vel * dt / dr(n)
    enddo

    call destroy(bpt)

  end subroutine make_w0

  subroutine make_w0_planar(n,vel,vel_old,Sbar_in,p0,rho0,gam1,eta,f,dt,dtold)

    use geometry, only: nr, dr
    use variables, only: rho_comp
    use bl_constants_module
    use probin_module, only: grav_const, evolve_base_state

    integer        , intent(in   ) :: n
    real(kind=dp_t), intent(  out) :: vel(0:)
    real(kind=dp_t), intent(in   ) :: vel_old(0:)
    real(kind=dp_t), intent(in   ) :: Sbar_in(0:)
    real(kind=dp_t), intent(in   ) :: p0(0:),rho0(0:),gam1(0:),eta(0:,:)
    real(kind=dp_t), intent(inout) ::   f(0:)
    real(kind=dp_t), intent(in   ) :: dt,dtold

    ! Local variables
    integer                      :: r
    real(kind=dp_t), allocatable :: vel_old_cen(:)
    real(kind=dp_t), allocatable :: vel_new_cen(:)
    real(kind=dp_t), allocatable ::   force(:)
    real(kind=dp_t), allocatable ::    edge(:)
    real(kind=dp_t), allocatable :: dpdroverrho(:)
    real(kind=dp_t)              :: eta_avg
    real(kind=dp_t)              :: vel_avg, div_avg, dt_avg

    ! Cell-centered
    allocate(vel_old_cen(0:nr(n)-1))
    allocate(vel_new_cen(0:nr(n)-1))
    allocate(      force(0:nr(n)-1))

    ! Cell-centered
    allocate(dpdroverrho(0:nr(n)-1))
 
    dpdroverrho(      0) = abs(grav_const)
    dpdroverrho(nr(n)-1) = abs(p0(nr(n)-1)-p0(nr(n)-2)) / dr(n) /  rho0(nr(n)-1)
    do r = 1, nr(n)-2
       dpdroverrho(r) = HALF*abs(p0(r+1)-p0(r-1)) / dr(n) /  rho0(r)
       dpdroverrho(r) = min(dpdroverrho(r), abs(grav_const))
    end do

    ! Initialize new velocity to zero.
    vel(0) = ZERO
    
    if (evolve_base_state) then
       do r = 1,nr(n)
          eta_avg = HALF * (eta(r,rho_comp)+eta(r-1,rho_comp))
          vel(r) = vel(r-1) + Sbar_in(r-1) * dr(n) - &
               ( eta_avg * dpdroverrho(r-1) / (gam1(r-1)*p0(r-1)) ) * dr(n)
       end do
    else
       do r = 1,nr(n)
          eta_avg = HALF * (eta(r,rho_comp)+eta(r-1,rho_comp))
          vel(r) = vel(r-1) + Sbar_in(r-1) * dr(n) - &
               ( eta_avg * abs(grav_const) / (gam1(r-1)*p0(r-1)) ) * dr(n)
       end do
    end if

    ! Compute the 1/rho0 grad pi0 term.
    dt_avg = HALF * (dt + dtold)
    do r = 0,nr(n)-1
       vel_old_cen(r) = HALF * (vel_old(r) + vel_old(r+1))
       vel_new_cen(r) = HALF * (vel    (r) + vel    (r+1))
       vel_avg = HALF * (dt *  vel_old_cen(r)           + dtold *  vel_new_cen(r)  ) / dt_avg
       div_avg = HALF * (dt * (vel_old(r+1)-vel_old(r)) + dtold * (vel(r+1)-vel(r))) / dt_avg
       f(r) = (vel_new_cen(r)-vel_old_cen(r)) / dt_avg + &
               vel_avg * div_avg / dr(n)
    end do

    deallocate(vel_old_cen,vel_new_cen,force)

  end subroutine make_w0_planar

  subroutine make_w0_spherical(n,vel,Sbar_in,p0,rho0,gam1)

    use geometry, only: base_cc_loc, nr, base_loedge_loc, dr
    use make_grav_module
    use cell_to_edge_module
    use bl_constants_module
    
    integer        , intent(in   ) :: n
    real(kind=dp_t), intent(  out) :: vel(0:)
    real(kind=dp_t), intent(in   ) :: p0(0:),rho0(0:),gam1(0:)
    real(kind=dp_t), intent(in   ) :: Sbar_in(0:)

    ! Local variables
    integer                      :: r
    real(kind=dp_t), allocatable :: c(:),d(:),e(:),u(:),rhs(:)
    real(kind=dp_t), allocatable :: m(:),grav_edge(:),rho0_edge(:)
    
    ! Cell-centered
    allocate(m(0:nr(n)-1))

    ! Edge-centered
    allocate(c(0:nr(n)),d(0:nr(n)),e(0:nr(n)),rhs(0:nr(n)),u(0:nr(n)))
    allocate(grav_edge(0:nr(n)),rho0_edge(0:nr(n)))

    c(:)   = ZERO
    d(:)   = ZERO
    e(:)   = ZERO
    rhs(:) = ZERO
    u(:)   = ZERO
   
    call make_grav_edge(n,grav_edge,rho0)

    do r = 1,nr(n)
       c(r) = gam1(r-1) * p0(r-1) * base_loedge_loc(n,r-1)**2 / base_cc_loc(n,r-1)**2
       c(r) = c(r) / dr(n)**2
    end do

    call cell_to_edge(n,rho0,rho0_edge)

    do r = 1,nr(n)-1

       d(r) = -( gam1(r-1) * p0(r-1) / base_cc_loc(n,r-1)**2 &
                +gam1(r  ) * p0(r  ) / base_cc_loc(n,r  )**2 ) &
                * (base_loedge_loc(n,r)**2/dr(n)**2) &
                - four * rho0_edge(r) * grav_edge(r) / base_loedge_loc(n,r)
    end do

    do r = 1,nr(n)-1
       rhs(r) = ( gam1(r  )*p0(r  )*Sbar_in(r) - gam1(r-1)*p0(r-1)*Sbar_in(r-1) ) 
       rhs(r) = rhs(r) / dr(n)
    end do

    do r = 0,nr(n)-1
       e(r) = gam1(r) * p0(r) * base_loedge_loc(n,r+1)**2 / base_cc_loc(n,r)**2
       e(r) = e(r) / dr(n)**2
    end do

    ! Lower boundary
       d(0) = one
       e(0) = zero
     rhs(0) = zero

    ! Upper boundary
       c(nr(n)) = zero
       d(nr(n)) = one
     rhs(nr(n)) = zero

    ! Call the tridiagonal solver
    call tridiag(c, d, e, rhs, u, nr(n)+1)

    do r = 0,nr(n)
       vel(r) = u(r)
    end do

    deallocate(c,d,e,rhs,u)
    deallocate(m,grav_edge,rho0_edge)

  end subroutine make_w0_spherical

  subroutine tridiag(a,b,c,r,u,n)

    use bl_error_module

    real(kind=dp_t), intent(in   ) :: a(:), b(:), c(:), r(:)
    real(kind=dp_t), intent(  out) :: u(:)
    integer, intent(in)            :: n

    real(kind=dp_t)              :: bet
    real(kind=dp_t), allocatable :: gam(:)
    integer                      :: j

    allocate(gam(n))
    
    if ( b(1) .eq. 0 ) call bl_error('tridiag: CANT HAVE B(1) = ZERO')
    
    bet = b(1)
    u(1) = r(1)/bet
    
    do j = 2,n
       gam(j) = c(j-1)/bet
       bet = b(j) - a(j)*gam(j)
       if ( bet .eq. 0 ) call bl_error('tridiag: TRIDIAG FAILED')
       u(j) = (r(j)-a(j)*u(j-1))/bet
    end do
    
    do j = n-1,1,-1
       u(j) = u(j) - gam(j+1)*u(j+1)
    end do
    
  end subroutine tridiag
  
end module make_w0_module
