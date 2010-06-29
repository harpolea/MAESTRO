! Adapted from the initdata.f90 file in the xrb/ directory

module init_module

  use bl_types
  use bl_constants_module
  use bc_module
  use multifab_physbc_module
  use define_bc_module
  use multifab_module
  use eos_module
  use variables
  use network
  use geometry, only: nr, spherical, dm, nlevs
  use ml_layout_module
  use ml_restriction_module
  use multifab_fill_ghost_module

  implicit none

  real(dp_t), save :: pert_height

  private
  public :: initscalardata, initscalardata_on_level, initveldata

contains

  subroutine initscalardata(s,s0_init,p0_init,dx,bc,mla)

    use probin_module, only: prob_lo, prob_hi, perturb_model, nova_pert_height

    type(multifab) , intent(inout) :: s(:)
    real(kind=dp_t), intent(in   ) :: s0_init(:,0:,:)
    real(kind=dp_t), intent(in   ) :: p0_init(:,0:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    type(bc_level) , intent(in   ) :: bc(:)
    type(ml_layout), intent(inout) :: mla

    real(kind=dp_t), pointer:: sop(:,:,:,:)
    integer :: lo(dm),hi(dm),ng
    integer :: i,n,r


    ng = s(1)%ng

    do n=1,nlevs

       do i = 1, s(n)%nboxes

          if ( multifab_remote(s(n),i) ) cycle

          sop => dataptr(s(n),i)
          lo =  lwb(get_box(s(n),i))
          hi =  upb(get_box(s(n),i))
          select case (dm)
          case (2)

             ! make sure our perturbation is within the domain
             if (perturb_model .and. &
                  (nova_pert_height < prob_lo(2) .or. &
                   nova_pert_height > prob_hi(2))) &
                   call bl_error("perturbation outside of the domain")

             call initscalardata_2d(sop(:,:,1,:), lo, hi, ng, dx(n,:), &
                  s0_init(n,:,:), p0_init(n,:))
          case (3)

             ! make sure our perturbation is within the domain
             if (perturb_model .and. &
                  (nova_pert_height < prob_lo(3) .or. &
                   nova_pert_height > prob_hi(3))) &
                   call bl_error("perturbation outside of the domain")

             call initscalardata_3d(sop(:,:,:,:), lo, hi, ng, dx(n,:), &
                  s0_init(n,:,:), p0_init(n,:))
          end select

       enddo

    enddo

    if (nlevs .eq. 1) then

       ! fill ghost cells for two adjacent grids at the same level
       ! this includes periodic domain boundary ghost cells
       call multifab_fill_boundary(s(nlevs))

       ! fill non-periodic domain boundary ghost cells
       call multifab_physbc(s(nlevs),rho_comp,dm+rho_comp,nscal,bc(nlevs))

    else

       ! the loop over nlevs must count backwards to make sure the finer grids
       ! are done first
       do n=nlevs,2,-1

          ! set level n-1 data to be the average of the level n data covering 
          ! it
          call ml_cc_restriction(s(n-1),s(n),mla%mba%rr(n-1,:))

          ! fill level n ghost cells using interpolation from level n-1 data
          ! note that multifab_fill_boundary and multifab_physbc are called for
          ! both levels n-1 and n
          call multifab_fill_ghost_cells(s(n),s(n-1),ng,mla%mba%rr(n-1,:), &
                                         bc(n-1),bc(n),rho_comp,dm+rho_comp, &
                                         nscal, fill_crse_input=.false.)

       enddo

    end if

  end subroutine initscalardata



  subroutine initscalardata_on_level(n,s,s0_init,p0_init,dx,bc)

    integer        , intent(in   ) :: n
    type(multifab) , intent(inout) :: s
    real(kind=dp_t), intent(in   ) :: s0_init(0:,:)
    real(kind=dp_t), intent(in   ) :: p0_init(0:)
    real(kind=dp_t), intent(in   ) :: dx(:)
    type(bc_level) , intent(in   ) :: bc

    ! local
    integer                    :: ng,i,r
    integer                    :: lo(dm),hi(dm)
    real(kind=dp_t), pointer   :: sop(:,:,:,:)


    ng = s%ng

    do i = 1, s%nboxes

       if ( multifab_remote(s,i) ) cycle

       sop => dataptr(s,i)
       lo =  lwb(get_box(s,i))
       hi =  upb(get_box(s,i))
       select case (dm)
       case (2)
          call initscalardata_2d(sop(:,:,1,:), lo, hi, ng, dx, s0_init, &
               p0_init)
       case (3)
          call initscalardata_3d(sop(:,:,:,:), lo, hi, ng, dx, s0_init, &
               p0_init)
       end select

    end do

    call multifab_fill_boundary(s)

    call multifab_physbc(s,rho_comp,dm+rho_comp,nscal,bc)

  end subroutine initscalardata_on_level



  subroutine initscalardata_2d(s,lo,hi,ng,dx,s0_init,p0_init)

    use probin_module, only: prob_lo, prob_hi, perturb_model, nova_pert_height

    integer           , intent(in   ) :: lo(:),hi(:),ng
    real (kind = dp_t), intent(inout) :: s(lo(1)-ng:,lo(2)-ng:,:)  
    real (kind = dp_t), intent(in   ) :: dx(:)
    real(kind=dp_t)   , intent(in   ) :: s0_init(0:,:)
    real(kind=dp_t)   , intent(in   ) :: p0_init(0:)

    ! Local variables
    integer         :: i,j
    real(kind=dp_t) :: x,y
    real(kind=dp_t) :: dens_pert, rhoh_pert, temp_pert
    real(kind=dp_t) :: rhoX_pert(nspec), trac_pert(ntrac)

    real(kind=dp_t) :: xcen, ycen, dist

    ! initial the domain with the base state
    s = ZERO

    ! initialize the scalars
    do j = lo(2), hi(2)

       do i = lo(1), hi(1)

          s(i,j,rho_comp)  = s0_init(j,rho_comp)
          s(i,j,rhoh_comp) = s0_init(j,rhoh_comp)
          s(i,j,temp_comp) = s0_init(j,temp_comp)
          s(i,j,spec_comp:spec_comp+nspec-1) = &
               s0_init(j,spec_comp:spec_comp+nspec-1)
          s(i,j,trac_comp:trac_comp+ntrac-1) = &
               s0_init(j,trac_comp:trac_comp+ntrac-1)

       enddo

    enddo
    
    ! add an optional perturbation
    if (perturb_model) then

       xcen = (prob_lo(1) + prob_hi(1)) / TWO
       ycen = nova_pert_height

       do j = lo(2), hi(2)
          y = prob_lo(2) + (dble(j)+HALF) * dx(2)
          
          do i = lo(1), hi(1)
             x = prob_lo(1) + (dble(i)+HALF) * dx(1)
          
             dist = sqrt((x-xcen)**2 + (y-ycen)**2)

             call perturb(dist, p0_init(j), s0_init(j,:), &
                          dens_pert, rhoh_pert, rhoX_pert, temp_pert, &
                          trac_pert)

             s(i,j,rho_comp) = dens_pert
             s(i,j,rhoh_comp) = rhoh_pert
             s(i,j,temp_comp) = temp_pert
             s(i,j,spec_comp:spec_comp+nspec-1) = rhoX_pert(1:)
             s(i,j,trac_comp:trac_comp+ntrac-1) = trac_pert(:)

          enddo

       enddo

    endif
    
  end subroutine initscalardata_2d




  subroutine initscalardata_3d(s,lo,hi,ng,dx,s0_init,p0_init)

    use probin_module, only: prob_lo, prob_hi, perturb_model, nova_pert_height
    
    integer           , intent(in   ) :: lo(:),hi(:),ng
    real (kind = dp_t), intent(inout) :: s(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)  
    real (kind = dp_t), intent(in   ) :: dx(:)
    real(kind=dp_t)   , intent(in   ) :: s0_init(0:,:)
    real(kind=dp_t)   , intent(in   ) :: p0_init(0:)

    ! local
    integer         :: i,j,k
    real(kind=dp_t) :: x, y, z
    real(kind=dp_t) :: dens_pert, rhoh_pert, temp_pert
    real(kind=dp_t) :: rhoX_pert(nspec), trac_pert(ntrac)

    real(kind=dp_t)            :: xcen, ycen, zcen, dist

    ! initial the domain with the base state
    s = ZERO
  
    if (spherical .eq. 1) then

       call bl_error('Error: initdata does not handle the spherical case')

    else 

       ! initialize the scalars
       do k = lo(3), hi(3)

          do j = lo(2), hi(2)

             do i = lo(1), hi(1)
                
                s(i,j,k,rho_comp)  = s0_init(k,rho_comp)
                s(i,j,k,rhoh_comp) = s0_init(k,rhoh_comp)
                s(i,j,k,temp_comp) = s0_init(k,temp_comp)
                s(i,j,k,spec_comp:spec_comp+nspec-1) = &
                     s0_init(k,spec_comp:spec_comp+nspec-1)
                s(i,j,k,trac_comp:trac_comp+ntrac-1) = &
                     s0_init(k,trac_comp:trac_comp+ntrac-1)

             enddo

          enddo

       enddo
       
       if (perturb_model) then

          xcen = (prob_lo(1) + prob_hi(1)) / TWO
          ycen = (prob_lo(2) + prob_hi(2)) / TWO
          zcen = nova_pert_height

          ! add an optional perturbation
          do k = lo(3), hi(3)
             z = prob_lo(3) + (dble(k)+HALF) * dx(3)
             
             do j = lo(2), hi(2)
                y = prob_lo(2) + (dble(j)+HALF) * dx(2)
                
                do i = lo(1), hi(1)
                   x = prob_lo(1) + (dble(i)+HALF) * dx(1)

                   dist = sqrt((x-xcen)**2 + (y-ycen)**2 + (z-zcen)**2)
                   
                   call perturb(dist, p0_init(k), s0_init(k,:), &
                                dens_pert, rhoh_pert, rhoX_pert, temp_pert, &
                                trac_pert)

                   s(i,j,k,rho_comp) = dens_pert
                   s(i,j,k,rhoh_comp) = rhoh_pert
                   s(i,j,k,temp_comp) = temp_pert
                   s(i,j,k,spec_comp:spec_comp+nspec-1) = rhoX_pert(:)
                   s(i,j,k,trac_comp:trac_comp+ntrac-1) = trac_pert(:)

                enddo

             enddo

          enddo

       endif

    end if
    
  end subroutine initscalardata_3d




  subroutine initveldata(u,s0_init,p0_init,dx,bc,mla)
    
    use probin_module, only: prob_lo, prob_hi, num_vortices

    type(multifab) , intent(inout) :: u(:)
    real(kind=dp_t), intent(in   ) :: s0_init(:,0:,:)
    real(kind=dp_t), intent(in   ) :: p0_init(:,0:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    type(bc_level) , intent(in   ) :: bc(:)
    type(ml_layout), intent(inout) :: mla

    real(kind=dp_t), pointer:: uop(:,:,:,:)
    integer :: lo(dm),hi(dm),ng
    integer :: i,n

    real(kind=dp_t) :: xloc_vortices(num_vortices)
    real(kind=dp_t) :: offset

    ! for now, this is calculated even if we don't use the velocity field in
    ! the initveldata_2d routine below
    offset = (prob_hi(1) - prob_lo(1)) / (num_vortices + 1)

    do i = 1, num_vortices
       xloc_vortices(i) = dble(i) * offset
    enddo
    
    
    ng = u(1)%ng

    do n=1,nlevs
       do i = 1, u(n)%nboxes
          if ( multifab_remote(u(n),i) ) cycle
          uop => dataptr(u(n),i)
          lo =  lwb(get_box(u(n),i))
          hi =  upb(get_box(u(n),i))
          select case (dm)
          case (2)
             call initveldata_2d(uop(:,:,1,:), lo, hi, ng, dx(n,:), &
                                 s0_init(n,:,:), p0_init(n,:), xloc_vortices)
          ! 3d doesn't currently have vortice information coded !
          case (3) 
             call initveldata_3d(uop(:,:,:,:), lo, hi, ng, dx(n,:), &
                                 s0_init(n,:,:), p0_init(n,:))
          end select
       end do
    enddo

    if (nlevs .eq. 1) then

       ! fill ghost cells for two adjacent grids at the same level
       ! this includes periodic domain boundary ghost cells
       call multifab_fill_boundary(u(nlevs))

       ! fill non-periodic domain boundary ghost cells
       call multifab_physbc(u(nlevs),1,1,dm,bc(nlevs))
    else
    
       ! the loop over nlevs must count backwards to make sure the finer grids
       ! are done first
       do n=nlevs,2,-1

          ! set level n-1 data to be the average of the level n data covering 
          ! it
          call ml_cc_restriction(u(n-1),u(n),mla%mba%rr(n-1,:))
          
          ! fill level n ghost cells using interpolation from level n-1 data
          ! note that multifab_fill_boundary and multifab_physbc are called for
          ! both levels n-1 and n
          call multifab_fill_ghost_cells(u(n),u(n-1),ng,mla%mba%rr(n-1,:), &
                                         bc(n-1),bc(n),1,1,dm, &
                                         fill_crse_input=.false.)

       enddo

    end if

  end subroutine initveldata




  subroutine initveldata_2d(u,lo,hi,ng,dx,s0_init,p0_init,xloc_vortices)

    use probin_module, only: prob_lo, apply_vel_field, velpert_scale, &
                             velpert_amplitude, velpert_height_loc

    integer           , intent(in   ) :: lo(:),hi(:),ng
    real (kind = dp_t), intent(  out) :: u(lo(1)-ng:,lo(2)-ng:,:)  
    real (kind = dp_t), intent(in   ) :: dx(:)
    real(kind=dp_t)   , intent(in   ) :: s0_init(0:,:)
    real(kind=dp_t)   , intent(in   ) :: p0_init(0:)
    real(kind=dp_t)   , intent(in   ) :: xloc_vortices(:)

    ! local variables
    real(kind=dp_t) :: xloc(2), upert(2)
    integer :: i, j, vortex
    real(kind=dp_t) :: xdist, ydist, r


    u = ZERO

    if (apply_vel_field) then

       do j = lo(2), hi(2)

          xloc(2) = prob_lo(2) + (dble(j)+HALF)*dx(2)

          ydist = xloc(2) - velpert_height_loc
          
          do i = lo(1), hi(1)

             upert = ZERO

             xloc(1) = prob_lo(1) + (dble(i)+HALF)*dx(1)

             ! loop over each vortex
             do vortex = 1, size(xloc_vortices, dim=1)

                xdist = xloc(1) - xloc_vortices(vortex)

                r = xdist**2 + ydist**2
                r = sqrt(r)

                ! e.g. Calder et al. ApJSS 143, 201-229 (2002)
                ! we set things up so that every other vortex has the same
                ! orientation
                upert(1) = upert(1) - ydist * &
                     velpert_amplitude * exp( -r**2/(TWO*velpert_scale)) &
                     * (-ONE)**vortex

                upert(2) = upert(2) + xdist * &
                     velpert_amplitude * exp(-r**2/(TWO*velpert_scale)) &
                     * (-ONE)**vortex
             enddo

             u(i,j,:) = u(i,j,:) + upert(:)

          enddo

       enddo
       
                
    endif

  end subroutine initveldata_2d




  subroutine initveldata_3d(u,lo,hi,ng,dx,s0_init,p0_init)

    use probin_module, only: apply_vel_field

    integer           , intent(in   ) :: lo(:), hi(:), ng
    real (kind = dp_t), intent(  out) :: u(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)  
    real (kind = dp_t), intent(in   ) :: dx(:)
    real(kind=dp_t)   , intent(in   ) :: s0_init(0:,:)
    real(kind=dp_t)   , intent(in   ) :: p0_init(0:)

    ! initial the velocity
    u = ZERO

    ! we don't support the vortices in 3d yet...
    if (apply_vel_field) call bl_error("apply_vel_field not supported for 3d")

    
  end subroutine initveldata_3d




  subroutine perturb(distance, p0_init, s0_init,  &
                     dens_pert, rhoh_pert, rhoX_pert, temp_pert, trac_pert)

    use probin_module, ONLY: nova_pert_size, nova_pert_factor, nova_pert_type
    ! apply an optional perturbation to the initial temperature field
    ! to see some bubbles

    real(kind=dp_t), intent(in ) :: distance
    real(kind=dp_t), intent(in ) :: p0_init, s0_init(:)
    real(kind=dp_t), intent(out) :: dens_pert, rhoh_pert, temp_pert
    real(kind=dp_t), intent(out) :: rhoX_pert(:)
    real(kind=dp_t), intent(out) :: trac_pert(:)

    real(kind=dp_t) :: temp,dens,t0,d0,rad_pert

    integer, parameter :: perturb_temp = 1, perturb_dens = 2
    integer :: eos_input_flag

    rad_pert = -nova_pert_size**2 / (FOUR*log(HALF))

    select case (nova_pert_type)
    case(perturb_temp)

       t0 = s0_init(temp_comp)

       temp = t0 * (ONE + nova_pert_factor * dexp(-distance**2 / rad_pert) )

       dens = s0_init(rho_comp)

       eos_input_flag = eos_input_tp

    case(perturb_dens)
          
       d0 = s0_init(rho_comp)
       
       dens = d0 * (ONE + nova_pert_factor * dexp(-distance**2 / rad_pert) )
       
       temp = s0_init(temp_comp)

       eos_input_flag = eos_input_rp

    end select

    temp_eos(1) = temp
    p_eos(1) = p0_init
    den_eos(1) = dens
    xn_eos(1,:) = s0_init(spec_comp:spec_comp+nspec-1)/s0_init(rho_comp)

    call eos(eos_input_flag, den_eos, temp_eos, &
             npts, &
             xn_eos, &
             p_eos, h_eos, e_eos, &
             cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
             dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
             dpdX_eos, dhdX_eos, &
             gam1_eos, cs_eos, s_eos, &
             dsdt_eos, dsdr_eos, &
             .false.)

    dens_pert = den_eos(1)
    rhoh_pert = den_eos(1)*h_eos(1)
    rhoX_pert(:) = dens_pert*xn_eos(1,:)

    temp_pert = temp_eos(1)
    
    trac_pert(:) = ZERO

  end subroutine perturb




end module init_module