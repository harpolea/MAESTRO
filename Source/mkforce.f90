module mk_vel_force_module

  ! Compute the force that appears in the velocity (or momentum)
  ! equations.  This is used both when predicting the interface
  ! states and in the final, conservative update.

  ! for the final conservative update of the velocity, we need to
  ! time-center the Coriolis term ( -2 omega x U ), which means we
  ! should use umac.  This is selected by setting is_final_update = T

  use multifab_module
  use define_bc_module
  use ml_layout_module

  implicit none

  private
  public :: mk_vel_force

contains

  subroutine mk_vel_force(vel_force,is_final_update, &
                          uold,umac,w0,w0mac,gpi,s,index_rho,index_Dh,normal, &
                          D0,Dh0,grav,dx,w0_force,w0_force_cart,the_bc_level,mla, &
                          do_add_utilde_force,u0,chrls,gam)

    ! index_rho refers to the index into s where the density lives.
    ! Usually s will be the full state array, and index_rho would
    ! be rho_comp, but sometimes we pass in only a single-variable
    ! multifab array, so index_rho may be different.

    use bl_prof_module
    use geometry, only: spherical, nr_fine, dr
    use bl_constants_module
    use ml_restrict_fill_module
    use probin_module, only: evolve_base_state
    use fill_3d_module, only : put_1d_array_on_cart
    use variables, only : foextrap_comp

    type(multifab) , intent(inout) :: vel_force(:)
    logical        , intent(in   ) :: is_final_update
    type(multifab) , intent(in   ) :: uold(:)
    type(multifab) , intent(in   ) :: normal(:)
    type(multifab) , intent(in   ) :: umac(:,:)
    real(kind=dp_t), intent(in   ) :: w0(:,0:)
    type(multifab) , intent(in   ) :: w0mac(:,:)
    type(multifab) , intent(in   ) :: gpi(:)
    type(multifab) , intent(in   ) :: s(:)
    integer                        :: index_rho
    integer                        :: index_Dh
    real(kind=dp_t), intent(in   ) :: D0(:,0:)
    real(kind=dp_t), intent(in   ) :: Dh0(:,0:)
    real(kind=dp_t), intent(in   ) :: grav(:,0:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    real(kind=dp_t), intent(in   ) :: w0_force(:,0:)
    type(multifab) , intent(in   ) :: w0_force_cart(:)
    type(bc_level) , intent(in   ) :: the_bc_level(:)
    type(ml_layout), intent(inout) :: mla
    logical        , intent(in   ) :: do_add_utilde_force
    type(multifab) , intent(in   ) :: u0(:)
    real(kind=dp_t), intent(in   ) :: chrls(:,0:,0:,0:,0:,0:,0:)
    type(multifab) , intent(in   ) :: gam(:)

    ! Local variables
    real(kind=dp_t), pointer ::  uop(:,:,:,:)
    real(kind=dp_t), pointer ::  ump(:,:,:,:)
    real(kind=dp_t), pointer ::  vmp(:,:,:,:)
    real(kind=dp_t), pointer ::  wmp(:,:,:,:)
    real(kind=dp_t), pointer ::   np(:,:,:,:)
    real(kind=dp_t), pointer :: w0cp(:,:,:,:)
    real(kind=dp_t), pointer :: w0xp(:,:,:,:)
    real(kind=dp_t), pointer :: w0yp(:,:,:,:)
    real(kind=dp_t), pointer ::  gpp(:,:,:,:)
    real(kind=dp_t), pointer ::   fp(:,:,:,:)
    real(kind=dp_t), pointer ::   rp(:,:,:,:)
    real(kind=dp_t), pointer ::  w0p(:,:,:,:)
    real(kind=dp_t), pointer :: gw0p(:,:,:,:)
    real(kind=dp_t), pointer ::  u0p(:,:,:,:)
    real(kind=dp_t), pointer ::  gamp(:,:,:,:)

    ! stuff for spherical only
    real(kind=dp_t) :: gradw0_rad(1,0:nr_fine-1)
    type(multifab)  :: gradw0_cart(mla%nlevel)

    integer                  :: i,r,lo(mla%dim),hi(mla%dim),dm,nlevs
    integer                  :: ng_s,ng_f,ng_gp,n,ng_uo,ng_um, ng_n, ng_u0

    type(multifab) :: w0_cart(mla%nlevel)
    integer :: ng_wc, ng_wm, ng_w, ng_gw

    type(bl_prof_timer), save :: bpt

    call build(bpt, "mk_vel_force")

    dm = mla%dim
    nlevs = mla%nlevel

    ng_s  = nghost(s(1))
    ng_f  = nghost(vel_force(1))
    ng_gp = nghost(gpi(1))


    ! put w0 and gradw0 on cell centers
    if (spherical .eq. 1) then
       do n=1,nlevs
          ! w0_cart will contain the cell-centered Cartesian components of w0,
          ! for use in computing the Coriolis term in the prediction
          ! w0mac is passed in and is used to compute the Coriolis term
          ! in the cell update
          call build(w0_cart(n),mla%la(n),dm,0)
          call setval(w0_cart(n), ZERO, all=.true.)

          call build(gradw0_cart(n),get_layout(vel_force(n)),1,1)
          call setval(gradw0_cart(n), ZERO, all=.true.)
       enddo

       if (evolve_base_state) then
          ! fill the all dm components of the cell-centered w0_cart
          call put_1d_array_on_cart(w0,w0_cart,1,.true.,.true.,dx, &
                                    the_bc_level,mla)

          if (do_add_utilde_force) then
             !$OMP PARALLEL DO PRIVATE(r)
             do r=0,nr_fine-1
                gradw0_rad(1,r) = (w0(1,r+1) - w0(1,r)) / dr(1)
             enddo
             !$OMP END PARALLEL DO

             call put_1d_array_on_cart(gradw0_rad,gradw0_cart,foextrap_comp, &
                                       .false.,.false.,dx,the_bc_level,mla)
          endif

       end if

    endif

    do n=1,nlevs
       do i=1,nfabs(s(n))
          fp  => dataptr(vel_force(n),i)
          gpp => dataptr(gpi(n),i)
          rp  => dataptr(s(n),i)

          lo = lwb(get_box(s(n),i))
          hi = upb(get_box(s(n),i))

          ump => dataptr(umac(n,1),i)
          ng_um = nghost(umac(1,1))

          select case (dm)
          case (1)
             call mk_vel_force_1d(fp(:,1,1,1),ng_f,gpp(:,1,1,1),ng_gp, &
                                  rp(:,1,1,index_rho),ng_s, &
                                  ump(:,1,1,1), ng_um, &
                                  D0(n,:),grav(n,:),w0(n,:),w0_force(n,:),lo,hi,n, &
                                  do_add_utilde_force)

          case (2)
             vmp => dataptr(umac(n,2),i)
             call mk_vel_force_2d(fp(:,:,1,:),ng_f,gpp(:,:,1,:),ng_gp, &
                                  rp(:,:,1,index_rho),ng_s, &
                                  vmp(:,:,1,1), ng_um, &
                                  D0(n,:),grav(n,:),w0(n,:),w0_force(n,:),lo,hi,n, &
                                  do_add_utilde_force)

          case (3)
             uop => dataptr(uold(n),i)
             vmp => dataptr(umac(n,2),i)
             wmp => dataptr(umac(n,3),i)

             ng_uo = nghost(uold(1))
             ng_u0 = nghost(u0(1))

             u0p => dataptr(u0(n),i)
             gamp => dataptr(gam(n),i)

             if (spherical .eq. 1) then
                w0cp  => dataptr(w0_cart(n), i)
                w0xp  => dataptr(w0mac(n,1),i)
                w0yp  => dataptr(w0mac(n,2),i)
                w0p   => dataptr(w0_force_cart(n), i)
                np    => dataptr(normal(n),i)
                gw0p   => dataptr(gradw0_cart(n),i)

                ng_wm = nghost(w0mac(1,1))
                ng_wc = nghost(w0_cart(1))
                ng_w  = nghost(w0_force_cart(1))
                ng_gw = nghost(gradw0_cart(1))
                ng_n  = nghost(normal(1))

                call mk_vel_force_3d_sphr(fp(:,:,:,:),ng_f,is_final_update, &
                                          uop(:,:,:,:),ng_uo,np(:,:,:,:),ng_n, &
                                          ump(:,:,:,1),vmp(:,:,:,1),wmp(:,:,:,1),ng_um, &
                                          w0cp(:,:,:,:),ng_wc,gw0p(:,:,:,1),ng_gw, &
                                          w0xp(:,:,:,1),w0yp(:,:,:,1),ng_wm, &
                                          gpp(:,:,:,:),ng_gp,rp(:,:,:,index_rho),ng_s, &
                                          D0(1,:),grav(1,:),w0p(:,:,:,:),ng_w,lo,hi,dx(n,:), &
                                          do_add_utilde_force)

             else
                call mk_vel_force_3d_cart(fp(:,:,:,:),ng_f,is_final_update, &
                                          uop(:,:,:,:),ng_uo, &
                                          ump(:,:,:,1),vmp(:,:,:,1),wmp(:,:,:,1),ng_um, &
                                          w0(n,:), &
                                          gpp(:,:,:,:),ng_gp,rp(:,:,:,index_rho),rp(:,:,:,index_Dh),ng_s, &
                                          D0(n,:),Dh0(n,:),grav(n,:),w0_force(n,:),lo,hi,n, &
                                          do_add_utilde_force,u0p(:,:,:,1), &
                                          ng_u0,chrls(n,:,:,:,:,:,:),gamp(:,:,:,:))
             end if
          end select
       end do
    enddo

    if (spherical .eq. 1) then
       do n=1,nlevs
          call destroy(w0_cart(n))
          call destroy(gradw0_cart(n))
       end do
    end if

    ! restrict data and fill all ghost cells
    call ml_restrict_and_fill(nlevs,vel_force,mla%mba%rr,the_bc_level, &
                              icomp=1, &
                              bcomp=1, &
                              nc=dm, &
                              ng=vel_force(1)%ng)

    call destroy(bpt)

  end subroutine mk_vel_force

  subroutine mk_vel_force_1d(vel_force,ng_f,gpi,ng_gp, &
                             D,ng_s, &
                             umac,ng_um, &
                             D0,grav,w0,w0_force,lo,hi,n, &
                             do_add_utilde_force)

    use geometry, only: nr, dr
    use bl_constants_module
    use probin_module, only: base_cutoff_density, buoyancy_cutoff_factor

    integer        , intent(in   ) ::  lo(:),hi(:),ng_f,ng_gp,ng_s,ng_um
    real(kind=dp_t), intent(inout) :: vel_force(lo(1)-ng_f :)
    real(kind=dp_t), intent(in   ) ::     gpi(lo(1)-ng_gp:)
    real(kind=dp_t), intent(in   ) ::     D(lo(1)-ng_s :)
    real(kind=dp_t), intent(in   ) ::    umac(lo(1)-ng_um:)
    real(kind=dp_t), intent(in   ) :: D0(0:)
    real(kind=dp_t), intent(in   ) :: grav(0:)
    real(kind=dp_t), intent(in   ) :: w0(0:), w0_force(0:)
    logical        , intent(in   ) :: do_add_utilde_force
    integer        , intent(in   ) :: n
    integer         :: i
    real(kind=dp_t) :: Dpert

    vel_force = ZERO

    do i = lo(1),hi(1)

       Dpert = D(i) - D0(i)

       ! cutoff the buoyancy term if we are outside of the star
       if (D(i) .lt. buoyancy_cutoff_factor*base_cutoff_density) then
          Dpert = 0.d0
       end if

       ! note: if use_alt_energy_fix = T, then gphi is already weighted
       ! by beta_0
       vel_force(i) =  Dpert / D(i) * grav(i) - gpi(i) / D(i) - w0_force(i)

    end do

    if (do_add_utilde_force) then
       do i=lo(1),hi(1)

          if (i .le. -1) then
             ! do not modify force since dw0/dr=0
          else if (i .ge. nr(n)) then
             ! do not modify force since dw0/dr=0
          else
             vel_force(i) = vel_force(i) &
                  - (umac(i+1)+umac(i))*(w0(i+1)-w0(i)) / (TWO*dr(n))
          end if

       enddo
    endif


  end subroutine mk_vel_force_1d

  subroutine mk_vel_force_2d(vel_force,ng_f,gpi,ng_gp, &
                             D,ng_s, &
                             vmac, ng_um, &
                             D0,grav,w0,w0_force,lo,hi,n, &
                             do_add_utilde_force)

    use geometry, only: nr, dr
    use bl_constants_module
    use probin_module, only: base_cutoff_density, buoyancy_cutoff_factor

    integer        , intent(in   ) ::  lo(:),hi(:),ng_f,ng_gp,ng_s,ng_um, n
    real(kind=dp_t), intent(inout) :: vel_force(lo(1)-ng_f :,lo(2)-ng_f :,:)
    real(kind=dp_t), intent(in   ) ::     gpi(lo(1)-ng_gp:,lo(2)-ng_gp:,:)
    real(kind=dp_t), intent(in   ) ::     D(lo(1)-ng_s :,lo(2)-ng_s :)
    real(kind=dp_t), intent(in   ) ::    vmac(lo(1)-ng_um:,lo(2)-ng_um:)
    real(kind=dp_t), intent(in   ) :: D0(0:)
    real(kind=dp_t), intent(in   ) :: grav(0:)
    real(kind=dp_t), intent(in   ) :: w0(0:),w0_force(0:)
    logical        , intent(in   ) :: do_add_utilde_force

    integer         :: i,j
    real(kind=dp_t) :: Dpert

    vel_force = ZERO

    do j = lo(2),hi(2)
       do i = lo(1),hi(1)

          Dpert = D(i,j) - D0(j)

          ! cutoff the buoyancy term if we are outside of the star
          if (D(i,j) .lt. buoyancy_cutoff_factor*base_cutoff_density) then
             Dpert = 0.d0
          end if

          ! note: if use_alt_energy_fix = T, then gphi is already weighted
          ! by beta_0
          vel_force(i,j,1) = - gpi(i,j,1) / D(i,j)
          vel_force(i,j,2) =  Dpert / D(i,j) * grav(j) &
               - gpi(i,j,2) / D(i,j) - w0_force(j)
       end do
    end do

    if (do_add_utilde_force) then

       do j=lo(2),hi(2)
          do i=lo(1),hi(1)

             if (j .le. -1) then
                ! do not modify force since dw0/dr=0
             else if (j .ge. nr(n)) then
                ! do not modify force since dw0/dr=0
             else
                vel_force(i,j,2) = vel_force(i,j,2) &
                     - (vmac(i,j+1)+vmac(i,j))*(w0(j+1)-w0(j)) / (TWO*dr(n))
             end if

          end do
       end do

    endif


  end subroutine mk_vel_force_2d

  subroutine mk_vel_force_3d_cart(vel_force,ng_f,is_final_update, &
                                  uold,ng_uo, &
                                  umac,vmac,wmac,ng_um, &
                                  w0, &
                                  gpi,ng_gp,D,Dh,ng_s, &
                                  D0,Dh0,grav,w0_force,lo,hi,n, &
                                  do_add_utilde_force,u0,ng_u0,chrls,gam)

    use geometry,  only: sin_theta, cos_theta, omega, nr, dr
    use bl_constants_module
    use probin_module, only: base_cutoff_density, buoyancy_cutoff_factor, &
                             rotation_radius

    integer        , intent(in   ) ::  lo(:),hi(:),ng_f,ng_gp,ng_s, ng_uo, ng_um, n, ng_u0
    real(kind=dp_t), intent(inout) :: vel_force(lo(1)-ng_f :,lo(2)-ng_f :,lo(3)-ng_f :,:)
    logical        , intent(in   ) :: is_final_update
    real(kind=dp_t), intent(in   ) ::      uold(lo(1)-ng_uo:,lo(2)-ng_uo:,lo(3)-ng_uo:,:)
    real(kind=dp_t), intent(in   ) ::      umac(lo(1)-ng_um:,lo(2)-ng_um:,lo(3)-ng_um:)
    real(kind=dp_t), intent(in   ) ::      vmac(lo(1)-ng_um:,lo(2)-ng_um:,lo(3)-ng_um:)
    real(kind=dp_t), intent(in   ) ::      wmac(lo(1)-ng_um:,lo(2)-ng_um:,lo(3)-ng_um:)
    real(kind=dp_t), intent(in   ) ::   w0(0:)
    real(kind=dp_t), intent(in   ) ::     gpi(lo(1)-ng_gp:,lo(2)-ng_gp:,lo(3)-ng_gp:,:)
    real(kind=dp_t), intent(in   ) ::       D(lo(1)-ng_s :,lo(2)-ng_s :,lo(3)-ng_s :)
    real(kind=dp_t), intent(in   ) ::       Dh(lo(1)-ng_s :,lo(2)-ng_s :,lo(3)-ng_s :)
    real(kind=dp_t), intent(in   ) :: D0(0:)
    real(kind=dp_t), intent(in   ) :: Dh0(0:)
    real(kind=dp_t), intent(in   ) :: grav(0:)
    real(kind=dp_t), intent(in   ) :: w0_force(0:)
    logical        , intent(in   ) :: do_add_utilde_force
    real(kind=dp_t), intent(in   ) ::      u0(lo(1)-ng_u0:,lo(2)-ng_u0:,lo(3)-ng_u0:)
    real(kind=dp_t), intent(in   ) :: chrls(0:,0:,0:,0:,0:,0:)
    real(kind=dp_t), intent(in   ) :: gam(lo(1)-ng_u0:,lo(2)-ng_u0:,lo(3)-ng_u0:,:)

    integer         :: i,j,k,m
    real(kind=dp_t) :: Dpert

    real(kind=dp_t) :: coriolis_term(3), centrifugal_term(3), gr_term(3)

    vel_force = ZERO

    ! CURRENTLY for rotation in plane-parallel, we make the (bad) assumption
    ! that all points within the patch have the same centrifugal forcing terms.
    !
    ! We assume the centrifugal term applies at a constant radius,
    ! rotation_radius, for the patch.  In otherwords, the patch lives on the
    ! surface of a sphere of radius rotation_radius.
    !
    ! Furthermore, we assume the patch lives at longitude = 0.
    !
    ! Then the orientation of the patch is such that e_z is in the
    ! outward radial direction of the star, e_x is in the co_latitude (polar)
    ! angle direction and e_y is in the global y-direction.
    !
    ! centrifugal_term = omega x (omega x r) = (omega dot r) * omega
    !                                          - omega^2 * r
    ! where omega = (-|omega| sin_theta) e_x + (|omega| cos_theta) e_z
    !           r = rotation_radius e_z
    !
    ! See docs/rotation for derivation and figures.
    !
    gr_term = ZERO

    centrifugal_term(1) = - omega**2 * rotation_radius * sin_theta * sin_theta
    centrifugal_term(2) = ZERO
    centrifugal_term(3) = omega**2 * rotation_radius * cos_theta * sin_theta &
                          - omega**2 * rotation_radius

    !$OMP PARALLEL DO PRIVATE(i,j,k,Dpert,coriolis_term)
    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)

             Dpert = D(i,j,k) - D0(k)

             ! cutoff the buoyancy term if we are outside of the star
             if (D(i,j,k) .lt. buoyancy_cutoff_factor*base_cutoff_density) then
                Dpert = 0.d0
             end if

             ! the coriolis term is:
             !    TWO * omega x U
             ! where omega is given above and U = (u, v, w) is the velocity

             if (is_final_update) then

                ! use umac so we are time-centered
                coriolis_term(1) = -TWO * omega * &
                     HALF*(vmac(i,j,k) + vmac(i,j+1,k)) * cos_theta

                coriolis_term(2) =  TWO * omega * &
                     (HALF*(wmac(i,j,k)   + w0(k) + &
                            wmac(i,j,k+1) + w0(k+1)) * sin_theta + &
                      HALF*(umac(i,j,k) + umac(i+1,j,k)) * cos_theta)

                coriolis_term(3) = -TWO * omega * &
                     HALF*(vmac(i,j,k) + vmac(i,j+1,k)) * sin_theta

             else
                coriolis_term(1) = -TWO * omega * uold(i,j,k,2) * cos_theta

                coriolis_term(2) =  TWO * omega * ((uold(i,j,k,3) + HALF*(w0(k) + w0(k+1))) * sin_theta + &
                                                   uold(i,j,k,1) * cos_theta)

                coriolis_term(3) = -TWO * omega * uold(i,j,k,2) * sin_theta
             endif

             ! NOTE: using the fact that the gamma 3-metric is diagonal
             ! Do 0 index separately as need time component of velocity, which for conserved variables is just 1
             m = 0
             gr_term(1) = gr_term(1) + gam(i,j,k,1) * &
              (chrls(m,1,1,i,j,k) * (umac(i,j,k) + umac(i+1,j,k)) + &
               chrls(m,2,1,i,j,k) * (vmac(i,j,k) + vmac(i+1,j,k)) + &
               chrls(m,3,1,i,j,k) * (wmac(i,j,k) + wmac(i+1,j,k)))
             gr_term(2) = gr_term(2) + gam(i,j,k,2) * &
              (chrls(m,1,2,i,j,k) * (umac(i,j,k) + umac(i,j+1,k)) + &
               chrls(m,2,2,i,j,k) * (vmac(i,j,k) + vmac(i,j+1,k)) + &
               chrls(m,3,2,i,j,k) * (wmac(i,j,k) + wmac(i,j+1,k)))
             gr_term(3) = gr_term(3) + gam(i,j,k,3) * &
              (chrls(m,1,3,i,j,k) * (umac(i,j,k) + umac(i,j,k+1)) + &
               chrls(m,2,3,i,j,k) * (vmac(i,j,k) + vmac(i,j,k+1)) + &
               chrls(m,3,3,i,j,k) * (wmac(i,j,k) + wmac(i,j,k+1)))

             do m = 1, 3
                gr_term(1) = gr_term(1) + gam(i,j,k,1) * uold(i,j,k,m) * &
                 (chrls(m,1,1,i,j,k) * (umac(i,j,k) + umac(i+1,j,k)) + &
                  chrls(m,2,1,i,j,k) * (vmac(i,j,k) + vmac(i+1,j,k)) + &
                  chrls(m,3,1,i,j,k) * (wmac(i,j,k) + wmac(i+1,j,k)))
                gr_term(2) = gr_term(2) + gam(i,j,k,2) * uold(i,j,k,m) * &
                 (chrls(m,1,2,i,j,k) * (umac(i,j,k) + umac(i,j+1,k)) + &
                  chrls(m,2,2,i,j,k) * (vmac(i,j,k) + vmac(i,j+1,k)) + &
                  chrls(m,3,2,i,j,k) * (wmac(i,j,k) + wmac(i,j+1,k)))
                gr_term(3) = gr_term(3) + gam(i,j,k,3) * uold(i,j,k,m) * &
                 (chrls(m,1,3,i,j,k) * (umac(i,j,k) + umac(i,j,k+1)) + &
                  chrls(m,2,3,i,j,k) * (vmac(i,j,k) + vmac(i,j,k+1)) + &
                  chrls(m,3,3,i,j,k) * (wmac(i,j,k) + wmac(i,j,k+1)))
             enddo


             ! note: if use_alt_energy_fix = T, then gphi is already
             ! weighted by beta_0
             vel_force(i,j,k,1) = -coriolis_term(1) - centrifugal_term(1) - &
                  gpi(i,j,k,1) / (Dh(i,j,k) * u0(i,j,k)) + gr_term(1)

             vel_force(i,j,k,2) = -coriolis_term(2) - centrifugal_term(2) - &
                  gpi(i,j,k,2) / (Dh(i,j,k) * u0(i,j,k)) + gr_term(2)

             vel_force(i,j,k,3) = -coriolis_term(3) - centrifugal_term(3) + &
                  ( Dpert * grav(k) - gpi(i,j,k,3) ) / (Dh(i,j,k) * u0(i,j,k)) &
                  - w0_force(k) + gr_term(3)

          end do
       end do
    end do
    !$OMP END PARALLEL DO


    if (do_add_utilde_force) then
       !$OMP PARALLEL DO PRIVATE(i,j,k)
       do k=lo(3),hi(3)
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)

                if (k .le. -1) then
                   ! do not modify force since dw0/dr=0
                else if (k .ge. nr(n)) then
                   ! do not modify force since dw0/dr=0
                else
                   vel_force(i,j,k,3) = vel_force(i,j,k,3) &
                        - (wmac(i,j,k+1)+wmac(i,j,k))*(w0(k+1)-w0(k)) / (TWO*dr(n))
                end if

             end do
          end do
       end do
       !$OMP END PARALLEL DO

    endif

  end subroutine mk_vel_force_3d_cart

  subroutine mk_vel_force_3d_sphr(vel_force,ng_f,is_final_update, &
                                  uold,ng_uo,normal,ng_n, &
                                  umac,vmac,wmac,ng_um, &
                                  w0_cart,ng_wc,gradw0_cart,ng_gw, &
                                  w0macx,w0macy,ng_wm, &
                                  gpi,ng_gp,D,ng_s, &
                                  D0,grav,w0_force_cart,ng_w,lo,hi,dx, &
                                  do_add_utilde_force)

    use fill_3d_module
    use bl_constants_module
    use geometry,  only: omega, center
    use probin_module, only: base_cutoff_density, buoyancy_cutoff_factor, prob_lo

    integer        , intent(in   ) :: lo(:),hi(:),ng_f,ng_gp,ng_s,ng_uo,ng_um,ng_wc,ng_wm,ng_w,ng_n,ng_gw
    real(kind=dp_t), intent(inout) :: vel_force(lo(1)-ng_f :,lo(2)-ng_f :,lo(3)-ng_f :,:)
    logical        , intent(in   ) :: is_final_update
    real(kind=dp_t), intent(in   ) ::       uold(lo(1)-ng_uo:,lo(2)-ng_uo:,lo(3)-ng_uo:,:)
    real(kind=dp_t), intent(in   ) ::     normal(lo(1)-ng_n :,lo(2)-ng_n :,lo(3)-ng_n :,:)
    real(kind=dp_t), intent(in   ) ::       umac(lo(1)-ng_um:,lo(2)-ng_um:,lo(3)-ng_um:)
    real(kind=dp_t), intent(in   ) ::       vmac(lo(1)-ng_um:,lo(2)-ng_um:,lo(3)-ng_um:)
    real(kind=dp_t), intent(in   ) ::       wmac(lo(1)-ng_um:,lo(2)-ng_um:,lo(3)-ng_um:)
    real(kind=dp_t), intent(in   ) ::    w0_cart(lo(1)-ng_wc:,lo(2)-ng_wc:,lo(3)-ng_wc:,:)
    real(kind=dp_t), intent(in   ) ::gradw0_cart(lo(1)-ng_wc:,lo(2)-ng_wc:,lo(3)-ng_wc:)
    real(kind=dp_t), intent(in   ) ::     w0macx(lo(1)-ng_wm:,lo(2)-ng_wm:,lo(3)-ng_wm:)
    real(kind=dp_t), intent(in   ) ::     w0macy(lo(1)-ng_wm:,lo(2)-ng_wm:,lo(3)-ng_wm:)
    real(kind=dp_t), intent(in   ) ::        gpi(lo(1)-ng_gp:,lo(2)-ng_gp:,lo(3)-ng_gp:,:)
    real(kind=dp_t), intent(in   ) ::        D(lo(1)-ng_s :,lo(2)-ng_s :,lo(3)-ng_s :)
    real(kind=dp_t), intent(in   ) :: w0_force_cart(lo(1)-ng_w:,lo(2)-ng_w:,lo(3)-ng_w:,:)
    real(kind=dp_t), intent(in   ) :: D0(0:)
    real(kind=dp_t), intent(in   ) :: grav(0:)
    real(kind=dp_t), intent(in   ) ::   dx(:)
    logical        , intent(in   ) :: do_add_utilde_force

    integer         :: i,j,k

    real(kind=dp_t), allocatable :: D0_cart(:,:,:,:)
    real(kind=dp_t), allocatable :: grav_cart(:,:,:,:)

    real(kind=dp_t) :: Dpert
    real(kind=dp_t) :: xx, yy, zz
    real(kind=dp_t) :: centrifugal_term(3), coriolis_term(3)

    real(kind=dp_t) :: Ut_dot_er

    allocate(D0_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1))
    allocate(grav_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),3))

    vel_force = ZERO

    call put_1d_array_on_cart_3d_sphr(.false.,.false.,D0,D0_cart,lo,hi,dx,0)
    call put_1d_array_on_cart_3d_sphr(.false.,.true.,grav,grav_cart,lo,hi,dx,0)

    !$OMP PARALLEL DO PRIVATE(i,j,k,xx,yy,zz,Dpert,centrifugal_term,coriolis_term)
    do k = lo(3),hi(3)
       zz = prob_lo(3) + (dble(k) + HALF)*dx(3) - center(3)
       do j = lo(2),hi(2)
          yy = prob_lo(2) + (dble(j) + HALF)*dx(2) - center(2)
          do i = lo(1),hi(1)
             xx = prob_lo(1) + (dble(i) + HALF)*dx(1) - center(1)

             Dpert = D(i,j,k) - D0_cart(i,j,k,1)

             ! cutoff the buoyancy term if we are outside of the star
             if (D(i,j,k) .lt. buoyancy_cutoff_factor*base_cutoff_density) then
                Dpert = 0.d0
             end if


             ! Coriolis and centrifugal forces.  We assume that the
             ! rotation axis is the z direction, with angular velocity
             ! omega

             ! omega x (omega x r ) = - omega^2 x e_x  - omega^2 y e_y
             ! (with omega = omega e_z)
             centrifugal_term(1) = -omega * omega * xx
             centrifugal_term(2) = -omega * omega * yy
             centrifugal_term(3) = ZERO

             ! cutoff the centrifugal term if we are outside the star
             if (D(i,j,k) .lt. buoyancy_cutoff_factor*base_cutoff_density) then
                centrifugal_term(:) = 0.d0
             end if


             ! 2 omega x U = - 2 omega v e_x  + 2 omega u e_y
             ! (with omega = omega e_z)
             if (is_final_update) then

                ! use umac so we are time-centered
                coriolis_term(1) = -TWO * omega * &
                     HALF*(vmac(i,j,k)   + w0macy(i,j,k) + &
                           vmac(i,j+1,k) + w0macy(i,j+1,k))

                coriolis_term(2) =  TWO * omega * &
                     HALF*(umac(i,j,k)   + w0macx(i,j,k) + &
                           umac(i+1,j,k) + w0macx(i+1,j,k))

                coriolis_term(3) = ZERO

             else
                coriolis_term(1) = -TWO * omega * (uold(i,j,k,2) + w0_cart(i,j,k,2))
                coriolis_term(2) =  TWO * omega * (uold(i,j,k,1) + w0_cart(i,j,k,1))
                coriolis_term(3) = ZERO
             endif


             ! F_Coriolis = -2 omega x U
             ! F_centrifugal = - omega x (omega x r)

             ! we just computed the absolute value of the forces above, so use
             ! the right sign here

             ! note: if use_alt_energy_fix = T, then gphi is already weighted
             ! by beta_0
             vel_force(i,j,k,1) = -coriolis_term(1) - centrifugal_term(1) + &
                  ( Dpert * grav_cart(i,j,k,1) - gpi(i,j,k,1) ) / D(i,j,k) &
                  - w0_force_cart(i,j,k,1)

             vel_force(i,j,k,2) = -coriolis_term(2) - centrifugal_term(2) + &
                  ( Dpert * grav_cart(i,j,k,2) - gpi(i,j,k,2) ) / D(i,j,k) &
                  - w0_force_cart(i,j,k,2)

             vel_force(i,j,k,3) = -coriolis_term(3) - centrifugal_term(3) + &
                  ( Dpert * grav_cart(i,j,k,3) - gpi(i,j,k,3) ) / D(i,j,k) &
                  - w0_force_cart(i,j,k,3)

          end do
       end do
    end do
    !$OMP END PARALLEL DO


    if (do_add_utilde_force) then

       !$OMP PARALLEL DO PRIVATE(i,j,k,Ut_dot_er)
       do k=lo(3),hi(3)
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)

                Ut_dot_er = &
                     HALF*(umac(i,j,k)+umac(i+1,j  ,k  ))*normal(i,j,k,1) + &
                     HALF*(vmac(i,j,k)+vmac(i  ,j+1,k  ))*normal(i,j,k,2) + &
                     HALF*(wmac(i,j,k)+wmac(i  ,j,  k+1))*normal(i,j,k,3)

                vel_force(i,j,k,1) = vel_force(i,j,k,1) - Ut_dot_er*gradw0_cart(i,j,k)*normal(i,j,k,1)
                vel_force(i,j,k,2) = vel_force(i,j,k,2) - Ut_dot_er*gradw0_cart(i,j,k)*normal(i,j,k,2)
                vel_force(i,j,k,3) = vel_force(i,j,k,3) - Ut_dot_er*gradw0_cart(i,j,k)*normal(i,j,k,3)

             end do
          end do
       end do
       !$OMP END PARALLEL DO

    endif

    deallocate(D0_cart,grav_cart)

  end subroutine mk_vel_force_3d_sphr

end module mk_vel_force_module
