module enthalpy_advance_module

  use bl_types
  use multifab_module
  use ml_layout_module
  use define_bc_module

  implicit none

  private

  public :: enthalpy_advance

contains

  subroutine enthalpy_advance(mla,which_step,uold,sold,snew,sedge,sflux,scal_force,&
                              thermal,umac,w0,w0mac, &
                              D0_old,Dh0_old,D0_new,Dh0_new,p0_old,p0_new, &
                              tempbar,psi,dx,dt,the_bc_level,u0_1d,alpha,beta,gam,u0)

    use bl_prof_module
    use bl_constants_module
    use make_edge_scal_module
    use bds_module
    use mkflux_module
    use mkscalforce_module
    use update_scal_module
    use addw0_module
    use define_bc_module
    use fill_3d_module
    use pert_form_module
    use cell_to_edge_module
    use rhoh_vs_t_module
    use geometry,      only: spherical, nr_fine, r_start_coord, r_end_coord, &
         numdisjointchunks, nlevs_radial
    use variables,     only: temp_comp, rho_comp, rhoh_comp, foextrap_comp, nscal
    use probin_module, only: enthalpy_pred_type, verbose, bds_type
    use pred_parameters
    use modify_scal_force_module
    use convert_rhoX_to_X_module
    use metric_module, only: cons_to_prim, prim_to_cons

    type(ml_layout), intent(inout) :: mla
    integer        , intent(in   ) :: which_step
    type(multifab) , intent(in   ) :: uold(:)
    type(multifab) , intent(inout) :: sold(:)
    type(multifab) , intent(inout) :: snew(:)
    type(multifab) , intent(inout) :: sedge(:,:)
    type(multifab) , intent(inout) :: sflux(:,:)
    type(multifab) , intent(inout) :: scal_force(:)
    type(multifab) , intent(in   ) :: thermal(:)
    type(multifab) , intent(inout) :: umac(:,:)
    real(kind=dp_t), intent(in   ) :: w0(:,0:)
    type(multifab) , intent(in   ) :: w0mac(:,:)
    real(kind=dp_t), intent(in   ) :: D0_old(:,0:)
    real(kind=dp_t), intent(in   ) :: Dh0_old(:,0:)
    real(kind=dp_t), intent(in   ) :: D0_new(:,0:)
    real(kind=dp_t), intent(in   ) :: Dh0_new(:,0:)
    real(kind=dp_t), intent(in   ) :: p0_old(:,0:)
    real(kind=dp_t), intent(in   ) :: p0_new(:,0:)
    real(kind=dp_t), intent(in   ) :: tempbar(:,0:)
    real(kind=dp_t), intent(in   ) :: psi(:,0:)
    real(kind=dp_t), intent(in   ) :: dx(:,:),dt
    type(bc_level) , intent(in   ) :: the_bc_level(:)
    real(kind=dp_t), intent(in   ) :: u0_1d(:,0:)
    type(multifab), intent(in) :: alpha(:)
    type(multifab), intent(in) :: beta(:)
    type(multifab), intent(in) :: gam(:)
    type(multifab),  intent(in) :: u0(:)

    type(multifab) :: Dh0_old_cart(mla%nlevel)
    type(multifab) :: p0_new_cart(mla%nlevel)

    type(multifab) :: D0mac_old(mla%nlevel,mla%dim)
    type(multifab) :: D0mac_new(mla%nlevel,mla%dim)
    type(multifab) :: Dh0mac_old(mla%nlevel,mla%dim)
    type(multifab) :: Dh0mac_new(mla%nlevel,mla%dim)
    type(multifab) :: h0mac_old(mla%nlevel,mla%dim)
    type(multifab) :: h0mac_new(mla%nlevel,mla%dim)
    type(multifab) :: s_prim_edge(mla%nlevel,mla%dim)
    type(multifab) :: s_prim(mla%nlevel)
    type(multifab) :: u_prim(mla%nlevel)

    integer    :: pred_comp,n,r,i,comp,dm,nlevs,ng_s
    logical    :: is_vel
    real(dp_t) :: smin,smax
    logical    :: is_prediction

    ! Create cell-centered base state quantity
    real(kind=dp_t), allocatable :: h0_old(:,:)
    real(kind=dp_t), allocatable :: h0_new(:,:)

    ! Create edge-centered base state quantities.
    ! Note: D0_edge_{old,new} and Dh0_edge_{old,new}
    ! contain edge-centered quantities created via spatial interpolation.
    real(kind=dp_t), allocatable ::  rho0_old(:,:)
    real(kind=dp_t), allocatable ::  rho0_new(:,:)
    real(kind=dp_t), allocatable :: rhoh0_old(:,:)
    real(kind=dp_t), allocatable :: rhoh0_new(:,:)
    real(kind=dp_t), allocatable ::  D0_edge_old(:,:)
    real(kind=dp_t), allocatable ::  D0_edge_new(:,:)
    real(kind=dp_t), allocatable :: Dh0_edge_old(:,:)
    real(kind=dp_t), allocatable :: Dh0_edge_new(:,:)
    real(kind=dp_t), allocatable ::    t0_edge_old(:,:)
    real(kind=dp_t), allocatable ::    t0_edge_new(:,:)
    real(kind=dp_t), allocatable ::  rho0_edge_old(:,:)
    real(kind=dp_t), allocatable ::  rho0_edge_new(:,:)
    real(kind=dp_t), allocatable :: rhoh0_edge_old(:,:)
    real(kind=dp_t), allocatable :: rhoh0_edge_new(:,:)


    type(bl_prof_timer), save :: bpt

    call build(bpt, "enthalpy_advance")

    allocate( h0_old(nlevs_radial,0:nr_fine-1))
    allocate( h0_new(nlevs_radial,0:nr_fine-1))
    allocate( rho0_old(nlevs_radial,0:nr_fine-1))
    allocate( rho0_new(nlevs_radial,0:nr_fine-1))
    allocate( rhoh0_old(nlevs_radial,0:nr_fine-1))
    allocate( rhoh0_new(nlevs_radial,0:nr_fine-1))
    allocate(  D0_edge_old(nlevs_radial,0:nr_fine))
    allocate(  D0_edge_new(nlevs_radial,0:nr_fine))
    allocate( Dh0_edge_old(nlevs_radial,0:nr_fine))
    allocate( Dh0_edge_new(nlevs_radial,0:nr_fine))
    allocate(  rho0_edge_old(nlevs_radial,0:nr_fine))
    allocate(  rho0_edge_new(nlevs_radial,0:nr_fine))
    allocate( rhoh0_edge_old(nlevs_radial,0:nr_fine))
    allocate( rhoh0_edge_new(nlevs_radial,0:nr_fine))
    allocate(    t0_edge_old(nlevs_radial,0:nr_fine))
    allocate(    t0_edge_new(nlevs_radial,0:nr_fine))

    is_vel  = .false.
    dm = mla%dim
    nlevs = mla%nlevel

    if (spherical .eq. 0) then
       call cell_to_edge(D0_old,D0_edge_old)
       call cell_to_edge(D0_new,D0_edge_new)
       call cell_to_edge(Dh0_old,Dh0_edge_old)
       call cell_to_edge(Dh0_new,Dh0_edge_new)
       call cell_to_edge(tempbar,t0_edge_old)
       call cell_to_edge(tempbar,t0_edge_new)
    end if

    if (enthalpy_pred_type .eq. predict_h .or. &
        enthalpy_pred_type .eq. predict_hprime) then
       ! convert (rho h) -> h
       call convert_rhoh_to_h(sold,.true.,mla,the_bc_level)
    end if

    !**************************************************************************
    ! Create scalar source term at time n
    !**************************************************************************

    do n = 1, nlevs
       call setval(scal_force(n),ZERO,all=.true.)
    end do

    ! compute forcing terms
    if (enthalpy_pred_type .eq. predict_rhohprime) then

       ! make force for (rho h)'
       is_prediction = .true.
       call mkrhohforce(mla,scal_force,is_prediction, &
                        thermal,umac,p0_old,p0_old,D0_old,D0_old,&
                        Dh0_old,Dh0_new, &
                        psi,dx,.true.,the_bc_level,u0_1d)

       do n=1,nlevs
          call build(Dh0_old_cart(n), get_layout(sold(n)), 1, 1)
       end do

       call put_1d_array_on_cart(Dh0_old,Dh0_old_cart,dm+rhoh_comp,.false., &
                                 .false.,dx,the_bc_level,mla)

       call modify_scal_force(scal_force,sold,umac,Dh0_old, &
                              Dh0_edge_old,w0,dx,Dh0_old_cart,rhoh_comp,mla,the_bc_level)

       do n=1,nlevs
          call destroy(Dh0_old_cart(n))
       end do

    else if (enthalpy_pred_type .eq. predict_h .or. &
             enthalpy_pred_type .eq. predict_rhoh) then

       is_prediction = .true.
       call mkrhohforce(mla,scal_force,is_prediction,&
                        thermal,umac,p0_old,p0_old,D0_old,D0_old,&
                        Dh0_old,Dh0_new, &
                        psi,dx,.true.,the_bc_level,u0_1d)

       if (enthalpy_pred_type .eq. predict_h) then
          ! make force for h by calling mkrhohforce then dividing by rho
          do n=1,nlevs
             call multifab_div_div_c(scal_force(n),rhoh_comp,sold(n),rho_comp,1,1)
          end do
       end if

    else if (enthalpy_pred_type .eq. predict_hprime) then

       ! first compute h0_old
       do n=1,nlevs_radial
          do i=1,numdisjointchunks(n)
             do r=r_start_coord(n,i),r_end_coord(n,i)
                h0_old(n,r) = Dh0_old(n,r) / D0_old(n,r)
             end do
          end do
       end do

       ! make force for hprime
       is_prediction = .true.
       call mkhprimeforce(mla,sold,sold,scal_force,is_prediction,thermal,umac,p0_old, &
                          p0_old,h0_old,h0_old,psi,dx,.true.,the_bc_level)

    else if ( (enthalpy_pred_type .eq. predict_T_then_rhohprime) .or. &
              (enthalpy_pred_type .eq. predict_T_then_h        ) .or. &
              (enthalpy_pred_type .eq. predict_Tprime_then_h) ) then

       ! make force for temperature
       ! FIXME: this should take the primitive variables, but as function
       ! isn't called in my code, not going to bother changing.
       call mktempforce(mla,scal_force,umac,sold,thermal,p0_old,p0_old,psi,dx,the_bc_level)

    end if

    !**************************************************************************
    !     Add w0 to MAC velocities (trans velocities already have w0).
    !**************************************************************************

    call addw0(umac,the_bc_level,mla,w0,w0mac,mult=ONE)

    !**************************************************************************
    !     Create the edge states of (rho h)' or h or T
    !**************************************************************************

    if (enthalpy_pred_type .eq. predict_rhohprime) then
       ! convert (rho h) -> (rho h)' or
       call put_in_pert_form(mla,sold,Dh0_old,dx,rhoh_comp,foextrap_comp,.true., &
                             the_bc_level)
    end if

    if (enthalpy_pred_type .eq. predict_hprime) then
       ! convert h -> h'
       call put_in_pert_form(mla,sold,h0_old,dx,rhoh_comp,foextrap_comp,.true.,the_bc_level)
    end if

    if (enthalpy_pred_type .eq. predict_Tprime_then_h) then
       ! convert T -> T'
       call put_in_pert_form(mla,sold,tempbar,dx,temp_comp,foextrap_comp,.true.,the_bc_level)
    end if

    ! predict either T, h, or (rho h)' at the edges
    if ( (enthalpy_pred_type .eq. predict_T_then_rhohprime) .or. &
         (enthalpy_pred_type .eq. predict_T_then_h        ) .or. &
         (enthalpy_pred_type .eq. predict_Tprime_then_h) ) then
       pred_comp = temp_comp
    else
       pred_comp = rhoh_comp
    end if

    if (enthalpy_pred_type .eq. predict_rhoh) then
       ! use the conservative form of the prediction
       if (bds_type .eq. 0) then
          call make_edge_scal(sold,sedge,umac,scal_force, &
                              dx,dt,is_vel,the_bc_level, &
                              pred_comp,dm+pred_comp,1,.true.,mla)
       else if (bds_type .eq. 1) then
          call bds(sold,sedge,umac,scal_force, &
                   dx,dt,is_vel,the_bc_level, &
                   pred_comp,dm+pred_comp,1,.true.,mla)
       end if
    else
       ! use the advective form of the prediction
       if (bds_type .eq. 0) then
          call make_edge_scal(sold,sedge,umac,scal_force, &
                              dx,dt,is_vel,the_bc_level, &
                              pred_comp,dm+pred_comp,1,.false.,mla)
       else if (bds_type .eq. 1) then
          call bds(sold,sedge,umac,scal_force, &
                   dx,dt,is_vel,the_bc_level, &
                   pred_comp,dm+pred_comp,1,.false.,mla)
       end if
    end if

    if (enthalpy_pred_type .eq. predict_rhohprime) then
       ! convert (rho h)' -> (rho h)
       call put_in_pert_form(mla,sold,Dh0_old,dx,rhoh_comp,dm+rhoh_comp,.false., &
                             the_bc_level)
    end if

    if (enthalpy_pred_type .eq. predict_hprime) then
       ! convert h' -> h
       call put_in_pert_form(mla,sold,h0_old,dx,rhoh_comp,foextrap_comp,.false.,the_bc_level)
    end if

    if (enthalpy_pred_type .eq. predict_Tprime_then_h) then
       ! convert T' -> T
       call put_in_pert_form(mla,sold,tempbar,dx,temp_comp,dm+temp_comp,.false., &
                             the_bc_level)
    end if

    if (enthalpy_pred_type .eq. predict_h .or. &
        enthalpy_pred_type .eq. predict_hprime) then
       ! convert h -> (rho h)
       call convert_rhoh_to_h(sold,.false.,mla,the_bc_level)
    end if

    ! FIXME: need primitive variables here
    ! It's going to be really hard to do this properly, so going to just do it
    ! in a really hacky way for now.
    rho0_old(:,:) = D0_old(:,:) / u0_1d(:,:)
    rho0_new(:,:) = D0_new(:,:) / u0_1d(:,:)
    rhoh0_old(:,:) = Dh0_old(:,:) / u0_1d(:,:)
    rhoh0_new(:,:) = Dh0_new(:,:) / u0_1d(:,:)

    if (spherical .eq. 0) then
       call cell_to_edge(rho0_old,rho0_edge_old)
       call cell_to_edge(rho0_new,rho0_edge_new)
       call cell_to_edge(rhoh0_old,rhoh0_edge_old)
       call cell_to_edge(rhoh0_new,rhoh0_edge_new)
    end if

    do n=1,nlevs
       ng_s = nghost(sold(n))
       do comp=1,dm
           call multifab_build_edge(s_prim_edge(n,comp),mla%la(n),nscal, ng_s,comp)
       enddo
       call multifab_build(s_prim(n), mla%la(n), nscal, ng_s)
       call multifab_build(u_prim(n), mla%la(n), dm, ng_s)
    end do
    call cons_to_prim(sold, uold, alpha, beta, gam, s_prim, u_prim, mla,the_bc_level)

    if (enthalpy_pred_type .eq. predict_rhoh) then
       ! use the conservative form of the prediction
       if (bds_type .eq. 0) then
          call make_edge_scal(s_prim,s_prim_edge,umac,scal_force, &
                              dx,dt,is_vel,the_bc_level, &
                              pred_comp,dm+pred_comp,1,.true.,mla)
       else if (bds_type .eq. 1) then
          call bds(s_prim,s_prim_edge,umac,scal_force, &
                   dx,dt,is_vel,the_bc_level, &
                   pred_comp,dm+pred_comp,1,.true.,mla)
       end if
    else
       ! use the advective form of the prediction
       if (bds_type .eq. 0) then
          call make_edge_scal(s_prim,s_prim_edge,umac,scal_force, &
                              dx,dt,is_vel,the_bc_level, &
                              pred_comp,dm+pred_comp,1,.false.,mla)
       else if (bds_type .eq. 1) then
          call bds(s_prim,s_prim_edge,umac,scal_force, &
                   dx,dt,is_vel,the_bc_level, &
                   pred_comp,dm+pred_comp,1,.false.,mla)
       end if
    end if

    ! Compute enthalpy edge states if we were predicting temperature.  This
    ! needs to be done after the state was returned to the full state.
    if ( (enthalpy_pred_type .eq. predict_T_then_rhohprime) .or. &
         (enthalpy_pred_type .eq. predict_T_then_h        ) .or. &
         (enthalpy_pred_type .eq. predict_Tprime_then_h) ) then
       call makeHfromRhoT_edge(u_prim,s_prim_edge,rho0_old,rhoh0_old,tempbar,rho0_edge_old, &
                               rhoh0_edge_old,t0_edge_old,rho0_new,rhoh0_new,tempbar, &
                               rho0_edge_new,rhoh0_edge_new,t0_edge_new,the_bc_level,dx,mla)
    end if

    !**************************************************************************
    !     Subtract w0 from MAC velocities.
    !**************************************************************************

    call addw0(umac,the_bc_level,mla,w0,w0mac,mult=-ONE)

    !**************************************************************************
    !     Compute fluxes
    !**************************************************************************

    ! for which_step .eq. 1, we pass in only the old base state quantities
    ! for which_step .eq. 2, we pass in the old and new for averaging within mkflux
    if (which_step .eq. 1) then

       if (spherical .eq. 1) then
          do n=1,nlevs
             do comp=1,dm
                call multifab_build_edge(D0mac_old(n,comp),mla%la(n),1,1,comp)
                call multifab_build_edge(Dh0mac_old(n,comp),mla%la(n),1,1,comp)
                call multifab_build_edge(h0mac_old(n,comp),mla%la(n),1,1,comp)
             end do
          end do

          do n=1,nlevs_radial
             do i=1,numdisjointchunks(n)
                do r=r_start_coord(n,i),r_end_coord(n,i)
                   h0_old(n,r) = Dh0_old(n,r) / D0_old(n,r)
                end do
             end do
          end do

          call make_s0mac(mla,D0_old,D0mac_old,dx,dm+rho_comp,the_bc_level)
          call make_s0mac(mla,Dh0_old,Dh0mac_old,dx,dm+rhoh_comp,the_bc_level)
          call make_s0mac(mla,h0_old,h0mac_old,dx,foextrap_comp,the_bc_level)
       end if

       ! compute enthalpy fluxes
       call mk_rhoh_flux(mla,sflux,sold,sedge,umac,w0,w0mac, &
                         D0_old,D0_edge_old,D0mac_old, &
                         D0_old,D0_edge_old,D0mac_old, &
                         Dh0_old,Dh0_edge_old,Dh0mac_old, &
                         Dh0_old,Dh0_edge_old,Dh0mac_old, &
                         h0mac_old,h0mac_old)

      if (spherical .eq. 1) then
          do n=1,nlevs
             do comp=1,dm
                call destroy(D0mac_old(n,comp))
                call destroy(Dh0mac_old(n,comp))
                call destroy(h0mac_old(n,comp))
                call destroy(s_prim_edge(n,comp))
            enddo
            call destroy(s_prim(n))
            call destroy(u_prim(n))
          end do
       end if

    else if (which_step .eq. 2) then

       if (spherical .eq. 1) then
          do n=1,nlevs
             do comp=1,dm
                call multifab_build_edge( D0mac_old(n,comp),mla%la(n),1,1,comp)
                call multifab_build_edge(Dh0mac_old(n,comp),mla%la(n),1,1,comp)
                call multifab_build_edge(   h0mac_old(n,comp),mla%la(n),1,1,comp)
                call multifab_build_edge( D0mac_new(n,comp),mla%la(n),1,1,comp)
                call multifab_build_edge(Dh0mac_new(n,comp),mla%la(n),1,1,comp)
                call multifab_build_edge(   h0mac_new(n,comp),mla%la(n),1,1,comp)
             end do
          end do

          do n=1,nlevs_radial
             do i=1,numdisjointchunks(n)
                do r=r_start_coord(n,i),r_end_coord(n,i)
                   h0_old(n,r) = Dh0_old(n,r) / D0_old(n,r)
                   h0_new(n,r) = Dh0_new(n,r) / D0_new(n,r)
                end do
             end do
          end do

          call make_s0mac(mla,D0_old,D0mac_old,dx,dm+rho_comp,the_bc_level)
          call make_s0mac(mla,Dh0_old,Dh0mac_old,dx,dm+rhoh_comp,the_bc_level)
          call make_s0mac(mla,h0_old,h0mac_old,dx,foextrap_comp,the_bc_level)
          call make_s0mac(mla,D0_new,D0mac_new,dx,dm+rho_comp,the_bc_level)
          call make_s0mac(mla,Dh0_new,Dh0mac_new,dx,dm+rhoh_comp,the_bc_level)
          call make_s0mac(mla,h0_new,h0mac_new,dx,foextrap_comp,the_bc_level)
       end if

       ! compute enthalpy fluxes
       call mk_rhoh_flux(mla,sflux,sold,sedge,umac,w0,w0mac, &
                         D0_old,D0_edge_old,D0mac_old, &
                         D0_new,D0_edge_new,D0mac_new, &
                         Dh0_old,Dh0_edge_old,Dh0mac_old, &
                         Dh0_new,Dh0_edge_new,Dh0mac_new, &
                         h0mac_old,h0mac_new)

      if (spherical .eq. 1) then
          do n=1,nlevs
             do comp=1,dm
                call destroy(D0mac_old(n,comp))
                call destroy(Dh0mac_old(n,comp))
                call destroy(h0mac_old(n,comp))
                call destroy(D0mac_new(n,comp))
                call destroy(Dh0mac_new(n,comp))
                call destroy(h0mac_new(n,comp))
                call destroy(s_prim_edge(n,comp))
            enddo
            call destroy(s_prim(n))
            call destroy(u_prim(n))
          end do
       end if

    end if

    do n=1,nlevs
       call setval(scal_force(n),ZERO,all=.true.)
    end do

    !**************************************************************************
    !     1) Create (rho h)' force at time n+1/2.
    !          (NOTE: we don't worry about filling ghost cells of the scal_force
    !                 because we only need them in valid regions...)
    !     2) Update (rho h) with conservative differencing.
    !**************************************************************************

    if (which_step .eq. 1) then
      ! Here just send p0_old and p0_old
       is_prediction = .false.
       call mkrhohforce(mla,scal_force,is_prediction, &
                        thermal,umac,p0_old,p0_old,D0_old,D0_old,&
                        Dh0_old,Dh0_new, &
                        psi,dx,.false.,the_bc_level,u0_1d)
    else
      ! Here send p0_old and p0_new
       is_prediction = .false.
       call mkrhohforce(mla,scal_force,is_prediction, &
                        thermal,umac,p0_old,p0_new,D0_old,D0_new,&
                        Dh0_old,Dh0_new, &
                        psi,dx,.false.,the_bc_level,u0_1d)
    end if

    if (spherical .eq. 1) then
       do n=1,nlevs
          call build(p0_new_cart(n), get_layout(sold(n)), 1, 1)
       end do

       call put_1d_array_on_cart(p0_new,p0_new_cart,foextrap_comp,.false., &
                                 .false.,dx,the_bc_level,mla)
    end if

    call update_scal(mla,rhoh_comp,rhoh_comp,sold,snew,sflux,scal_force, &
                     p0_new,p0_new_cart,dx,dt,the_bc_level,uold,alpha,beta,gam,u0)

    if (spherical .eq. 1) then
       do n=1,nlevs
          call destroy(p0_new_cart(n))
       end do
    end if


    if ( verbose .ge. 1 ) then
       do n=1,nlevs
          smin = multifab_min_c(snew(n),rhoh_comp)
          smax = multifab_max_c(snew(n),rhoh_comp)
          if (parallel_IOProcessor()) then
             write (6,1999) n
             write (6,2001) smin,smax
          end if
       end do
    end if

    if (parallel_IOProcessor()) write(6,2004)

    call destroy(bpt)

1999 format('... Level ', i1, ' update:')
2001 format('... new min/max : rho * H           ',e17.10,2x,e17.10)
2004 format(' ')

  end subroutine enthalpy_advance

end module enthalpy_advance_module
