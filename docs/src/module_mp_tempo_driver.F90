module module_mp_tempo_driver
  !! tempo driver that converts 3d model input to 1d arrays used by the main code
  !! also allocates and fills diagnostic arrays
  use module_mp_tempo_params, only : wp, sp, dp, tempo_cfgs
  use module_mp_tempo_main, only : tempo_main, ty_tempo_main_diags
  implicit none
  private

  public :: tempo_driver, ty_tempo_driver_diags, tempo_aerosol_surface_emissions

  type :: ty_tempo_driver_diags
    real(wp), dimension(:,:), allocatable :: rain_precip
    real(wp), dimension(:,:), allocatable :: ice_liquid_equiv_precip
    real(wp), dimension(:,:), allocatable :: snow_liquid_equiv_precip
    real(wp), dimension(:,:), allocatable :: graupel_liquid_equiv_precip
    real(wp), dimension(:,:), allocatable :: frozen_fraction
    real(wp), dimension(:,:,:), allocatable :: refl10cm
    real(wp), dimension(:,:,:), allocatable :: re_cloud
    real(wp), dimension(:,:,:), allocatable :: re_ice
    real(wp), dimension(:,:,:), allocatable :: re_snow
    real(wp), dimension(:,:,:), allocatable :: max_hail_diameter
    real(wp), dimension(:,:,:), allocatable :: rain_med_vol_diam
    real(wp), dimension(:,:,:), allocatable :: graupel_med_vol_diam
  end type

  contains

  subroutine tempo_driver(dt, itimestep, &
    t, th, pii, p, w, dz, &
    qv, qc, qr, qi, qs, qg, ni, nr, &
    nc, nwfa, nifa, ng, qb, &
    ids, ide, jds, jde, kds, kde, &
    ims, ime, jms, jme, kms, kme, &
    its, ite, jts, jte, kts, kte, tempo_diags)

    real(wp), intent(in) :: dt !! timestep \([s]]\)
    integer, intent(in) :: itimestep !! integer timestep = integration time / dt
    integer, intent(in) :: ids, ide, jds, jde, kds, kde !! domain locations
    integer, intent(in) :: ims, ime, jms, jme, kms, kme !! memory locations
    integer, intent(in) :: its, ite, jts, jte, kts, kte !! tile locations

    real(wp), dimension(ims:ime, kms:kme, jms:jme), intent(inout), optional :: t !! temperature \([K]\)
    real(wp), dimension(ims:ime, kms:kme, jms:jme), intent(inout), optional :: th !! theta \([K]\)

    real(wp), dimension(ims:ime, kms:kme, jms:jme), intent(in) :: p !! pressure \([Pa]\)
    real(wp), dimension(ims:ime, kms:kme, jms:jme), intent(in) :: w !! vertical velocity \([m\; s^{-1}]\)
    real(wp), dimension(ims:ime, kms:kme, jms:jme), intent(in) :: dz !! vertical grid spacing \([m]\)
    real(wp), dimension(ims:ime, kms:kme, jms:jme), intent(in), optional :: pii !! exner function

    real(wp), dimension(ims:ime, kms:kme, jms:jme), intent(inout) :: qv !! 3D water vapor mixing ratio \([kg\; kg^{-1}]\)
    real(wp), dimension(ims:ime, kms:kme, jms:jme), intent(inout) :: qc !! 3D cloud water mass mixing ratio \([kg\; kg^{-1}]\)
    real(wp), dimension(ims:ime, kms:kme, jms:jme), intent(inout) :: qr !! 3D rain water mass mixing ratio \([kg\; kg^{-1}]\)
    real(wp), dimension(ims:ime, kms:kme, jms:jme), intent(inout) :: qi !! 3D cloud ice mass mixing ratio \([kg\; kg^{-1}]\)
    real(wp), dimension(ims:ime, kms:kme, jms:jme), intent(inout) :: qs !! 3D snow mass mixing ratio \([kg\; kg^{-1}]\)
    real(wp), dimension(ims:ime, kms:kme, jms:jme), intent(inout) :: qg !! 3D graupel mass mixing ratio \([kg\; kg^{-1}]\)
    real(wp), dimension(ims:ime, kms:kme, jms:jme), intent(inout) :: ni !! 3D cloud ice number mixing ratio \([kg^{-1}]\)
    real(wp), dimension(ims:ime, kms:kme, jms:jme), intent(inout) :: nr !! 3D rain water number mixing ratio \([kg^{-1}]\)

    real(wp), dimension(ims:ime, kms:kme, jms:jme), intent(inout), optional :: nc !! 3D cloud water number mixing ratio \([kg^{-1}]\) (aerosol-aware)
    real(wp), dimension(ims:ime, kms:kme, jms:jme), intent(inout), optional :: nwfa !! 3D water-friendly aerosol number mixing ratio \([kg^{-1}]\) (aerosol-aware)
    real(wp), dimension(ims:ime, kms:kme, jms:jme), intent(inout), optional :: nifa !! 3D ice-friendly aerosol number mixing ratio \([kg^{-1}]\) (aerosol-aware)
    real(wp), dimension(ims:ime, kms:kme, jms:jme), intent(inout), optional :: qb !! 3D graupel volume mixing ratio \([m^{-3}\; kg^{-1}]\) (hail-aware)
    real(wp), dimension(ims:ime, kms:kme, jms:jme), intent(inout), optional :: ng !! 3D graupel number mixing ratio \([kg^{-1}]\) (hail-aware)

    real(wp), dimension(kts:kte) :: t1d !! 1D temperature \([K]\)
    real(wp), dimension(kts:kte) :: p1d !! 1D pressure \([Pa]\)
    real(wp), dimension(kts:kte) :: qv1d !! 1D water vapor mixing ratio \([kg\; kg^{-1}]\)
    real(wp), dimension(kts:kte) :: qc1d !! 1D cloud water mass mixing ratio \([kg\; kg^{-1}]\)
    real(wp), dimension(kts:kte) :: qr1d !! 1D rain water mass mixing ratio \([kg\; kg^{-1}]\)
    real(wp), dimension(kts:kte) :: qi1d !! 1D cloud ice mass mixing ratio \([kg\; kg^{-1}]\)
    real(wp), dimension(kts:kte) :: qs1d !! 1D snow mass mixing ratio \([kg\; kg^{-1}]\)
    real(wp), dimension(kts:kte) :: qg1d !! 1D graupel mass mixing ratio \([kg\; kg^{-1}]\)
    real(wp), dimension(kts:kte) :: ni1d !! 1D cloud ice number mixing ratio \([kg^{-1}]\)
    real(wp), dimension(kts:kte) :: nr1d !! 1D rain water number mixing ratio \([kg^{-1}]\)
    real(wp), dimension(kts:kte) :: w1d !! 1D vertical velocity \(m\; s^{-1}]\)
    real(wp), dimension(kts:kte) :: dz1d !! 1D vertical grid spacing \([m]\)

    real(wp), dimension(:), allocatable :: nc1d !! 1D cloud water number mixing ratio \([kg^{-1}]\) (aerosol-aware)
    real(wp), dimension(:), allocatable :: nwfa1d !! 1D water-friendly aerosol number mixing ratio \([kg^{-1}]\) (aerosol-aware)
    real(wp), dimension(:), allocatable :: nifa1d !! 1D ice-friendly aerosol number mixing ratio \([kg^{-1}]\) (aerosol-aware)
    real(wp), dimension(:), allocatable :: qb1d !! 1D graupel volume mixing ratio \([m^{-3}\; kg^{-1}]\) (hail-aware)
    real(wp), dimension(:), allocatable :: ng1d !! 1D graupel number mixing ratio \([kg^{-1}]\) (hail-aware)
    integer :: i, j, k, nz
    logical :: use_temperature 

    type(ty_tempo_main_diags) :: tempo_main_diags
    type(ty_tempo_driver_diags), intent(out) :: tempo_diags

    nz = kte - kts + 1
    ! allocate 1d arrays if 3d arrays are present
    if (present(nwfa)) allocate(nwfa1d(nz), source=0._wp)
    if (present(nifa)) allocate(nifa1d(nz), source=0._wp)
    if (present(nc)) allocate(nc1d(nz), source=0._wp)
    if (present(ng)) allocate(ng1d(nz), source=0._wp)
    if (present(qb)) allocate(qb1d(nz), source=0._wp)

    ! allocate diagnostics
    ! 3d diagnostics have configuration flags
    if (tempo_cfgs%rain_med_vol_diam) allocate(tempo_diags%rain_med_vol_diam(its:ite, kts:kte, jts:jte), source=0._wp)
    if (tempo_cfgs%graupel_med_vol_diam) allocate(tempo_diags%graupel_med_vol_diam(its:ite, kts:kte, jts:jte), source=0._wp)
    if (tempo_cfgs%refl10cm) allocate(tempo_diags%refl10cm(its:ite, kts:kte, jts:jte), source=-35._wp)
    if (tempo_cfgs%re_cloud) allocate(tempo_diags%re_cloud(its:ite, kts:kte, jts:jte), source=0._wp)
    if (tempo_cfgs%re_ice) allocate(tempo_diags%re_ice(its:ite, kts:kte, jts:jte), source=0._wp)
    if (tempo_cfgs%re_snow) allocate(tempo_diags%re_snow(its:ite, kts:kte, jts:jte), source=0._wp)
    if (tempo_cfgs%max_hail_diameter) allocate(tempo_diags%max_hail_diameter(its:ite, kts:kte, jts:jte), source=0._wp)

    ! precipitation
    allocate(tempo_diags%rain_precip(its:ite, jts:jte), source=0._wp)
    allocate(tempo_diags%ice_liquid_equiv_precip(its:ite, jts:jte), source=0._wp)
    allocate(tempo_diags%snow_liquid_equiv_precip(its:ite, jts:jte), source=0._wp)
    allocate(tempo_diags%graupel_liquid_equiv_precip(its:ite, jts:jte), source=0._wp)
    allocate(tempo_diags%frozen_fraction(its:ite, jts:jte), source=0._wp)

    ! temperature or theta and exner
    if (present(t)) then
      use_temperature = .true.
    elseif (present(th) .and. present(pii)) then
      use_temperature = .false.
    else  
      error stop "tempo_driver() --- requires either temperature or theta and Exner function"
    endif 

    ! tempo driver code
    do j = jts, jte
      do i = its, ite
        do k = kts, kte
          if (use_temperature) then
            t1d(k) = t(i,k,j)
          else  
            t1d(k) = th(i,k,j) * pii(i,k,j)
          endif 
          p1d(k) = p(i,k,j)
          w1d(k) = w(i,k,j)
          dz1d(k) = dz(i,k,j)
          qv1d(k) = qv(i,k,j)
          qc1d(k) = qc(i,k,j)
          qi1d(k) = qi(i,k,j)
          qr1d(k) = qr(i,k,j)
          qs1d(k) = qs(i,k,j)
          qg1d(k) = qg(i,k,j)
          ni1d(k) = ni(i,k,j)
          nr1d(k) = nr(i,k,j)

          ! nwfa, nifa, and nc are optional aerosol-aware variables
          if (present(nwfa)) nwfa1d(k) = nwfa(i,k,j)
          if (present(nifa)) nifa1d(k) = nifa(i,k,j)
          if (present(nc)) nc1d(k) = nc(i,k,j)

          ! ng and qb are optional hail-aware variables
          if ((present(ng)) .and. (present(qb))) then
            ng1d(k) = ng(i,k,j)
            qb1d(k) = qb(i,k,j)
          endif 
        enddo
          
        ! main call to the 1d tempo microphysics
        call tempo_main(qv1d=qv1d, qc1d=qc1d, qi1d=qi1d, qr1d=qr1d, qs1d=qs1d, qg1d=qg1d, qb1d=qb1d, &
          ni1d=ni1d, nr1d=nr1d, nc1d=nc1d, ng1d=ng1d, nwfa1d=nwfa1d, nifa1d=nifa1d, t1d=t1d, p1d=p1d, &
          w1d=w1d, dz1d=dz1d, kts=kts, kte=kte, dt=dt, ii=i, jj=j, tempo_main_diags=tempo_main_diags)
          
        ! precipitation
        tempo_diags%rain_precip(i,j) = tempo_main_diags%rain_precip
        tempo_diags%ice_liquid_equiv_precip(i,j) = tempo_main_diags%ice_liquid_equiv_precip
        tempo_diags%snow_liquid_equiv_precip(i,j) = tempo_main_diags%snow_liquid_equiv_precip
        tempo_diags%graupel_liquid_equiv_precip(i,j) = tempo_main_diags%graupel_liquid_equiv_precip
        tempo_diags%frozen_fraction(i,j) = tempo_main_diags%frozen_fraction

        ! 3d diagnostics
        if (allocated(tempo_diags%rain_med_vol_diam) .and. allocated(tempo_main_diags%rain_med_vol_diam)) then
          tempo_diags%rain_med_vol_diam(i,:,j) = tempo_main_diags%rain_med_vol_diam
        endif 
        if (allocated(tempo_diags%graupel_med_vol_diam) .and. allocated(tempo_main_diags%graupel_med_vol_diam)) then
          tempo_diags%graupel_med_vol_diam(i,:,j) = tempo_main_diags%graupel_med_vol_diam
        endif 
        if (allocated(tempo_diags%re_cloud) .and. allocated(tempo_main_diags%re_cloud)) then
          tempo_diags%re_cloud(i,:,j) = tempo_main_diags%re_cloud
        endif 
        if (allocated(tempo_diags%re_ice) .and. allocated(tempo_main_diags%re_ice)) then
          tempo_diags%re_ice(i,:,j) = tempo_main_diags%re_ice
        endif 
        if (allocated(tempo_diags%re_snow) .and. allocated(tempo_main_diags%re_snow)) then
          tempo_diags%re_snow(i,:,j) = tempo_main_diags%re_snow
        endif 
        if (allocated(tempo_diags%refl10cm) .and. allocated(tempo_main_diags%refl10cm)) then
          tempo_diags%refl10cm(i,:,j) = tempo_main_diags%refl10cm
        endif
        if (allocated(tempo_diags%max_hail_diameter) .and. allocated(tempo_main_diags%max_hail_diameter)) then
          tempo_diags%max_hail_diameter(i,:,j) = tempo_main_diags%max_hail_diameter
        endif 
        ! return variables to model
        do k = kts, kte
          if (present(nc)) nc(i,k,j) = nc1d(k)
          if (present(nwfa)) nwfa(i,k,j) = nwfa1d(k)
          if (present(nifa)) nifa(i,k,j) = nifa1d(k)
          if ((present(ng)) .and. (present(qb))) then
            ng(i,k,j) = ng1d(k)
            qb(i,k,j) = qb1d(k)
          endif 
          qv(i,k,j) = qv1d(k)
          qc(i,k,j) = qc1d(k)
          qi(i,k,j) = qi1d(k)
          qr(i,k,j) = qr1d(k)
          qs(i,k,j) = qs1d(k)
          qg(i,k,j) = qg1d(k)
          ni(i,k,j) = ni1d(k)
          nr(i,k,j) = nr1d(k)
          
          ! tempo main returns temperature (t1d), so convert to theta if needed
          if (use_temperature) then
            t(i,k,j) = t1d(k)
          else  
            th(i,k,j) = t1d(k) / pii(i,k,j)
          endif 
        enddo
      enddo
    enddo
  end subroutine tempo_driver


  subroutine tempo_aerosol_surface_emissions(dt, nwfa, nwfa2d, ims, ime, jms, jme, kms, kme, kts)
    !! adds aerosol surface emissions to the 3D field
    real(wp), intent(in) :: dt
    real(wp), dimension(ims:ime, kms:kme, jms:jme), intent(inout) :: nwfa 
    real(wp), dimension(ims:ime, jms:jme), intent(in) :: nwfa2d
    integer, intent(in) :: ims, ime, jms, jme, kms, kme, kts
    integer :: i, j, k

    do j = jms, jme
      do i = ims, ime
        nwfa(i,kts,j) = nwfa(i,kts,j) + nwfa2d(i,j) * dt
      enddo
    enddo
  end subroutine tempo_aerosol_surface_emissions

end module module_mp_tempo_driver
