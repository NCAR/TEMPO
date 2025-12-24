module module_mp_tempo_driver

  use module_mp_tempo_params, only : wp, sp, dp
  use module_mp_tempo_main, only : tempo_main, ty_tempo_main_diags
  use module_mp_tempo_ml, only : predict_number_sub

  implicit none
  private

  public :: tempo_driver, ty_tempo_driver_cfgs, ty_tempo_driver_diags
  
  ! tempo configuration flags for the driver
  type :: ty_tempo_driver_cfgs
    logical :: semi_sedi = .false. !! flag for semi-lagrangian sedimentation
  end type

  type :: ty_tempo_driver_diags
    real(wp), allocatable, dimension(:,:) :: rain_precipitation
    real(wp), allocatable, dimension(:,:) :: snow_liquid_equiv_precipitation
    real(wp), allocatable, dimension(:,:,:) :: reflectivity
  end type

  contains

  subroutine tempo_driver(dt, itimestep, &
    t, th, pii, p, w, dz, &
    qv, qc, qr, qi, qs, qg, ni, nr, &
    nc, nwfa, nifa, ng, qb, &
    nwfa2d, nifa2d, &
    ids, ide, jds, jde, kds, kde, &
    ims, ime, jms, jme, kms, kme, &
    its, ite, jts, jte, kts, kte, tempo_driver_cfgs, tempo_driver_diags)

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
    real, dimension(ims:ime, jms:jme), intent(in), optional :: nwfa2d, nifa2d

    !real, dimension(ims:ime, kms:kme, jms:jme), intent(inout) :: re_cloud, re_ice, re_snow
    !real, dimension(ims:ime, jms:jme), intent(inout) :: rainnc, rainncv, sr
    !real, optional, dimension(ims:ime,jms:jme), intent(inout) :: frainnc, max_hail_diameter_column, max_hail_diameter_sfc
    !real, dimension(ims:ime, kms:kme, jms:jme), intent(inout), optional :: refl_10cm
    !real, dimension(ims:ime, kms:kme, jms:jme), intent(in), optional :: qcbl, cldfrac
    !real, dimension(ims:ime, jms:jme), intent(inout), optional :: snownc, snowncv, graupelnc, graupelncv

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

    real(wp), dimension(:), allocatable :: nc1d !! 1D cloud water number mixing ratio \([kg^{-1}]\)
    real(wp), dimension(:), allocatable :: nwfa1d !! 1D water-friendly aerosol number mixing ratio \([kg^{-1}]\)
    real(wp), dimension(:), allocatable :: nifa1d !! 1D ice-friendly aerosol number mixing ratio \([kg^{-1}]\)
    real(wp), dimension(:), allocatable :: qb1d !! 1D graupel volume mixing ratio \([m^{-3}\; kg^{-1}]\)
    real(wp), dimension(:), allocatable :: ng1d !! 1D graupel number mixing ratio \([kg^{-1}]\)

    
    ! Local (1d) variables
    ! real, dimension(kts:kte) :: rho, dbz, qcbl1d, cldfrac1d, qg_max_diam1d
    ! real, dimension(kts:kte) :: re_qc1d, re_qi1d, re_qs1d
    ! double precision, dimension(kts:kte) :: ncbl1d
    ! real, dimension(its:ite, jts:jte) :: pcp_ra, pcp_sn, pcp_gr, pcp_ic, frain
    ! real :: pptrain, pptsnow, pptgraul, pptice
    !real :: nwfa1
    !real :: ygra1, zans1
    !real :: graupel_vol
    !real :: tmprc, tmpnc, xDc
    !integer :: nu_c
    !logical, dimension(kts:kte) :: sgs_clouds
    !double precision :: lamg, lam_exp, lamr, n0_min, n0_exp, lamc
    integer :: i, j, k, nz
    !logical, optional, intent(in) :: diagflag
    !integer, optional, intent(in) :: do_radar_ref
    !character(len=132) :: message

    type(ty_tempo_main_diags) :: tempo_main_diags
    type(ty_tempo_driver_diags), intent(out) :: tempo_driver_diags
    type(ty_tempo_driver_cfgs), intent(in) :: tempo_driver_cfgs

    nz = kte - kts + 1
    allocate(tempo_driver_diags%reflectivity(its:ite, kts:kte, jts:jte), source=-35._wp)
    allocate(tempo_driver_diags%rain_precipitation(its:ite, jts:jte), source=0._wp)
    allocate(tempo_driver_diags%snow_liquid_equiv_precipitation(its:ite, jts:jte), source=0._wp)

    if (present(nwfa)) allocate(nwfa1d(nz), source=0._wp)
    if (present(nifa)) allocate(nifa1d(nz), source=0._wp)
    if (present(nc)) allocate(nc1d(nz), source=0._wp)
    if (present(ng)) allocate(ng1d(nz), source=0._wp)
    if (present(qb)) allocate(qb1d(nz), source=0._wp)
    do j = jts, jte
      do i = its, ite
        ! Begin k loop
        do k = kts, kte
          t1d(k) = th(i,k,j) * pii(i,k,j)
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
          if (present(nwfa)) then
            if (present(nwfa2d)) then
              if (k == kts) then
                nwfa(i,k,j) = nwfa(i,k,j) + nwfa2d(i,j) * dt
              endif
            endif
            nwfa1d(k) = nwfa(i,k,j)
          endif

          if (present(nifa)) nifa1d(k) = nifa(i,k,j)
          if (present(nc)) nc1d(k) = nc(i,k,j)

          ! ng and qb are optional hail-aware variables
          if ((present(ng)) .and. (present(qb))) then
            ng1d(k) = ng(i,k,j)
            qb1d(k) = qb(i,k,j)
          endif 
        enddo
          
        ! Main call to the 1D microphysics
        call tempo_main(qv1d=qv1d, qc1d=qc1d, qi1d=qi1d, qr1d=qr1d, qs1d=qs1d, &
          qg1d=qg1d, qb1d=qb1d, &
          ni1d=ni1d, nr1d=nr1d, nc1d=nc1d, ng1d=ng1d, nwfa1d=nwfa1d, nifa1d=nifa1d, t1d=t1d, p1d=p1d, w1d=w1d, dz1d=dz1d, kts=kts, kte=kte, dt=dt, ii=i, jj=j, tempo_main_diags=tempo_main_diags, itimestep=itimestep)
      
        tempo_driver_diags%rain_precipitation(i,j) = tempo_main_diags%rain_precipitation
        tempo_driver_diags%snow_liquid_equiv_precipitation(i,j) = &
          tempo_main_diags%snow_liquid_equiv_precipitation

        if (allocated(tempo_main_diags%reflectivity)) then
          tempo_driver_diags%reflectivity(i,:,j) = tempo_main_diags%reflectivity
        endif

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
          th(i,k,j) = t1d(k) / pii(i,k,j)
        enddo
      enddo
    enddo 
  end subroutine tempo_driver

end module module_mp_tempo_driver
