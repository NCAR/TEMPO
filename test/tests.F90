module tests
  !! Tests for TEMPO microphysics
  use module_mp_tempo_params, only : tempo_init_cfgs
  use module_mp_tempo_init, only : tempo_init
  use module_mp_tempo_driver, only : tempo_driver, ty_tempo_driver_diags !tempo_driver_cfgs
  implicit none
  private

  public :: test_tempo_init, test_tempo_driver

  contains

  subroutine test_tempo_init()
    !! test tempo initialization procedure
    call tempo_init()
  end subroutine test_tempo_init


  subroutine test_tempo_driver()
    use module_mp_tempo_params, only : wp, sp, dp

    integer :: itimestep
    real(wp) :: dt
    integer, parameter :: nz = 100
    integer, parameter :: ids = 1, ide = 1, ims = 1, ime = 1, its = 1, ite = 1
    integer, parameter :: jds = 1, jde = 1, jms = 1, jme = 1, jts = 1, jte = 1
    integer, parameter :: kds = 1, kde = nz, kms = 1, kme = nz, kts = 1, kte = nz
    real(wp), dimension(its:ite, kts:kte, jts:jte) :: qv, t, th, pii, p, w, dz, &
      qc, qr, qi, qs, qg, ni, nr, re_cloud, re_ice, re_snow
    real(wp), dimension(its:ite, jts:jte) :: rainnc, rainncv, sr
    type(ty_tempo_driver_diags) :: tempo_driver_diags
    integer :: k

    dt = 20.
    qv = 10.e-3_wp
    t = 273.15_wp
    th = 273.15
    pii = 1.
    p = 100000._wp
    w = 0.
    dz = 100.
    qc = 1.e-3_wp
    qr = 0.
    nr = 0.
    qi = 0.
    qs = 0.
    qg = 0.
    ni = 0.
    re_cloud = 0.
    re_snow = 0.
    re_ice = 0.
    rainnc = 0.
    rainncv = 0.
    sr = 0.
    
    do k = 1, nz
      if (k < 50) then
        !qr(:,k,:) = 0.
        !nr(:,k,:) = 0.
        th(:,k,:) = 276.
      else
        qs(:,k,:) = 3.e-3
        !qr(:,k,:) = 1.e-3
        !nr(:,k,:) = 2000.
        th(:,k,:) = 273.
      endif 
      ! if ((k >= 50) .and. (k < 60)) then
      !   th(:,k,:) = 235.15
      ! endif 
    enddo 
    
    do k = 1, nz
      write(96,*) k, qs(1,k,1)
    enddo 

    do itimestep = 1, 100
      call  tempo_driver(itimestep=itimestep, dt=dt, &
                        ids=ids, ide=ide, ims=ims, ime=ime, its=its, ite=ite, &
                        jds=jds, jde=jde, jms=jms, jme=jme, jts=jts, jte=jte, &
                        kds=kds, kde=kde, kms=kms, kme=kme, kts=kts, kte=kte, &
                        t=t, p=p, w=w, dz=dz, qv=qv, th=th, pii=pii, &
                        qc=qc, qr=qr, qi=qi, qs=qs, qg=qg, ni=ni, nr=nr, &
                        rainnc=rainnc, rainncv=rainncv, sr=sr, re_cloud=re_cloud, &
                        re_ice=re_ice, re_snow=re_snow, tempo_driver_diags=tempo_driver_diags)
    enddo
    do k = 1, nz
      write(19,*) k, tempo_driver_diags%reflectivity3d(1,k,1) ! qr(1,k,1), nr(1,k,1), qs(1,k,1)
    enddo 
  end subroutine test_tempo_driver
end module tests