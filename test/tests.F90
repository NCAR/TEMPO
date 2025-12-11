module tests
  !! Tests for TEMPO microphysics
  use module_mp_tempo_params, only : tempo_init_cfgs
  use module_mp_tempo_init, only : tempo_init
  ! use module_mp_tempo_driver, only : tempo_driver, tempo_driver_cfgs
  implicit none
  private

  public :: test_tempo_init !, test_tempo_driver

  contains

  subroutine test_tempo_init()
    !! test tempo initialization procedure
    call tempo_init()
  end subroutine test_tempo_init


  ! subroutine test_tempo_driver()
  !   use module_mp_tempo_params, only : wp, sp, dp

  !   integer :: itimestep
  !   real(wp) :: dt
  !   integer, parameter :: nz = 50
  !   integer, parameter :: ids = 1, ide = 1, ims = 1, ime = 1, its = 1, ite = 1
  !   integer, parameter :: jds = 1, jde = 1, jms = 1, jme = 1, jts = 1, jte = 1
  !   integer, parameter :: kds = 1, kde = nz, kms = 1, kme = nz, kts = 1, kte = nz
  !   real(wp), dimension(its:ite, kts:kte, jts:jte) :: qv, t, th, pii, p, w, dz, &
  !     qc, qr, qi, qs, qg, ni, nr
        
  !   dt = 1.
  !   qv = 1.e-3_wp
  !   t = 273.15_wp
  !   p = 100000._wp
  !   w = 0.
  !   dz = 100.
  !   qc = 1.e-3_wp
  !   qr = 0.
  !   qi = 0.
  !   qs = 0.
  !   qg = 0.
  !   ni = 0.
  !   nr = 0.

  !   do itimestep = 1, 100

  !     call  tempo_driver(itimestep=itimestep, dt=dt, &
  !                       ids=ids, ide=ide, ims=ims, ime=ime, its=its, ite=ite, &
  !                       jds=jds, jde=jde, jms=jms, jme=jme, jts=jts, jte=jte, &
  !                       kds=kds, kde=kde, kms=kms, kme=kme, kts=kts, kte=kte, &
  !                       t=t, p=p, w=w, dz=dz, qv=qv, &
  !                       qc=qc, qr=qr, qi=qi, qs=qs, qg=qg, ni=ni, nr=nr)
  !   enddo
  ! end subroutine test_tempo_driver
end module tests
