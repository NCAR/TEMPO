module tests
  !! Tests for TEMPO microphysics
  use module_mp_tempo_params, only : tempo_cfgs
  use module_mp_tempo_init, only : tempo_init
  use module_mp_tempo_driver, only : tempo_driver, ty_tempo_driver_diags
  implicit none
  private

  public :: test_tempo_init, test_tempo_driver, test_rain_sedimentation, &
    test_graupel_sedimentation

  contains

  subroutine test_tempo_init()
    !! test tempo initialization procedure
    call tempo_init()
  end subroutine test_tempo_init


 subroutine test_graupel_sedimentation(dt)
    use module_mp_tempo_params, only : wp, sp, dp

    integer :: itimestep
    real(wp), intent(in) :: dt
    integer, parameter :: nz = 100
    integer, parameter :: ids = 1, ide = 1, ims = 1, ime = 1, its = 1, ite = 1
    integer, parameter :: jds = 1, jde = 1, jms = 1, jme = 1, jts = 1, jte = 1
    integer, parameter :: kds = 1, kde = nz, kms = 1, kme = nz, kts = 1, kte = nz
    real(wp), dimension(its:ite, kts:kte, jts:jte) :: qv, t, th, pii, p, w, dz, &
      qc, qr, qi, qs, qg, ni, nr, ng, qb
    type(ty_tempo_driver_diags) :: tempo_driver_diags
    real(wp) :: zsum, psum
    character(len=20) :: string
    integer :: io, io3, k, total_timesteps
    write(string, '(I5)') int(dt)

    qv = 10.e-3_wp
    t = 273.15_wp
    p = 100000._wp
    w = 0._wp
    ! dz = 70._wp
    qc = 0._wp
    qr = 0._wp
    nr = 0._wp
    qi = 0._wp
    qs = 0._wp
    qg = 0._wp
    ni = 0._wp
    qb = 0._wp

    zsum = 0.
    do k = 1, nz
      dz(:,k,:) = k*10.
      zsum = zsum + dz(1,k,1)
      if ((k > 30) .and. (k < 60)) then
        qg(:,k,:) = 1.e-3
        ng(:,k,:) = 2000.
        qb(:,k,:)= 1000.*qg(:,k,:)/(real(k)*15.)
        write(*,*), k, (600. - (real(k)-30.)*10. )
      endif 
    enddo 
    
    tempo_cfgs%all_processes_off = .true.

    total_timesteps = int(2400/dt)
    !total_timesteps = int(1200/dt)
    !total_timesteps = int(240/dt)

    psum = 0.
    open(newunit=io3, file="precip_"//trim(adjustl(string))//".txt", status="new", action="write")

    do itimestep = 1, total_timesteps

      call  tempo_driver(itimestep=itimestep, dt=dt, &
                        ids=ids, ide=ide, ims=ims, ime=ime, its=its, ite=ite, &
                        jds=jds, jde=jde, jms=jms, jme=jme, jts=jts, jte=jte, &
                        kds=kds, kde=kde, kms=kms, kme=kme, kts=kts, kte=kte, &
                        t=t, p=p, w=w, dz=dz, qv=qv, th=th, pii=pii, &
                        qc=qc, qr=qr, qi=qi, qs=qs, qg=qg, ng=ng, qb=qb, ni=ni, nr=nr, &
                        tempo_driver_diags=tempo_driver_diags)
      psum = psum + tempo_driver_diags%graupel_liquid_equiv_precipitation(1,1)
      write(io3,*) itimestep*dt, psum

      if (itimestep == 1) then
        open(newunit=io, file="test_graupel_sedimentation_dt_"//trim(adjustl(string))//"_t0.txt", status="new", action="write")
        do k = 1, nz
          write(io,'(I5, 3E12.4)') k, qg(1,k,1), ng(1,k,1), 1000.*qg(1,k,1)/qb(1,k,1)
        enddo 
      endif 
    enddo

    open(newunit=io, file="test_graupel_sedimentation_dt_"//trim(adjustl(string))//"_t1.txt", status="new", action="write")
    do k = 1, nz
      write(io,'(I5, 3E12.4)') k, qg(1,k,1), ng(1,k,1), 1000.*qg(1,k,1)/qb(1,k,1)
    enddo 
  end subroutine test_graupel_sedimentation


 subroutine test_rain_sedimentation(dt)
    use module_mp_tempo_params, only : wp, sp, dp

    integer :: itimestep
    real(wp), intent(in) :: dt
    integer, parameter :: nz = 100
    integer, parameter :: ids = 1, ide = 1, ims = 1, ime = 1, its = 1, ite = 1
    integer, parameter :: jds = 1, jde = 1, jms = 1, jme = 1, jts = 1, jte = 1
    integer, parameter :: kds = 1, kde = nz, kms = 1, kme = nz, kts = 1, kte = nz
    real(wp), dimension(its:ite, kts:kte, jts:jte) :: qv, t, th, pii, p, w, dz, &
      qc, qr, qi, qs, qg, ni, nr
    type(ty_tempo_driver_diags) :: tempo_driver_diags
    real(wp) :: zsum, psum
    character(len=20) :: string
    integer :: io, io3, k, total_timesteps
    write(string, '(I5)') int(dt)

    qv = 10.e-3_wp
    t = 273.15_wp
    p = 100000._wp
    w = 0._wp
    ! dz = 70._wp
    qc = 0._wp
    qr = 0._wp
    nr = 0._wp
    qi = 0._wp
    qs = 0._wp
    qg = 0._wp
    ni = 0._wp
      
    zsum = 0.
    do k = 1, nz
      dz(:,k,:) = k*10.
      zsum = zsum + dz(1,k,1)
      ! write(*,*) k, dz(1,k,1), zsum
      if ((k > 30) .and. (k < 60)) then
        qr(:,k,:) = 1.e-3
        nr(:,k,:) = 2000.
      endif 
    enddo 
    
    tempo_cfgs%all_processes_off = .true.

    total_timesteps = int(2400/dt)
    ! total_timesteps = int(1200/dt)
    psum = 0.
    open(newunit=io3, file="precip_"//trim(adjustl(string))//".txt", status="new", action="write")

    do itimestep = 1, total_timesteps

      call  tempo_driver(itimestep=itimestep, dt=dt, &
                        ids=ids, ide=ide, ims=ims, ime=ime, its=its, ite=ite, &
                        jds=jds, jde=jde, jms=jms, jme=jme, jts=jts, jte=jte, &
                        kds=kds, kde=kde, kms=kms, kme=kme, kts=kts, kte=kte, &
                        t=t, p=p, w=w, dz=dz, qv=qv, th=th, pii=pii, &
                        qc=qc, qr=qr, qi=qi, qs=qs, qg=qg, ni=ni, nr=nr, &
                        tempo_driver_diags=tempo_driver_diags)
      psum = psum + tempo_driver_diags%rain_precipitation(1,1)
      write(io3,*) itimestep*dt, psum

      if (itimestep == 1) then
        open(newunit=io, file="test_rain_sedimentation_dt_"//trim(adjustl(string))//"_t0.txt", status="new", action="write")
        do k = 1, nz
          write(io,'(I5, 3E12.4)') k, qr(1,k,1), nr(1,k,1), tempo_driver_diags%mvd_r(1,k,1)
        enddo 
      endif 
    enddo

    open(newunit=io, file="test_rain_sedimentation_dt_"//trim(adjustl(string))//"_t1.txt", status="new", action="write")
    do k = 1, nz
      write(io,'(I5, 3E12.4)') k, qr(1,k,1), nr(1,k,1), tempo_driver_diags%mvd_r(1,k,1)
    enddo 
  end subroutine test_rain_sedimentation


  subroutine test_tempo_driver()
    use module_mp_tempo_params, only : wp, sp, dp

    integer :: itimestep
    real(wp) :: dt
    integer, parameter :: nz = 100
    integer, parameter :: ids = 1, ide = 1, ims = 1, ime = 1, its = 1, ite = 1
    integer, parameter :: jds = 1, jde = 1, jms = 1, jme = 1, jts = 1, jte = 1
    integer, parameter :: kds = 1, kde = nz, kms = 1, kme = nz, kts = 1, kte = nz
    real(wp), dimension(its:ite, kts:kte, jts:jte) :: qv, t, th, pii, p, w, dz, &
      qc, qr, qi, qs, qg, ni, nr
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
                        tempo_driver_diags=tempo_driver_diags)
    enddo
    do k = 1, nz
      write(19,*) k, tempo_driver_diags%reflectivity(1,k,1) ! qr(1,k,1), nr(1,k,1), qs(1,k,1)
    enddo 

    !  do k = 1, nz
    !   if (qsave(k) <= r1) then
    !     if (qcten(k) == 0._wp) then
    !       if (qc1d(k) > 0._wp) then
    !         write(*,*) 'ERROR', k, qsave(k), qc1d(k), qcten(k)
    !       endif 
    !     endif 
    !   endif 
    !   if ((qsave(k) <= r1) .and. (qc1d(k) > r1)) then
    !     if (qcten(k) == 0._wp) then
    !       write(*,*) 'ERROR', k, qsave(k), qc1d(k), qcten(k)
    !     endif 
    !   endif 
    ! enddo 
  end subroutine test_tempo_driver
end module tests