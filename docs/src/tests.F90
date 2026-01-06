module tests
  !! TEMPO tests
  use module_mp_tempo_cfgs, only : ty_tempo_cfgs
  use module_mp_tempo_init, only : tempo_init
  use module_mp_tempo_driver, only : tempo_driver, ty_tempo_driver_diags
  implicit none
  private

  public :: test_tempo_init, test_graupel_sedimentation, test_snow_sedimentation

  contains

  subroutine test_tempo_init()
    !! test tempo initialization procedure
    call tempo_init()
  end subroutine test_tempo_init


 subroutine test_graupel_sedimentation(dt, semi_sedi)
    !! test graupel sedimentation
    use module_mp_tempo_params, only : wp, sp, dp

    integer :: itimestep
    real(wp), intent(in) :: dt
    logical, intent(in) :: semi_sedi
    integer, parameter :: nz = 59
    integer, parameter :: ids = 1, ide = 1, ims = 1, ime = 1, its = 1, ite = 1
    integer, parameter :: jds = 1, jde = 1, jms = 1, jme = 1, jts = 1, jte = 1
    integer, parameter :: kds = 1, kde = nz, kms = 1, kme = nz, kts = 1, kte = nz
    real(wp), dimension(its:ite, kts:kte, jts:jte) :: qv, t, th, pii, p, w, dz, &
      qc, qr, qi, qs, qg, ni, nr, ng, qb
    real(wp), dimension(nz) :: klevs_in, qv_in, qc_in, qr_in, qi_in, qs_in, &
      qg_in, ni_in, nr_in, nc_in, nwfa_in, nifa_in, theta_in, ng_in, volg_in, &
      pressure_in, w_in, dz_in
    real(wp) :: precip_sum
    character(len=20) :: dt_string, tt_string
    character(len=20) :: semi_sedi_string
    type(ty_tempo_driver_diags) :: tempo_driver_diags
    type(ty_tempo_cfgs) :: tempo_cfgs
    integer :: io1, io2, io3, io4, k, total_timesteps
    integer, parameter :: integration_time = 1200

    write(dt_string, '(I7)') int(dt)
    
    if (semi_sedi) then
      write(semi_sedi_string, '(A)') '_semi_sedi'
    else
      write(semi_sedi_string, '(A)') ''
    endif

    ! read input file
    io1 = 11
    open (io1, file='./test/data/mpas_59lev_test.txt', status='old')
    read(io1,*) ! header
    do k = 1, nz
      read(io1,*) klevs_in(k), qv_in(k), qc_in(k), qr_in(k), qi_in(k), qs_in(k), &
          qg_in(k), ni_in(k), nr_in(k), nc_in(k), nwfa_in(k), nifa_in(k), &
          theta_in(k), ng_in(k), volg_in(k), pressure_in(k), w_in(k), dz_in(k)
    end do
    close(io1)

    qv(1,:,1) = qv_in
    t(1,:,1) = theta_in * (pressure_in/100000.)**0.286
    p(1,:,1) = pressure_in
    dz(1,:,1) = dz_in
    w = 0._wp
    qc = 0._wp
    qr = 0._wp
    nr = 0._wp
    qi = 0._wp
    qs = 0._wp
    ni = 0._wp
    qg(1,:,1) = qg_in
    ng(1,:,1) = ng_in
    qb(1,:,1) = volg_in * 1000. ! convert meters^3 -> liters

    ! set configs
    tempo_cfgs%turn_off_micro_flag = .true.
    tempo_cfgs%graupel_med_vol_diam_flag = .true.
    tempo_cfgs%semi_sedi_flag = semi_sedi

    total_timesteps = int(integration_time/dt)
    precip_sum = 0._wp
    open(newunit=io2, file="graupel_precip_dt_"//trim(adjustl(dt_string))//""//trim(semi_sedi_string)//".txt", &
      status="new", action="write")

    do itimestep = 1, total_timesteps
      call  tempo_driver(tempo_cfgs=tempo_cfgs, itimestep=itimestep, dt=dt, &
                        ids=ids, ide=ide, ims=ims, ime=ime, its=its, ite=ite, &
                        jds=jds, jde=jde, jms=jms, jme=jme, jts=jts, jte=jte, &
                        kds=kds, kde=kde, kms=kms, kme=kme, kts=kts, kte=kte, &
                        t=t, p=p, w=w, dz=dz, qv=qv, th=th, pii=pii, &
                        qc=qc, qr=qr, qi=qi, qs=qs, qg=qg,  ni=ni, nr=nr, ng=ng, qb=qb, &
                        tempo_diags=tempo_driver_diags)
      precip_sum = precip_sum + tempo_driver_diags%graupel_liquid_equiv_precip(1,1)
      write(io2,'(I7, 1E12.4)') int(itimestep*dt), precip_sum

      if (dt == 1. .and. .not. semi_sedi .and. itimestep == 1) then
        open(newunit=io3, file="graupel_sedi_init.txt", status="new", action="write")
        do k = 1, nz
          write(io3,'(I5, 4E12.4)') k, qg(1,k,1), ng(1,k,1), 1000.*qg(1,k,1)/qb(1,k,1), &
            tempo_driver_diags%graupel_med_vol_diam(1,k,1)
        enddo 
      endif 
    enddo

    write(tt_string, '(I7)') integration_time
    open(newunit=io4, file="graupel_sedi_dt_"//trim(adjustl(dt_string))//""//trim(semi_sedi_string)//&
      "_runtime"//trim(adjustl(tt_string))//".txt", status="new", action="write")
    do k = 1, nz
      write(io4,'(I5, 4E12.4)') k, qg(1,k,1), ng(1,k,1), 1000.*qg(1,k,1)/qb(1,k,1), &
        tempo_driver_diags%graupel_med_vol_diam(1,k,1)
    enddo 
  end subroutine test_graupel_sedimentation


  subroutine test_snow_sedimentation(dt)
    use module_mp_tempo_params, only : wp, sp, dp

    integer :: itimestep
    real(wp), intent(in) :: dt
    integer, parameter :: nz = 59
    integer, parameter :: ids = 1, ide = 1, ims = 1, ime = 1, its = 1, ite = 1
    integer, parameter :: jds = 1, jde = 1, jms = 1, jme = 1, jts = 1, jte = 1
    integer, parameter :: kds = 1, kde = nz, kms = 1, kme = nz, kts = 1, kte = nz
    real(wp), dimension(its:ite, kts:kte, jts:jte) :: qv, t, th, pii, p, w, dz, &
      qc, qr, qi, qs, qg, ni, nr, ng, qb
    real(wp), dimension(nz) :: klevs_in, qv_in, qc_in, qr_in, qi_in, qs_in, &
      qg_in, ni_in, nr_in, nc_in, nwfa_in, nifa_in, theta_in, ng_in, volg_in, &
      pressure_in, w_in, dz_in
    real(wp) :: precip_sum
    character(len=20) :: dt_string, tt_string
    type(ty_tempo_driver_diags) :: tempo_driver_diags
    type(ty_tempo_cfgs) :: tempo_cfgs
    integer :: io1, io2, io3, io4, k, total_timesteps
    integer, parameter :: integration_time = 1200

    write(dt_string, '(I7)') int(dt)


    ! read input file
    io1 = 11
    open (io1, file='./test/data/mpas_59lev_test.txt', status='old')
    read(io1,*) ! header
    do k = 1, nz
      read(io1,*) klevs_in(k), qv_in(k), qc_in(k), qr_in(k), qi_in(k), qs_in(k), &
          qg_in(k), ni_in(k), nr_in(k), nc_in(k), nwfa_in(k), nifa_in(k), &
          theta_in(k), ng_in(k), volg_in(k), pressure_in(k), w_in(k), dz_in(k)
    end do
    close(io1)

    qv(1,:,1) = qv_in
    t(1,:,1) = theta_in * (pressure_in/100000.)**0.286
    p(1,:,1) = pressure_in
    dz(1,:,1) = dz_in
    w = 0._wp
    qc = 0._wp
    qr = 0._wp
    nr = 0._wp
    qi = 0._wp
    qs(1,:,1) = qs_in
    ni = 0._wp
    qg = 0._wp
    ng = 0._wp
    qb = 0._wp

    ! set configs
    tempo_cfgs%turn_off_micro_flag= .true.

    total_timesteps = int(integration_time/dt)
    precip_sum = 0._wp
    open(newunit=io2, file="snow_precip_dt_"//trim(adjustl(dt_string))//".txt", &
      status="new", action="write")

    do itimestep = 1, total_timesteps
      call  tempo_driver(tempo_cfgs=tempo_cfgs, itimestep=itimestep, dt=dt, &
                        ids=ids, ide=ide, ims=ims, ime=ime, its=its, ite=ite, &
                        jds=jds, jde=jde, jms=jms, jme=jme, jts=jts, jte=jte, &
                        kds=kds, kde=kde, kms=kms, kme=kme, kts=kts, kte=kte, &
                        t=t, p=p, w=w, dz=dz, qv=qv, th=th, pii=pii, &
                        qc=qc, qr=qr, qi=qi, qs=qs, qg=qg,  ni=ni, nr=nr, ng=ng, qb=qb, &
                        tempo_diags=tempo_driver_diags)
      precip_sum = precip_sum + tempo_driver_diags%snow_liquid_equiv_precip(1,1)
      write(io2,'(I7, 1E12.4)') int(itimestep*dt), precip_sum

      if (dt == 1. .and. itimestep == 1) then
        open(newunit=io3, file="snow_sedi_init.txt", status="new", action="write")
        do k = 1, nz
          write(io3,'(I5, 4E12.4)') k, qs(1,k,1)
        enddo 
      endif 
    enddo

    write(tt_string, '(I7)') integration_time
    open(newunit=io4, file="snow_sedi_dt_"//trim(adjustl(dt_string))// &
      "_runtime"//trim(adjustl(tt_string))//".txt", status="new", action="write")
    do k = 1, nz
      write(io4,'(I5, 4E12.4)') k, qs(1,k,1)
    enddo 
  end subroutine test_snow_sedimentation

end module tests