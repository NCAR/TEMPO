module tests
  !! TEMPO tests
  use module_mp_tempo_cfgs, only : ty_tempo_cfgs
  use module_mp_tempo_driver, only : tempo_init, tempo_run, ty_tempo_driver_diags, &
    tempo_run_enter_data, tempo_run_exit_data, tempo_run_sync_host_diags, tempo_run_sync_host_fields
  use module_mp_tempo_main, only : tempo_tend_init, tempo_tend_finalize
  use module_mp_nvtx, only : nvtx_range_push, nvtx_range_pop
  implicit none
  private

  public :: test_tempo_init, test_graupel_sedimentation, test_snow_sedimentation, &
    test_cloud_number_aerosolaware, test_cloud_number_non_aerosolaware, &
    test_cloud_number_ml, test_ml_cloud_effective_radius, set_ncells, set_nout_values, set_stride

  integer :: nCells = 1  !! number of cells (set via -c command line); ide, ime, ite, jde, jme, jte = nCells
  integer :: nOutValues = 3  !! number of (i,j) pairs written in graupel sedimentation profiles (see set_nout_values)
  integer :: nStride = -1  !! hybrid horizontal block size (set via -s command line); <=0 => full plane (=nCells)

  type(ty_tempo_cfgs) :: tempo_cfgs
  contains

  subroutine set_ncells(n)
    integer, intent(in) :: n
    nCells = max(1, n)
  end subroutine set_ncells

  subroutine set_nout_values(n)
    integer, intent(in) :: n
    nOutValues = max(1, n)
  end subroutine set_nout_values

  subroutine set_stride(n)
    integer, intent(in) :: n
    nStride = n
  end subroutine set_stride

  !! Resolve the effective hybrid stride for a tile of n columns:
  !! unset / <=0 / > n  => full plane (n); otherwise the requested block size.
  integer function effective_stride(n) result(s)
    integer, intent(in) :: n
    if (nStride <= 0 .or. nStride > n) then
      s = n
    else
      s = nStride
    endif
  end function effective_stride

  subroutine test_tempo_init()
    !! test tempo initialization procedure
    !! use ml_for_bl_nc_flag = .true. to initialize ml data
    !! which is used for a few tests
    !! test specific flags are then set in each test
    call nvtx_range_push('test_tempo_init')
    call tempo_init(ml_for_bl_nc_flag = .true., tempo_cfgs=tempo_cfgs)
    call nvtx_range_pop()
  end subroutine test_tempo_init


 subroutine test_graupel_sedimentation(dt, semi_sedi)
    !! test graupel sedimentation
    use module_mp_tempo_params, only : wp, sp, dp

    integer :: itimestep
    real(wp), intent(in) :: dt
    logical, intent(in) :: semi_sedi
    integer, parameter :: nz = 59
    integer :: ids, ide, ims, ime, its, ite, jds, jde, jms, jme, jts, jte
    integer, parameter :: kds = 1, kde = nz, kms = 1, kme = nz, kts = 1, kte = nz
    real(wp), dimension(:,:,:), allocatable :: qv, t, th, pii, p, w, dz, &
      qc, qr, qi, qs, qg, ni, nr, ng, qb
    real(wp), dimension(nz) :: klevs_in, qv_in, qc_in, qr_in, qi_in, qs_in, &
      qg_in, ni_in, nr_in, nc_in, nwfa_in, nifa_in, theta_in, ng_in, volg_in, &
      pressure_in, w_in, dz_in
    real(wp) :: precip_sum
    real(wp) :: elapsed_time
    character(len=20) :: dt_string, tt_string
    character(len=20) :: semi_sedi_string
    type(ty_tempo_driver_diags) :: tempo_driver_diags
    type(ty_tempo_cfgs) :: tempo_cfgs
    integer :: io1, io2, io3, io4, k, i, j, total_timesteps
    integer(kind=8) :: count_start, count_end, count_rate
    integer :: ip, ij_samp, n_samp
    integer, parameter :: integration_time = 1200

    ids = 1; ide = nCells; ims = 1; ime = nCells; its = 1; ite = nCells
    jds = 1; jde = nCells; jms = 1; jme = nCells; jts = 1; jte = nCells
    allocate(qv(kts:kte, its:ite, jts:jte), t(kts:kte, its:ite, jts:jte), th(kts:kte, its:ite, jts:jte), &
      pii(kts:kte, its:ite, jts:jte), p(kts:kte, its:ite, jts:jte), w(kts:kte, its:ite, jts:jte), &
      dz(kts:kte, its:ite, jts:jte), qc(kts:kte, its:ite, jts:jte), qr(kts:kte, its:ite, jts:jte), &
      qi(kts:kte, its:ite, jts:jte), qs(kts:kte, its:ite, jts:jte), qg(kts:kte, its:ite, jts:jte), &
      ni(kts:kte, its:ite, jts:jte), nr(kts:kte, its:ite, jts:jte), ng(kts:kte, its:ite, jts:jte), &
      qb(kts:kte, its:ite, jts:jte))

    call nvtx_range_push('test_graupel_sedimentation')

    write(dt_string, '(I7)') int(dt)
    
    if (semi_sedi) then
      write(semi_sedi_string, '(A)') '_semi_sedi'
    else
      write(semi_sedi_string, '(A)') ''
    endif

    ! read input file
    call nvtx_range_push('test_graupel_sedimentation_read_input')
    io1 = 11
    open (io1, file='test/data/mpas_59lev_test.txt', status='old')
    read(io1,*) ! header
    do k = 1, nz
      read(io1,*) klevs_in(k), qv_in(k), qc_in(k), qr_in(k), qi_in(k), qs_in(k), &
          qg_in(k), ni_in(k), nr_in(k), nc_in(k), nwfa_in(k), nifa_in(k), &
          theta_in(k), ng_in(k), volg_in(k), pressure_in(k), w_in(k), dz_in(k)
    end do
    close(io1)
    call nvtx_range_pop()

    call nvtx_range_push('test_graupel_sedimentation_fill_tile')
    do j = jts, jte
      do i = its, ite
        qv(:,i,j) = qv_in
        t(:,i,j) = theta_in * (pressure_in/100000.)**0.286
        p(:,i,j) = pressure_in
        dz(:,i,j) = dz_in
        qg(:,i,j) = qg_in
        ng(:,i,j) = ng_in
        qb(:,i,j) = volg_in * 1000. ! convert meters^3 -> liters
      end do
    end do
    w = 0._wp
    qc = 0._wp
    qr = 0._wp
    nr = 0._wp
    qi = 0._wp
    qs = 0._wp
    ni = 0._wp
    call nvtx_range_pop()

    ! set configs
    tempo_cfgs%turn_off_micro_flag = .true.
    tempo_cfgs%graupel_med_vol_diam_flag = .true.
    tempo_cfgs%semi_sedi_flag = semi_sedi

    total_timesteps = int(integration_time/dt)
    precip_sum = 0._wp
    open(newunit=io2, file="graupel_precip_dt_"//trim(adjustl(dt_string))//""//trim(semi_sedi_string)//".txt", &
      status="new", action="write")

    call tempo_tend_init(kts, kte, its, ite, jts, jte)

    call nvtx_range_push('tempo_run_enter_data')
    call tempo_run_enter_data(tempo_cfgs=tempo_cfgs, tempo_diags=tempo_driver_diags, use_temperature=.true., &
      t=t, p=p, w=w, dz=dz, qv=qv, th=th, pii=pii, qc=qc, qr=qr, qi=qi, qs=qs, qg=qg, ni=ni, nr=nr, &
      ng=ng, qb=qb, its=its, ite=ite, jts=jts, jte=jte, kts=kts, kte=kte)
    call nvtx_range_pop()

    call nvtx_range_push('test_graupel_sedimentation_timestep_integrate')
    elapsed_time = 0._wp
    do itimestep = 1, total_timesteps
      call nvtx_range_push('test_graupel_sedimentation_timestep')
      call system_clock(count_start, count_rate)
      call nvtx_range_push('test_graupel_sedimentation_tempo_run')
      call  tempo_run(tempo_cfgs=tempo_cfgs, itimestep=itimestep, dt=dt, &
                      ids=ids, ide=ide, ims=ims, ime=ime, its=its, ite=ite, &
                      jds=jds, jde=jde, jms=jms, jme=jme, jts=jts, jte=jte, &
                      kds=kds, kde=kde, kms=kms, kme=kme, kts=kts, kte=kte, &
                      t=t, p=p, w=w, dz=dz, qv=qv, th=th, pii=pii, &
                      qc=qc, qr=qr, qi=qi, qs=qs, qg=qg,  ni=ni, nr=nr, ng=ng, qb=qb, &
                      tempo_diags=tempo_driver_diags, stride=effective_stride(nCells), arguments_on_device=.true.)
      call nvtx_range_pop()
      call system_clock(count_end)
      if (count_rate > 0) elapsed_time = elapsed_time + real(count_end - count_start, kind=wp) / real(count_rate, wp)
      call tempo_run_sync_host_diags(tempo_driver_diags)
      precip_sum = precip_sum + tempo_driver_diags%graupel_liquid_equiv_precip(1,1)
      write(io2,'(I7, 1E12.4)') int(itimestep*dt), precip_sum

      if (dt == 1. .and. .not. semi_sedi .and. itimestep == 1) then
        call tempo_run_sync_host_fields(use_temperature=.true., t=t, qv=qv, qc=qc, qr=qr, qi=qi, qs=qs, qg=qg, &
          ni=ni, nr=nr, ng=ng, qb=qb)
        open(newunit=io3, file="graupel_sedi_init.txt", status="new", action="write")
        write(io3, '(5A)') 'k ', 'mass ', 'number ', 'density ' , 'mvd'
        n_samp = max(1, nOutValues)
        do ip = 1, n_samp
          if (n_samp == 1) then
            ij_samp = 1
          else
            ij_samp = 1 + nint(real((ip - 1) * (nCells - 1), wp) / real(n_samp - 1, wp))
          endif
          ij_samp = max(1, min(nCells, ij_samp))
          i = ij_samp
          j = ij_samp
          write(io3,'(A,2I8)') '# i j ', i, j
          do k = 1, nz
            if (qb(k,i,j) > 0._wp) then
              write(io3,'(I5, 4E12.4)') k, qg(k,i,j), ng(k,i,j), 1000._wp*qg(k,i,j)/qb(k,i,j), &
                tempo_driver_diags%graupel_med_vol_diam(k,i,j)
            else
              write(io3,'(I5, 4E12.4)') k, qg(k,i,j), ng(k,i,j), 0._wp, &
                tempo_driver_diags%graupel_med_vol_diam(k,i,j)
            endif
          enddo
        enddo
      endif 
      call nvtx_range_pop()
    enddo

    call nvtx_range_pop()

    call nvtx_range_push('tempo_run_exit_data')
    call tempo_run_exit_data(tempo_cfgs=tempo_cfgs, tempo_diags=tempo_driver_diags, use_temperature=.true., &
      t=t, p=p, w=w, dz=dz, qv=qv, th=th, pii=pii, qc=qc, qr=qr, qi=qi, qs=qs, qg=qg, ni=ni, nr=nr, &
      ng=ng, qb=qb, kts=kts, kte=kte, its=its, ite=ite, jts=jts, jte=jte)
    call nvtx_range_pop()

    call tempo_tend_finalize()

    write(io2,'(A,F12.4,A)') 'test_graupel_sedimentation compute time: ', elapsed_time, ' s'
    write(*,'(A,F12.4,A)') 'test_graupel_sedimentation compute time: ', elapsed_time, ' s'

    write(tt_string, '(I7)') integration_time
    open(newunit=io4, file="graupel_sedi_dt_"//trim(adjustl(dt_string))//""//trim(semi_sedi_string)//&
      "_runtime"//trim(adjustl(tt_string))//".txt", status="new", action="write")
    write(io4, '(5A)') 'k ', 'mass ', 'number ', 'density ' , 'mvd'
    n_samp = max(1, nOutValues)
    do ip = 1, n_samp
      if (n_samp == 1) then
        ij_samp = 1
      else
        ij_samp = 1 + nint(real((ip - 1) * (nCells - 1), wp) / real(n_samp - 1, wp))
      endif
      ij_samp = max(1, min(nCells, ij_samp))
      i = ij_samp
      j = ij_samp
      write(io4,'(A,2I8)') '# i j ', i, j
      do k = 1, nz
        if (qb(k,i,j) > 0._wp) then
          write(io4,'(I5, 4E12.4)') k, qg(k,i,j), ng(k,i,j), 1000._wp*qg(k,i,j)/qb(k,i,j), &
            tempo_driver_diags%graupel_med_vol_diam(k,i,j)
        else
          write(io4,'(I5, 4E12.4)') k, qg(k,i,j), ng(k,i,j), 0._wp, &
            tempo_driver_diags%graupel_med_vol_diam(k,i,j)
        endif
      enddo
    enddo
    call nvtx_range_pop()
    deallocate(qv, t, th, pii, p, w, dz, qc, qr, qi, qs, qg, ni, nr, ng, qb)
  end subroutine test_graupel_sedimentation


  subroutine test_snow_sedimentation(dt)
    use module_mp_tempo_params, only : wp, sp, dp

    integer :: itimestep
    real(wp), intent(in) :: dt
    integer, parameter :: nz = 59
    integer :: ids, ide, ims, ime, its, ite, jds, jde, jms, jme, jts, jte
    integer, parameter :: kds = 1, kde = nz, kms = 1, kme = nz, kts = 1, kte = nz
    real(wp), dimension(:,:,:), allocatable :: qv, t, th, pii, p, w, dz, &
      qc, qr, qi, qs, qg, ni, nr, ng, qb
    real(wp), dimension(nz) :: klevs_in, qv_in, qc_in, qr_in, qi_in, qs_in, &
      qg_in, ni_in, nr_in, nc_in, nwfa_in, nifa_in, theta_in, ng_in, volg_in, &
      pressure_in, w_in, dz_in
    real(wp) :: precip_sum
    real(wp) :: elapsed_time
    character(len=20) :: dt_string, tt_string
    type(ty_tempo_driver_diags) :: tempo_driver_diags
    type(ty_tempo_cfgs) :: tempo_cfgs
    integer :: io1, io2, io3, io4, k, i, j, total_timesteps
    integer(kind=8) :: count_start, count_end, count_rate
    integer :: ip, ij_samp, n_samp
    integer, parameter :: integration_time = 1200

    ids = 1; ide = nCells; ims = 1; ime = nCells; its = 1; ite = nCells
    jds = 1; jde = nCells; jms = 1; jme = nCells; jts = 1; jte = nCells
    allocate(qv(kts:kte, its:ite, jts:jte), t(kts:kte, its:ite, jts:jte), th(kts:kte, its:ite, jts:jte), &
      pii(kts:kte, its:ite, jts:jte), p(kts:kte, its:ite, jts:jte), w(kts:kte, its:ite, jts:jte), &
      dz(kts:kte, its:ite, jts:jte), qc(kts:kte, its:ite, jts:jte), qr(kts:kte, its:ite, jts:jte), &
      qi(kts:kte, its:ite, jts:jte), qs(kts:kte, its:ite, jts:jte), qg(kts:kte, its:ite, jts:jte), &
      ni(kts:kte, its:ite, jts:jte), nr(kts:kte, its:ite, jts:jte), ng(kts:kte, its:ite, jts:jte), &
      qb(kts:kte, its:ite, jts:jte))

    call nvtx_range_push('test_snow_sedimentation')

    write(dt_string, '(I7)') int(dt)

    ! read input file
    call nvtx_range_push('test_snow_sedimentation_read_input')
    io1 = 11
    open (io1, file='test/data/mpas_59lev_test.txt', status='old')
    read(io1,*) ! header
    do k = 1, nz
      read(io1,*) klevs_in(k), qv_in(k), qc_in(k), qr_in(k), qi_in(k), qs_in(k), &
          qg_in(k), ni_in(k), nr_in(k), nc_in(k), nwfa_in(k), nifa_in(k), &
          theta_in(k), ng_in(k), volg_in(k), pressure_in(k), w_in(k), dz_in(k)
    end do
    close(io1)
    call nvtx_range_pop()

    call nvtx_range_push('test_snow_sedimentation_fill_tile')
    do j = jts, jte
      do i = its, ite
        qv(:,i,j) = qv_in
        t(:,i,j) = theta_in * (pressure_in/100000.)**0.286
        p(:,i,j) = pressure_in
        dz(:,i,j) = dz_in
        qs(:,i,j) = qs_in
      end do
    end do
    w = 0._wp
    qc = 0._wp
    qr = 0._wp
    nr = 0._wp
    qi = 0._wp
    ni = 0._wp
    qg = 0._wp
    ng = 0._wp
    qb = 0._wp
    call nvtx_range_pop()

    ! set configs
    tempo_cfgs%turn_off_micro_flag= .true.

    total_timesteps = int(integration_time/dt)
    precip_sum = 0._wp
    open(newunit=io2, file="snow_precip_dt_"//trim(adjustl(dt_string))//".txt", &
      status="new", action="write")

    call tempo_tend_init(kts, kte, its, ite, jts, jte)

    call nvtx_range_push('tempo_run_enter_data')
    call tempo_run_enter_data(tempo_cfgs=tempo_cfgs, tempo_diags=tempo_driver_diags, use_temperature=.true., &
      t=t, p=p, w=w, dz=dz, qv=qv, th=th, pii=pii, qc=qc, qr=qr, qi=qi, qs=qs, qg=qg, ni=ni, nr=nr, &
      ng=ng, qb=qb, its=its, ite=ite, jts=jts, jte=jte, kts=kts, kte=kte)
    call nvtx_range_pop()

    call nvtx_range_push('test_snow_sedimentation_timestep_integrate')
    elapsed_time = 0._wp
    do itimestep = 1, total_timesteps
      call nvtx_range_push('test_snow_sedimentation_timestep')
      call system_clock(count_start, count_rate)
      call nvtx_range_push('test_snow_sedimentation_tempo_run')
      call  tempo_run(tempo_cfgs=tempo_cfgs, itimestep=itimestep, dt=dt, &
                      ids=ids, ide=ide, ims=ims, ime=ime, its=its, ite=ite, &
                      jds=jds, jde=jde, jms=jms, jme=jme, jts=jts, jte=jte, &
                      kds=kds, kde=kde, kms=kms, kme=kme, kts=kts, kte=kte, &
                      t=t, p=p, w=w, dz=dz, qv=qv, th=th, pii=pii, &
                      qc=qc, qr=qr, qi=qi, qs=qs, qg=qg,  ni=ni, nr=nr, ng=ng, qb=qb, &
                      tempo_diags=tempo_driver_diags, stride=effective_stride(nCells), arguments_on_device=.true.)
      call nvtx_range_pop()
      call system_clock(count_end)
      if (count_rate > 0) elapsed_time = elapsed_time + real(count_end - count_start, kind=wp) / real(count_rate, wp)
      call tempo_run_sync_host_diags(tempo_driver_diags)
      precip_sum = precip_sum + tempo_driver_diags%snow_liquid_equiv_precip(1,1)
      write(io2,'(I7, 1E12.4)') int(itimestep*dt), precip_sum

      if (dt == 1. .and. itimestep == 1) then
        call tempo_run_sync_host_fields(use_temperature=.true., t=t, qv=qv, qc=qc, qr=qr, qi=qi, qs=qs, qg=qg, &
          ni=ni, nr=nr, ng=ng, qb=qb)
        open(newunit=io3, file="snow_sedi_init.txt", status="new", action="write")
        n_samp = max(1, nOutValues)
        do ip = 1, n_samp
          if (n_samp == 1) then
            ij_samp = 1
          else
            ij_samp = 1 + nint(real((ip - 1) * (nCells - 1), wp) / real(n_samp - 1, wp))
          endif
          ij_samp = max(1, min(nCells, ij_samp))
          i = ij_samp
          j = ij_samp
          write(io3,'(A,2I8)') '# i j ', i, j
          do k = 1, nz
            write(io3,'(I5, 4E12.4)') k, qs(k,i,j)
          enddo
        enddo
      endif 
      call nvtx_range_pop()
    enddo

    call nvtx_range_pop()

    call nvtx_range_push('tempo_run_exit_data')
    call tempo_run_exit_data(tempo_cfgs=tempo_cfgs, tempo_diags=tempo_driver_diags, use_temperature=.true., &
      t=t, p=p, w=w, dz=dz, qv=qv, th=th, pii=pii, qc=qc, qr=qr, qi=qi, qs=qs, qg=qg, ni=ni, nr=nr, &
      ng=ng, qb=qb, kts=kts, kte=kte, its=its, ite=ite, jts=jts, jte=jte)
    call nvtx_range_pop()

    call tempo_tend_finalize()

    write(io2,'(A,F12.4,A)') 'test_snow_sedimentation compute time: ', elapsed_time, ' s'
    write(*,'(A,F12.4,A)') 'test_snow_sedimentation compute time: ', elapsed_time, ' s'

    write(tt_string, '(I7)') integration_time
    open(newunit=io4, file="snow_sedi_dt_"//trim(adjustl(dt_string))// &
      "_runtime"//trim(adjustl(tt_string))//".txt", status="new", action="write")
    n_samp = max(1, nOutValues)
    do ip = 1, n_samp
      if (n_samp == 1) then
        ij_samp = 1
      else
        ij_samp = 1 + nint(real((ip - 1) * (nCells - 1), wp) / real(n_samp - 1, wp))
      endif
      ij_samp = max(1, min(nCells, ij_samp))
      i = ij_samp
      j = ij_samp
      write(io4,'(A,2I8)') '# i j ', i, j
      do k = 1, nz
        write(io4,'(I5, 4E12.4)') k, qs(k,i,j)
      enddo
    enddo
    call nvtx_range_pop()
    deallocate(qv, t, th, pii, p, w, dz, qc, qr, qi, qs, qg, ni, nr, ng, qb)
  end subroutine test_snow_sedimentation


  subroutine test_cloud_number_aerosolaware(dt)
    use module_mp_tempo_params, only : wp, sp, dp

    integer :: itimestep
    real(wp), intent(in) :: dt
    integer, parameter :: nz = 59
    integer :: ids, ide, ims, ime, its, ite, jds, jde, jms, jme, jts, jte
    integer, parameter :: kds = 1, kde = nz, kms = 1, kme = nz, kts = 1, kte = nz
    real(wp), dimension(:,:,:), allocatable :: qv, t, th, pii, p, w, dz, &
      qc, qr, qi, qs, qg, ni, nr, ng, qb, nc
    real(wp), dimension(nz) :: klevs_in, qv_in, qc_in, qr_in, qi_in, qs_in, &
      qg_in, ni_in, nr_in, nc_in, nwfa_in, nifa_in, theta_in, ng_in, volg_in, &
      pressure_in, w_in, dz_in
    real(wp) :: precip_sum
    real(wp) :: elapsed_time
    character(len=20) :: dt_string, tt_string
    type(ty_tempo_driver_diags) :: tempo_driver_diags
    type(ty_tempo_cfgs) :: tempo_cfgs
    integer :: io1, io2, io3, io4, k, i, j, total_timesteps
    integer(kind=8) :: count_start, count_end, count_rate
    integer :: ip, ij_samp, n_samp
    integer, parameter :: integration_time = 1200

    ids = 1; ide = nCells; ims = 1; ime = nCells; its = 1; ite = nCells
    jds = 1; jde = nCells; jms = 1; jme = nCells; jts = 1; jte = nCells
    allocate(qv(kts:kte, its:ite, jts:jte), t(kts:kte, its:ite, jts:jte), th(kts:kte, its:ite, jts:jte), &
      pii(kts:kte, its:ite, jts:jte), p(kts:kte, its:ite, jts:jte), w(kts:kte, its:ite, jts:jte), &
      dz(kts:kte, its:ite, jts:jte), qc(kts:kte, its:ite, jts:jte), qr(kts:kte, its:ite, jts:jte), &
      qi(kts:kte, its:ite, jts:jte), qs(kts:kte, its:ite, jts:jte), qg(kts:kte, its:ite, jts:jte), &
      ni(kts:kte, its:ite, jts:jte), nr(kts:kte, its:ite, jts:jte), ng(kts:kte, its:ite, jts:jte), &
      qb(kts:kte, its:ite, jts:jte), nc(kts:kte, its:ite, jts:jte))

    call nvtx_range_push('test_cloud_number_aerosolaware')

    write(dt_string, '(I7)') int(dt)

    ! read input file
    call nvtx_range_push('test_cloud_number_aerosolaware_read_input')
    io1 = 11
    open (io1, file='test/data/mpas_59lev_test.txt', status='old')
    read(io1,*) ! header
    do k = 1, nz
      read(io1,*) klevs_in(k), qv_in(k), qc_in(k), qr_in(k), qi_in(k), qs_in(k), &
          qg_in(k), ni_in(k), nr_in(k), nc_in(k), nwfa_in(k), nifa_in(k), &
          theta_in(k), ng_in(k), volg_in(k), pressure_in(k), w_in(k), dz_in(k)
    end do
    close(io1)
    call nvtx_range_pop()

    call nvtx_range_push('test_cloud_number_aerosolaware_fill_tile')
    do j = jts, jte
      do i = its, ite
        qv(:,i,j) = qv_in
        t(:,i,j) = theta_in * (pressure_in/100000.)**0.286
        p(:,i,j) = pressure_in
        dz(:,i,j) = dz_in
        w(:,i,j) = w_in
        qc(:,i,j) = qc_in
        nc(:,i,j) = nc_in
        qr(:,i,j) = qr_in
        nr(:,i,j) = nr_in
        qi(:,i,j) =  qi_in
        qs(:,i,j) = qs_in
        ni(:,i,j) = ni_in
      end do
    end do
    qg = 0._wp
    ng = 0._wp
    qb = 0._wp
    call nvtx_range_pop()

    ! set configs
    tempo_cfgs%turn_off_micro_flag= .true.

    total_timesteps = int(integration_time/dt)

    call tempo_tend_init(kts, kte, its, ite, jts, jte)

    call nvtx_range_push('tempo_run_enter_data')
    call tempo_run_enter_data(tempo_cfgs=tempo_cfgs, tempo_diags=tempo_driver_diags, use_temperature=.true., &
      t=t, p=p, w=w, dz=dz, qv=qv, th=th, pii=pii, qc=qc, nc=nc, qr=qr, qi=qi, qs=qs, qg=qg, ni=ni, nr=nr, &
      ng=ng, qb=qb, its=its, ite=ite, jts=jts, jte=jte, kts=kts, kte=kte)
    call nvtx_range_pop()

    call nvtx_range_push('test_cloud_number_aerosolaware_timestep_integrate')
    elapsed_time = 0._wp
    do itimestep = 1, total_timesteps
      call nvtx_range_push('test_cloud_number_aerosolaware_timestep')
      call system_clock(count_start, count_rate)
      call nvtx_range_push('test_cloud_number_aerosolaware_tempo_run')
      call  tempo_run(tempo_cfgs=tempo_cfgs, itimestep=itimestep, dt=dt, &
                      ids=ids, ide=ide, ims=ims, ime=ime, its=its, ite=ite, &
                      jds=jds, jde=jde, jms=jms, jme=jme, jts=jts, jte=jte, &
                      kds=kds, kde=kde, kms=kms, kme=kme, kts=kts, kte=kte, &
                      t=t, p=p, w=w, dz=dz, qv=qv, th=th, pii=pii, &
                      qc=qc, nc=nc, qr=qr, qi=qi, qs=qs, qg=qg, ni=ni, nr=nr, &
                      tempo_diags=tempo_driver_diags, stride=effective_stride(nCells), arguments_on_device=.true.)
      call nvtx_range_pop()
      call system_clock(count_end)
      if (count_rate > 0) elapsed_time = elapsed_time + real(count_end - count_start, kind=wp) / real(count_rate, wp)

      if (itimestep == 1) then
        call tempo_run_sync_host_fields(use_temperature=.true., t=t, qv=qv, qc=qc, qr=qr, qi=qi, qs=qs, qg=qg, &
          ni=ni, nr=nr, ng=ng, qb=qb, nc=nc)
        open(newunit=io3, file="cloud_number_init.txt", status="new", action="write")
        n_samp = max(1, nOutValues)
        do ip = 1, n_samp
          if (n_samp == 1) then
            ij_samp = 1
          else
            ij_samp = 1 + nint(real((ip - 1) * (nCells - 1), wp) / real(n_samp - 1, wp))
          endif
          ij_samp = max(1, min(nCells, ij_samp))
          i = ij_samp
          j = ij_samp
          write(io3,'(A,2I8)') '# i j ', i, j
          do k = 1, nz
            write(io3,'(I5, 4E12.4)') k, qc(k,i,j), nc(k,i,j)
          enddo
        enddo
      endif 
      call nvtx_range_pop()
    enddo

    call nvtx_range_pop()

    call nvtx_range_push('tempo_run_exit_data')
    call tempo_run_exit_data(tempo_cfgs=tempo_cfgs, tempo_diags=tempo_driver_diags, use_temperature=.true., &
      t=t, p=p, w=w, dz=dz, qv=qv, th=th, pii=pii, qc=qc, nc=nc, qr=qr, qi=qi, qs=qs, qg=qg, ni=ni, nr=nr, &
      ng=ng, qb=qb, kts=kts, kte=kte, its=its, ite=ite, jts=jts, jte=jte)
    call nvtx_range_pop()

    call tempo_tend_finalize()

    write(tt_string, '(I7)') integration_time
    open(newunit=io4, file="cloud_number_dt_"//trim(adjustl(dt_string))// &
      "_runtime"//trim(adjustl(tt_string))//".txt", status="new", action="write")
    n_samp = max(1, nOutValues)
    do ip = 1, n_samp
      if (n_samp == 1) then
        ij_samp = 1
      else
        ij_samp = 1 + nint(real((ip - 1) * (nCells - 1), wp) / real(n_samp - 1, wp))
      endif
      ij_samp = max(1, min(nCells, ij_samp))
      i = ij_samp
      j = ij_samp
      write(io4,'(A,2I8)') '# i j ', i, j
      do k = 1, nz
        write(io4,'(I5, 4E12.4)') k, qc(k,i,j), nc(k,i,j)
      enddo
    enddo
    write(io4,'(A,F12.4,A)') 'test_cloud_number_aerosolaware compute time: ', elapsed_time, ' s'
    write(*,'(A,F12.4,A)') 'test_cloud_number_aerosolaware compute time: ', elapsed_time, ' s'

    call nvtx_range_pop()
    deallocate(qv, t, th, pii, p, w, dz, qc, qr, qi, qs, qg, ni, nr, ng, qb, nc)
  end subroutine test_cloud_number_aerosolaware


  subroutine test_cloud_number_non_aerosolaware(dt)
    use module_mp_tempo_params, only : wp, sp, dp

    integer :: itimestep
    real(wp), intent(in) :: dt
    integer, parameter :: nz = 59
    integer :: ids, ide, ims, ime, its, ite, jds, jde, jms, jme, jts, jte
    integer, parameter :: kds = 1, kde = nz, kms = 1, kme = nz, kts = 1, kte = nz
    real(wp), dimension(:,:,:), allocatable :: qv, t, th, pii, p, w, dz, &
      qc, qr, qi, qs, qg, ni, nr, ng, qb, nc
    real(wp), dimension(nz) :: klevs_in, qv_in, qc_in, qr_in, qi_in, qs_in, &
      qg_in, ni_in, nr_in, nc_in, nwfa_in, nifa_in, theta_in, ng_in, volg_in, &
      pressure_in, w_in, dz_in
    real(wp) :: precip_sum
    real(wp) :: elapsed_time
    character(len=20) :: dt_string, tt_string
    type(ty_tempo_driver_diags) :: tempo_driver_diags
    type(ty_tempo_cfgs) :: tempo_cfgs
    integer :: io1, io2, io3, io4, k, i, j, total_timesteps
    integer(kind=8) :: count_start, count_end, count_rate
    integer :: ip, ij_samp, n_samp
    integer, parameter :: integration_time = 1200

    ids = 1; ide = nCells; ims = 1; ime = nCells; its = 1; ite = nCells
    jds = 1; jde = nCells; jms = 1; jme = nCells; jts = 1; jte = nCells
    allocate(qv(kts:kte, its:ite, jts:jte), t(kts:kte, its:ite, jts:jte), th(kts:kte, its:ite, jts:jte), &
      pii(kts:kte, its:ite, jts:jte), p(kts:kte, its:ite, jts:jte), w(kts:kte, its:ite, jts:jte), &
      dz(kts:kte, its:ite, jts:jte), qc(kts:kte, its:ite, jts:jte), qr(kts:kte, its:ite, jts:jte), &
      qi(kts:kte, its:ite, jts:jte), qs(kts:kte, its:ite, jts:jte), qg(kts:kte, its:ite, jts:jte), &
      ni(kts:kte, its:ite, jts:jte), nr(kts:kte, its:ite, jts:jte), ng(kts:kte, its:ite, jts:jte), &
      qb(kts:kte, its:ite, jts:jte), nc(kts:kte, its:ite, jts:jte))

    call nvtx_range_push('test_cloud_number_non_aerosolaware')

    write(dt_string, '(I7)') int(dt)

    ! read input file
    call nvtx_range_push('test_cloud_number_non_aerosolaware_read_input')
    io1 = 11
    open (io1, file='test/data/mpas_59lev_test.txt', status='old')
    read(io1,*) ! header
    do k = 1, nz
      read(io1,*) klevs_in(k), qv_in(k), qc_in(k), qr_in(k), qi_in(k), qs_in(k), &
          qg_in(k), ni_in(k), nr_in(k), nc_in(k), nwfa_in(k), nifa_in(k), &
          theta_in(k), ng_in(k), volg_in(k), pressure_in(k), w_in(k), dz_in(k)
    end do
    close(io1)
    call nvtx_range_pop()

    call nvtx_range_push('test_cloud_number_non_aerosolaware_fill_tile')
    do j = jts, jte
      do i = its, ite
        qv(:,i,j) = qv_in
        t(:,i,j) = theta_in * (pressure_in/100000.)**0.286
        p(:,i,j) = pressure_in
        dz(:,i,j) = dz_in
        w(:,i,j) = w_in
        qc(:,i,j) = qc_in
        nc(:,i,j) = 0._wp
        qr(:,i,j) = qr_in
        nr(:,i,j) = nr_in
        qi(:,i,j) =  qi_in
        qs(:,i,j) = qs_in
        ni(:,i,j) = ni_in
      end do
    end do
    qg = 0._wp
    ng = 0._wp
    qb = 0._wp
    call nvtx_range_pop()

    ! set configs
    tempo_cfgs%turn_off_micro_flag= .true.
    tempo_cfgs%cloud_number_mixing_ratio_flag= .true.
    tempo_cfgs%ml_for_nc_flag = .false.

    total_timesteps = int(integration_time/dt)

    call tempo_tend_init(kts, kte, its, ite, jts, jte)

    call nvtx_range_push('tempo_run_enter_data')
    call tempo_run_enter_data(tempo_cfgs=tempo_cfgs, tempo_diags=tempo_driver_diags, use_temperature=.true., &
      t=t, p=p, w=w, dz=dz, qv=qv, th=th, pii=pii, qc=qc, qr=qr, qi=qi, qs=qs, qg=qg, ni=ni, nr=nr, &
      ng=ng, qb=qb, its=its, ite=ite, jts=jts, jte=jte, kts=kts, kte=kte)
    call nvtx_range_pop()

    call nvtx_range_push('test_cloud_number_non_aerosolaware_timestep_integrate')
    elapsed_time = 0._wp
    do itimestep = 1, total_timesteps
      call nvtx_range_push('test_cloud_number_non_aerosolaware_timestep')
      call system_clock(count_start, count_rate)
      call nvtx_range_push('test_cloud_number_non_aerosolaware_tempo_run')
      call  tempo_run(tempo_cfgs=tempo_cfgs, itimestep=itimestep, dt=dt, &
                      ids=ids, ide=ide, ims=ims, ime=ime, its=its, ite=ite, &
                      jds=jds, jde=jde, jms=jms, jme=jme, jts=jts, jte=jte, &
                      kds=kds, kde=kde, kms=kms, kme=kme, kts=kts, kte=kte, &
                      t=t, p=p, w=w, dz=dz, qv=qv, th=th, pii=pii, &
                      qc=qc, qr=qr, qi=qi, qs=qs, qg=qg, ni=ni, nr=nr, &
                      tempo_diags=tempo_driver_diags, stride=effective_stride(nCells), arguments_on_device=.true.)
      call nvtx_range_pop()
      call system_clock(count_end)
      if (count_rate > 0) elapsed_time = elapsed_time + real(count_end - count_start, kind=wp) / real(count_rate, wp)

      if (itimestep == 1) then
        call tempo_run_sync_host_diags(tempo_driver_diags)
        call tempo_run_sync_host_fields(use_temperature=.true., t=t, qv=qv, qc=qc, qr=qr, qi=qi, qs=qs, qg=qg, &
          ni=ni, nr=nr, ng=ng, qb=qb)
        open(newunit=io3, file="cloud_number_constant_init.txt", status="new", action="write")
        n_samp = max(1, nOutValues)
        do ip = 1, n_samp
          if (n_samp == 1) then
            ij_samp = 1
          else
            ij_samp = 1 + nint(real((ip - 1) * (nCells - 1), wp) / real(n_samp - 1, wp))
          endif
          ij_samp = max(1, min(nCells, ij_samp))
          i = ij_samp
          j = ij_samp
          write(io3,'(A,2I8)') '# i j ', i, j
          do k = 1, nz
            write(io3,'(I5, 4E12.4)') k, qc(k,i,j), tempo_driver_diags%cloud_number_mixing_ratio(k,i,j)
          enddo
        enddo
      endif 
      call nvtx_range_pop()
    enddo

    call nvtx_range_pop()

    call nvtx_range_push('tempo_run_exit_data')
    call tempo_run_exit_data(tempo_cfgs=tempo_cfgs, tempo_diags=tempo_driver_diags, use_temperature=.true., &
      t=t, p=p, w=w, dz=dz, qv=qv, th=th, pii=pii, qc=qc, qr=qr, qi=qi, qs=qs, qg=qg, ni=ni, nr=nr, &
      ng=ng, qb=qb, kts=kts, kte=kte, its=its, ite=ite, jts=jts, jte=jte)
    call nvtx_range_pop()

    call tempo_tend_finalize()

    write(tt_string, '(I7)') integration_time
    open(newunit=io4, file="cloud_number_constant_dt_"//trim(adjustl(dt_string))// &
      "_runtime"//trim(adjustl(tt_string))//".txt", status="new", action="write")
    n_samp = max(1, nOutValues)
    do ip = 1, n_samp
      if (n_samp == 1) then
        ij_samp = 1
      else
        ij_samp = 1 + nint(real((ip - 1) * (nCells - 1), wp) / real(n_samp - 1, wp))
      endif
      ij_samp = max(1, min(nCells, ij_samp))
      i = ij_samp
      j = ij_samp
      write(io4,'(A,2I8)') '# i j ', i, j
      do k = 1, nz
        write(io4,'(I5, 4E12.4)') k, qc(k,i,j), tempo_driver_diags%cloud_number_mixing_ratio(k,i,j)
      enddo
    enddo
    write(io4,'(A,F12.4,A)') 'test_cloud_number_non_aerosolaware compute time: ', elapsed_time, ' s'
    write(*,'(A,F12.4,A)') 'test_cloud_number_non_aerosolaware compute time: ', elapsed_time, ' s'

    call nvtx_range_pop()
    deallocate(qv, t, th, pii, p, w, dz, qc, qr, qi, qs, qg, ni, nr, ng, qb, nc)
  end subroutine test_cloud_number_non_aerosolaware


  subroutine test_cloud_number_ml(dt)
    use module_mp_tempo_params, only : wp, sp, dp

    integer :: itimestep
    real(wp), intent(in) :: dt
    integer, parameter :: nz = 59
    integer :: ids, ide, ims, ime, its, ite, jds, jde, jms, jme, jts, jte
    integer, parameter :: kds = 1, kde = nz, kms = 1, kme = nz, kts = 1, kte = nz
    real(wp), dimension(:,:,:), allocatable :: qv, t, th, pii, p, w, dz, &
      qc, qr, qi, qs, qg, ni, nr, ng, qb, nc
    real(wp), dimension(nz) :: klevs_in, qv_in, qc_in, qr_in, qi_in, qs_in, &
      qg_in, ni_in, nr_in, nc_in, nwfa_in, nifa_in, theta_in, ng_in, volg_in, &
      pressure_in, w_in, dz_in
    real(wp) :: precip_sum
    real(wp) :: elapsed_time
    character(len=20) :: dt_string, tt_string
    type(ty_tempo_driver_diags) :: tempo_driver_diags
    type(ty_tempo_cfgs) :: tempo_cfgs
    integer :: io1, io2, io3, io4, k, i, j, total_timesteps
    integer(kind=8) :: count_start, count_end, count_rate
    integer :: ip, ij_samp, n_samp
    integer, parameter :: integration_time = 1200

    ids = 1; ide = nCells; ims = 1; ime = nCells; its = 1; ite = nCells
    jds = 1; jde = nCells; jms = 1; jme = nCells; jts = 1; jte = nCells
    allocate(qv(kts:kte, its:ite, jts:jte), t(kts:kte, its:ite, jts:jte), th(kts:kte, its:ite, jts:jte), &
      pii(kts:kte, its:ite, jts:jte), p(kts:kte, its:ite, jts:jte), w(kts:kte, its:ite, jts:jte), &
      dz(kts:kte, its:ite, jts:jte), qc(kts:kte, its:ite, jts:jte), qr(kts:kte, its:ite, jts:jte), &
      qi(kts:kte, its:ite, jts:jte), qs(kts:kte, its:ite, jts:jte), qg(kts:kte, its:ite, jts:jte), &
      ni(kts:kte, its:ite, jts:jte), nr(kts:kte, its:ite, jts:jte), ng(kts:kte, its:ite, jts:jte), &
      qb(kts:kte, its:ite, jts:jte), nc(kts:kte, its:ite, jts:jte))

    call nvtx_range_push('test_cloud_number_ml')

    write(dt_string, '(I7)') int(dt)

    ! read input file
    call nvtx_range_push('test_cloud_number_ml_read_input')
    io1 = 11
    open (io1, file='test/data/mpas_59lev_test.txt', status='old')
    read(io1,*) ! header
    do k = 1, nz
      read(io1,*) klevs_in(k), qv_in(k), qc_in(k), qr_in(k), qi_in(k), qs_in(k), &
          qg_in(k), ni_in(k), nr_in(k), nc_in(k), nwfa_in(k), nifa_in(k), &
          theta_in(k), ng_in(k), volg_in(k), pressure_in(k), w_in(k), dz_in(k)
    end do
    close(io1)
    call nvtx_range_pop()

    call nvtx_range_push('test_cloud_number_ml_fill_tile')
    do j = jts, jte
      do i = its, ite
        qv(:,i,j) = qv_in
        t(:,i,j) = theta_in * (pressure_in/100000.)**0.286
        p(:,i,j) = pressure_in
        dz(:,i,j) = dz_in
        w(:,i,j) = w_in
        qc(:,i,j) = qc_in
        nc(:,i,j) = 0._wp
        qr(:,i,j) = qr_in
        nr(:,i,j) = nr_in
        qi(:,i,j) =  qi_in
        qs(:,i,j) = qs_in
        ni(:,i,j) = ni_in
      end do
    end do
    qg = 0._wp
    ng = 0._wp
    qb = 0._wp
    call nvtx_range_pop()

    ! set configs
    tempo_cfgs%turn_off_micro_flag = .true.
    tempo_cfgs%cloud_number_mixing_ratio_flag = .true.
    tempo_cfgs%ml_for_nc_flag = .true.

    total_timesteps = int(integration_time/dt)

    call tempo_tend_init(kts, kte, its, ite, jts, jte)

    call nvtx_range_push('tempo_run_enter_data')
    call tempo_run_enter_data(tempo_cfgs=tempo_cfgs, tempo_diags=tempo_driver_diags, use_temperature=.true., &
      t=t, p=p, w=w, dz=dz, qv=qv, th=th, pii=pii, qc=qc, qr=qr, qi=qi, qs=qs, qg=qg, ni=ni, nr=nr, &
      ng=ng, qb=qb, its=its, ite=ite, jts=jts, jte=jte, kts=kts, kte=kte)
    call nvtx_range_pop()

    call nvtx_range_push('test_cloud_number_ml_timestep_integrate')
    elapsed_time = 0._wp
    do itimestep = 1, total_timesteps
      call nvtx_range_push('test_cloud_number_ml_timestep')
      call system_clock(count_start, count_rate)
      call nvtx_range_push('test_cloud_number_ml_tempo_run')
      call  tempo_run(tempo_cfgs=tempo_cfgs, itimestep=itimestep, dt=dt, &
                      ids=ids, ide=ide, ims=ims, ime=ime, its=its, ite=ite, &
                      jds=jds, jde=jde, jms=jms, jme=jme, jts=jts, jte=jte, &
                      kds=kds, kde=kde, kms=kms, kme=kme, kts=kts, kte=kte, &
                      t=t, p=p, w=w, dz=dz, qv=qv, th=th, pii=pii, &
                      qc=qc, qr=qr, qi=qi, qs=qs, qg=qg, ni=ni, nr=nr, &
                      tempo_diags=tempo_driver_diags, stride=effective_stride(nCells), arguments_on_device=.true.)
      call nvtx_range_pop()
      call system_clock(count_end)
      if (count_rate > 0) elapsed_time = elapsed_time + real(count_end - count_start, kind=wp) / real(count_rate, wp)

      if (itimestep == 1) then
        call tempo_run_sync_host_diags(tempo_driver_diags)
        call tempo_run_sync_host_fields(use_temperature=.true., t=t, qv=qv, qc=qc, qr=qr, qi=qi, qs=qs, qg=qg, &
          ni=ni, nr=nr, ng=ng, qb=qb)
        open(newunit=io3, file="cloud_number_ml_init.txt", status="new", action="write")
        n_samp = max(1, nOutValues)
        do ip = 1, n_samp
          if (n_samp == 1) then
            ij_samp = 1
          else
            ij_samp = 1 + nint(real((ip - 1) * (nCells - 1), wp) / real(n_samp - 1, wp))
          endif
          ij_samp = max(1, min(nCells, ij_samp))
          i = ij_samp
          j = ij_samp
          write(io3,'(A,2I8)') '# i j ', i, j
          do k = 1, nz
            write(io3,'(I5, 4E12.4)') k, qc(k,i,j), tempo_driver_diags%cloud_number_mixing_ratio(k,i,j)
          enddo
        enddo
      endif 
      call nvtx_range_pop()
    enddo

    call nvtx_range_pop()

    call nvtx_range_push('tempo_run_exit_data')
    call tempo_run_exit_data(tempo_cfgs=tempo_cfgs, tempo_diags=tempo_driver_diags, use_temperature=.true., &
      t=t, p=p, w=w, dz=dz, qv=qv, th=th, pii=pii, qc=qc, qr=qr, qi=qi, qs=qs, qg=qg, ni=ni, nr=nr, &
      ng=ng, qb=qb, kts=kts, kte=kte, its=its, ite=ite, jts=jts, jte=jte)
    call nvtx_range_pop()

    call tempo_tend_finalize()

    write(tt_string, '(I7)') integration_time
    open(newunit=io4, file="cloud_number_ml_dt_"//trim(adjustl(dt_string))// &
      "_runtime"//trim(adjustl(tt_string))//".txt", status="new", action="write")
    n_samp = max(1, nOutValues)
    do ip = 1, n_samp
      if (n_samp == 1) then
        ij_samp = 1
      else
        ij_samp = 1 + nint(real((ip - 1) * (nCells - 1), wp) / real(n_samp - 1, wp))
      endif
      ij_samp = max(1, min(nCells, ij_samp))
      i = ij_samp
      j = ij_samp
      write(io4,'(A,2I8)') '# i j ', i, j
      do k = 1, nz
        write(io4,'(I5, 4E12.4)') k, qc(k,i,j), tempo_driver_diags%cloud_number_mixing_ratio(k,i,j)
      enddo
    enddo
    write(io4,'(A,F12.4,A)') 'test_cloud_number_ml compute time: ', elapsed_time, ' s'
    write(*,'(A,F12.4,A)') 'test_cloud_number_ml compute time: ', elapsed_time, ' s'

    call nvtx_range_pop()
    deallocate(qv, t, th, pii, p, w, dz, qc, qr, qi, qs, qg, ni, nr, ng, qb, nc)
  end subroutine test_cloud_number_ml


  subroutine test_ml_cloud_effective_radius(dt)
    use module_mp_tempo_params, only : wp, sp, dp

    integer :: itimestep
    real(wp), intent(in) :: dt
    integer, parameter :: nz = 59
    integer :: ids, ide, ims, ime, its, ite, jds, jde, jms, jme, jts, jte
    integer, parameter :: kds = 1, kde = nz, kms = 1, kme = nz, kts = 1, kte = nz
    real(wp), dimension(:,:,:), allocatable :: qv, t, th, pii, p, w, dz, &
      qc, qr, qi, qs, qg, ni, nr, ng, qb, nc, qc_bl, qcfrac_bl
    real(wp), dimension(nz) :: klevs_in, qv_in, qc_in, qr_in, qi_in, qs_in, &
      qg_in, ni_in, nr_in, nc_in, nwfa_in, nifa_in, theta_in, ng_in, volg_in, &
      pressure_in, w_in, dz_in
    real(wp) :: precip_sum
    real(wp) :: elapsed_time
    character(len=20) :: dt_string, tt_string
    type(ty_tempo_driver_diags) :: tempo_driver_diags
    type(ty_tempo_cfgs) :: tempo_cfgs
    integer :: io1, io2, io3, io4, k, i, j, total_timesteps
    integer(kind=8) :: count_start, count_end, count_rate
    integer :: ip, ij_samp, n_samp
    integer, parameter :: integration_time = 1200

    ids = 1; ide = nCells; ims = 1; ime = nCells; its = 1; ite = nCells
    jds = 1; jde = nCells; jms = 1; jme = nCells; jts = 1; jte = nCells
    allocate(qv(kts:kte, its:ite, jts:jte), t(kts:kte, its:ite, jts:jte), th(kts:kte, its:ite, jts:jte), &
      pii(kts:kte, its:ite, jts:jte), p(kts:kte, its:ite, jts:jte), w(kts:kte, its:ite, jts:jte), &
      dz(kts:kte, its:ite, jts:jte), qc(kts:kte, its:ite, jts:jte), qr(kts:kte, its:ite, jts:jte), &
      qi(kts:kte, its:ite, jts:jte), qs(kts:kte, its:ite, jts:jte), qg(kts:kte, its:ite, jts:jte), &
      ni(kts:kte, its:ite, jts:jte), nr(kts:kte, its:ite, jts:jte), ng(kts:kte, its:ite, jts:jte), &
      qb(kts:kte, its:ite, jts:jte), nc(kts:kte, its:ite, jts:jte), &
      qc_bl(kts:kte, its:ite, jts:jte), qcfrac_bl(kts:kte, its:ite, jts:jte))

    call nvtx_range_push('test_ml_cloud_effective_radius')

    write(dt_string, '(I7)') int(dt)

    ! read input file
    call nvtx_range_push('test_ml_cloud_effective_radius_read_input')
    io1 = 11
    open (io1, file='test/data/mpas_59lev_test.txt', status='old')
    read(io1,*) ! header
    do k = 1, nz
      read(io1,*) klevs_in(k), qv_in(k), qc_in(k), qr_in(k), qi_in(k), qs_in(k), &
          qg_in(k), ni_in(k), nr_in(k), nc_in(k), nwfa_in(k), nifa_in(k), &
          theta_in(k), ng_in(k), volg_in(k), pressure_in(k), w_in(k), dz_in(k)
    end do
    close(io1)
    call nvtx_range_pop()

    call nvtx_range_push('test_ml_cloud_effective_radius_fill_tile')
    do j = jts, jte
      do i = its, ite
        qv(:,i,j) = qv_in
        t(:,i,j) = theta_in * (pressure_in/100000.)**0.286
        p(:,i,j) = pressure_in
        dz(:,i,j) = dz_in
        w(:,i,j) = w_in
        qc(:,i,j) = qc_in
        nc(:,i,j) = nc_in
        qr(:,i,j) = qr_in
        nr(:,i,j) = nr_in
        qi(:,i,j) =  qi_in
        qs(:,i,j) = qs_in
        ni(:,i,j) = ni_in
      end do
    end do
    call nvtx_range_pop()

    call nvtx_range_push('test_ml_cloud_effective_radius_bl_qc')
    do j = jts, jte
      do i = its, ite
        qc_bl(:,i,j) = 0._wp
        qcfrac_bl(:,i,j) = 0._wp
        do k = 1, nz
          if (k < 20) then
            qc_bl(k,i,j) = qc_in(k+15)*0.2_wp
          endif
          if (qc_bl(k,i,j) > 0._wp) then
            qcfrac_bl(k,i,j) = 0.35_wp
          endif
        enddo
      end do
    end do
    call nvtx_range_pop()

    qg = 0._wp
    ng = 0._wp
    qb = 0._wp

    ! set configs
    tempo_cfgs%turn_off_micro_flag = .true.
    tempo_cfgs%ml_for_bl_nc_flag = .true.
    tempo_cfgs%cloud_number_mixing_ratio_flag = .true.

    total_timesteps = int(integration_time/dt)

    call tempo_tend_init(kts, kte, its, ite, jts, jte)

    call nvtx_range_push('tempo_run_enter_data')
    call tempo_run_enter_data(tempo_cfgs=tempo_cfgs, tempo_diags=tempo_driver_diags, use_temperature=.true., &
      t=t, p=p, w=w, dz=dz, qv=qv, th=th, pii=pii, qc=qc, nc=nc, qr=qr, qi=qi, qs=qs, qg=qg, ni=ni, nr=nr, &
      qc_bl=qc_bl, qcfrac_bl=qcfrac_bl, ng=ng, qb=qb, its=its, ite=ite, jts=jts, jte=jte, kts=kts, kte=kte)
    call nvtx_range_pop()

    call nvtx_range_push('test_ml_cloud_effective_radius_timestep_integrate')
    elapsed_time = 0._wp
    do itimestep = 1, total_timesteps
      call nvtx_range_push('test_ml_cloud_effective_radius_timestep')
      call system_clock(count_start, count_rate)
      call nvtx_range_push('test_ml_cloud_effective_radius_tempo_run')
      call  tempo_run(tempo_cfgs=tempo_cfgs, itimestep=itimestep, dt=dt, &
                      ids=ids, ide=ide, ims=ims, ime=ime, its=its, ite=ite, &
                      jds=jds, jde=jde, jms=jms, jme=jme, jts=jts, jte=jte, &
                      kds=kds, kde=kde, kms=kms, kme=kme, kts=kts, kte=kte, &
                      t=t, p=p, w=w, dz=dz, qv=qv, th=th, pii=pii, &
                      qc_bl=qc_bl, qcfrac_bl=qcfrac_bl, &
                      qc=qc, nc=nc, qr=qr, qi=qi, qs=qs, qg=qg, ni=ni, nr=nr, &
                      tempo_diags=tempo_driver_diags, stride=effective_stride(nCells), arguments_on_device=.true.)
      call nvtx_range_pop()
      call system_clock(count_end)
      if (count_rate > 0) elapsed_time = elapsed_time + real(count_end - count_start, kind=wp) / real(count_rate, wp)

      if (itimestep == 1) then
        call tempo_run_sync_host_diags(tempo_driver_diags)
        call tempo_run_sync_host_fields(use_temperature=.true., t=t, qv=qv, qc=qc, qr=qr, qi=qi, qs=qs, qg=qg, &
          ni=ni, nr=nr, ng=ng, qb=qb, nc=nc)
        open(newunit=io3, file="cloud_re_init.txt", status="new", action="write")
        n_samp = max(1, nOutValues)
        do ip = 1, n_samp
          if (n_samp == 1) then
            ij_samp = 1
          else
            ij_samp = 1 + nint(real((ip - 1) * (nCells - 1), wp) / real(n_samp - 1, wp))
          endif
          ij_samp = max(1, min(nCells, ij_samp))
          i = ij_samp
          j = ij_samp
          write(io3,'(A,2I8)') '# i j ', i, j
          do k = 1, nz
            write(io3,'(I5, 5E12.4)') k, qc(k,i,j), qc_bl(k,i,j), qcfrac_bl(k,i,j), &
              tempo_driver_diags%cloud_number_mixing_ratio(k,i,j), tempo_driver_diags%re_cloud(k,i,j)*1.e6_wp
          enddo
        enddo
      endif 
      call nvtx_range_pop()
    enddo

    call nvtx_range_pop()

    call nvtx_range_push('tempo_run_exit_data')
    call tempo_run_exit_data(tempo_cfgs=tempo_cfgs, tempo_diags=tempo_driver_diags, use_temperature=.true., &
      t=t, p=p, w=w, dz=dz, qv=qv, th=th, pii=pii, qc=qc, nc=nc, qr=qr, qi=qi, qs=qs, qg=qg, ni=ni, nr=nr, &
      qc_bl=qc_bl, qcfrac_bl=qcfrac_bl, ng=ng, qb=qb, kts=kts, kte=kte, its=its, ite=ite, jts=jts, jte=jte)
    call nvtx_range_pop()

    call tempo_tend_finalize()

    write(tt_string, '(I7)') integration_time
    open(newunit=io4, file="cloud_re_dt_"//trim(adjustl(dt_string))// &
      "_runtime"//trim(adjustl(tt_string))//".txt", status="new", action="write")
    n_samp = max(1, nOutValues)
    do ip = 1, n_samp
      if (n_samp == 1) then
        ij_samp = 1
      else
        ij_samp = 1 + nint(real((ip - 1) * (nCells - 1), wp) / real(n_samp - 1, wp))
      endif
      ij_samp = max(1, min(nCells, ij_samp))
      i = ij_samp
      j = ij_samp
      write(io4,'(A,2I8)') '# i j ', i, j
      do k = 1, nz
        write(io4,'(I5, 5E12.4)') k, qc(k,i,j), qc_bl(k,i,j), qcfrac_bl(k,i,j), &
          tempo_driver_diags%cloud_number_mixing_ratio(k,i,j), tempo_driver_diags%re_cloud(k,i,j)*1.e6_wp
      enddo
    enddo
    write(io4,'(A,F12.4,A)') 'test_ml_cloud_effective_radius compute time: ', elapsed_time, ' s'
    write(*,'(A,F12.4,A)') 'test_ml_cloud_effective_radius compute time: ', elapsed_time, ' s'

    call nvtx_range_pop()
    deallocate(qv, t, th, pii, p, w, dz, qc, qr, qi, qs, qg, ni, nr, ng, qb, nc, qc_bl, qcfrac_bl)
  end subroutine test_ml_cloud_effective_radius

end module tests