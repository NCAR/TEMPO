module module_mp_tempo_init
  !! initialize variables for tempo microphysics
  !! includes a procedure to build and save lookup tables
  use module_mp_tempo_params, only : wp, sp, dp, tempo_init_cfgs, tempo_table_cfgs
  use module_mp_tempo_utils, only : snow_moments, calc_gamma_p, get_nuc
  ! use module_mp_radar

#ifdef build_tables_with_mpi
  use mpi_f08 
#endif

  implicit none
  private

#ifdef unit_testing
  public :: compute_efrw, compute_efsw
#endif

  public :: tempo_init, tempo_build_tables

  ! type(ty_tempo_init_cfgs) :: tempo_init_cfgs
  ! type(ty_tempo_table_cfgs) :: tempo_table_cfgs

  contains

  subroutine tempo_init(aerosolaware_flag, hailaware_flag, restart_flag)
    !! public procedure called to initialize tempo microphysics
    use module_mp_tempo_params, only : tempo_version, t_efrw, &
      initialize_graupel_vars, initialize_parameters, initialize_bins_for_tables, &
      initialize_array_efrw, initialize_array_efsw, initialize_arrays_drop_evap, initialize_arrays_ccn, initialize_arrays_qi_aut_qs, &
      initialize_arrays_qr_acr_qs, initialize_arrays_qr_acr_qg, initialize_arrays_freezewater, &
      am_r, bm_r, mu_r, am_s, bm_s, mu_s, am_g, bm_g, mu_g, idx_bg1 ! For radar

    logical, intent(in), optional :: aerosolaware_flag, hailaware_flag, restart_flag

    character(len=100) :: table_filename
    integer :: table_size
    logical :: initialize_mp_vars

    ! get tempo version from readme file
    call get_version(tempo_version) 

    ! check an allocatable array (t_efrw) to see if initialization can be skipped
    initialize_mp_vars = .true.
    if (allocated(t_efrw)) initialize_mp_vars = .false.

    if (initialize_mp_vars) then
      if (present(aerosolaware_flag)) tempo_init_cfgs%aerosolaware_flag = aerosolaware_flag
      if (present(hailaware_flag)) tempo_init_cfgs%hailaware_flag = hailaware_flag
      if (present(restart_flag)) tempo_init_cfgs%restart_flag = restart_flag

      write(*,'(A)') 'tempo_init() --- TEMPO microphysics configuration options: '
      write(*,'(A,L)') 'tempo_init() --- aerosol aware = ', tempo_init_cfgs%aerosolaware_flag
      write(*,'(A,L)') 'tempo_init() --- hail aware = ', tempo_init_cfgs%hailaware_flag
      write(*,'(A,L)') 'tempo_init() --- restart = ', tempo_init_cfgs%restart_flag

      ! set graupel variables from hail_aware_flag
      call initialize_graupel_vars(tempo_init_cfgs%hailaware_flag) 
      write(*,'(A,L)') 'tempo_init() --- initialized graupel variables using hail aware = ', tempo_init_cfgs%hailaware_flag

      ! set parameters that can depend on the host model
      call initialize_parameters() 
      write(*,'(A)') 'tempo_init() --- initialized parameters'

      ! creates log-spaced bins of hydrometers for tables
      call initialize_bins_for_tables() 
      write(*,'(A)') 'tempo_init() --- initialized bins for lookup tables'

      ! collision efficiencies between rain/snow and cloud water.
      call initialize_array_efrw()
      call compute_efrw()
      write(*,'(A)') 'tempo_init() --- initialized collision efficiency data for rain collecting cloud water'
      call initialize_array_efsw()
      call compute_efsw()
      write(*,'(A)') 'tempo_init() --- initialized collision efficiency data for snow collecting cloud water'

      ! drop evaporation
      call initialize_arrays_drop_evap()
      call compute_drop_evap()
      write(*,'(A)') 'tempo_init() --- initialized drop evaporation data'

      ! cloud ice to snow and depositional growth
      ! consider explicitly passing arrays at init
      call initialize_arrays_qi_aut_qs()
      call qi_aut_qs()

      ! CCN activation table
      table_filename = tempo_table_cfgs%ccn_table_name
      call initialize_arrays_ccn(table_size)
      call read_table_ccn(trim(table_filename), table_size)
      write(*,'(A)') 'tempo_init() --- initialized data for ccn lookup table'

      ! freeze water collection lookup table
      table_filename = tempo_table_cfgs%freezewater_table_name
      call initialize_arrays_freezewater(table_size)
      call read_table_freezewater(trim(table_filename), table_size)
      write(*,'(A)') 'tempo_init() --- initialized data for frozen cloud water and rain lookup table'

      ! rain-snow collection lookup table
      table_filename = tempo_table_cfgs%qrqs_table_name
      call initialize_arrays_qr_acr_qs(table_size)
      call read_table_qr_acr_qs(trim(table_filename), table_size)
      write(*,'(A)') 'tempo_init() --- initialized data for rain-snow collection lookup table'

      ! rain-graupel collection lookup table
      table_filename = tempo_table_cfgs%qrqg_table_name
      call initialize_arrays_qr_acr_qg(table_size)
      call read_table_qr_acr_qg(trim(table_filename), table_size)
      write(*,'(A)') 'tempo_init() --- initialized data for rain-snow collection lookup table'
     
      ! Initialize various constants for computing radar reflectivity.
      ! xam_r = am_r
      ! xbm_r = bm_r
      ! xmu_r = mu_r
      ! xam_s = am_s
      ! xbm_s = bm_s
      ! xmu_s = mu_s
      ! xam_g = am_g(idx_bg1)
      ! xbm_g = bm_g
      ! xmu_g = mu_g
    endif
  end subroutine tempo_init


  subroutine tempo_build_tables(build_tables_rank, build_tables_num_proc)
    !! public procedure to build 3 lookup tables for tempo microphysics
    use module_mp_tempo_params, only : tempo_version, &
      initialize_graupel_vars, initialize_parameters, initialize_bins_for_tables

    integer, intent(in) :: build_tables_rank, build_tables_num_proc

    character(len=100) :: table_filename
    logical, parameter :: build_table_hail_flag = .true.

    ! MPI to speed up the table building process
    integer :: rank, num_proc

    ! set global variables rank and num_proc from MPI_Comm in build_tables program 
    rank = build_tables_rank
    num_proc = build_tables_num_proc

    ! get tempo version from readme file
    call get_version(tempo_version) 

#ifdef build_tables_with_mpi
      write(*,'(A,I4,A)') 'tempo_build_tables() --- building lookup tables with MPI and', num_proc, ' process(es)' 
#else
      write(*,'(A,I4,A)') 'tempo_build_tables() --- building lookup tables with', num_proc, ' process' 
#endif

      ! hard-code hail aware = true to build lookup tables
      call initialize_graupel_vars(build_table_hail_flag) 
      write(*,'(A,L)') 'tempo_build_tables() --- initialized graupel variables using hail aware = ', build_table_hail_flag

      ! set parameters that can depend on the host model
      call initialize_parameters() 
      write(*,'(A)') 'tempo_build_tables() --- initialized parameters'

      ! creates log-spaced bins of hydrometers for tables
      call initialize_bins_for_tables() 
      write(*,'(A)') 'tempo_build_tables() --- initialized bins for lookup tables'

      ! freeze water collection lookup table
      table_filename = tempo_table_cfgs%freezewater_table_name
      write(*,'(2A)') 'tempo_build_tables() --- building table ', trim(table_filename)
      call build_table_freezewater()
      if (rank == 0) call write_table_freezewater(trim(table_filename))

      ! rain-snow collection lookup table
      table_filename = tempo_table_cfgs%qrqs_table_name
      write(*,'(2A)') 'tempo_build_tables() --- building table ', trim(table_filename)
      call build_table_qr_acr_qs(rank, num_proc)
      if (rank == 0) call write_table_qr_acr_qs(trim(table_filename))

      ! rain-graupel collection lookup table
      table_filename = tempo_table_cfgs%qrqg_table_name
      write(*,'(2A)') 'tempo_build_tables() --- building table ', trim(table_filename)
      call build_table_qr_acr_qg(rank, num_proc)
      if (rank == 0) call write_table_qr_acr_qg(trim(table_filename))
  end subroutine tempo_build_tables


  subroutine build_table_freezewater()
    !! build lookup table data for frozen cloud and rain water
    use module_mp_tempo_params, only : initialize_arrays_freezewater

    real(wp) :: timing_start, timing_end

    call initialize_arrays_freezewater()
    call cpu_time(timing_start)
    call freezewater()
    call cpu_time(timing_end)
    write(*,'(A,I5,A)') 'build_table_freezewater() --- time to build table: ', int(timing_end - timing_start), ' s'
  end subroutine build_table_freezewater 


 subroutine write_table_freezewater(filename)
    !! write data for frozen cloud water and rain to a file
    use module_mp_tempo_params, only : tpi_qrfz, tni_qrfz, tpg_qrfz, &
      tnr_qrfz, tpi_qcfz, tni_qcfz

    character(len=*), intent(in) :: filename
    integer :: mp_unit, istat
    logical :: fileexists

    inquire(file=filename, exist=fileexists)
    if (fileexists) call rename_file_if_exists(filename)

    mp_unit = 11
    open(unit=mp_unit, file=filename, form='unformatted', status='new', &
      iostat=istat, convert='big_endian')
    write(mp_unit) tpi_qrfz
    write(mp_unit) tni_qrfz
    write(mp_unit) tpg_qrfz
    write(mp_unit) tnr_qrfz
    write(mp_unit) tpi_qcfz
    write(mp_unit) tni_qcfz
    close(unit=mp_unit)
  end subroutine write_table_freezewater


  subroutine read_table_freezewater(filename, table_size)
    !! read lookup table for rain-graupel collection
    use module_mp_tempo_params, only : tpi_qrfz, tni_qrfz, &
      tpg_qrfz, tnr_qrfz, tpi_qcfz, tni_qcfz

    character(len=*), intent(in) :: filename
    integer, intent(in) :: table_size
    
    integer :: mp_unit, istat

    mp_unit = 11
    call check_before_table_read(filename, table_size)
    open(unit=mp_unit, file=filename, form='unformatted', status='old', &
      action='read', iostat=istat, convert='big_endian')
    read(mp_unit) tpi_qrfz
    read(mp_unit) tni_qrfz
    read(mp_unit) tpg_qrfz
    read(mp_unit) tnr_qrfz
    read(mp_unit) tpi_qcfz
    read(mp_unit) tni_qcfz
    close(unit=mp_unit)
  end subroutine read_table_freezewater


  subroutine build_table_qr_acr_qs(rank, num_proc)
    !! build lookup table data for rain-snow collection
    use module_mp_tempo_params, only : table_dp, initialize_arrays_qr_acr_qs, &
      tcs_racs1, tmr_racs1, tcs_racs2, tmr_racs2, tcr_sacr1, tms_sacr1, & ! data arrays
      tcr_sacr2, tms_sacr2, tnr_racs1, tnr_racs2, tnr_sacr1, tnr_sacr2, & ! data arrays
      ntb_s, ntb_t, ntb_r1, ntb_r ! dimensions
  
    integer, intent(in) :: rank, num_proc
#ifdef build_tables_with_mpi
    integer :: ierror
#endif

    real(wp) :: timing_start, timing_end
    integer :: start_idx, end_idx, local_dim_size, local_flat_size
    integer, allocatable, dimension(:) :: sendcounts, displacements
    real(table_dp), allocatable, dimension(:) :: tcs_racs1_flat, tmr_racs1_flat, &
      tcs_racs2_flat, tmr_racs2_flat, tcr_sacr1_flat, tms_sacr1_flat, &
      tcr_sacr2_flat, tms_sacr2_flat, tnr_racs1_flat, tnr_racs2_flat, &
      tnr_sacr1_flat, tnr_sacr2_flat
    real(table_dp), allocatable, dimension(:,:,:,:) :: tcs_racs1_, tmr_racs1_, &
      tcs_racs2_, tmr_racs2_, tcr_sacr1_, tms_sacr1_, tcr_sacr2_, &
      tms_sacr2_, tnr_racs1_, tnr_racs2_, tnr_sacr1_, tnr_sacr2_

    if (rank == 0) then
      ! initialize lookup table arrays and flatten for MPI
      call initialize_arrays_qr_acr_qs()
      allocate(tcs_racs1_flat(size(tcs_racs1)))
      allocate(tmr_racs1_flat(size(tmr_racs1)))
      allocate(tcs_racs2_flat(size(tcs_racs2)))
      allocate(tmr_racs2_flat(size(tmr_racs2)))
      allocate(tcr_sacr1_flat(size(tcr_sacr1)))
      allocate(tms_sacr1_flat(size(tms_sacr1)))
      allocate(tcr_sacr2_flat(size(tcr_sacr2)))
      allocate(tms_sacr2_flat(size(tms_sacr2)))
      allocate(tnr_racs1_flat(size(tnr_racs1)))
      allocate(tnr_racs2_flat(size(tnr_racs2)))
      allocate(tnr_sacr1_flat(size(tnr_sacr1)))
      allocate(tnr_sacr2_flat(size(tnr_sacr2)))
    endif 

    ! split over the last dimension, nrb_r
    call get_index_for_rank(ntb_r, rank, num_proc, start_idx, end_idx)    
    local_dim_size = end_idx - start_idx + 1
    local_flat_size = local_dim_size * ntb_s*ntb_t*ntb_r1

    ! local arrays for MPI
    allocate(tcs_racs1_(ntb_s,ntb_t,ntb_r1,local_dim_size))
    allocate(tmr_racs1_(ntb_s,ntb_t,ntb_r1,local_dim_size))
    allocate(tcs_racs2_(ntb_s,ntb_t,ntb_r1,local_dim_size))
    allocate(tmr_racs2_(ntb_s,ntb_t,ntb_r1,local_dim_size))
    allocate(tcr_sacr1_(ntb_s,ntb_t,ntb_r1,local_dim_size))
    allocate(tms_sacr1_(ntb_s,ntb_t,ntb_r1,local_dim_size))
    allocate(tcr_sacr2_(ntb_s,ntb_t,ntb_r1,local_dim_size))
    allocate(tms_sacr2_(ntb_s,ntb_t,ntb_r1,local_dim_size))
    allocate(tnr_racs1_(ntb_s,ntb_t,ntb_r1,local_dim_size))
    allocate(tnr_racs2_(ntb_s,ntb_t,ntb_r1,local_dim_size))
    allocate(tnr_sacr1_(ntb_s,ntb_t,ntb_r1,local_dim_size))
    allocate(tnr_sacr2_(ntb_s,ntb_t,ntb_r1,local_dim_size))

    allocate(sendcounts(num_proc), displacements(num_proc))
    if (num_proc == 1) then
      sendcounts(1) = local_flat_size
      displacements(1) = 0 ! MPI displacements start at zero
    endif 

#ifdef build_tables_with_mpi
    call MPI_Allgather(local_flat_size, 1, MPI_INTEGER, sendcounts, 1, MPI_INTEGER, MPI_COMM_WORLD, ierror)
    ! MPI displacements start at zero, i.e., start_idx-1 below
    call MPI_Allgather((start_idx-1)*ntb_s*ntb_t*ntb_r1, 1, MPI_INTEGER, displacements, 1, MPI_INTEGER, MPI_COMM_WORLD, ierror)
#endif

#ifdef build_tables_with_mpi
    timing_start = MPI_Wtime()
#else
    call cpu_time(timing_start)
#endif

    call qr_acr_qs(start_idx, end_idx, &
      tcs_racs1_, tmr_racs1_, tcs_racs2_, tmr_racs2_, &
      tcr_sacr1_, tms_sacr1_, tcr_sacr2_, tms_sacr2_, &
      tnr_racs1_, tnr_racs2_, tnr_sacr1_, tnr_sacr2_)

#ifdef build_tables_with_mpi
    timing_end = MPI_Wtime()
#else
    call cpu_time(timing_end)
#endif
    if (rank == 0) write(*,'(A,I5,A)') 'build_table_qr_acr_qs() --- time to build table: ', int(timing_end - timing_start), ' s'

#ifdef build_tables_with_mpi
    call MPI_Barrier(MPI_COMM_WORLD, ierror)
    call gather(reshape(tcs_racs1_, (/local_flat_size/)), tcs_racs1_flat, sendcounts, displacements, ierror)
    call gather(reshape(tmr_racs1_, (/local_flat_size/)), tmr_racs1_flat, sendcounts, displacements, ierror)
    call gather(reshape(tcs_racs2_, (/local_flat_size/)), tcs_racs2_flat, sendcounts, displacements, ierror)
    call gather(reshape(tmr_racs2_, (/local_flat_size/)), tmr_racs2_flat, sendcounts, displacements, ierror)
    call gather(reshape(tcr_sacr1_, (/local_flat_size/)), tcr_sacr1_flat, sendcounts, displacements, ierror)
    call gather(reshape(tms_sacr1_, (/local_flat_size/)), tms_sacr1_flat, sendcounts, displacements, ierror)
    call gather(reshape(tcr_sacr2_, (/local_flat_size/)), tcr_sacr2_flat, sendcounts, displacements, ierror)
    call gather(reshape(tms_sacr2_, (/local_flat_size/)), tms_sacr2_flat, sendcounts, displacements, ierror)
    call gather(reshape(tnr_racs1_, (/local_flat_size/)), tnr_racs1_flat, sendcounts, displacements, ierror)
    call gather(reshape(tnr_racs2_, (/local_flat_size/)), tnr_racs2_flat, sendcounts, displacements, ierror)
    call gather(reshape(tnr_sacr1_, (/local_flat_size/)), tnr_sacr1_flat, sendcounts, displacements, ierror)
    call gather(reshape(tnr_sacr2_, (/local_flat_size/)), tnr_sacr2_flat, sendcounts, displacements, ierror)

    if (rank == 0) then
      tcs_racs1 = reshape(tcs_racs1_flat, (/ntb_s,ntb_t,ntb_r1,ntb_r/))
      tmr_racs1 = reshape(tmr_racs1_flat, (/ntb_s,ntb_t,ntb_r1,ntb_r/))
      tcs_racs2 = reshape(tcs_racs2_flat, (/ntb_s,ntb_t,ntb_r1,ntb_r/))
      tmr_racs2 = reshape(tmr_racs2_flat, (/ntb_s,ntb_t,ntb_r1,ntb_r/))
      tcr_sacr1 = reshape(tcr_sacr1_flat, (/ntb_s,ntb_t,ntb_r1,ntb_r/))
      tms_sacr1 = reshape(tms_sacr1_flat, (/ntb_s,ntb_t,ntb_r1,ntb_r/))
      tcr_sacr2 = reshape(tcr_sacr2_flat, (/ntb_s,ntb_t,ntb_r1,ntb_r/))
      tms_sacr2 = reshape(tms_sacr2_flat, (/ntb_s,ntb_t,ntb_r1,ntb_r/))
      tnr_racs1 = reshape(tnr_racs1_flat, (/ntb_s,ntb_t,ntb_r1,ntb_r/))
      tnr_racs2 = reshape(tnr_racs2_flat, (/ntb_s,ntb_t,ntb_r1,ntb_r/))
      tnr_sacr1 = reshape(tnr_sacr1_flat, (/ntb_s,ntb_t,ntb_r1,ntb_r/))
      tnr_sacr2 = reshape(tnr_sacr2_flat, (/ntb_s,ntb_t,ntb_r1,ntb_r/))
    endif 
#else
    tcs_racs1 = tcs_racs1_
    tmr_racs1 = tmr_racs1_
    tcs_racs2 = tcs_racs2_
    tmr_racs2 = tmr_racs2_
    tcr_sacr1 = tcr_sacr1_
    tms_sacr1 = tms_sacr1_
    tcr_sacr2 = tcr_sacr2_
    tms_sacr2 = tms_sacr2_
    tnr_racs1 = tnr_racs1_
    tnr_racs2 = tnr_racs2_
    tnr_sacr1 = tnr_sacr1_
    tnr_sacr2 = tnr_sacr2_
#endif                                                                                       

#if build_tables_with_mpi
    call MPI_Barrier(MPI_COMM_WORLD, ierror)
#endif
  end subroutine build_table_qr_acr_qs


  subroutine write_table_qr_acr_qs(filename)
    !! write data for rain-snow collection to a file
    use module_mp_tempo_params, only : tcs_racs1, tmr_racs1, tcs_racs2, &
      tmr_racs2, tcr_sacr1, tms_sacr1, tcr_sacr2, tms_sacr2, tnr_racs1, &
      tnr_racs2, tnr_sacr1, tnr_sacr2

    character(len=*), intent(in) :: filename
    integer :: mp_unit, istat
    logical :: fileexists

    inquire(file=filename, exist=fileexists)
    if (fileexists) call rename_file_if_exists(filename)

    mp_unit = 11
    open(unit=mp_unit, file=filename, form='unformatted', status='new', &
      iostat=istat, convert='big_endian')
    write(mp_unit) tcs_racs1
    write(mp_unit) tmr_racs1
    write(mp_unit) tcs_racs2
    write(mp_unit) tmr_racs2
    write(mp_unit) tcr_sacr1
    write(mp_unit) tms_sacr1
    write(mp_unit) tcr_sacr2
    write(mp_unit) tms_sacr2
    write(mp_unit) tnr_racs1
    write(mp_unit) tnr_racs2
    write(mp_unit) tnr_sacr1
    write(mp_unit) tnr_sacr2
    close(unit=mp_unit)
  end subroutine write_table_qr_acr_qs


  subroutine read_table_qr_acr_qs(filename, table_size)
    !! read lookup table for rain-snow collection
    use module_mp_tempo_params, only : tcs_racs1, tmr_racs1, &
      tcs_racs2, tmr_racs2, tcr_sacr1, tms_sacr1, tcr_sacr2, &
      tms_sacr2, tnr_racs1, tnr_racs2, tnr_sacr1, tnr_sacr2

    character(len=*), intent(in) :: filename
    integer, intent(in) :: table_size
    
    integer :: mp_unit, istat

    mp_unit = 11
    call check_before_table_read(filename, table_size)
    open(unit=mp_unit, file=filename, form='unformatted', status='old', &
      action='read', iostat=istat, convert='big_endian')
    read(mp_unit) tcs_racs1
    read(mp_unit) tmr_racs1
    read(mp_unit) tcs_racs2
    read(mp_unit) tmr_racs2
    read(mp_unit) tcr_sacr1
    read(mp_unit) tms_sacr1
    read(mp_unit) tcr_sacr2
    read(mp_unit) tms_sacr2
    read(mp_unit) tnr_racs1
    read(mp_unit) tnr_racs2
    read(mp_unit) tnr_sacr1
    read(mp_unit) tnr_sacr2
    close(unit=mp_unit)
  end subroutine read_table_qr_acr_qs


  subroutine build_table_qr_acr_qg(rank, num_proc)
    !! build lookup table data for rain-graupel collection
    use module_mp_tempo_params, only : table_dp, initialize_arrays_qr_acr_qg, &
      tcg_racg, tmr_racg, tcr_gacr, tnr_racg, tnr_gacr, &
      ntb_g1, ntb_g, nrhg, ntb_r1, ntb_r
  
    integer, intent(in) :: rank, num_proc
#ifdef build_tables_with_mpi
    integer :: ierror
#endif

    real(wp) :: timing_start, timing_end
    integer :: start_idx, end_idx, local_dim_size, local_flat_size
    integer, allocatable, dimension(:) :: sendcounts, displacements
    real(table_dp), allocatable, dimension(:) ::  tcg_racg_flat, tmr_racg_flat, &
      tcr_gacr_flat, tnr_racg_flat, tnr_gacr_flat
    real(table_dp), allocatable, dimension(:,:,:,:,:) :: tcg_racg_, &
      tmr_racg_, tcr_gacr_, tnr_racg_, tnr_gacr_

    if (rank == 0) then
      ! initialize lookup table arrays and flatten for MPI
      call initialize_arrays_qr_acr_qg()
      allocate(tcg_racg_flat(size(tcg_racg)))
      allocate(tmr_racg_flat(size(tmr_racg)))
      allocate(tcr_gacr_flat(size(tcr_gacr)))
      allocate(tnr_racg_flat(size(tnr_racg)))
      allocate(tnr_gacr_flat(size(tnr_gacr)))
    endif 

    ! split over the last dimension, nrb_r
    call get_index_for_rank(ntb_r, rank, num_proc, start_idx, end_idx)    
    local_dim_size = end_idx - start_idx + 1
    local_flat_size = local_dim_size * ntb_g1*ntb_g*nrhg*ntb_r1

    ! local arrays for MPI
    allocate(tcg_racg_(ntb_g1,ntb_g,nrhg,ntb_r1,local_dim_size))
    allocate(tmr_racg_(ntb_g1,ntb_g,nrhg,ntb_r1,local_dim_size))
    allocate(tcr_gacr_(ntb_g1,ntb_g,nrhg,ntb_r1,local_dim_size))
    allocate(tnr_racg_(ntb_g1,ntb_g,nrhg,ntb_r1,local_dim_size))
    allocate(tnr_gacr_(ntb_g1,ntb_g,nrhg,ntb_r1,local_dim_size))

    allocate(sendcounts(num_proc), displacements(num_proc))
    if (num_proc == 1) then
      sendcounts(1) = local_flat_size
      displacements(1) = 0 ! MPI displacements start at zero
    endif 

#ifdef build_tables_with_mpi
    call MPI_Allgather(local_flat_size, 1, MPI_INTEGER, sendcounts, 1, MPI_INTEGER, MPI_COMM_WORLD, ierror)
    ! MPI displacements start at zero, i.e., start_idx-1 below
    call MPI_Allgather((start_idx-1)*ntb_g1*ntb_g*nrhg*ntb_r1, 1, MPI_INTEGER, displacements, 1, MPI_INTEGER, MPI_COMM_WORLD, ierror)
#endif

#ifdef build_tables_with_mpi
    timing_start = MPI_Wtime()
#else
    call cpu_time(timing_start)
#endif

    call qr_acr_qg(start_idx, end_idx, &
      tcg_racg_, tmr_racg_, tcr_gacr_, tnr_racg_, tnr_gacr_)

#ifdef build_tables_with_mpi
    timing_end = MPI_Wtime()
#else
    call cpu_time(timing_end)
#endif
    if (rank == 0) write(*,'(A,I5,A)') 'build_table_qr_acr_qg() --- time to build table: ', int(timing_end - timing_start), ' s'

#ifdef build_tables_with_mpi
    call MPI_Barrier(MPI_COMM_WORLD, ierror)
    call gather(reshape(tcg_racg_, (/local_flat_size/)), tcg_racg_flat, sendcounts, displacements, ierror)
    call gather(reshape(tmr_racg_, (/local_flat_size/)), tmr_racg_flat, sendcounts, displacements, ierror)
    call gather(reshape(tcr_gacr_, (/local_flat_size/)), tcr_gacr_flat, sendcounts, displacements, ierror)
    call gather(reshape(tnr_racg_, (/local_flat_size/)), tnr_racg_flat, sendcounts, displacements, ierror)
    call gather(reshape(tnr_gacr_, (/local_flat_size/)), tnr_gacr_flat, sendcounts, displacements, ierror)

    if (rank == 0) then
      tcg_racg = reshape(tcg_racg_flat, (/ntb_g1,ntb_g,nrhg,ntb_r1,ntb_r/))
      tmr_racg = reshape(tmr_racg_flat, (/ntb_g1,ntb_g,nrhg,ntb_r1,ntb_r/))
      tcr_gacr = reshape(tcr_gacr_flat, (/ntb_g1,ntb_g,nrhg,ntb_r1,ntb_r/))
      tnr_racg = reshape(tnr_racg_flat, (/ntb_g1,ntb_g,nrhg,ntb_r1,ntb_r/))
      tnr_gacr = reshape(tnr_gacr_flat, (/ntb_g1,ntb_g,nrhg,ntb_r1,ntb_r/))
    endif 
#else
    tcg_racg = tcg_racg_
    tmr_racg = tmr_racg_
    tcr_gacr = tcr_gacr_
    tnr_racg = tnr_racg_
    tnr_gacr = tnr_gacr_
#endif                                                                                       

#if build_tables_with_mpi
    call MPI_Barrier(MPI_COMM_WORLD, ierror)
#endif
  end subroutine build_table_qr_acr_qg


  subroutine write_table_qr_acr_qg(filename)
    !! write data for rain-graupel collection to a file
    use module_mp_tempo_params, only : tcg_racg, tmr_racg, &
      tcr_gacr, tnr_racg, tnr_gacr

    character(len=*), intent(in) :: filename
    integer :: mp_unit, istat
    logical :: fileexists

    inquire(file=filename, exist=fileexists)
    if (fileexists) call rename_file_if_exists(filename)

    mp_unit = 11
    open(unit=mp_unit, file=filename, form='unformatted', status='new', &
      iostat=istat, convert='big_endian')
    write(mp_unit) tcg_racg
    write(mp_unit) tmr_racg
    write(mp_unit) tcr_gacr
    write(mp_unit) tnr_racg
    write(mp_unit) tnr_gacr
    close(unit=mp_unit)
  end subroutine write_table_qr_acr_qg


  subroutine read_table_qr_acr_qg(filename, table_size)
    !! read lookup table for rain-graupel collection
    use module_mp_tempo_params, only : tcg_racg, tmr_racg, &
      tcr_gacr, tnr_racg, tnr_gacr

    character(len=*), intent(in) :: filename
    integer, intent(in) :: table_size
    
    integer :: mp_unit, istat

    mp_unit = 11
    call check_before_table_read(filename, table_size)
    open(unit=mp_unit, file=filename, form='unformatted', status='old', &
      action='read', iostat=istat, convert='big_endian')
    read(mp_unit) tcg_racg
    read(mp_unit) tmr_racg
    read(mp_unit) tcr_gacr
    read(mp_unit) tnr_racg
    read(mp_unit) tnr_gacr
    close(unit=mp_unit)
  end subroutine read_table_qr_acr_qg

#ifdef build_tables_with_mpi
  subroutine gather(local_flat, global_flat, sendcounts, displacements, ierror)
    !! wrapper to simplify MPI_Gatherv
    use module_mp_tempo_params, only : table_dp

    real(table_dp), dimension(:), intent(in) :: local_flat
    real(table_dp), dimension(:), intent(out) :: global_flat
    integer, intent(inout) :: ierror
    integer, dimension(:), intent(in) :: sendcounts, displacements
    integer :: local_size

    local_size = size(local_flat)
    call MPI_Gatherv(local_flat, local_size, MPI_DOUBLE_PRECISION, &
      global_flat, sendcounts, displacements, &
      MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
  end subroutine gather
#endif

  subroutine get_index_for_rank(idx, rank, num_proc, start_idx, end_idx)
    !! returns start and end index values for an array dimension
    !! of size idx distributed over num_procs
    
    integer, intent(in) :: idx, rank, num_proc
    integer, intent(out) :: start_idx, end_idx
    integer :: values_per_proc

    if (num_proc > idx) then
      write(*,'(A,I4,A,I4)') 'num_proc', num_proc, 'cannot be larger than idx', idx
      error stop '--- reduce the number of processes'
    endif

    values_per_proc = idx/num_proc

    if(rank < mod(idx, num_proc)) then
      values_per_proc = values_per_proc + 1
    endif

    if(rank < mod(idx, num_proc)) then
      start_idx = rank * values_per_proc + 1
      end_idx = (rank+1) * values_per_proc
    else
      start_idx = mod(idx, num_proc) + rank*values_per_proc + 1
      end_idx = mod(idx, num_proc) + (rank+1) * values_per_proc
    endif
  end subroutine get_index_for_rank


  subroutine get_version(version)
    !! returns the tempo version string from the README.md file or returns empty string if not found
  
    character(len=*), intent(inout) :: version
    character(len=100) :: first_line, filename
    integer :: io_unit
    logical :: fileexists

    filename = 'README.md'
    inquire(file=trim(filename), exist=fileexists)
    if (.not. fileexists) then
      version = ''
      write(*,'(A)') 'Unable to determine TEMPO Microphysics Version'
      return
    endif
  
    open(newunit=io_unit, file=filename, status='old', action='read')
    read(io_unit, '(A)') first_line
    close(io_unit)

    ! format is tempo-vX.X.X
    version = trim(first_line(8:))
    write(*,'(A)') 'TEMPO Microphysics Version: '//trim(version)
  end subroutine get_version


  subroutine read_table_ccn(filename, table_size)
    !! read static file containing CCN activation of aerosols;
    !! The data were created from a parcel model by Feingold & Heymsfield
    !! with further changes by Eidhammer and Kriedenweis
    use module_mp_tempo_params, only : tnccn_act
  
    character(len=*), intent(in) :: filename
    integer, intent(in) :: table_size
    
    integer :: mp_unit, istat

    call check_before_table_read(filename=filename, table_size=table_size)

    mp_unit = 11
    open(unit=mp_unit, file=filename, form='unformatted', status='old', &
      action='read', iostat=istat, convert='big_endian')
    read(mp_unit) tnccn_act
    close(unit=mp_unit)
  end subroutine read_table_ccn

  
  subroutine check_before_table_read(filename, table_size)
    !! checks that lookup tables exist and are the correct size before attempting to read them
    character(len=*), intent(in) :: filename
    integer, intent(in) :: table_size

    logical :: fileexists
    integer :: filesize
    character(len=100) :: int_to_str1, int_to_str2

    inquire(file=filename, size=filesize, exist=fileexists)
    if (.not. fileexists) then
      write(*,'(3A)') 'tempo_init() --- *** FATAL *** file "', filename, &
        '" was not found in this directory.'
      write(*,'(A)') ''
      write(*,'(A)') 'How to fix issues with tables (datasets stored in files):'
      write(*,'(3A)') '   (1) The table ', trim(tempo_table_cfgs%ccn_table_name), &
        ' is located in the TEMPO/tables/ directory. Copy this file to the directory where the model executable is located.'
      write(*,'(8A)') '   (2) Three tables:', trim(tempo_table_cfgs%qrqs_table_name), ', ', trim(tempo_table_cfgs%qrqg_table_name), ', and ', &
        trim(tempo_table_cfgs%freezewater_table_name), &
        ' can be build by compiling and running the executable "build_tables" in the main TEMPO directory. ', &
        'Then copy the file to the directory where the model executable is located.'
      write(*,'(A)') '   (3) Ask the developers for tables. They are willing to share.' 
      write(*,'(A)') ''      
      error stop '--- file "' // filename // '" needed for TEMPO microphysics was not found.'
    endif

    if (filesize /= table_size) then
      write(int_to_str1, '(I0)') filesize
      write(int_to_str2, '(I0)') table_size
      write(*,'(7A)') 'tempo_init() --- *** FATAL *** file "', filename, '" has a size of ', &
        trim(int_to_str1), ' bytes but the array allocated to hold the data expects a file size of ', &
          trim(int_to_str2), ' bytes.'
      write(*,'(A)') ''
      write(*,'(A)') 'How to fix issues with tables (datasets stored in files):'
      write(*,'(3A)') '   (1) The table ', trim(tempo_table_cfgs%ccn_table_name), &
        ' is located in the TEMPO/tables/ directory. Copy this file to the directory where the model executable is located.'
      write(*,'(8A)') '   (2) Three tables: ', trim(tempo_table_cfgs%qrqs_table_name), ', ', trim(tempo_table_cfgs%qrqg_table_name), ', and ', &
        trim(tempo_table_cfgs%freezewater_table_name), &
        ' can be build by compiling and running the executable "build_tables" in the main TEMPO directory. ', &
        'Then copy the file to the directory where the model executable is located.'
      write(*,'(A)') '   (3) Ask the developers for tables. They are willing to share.' 
      write(*,'(A)') ''
      error stop '--- size of file "' // filename // '" needed for TEMPO microphysics is inconsistent with expected size.'
    endif
  end subroutine check_before_table_read


  subroutine rename_file_if_exists(oldfilename)
    !! rename a lookup table if attempting to write to that file and it already exists
    character(len=*), intent(in) :: oldfilename
    character(len=100) :: newfilename
    character(8) :: date
    character(6) :: time
    integer :: renamestat
    character(len=20) :: fileappend

    call date_and_time(date, time)
    fileappend = '_old_' // trim(date) // '_' // trim(time)
    newfilename = trim(oldfilename // trim(fileappend))
    renamestat = rename(oldfilename, newfilename)
    if (renamestat /= 0) then
      write(*,'(3A)') 'rename_file_if_exists() --- *** FATAL *** unable to rename existing table file "', &
        oldfilename, '" to "', newfilename, '"'
      error stop '--- file rename failed.'
    else
    write(*,'(5A)') 'rename_file_if_exists() --- existing table file "', oldfilename, &
      '" has been renamed to "', newfilename, '"'
    endif
  end subroutine rename_file_if_exists


  subroutine compute_efrw()
    !! collision efficiency for rain collecting cloud water from Beard and Grover (1974)
    !! if a/A < 0.25, otherwise uses polynomials to get close match of Pruppacher & Klett Fig. 14-9
    use module_mp_tempo_params, only : nbc, nbr, dc, dr, t_efrw, rho_w, pi

    real(dp) :: vtr, stokes, reynolds, ef_rw
    real(dp) :: p, yc0, f, g, h, z, k0, x
    integer :: i, j

    do j = 1, nbc
      do i = 1, nbr
        ef_rw = 0.0_dp
        p = dc(j) / dr(i)
        if (dr(i) < 50.e-6_dp .or. dc(j) < 3.e-6_dp) then
          t_efrw(i,j) = 0.0_dp
        elseif (p > 0.25_dp) then
          x = dc(j) * 1.e6_dp
          if (dr(i) < 75.e-6_dp) then
            ef_rw = 0.026794_dp*x - 0.20604_dp
          elseif (dr(i) < 125.e-6_dp) then
            ef_rw = -0.00066842_dp*x*x + 0.061542_dp*x - 0.37089_dp
          elseif (dr(i) < 175.e-6_dp) then
            ef_rw = 4.091e-06_dp*x*x*x*x - 0.00030908_dp*x*x*x + &
              0.0066237_dp*x*x - 0.0013687_dp*x - 0.073022_dp
          elseif (dr(i) < 250.e-6_dp) then
            ef_rw = 9.6719e-5_dp*x*x*x - 0.0068901_dp*x*x + 0.17305_dp*x - 0.65988_dp
          elseif (dr(i) < 350.e-6_dp) then
            ef_rw = 9.0488e-5_dp*x*x*x - 0.006585_dp*x*x + 0.16606_dp*x - 0.56125_dp
          else
            ef_rw = 0.00010721_dp*x*x*x - 0.0072962_dp*x*x + 0.1704_dp*x - 0.46929_dp
          endif
        else
          vtr = -0.1021_dp + 4.932e3_dp*dr(i) - 0.9551e6_dp*dr(i)*dr(i) + 0.07934e9_dp*dr(i)*dr(i)*dr(i) &
            - 0.002362e12_dp*dr(i)*dr(i)*dr(i)*dr(i)
          stokes = dc(j) * dc(j) * vtr * rho_w / (9._dp*1.718e-5_dp*dr(i))
          reynolds = 9._dp * stokes / (p*p*rho_w)

          f = log(reynolds)
          g = -0.1007_dp - 0.358_dp*f + 0.0261_dp*f*f
          k0 = exp(g)
          z = log(stokes / (k0+1.e-15_dp))
          H = 0.1465_dp + 1.302_dp*z - 0.607_dp*z*z + 0.293_dp*z*z*z
          yc0 = 2.0_dp / pi * atan(h)
          ef_rw = (yc0+p)*(yc0+p) / ((1.+p)*(1.+p))
        endif
        t_efrw(i,j) = max(0.0_dp, min(ef_rw, 0.95_dp))
      enddo
    enddo
  end subroutine compute_efrw


  subroutine compute_efsw()
    !! collision efficiency for snow collecting cloud water from Wang and Ji (2000) except
    !! equate melted snow diameter to their "effective collision cross-section."
    use module_mp_tempo_params, only : wp, sp, dp, &
      nbc, dc, av_s, bv_s, ds, nbs, am_s, bm_s, am_r, obmr, fv_s, &
      t_efsw, d0s, rho_w, pi

    real(dp) :: ds_m, vts, vtc, stokes, reynolds, ef_sw
    real(dp) :: p, yc0, f, g, h, z, k0
    integer :: i, j

    do j = 1, nbc
      vtc = 1.19e4_dp * (1.0e4_dp*dc(j)*dc(j)*0.25_dp)
      do i = 1, nbs
        vts = av_s*ds(i)**bv_s * exp(real(-fv_s*ds(i), kind=dp)) - vtc
        ds_m = (am_s*ds(i)**bm_s / am_r)**obmr
        p = dc(j) / ds_m

        if (p > 0.25_dp .or. ds(i) < d0s .or. dc(j) < 6.e-6_dp .or. vts < 1.e-3_dp) then
          t_efsw(i,j) = 0.0_dp
        else
          stokes = dc(j) * dc(j) * vts * rho_w / (9.*1.718e-5_dp*ds_m)
          reynolds = 9._dp * stokes / (p*p*rho_w)

          f = log(reynolds)
          g = -0.1007_dp - 0.358_dp*f + 0.0261_dp*f*f
          k0 = exp(g)
          z = log(stokes / (k0+1.e-15_dp))
          h = 0.1465_dp + 1.302_dp*z - 0.607_dp*z*z + 0.293_dp*z*z*z
          yc0 = 2.0_dp / pi * atan(h)
          ef_sw = (yc0+p)*(yc0+p) / ((1.+p)*(1.+p))

          t_efsw(i,j) = max(0.0_dp, min(ef_sw, 0.95_dp))
        endif
      enddo
    enddo
  end subroutine compute_efsw


 subroutine compute_drop_evap()
    !! calculates droplet evaporation data
    use module_mp_tempo_params, only: wp, sp, dp, &
      nbc, am_r, dc, bm_r, nu_c_scale, nu_c_max, nu_c_min, &
      t_nc, ntb_c, cce, ccg, ocg1, r_c, obmr, dtc, &
      tpc_wev, tnc_wev

    integer :: i, j, k, n
    real(dp), dimension(nbc) :: n_c, massc
    real(dp) :: summ, summ2, lamc, n0_c
    integer :: nu_c

    do n = 1, nbc
      massc(n) = am_r*dc(n)**bm_r
    enddo

    do k = 1, nbc
      nu_c = get_nuc(t_nc(k))
      do j = 1, ntb_c
        lamc = (t_nc(k)*am_r* ccg(2,nu_c)*ocg1(nu_c) / r_c(j))**obmr
        n0_c = t_nc(k)*ocg1(nu_c) * lamc**cce(1,nu_c)
        do i = 1, nbc
          n_c(i) = n0_c* dc(i)**nu_c*exp(-lamc*dc(i))*dtc(i)
          summ = 0._dp
          summ2 = 0._dp
          do n = 1, i
            summ = summ + massc(n)*n_c(n)
            summ2 = summ2 + n_c(n)
          enddo
          tpc_wev(i,j,k) = summ
          tnc_wev(i,j,k) = summ2
        enddo
      enddo
    enddo
  end subroutine compute_drop_evap


  subroutine qr_acr_qs(local_start, local_end, &
    ltcs_racs1, ltmr_racs1, ltcs_racs2, ltmr_racs2, ltcr_sacr1, ltms_sacr1, &
    ltcr_sacr2, ltms_sacr2, ltnr_racs1, ltnr_racs2, ltnr_sacr1, ltnr_sacr2)
    !! calculate rain collecting snow (and inverse)
    use module_mp_tempo_params, only : table_dp, &
      nbr, nbs, dr, av_s, bv_s, ds, fv_s, &
      ntb_r, ntb_r1, n0r_exp, am_r, cre, crg, ore1, r_r, &
      org1, org2, obmr, mu_r, dtr, ntb_t, ntb_s, r_s, &
      sa, sb, tc, bm_s, mu_s, lam0, lam1, kap0, kap1, dts, &
      bm_r, am_s, pi, ef_rs, table_dp

    integer, intent(in) :: local_start, local_end
    real(table_dp), intent(out), dimension(:,:,:,:) :: &
      ltcs_racs1, ltmr_racs1, ltcs_racs2, ltmr_racs2, &
      ltcr_sacr1, ltms_sacr1, ltcr_sacr2, ltms_sacr2, &
      ltnr_racs1, ltnr_racs2, ltnr_sacr1, ltnr_sacr2

    integer :: i, j, k, m, n, n2
    real(dp), dimension(nbr) :: vr, d1, n_r
    real(dp), dimension(nbs) :: vs, n_s
    real(dp) :: m0, m2, m3, mrat, om3
    real(dp) :: n0_r, lam_exp, lamr, slam1, slam2
    real(dp) :: dvs, dvr, masss, massr
    real(dp) :: t1, t2, t3, t4, z1, z2, z3, z4
    real(dp) :: y1, y2, y3, y4

    do n2 = 1, nbr
      vr(n2) = -0.1021_dp + 4.932e3_dp*dr(n2) - 0.9551e6_dp*dr(n2)*dr(n2) &
        + 0.07934e9_dp*dr(n2)*dr(n2)*dr(n2) &
        - 0.002362e12_dp*dr(n2)*dr(n2)*dr(n2)*dr(n2)
      d1(n2) = (vr(n2)/av_s)**(1._dp/bv_s)
    enddo
    do n = 1, nbs
      vs(n) = 1.5_dp*av_s*ds(n)**bv_s * exp(real(-fv_s*ds(n), kind=dp))
    enddo

    do m = local_start, local_end
      do k = 1, ntb_r1
        lam_exp = (n0r_exp(k)*am_r*crg(1)/r_r(m))**ore1
        lamr = lam_exp * (crg(3)*org2*org1)**obmr
        n0_r = n0r_exp(k)/(crg(2)*lam_exp) * lamr**cre(2)
        do n2 = 1, nbr
          n_r(n2) = n0_r*dr(n2)**mu_r * exp(-lamr*dr(n2))*dtr(n2)
        enddo

        do j = 1, ntb_t
          do i = 1, ntb_s
            call snow_moments(rs=r_s(i), tc=tc(j), smob=m2, smoc=m3)

            om3 = 1._wp/m3
            mrat = m2*(m2*om3)*(m2*om3)*(m2*om3)
            m0   = (m2*om3)**mu_s
            slam1 = m2 * om3 * lam0
            slam2 = m2 * om3 * lam1

            do n = 1, nbs
              n_s(n) = mrat*(kap0*exp(-slam1*ds(n)) &
                + kap1*m0*ds(n)**mu_s * exp(-slam2*ds(n)))*dts(n)
            enddo

            t1 = 0._dp
            t2 = 0._dp
            t3 = 0._dp
            t4 = 0._dp
            z1 = 0._dp
            z2 = 0._dp
            z3 = 0._dp
            z4 = 0._dp
            y1 = 0._dp
            y2 = 0._dp
            y3 = 0._dp
            y4 = 0._dp
            do n2 = 1, nbr
              massr = am_r * dr(n2)**bm_r
              do n = 1, nbs
                masss = am_s * ds(n)**bm_s

                dvs = 0.5_dp*((vr(n2) - vs(n)) + abs(vr(n2)-vs(n)))
                dvr = 0.5_dp*((vs(n) - vr(n2)) + abs(vs(n)-vr(n2)))
                if (massr > 1.5*masss) then
                  t1 = t1+ pi*.25_dp*ef_rs*(ds(n)+dr(n2))*(ds(n)+dr(n2)) &
                      *dvs*masss * n_s(n)* n_r(n2)
                  z1 = z1+ pi*.25_dp*ef_rs*(ds(n)+dr(n2))*(ds(n)+dr(n2)) &
                      *dvs*massr * n_s(n)* n_r(n2)
                  y1 = y1+ pi*.25_dp*ef_rs*(ds(n)+dr(n2))*(ds(n)+dr(n2)) &
                      *dvs       * n_s(n)* n_r(n2)
                else
                  t3 = t3+ pi*.25_dp*ef_rs*(ds(n)+dr(n2))*(ds(n)+dr(n2)) &
                      *dvs*masss * n_s(n)* n_r(n2)
                  z3 = z3+ pi*.25_dp*ef_rs*(ds(n)+dr(n2))*(ds(n)+dr(n2)) &
                      *dvs*massr * n_s(n)* n_r(n2)
                  y3 = y3+ pi*.25_dp*ef_rs*(ds(n)+dr(n2))*(ds(n)+dr(n2)) &
                      *dvs       * n_s(n)* n_r(n2)
                endif

                if (massr > 1.5_dp*masss) then
                  t2 = t2+ pi*.25_dp*ef_rs*(ds(n)+dr(n2))*(ds(n)+dr(n2)) &
                      *dvr*massr * n_s(n)* n_r(n2)
                  y2 = y2+ pi*.25_dp*ef_rs*(ds(n)+dr(n2))*(ds(n)+dr(n2)) &
                      *dvr       * n_s(n)* n_r(n2)
                  z2 = z2+ pi*.25_dp*ef_rs*(ds(n)+dr(n2))*(ds(n)+dr(n2)) &
                      *dvr*masss * n_s(n)* n_r(n2)
                else
                  t4 = t4+ pi*.25_dp*ef_rs*(ds(n)+dr(n2))*(ds(n)+dr(n2)) &
                      *dvr*massr * n_s(n)* n_r(n2)
                  y4 = y4+ pi*.25_dp*ef_rs*(ds(n)+dr(n2))*(ds(n)+dr(n2)) &
                      *dvr       * n_s(n)* n_r(n2)
                  z4 = z4+ pi*.25_dp*ef_rs*(ds(n)+dr(n2))*(ds(n)+dr(n2)) &
                      *dvr*masss * n_s(n)* n_r(n2)
                endif
              enddo
            enddo
            ltcs_racs1(i,j,k,m-local_start+1) = t1
            ltmr_racs1(i,j,k,m-local_start+1) = min(z1, real(r_r(m), kind=dp))
            ltcs_racs2(i,j,k,m-local_start+1) = t3
            ltmr_racs2(i,j,k,m-local_start+1) = z3
            ltcr_sacr1(i,j,k,m-local_start+1) = t2
            ltms_sacr1(i,j,k,m-local_start+1) = z2
            ltcr_sacr2(i,j,k,m-local_start+1) = t4
            ltms_sacr2(i,j,k,m-local_start+1) = z4
            ltnr_racs1(i,j,k,m-local_start+1) = y1
            ltnr_racs2(i,j,k,m-local_start+1) = y3
            ltnr_sacr1(i,j,k,m-local_start+1) = y2
            ltnr_sacr2(i,j,k,m-local_start+1) = y4
          enddo
        enddo
      enddo
    enddo
  end subroutine qr_acr_qs


  subroutine qr_acr_qg(local_start, local_end, &
    ltcg_racg, ltmr_racg, ltcr_gacr, ltnr_racg, ltnr_gacr)
    !! Rain collecting graupel (and inverse).  Explicit CE integration.
    use module_mp_tempo_params, only : table_dp, nrhg, nbg, nbr, dr, &
      av_g, dg, bv_g, ntb_r, ntb_r1, &
      n0r_exp, am_r, cre, crg, r_r, ore1, org1, org2, &
      obmr, mu_r, dtr, ntb_g, ntb_g1, n0g_exp, am_g, cge, cgg, &
      r_g, oge1, ogg1, ogg2, obmg, mu_g, dtg, bm_r, bm_g, pi, ef_rg

    integer, intent(in) :: local_start, local_end
    real(table_dp), intent(out), dimension(:,:,:,:,:) :: &
      ltcg_racg, ltmr_racg, ltcr_gacr, ltnr_racg, ltnr_gacr

    integer :: i, j, k, m, n, n2, n3
    real(dp), dimension(nbg) :: n_g
    real(dp), dimension(nbg, nrhg) :: vg
    real(dp), dimension(nbr):: vr, n_r
    real(dp) :: n0_r, n0_g, lam_exp, lamg, lamr
    real(dp) :: massg, massr, dvg, dvr, t1, t2, z1, z2, y1, y2

    do n2 = 1, nbr
      vr(n2) = -0.1021_dp + 4.932e3_dp*dr(n2) - 0.9551e6_dp*dr(n2)*dr(n2) &
        + 0.07934e9_dp*dr(n2)*dr(n2)*dr(n2) &
        - 0.002362e12_dp*dr(n2)*dr(n2)*dr(n2)*dr(n2)
    enddo
    do n3 = 1, nrhg
      do n = 1, nbg
        vg(n,n3) = av_g(n3)*dg(n)**bv_g(n3)
      enddo
    enddo

    do m = local_start, local_end
      do k = 1, ntb_r1

        lam_exp = (n0r_exp(k)*am_r*crg(1)/r_r(m))**ore1
        lamr = lam_exp * (crg(3)*org2*org1)**obmr
        n0_r = n0r_exp(k)/(crg(2)*lam_exp) * lamr**cre(2)
        do n2 = 1, nbr
          n_r(n2) = n0_r*dr(n2)**mu_r *exp(-lamr*dr(n2))*dtr(n2)
        enddo

        do n3 = 1, nrhg
          do j = 1, ntb_g
            do i = 1, ntb_g1
              lam_exp = (n0g_exp(i)*am_g(n3)*cgg(1,1)/r_g(j))**oge1
              lamg = lam_exp * (cgg(3,1)*ogg2*ogg1)**obmg
              n0_g = n0g_exp(i)/(cgg(2,1)*lam_exp) * lamg**cge(2,1)
              do n = 1, nbg
                n_g(n) = n0_g*dg(n)**mu_g * exp(-lamg*dg(n))*dtg(n)
              enddo

              t1 = 0._dp
              t2 = 0._dp
              z1 = 0._dp
              z2 = 0._dp
              y1 = 0._dp
              y2 = 0._dp
              do n2 = 1, nbr
                massr = am_r * dr(n2)**bm_r
                do n = 1, nbg
                  massg = am_g(n3) * dg(n)**bm_g

                  dvg = 0.5_dp*((vr(n2) - vg(n,n3)) + abs(vr(n2)-vg(n,n3)))
                  dvr = 0.5_dp*((vg(n,n3) - vr(n2)) + abs(vg(n,n3)-vr(n2)))

                  t1 = t1+ pi*.25_wp*ef_rg*(dg(n)+dr(n2))*(dg(n)+dr(n2)) &
                      *dvg*massg * n_g(n)* n_r(n2)
                  z1 = z1+ pi*.25_wp*ef_rg*(dg(n)+dr(n2))*(dg(n)+dr(n2)) &
                      *dvg*massr * n_g(n)* n_r(n2)
                  y1 = y1+ pi*.25_wp*ef_rg*(dg(n)+dr(n2))*(dg(n)+dr(n2)) &
                      *dvg       * n_g(n)* n_r(n2)

                  t2 = t2+ pi*.25_wp*ef_rg*(dg(n)+dr(n2))*(dg(n)+dr(n2)) &
                      *dvr*massr * n_g(n)* n_r(n2)
                  y2 = y2+ pi*.25_wp*ef_rg*(dg(n)+dr(n2))*(dg(n)+dr(n2)) &
                      *dvr       * n_g(n)* n_r(n2)
                  z2 = z2+ pi*.25_wp*ef_rg*(dg(n)+dr(n2))*(dg(n)+dr(n2)) &
                      *dvr*massg * n_g(n)* n_r(n2)
                enddo
              enddo
              ltcg_racg(i,j,n3,k,m-local_start+1) = t1
              ltmr_racg(i,j,n3,k,m-local_start+1) = min(z1, real(r_r(m), kind=dp))
              ltcr_gacr(i,j,n3,k,m-local_start+1) = t2
              ltnr_racg(i,j,n3,k,m-local_start+1) = y1
              ltnr_gacr(i,j,n3,k,m-local_start+1) = y2
            enddo
          enddo
        enddo
      enddo
    enddo 
  end subroutine qr_acr_qg


  subroutine freezewater()
    !! freeze water from  Bigg (1953), calculating the probability of drops of a particular volume freezing
    use module_mp_tempo_params, only : nbc, nbr, rho_w, &
      am_r, dr, dtr, bm_r, dc, dtc, ntb_in, ntb_r1, ntb_r, nt_in, &
      n0r_exp, cre, crg, r_r, ore1, org2, org1, obmr, mu_r, &
      nu_c_scale, xm0g, t_nc, cce, ccg, ocg1, r_c, ntb_c, &
      tpi_qrfz, tni_qrfz, tpg_qrfz, tnr_qrfz, tpi_qcfz, tni_qcfz ! data arrays

    integer :: i, j, k, m, n, n2
    real(dp) :: n_r, n_c
    real(dp), dimension(nbr):: massr
    real(dp), dimension(nbc):: massc
    real(dp) :: sum1, sum2, sumn1, sumn2, &
        prob, vol, texp, orho_w, &
        lam_exp, lamr, N0_r, lamc, n0_c
    integer :: nu_c
    real(wp) :: t_adjust

    orho_w = 1._wp/rho_w

    do n2 = 1, nbr
      massr(n2) = am_r*dr(n2)**bm_r
    enddo
    do n = 1, nbc
      massc(n) = am_r*dc(n)**bm_r
    enddo

    ! the smallest drops become cloud ice, otherwise graupel
    do m = 1, ntb_in
      t_adjust = max(-3.0_wp, min(3.0_wp - log10(nt_in(m)), 3.0_wp))
      do k = 1, 45
        texp = exp(real(k, kind=dp) - real(t_adjust, kind=dp)) - 1.0_dp
        do j = 1, ntb_r1
          do i = 1, ntb_r
            lam_exp = (n0r_exp(j)*am_r*crg(1)/r_r(i))**ore1
            lamr = lam_exp * (crg(3)*org2*org1)**obmr
            n0_r = n0r_exp(j)/(crg(2)*lam_exp) * lamr**cre(2)
            sum1 = 0._dp
            sum2 = 0._dp
            sumn1 = 0._dp
            sumn2 = 0._dp
            do n2 = nbr, 1, -1
              n_r = n0_r*dr(n2)**mu_r*exp(-lamr*dr(n2))*dtr(n2)
              vol = massr(n2)*orho_w
              prob = max(0.0_dp, 1.0_dp - exp(-120.0_dp*vol*5.2d-4 * texp))
              !!@note
              !! Graupel tuning parameter: In this table, if the frozen raindrop mass is >= xm0g,
              !! initialization to graupel occurs. 
              !!@endnote
              if (massr(n2) < xm0g) then
                sumn1 = sumn1 + prob*n_r
                sum1 = sum1 + prob*n_r*massr(n2)
              else
                sumn2 = sumn2 + prob*n_r
                sum2 = sum2 + prob*n_r*massr(n2)
              endif
              if ((sum1+sum2) >= r_r(i)) exit
            enddo
            tpi_qrfz(i,j,k,m) = sum1
            tni_qrfz(i,j,k,m) = sumn1
            tpg_qrfz(i,j,k,m) = sum2
            tnr_qrfz(i,j,k,m) = sumn2
          enddo
        enddo

        do j = 1, nbc
          nu_c = min(15, nint(nu_c_scale/t_nc(j)) + 2)
          do i = 1, ntb_c
            lamc = (t_nc(j)*am_r* ccg(2,nu_c) * ocg1(nu_c) / r_c(i))**obmr
            n0_c = t_nc(j)*ocg1(nu_c) * lamc**cce(1,nu_c)
            sum1 = 0._dp
            sumn2 = 0._dp
            do n = nbc, 1, -1
              vol = massc(n)*orho_w
              prob = max(0.0_dp, 1.0_dp - exp(-120.0_dp*vol*5.2e-4_dp * texp))
              n_c = n0_c*dc(n)**nu_c*exp(-lamc*dc(n))*dtc(n)
              sumn2 = min(t_nc(j), sumn2 + prob*n_c)
              sum1 = sum1 + prob*n_c*massc(n)
              if (sum1 >= r_c(i)) exit
            enddo
            tpi_qcfz(i,j,k,m) = sum1
            tni_qcfz(i,j,k,m) = sumn2
          enddo
        enddo
      enddo
    enddo
  end subroutine freezewater


  subroutine qi_aut_qs()
    !! fills data array that contain values of cloud ice conversion to snow
    !! and depositional growth by binning cloud ice distributions and 
    !! determining both the size bins > d0s (that are converted to snow) and the 
    !! depositional growth up to d0s (for cloud ice) and > d0s for snow 
    !! following Harrington et al. (1995)
    use module_mp_tempo_params, only : wp, sp, dp, &
      nbi, ntb_i, ntb_i1, &
      am_i, cie, cig, oig1, nt_i, r_i, obmi, bm_i, mu_i, d0s, d0i, &
      tpi_ide, tps_iaus, tni_iaus, di, dti

    integer :: i, j, n2
    real(dp), dimension(nbi) :: n_i
    real(dp) :: n0_i, lami, di_mean, t1, t2
    real(wp) :: xlimit_intg

    do j = 1, ntb_i1
      do i = 1, ntb_i
        lami = (am_i*cig(2)*oig1*nt_i(j)/r_i(i))**obmi
        di_mean = (bm_i + mu_i + 1.) / lami
        n0_i = nt_i(j)*oig1 * lami**cie(1)
        t1 = 0.
        t2 = 0.
        if (real(di_mean, kind=wp) > 5.*d0s) then
          t1 = r_i(i)
          t2 = nt_i(j)
          tpi_ide(i,j) = 0.
        elseif (real(di_mean, kind=wp) < d0i) then
          t1 = 0.
          t2 = 0.
          tpi_ide(i,j) = 1.
        else
          xlimit_intg = lami*d0s
          tpi_ide(i,j) = real(calc_gamma_p(mu_i+2.0, xlimit_intg), kind=dp)
          do n2 = 1, nbi
            n_i(n2) = n0_i*di(n2)**mu_i * exp(-lami*di(n2))*dti(n2)
            if (di(n2) >= d0s) then
              t1 = t1 + n_i(n2) * am_i*di(n2)**bm_i
              t2 = t2 + n_i(n2)
            endif
          enddo
        endif
        tps_iaus(i,j) = t1
        tni_iaus(i,j) = t2
      enddo
    enddo
  end subroutine qi_aut_qs

end module module_mp_tempo_init