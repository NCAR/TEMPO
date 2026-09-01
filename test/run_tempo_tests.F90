program run_tempo_tests
  !! runs tempo tests
  use tests, only : test_tempo_init, test_graupel_sedimentation, &
    test_snow_sedimentation, test_cloud_number_aerosolaware, &
    test_cloud_number_non_aerosolaware, test_cloud_number_ml, &
    test_ml_cloud_effective_radius, set_ncells, set_nout_values, set_stride
  use module_mp_nvtx, only : nvtx_range_push, nvtx_range_pop

  implicit none

  real, dimension(7) :: sedi_tests = &
    [1., 10., 20., 60., 120., 300., 600.]
  real, dimension(1) :: sedi_tests = &
    [20.]
  integer :: t, ncells, nargs, i, nout, stride
  character(len=32) :: arg
  character(len=160) :: nvtx_lbl
  logical :: have_nout, have_stride

  ! Parse command line: -c <n> sets number of cells (nCells); ide, ime, ite, jde, jme, jte = ncells
  ! -s <n> sets the hybrid horizontal block size (stride); 1=per-column (CPU/column form),
  !        ncells=full plane (GPU/plane form), <=0 or unset => default (full plane).
  ! -o <n> sets nOutValues for graupel sedimentation profile output (equally spaced (i,j) on diagonal)
  ncells = 1
  have_nout = .false.
  nout = 1
  have_stride = .false.
  stride = -1
  nargs = command_argument_count()
  i = 1
  do while (i <= nargs)
    call get_command_argument(i, arg)
    if (arg == '-c' .and. i < nargs) then
      call get_command_argument(i + 1, arg)
      read(arg, *, err=1) ncells
      ncells = max(1, ncells)
      i = i + 2
      cycle
    1 ncells = 1
      i = i + 2
      cycle
    endif
    if (arg == '-s' .and. i < nargs) then
      call get_command_argument(i + 1, arg)
      read(arg, *, err=3) stride
      have_stride = .true.
      i = i + 2
      cycle
    3 stride = -1
      have_stride = .false.
      i = i + 2
      cycle
    endif
    if (arg == '-o' .and. i < nargs) then
      call get_command_argument(i + 1, arg)
      read(arg, *, err=2) nout
      nout = max(1, nout)
      have_nout = .true.
      i = i + 2
      cycle
    2 nout = 1
      have_nout = .true.
      i = i + 2
      cycle
    endif
    i = i + 1
  enddo
  call set_ncells(ncells)
  call set_stride(stride)
  if (have_nout) call set_nout_values(nout)
  write(*, '(A,I0)') 'Running TEMPO tests with nCells = ', ncells
  if (have_stride) then
    write(*, '(A,I0)') 'Hybrid stride (horizontal block size) = ', stride
  else
    write(*, '(A)') 'Hybrid stride (horizontal block size) = default (full plane = nCells)'
  endif
  if (have_nout) then
    write(*, '(A,I0)') 'Graupel sedimentation profile output pairs (nOutValues) = ', nout
  endif

  ! tempo init
  call nvtx_range_push('test_tempo_init')
  call test_tempo_init()
  call nvtx_range_pop()

  ! ml cloud effective radius
  call nvtx_range_push('test_ml_cloud_effective_radius')
  call test_ml_cloud_effective_radius(dt=20.)
  call nvtx_range_pop()

  ! test cloud number concentration
  call nvtx_range_push('test_cloud_number_aerosolaware')
  call test_cloud_number_aerosolaware(dt=20.)
  call nvtx_range_pop()

  call nvtx_range_push('test_cloud_number_non_aerosolaware')
  call test_cloud_number_non_aerosolaware(dt=20.)
  call nvtx_range_pop()

  call nvtx_range_push('test_cloud_number_ml')
  call test_cloud_number_ml(dt=20.)
  call nvtx_range_pop()

  ! graupel sedimentation
  do t = 1, size(sedi_tests)
    write(nvtx_lbl, '(A,F0.2,A)') 'test_graupel_sedimentation dt=', sedi_tests(t), ' semi=F'
    call nvtx_range_push(trim(nvtx_lbl))
    call test_graupel_sedimentation(dt=sedi_tests(t), semi_sedi=.false.)
    call nvtx_range_pop()
  enddo
  do t = 1, size(sedi_tests)
    write(nvtx_lbl, '(A,F0.2,A)') 'test_graupel_sedimentation dt=', sedi_tests(t), ' semi=T'
    call nvtx_range_push(trim(nvtx_lbl))
    call test_graupel_sedimentation(dt=sedi_tests(t), semi_sedi=.true.)
    call nvtx_range_pop()
  enddo

  ! snow sedimentation
  do t = 1, size(sedi_tests)
    write(nvtx_lbl, '(A,F0.2,A)') 'test_snow_sedimentation dt=', sedi_tests(t)
    call nvtx_range_push(trim(nvtx_lbl))
    call test_snow_sedimentation(dt=sedi_tests(t))
    call nvtx_range_pop()
  enddo

end program run_tempo_tests
