program run_tempo_tests
  !! This program runs TEMPO tests
  use tests, only : test_tempo_init, test_graupel_sedimentation, &
    test_snow_sedimentation
  implicit none

  real, dimension(7) :: sedi_tests = &
    [1., 10., 20., 60., 120., 300., 600.]

  integer :: t

  call test_tempo_init()
  
  ! graupel sedimentation
  do t = 1, size(sedi_tests)
    call test_graupel_sedimentation(dt=sedi_tests(t), semi_sedi=.false.)
  enddo
  do t = 1, size(sedi_tests)
    call test_graupel_sedimentation(dt=sedi_tests(t), semi_sedi=.true.)
  enddo

  ! snow sedimentation
  do t = 1, size(sedi_tests)
    call test_snow_sedimentation(dt=sedi_tests(t))
  enddo
  ! test_tempo_driver()

end program run_tempo_tests
