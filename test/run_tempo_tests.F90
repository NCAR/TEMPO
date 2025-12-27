program run_tempo_tests
  !! This program runs TEMPO tests
  use tests, only : test_tempo_init, test_tempo_driver, &
    test_rain_sedimentation
  implicit none

  call test_tempo_init()
  !call test_tempo_driver()
  call test_rain_sedimentation(dt=1.)
  call test_rain_sedimentation(dt=10.)
  call test_rain_sedimentation(dt=20.)
  call test_rain_sedimentation(dt=60.)
  call test_rain_sedimentation(dt=120.)
  call test_rain_sedimentation(dt=300.)
  call test_rain_sedimentation(dt=600.)
end program run_tempo_tests
