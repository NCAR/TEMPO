program run_tempo_tests
  !! This program runs TEMPO tests
  use tests, only : test_tempo_init !, test_tempo_driver
  implicit none

  call test_tempo_init()
  ! call test_tempo_driver()
  
end program run_tempo_tests
