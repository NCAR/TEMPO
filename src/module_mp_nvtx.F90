module module_mp_nvtx
  !! NVTX ranges for nsys (same pattern as MPAS nvtxTimer.F + nvtxMakefile: use nvtx, nvtxStartRange / nvtxEndRange,
  !! link with -lnvhpcwrapnvtx when built with -DTEMPO_NVTX_RANGES).
#ifdef TEMPO_NVTX_RANGES
  use nvtx
#endif
  implicit none
  private
  public :: nvtx_range_push, nvtx_range_pop

contains

  subroutine nvtx_range_push(name)
    character(len=*), intent(in) :: name
#ifdef TEMPO_NVTX_RANGES
    call nvtxStartRange(trim(name))
#endif
  end subroutine nvtx_range_push

  subroutine nvtx_range_pop()
#ifdef TEMPO_NVTX_RANGES
    call nvtxEndRange()
#endif
  end subroutine nvtx_range_pop

end module module_mp_nvtx
