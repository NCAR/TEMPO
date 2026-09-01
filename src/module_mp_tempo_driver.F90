module module_mp_tempo_driver
  !! tempo driver: 3D tile workspace passed to tempo_main (horizontal loop is inside tempo_main)
  !! also allocates and fills diagnostic arrays
  use module_mp_tempo_cfgs, only : ty_tempo_cfgs, ty_tempo_table_cfgs
  use module_mp_tempo_params, only : wp, sp, dp
  use module_mp_nvtx, only : nvtx_range_push, nvtx_range_pop
  use module_mp_tempo_main, only : tempo_main, ty_tempo_main_diags
  use module_mp_tempo_utils, only : compute_efrw, compute_efsw, compute_drop_evap, qi_aut_qs
  use module_mp_tempo_ml, only : ty_tempo_ml_data, nc_ml_nodes, nc_ml_input, nc_ml_output, &
    nc_ml_trans_mean, nc_ml_trans_var, nc_ml_w00, nc_ml_w01, nc_ml_b00, nc_ml_b01, save_or_read_ml_data

  implicit none
  private

  public :: tempo_init, tempo_run, ty_tempo_driver_diags, tempo_run_enter_data, tempo_run_exit_data, &
    tempo_run_sync_host_diags, tempo_run_sync_host_fields, tempo_aerosol_surface_emissions
  
  type(ty_tempo_table_cfgs) :: tempo_table_cfgs

  type :: ty_tempo_driver_diags
    real(wp), dimension(:,:), allocatable :: rain_precip
    real(wp), dimension(:,:), allocatable :: ice_liquid_equiv_precip
    real(wp), dimension(:,:), allocatable :: snow_liquid_equiv_precip
    real(wp), dimension(:,:), allocatable :: graupel_liquid_equiv_precip
    real(wp), dimension(:,:), allocatable :: frozen_fraction
    real(wp), dimension(:,:), allocatable :: frz_rain_precip
    real(wp), dimension(:,:), allocatable :: max_hail_diameter_sfc
    real(wp), dimension(:,:), allocatable :: max_hail_diameter_column
    real(wp), dimension(:,:,:), allocatable :: refl10cm
    real(wp), dimension(:,:,:), allocatable :: re_cloud
    real(wp), dimension(:,:,:), allocatable :: re_ice
    real(wp), dimension(:,:,:), allocatable :: re_snow
    real(wp), dimension(:,:,:), allocatable :: rain_med_vol_diam
    real(wp), dimension(:,:,:), allocatable :: graupel_med_vol_diam
    real(wp), dimension(:,:,:), allocatable :: cloud_number_mixing_ratio
  end type

  contains

!> initialize tempo microphysics
!! \section arg_table_tempo_init Argument Table
!! \htmlinclude tempo_init.html
!!
  subroutine tempo_init(aerosolaware_flag, hailaware_flag, semi_sedi_flag, cloud_condensation_flag, &
    refl10cm_from_melting_flag, ml_for_bl_nc_flag, ml_for_nc_flag, force_init_flag, tempo_cfgs)
    !! initialize tempo microphysics
    use module_mp_tempo_params, only : get_version, tempo_version, t_efrw, &
      initialize_graupel_vars, initialize_parameters, initialize_bins_for_tables, &
      initialize_array_efrw, initialize_array_efsw, initialize_arrays_drop_evap, &
      initialize_arrays_ccn, initialize_arrays_qi_aut_qs, &
      initialize_arrays_qr_acr_qs, initialize_arrays_qr_acr_qg, initialize_arrays_freezewater, &
      initialize_bins_for_hail_size, initialize_bins_for_radar

    logical, intent(in), optional :: aerosolaware_flag, hailaware_flag, refl10cm_from_melting_flag, &
      ml_for_bl_nc_flag, ml_for_nc_flag, force_init_flag, semi_sedi_flag, cloud_condensation_flag
    type(ty_tempo_cfgs), intent(out) :: tempo_cfgs

    character(len=100) :: table_filename
    integer :: table_size
    logical :: initialize_mp_vars, force_init

    call nvtx_range_push('tempo_init')

    ! get tempo version from readme file
    call get_version(tempo_version) 

    ! check an allocatable array (t_efrw) to see if initialization can be skipped
    ! but allow for force initialization useful for testing
    force_init = .false.
    if (present(force_init_flag)) force_init = force_init_flag

    initialize_mp_vars = .true.
    if (allocated(t_efrw)) initialize_mp_vars = .false.
    if (force_init) initialize_mp_vars = .true.

    if (initialize_mp_vars) then
      call nvtx_range_push('tempo_init_initialize_mp_vars')
      if (present(aerosolaware_flag)) tempo_cfgs%aerosolaware_flag = aerosolaware_flag
      if (present(hailaware_flag)) tempo_cfgs%hailaware_flag = hailaware_flag
      if (present(ml_for_bl_nc_flag)) tempo_cfgs%ml_for_bl_nc_flag = ml_for_bl_nc_flag
      if (present(ml_for_nc_flag)) tempo_cfgs%ml_for_nc_flag = ml_for_nc_flag
      if (present(semi_sedi_flag)) tempo_cfgs%semi_sedi_flag = semi_sedi_flag
      if (present(cloud_condensation_flag)) tempo_cfgs%cloud_condensation_flag = cloud_condensation_flag
      if (present(refl10cm_from_melting_flag)) tempo_cfgs%refl10cm_from_melting_flag = refl10cm_from_melting_flag

      if (tempo_cfgs%verbose) then
        write(*,'(A)') 'tempo_init() --- TEMPO microphysics configuration options: '
        write(*,'(A,L)') 'tempo_init() --- aerosol aware = ', tempo_cfgs%aerosolaware_flag
        write(*,'(A,L)') 'tempo_init() --- hail aware = ', tempo_cfgs%hailaware_flag
        write(*,'(A,L)') 'tempo_init() --- ML for subgrid cloud number = ', tempo_cfgs%ml_for_bl_nc_flag
        write(*,'(A,L)') 'tempo_init() --- ML for cloud number = ', tempo_cfgs%ml_for_nc_flag
        write(*,'(A,L)') 'tempo_init() --- reflectivity from melting snow/graupel = ', tempo_cfgs%refl10cm_from_melting_flag
        write(*,'(A,L)') 'tempo_init() --- semi-lagrangian sedimentation = ', tempo_cfgs%semi_sedi_flag
      endif 

      ! set graupel variables from hail_aware_flag
      call initialize_graupel_vars(tempo_cfgs%hailaware_flag) 
      if (tempo_cfgs%verbose) then
        write(*,'(A,L)') 'tempo_init() --- initialized graupel variables using hail aware = ', tempo_cfgs%hailaware_flag
      endif 

      ! set parameters that can depend on the host model
      call initialize_parameters() 
      if (tempo_cfgs%verbose) write(*,'(A)') 'tempo_init() --- initialized parameters'
      
      ! creates log-spaced bins of hydrometers for tables
      call initialize_bins_for_tables() 
      if (tempo_cfgs%verbose) write(*,'(A)') 'tempo_init() --- initialized bins for lookup tables'

      ! collision efficiencies between rain/snow and cloud water.
      call initialize_array_efrw()
      call compute_efrw()
      if (tempo_cfgs%verbose) then
        write(*,'(A)') 'tempo_init() --- initialized collision efficiency data for rain collecting cloud water'
      endif 
      call initialize_array_efsw()
      call compute_efsw()
      if (tempo_cfgs%verbose) then
        write(*,'(A)') 'tempo_init() --- initialized collision efficiency data for snow collecting cloud water'
      endif 

      ! drop evaporation
      call initialize_arrays_drop_evap()
      call compute_drop_evap()
      if (tempo_cfgs%verbose) write(*,'(A)') 'tempo_init() --- initialized drop evaporation data'

      ! cloud ice to snow and depositional growth
      call initialize_arrays_qi_aut_qs()
      call qi_aut_qs()

      ! CCN activation table
      table_filename = tempo_table_cfgs%ccn_table_name
      call initialize_arrays_ccn(table_size)
      call read_table_ccn(trim(table_filename), table_size)
      if (tempo_cfgs%verbose) write(*,'(A)') 'tempo_init() --- initialized data for ccn lookup table'

      ! freeze water collection lookup table
      table_filename = tempo_table_cfgs%freezewater_table_name
      call initialize_arrays_freezewater(table_size)
      call read_table_freezewater(trim(table_filename), table_size)
      if (tempo_cfgs%verbose) then
        write(*,'(A)') 'tempo_init() --- initialized data for frozen cloud water and rain lookup table'
      endif 

      ! rain-snow collection lookup table
      table_filename = tempo_table_cfgs%qrqs_table_name
      call initialize_arrays_qr_acr_qs(table_size)
      call read_table_qr_acr_qs(trim(table_filename), table_size)
      if (tempo_cfgs%verbose) then
        write(*,'(A)') 'tempo_init() --- initialized data for rain-snow collection lookup table'
      endif 

      ! rain-graupel collection lookup table
      table_filename = tempo_table_cfgs%qrqg_table_name
      call initialize_arrays_qr_acr_qg(table_size)
      call read_table_qr_acr_qg(trim(table_filename), table_size)
      if (tempo_cfgs%verbose) then
        write(*,'(A)') 'tempo_init() --- initialized data for rain-graupel collection lookup table'
      endif 

      ! bins used for optional refl10cm calculation with melting
      if (tempo_cfgs%refl10cm_from_melting_flag) then
        call initialize_bins_for_radar()
        if (tempo_cfgs%verbose) then
          write(*,'(A,L)') 'tempo_init() ---  flag to calcuate reflectivity with contributions from melting snow and graupel = ', &
            tempo_cfgs%refl10cm_from_melting_flag
          write(*,'(A)') 'tempo_init() --- initialized bins for reflectivity calcuation with meting snow and graupel'
        endif 
      endif

      ! bins used for optional hail size calculation
      if (tempo_cfgs%max_hail_diameter_flag) then
        call initialize_bins_for_hail_size()
        if (tempo_cfgs%verbose) then
          write(*,'(A,L)') 'tempo_init() ---  flag to calculate max hail diameter = ', &
            tempo_cfgs%max_hail_diameter_flag
          write(*,'(A)') 'tempo_init() --- initialized bins for hail size calculation'
        endif
      endif

      ! data for machine learning
      if(tempo_cfgs%ml_for_bl_nc_flag .or. tempo_cfgs%ml_for_nc_flag) then
        call init_ml_data()
        if (tempo_cfgs%verbose) write(*,'(A)') 'tempo_init() --- initialized data for cloud number machine learning'
      endif 
      call nvtx_range_pop()
    endif
    call nvtx_range_pop()
  end subroutine tempo_init

  subroutine tempo_run_setup_diags(tempo_cfgs, tempo_diags, kts, kte, its, ite, jts, jte)
    !! Allocate driver diagnostic arrays on the host. Initialization is handled separately.
    type(ty_tempo_cfgs), intent(in) :: tempo_cfgs
    type(ty_tempo_driver_diags), intent(inout) :: tempo_diags
    integer, intent(in) :: kts, kte, its, ite, jts, jte

    if (tempo_cfgs%cloud_number_mixing_ratio_flag) then
      if (.not. allocated(tempo_diags%cloud_number_mixing_ratio)) then
        allocate(tempo_diags%cloud_number_mixing_ratio(kts:kte, its:ite, jts:jte))
      endif
    endif

    if (tempo_cfgs%rain_med_vol_diam_flag) then
      if (.not. allocated(tempo_diags%rain_med_vol_diam)) then
        allocate(tempo_diags%rain_med_vol_diam(kts:kte, its:ite, jts:jte))
      endif
    endif

    if (tempo_cfgs%graupel_med_vol_diam_flag) then
      if (.not. allocated(tempo_diags%graupel_med_vol_diam)) then
        allocate(tempo_diags%graupel_med_vol_diam(kts:kte, its:ite, jts:jte))
      endif
    endif

    if (tempo_cfgs%refl10cm_flag) then
      if (.not. allocated(tempo_diags%refl10cm)) then
        allocate(tempo_diags%refl10cm(kts:kte, its:ite, jts:jte))
      endif
    endif

    if (tempo_cfgs%re_cloud_flag) then
      if (.not. allocated(tempo_diags%re_cloud)) then
        allocate(tempo_diags%re_cloud(kts:kte, its:ite, jts:jte))
      endif
    endif

    if (tempo_cfgs%re_ice_flag) then
      if (.not. allocated(tempo_diags%re_ice)) then
        allocate(tempo_diags%re_ice(kts:kte, its:ite, jts:jte))
      endif
    endif

    if (tempo_cfgs%re_snow_flag) then
      if (.not. allocated(tempo_diags%re_snow)) then
        allocate(tempo_diags%re_snow(kts:kte, its:ite, jts:jte))
      endif
    endif

    if (tempo_cfgs%max_hail_diameter_flag) then
      if (.not. allocated(tempo_diags%max_hail_diameter_sfc)) then
        allocate(tempo_diags%max_hail_diameter_sfc(its:ite, jts:jte))
        allocate(tempo_diags%max_hail_diameter_column(its:ite, jts:jte))
      endif
    endif

    if (.not. allocated(tempo_diags%rain_precip)) then
      allocate(tempo_diags%rain_precip(its:ite, jts:jte))
      allocate(tempo_diags%ice_liquid_equiv_precip(its:ite, jts:jte))
      allocate(tempo_diags%snow_liquid_equiv_precip(its:ite, jts:jte))
      allocate(tempo_diags%graupel_liquid_equiv_precip(its:ite, jts:jte))
      allocate(tempo_diags%frozen_fraction(its:ite, jts:jte))
      allocate(tempo_diags%frz_rain_precip(its:ite, jts:jte))
    endif
  end subroutine tempo_run_setup_diags


  subroutine tempo_run_reset_diags_device(tempo_diags, kts, kte, its, ite, jts, jte)
    !! Initialize/reset driver diagnostic arrays where they are mapped.
    type(ty_tempo_driver_diags), intent(inout) :: tempo_diags
    integer, intent(in) :: kts, kte, its, ite, jts, jte
    integer :: i, j, k

    if (allocated(tempo_diags%rain_precip)) then
      !$acc parallel
      !$acc loop gang vector collapse(2)
      do j = jts, jte
        do i = its, ite
          tempo_diags%rain_precip(i,j) = 0._wp
          tempo_diags%ice_liquid_equiv_precip(i,j) = 0._wp
          tempo_diags%snow_liquid_equiv_precip(i,j) = 0._wp
          tempo_diags%graupel_liquid_equiv_precip(i,j) = 0._wp
          tempo_diags%frozen_fraction(i,j) = 0._wp
          tempo_diags%frz_rain_precip(i,j) = 0._wp
        enddo
      enddo
      !$acc end parallel
    endif

    if (allocated(tempo_diags%max_hail_diameter_sfc)) then
      !$acc parallel
      !$acc loop gang vector collapse(2)
      do j = jts, jte
        do i = its, ite
          tempo_diags%max_hail_diameter_sfc(i,j) = 0._wp
          tempo_diags%max_hail_diameter_column(i,j) = 0._wp
        enddo
      enddo
      !$acc end parallel
    endif

    if (allocated(tempo_diags%cloud_number_mixing_ratio)) then
      !$acc parallel
      !$acc loop gang vector collapse(3)
      do j = jts, jte
        do i = its, ite
          do k = kts, kte
            tempo_diags%cloud_number_mixing_ratio(k,i,j) = 0._wp
          enddo
        enddo
      enddo
      !$acc end parallel
    endif

    if (allocated(tempo_diags%rain_med_vol_diam)) then
      !$acc parallel
      !$acc loop gang vector collapse(3)
      do j = jts, jte
        do i = its, ite
          do k = kts, kte
            tempo_diags%rain_med_vol_diam(k,i,j) = 0._wp
          enddo
        enddo
      enddo
      !$acc end parallel
    endif

    if (allocated(tempo_diags%graupel_med_vol_diam)) then
      !$acc parallel
      !$acc loop gang vector collapse(3)
      do j = jts, jte
        do i = its, ite
          do k = kts, kte
            tempo_diags%graupel_med_vol_diam(k,i,j) = 0._wp
          enddo
        enddo
      enddo
      !$acc end parallel
    endif

    if (allocated(tempo_diags%refl10cm)) then
      !$acc parallel
      !$acc loop gang vector collapse(3)
      do j = jts, jte
        do i = its, ite
          do k = kts, kte
            tempo_diags%refl10cm(k,i,j) = -35._wp
          enddo
        enddo
      enddo
      !$acc end parallel
    endif

    if (allocated(tempo_diags%re_cloud) .and. allocated(tempo_diags%re_ice) .and. allocated(tempo_diags%re_snow)) then
      !$acc parallel
      !$acc loop gang vector collapse(3)
      do j = jts, jte
        do i = its, ite
          do k = kts, kte
            tempo_diags%re_cloud(k,i,j) = 0._wp
            tempo_diags%re_ice(k,i,j) = 0._wp
            tempo_diags%re_snow(k,i,j) = 0._wp
          enddo
        enddo
      enddo
      !$acc end parallel
    endif
  end subroutine tempo_run_reset_diags_device


  subroutine tempo_run_setup_main_diags(tempo_cfgs, tempo_main_diags, kts, kte, its, ite, jts, jte)
    type(ty_tempo_cfgs), intent(in) :: tempo_cfgs
    type(ty_tempo_main_diags), intent(inout) :: tempo_main_diags
    integer, intent(in) :: kts, kte, its, ite, jts, jte

    if (.not. allocated(tempo_main_diags%rain_precip)) then
      allocate(tempo_main_diags%rain_precip(its:ite, jts:jte), source=0._wp)
      allocate(tempo_main_diags%cloud_precip(its:ite, jts:jte), source=0._wp)
      allocate(tempo_main_diags%ice_liquid_equiv_precip(its:ite, jts:jte), source=0._wp)
      allocate(tempo_main_diags%snow_liquid_equiv_precip(its:ite, jts:jte), source=0._wp)
      allocate(tempo_main_diags%graupel_liquid_equiv_precip(its:ite, jts:jte), source=0._wp)
      allocate(tempo_main_diags%frozen_fraction(its:ite, jts:jte), source=0._wp)
      allocate(tempo_main_diags%frz_rain_precip(its:ite, jts:jte), source=0._wp)
    endif
    if (tempo_cfgs%cloud_number_mixing_ratio_flag .and. .not. allocated(tempo_main_diags%cloud_number_mixing_ratio)) then
      allocate(tempo_main_diags%cloud_number_mixing_ratio(kts:kte, its:ite, jts:jte), source=0._wp)
    endif
    if (tempo_cfgs%rain_med_vol_diam_flag .and. .not. allocated(tempo_main_diags%rain_med_vol_diam)) then
      allocate(tempo_main_diags%rain_med_vol_diam(kts:kte, its:ite, jts:jte), source=0._wp)
    endif
    if (tempo_cfgs%graupel_med_vol_diam_flag .and. .not. allocated(tempo_main_diags%graupel_med_vol_diam)) then
      allocate(tempo_main_diags%graupel_med_vol_diam(kts:kte, its:ite, jts:jte), source=0._wp)
    endif
    if (tempo_cfgs%max_hail_diameter_flag .and. .not. allocated(tempo_main_diags%max_hail_diameter)) then
      allocate(tempo_main_diags%max_hail_diameter(kts:kte, its:ite, jts:jte), source=0._wp)
    endif
    if (tempo_cfgs%refl10cm_flag .and. .not. allocated(tempo_main_diags%refl10cm)) then
      allocate(tempo_main_diags%refl10cm(kts:kte, its:ite, jts:jte), source=-35._wp)
    endif
    if ((tempo_cfgs%re_cloud_flag .and. tempo_cfgs%re_ice_flag .and. tempo_cfgs%re_snow_flag) &
        .and. .not. allocated(tempo_main_diags%re_cloud)) then
      allocate(tempo_main_diags%re_cloud(kts:kte, its:ite, jts:jte), source=0._wp)
      allocate(tempo_main_diags%re_ice(kts:kte, its:ite, jts:jte), source=0._wp)
      allocate(tempo_main_diags%re_snow(kts:kte, its:ite, jts:jte), source=0._wp)
    endif
  end subroutine tempo_run_setup_main_diags

  subroutine tempo_run_reset_block_diags(block_diags, kts, kte, bi, bj)
    !! reset the (1:bi,1:bj) region of the per-block diagnostics to their initial values, matching
    !! the source= initialization the full-tile tempo_main_diags receives from a fresh allocation
    !! each timestep. Required because block_diags is reused across the hybrid block loop.
    type(ty_tempo_main_diags), intent(inout) :: block_diags
    integer, intent(in) :: kts, kte, bi, bj
    integer :: i, j, k

    if (allocated(block_diags%rain_precip)) then
      !$acc parallel
      !$acc loop gang vector collapse(2)
      do j = 1, bj
        do i = 1, bi
          block_diags%rain_precip(i,j) = 0._wp
          block_diags%cloud_precip(i,j) = 0._wp
          block_diags%ice_liquid_equiv_precip(i,j) = 0._wp
          block_diags%snow_liquid_equiv_precip(i,j) = 0._wp
          block_diags%graupel_liquid_equiv_precip(i,j) = 0._wp
          block_diags%frozen_fraction(i,j) = 0._wp
          block_diags%frz_rain_precip(i,j) = 0._wp
        enddo
      enddo
      !$acc end parallel
    endif
    if (allocated(block_diags%cloud_number_mixing_ratio)) then
      !$acc parallel
      !$acc loop gang vector collapse(3)
      do j = 1, bj
        do i = 1, bi
          do k = kts, kte
            block_diags%cloud_number_mixing_ratio(k,i,j) = 0._wp
          enddo
        enddo
      enddo
      !$acc end parallel
    endif
    if (allocated(block_diags%rain_med_vol_diam)) then
      !$acc parallel
      !$acc loop gang vector collapse(3)
      do j = 1, bj
        do i = 1, bi
          do k = kts, kte
            block_diags%rain_med_vol_diam(k,i,j) = 0._wp
          enddo
        enddo
      enddo
      !$acc end parallel
    endif
    if (allocated(block_diags%graupel_med_vol_diam)) then
      !$acc parallel
      !$acc loop gang vector collapse(3)
      do j = 1, bj
        do i = 1, bi
          do k = kts, kte
            block_diags%graupel_med_vol_diam(k,i,j) = 0._wp
          enddo
        enddo
      enddo
      !$acc end parallel
    endif
    if (allocated(block_diags%max_hail_diameter)) then
      !$acc parallel
      !$acc loop gang vector collapse(3)
      do j = 1, bj
        do i = 1, bi
          do k = kts, kte
            block_diags%max_hail_diameter(k,i,j) = 0._wp
          enddo
        enddo
      enddo
      !$acc end parallel
    endif
    if (allocated(block_diags%refl10cm)) then
      !$acc parallel
      !$acc loop gang vector collapse(3)
      do j = 1, bj
        do i = 1, bi
          do k = kts, kte
            block_diags%refl10cm(k,i,j) = -35._wp
          enddo
        enddo
      enddo
      !$acc end parallel
    endif
    if (allocated(block_diags%re_cloud)) then
      !$acc parallel
      !$acc loop gang vector collapse(3)
      do j = 1, bj
        do i = 1, bi
          do k = kts, kte
            block_diags%re_cloud(k,i,j) = 0._wp
            block_diags%re_ice(k,i,j) = 0._wp
            block_diags%re_snow(k,i,j) = 0._wp
          enddo
        enddo
      enddo
      !$acc end parallel
    endif
  end subroutine tempo_run_reset_block_diags


  subroutine tempo_run(tempo_cfgs, dt, itimestep, &
    t, th, pii, p, w, dz, &
    qv, qc, qr, qi, qs, qg, ni, nr, &
    nc, nwfa, nifa, ng, qb, &
    qc_bl, qcfrac_bl, &
    qcfrac, qifrac, &
    thten_bl, qvten_bl, qcten_bl, qiten_bl, &
    thten_lwrad, thten_swrad, &
    ids, ide, jds, jde, kds, kde, &
    ims, ime, jms, jme, kms, kme, &
    its, ite, jts, jte, kts, kte, tempo_diags, arguments_on_device, stride)

    type(ty_tempo_cfgs), intent(in) :: tempo_cfgs
    real(wp), intent(in) :: dt !! timestep \([s]]\)
    integer, intent(in) :: itimestep !! integer timestep = integration time / dt
    integer, intent(in) :: ids, ide, jds, jde, kds, kde !! domain locations
    integer, intent(in) :: ims, ime, jms, jme, kms, kme !! memory locations
    integer, intent(in) :: its, ite, jts, jte, kts, kte !! tile locations

    real(wp), dimension(kms:kme, ims:ime, jms:jme), intent(inout), optional :: t !! temperature \([K]\)
    real(wp), dimension(kms:kme, ims:ime, jms:jme), intent(inout), optional :: th !! theta \([K]\)

    real(wp), dimension(kms:kme, ims:ime, jms:jme), intent(in) :: p !! pressure \([Pa]\)
    real(wp), dimension(kms:kme, ims:ime, jms:jme), intent(in) :: w !! vertical velocity \([m\; s^{-1}]\)
    real(wp), dimension(kms:kme, ims:ime, jms:jme), intent(in) :: dz !! vertical grid spacing \([m]\)
    real(wp), dimension(kms:kme, ims:ime, jms:jme), intent(in), optional :: pii !! exner function

    real(wp), dimension(kms:kme, ims:ime, jms:jme), intent(inout) :: qv !! 3D water vapor mixing ratio \([kg\; kg^{-1}]\)
    real(wp), dimension(kms:kme, ims:ime, jms:jme), intent(inout) :: qc !! 3D cloud water mass mixing ratio \([kg\; kg^{-1}]\)
    real(wp), dimension(kms:kme, ims:ime, jms:jme), intent(inout) :: qr !! 3D rain water mass mixing ratio \([kg\; kg^{-1}]\)
    real(wp), dimension(kms:kme, ims:ime, jms:jme), intent(inout) :: qi !! 3D cloud ice mass mixing ratio \([kg\; kg^{-1}]\)
    real(wp), dimension(kms:kme, ims:ime, jms:jme), intent(inout) :: qs !! 3D snow mass mixing ratio \([kg\; kg^{-1}]\)
    real(wp), dimension(kms:kme, ims:ime, jms:jme), intent(inout) :: qg !! 3D graupel mass mixing ratio \([kg\; kg^{-1}]\)
    real(wp), dimension(kms:kme, ims:ime, jms:jme), intent(inout) :: ni !! 3D cloud ice number mixing ratio \([kg^{-1}]\)
    real(wp), dimension(kms:kme, ims:ime, jms:jme), intent(inout) :: nr !! 3D rain water number mixing ratio \([kg^{-1}]\)

    real(wp), dimension(kms:kme, ims:ime, jms:jme), intent(inout), optional :: nc !! 3D cloud water number mixing ratio \([kg^{-1}]\) (aerosol-aware)
    real(wp), dimension(kms:kme, ims:ime, jms:jme), intent(inout), optional :: nwfa !! 3D water-friendly aerosol number mixing ratio \([kg^{-1}]\) (aerosol-aware)
    real(wp), dimension(kms:kme, ims:ime, jms:jme), intent(inout), optional :: nifa !! 3D ice-friendly aerosol number mixing ratio \([kg^{-1}]\) (aerosol-aware)
    real(wp), dimension(kms:kme, ims:ime, jms:jme), intent(inout), optional :: qb !! 3D graupel volume mixing ratio \([m^{-3}\; kg^{-1}]\) (hail-aware)
    real(wp), dimension(kms:kme, ims:ime, jms:jme), intent(inout), optional :: ng !! 3D graupel number mixing ratio \([kg^{-1}]\) (hail-aware)

    ! additional optional arguments
    real(wp), dimension(kms:kme, ims:ime, jms:jme), intent(inout), optional :: qcfrac
    real(wp), dimension(kms:kme, ims:ime, jms:jme), intent(inout), optional :: qifrac
    real(wp), dimension(kms:kme, ims:ime, jms:jme), intent(in), optional :: qc_bl
    real(wp), dimension(kms:kme, ims:ime, jms:jme), intent(in), optional :: qcfrac_bl
    real(wp), dimension(kms:kme, ims:ime, jms:jme), intent(in), optional :: thten_bl
    real(wp), dimension(kms:kme, ims:ime, jms:jme), intent(in), optional :: qvten_bl
    real(wp), dimension(kms:kme, ims:ime, jms:jme), intent(in), optional :: qcten_bl
    real(wp), dimension(kms:kme, ims:ime, jms:jme), intent(in), optional :: qiten_bl
    real(wp), dimension(kms:kme, ims:ime, jms:jme), intent(in), optional :: thten_lwrad
    real(wp), dimension(kms:kme, ims:ime, jms:jme), intent(in), optional :: thten_swrad

    !! 1-based block workspace, sized block_stride x block_stride. Host fields are packed
    !! directly into these per sub-tile (the pack plays the role of the gather), so no separate
    !! full-tile workspace is needed.
    real(wp), dimension(:,:,:), allocatable :: t3d  !! block workspace temperature \([K]\)
    real(wp), dimension(:,:,:), allocatable :: p3d  !! block workspace pressure \([Pa]\)
    real(wp), dimension(:,:,:), allocatable :: qv3d !! block workspace water vapor mixing ratio \([kg\; kg^{-1}]\)
    real(wp), dimension(:,:,:), allocatable :: qc3d !! block workspace cloud water mass mixing ratio \([kg\; kg^{-1}]\)
    real(wp), dimension(:,:,:), allocatable :: qr3d !! block workspace rain water mass mixing ratio \([kg\; kg^{-1}]\)
    real(wp), dimension(:,:,:), allocatable :: qi3d !! block workspace cloud ice mass mixing ratio \([kg\; kg^{-1}]\)
    real(wp), dimension(:,:,:), allocatable :: qs3d !! block workspace snow mass mixing ratio \([kg\; kg^{-1}]\)
    real(wp), dimension(:,:,:), allocatable :: qg3d !! block workspace graupel mass mixing ratio \([kg\; kg^{-1}]\)
    real(wp), dimension(:,:,:), allocatable :: ni3d !! block workspace cloud ice number mixing ratio \([kg^{-1}]\)
    real(wp), dimension(:,:,:), allocatable :: nr3d !! block workspace rain water number mixing ratio \([kg^{-1}]\)
    real(wp), dimension(:,:,:), allocatable :: w3d  !! block workspace vertical velocity \(m\; s^{-1}]\)
    real(wp), dimension(:,:,:), allocatable :: dz3d !! block workspace vertical grid spacing \([m]\)

    real(wp), dimension(:,:,:), allocatable, target :: nc3d !! 3D workspace cloud water number (aerosol-aware)
    real(wp), dimension(:,:,:), allocatable, target :: nwfa3d !! 3D workspace water-friendly aerosol number (aerosol-aware)
    real(wp), dimension(:,:,:), allocatable, target :: nifa3d !! 3D workspace ice-friendly aerosol number (aerosol-aware)
    real(wp), dimension(:,:,:), allocatable, target :: qb3d !! 3D workspace graupel volume (hail-aware)
    real(wp), dimension(:,:,:), allocatable, target :: ng3d !! 3D workspace graupel number (hail-aware)

    real(wp), dimension(:,:,:), allocatable, target :: qcfrac3d
    real(wp), dimension(:,:,:), allocatable, target :: qifrac3d
    real(wp), dimension(:,:,:), allocatable, target :: qc_bl3d
    real(wp), dimension(:,:,:), allocatable, target :: qcfrac_bl3d
    real(wp), dimension(:,:,:), allocatable, target :: thten_bl3d
    real(wp), dimension(:,:,:), allocatable, target :: qvten_bl3d
    real(wp), dimension(:,:,:), allocatable, target :: qcten_bl3d
    real(wp), dimension(:,:,:), allocatable, target :: qiten_bl3d
    real(wp), dimension(:,:,:), allocatable, target :: thten_lwrad3d
    real(wp), dimension(:,:,:), allocatable, target :: thten_swrad3d

    integer :: i, j, k, kk
    integer :: i_out, j_out, i0, i1, j0, j1, ii, jj, bi, bj, bs !! hybrid stride-block bounds (packed bs x bs sub-tile)
    real(wp) :: hail_col_max
    logical :: use_temperature

    type(ty_tempo_main_diags) :: tempo_main_diags
    type(ty_tempo_driver_diags), intent(inout) :: tempo_diags
    logical, intent(in), optional :: arguments_on_device
    integer, intent(in), optional :: stride !! hybrid horizontal block size (1=per-column, ncells=full plane); <=0/absent => full tile

    logical :: xfer_arguments
    integer :: block_stride !! effective horizontal block size used to tile the tempo_main calls
    integer :: ncols_i, ncols_j
    character(len=96) :: nvtx_lbl_run

    write(nvtx_lbl_run, '(A,I0,A,ES11.4)') 'tempo_run it=', itimestep, ' dt=', dt
    call nvtx_range_push(trim(nvtx_lbl_run))

    ! Resolve the hybrid horizontal block size (stride).
    ! stride=1 => per-column (column/CPU form); stride>=tile extent => full plane (plane/GPU form).
    ! Absent or <=0 defaults to the full tile (preserving the plane behavior).
    ncols_i = ite - its + 1
    ncols_j = jte - jts + 1
    block_stride = max(ncols_i, ncols_j)
    if (present(stride)) then
      if (stride > 0) block_stride = min(stride, max(ncols_i, ncols_j))
    endif
#ifdef TEMPO_STRIDE1
    ! Compile-time per-column build: tempo_main (and its block-context helpers) are
    ! specialized to a hardcoded horizontal extent of 1, so force the per-column block
    ! regardless of the requested stride. Every invocation runs the per-column form.
    block_stride = 1
#endif
    bs = block_stride

    ! allocate the 1-based block workspace (block_stride x block_stride); optional model
    ! fields only when present so an unallocated allocatable is cleanly "not present".
    allocate(t3d(kts:kte,bs,bs), p3d(kts:kte,bs,bs), qv3d(kts:kte,bs,bs), qc3d(kts:kte,bs,bs), &
      qr3d(kts:kte,bs,bs), qi3d(kts:kte,bs,bs), qs3d(kts:kte,bs,bs), qg3d(kts:kte,bs,bs), &
      ni3d(kts:kte,bs,bs), nr3d(kts:kte,bs,bs), w3d(kts:kte,bs,bs), dz3d(kts:kte,bs,bs))
    if (present(nwfa)) allocate(nwfa3d(kts:kte,bs,bs), source=0._wp)
    if (present(nifa)) allocate(nifa3d(kts:kte,bs,bs), source=0._wp)
    if (present(nc)) allocate(nc3d(kts:kte,bs,bs), source=0._wp)
    if (present(ng)) allocate(ng3d(kts:kte,bs,bs), source=0._wp)
    if (present(qb)) allocate(qb3d(kts:kte,bs,bs), source=0._wp)

    if (present(qcfrac)) allocate(qcfrac3d(kts:kte,bs,bs), source=0._wp)
    if (present(qifrac)) allocate(qifrac3d(kts:kte,bs,bs), source=0._wp)
    if (present(qc_bl)) allocate(qc_bl3d(kts:kte,bs,bs), source=0._wp)
    if (present(qcfrac_bl)) allocate(qcfrac_bl3d(kts:kte,bs,bs), source=0._wp)
    if (present(thten_bl)) allocate(thten_bl3d(kts:kte,bs,bs), source=0._wp)
    if (present(qvten_bl)) allocate(qvten_bl3d(kts:kte,bs,bs), source=0._wp)
    if (present(qcten_bl)) allocate(qcten_bl3d(kts:kte,bs,bs), source=0._wp)
    if (present(qiten_bl)) allocate(qiten_bl3d(kts:kte,bs,bs), source=0._wp)
    if (present(thten_lwrad)) allocate(thten_lwrad3d(kts:kte,bs,bs), source=0._wp)
    if (present(thten_swrad)) allocate(thten_swrad3d(kts:kte,bs,bs), source=0._wp)

    ! per-block diagnostics (block-sized; scattered straight into tempo_diags below)
    call tempo_run_setup_main_diags(tempo_cfgs, tempo_main_diags, kts, kte, 1, bs, 1, bs)

    ! temperature or theta and exner
    if (present(t)) then
      use_temperature = .true.
    elseif (present(th) .and. present(pii)) then
      use_temperature = .false.
    else
      error stop "tempo_run() --- requires either temperature or theta and Exner function"
    endif

    xfer_arguments = .true.
    if (present(arguments_on_device)) xfer_arguments = .not. arguments_on_device

    !! ============================================================
    !! Explicit OpenACC data movement (replaces -gpu=mem:managed)
    !!
    !! Required host inputs/outputs + 3D workspace. Optional dummy arguments cannot be passed
    !! unconditionally to a data clause when absent, so they are handled in dedicated
    !! present(...)-guarded enter/exit data blocks below. Matching exits run before the routine
    !! returns to copy results out and release device images.
    !! ============================================================
    call nvtx_range_push('tempo_run_enter_data')
    if (xfer_arguments) then
      call tempo_run_enter_data(tempo_cfgs, tempo_diags, use_temperature, &
        t=t, th=th, pii=pii, p=p, w=w, dz=dz, qv=qv, qc=qc, qr=qr, qi=qi, qs=qs, qg=qg, ni=ni, nr=nr, &
        nc=nc, nwfa=nwfa, nifa=nifa, ng=ng, qb=qb, qcfrac=qcfrac, qifrac=qifrac, qc_bl=qc_bl, &
        qcfrac_bl=qcfrac_bl, thten_bl=thten_bl, qvten_bl=qvten_bl, qcten_bl=qcten_bl, qiten_bl=qiten_bl, &
        thten_lwrad=thten_lwrad, thten_swrad=thten_swrad, kts=kts, kte=kte, its=its, ite=ite, jts=jts, jte=jte)
    endif
    !$acc enter data create(t3d, p3d, qv3d, qc3d, qr3d, qi3d, qs3d, qg3d, ni3d, nr3d, w3d, dz3d)

    if (present(nc)) then
      if (xfer_arguments) then
        !$acc enter data copyin(nc) create(nc3d)
      else
        !$acc enter data create(nc3d)
      endif
    endif
    if (present(nwfa)) then
      if (xfer_arguments) then
        !$acc enter data copyin(nwfa) create(nwfa3d)
      else
        !$acc enter data create(nwfa3d)
      endif
    endif
    if (present(nifa)) then
      if (xfer_arguments) then
        !$acc enter data copyin(nifa) create(nifa3d)
      else
        !$acc enter data create(nifa3d)
      endif
    endif
    if (present(ng)) then
      if (xfer_arguments) then
        !$acc enter data copyin(ng) create(ng3d)
      else
        !$acc enter data create(ng3d)
      endif
    endif
    if (present(qb)) then
      if (xfer_arguments) then
        !$acc enter data copyin(qb) create(qb3d)
      else
        !$acc enter data create(qb3d)
      endif
    endif
    if (present(qcfrac)) then
      if (xfer_arguments) then
        !$acc enter data copyin(qcfrac) create(qcfrac3d)
      else
        !$acc enter data create(qcfrac3d)
      endif
    endif
    if (present(qifrac)) then
      if (xfer_arguments) then
        !$acc enter data copyin(qifrac) create(qifrac3d)
      else
        !$acc enter data create(qifrac3d)
      endif
    endif
    if (present(qc_bl)) then
      if (xfer_arguments) then
        !$acc enter data copyin(qc_bl) create(qc_bl3d)
      else
        !$acc enter data create(qc_bl3d)
      endif
    endif
    if (present(qcfrac_bl)) then
      if (xfer_arguments) then
        !$acc enter data copyin(qcfrac_bl) create(qcfrac_bl3d)
      else
        !$acc enter data create(qcfrac_bl3d)
      endif
    endif
    if (present(thten_bl)) then
      if (xfer_arguments) then
        !$acc enter data copyin(thten_bl) create(thten_bl3d)
      else
        !$acc enter data create(thten_bl3d)
      endif
    endif
    if (present(qvten_bl)) then
      if (xfer_arguments) then
        !$acc enter data copyin(qvten_bl) create(qvten_bl3d)
      else
        !$acc enter data create(qvten_bl3d)
      endif
    endif
    if (present(qcten_bl)) then
      if (xfer_arguments) then
        !$acc enter data copyin(qcten_bl) create(qcten_bl3d)
      else
        !$acc enter data create(qcten_bl3d)
      endif
    endif
    if (present(qiten_bl)) then
      if (xfer_arguments) then
        !$acc enter data copyin(qiten_bl) create(qiten_bl3d)
      else
        !$acc enter data create(qiten_bl3d)
      endif
    endif
    if (present(thten_lwrad)) then
      if (xfer_arguments) then
        !$acc enter data copyin(thten_lwrad) create(thten_lwrad3d)
      else
        !$acc enter data create(thten_lwrad3d)
      endif
    endif
    if (present(thten_swrad)) then
      if (xfer_arguments) then
        !$acc enter data copyin(thten_swrad) create(thten_swrad3d)
      else
        !$acc enter data create(thten_swrad3d)
      endif
    endif

    !$acc enter data copyin(tempo_main_diags)
    !$acc enter data create(tempo_main_diags%rain_precip, tempo_main_diags%cloud_precip, &
    !$acc                   tempo_main_diags%ice_liquid_equiv_precip, tempo_main_diags%snow_liquid_equiv_precip, &
    !$acc                   tempo_main_diags%graupel_liquid_equiv_precip, tempo_main_diags%frozen_fraction, &
    !$acc                   tempo_main_diags%frz_rain_precip)
    if (allocated(tempo_main_diags%cloud_number_mixing_ratio)) then
      !$acc enter data create(tempo_main_diags%cloud_number_mixing_ratio)
    endif
    if (allocated(tempo_main_diags%rain_med_vol_diam)) then
      !$acc enter data create(tempo_main_diags%rain_med_vol_diam)
    endif
    if (allocated(tempo_main_diags%graupel_med_vol_diam)) then
      !$acc enter data create(tempo_main_diags%graupel_med_vol_diam)
    endif
    if (allocated(tempo_main_diags%max_hail_diameter)) then
      !$acc enter data create(tempo_main_diags%max_hail_diameter)
    endif
    if (allocated(tempo_main_diags%refl10cm)) then
      !$acc enter data create(tempo_main_diags%refl10cm)
    endif
    if (allocated(tempo_main_diags%re_cloud)) then
      !$acc enter data create(tempo_main_diags%re_cloud, tempo_main_diags%re_ice, tempo_main_diags%re_snow)
    endif
    call nvtx_range_pop()


    if (present(arguments_on_device)) then
      if (arguments_on_device) then
        call tempo_run_setup_diags(tempo_cfgs, tempo_diags, kts, kte, its, ite, jts, jte)
        ! Device diagnostics are overwritten below by the per-block scatter kernels, so avoid
        ! resetting them on every timestep.
      endif
    endif

    ! (2) Microphysics over the tile, processed in block_stride x block_stride sub-tiles.
    !! For each sub-tile: pack the host fields straight into the 1-based block workspace (the pack
    !! plays the role of the gather), run tempo_main on the dense block (its=1..bi, jts=1..bj), then
    !! unpack the updated fields back to the host model arrays and scatter the per-block diagnostics
    !! into tempo_diags, all at absolute tile positions. tempo_main_diags is the block-sized
    !! diagnostic workspace, reset each block. Columns are horizontally independent, so any block
    !! size yields identical numerics; it only trades kernel granularity vs launch count.
    !!   block_stride == 1          -> one tempo_main call per column (column/CPU form).
    !!   block_stride >= tile extent -> a single block spanning the whole tile (plane/GPU form).
    call nvtx_range_push('tempo_main')
      do j_out = jts, jte, bs
        j0 = j_out
        j1 = min(j_out + bs - 1, jte)
        bj = j1 - j0 + 1
        do i_out = its, ite, bs
          i0 = i_out
          i1 = min(i_out + bs - 1, ite)
          bi = i1 - i0 + 1

          ! pack host sub-tile -> 1-based block workspace (the pack is the gather; optionals guarded)
          call nvtx_range_push('tempo_run_pack_host_workspace')
          !$acc parallel
          !$acc loop gang vector collapse(3)
          do jj = 1, bj
            do ii = 1, bi
              do k = kts, kte
                if (use_temperature) then
                  t3d(k,ii,jj) = t(k,i0+ii-1,j0+jj-1)
                else
                  t3d(k,ii,jj) = th(k,i0+ii-1,j0+jj-1) * pii(k,i0+ii-1,j0+jj-1)
                endif
                p3d(k,ii,jj)  = p(k,i0+ii-1,j0+jj-1)
                w3d(k,ii,jj)  = w(k,i0+ii-1,j0+jj-1)
                dz3d(k,ii,jj) = dz(k,i0+ii-1,j0+jj-1)
                qv3d(k,ii,jj) = qv(k,i0+ii-1,j0+jj-1)
                qc3d(k,ii,jj) = qc(k,i0+ii-1,j0+jj-1)
                qi3d(k,ii,jj) = qi(k,i0+ii-1,j0+jj-1)
                qr3d(k,ii,jj) = qr(k,i0+ii-1,j0+jj-1)
                qs3d(k,ii,jj) = qs(k,i0+ii-1,j0+jj-1)
                qg3d(k,ii,jj) = qg(k,i0+ii-1,j0+jj-1)
                ni3d(k,ii,jj) = ni(k,i0+ii-1,j0+jj-1)
                nr3d(k,ii,jj) = nr(k,i0+ii-1,j0+jj-1)
                if (present(nc))   nc3d(k,ii,jj)   = nc(k,i0+ii-1,j0+jj-1)
                if (present(nwfa)) nwfa3d(k,ii,jj) = nwfa(k,i0+ii-1,j0+jj-1)
                if (present(nifa)) nifa3d(k,ii,jj) = nifa(k,i0+ii-1,j0+jj-1)
                if (present(ng) .and. present(qb)) then
                  ng3d(k,ii,jj) = ng(k,i0+ii-1,j0+jj-1)
                  qb3d(k,ii,jj) = qb(k,i0+ii-1,j0+jj-1)
                endif
                if (present(qc_bl) .and. present(qcfrac_bl)) then
                  qc_bl3d(k,ii,jj)     = qc_bl(k,i0+ii-1,j0+jj-1)
                  qcfrac_bl3d(k,ii,jj) = qcfrac_bl(k,i0+ii-1,j0+jj-1)
                endif
              enddo
            enddo
          enddo
          !$acc end parallel
          call nvtx_range_pop()

          ! reset per-block diagnostics to their initial values (matches the fresh source=
          ! initialization the diagnostic allocation gets each timestep; block reuse must match it)
          call tempo_run_reset_block_diags(tempo_main_diags, kts, kte, bi, bj)

          call tempo_main(tempo_cfgs=tempo_cfgs, &
            qv3d=qv3d, qc3d=qc3d, qi3d=qi3d, qr3d=qr3d, qs3d=qs3d, qg3d=qg3d, qb3d=qb3d, &
            ni3d=ni3d, nr3d=nr3d, nc3d=nc3d, ng3d=ng3d, nwfa3d=nwfa3d, nifa3d=nifa3d, &
            t3d=t3d, p3d=p3d, w3d=w3d, dz3d=dz3d, &
            qcfrac3d=qcfrac3d, qifrac3d=qifrac3d, qc_bl3d=qc_bl3d, qcfrac_bl3d=qcfrac_bl3d, &
            thten_bl3d=thten_bl3d, qvten_bl3d=qvten_bl3d, qcten_bl3d=qcten_bl3d, qiten_bl3d=qiten_bl3d, &
            thten_lwrad3d=thten_lwrad3d, thten_swrad3d=thten_swrad3d, &
            kts=kts, kte=kte, dt=dt, its=1, ite=bi, jts=1, jte=bj, tempo_main_diags=tempo_main_diags)

          ! unpack updated block fields back to the host model arrays at absolute positions.
          ! Split on use_temperature so the kernel only touches device-mapped fields (t XOR th,pii).
          call nvtx_range_push('tempo_run_unpack_host_workspace')
          if (use_temperature) then
            !$acc parallel
            !$acc loop gang vector collapse(3)
            do jj = 1, bj
              do ii = 1, bi
                do k = kts, kte
                  if (present(nc))   nc(k,i0+ii-1,j0+jj-1)   = nc3d(k,ii,jj)
                  if (present(nwfa)) nwfa(k,i0+ii-1,j0+jj-1) = nwfa3d(k,ii,jj)
                  if (present(nifa)) nifa(k,i0+ii-1,j0+jj-1) = nifa3d(k,ii,jj)
                  if (present(ng) .and. present(qb)) then
                    ng(k,i0+ii-1,j0+jj-1) = ng3d(k,ii,jj)
                    qb(k,i0+ii-1,j0+jj-1) = qb3d(k,ii,jj)
                  endif
                  qv(k,i0+ii-1,j0+jj-1) = qv3d(k,ii,jj)
                  qc(k,i0+ii-1,j0+jj-1) = qc3d(k,ii,jj)
                  qi(k,i0+ii-1,j0+jj-1) = qi3d(k,ii,jj)
                  qr(k,i0+ii-1,j0+jj-1) = qr3d(k,ii,jj)
                  qs(k,i0+ii-1,j0+jj-1) = qs3d(k,ii,jj)
                  qg(k,i0+ii-1,j0+jj-1) = qg3d(k,ii,jj)
                  ni(k,i0+ii-1,j0+jj-1) = ni3d(k,ii,jj)
                  nr(k,i0+ii-1,j0+jj-1) = nr3d(k,ii,jj)
                  t(k,i0+ii-1,j0+jj-1)  = t3d(k,ii,jj)
                enddo
              enddo
            enddo
            !$acc end parallel
          else
            !$acc parallel
            !$acc loop gang vector collapse(3)
            do jj = 1, bj
              do ii = 1, bi
                do k = kts, kte
                  if (present(nc))   nc(k,i0+ii-1,j0+jj-1)   = nc3d(k,ii,jj)
                  if (present(nwfa)) nwfa(k,i0+ii-1,j0+jj-1) = nwfa3d(k,ii,jj)
                  if (present(nifa)) nifa(k,i0+ii-1,j0+jj-1) = nifa3d(k,ii,jj)
                  if (present(ng) .and. present(qb)) then
                    ng(k,i0+ii-1,j0+jj-1) = ng3d(k,ii,jj)
                    qb(k,i0+ii-1,j0+jj-1) = qb3d(k,ii,jj)
                  endif
                  qv(k,i0+ii-1,j0+jj-1) = qv3d(k,ii,jj)
                  qc(k,i0+ii-1,j0+jj-1) = qc3d(k,ii,jj)
                  qi(k,i0+ii-1,j0+jj-1) = qi3d(k,ii,jj)
                  qr(k,i0+ii-1,j0+jj-1) = qr3d(k,ii,jj)
                  qs(k,i0+ii-1,j0+jj-1) = qs3d(k,ii,jj)
                  qg(k,i0+ii-1,j0+jj-1) = qg3d(k,ii,jj)
                  ni(k,i0+ii-1,j0+jj-1) = ni3d(k,ii,jj)
                  nr(k,i0+ii-1,j0+jj-1) = nr3d(k,ii,jj)
                  th(k,i0+ii-1,j0+jj-1) = t3d(k,ii,jj) / pii(k,i0+ii-1,j0+jj-1)
                enddo
              enddo
            enddo
            !$acc end parallel
          endif
          call nvtx_range_pop()

          ! scatter the per-block diagnostics straight into tempo_diags at absolute positions
          call nvtx_range_push('tempo_run_copy_driver_diags')
          if (allocated(tempo_diags%rain_precip)) then
            !$acc parallel
            !$acc loop gang vector collapse(2)
            do jj = 1, bj
              do ii = 1, bi
                tempo_diags%rain_precip(i0+ii-1,j0+jj-1)                 = tempo_main_diags%rain_precip(ii,jj)
                tempo_diags%ice_liquid_equiv_precip(i0+ii-1,j0+jj-1)     = tempo_main_diags%ice_liquid_equiv_precip(ii,jj)
                tempo_diags%snow_liquid_equiv_precip(i0+ii-1,j0+jj-1)    = tempo_main_diags%snow_liquid_equiv_precip(ii,jj)
                tempo_diags%graupel_liquid_equiv_precip(i0+ii-1,j0+jj-1) = tempo_main_diags%graupel_liquid_equiv_precip(ii,jj)
                tempo_diags%frozen_fraction(i0+ii-1,j0+jj-1)            = tempo_main_diags%frozen_fraction(ii,jj)
                tempo_diags%frz_rain_precip(i0+ii-1,j0+jj-1)            = tempo_main_diags%frz_rain_precip(ii,jj)
              enddo
            enddo
            !$acc end parallel
          endif

          if (allocated(tempo_diags%cloud_number_mixing_ratio) .and. allocated(tempo_main_diags%cloud_number_mixing_ratio)) then
            !$acc parallel
            !$acc loop gang vector collapse(3)
            do jj = 1, bj
              do ii = 1, bi
                do k = kts, kte
                  tempo_diags%cloud_number_mixing_ratio(k,i0+ii-1,j0+jj-1) = tempo_main_diags%cloud_number_mixing_ratio(k,ii,jj)
                enddo
              enddo
            enddo
            !$acc end parallel
          endif
          if (allocated(tempo_diags%rain_med_vol_diam) .and. allocated(tempo_main_diags%rain_med_vol_diam)) then
            !$acc parallel
            !$acc loop gang vector collapse(3)
            do jj = 1, bj
              do ii = 1, bi
                do k = kts, kte
                  tempo_diags%rain_med_vol_diam(k,i0+ii-1,j0+jj-1) = tempo_main_diags%rain_med_vol_diam(k,ii,jj)
                enddo
              enddo
            enddo
            !$acc end parallel
          endif
          if (allocated(tempo_diags%graupel_med_vol_diam) .and. allocated(tempo_main_diags%graupel_med_vol_diam)) then
            !$acc parallel
            !$acc loop gang vector collapse(3)
            do jj = 1, bj
              do ii = 1, bi
                do k = kts, kte
                  tempo_diags%graupel_med_vol_diam(k,i0+ii-1,j0+jj-1) = tempo_main_diags%graupel_med_vol_diam(k,ii,jj)
                enddo
              enddo
            enddo
            !$acc end parallel
          endif
          if (allocated(tempo_diags%re_cloud) .and. allocated(tempo_main_diags%re_cloud)) then
            !$acc parallel
            !$acc loop gang vector collapse(3)
            do jj = 1, bj
              do ii = 1, bi
                do k = kts, kte
                  tempo_diags%re_cloud(k,i0+ii-1,j0+jj-1) = tempo_main_diags%re_cloud(k,ii,jj)
                  tempo_diags%re_ice(k,i0+ii-1,j0+jj-1)   = tempo_main_diags%re_ice(k,ii,jj)
                  tempo_diags%re_snow(k,i0+ii-1,j0+jj-1)  = tempo_main_diags%re_snow(k,ii,jj)
                enddo
              enddo
            enddo
            !$acc end parallel
          endif
          if (allocated(tempo_diags%refl10cm) .and. allocated(tempo_main_diags%refl10cm)) then
            !$acc parallel
            !$acc loop gang vector collapse(3)
            do jj = 1, bj
              do ii = 1, bi
                do k = kts, kte
                  tempo_diags%refl10cm(k,i0+ii-1,j0+jj-1) = tempo_main_diags%refl10cm(k,ii,jj)
                enddo
              enddo
            enddo
            !$acc end parallel
          endif
          if (allocated(tempo_diags%max_hail_diameter_sfc) .and. allocated(tempo_diags%max_hail_diameter_column) &
              .and. allocated(tempo_main_diags%max_hail_diameter)) then
            !$acc parallel
            !$acc loop gang vector collapse(2) private(hail_col_max)
            do jj = 1, bj
              do ii = 1, bi
                tempo_diags%max_hail_diameter_sfc(i0+ii-1,j0+jj-1) = tempo_main_diags%max_hail_diameter(kts,ii,jj)
                hail_col_max = tempo_main_diags%max_hail_diameter(kts,ii,jj)
                !$acc loop seq
                do kk = kts + 1, kte
                  hail_col_max = max(hail_col_max, tempo_main_diags%max_hail_diameter(kk,ii,jj))
                enddo
                tempo_diags%max_hail_diameter_column(i0+ii-1,j0+jj-1) = hail_col_max
              enddo
            enddo
            !$acc end parallel
          endif
          call nvtx_range_pop()
        enddo
      enddo

    call nvtx_range_pop()

    !! ============================================================
    !! Explicit OpenACC data movement: exit / copyout
    !! Mirror the enter-data block above. tempo_main_diags%* are device-only working
    !! buffers (delete, no copyout). tempo_diags%* + host inout fields are copied back.
    !! Order: most-derived first so device images outlive any final references.
    !! ============================================================
    call nvtx_range_push('tempo_run_exit_data')

    if (allocated(tempo_main_diags%re_cloud)) then
      !$acc exit data delete(tempo_main_diags%re_cloud, tempo_main_diags%re_ice, tempo_main_diags%re_snow)
    endif
    if (allocated(tempo_main_diags%refl10cm)) then
      !$acc exit data delete(tempo_main_diags%refl10cm)
    endif
    if (allocated(tempo_main_diags%max_hail_diameter)) then
      !$acc exit data delete(tempo_main_diags%max_hail_diameter)
    endif
    if (allocated(tempo_main_diags%graupel_med_vol_diam)) then
      !$acc exit data delete(tempo_main_diags%graupel_med_vol_diam)
    endif
    if (allocated(tempo_main_diags%rain_med_vol_diam)) then
      !$acc exit data delete(tempo_main_diags%rain_med_vol_diam)
    endif
    if (allocated(tempo_main_diags%cloud_number_mixing_ratio)) then
      !$acc exit data delete(tempo_main_diags%cloud_number_mixing_ratio)
    endif
    !$acc exit data delete(tempo_main_diags%rain_precip, tempo_main_diags%cloud_precip, &
    !$acc                  tempo_main_diags%ice_liquid_equiv_precip, tempo_main_diags%snow_liquid_equiv_precip, &
    !$acc                  tempo_main_diags%graupel_liquid_equiv_precip, tempo_main_diags%frozen_fraction, &
    !$acc                  tempo_main_diags%frz_rain_precip)
    !$acc exit data delete(tempo_main_diags)

    if (present(thten_swrad)) then
      !$acc exit data delete(thten_swrad3d)
      if (xfer_arguments) then
        !$acc exit data delete(thten_swrad)
      endif
    endif
    if (present(thten_lwrad)) then
      !$acc exit data delete(thten_lwrad3d)
      if (xfer_arguments) then
        !$acc exit data delete(thten_lwrad)
      endif
    endif
    if (present(qiten_bl)) then
      !$acc exit data delete(qiten_bl3d)
      if (xfer_arguments) then
        !$acc exit data delete(qiten_bl)
      endif
    endif
    if (present(qcten_bl)) then
      !$acc exit data delete(qcten_bl3d)
      if (xfer_arguments) then
        !$acc exit data delete(qcten_bl)
      endif
    endif
    if (present(qvten_bl)) then
      !$acc exit data delete(qvten_bl3d)
      if (xfer_arguments) then
        !$acc exit data delete(qvten_bl)
      endif
    endif
    if (present(thten_bl)) then
      !$acc exit data delete(thten_bl3d)
      if (xfer_arguments) then
        !$acc exit data delete(thten_bl)
      endif
    endif
    if (present(qcfrac_bl)) then
      !$acc exit data delete(qcfrac_bl3d)
      if (xfer_arguments) then
        !$acc exit data delete(qcfrac_bl)
      endif
    endif
    if (present(qc_bl)) then
      !$acc exit data delete(qc_bl3d)
      if (xfer_arguments) then
        !$acc exit data delete(qc_bl)
      endif
    endif
    if (present(qifrac)) then
      !$acc exit data delete(qifrac3d)
      if (xfer_arguments) then
        !$acc exit data copyout(qifrac)
      endif
    endif
    if (present(qcfrac)) then
      !$acc exit data delete(qcfrac3d)
      if (xfer_arguments) then
        !$acc exit data copyout(qcfrac)
      endif
    endif
    if (present(qb)) then
      !$acc exit data delete(qb3d)
      if (xfer_arguments) then
        !$acc exit data copyout(qb)
      endif
    endif
    if (present(ng)) then
      !$acc exit data delete(ng3d)
      if (xfer_arguments) then
        !$acc exit data copyout(ng)
      endif
    endif
    if (present(nifa)) then
      !$acc exit data delete(nifa3d)
      if (xfer_arguments) then
        !$acc exit data copyout(nifa)
      endif
    endif
    if (present(nwfa)) then
      !$acc exit data delete(nwfa3d)
      if (xfer_arguments) then
        !$acc exit data copyout(nwfa)
      endif
    endif
    if (present(nc)) then
      !$acc exit data delete(nc3d)
      if (xfer_arguments) then
        !$acc exit data copyout(nc)
      endif
    endif

    !$acc exit data delete(t3d, p3d, qv3d, qc3d, qr3d, qi3d, qs3d, qg3d, ni3d, nr3d, w3d, dz3d)

    if (xfer_arguments) then
      call tempo_run_exit_data(tempo_cfgs, tempo_diags, use_temperature, &
        t=t, th=th, pii=pii, p=p, w=w, dz=dz, qv=qv, qc=qc, qr=qr, qi=qi, qs=qs, qg=qg, ni=ni, nr=nr, &
        nc=nc, nwfa=nwfa, nifa=nifa, ng=ng, qb=qb, qcfrac=qcfrac, qifrac=qifrac, qc_bl=qc_bl, &
        qcfrac_bl=qcfrac_bl, thten_bl=thten_bl, qvten_bl=qvten_bl, qcten_bl=qcten_bl, qiten_bl=qiten_bl, &
        thten_lwrad=thten_lwrad, thten_swrad=thten_swrad, kts=kts, kte=kte, its=its, ite=ite, jts=jts, jte=jte)
    endif
    call nvtx_range_pop()

    call nvtx_range_pop()

  end subroutine tempo_run

  subroutine tempo_run_enter_data(tempo_cfgs, tempo_diags, use_temperature, &
      t, th, pii, p, w, dz, qv, qc, qr, qi, qs, qg, ni, nr, &
      nc, nwfa, nifa, ng, qb, qcfrac, qifrac, qc_bl, qcfrac_bl, &
      thten_bl, qvten_bl, qcten_bl, qiten_bl, thten_lwrad, thten_swrad, &
      kts, kte, its, ite, jts, jte)
    !! Device data setup for tempo_run dummy arguments (call before the timestep loop).
    type(ty_tempo_cfgs), intent(in) :: tempo_cfgs
    type(ty_tempo_driver_diags), intent(inout) :: tempo_diags
    logical, intent(in) :: use_temperature
    integer, intent(in) :: kts, kte, its, ite, jts, jte
    real(wp), dimension(:,:,:), intent(inout), optional :: t, th
    real(wp), dimension(:,:,:), intent(in), optional :: pii
    real(wp), dimension(:,:,:), intent(in) :: p, w, dz
    real(wp), dimension(:,:,:), intent(inout) :: qv, qc, qr, qi, qs, qg, ni, nr
    real(wp), dimension(:,:,:), intent(inout), optional :: nc, nwfa, nifa, ng, qb
    real(wp), dimension(:,:,:), intent(inout), optional :: qcfrac, qifrac
    real(wp), dimension(:,:,:), intent(in), optional :: qc_bl, qcfrac_bl
    real(wp), dimension(:,:,:), intent(in), optional :: thten_bl, qvten_bl, qcten_bl, qiten_bl
    real(wp), dimension(:,:,:), intent(in), optional :: thten_lwrad, thten_swrad

    call tempo_run_setup_diags(tempo_cfgs, tempo_diags, kts, kte, its, ite, jts, jte)

    !$acc enter data copyin(tempo_cfgs)
    !$acc enter data copyin(p, w, dz) copyin(qv, qc, qr, qi, qs, qg, ni, nr)
    if (use_temperature) then
      !$acc enter data copyin(t)
    else
      !$acc enter data copyin(th) copyin(pii)
    endif
    if (present(nc)) then
      !$acc enter data copyin(nc)
    endif
    if (present(nwfa)) then
      !$acc enter data copyin(nwfa)
    endif
    if (present(nifa)) then
      !$acc enter data copyin(nifa)
    endif
    if (present(ng)) then
      !$acc enter data copyin(ng)
    endif
    if (present(qb)) then
      !$acc enter data copyin(qb)
    endif
    if (present(qcfrac)) then
      !$acc enter data copyin(qcfrac)
    endif
    if (present(qifrac)) then
      !$acc enter data copyin(qifrac)
    endif
    if (present(qc_bl)) then
      !$acc enter data copyin(qc_bl)
    endif
    if (present(qcfrac_bl)) then
      !$acc enter data copyin(qcfrac_bl)
    endif
    if (present(thten_bl)) then
      !$acc enter data copyin(thten_bl)
    endif
    if (present(qvten_bl)) then
      !$acc enter data copyin(qvten_bl)
    endif
    if (present(qcten_bl)) then
      !$acc enter data copyin(qcten_bl)
    endif
    if (present(qiten_bl)) then
      !$acc enter data copyin(qiten_bl)
    endif
    if (present(thten_lwrad)) then
      !$acc enter data copyin(thten_lwrad)
    endif
    if (present(thten_swrad)) then
      !$acc enter data copyin(thten_swrad)
    endif

    !$acc enter data copyin(tempo_diags)
    if (allocated(tempo_diags%rain_precip)) then
      !$acc enter data create(tempo_diags%rain_precip, tempo_diags%ice_liquid_equiv_precip, &
      !$acc            tempo_diags%snow_liquid_equiv_precip, tempo_diags%graupel_liquid_equiv_precip, &
      !$acc            tempo_diags%frozen_fraction, tempo_diags%frz_rain_precip)
    endif
    if (allocated(tempo_diags%cloud_number_mixing_ratio)) then
      !$acc enter data create(tempo_diags%cloud_number_mixing_ratio)
    endif
    if (allocated(tempo_diags%rain_med_vol_diam)) then
      !$acc enter data create(tempo_diags%rain_med_vol_diam)
    endif
    if (allocated(tempo_diags%graupel_med_vol_diam)) then
      !$acc enter data create(tempo_diags%graupel_med_vol_diam)
    endif
    if (allocated(tempo_diags%refl10cm)) then
      !$acc enter data create(tempo_diags%refl10cm)
    endif
    if (allocated(tempo_diags%re_cloud)) then
      !$acc enter data create(tempo_diags%re_cloud, tempo_diags%re_ice, tempo_diags%re_snow)
    endif
    if (allocated(tempo_diags%max_hail_diameter_sfc)) then
      !$acc enter data create(tempo_diags%max_hail_diameter_sfc, tempo_diags%max_hail_diameter_column)
    endif
    call tempo_run_reset_diags_device(tempo_diags, kts, kte, its, ite, jts, jte)
  end subroutine tempo_run_enter_data


  subroutine tempo_run_exit_data(tempo_cfgs, tempo_diags, use_temperature, &
      t, th, pii, p, w, dz, qv, qc, qr, qi, qs, qg, ni, nr, &
      nc, nwfa, nifa, ng, qb, qcfrac, qifrac, qc_bl, qcfrac_bl, &
      thten_bl, qvten_bl, qcten_bl, qiten_bl, thten_lwrad, thten_swrad, &
      kts, kte, its, ite, jts, jte)
    !! Device-to-host transfer for tempo_run dummy arguments (call after the timestep loop).
    type(ty_tempo_cfgs), intent(in) :: tempo_cfgs
    type(ty_tempo_driver_diags), intent(inout) :: tempo_diags
    logical, intent(in) :: use_temperature
    integer, intent(in) :: kts, kte, its, ite, jts, jte
    real(wp), dimension(:,:,:), intent(inout), optional :: t, th
    real(wp), dimension(:,:,:), intent(in), optional :: pii
    real(wp), dimension(:,:,:), intent(in) :: p, w, dz
    real(wp), dimension(:,:,:), intent(inout) :: qv, qc, qr, qi, qs, qg, ni, nr
    real(wp), dimension(:,:,:), intent(inout), optional :: nc, nwfa, nifa, ng, qb
    real(wp), dimension(:,:,:), intent(inout), optional :: qcfrac, qifrac
    real(wp), dimension(:,:,:), intent(in), optional :: qc_bl, qcfrac_bl
    real(wp), dimension(:,:,:), intent(in), optional :: thten_bl, qvten_bl, qcten_bl, qiten_bl
    real(wp), dimension(:,:,:), intent(in), optional :: thten_lwrad, thten_swrad

    if (allocated(tempo_diags%max_hail_diameter_sfc)) then
      !$acc exit data copyout(tempo_diags%max_hail_diameter_sfc, tempo_diags%max_hail_diameter_column)
    endif
    if (allocated(tempo_diags%re_cloud)) then
      !$acc exit data copyout(tempo_diags%re_cloud, tempo_diags%re_ice, tempo_diags%re_snow)
    endif
    if (allocated(tempo_diags%refl10cm)) then
      !$acc exit data copyout(tempo_diags%refl10cm)
    endif
    if (allocated(tempo_diags%graupel_med_vol_diam)) then
      !$acc exit data copyout(tempo_diags%graupel_med_vol_diam)
    endif
    if (allocated(tempo_diags%rain_med_vol_diam)) then
      !$acc exit data copyout(tempo_diags%rain_med_vol_diam)
    endif
    if (allocated(tempo_diags%cloud_number_mixing_ratio)) then
      !$acc exit data copyout(tempo_diags%cloud_number_mixing_ratio)
    endif
    if (allocated(tempo_diags%rain_precip)) then
      !$acc exit data copyout(tempo_diags%rain_precip, tempo_diags%ice_liquid_equiv_precip, &
      !$acc               tempo_diags%snow_liquid_equiv_precip, tempo_diags%graupel_liquid_equiv_precip, &
      !$acc               tempo_diags%frozen_fraction, tempo_diags%frz_rain_precip)
    endif
    !$acc exit data delete(tempo_diags)

    if (present(thten_swrad)) then
      !$acc exit data delete(thten_swrad)
    endif
    if (present(thten_lwrad)) then
      !$acc exit data delete(thten_lwrad)
    endif
    if (present(qiten_bl)) then
      !$acc exit data delete(qiten_bl)
    endif
    if (present(qcten_bl)) then
      !$acc exit data delete(qcten_bl)
    endif
    if (present(qvten_bl)) then
      !$acc exit data delete(qvten_bl)
    endif
    if (present(thten_bl)) then
      !$acc exit data delete(thten_bl)
    endif
    if (present(qcfrac_bl)) then
      !$acc exit data delete(qcfrac_bl)
    endif
    if (present(qc_bl)) then
      !$acc exit data delete(qc_bl)
    endif
    if (present(qifrac)) then
      !$acc exit data copyout(qifrac)
    endif
    if (present(qcfrac)) then
      !$acc exit data copyout(qcfrac)
    endif
    if (present(qb)) then
      !$acc exit data copyout(qb)
    endif
    if (present(ng)) then
      !$acc exit data copyout(ng)
    endif
    if (present(nifa)) then
      !$acc exit data copyout(nifa)
    endif
    if (present(nwfa)) then
      !$acc exit data copyout(nwfa)
    endif
    if (present(nc)) then
      !$acc exit data copyout(nc)
    endif
    if (use_temperature) then
      !$acc exit data copyout(t)
    else
      !$acc exit data copyout(th) delete(pii)
    endif
    !$acc exit data copyout(qv, qc, qr, qi, qs, qg, ni, nr) delete(p, w, dz)
    !$acc exit data delete(tempo_cfgs)
  end subroutine tempo_run_exit_data


  subroutine tempo_run_sync_host_diags(tempo_diags)
    type(ty_tempo_driver_diags), intent(inout) :: tempo_diags
    if (allocated(tempo_diags%max_hail_diameter_sfc)) then
      !$acc update host(tempo_diags%max_hail_diameter_sfc, tempo_diags%max_hail_diameter_column)
    endif
    if (allocated(tempo_diags%re_cloud)) then
      !$acc update host(tempo_diags%re_cloud, tempo_diags%re_ice, tempo_diags%re_snow)
    endif
    if (allocated(tempo_diags%refl10cm)) then
      !$acc update host(tempo_diags%refl10cm)
    endif
    if (allocated(tempo_diags%graupel_med_vol_diam)) then
      !$acc update host(tempo_diags%graupel_med_vol_diam)
    endif
    if (allocated(tempo_diags%rain_med_vol_diam)) then
      !$acc update host(tempo_diags%rain_med_vol_diam)
    endif
    if (allocated(tempo_diags%cloud_number_mixing_ratio)) then
      !$acc update host(tempo_diags%cloud_number_mixing_ratio)
    endif
    if (allocated(tempo_diags%rain_precip)) then
      !$acc update host(tempo_diags%rain_precip, tempo_diags%ice_liquid_equiv_precip, &
      !$acc            tempo_diags%snow_liquid_equiv_precip, tempo_diags%graupel_liquid_equiv_precip, &
      !$acc            tempo_diags%frozen_fraction, tempo_diags%frz_rain_precip)
    endif
  end subroutine tempo_run_sync_host_diags


  subroutine tempo_run_sync_host_fields(use_temperature, t, th, qv, qc, qr, qi, qs, qg, ni, nr, nc, nwfa, nifa, ng, qb)
    logical, intent(in) :: use_temperature
    real(wp), dimension(:,:,:), intent(inout), optional :: t, th
    real(wp), dimension(:,:,:), intent(inout) :: qv, qc, qr, qi, qs, qg, ni, nr
    real(wp), dimension(:,:,:), intent(inout), optional :: nc, nwfa, nifa, ng, qb
    if (use_temperature) then
      if (present(t)) then
        !$acc update host(t)
      endif
    else
      if (present(th)) then
        !$acc update host(th)
      endif
    endif
    !$acc update host(qv, qc, qr, qi, qs, qg, ni, nr)
    if (present(nc)) then
      !$acc update host(nc)
    endif
    if (present(nwfa)) then
      !$acc update host(nwfa)
    endif
    if (present(nifa)) then
      !$acc update host(nifa)
    endif
    if (present(ng)) then
      !$acc update host(ng)
    endif
    if (present(qb)) then
      !$acc update host(qb)
    endif
  end subroutine tempo_run_sync_host_fields



  subroutine tempo_aerosol_surface_emissions(dt, nwfa, nwfa2d, ims, ime, jms, jme, kms, kme, kts)
    !! adds aerosol surface emissions to the 3D field
    real(wp), intent(in) :: dt
    real(wp), dimension(kms:kme, ims:ime, jms:jme), intent(inout) :: nwfa 
    real(wp), dimension(ims:ime, jms:jme), intent(in) :: nwfa2d
    integer, intent(in) :: ims, ime, jms, jme, kms, kme, kts
    integer :: i, j

    do j = jms, jme
      do i = ims, ime
        nwfa(kts,i,j) = nwfa(kts,i,j) + nwfa2d(i,j) * dt
      enddo
    enddo
  end subroutine tempo_aerosol_surface_emissions


  subroutine read_table_freezewater(filename, table_size)
    !! read lookup table for frozen cloud and rain water
    use module_mp_tempo_params, only : tpi_qrfz, tni_qrfz, &
      tpg_qrfz, tnr_qrfz, tpi_qcfz, tni_qcfz

    character(len=*), intent(in) :: filename
    integer, intent(in) :: table_size

    integer :: mp_unit, istat

    mp_unit = 11
    call check_before_table_read(filename, table_size)
    open(unit=mp_unit, file=filename, form='unformatted', status='old', access='stream', &
      action='read', iostat=istat, convert='big_endian')
    read(mp_unit) tpi_qrfz
    read(mp_unit) tni_qrfz
    read(mp_unit) tpg_qrfz
    read(mp_unit) tnr_qrfz
    read(mp_unit) tpi_qcfz
    read(mp_unit) tni_qcfz
    close(unit=mp_unit)
    !! push freezewater tables to device (host read; no managed memory)
    !$acc update device(tpi_qrfz, tni_qrfz, tpg_qrfz, tnr_qrfz, tpi_qcfz, tni_qcfz)
  end subroutine read_table_freezewater


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
    open(unit=mp_unit, file=filename, form='unformatted', status='old', access='stream', &
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
    !! push rain-snow collection tables to device (host read; no managed memory)
    !$acc update device(tcs_racs1, tmr_racs1, tcs_racs2, tmr_racs2, &
    !$acc               tcr_sacr1, tms_sacr1, tcr_sacr2, tms_sacr2, &
    !$acc               tnr_racs1, tnr_racs2, tnr_sacr1, tnr_sacr2)
  end subroutine read_table_qr_acr_qs


  subroutine read_table_qr_acr_qg(filename, table_size)
    !! read lookup table for rain-graupel collection
    use module_mp_tempo_params, only : tcg_racg, tmr_racg, &
      tcr_gacr, tnr_racg, tnr_gacr

    character(len=*), intent(in) :: filename
    integer, intent(in) :: table_size
    
    integer :: mp_unit, istat

    mp_unit = 11
    call check_before_table_read(filename, table_size)
    open(unit=mp_unit, file=filename, form='unformatted', status='old', access='stream', &
      action='read', iostat=istat, convert='big_endian')
    read(mp_unit) tcg_racg
    read(mp_unit) tmr_racg
    read(mp_unit) tcr_gacr
    read(mp_unit) tnr_racg
    read(mp_unit) tnr_gacr
    close(unit=mp_unit)
    !! push rain-graupel collection tables to device (host read; no managed memory)
    !$acc update device(tcg_racg, tmr_racg, tcr_gacr, tnr_racg, tnr_gacr)
  end subroutine read_table_qr_acr_qg


  subroutine read_table_ccn(filename, table_size)
    !! read static file containing CCN activation of aerosols;
    !! the data were created from a parcel model by Feingold and Heymsfield (1992)
    !! https://doi.org/10.1175/1520-0469(1992)049<2325:POCGOD>2.0.CO;2
    !! with further changes by Eidhammer and Kreidenweis
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
    !$acc update device(tnccn_act)
  end subroutine read_table_ccn


  subroutine check_before_table_read(filename, table_size)
    !! checks that lookup tables exist and are the correct size
    !! before attempting to read them

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


  subroutine init_ml_data()
    !! initialize machine learning data for tempo microphysics
    type(ty_tempo_ml_data) :: tempo_ml_data

    call nvtx_range_push('init_ml_data')

    ! cloud water
    tempo_ml_data%input_size = nc_ml_input
    tempo_ml_data%node_size = nc_ml_nodes
    tempo_ml_data%output_size = nc_ml_output

    if (.not.allocated(tempo_ml_data%transform_mean)) allocate(tempo_ml_data%transform_mean(nc_ml_input))
    if (.not.allocated(tempo_ml_data%transform_var)) allocate(tempo_ml_data%transform_var(nc_ml_input))

    tempo_ml_data%transform_mean = nc_ml_trans_mean
    tempo_ml_data%transform_var = nc_ml_trans_var

    if (.not.allocated(tempo_ml_data%weights00)) allocate(tempo_ml_data%weights00(nc_ml_nodes,nc_ml_input))
    if (.not.allocated(tempo_ml_data%weights01)) allocate(tempo_ml_data%weights01(nc_ml_output,nc_ml_nodes))
    if (.not.allocated(tempo_ml_data%bias00)) allocate(tempo_ml_data%bias00(nc_ml_nodes))
    if (.not.allocated(tempo_ml_data%bias01)) allocate(tempo_ml_data%bias01(nc_ml_output))

    tempo_ml_data%weights00 = reshape(nc_ml_w00, (/nc_ml_nodes, nc_ml_input/))
    tempo_ml_data%weights01 = reshape(nc_ml_w01, (/nc_ml_output, nc_ml_nodes/))
    tempo_ml_data%bias00 = nc_ml_b00
    tempo_ml_data%bias01 = nc_ml_b01

    ! save neural network
    call save_or_read_ml_data(ml_data_in=tempo_ml_data)
    call nvtx_range_pop()
  end subroutine init_ml_data

end module module_mp_tempo_driver
