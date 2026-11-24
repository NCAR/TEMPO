module module_mp_tempo_cfgs
  !! tempo configs

  use module_mp_tempo_params, only: wp

  implicit none
  private

  public :: ty_tempo_cfgs, ty_tempo_table_cfgs, ty_tempo_spp

!> \section arg_table_ty_tempo_cfgs Argument Table
!! \htmlinclude ty_tempo_cfgs.html
!!
  
  ! tempo configuration flags for init
  type :: ty_tempo_cfgs
    logical :: aerosolaware_flag = .true. !! flag to run aerosol-aware microphysics
    logical :: diagnostic_aerosols_flag = .false. !! flag to use diagnostic aerosols (no update) in the microphysics
    logical :: hailaware_flag = .true. !! flag to run hail-aware microphysics
    logical :: hailhyperaware_flag = .true. !! flag to run hail-hyperaware microphysics (splits graupel/hail and performs simple hail physics)
    logical :: ml_for_bl_nc_flag = .false. !! flag to run machine-learning prediction for subgrid cloud number concentration
    logical :: ml_for_nc_flag = .false. !! flag to run machine-learning prediction for tempo cloud number concentration
    logical :: semi_sedi_flag = .false. !! flag for semi-lagrangian sedimentation
    logical :: refl10cm_from_melting_flag = .false. !! flag to calculate reflectivity for melting snow and graupel
    logical :: turn_off_micro_flag = .false. !! flag to turn off all microphysical processes
    logical :: cloud_condensation_flag = .true. !! flag to control cloud condensation
    logical :: verbose = .false. !! flag to turn on verbose print statements
    logical :: check_tables = .false. !! flag to turn on check for lookup tables
    character(len=32) :: single_moment_nc_opt = 'land' !! option for single moment cloud number concentration when land input is not present (options include 'land' that uses the parameter nt_c_l, 'ocean' that uses nt_c_o, or string value in m^-3, e.g. '10.e6')
    ! flags to turn on/off diagnostic output
    logical :: refl10cm_flag = .true. !! flag to output 10cm reflectivity
    logical :: re_cloud_flag = .true. !! flag to output cloud effective radius
    logical :: re_ice_flag = .true. !! flag to output ice effective radius
    logical :: re_snow_flag = .true. !! flag to output snow effective radius
    logical :: max_hail_diameter_flag = .true. !! flag to output maximum hail diameter
    logical :: rain_med_vol_diam_flag = .false. !! flag to output median volume diameter for rain
    logical :: graupel_med_vol_diam_flag = .false. !! flag to output median volume diameter for graupel
    logical :: cloud_number_mixing_ratio_flag = .false. !! flag to output cloud number mixing ratio

    ! stochastically perturbed parameters (SPP); each column is perturbed by one value of a
    ! random pattern supplied by the host model (spp_pattern in tempo_run)
    logical :: spp_graupel_intercept_flag = .false. !! perturb the one-moment graupel intercept parameter (log10 N0 += pattern)
    logical :: spp_cloud_shape_flag = .false. !! perturb the cloud droplet gamma shape parameter (nu_c += nint(2 * pattern))
    logical :: spp_nucleation_flag = .false. !! perturb droplet activation and dust ice nucleation (aerosol-aware only)
    real(wp) :: spp_pattern_max = 0._wp !! largest magnitude of the SPP pattern (perturbation amplitude * standard deviation cutoff)

    contains
    procedure :: get_nc_val => resolve_nc_value
    procedure :: set_spp_options
    procedure :: spp_perturbations
  end type

  ! perturbations applied to one column, computed from the SPP pattern value for that column
  type :: ty_tempo_spp
    real(wp) :: graupel_log10_n0_offset = 0._wp !! added to log10 of the one-moment graupel intercept parameter
    integer :: cloud_shape_offset = 0 !! added to the cloud droplet gamma shape parameter nu_c (kept within 2-15)
    real(wp) :: nucleation_offset = 0._wp !! >= 0; added to w for droplet activation [m s-1] and scales dust ice nuclei by (1 + 50 * offset)
  end type

  ! tempo lookup table filenames
  type :: ty_tempo_table_cfgs
    character(len=100) :: ccn_table_name = 'ccn_activate.bin' !! ccn table name
    character(len=100) :: qrqg_table_name = 'qr_acr_qg_data_tempo_v3' !! rain-graupel collection table name
    character(len=100) :: qrqs_table_name = 'qr_acr_qs_data_tempo_v3' !! rain-snow collection table name
    character(len=100) :: freezewater_table_name = 'freeze_water_data_tempo_v3' !! freeze water collection table name
  end type

  contains

  real(wp) function resolve_nc_value(this, val_land, val_ocean)
    !! resolve the cloud number concentration value based on the input option

    class(ty_tempo_cfgs), intent(in) :: this
    real(wp), intent(in) :: val_land, val_ocean
    integer :: read_status

    select case (trim(adjustl(this%single_moment_nc_opt)))
    case ('land')
      resolve_nc_value = val_land
    case ('ocean')
      resolve_nc_value = val_ocean
    case default
      ! try to extract a real number from the string
      read(this%single_moment_nc_opt, *, iostat=read_status) resolve_nc_value

      ! catch typos
      if (read_status /= 0) then
        resolve_nc_value = val_land
        ! write(*,*) "invalid input: ", trim(this%single_moment_nc_opt)
        ! error stop "config error: single_moment_nc_opt must be 'land', 'ocean', or a number."
      end if
    end select
  end function resolve_nc_value


  subroutine set_spp_options(this, spp_mp, spp_amplitude, spp_stddev_cutoff)
    !! sets the SPP options from the bitmask used by host models (as in Thompson microphysics):
    !! 1 = graupel intercept, 2 = cloud droplet shape parameter, 4 = droplet activation and
    !! ice nucleation; values add, e.g. 3 = 1 + 2 and 7 = all
    !!
    !! spp_amplitude and spp_stddev_cutoff are the host's SPP perturbation magnitude and standard
    !! deviation cutoff for microphysics; their product bounds the pattern and keeps the
    !! nucleation perturbation non-negative

    class(ty_tempo_cfgs), intent(inout) :: this
    integer, intent(in) :: spp_mp
    real(wp), intent(in), optional :: spp_amplitude, spp_stddev_cutoff

    this%spp_graupel_intercept_flag = btest(spp_mp, 0)
    this%spp_cloud_shape_flag = btest(spp_mp, 1)
    this%spp_nucleation_flag = btest(spp_mp, 2)

    this%spp_pattern_max = 0._wp
    if (present(spp_amplitude) .and. present(spp_stddev_cutoff)) then
      this%spp_pattern_max = abs(spp_amplitude * spp_stddev_cutoff)
    endif
  end subroutine set_spp_options


  function spp_perturbations(this, pattern) result(spp)
    !! perturbations for one column from its SPP pattern value
    !!
    !! the graupel intercept offset is the pattern value; the cloud droplet shape offset is
    !! twice the pattern value rounded to an integer; the nucleation offset is shifted by
    !! spp_pattern_max so it is non-negative (0 to spp_pattern_max / 2 for a pattern within
    !! +/- spp_pattern_max)

    class(ty_tempo_cfgs), intent(in) :: this
    real(wp), intent(in) :: pattern
    type(ty_tempo_spp) :: spp

    if (this%spp_graupel_intercept_flag) spp%graupel_log10_n0_offset = pattern
    if (this%spp_cloud_shape_flag) spp%cloud_shape_offset = nint(2._wp * pattern)
    if (this%spp_nucleation_flag) spp%nucleation_offset = max(0._wp, 0.25_wp * (pattern + this%spp_pattern_max))
  end function spp_perturbations

end module module_mp_tempo_cfgs
