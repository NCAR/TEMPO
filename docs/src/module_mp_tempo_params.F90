module module_mp_tempo_params
  !! parameters and variables used in tempo microphysics

! define machine precision
#if defined(mpas)
  use mpas_kind_types, only: wp => RKIND, sp => R4KIND, dp => R8KIND
#elif defined(ccpp)
  use machine, only: wp => kind_phys, sp => kind_sngl_prec, dp => kind_dbl_prec
#else
  use machine, only: wp => kind_phys, sp => kind_sngl_prec, dp => kind_dbl_prec
#endif
  use iso_fortran_env, only : real32, real64 ! for machine-independent lookup table precisions
  implicit none

  public

  character(len=11) :: tempo_version !! tempo version string (max is xxx.xxx.xxx)

  ! tempo configuration flags for init
  type :: ty_tempo_init_cfgs
    logical :: aerosolaware_flag = .true. !! flag to run aerosol-aware microphysics
    logical :: hailaware_flag = .true. !! flag to run hail-aware microphysics
    logical :: restart_flag = .false. !! flag for restart or DA cycling
    character(len=4) :: model_flag !! flag for model
  end type

  ! tempo configuration flags for the driver
  type :: ty_tempo_driver_cfgs
    logical :: sedi_semi = .false. !! flag for semi-lagrangian sedimentation
  end type

  ! tempo lookup table filenames
  type :: ty_tempo_table_cfgs
    character(len=100) :: ccn_table_name = 'ccn_activate.bin' !! ccn table name
    character(len=100) :: qrqg_table_name = 'qr_acr_qg_data_tempo_v3' !! rain-graupel collection table name
    character(len=100) :: qrqs_table_name = 'qr_acr_qs_data_tempo_v3' !! rain-snow collection table name
    character(len=100) :: freezewater_table_name = 'freeze_water_data_tempo_v3' !! freeze water collection table name
  end type

  ! tempo configurations
  type(ty_tempo_init_cfgs) :: tempo_init_cfgs
  type(ty_tempo_driver_cfgs) :: tempo_driver_cfgs
  type(ty_tempo_table_cfgs) :: tempo_table_cfgs

  ! parameters that can be changed ------------------------------------------------------------------------
  integer, parameter :: idx_bg1 = 6 !! index from rho_g when hail_aware = false: density = 500 \(kg\, m^{-3}\)
  
  real(wp), parameter :: av_r = 4854._wp !! rain fallspeed power-law coefficient
  real(wp), parameter :: bv_r = 1.0_wp !! rain fallspeed power-law coefficient
    !! @note
    !! fallspeed power law relations are 
    !! \[ v =(a_{v}D^{b_{v}})\exp\left(-f_{v}D\right), \textrm{where}\,fv = 0\, \textrm{for graupel/ice} \]
    !! and coefficients are from from Ferrier (1994) for rain and 
    !! Thompson et al. (2008) for ice, snow, and graupel
    !! @endnote
  real(wp), parameter :: av_s = 40._wp !! snow fallspeed power-law coefficient
  real(wp), parameter :: bv_s = 0.55_wp !! snow fallspeed power-law coefficient
  real(wp), parameter :: fv_s = 100._wp !! snow fallspeed power-law coefficient
  real(wp), parameter :: bv_c = 2.0_wp !! cloud fallspeed power-law coefficient
  real(wp), parameter :: bv_i = 1.0_wp !! ice fallspeed power-law coefficient
  real(wp), parameter :: av_g_old = 442._wp !! graupel fallspeed power-law coefficient (hail_aware = false)
  real(wp), parameter :: bv_g_old = 0.89_wp !! graupel fallspeed power-law coefficient (hail_aware = false)

  real(wp), parameter :: am_s = 0.069_wp !! snow mass power-law coefficient
  real(wp), parameter :: bm_s = 2.0_wp !! snow mass power-law coefficient
    !! @note
    !! mass power law relations are
    !! \[ m = a_{m}D^{b_{m}} \]
    !! and coefficients for snow are from Field et al. (2005) and 
    !! others assume a spherical form
    !! @endnote
  real(wp), parameter :: bm_g = 3.0_wp !! graupel mass power-law coefficient
  real(wp), parameter :: bm_i = 3.0_wp !! ice mass power-law coefficient
  real(wp), parameter :: bm_r = 3.0_wp !! rain mass power-law coefficient

  real(wp), parameter :: rho_i = 890._wp !! density of cloud ice \([kg\, m^{-3}]\)
  real(wp), parameter :: xm0i = 1.e-12_wp !! ice initiates with this mass \([kg]\)
  real(wp), parameter :: d0c = 1.e-6_wp !! minimum diameter of cloud droplets \([m]\)
  real(wp), parameter :: d0r = 50.e-6_wp !! minimum diameter of raindrops \([m]\)
  real(wp), parameter :: d0s = 300.e-6_wp !! minimum diameter of snow \([m]\)
  real(wp), parameter :: d0g = 350.e-6_wp !! minimum diameter of graupel \([m]\)

  real(wp), parameter :: c_cube = 0.5_wp !! capacitance of a sphere \(\left(D^{3}\right)\)
  real(wp), parameter :: c_sqrd = 0.15_wp !! capacitance of plates/aggregates \(\left(D^{2}\right)\)

  real(wp), parameter :: mu_r = 0.0_wp !! shape parameter for rain
  real(wp), parameter :: mu_s = 0.6357_wp !! shape parameter for snow
  real(wp), parameter :: mu_g = 0.0_wp !! shape parameter for graupel
  real(wp), parameter :: mu_i = 0.0_wp !! shape parameter for cloud ice
    !! @note
    !! generalized gamma distributions for rain, graupel and cloud ice have the form
    !! \[ n\left(D\right) = n_{0} D^{\mu}\exp(-\lambda D) \]
    !! \[ \textrm{where}\\ \mu = 0 \\ \textrm{is exponential} \]
    !! @endnote
  real(wp), parameter :: nu_c_scale = 1000.e6_wp !! scaling parameter for nu_c
  integer, parameter :: nu_c_max = 15 !! maximum value for nu_c
  integer, parameter :: nu_c_min = 2 !! minimum value for nu_c

  ! parameters that should NOT be changed -----------------------------------------------------------------
  integer, parameter :: table_sp = real32 !! precision for lookup tables (machine independent)
  integer, parameter :: table_dp = real64 !! precision for lookup tables (machine independent)

  integer, parameter :: nrhg = 9 !! graupel density array size when hail_aware = true
  integer, parameter :: nrhg1 = 1 !! graupel density array size when hail_aware = false

  real(wp), parameter :: rho_w = 1000._wp !! density of liquid water \([kg\, m^{-3}]\)
  real(wp), dimension(nrhg), parameter :: rho_g = [50._wp, 100._wp, 200._wp, 300._wp, 400._wp, &
    500._wp, 600._wp, 700._wp, 800._wp] !! !! densities of graupel when hail_aware = true \([kg\, m^{-3}]\)

  real(wp), parameter :: sc = 0.632_wp !! schmidt number

  ! these can be overwritten by a host model and don't have a parameter attribute
  real(wp) :: pi = 3.1415926536_wp !! pi is approximately 355/113
  real(wp) :: lsub = 2.834e6_wp !! enthalpy of sublimation \([J\, kg^{-1}]\)
  real(wp) :: lvap0 = 2.5e6_wp !! enthalpy of vaporization \([J\, kg^{-1}]\)
  real(wp) :: rv = 461.5_wp !! gas constant for water vapor \([J\, K^{-1}\, kg^{-1}]\)
  real(wp) :: r = 287.04_wp !! gas constant for dry air \([J\, K^{-1}\, kg^{-1}]\)
  
  ! lookup table dimensions
  integer, parameter :: nbins = 100 !! lookup table dimension (number of bins)
  integer, parameter :: nbc = nbins !! lookup table dimension for cloud water
  integer, parameter :: nbr = nbins !! lookup table dimension for rain
  integer, parameter :: nbs = nbins !! lookup table dimension for snow
  integer, parameter :: nbi = nbins !! lookup table dimension
  integer, parameter :: nbg = nbins !! lookup table dimension
  integer, parameter :: ntb_i = 64 !! lookup table dimension for cloud ice
  integer, parameter :: ntb_i1 = 55 !! lookup table dimension for cloud ice
  integer, parameter :: ntb_c = 37 !! lookup table dimension for cloud water
  integer, parameter :: ntb_t = 9 !! lookup table dimension for temperature
  integer, parameter :: ntb_g1 = 37 !! lookup table dimension for graupel
  integer, parameter :: ntb_s = 37 !! lookup table dimension for snow
  integer, parameter :: ntb_g = 37 !! lookup table dimension for graupel
  integer, parameter :: ntb_r = 37 !! lookup table dimension for rain
  integer, parameter :: ntb_r1 = 37 !! lookup table dimension for rain
  integer, parameter :: ntb_t1 = 45 !! lookup table dimension for temperature
  integer, parameter :: ntb_in = 55 !! lookup table dimension for IN
  integer, parameter :: ntb_arc = 7 !! lookup table dimension for CCN activation
  integer, parameter :: ntb_arw = 9 !! lookup table dimension for CCN activation
  integer, parameter :: ntb_art = 7 !! lookup table dimension for CCN activation
  integer, parameter :: ntb_arr = 5 !! lookup table dimension for CCN activation
  integer, parameter :: ntb_ark = 4 !! lookup table dimension for CCN activation

  ! lookup table data
  real(wp), dimension(ntb_c), parameter :: &
    r_c = [1.e-6_wp,2.e-6_wp,3.e-6_wp,4.e-6_wp,5.e-6_wp,6.e-6_wp,7.e-6_wp,8.e-6_wp,9.e-6_wp, &
    1.e-5_wp,2.e-5_wp,3.e-5_wp,4.e-5_wp,5.e-5_wp,6.e-5_wp,7.e-5_wp,8.e-5_wp,9.e-5_wp, &
    1.e-4_wp,2.e-4_wp,3.e-4_wp,4.e-4_wp,5.e-4_wp,6.e-4_wp,7.e-4_wp,8.e-4_wp,9.e-4_wp, &
    1.e-3_wp,2.e-3_wp,3.e-3_wp,4.e-3_wp,5.e-3_wp,6.e-3_wp,7.e-3_wp,8.e-3_wp,9.e-3_wp, &
    1.e-2_wp] !! mass bins for cloud water \([kg\, m^{-3}]\)

  real(wp), dimension(ntb_i), parameter :: &
    r_i = [1.e-10_wp,2.e-10_wp,3.e-10_wp,4.e-10_wp, &
    5.e-10_wp,6.e-10_wp,7.e-10_wp,8.e-10_wp,9.e-10_wp, &
    1.e-9_wp,2.e-9_wp,3.e-9_wp,4.e-9_wp,5.e-9_wp,6.e-9_wp,7.e-9_wp,8.e-9_wp,9.e-9_wp, &
    1.e-8_wp,2.e-8_wp,3.e-8_wp,4.e-8_wp,5.e-8_wp,6.e-8_wp,7.e-8_wp,8.e-8_wp,9.e-8_wp, &
    1.e-7_wp,2.e-7_wp,3.e-7_wp,4.e-7_wp,5.e-7_wp,6.e-7_wp,7.e-7_wp,8.e-7_wp,9.e-7_wp, &
    1.e-6_wp,2.e-6_wp,3.e-6_wp,4.e-6_wp,5.e-6_wp,6.e-6_wp,7.e-6_wp,8.e-6_wp,9.e-6_wp, &
    1.e-5_wp,2.e-5_wp,3.e-5_wp,4.e-5_wp,5.e-5_wp,6.e-5_wp,7.e-5_wp,8.e-5_wp,9.e-5_wp, &
    1.e-4_wp,2.e-4_wp,3.e-4_wp,4.e-4_wp,5.e-4_wp,6.e-4_wp,7.e-4_wp,8.e-4_wp,9.e-4_wp, &
    1.e-3_wp] !! mass bins for ice water \([kg\, m^{-3}]\)

  real(wp), dimension(ntb_r), parameter :: &
    r_r = [1.e-6_wp,2.e-6_wp,3.e-6_wp,4.e-6_wp,5.e-6_wp,6.e-6_wp,7.e-6_wp,8.e-6_wp,9.e-6_wp, &
    1.e-5_wp,2.e-5_wp,3.e-5_wp,4.e-5_wp,5.e-5_wp,6.e-5_wp,7.e-5_wp,8.e-5_wp,9.e-5_wp, &
    1.e-4_wp,2.e-4_wp,3.e-4_wp,4.e-4_wp,5.e-4_wp,6.e-4_wp,7.e-4_wp,8.e-4_wp,9.e-4_wp, &
    1.e-3_wp,2.e-3_wp,3.e-3_wp,4.e-3_wp,5.e-3_wp,6.e-3_wp,7.e-3_wp,8.e-3_wp,9.e-3_wp, &
    1.e-2_wp] !! mass bins for rain \([kg\, m^{-3}]\)

  real(wp), dimension(ntb_s), parameter :: &
    r_s = [1.e-6_wp,2.e-6_wp,3.e-6_wp,4.e-6_wp,5.e-6_wp,6.e-6_wp,7.e-6_wp,8.e-6_wp,9.e-6_wp, &
    1.e-5_wp,2.e-5_wp,3.e-5_wp,4.e-5_wp,5.e-5_wp,6.e-5_wp,7.e-5_wp,8.e-5_wp,9.e-5_wp, &
    1.e-4_wp,2.e-4_wp,3.e-4_wp,4.e-4_wp,5.e-4_wp,6.e-4_wp,7.e-4_wp,8.e-4_wp,9.e-4_wp, &
    1.e-3_wp,2.e-3_wp,3.e-3_wp,4.e-3_wp,5.e-3_wp,6.e-3_wp,7.e-3_wp,8.e-3_wp,9.e-3_wp, &
    1.e-2_wp] !! mass bins for snow \([kg\, m^{-3}]\)

  real(wp), dimension(ntb_g), parameter :: &
    r_g = [1.e-6_wp,2.e-6_wp,3.e-6_wp,4.e-6_wp,5.e-6_wp,6.e-6_wp,7.e-6_wp,8.e-6_wp,9.e-6_wp, &
    1.e-5_wp,2.e-5_wp,3.e-5_wp,4.e-5_wp,5.e-5_wp,6.e-5_wp,7.e-5_wp,8.e-5_wp,9.e-5_wp, &
    1.e-4_wp,2.e-4_wp,3.e-4_wp,4.e-4_wp,5.e-4_wp,6.e-4_wp,7.e-4_wp,8.e-4_wp,9.e-4_wp, &
    1.e-3_wp,2.e-3_wp,3.e-3_wp,4.e-3_wp,5.e-3_wp,6.e-3_wp,7.e-3_wp,8.e-3_wp,9.e-3_wp, &
    1.e-2_wp] !! mass bins for graupel \([kg\, m^{-3}]\)

  real(wp), dimension(ntb_r1), parameter :: &
    n0r_exp = [1.e6_wp,2.e6_wp,3.e6_wp,4.e6_wp,5.e6_wp,6.e6_wp,7.e6_wp,8.e6_wp,9.e6_wp, &
    1.e7_wp,2.e7_wp,3.e7_wp,4.e7_wp,5.e7_wp,6.e7_wp,7.e7_wp,8.e7_wp,9.e7_wp, &
    1.e8_wp,2.e8_wp,3.e8_wp,4.e8_wp,5.e8_wp,6.e8_wp,7.e8_wp,8.e8_wp,9.e8_wp, &
    1.e9_wp,2.e9_wp,3.e9_wp,4.e9_wp,5.e9_wp,6.e9_wp,7.e9_wp,8.e9_wp,9.e9_wp, &
    1.e10_wp] !! y-intercept bins for rain \([m^{-4}]\)

  real(wp), dimension(ntb_g1), parameter :: &
    n0g_exp = [1.e2_wp,2.e2_wp,3.e2_wp,4.e2_wp,5.e2_wp,6.e2_wp,7.e2_wp,8.e2_wp,9.e2_wp, &
    1.e3_wp,2.e3_wp,3.e3_wp,4.e3_wp,5.e3_wp,6.e3_wp,7.e3_wp,8.e3_wp,9.e3_wp, &
    1.e4_wp,2.e4_wp,3.e4_wp,4.e4_wp,5.e4_wp,6.e4_wp,7.e4_wp,8.e4_wp,9.e4_wp, &
    1.e5_wp,2.e5_wp,3.e5_wp,4.e5_wp,5.e5_wp,6.e5_wp,7.e5_wp,8.e5_wp,9.e5_wp, &
    1.e6_wp] !! y-intercept bins for graupel \([m^{-4}]\)

  real(wp), dimension(ntb_i1), parameter :: &
    nt_i = [1.0_wp,2.0_wp,3.0_wp,4.0_wp,5.0_wp,6.0_wp,7.0_wp,8.0_wp,9.0_wp, &
    1.e1_wp,2.e1_wp,3.e1_wp,4.e1_wp,5.e1_wp,6.e1_wp,7.e1_wp,8.e1_wp,9.e1_wp, &
    1.e2_wp,2.e2_wp,3.e2_wp,4.e2_wp,5.e2_wp,6.e2_wp,7.e2_wp,8.e2_wp,9.e2_wp, &
    1.e3_wp,2.e3_wp,3.e3_wp,4.e3_wp,5.e3_wp,6.e3_wp,7.e3_wp,8.e3_wp,9.e3_wp, &
    1.e4_wp,2.e4_wp,3.e4_wp,4.e4_wp,5.e4_wp,6.e4_wp,7.e4_wp,8.e4_wp,9.e4_wp, &
    1.e5_wp,2.e5_wp,3.e5_wp,4.e5_wp,5.e5_wp,6.e5_wp,7.e5_wp,8.e5_wp,9.e5_wp, &
    1.e6_wp] !! number bins for ice \([m^{-3}]\)

  real(wp), dimension(ntb_in), parameter :: &
    nt_in = [1.0_wp,2.0_wp,3.0_wp,4.0_wp,5.0_wp,6.0_wp,7.0_wp,8.0_wp,9.0_wp, &
    1.e1_wp,2.e1_wp,3.e1_wp,4.e1_wp,5.e1_wp,6.e1_wp,7.e1_wp,8.e1_wp,9.e1_wp, &
    1.e2_wp,2.e2_wp,3.e2_wp,4.e2_wp,5.e2_wp,6.e2_wp,7.e2_wp,8.e2_wp,9.e2_wp, &
    1.e3_wp,2.e3_wp,3.e3_wp,4.e3_wp,5.e3_wp,6.e3_wp,7.e3_wp,8.e3_wp,9.e3_wp, &
    1.e4_wp,2.e4_wp,3.e4_wp,4.e4_wp,5.e4_wp,6.e4_wp,7.e4_wp,8.e4_wp,9.e4_wp, &
    1.e5_wp,2.e5_wp,3.e5_wp,4.e5_wp,5.e5_wp,6.e5_wp,7.e5_wp,8.e5_wp,9.e5_wp, &
    1.e6_wp] !! number bins for IN concentration from \(0.001-1000\, L^{-1}\) \([m^{-3}]\)
  
  ! variables ---------------------------------------------------------------------------------------------
  integer :: dim_nrhg !! number of dimensions for graupel density

  real(wp), dimension(nrhg) :: av_g = [45.9173813_wp, 67.0867386_wp, 98.0158463_wp, &
    122.353378_wp, 143.204224_wp, 161.794724_wp, &
    178.762115_wp, 194.488785_wp, 209.225876_wp] !! graupel fallspeed power-law coefficients (hail_aware = true)
  real(wp), dimension(nrhg) :: bv_g = [0.640961647_wp, 0.640961647_wp, 0.640961647_wp, &
    0.640961647_wp, 0.640961647_wp, 0.640961647_wp, &
    0.640961647_wp, 0.640961647_wp, 0.640961647_wp] !! graupel fallspeed power-law coefficients (hail_aware = true)
    !! @note
    !! av_g and bv_g values from A. Heymsfield: Best - Reynolds relationship
    !! @endnote

  real(wp) :: am_i !! ice mass-diameter power-law coefficient
  real(wp) :: am_r !! rain mass-diameter power-law coefficient
  real(wp), dimension (nrhg) :: am_g !! graupel mass-diameter power-law coefficient
  real(wp) :: lfus !! enthalpy of fusion \([J\, kg^{-1}]\)
  real(wp) :: olfus !! 1 / lfus \([kg\, J^{-1}]\)
  real(wp) :: orv !! 1 / rv \([K\, kg\, J^{-1}]\)

  real(wp) :: sc3 !! schmidt number to the 1/3 power

  real(wp) :: d0i !! minimum diameter of cloud ice \([m]\)
  real(wp) :: xm0s !! minimum mass of snow \([kg]\)
  real(wp) :: xm0g !! minimum mass of graupel \([kg]\)
  real(wp) :: obmi !! 1 / bm_i
  real(wp) :: obmr !! 1 / bm_r
  real(wp) :: oams !! 1 / am_s
  real(wp) :: obms !! 1 / bm_s
  real(wp) :: ocms !! oams ^ obms
  real(wp), dimension(nrhg) :: oamg !! 1 / am_g
  real(wp), dimension(nrhg) :: ocmg !! oamg ^ obmg
  real(wp) :: obmg !! 1 / bm_g

  ! various gamma calculations used throughout tempo
  real(wp), dimension(5,15) :: cce, ccg !! for \(ccg = \Gamma(x)\), cce is x for cloud water
  real(wp), dimension(15) :: ocg1, ocg2 !! inverse of specific ccg values
  real(wp), dimension(7) :: cie, cig !! for \(cig = \Gamma(x)\), cie is x for cloud ice
  real(wp) :: oig1, oig2 !! inverse of specific cig values
  real(wp), dimension(13) :: cre, crg !! for \(crg = \Gamma(x)\), cre is x for rain
  real(wp) :: ore1, org1, org2, org3 !! inverse of specific cre and crg values
  real(wp), dimension(17) :: cse, csg !! for \(csg = \Gamma(x)\), cse is x for snow
  real(wp), dimension(12,nrhg) :: cge, cgg !! for \(cgg = \Gamma(x)\), cge is x for graupel
  real(wp) :: oge1, ogg1, ogg2, ogg3 !! inverse of specific cge and cgg values

  ! precomputed constants in various rate equations
  real(wp) :: t1_qr_qc, t1_qr_qi, t2_qr_qi !! terms for rain collecting cloud water and cloud ice equations
  real(wp) :: t1_qs_qc, t1_qs_qi !! terms for snow collecting cloud water and cloud ice equations
  real(wp) :: t1_qr_ev, t2_qr_ev !! terms for rain evaporation equation
  real(wp) :: t1_qs_sd, t2_qs_sd !! terms for deposition/sublimation of snow equation
  real(wp) :: t1_qs_me, t2_qs_me !! terms for melting snow equation
  real(wp) :: t1_qg_sd !! term for deposition/sublimation of graupel equation
  real(wp) :: t1_qg_me !! term for melting graupel equation

  integer :: nic2, nii2, nii3, nir2, nir3, nis2, nig2, nig3, niin2 !! lookup table indexes
  integer :: nic1 !! used for cloud droplet number concentration lookup table
  !! @bug
  !! nic1 should be real(dp)
  !! @endbug

  real(dp), dimension(nbc) :: dc, dtc !! diameter and bin space for cloud water bins \([m]\)
  real(dp), dimension(nbi) :: di, dti !! diameter and bin space for ice bins \([m]\)
  real(dp), dimension(nbr) :: dr, dtr !! diameter and bin space for rain bins \([m]\)
  real(dp), dimension(nbs) :: ds, dts !! diameter and bin space for snow bins \([m]\)
  real(dp), dimension(nbg) :: dg, dtg !! diameter and bin space for graupel bins \([m]\)
  real(dp), dimension(nbc) :: t_nc !! cloud droplet number concentration bins \([cm^{-3}]\)

  ! lookup table data 
  real(dp), allocatable, dimension(:,:) :: t_efrw, t_efsw !! collection efficiency data arrays
  real(dp), allocatable, dimension(:,:,:) :: tpc_wev, tnc_wev !! evaporation data arrays
  real(table_sp), allocatable, dimension(:,:,:,:,:) :: tnccn_act !! cloud condensation nuclei data arrays
  real(table_dp), allocatable, dimension(:,:,:,:,:) :: tcg_racg, tmr_racg, tcr_gacr, &
    tnr_racg, tnr_gacr !! rain-graupel collection data arrays
  real(table_dp), allocatable, dimension(:,:,:,:) :: tcs_racs1, tmr_racs1, tcs_racs2, &
    tmr_racs2, tcr_sacr1, tms_sacr1, tcr_sacr2, tms_sacr2, &
    tnr_racs1, tnr_racs2, tnr_sacr1, tnr_sacr2 !! rain-snow collection data arrays
  real(table_dp), allocatable, dimension(:,:,:,:) :: tpi_qcfz, tni_qcfz !! cloud droplet freezing data arrays
  real(table_dp), allocatable, dimension(:,:,:,:) :: tpi_qrfz, tpg_qrfz, tni_qrfz, tnr_qrfz !! rain freezing data arrays
  real(table_dp), allocatable, dimension(:,:) :: tps_iaus, tni_iaus, tpi_ide !! cloud ice depositional growth and conversion to snow data arrays

    ! vars for the microphysics driver ------------------------------------------------------------------------

    real(wp)            :: RoverRv = 0.622
    real(wp)            :: Cp2 = 1004.0 ! AAJ change to Cp2

    ! ======================================================================
    ! Minimum microphys values
    ! R1 value, 1.e-12, cannot be set lower because of numerical
    ! problems with Paul Field's moments and should not be set larger
    ! because of truncation problems in snow/ice growth.
    real(wp), parameter :: R1 = 1.e-12
    real(wp), parameter :: R2 = 1.e-6
    real(wp), parameter :: eps = 1.e-15

    !    logical :: is_aerosol_aware = .true.
    logical :: merra2_aerosol_aware = .false.
    ! logical :: sedi_semi = .false.

    ! Hail-aware microphysics options
    integer :: dimNRHG

    ! Densities of rain, graupel, and cloud ice.
    real(wp), parameter :: rho_w2 = 1000.0 ! Change to rho_w2 to solve MPAS same var name conflict

    real(wp) :: t1_qg_qc, t2_qg_sd, t2_qg_me

    !=================================================================================================================
    ! Parameters needed by the microphysics driver

    ! Prescribed number of cloud droplets.  Set according to known data or
    ! roughly 100 per cc (100.e6 m^-3) for Maritime cases and
    ! 300 per cc (300.e6 m^-3) for Continental.  Gamma shape parameter,
    ! mu_c, calculated based on Nt_c is important in autoconversion
    ! scheme.  In 2-moment cloud water, Nt_c represents a maximum of
    ! droplet concentration and nu_c is also variable depending on local
    ! droplet number concentration.
    real(wp), parameter :: Nt_c_o = 50.e6
    real(wp), parameter :: Nt_c_l = 100.e6
    real(wp), parameter :: Nt_c_max = 1999.e6
    real(wp) :: Nt_c, mu_c
    real(wp) :: mu_c_o, mu_c_l

    real(wp) :: min_qv = 1.e-10
#if defined(ccpp_default)
    real(wp), parameter :: demott_nuc_ssati = 0.15 ! 0.15 for CCPP
#else
    real(wp), parameter :: demott_nuc_ssati = 0.25
#endif
    ! Declaration of constants for assumed CCN/IN aerosols when none in
    ! the input data.  Look inside the init routine for modifications
    ! due to surface land-sea points or vegetation characteristics.
    real(wp), parameter :: nwfa_default = 11.1e6
    real(wp), parameter :: naIN1 = 0.5e6
    real(wp), parameter :: nifa_default = naIN1*0.01
    real(wp), parameter :: aero_max = 9999.e6
    real(dp), parameter :: max_ni = 4999.e3
    real(wp), parameter :: icenuc_max = 1000.e3
    
    real(wp), parameter :: rime_threshold = 2.0 ! For MPAS
    real(wp), parameter :: rime_conversion = 0.95 ! For MPAS

    real(wp), parameter :: fv_r = 195.0
    real(wp) :: rho_s2 = 100.0 ! AAJ change to rho_s2 to solve MPAS same var name conflict
    real(wp), parameter :: av_c = 0.316946e8

    logical, parameter :: iiwarm = .false.
    logical, parameter :: dustyIce = .true.
    logical, parameter :: homogIce = .true.

    integer, parameter :: IFDRY = 0
    real(wp)           :: T_0 = 273.15

    real(wp), parameter :: naIN0 = 1.5e6
    real(wp), parameter :: naCCN0 = 300.0e6
    real(wp), parameter :: naCCN1 = 50.0e6

    ! Sum of two gamma distrib for snow (Field et al. 2005).
    ! N(D) = M2**4/M3**3 * [Kap0*exp(-M2*Lam0*D/M3)
    !      + Kap1*(M2/M3)**mu_s * D**mu_s * exp(-M2*Lam1*D/M3)]
    ! M2 and M3 are the (bm_s)th and (bm_s+1)th moments respectively
    ! calculated as function of ice water content and temperature.
    real(wp), parameter :: Kap0 = 490.6
    real(wp), parameter :: Kap1 = 17.46
    real(wp), parameter :: Lam0 = 20.78
    real(wp), parameter :: Lam1 = 3.29

    ! Y-intercept parameter for graupel is not constant and depends on
    ! mixing ratio.  Also, when mu_g is non-zero, these become equiv
    ! y-intercept for an exponential distrib and proper values are
    ! computed based on same mixing ratio and total number concentration.
    real(dp), parameter :: gonv_min = 1.e2
    real(dp), parameter :: gonv_max = 1.e6

    real(wp), parameter :: a_coeff = 0.47244157
    real(wp), parameter :: b_coeff = 0.54698726

#if defined(ccpp_default)
    real(wp) :: av_i
#else
    real(wp), parameter :: av_i = 1493.9
#endif

    ! Collection efficiencies.  Rain/snow/graupel collection of cloud
    ! droplets use variables (Ef_rw, Ef_sw, Ef_gw respectively) and
    ! get computed elsewhere because they are dependent on stokes
    ! number.
    real(wp), parameter :: Ef_si = 0.05
    real(wp), parameter :: Ef_rs = 0.95
    real(wp), parameter :: Ef_rg = 0.75
    real(wp), parameter :: Ef_ri = 0.95

    ! Constants in Cooper curve relation for cloud ice number.
    real(wp), parameter :: TNO = 5.0
    real(wp), parameter :: ATO = 0.304

    ! Rho_not used in fallspeed relations (rho_not/rho)**.5 adjustment.
    real(wp), parameter :: rho_not = 101325.0 / (287.05*298.0)

    ! Homogeneous freezing temperature
    real(wp), parameter :: HGFR = 235.16

    real(wp)            :: R_uni = 8.314                           ! J (mol K)-1

    real(dp)            :: k_b = 1.38065e-23                ! Boltzmann constant [J/K]
    real(dp)            :: M_w = 18.01528e-3                ! molecular mass of water [kg/mol]
    real(dp)            :: M_a = 28.96e-3                   ! molecular mass of air [kg/mol]
    real(dp)            :: N_avo = 6.022e23                 ! Avogadro number [1/mol]
    real(dp)            :: ma_w                             ! mass of water molecule [kg] (= M_w / N_avo, set in mp_tempo_params_init)
    real(wp)            :: ar_volume                        ! assume radius of 0.025 micrometer, 2.5e-6 cm (= 4.0 / 3.0 * PI * (2.5e-6)**3, set in mp_tempo_params_init)

    ! Aerosol table parameter: Number of available aerosols, vertical
    ! velocity, temperature, aerosol mean radius, and hygroscopicity.
    real(wp), dimension(ntb_arc), parameter :: &
        ta_Na = (/10.0, 31.6, 100.0, 316.0, 1000.0, 3160.0, 10000.0/)
    real(wp), dimension(ntb_arw), parameter :: &
        ta_Ww = (/0.01, 0.0316, 0.1, 0.316, 1.0, 3.16, 10.0, 31.6, 100.0/)
    real(wp), dimension(ntb_art), parameter :: &
        ta_Tk = (/243.15, 253.15, 263.15, 273.15, 283.15, 293.15, 303.15/)
    real(wp), dimension(ntb_arr), parameter :: &
        ta_Ra = (/0.01, 0.02, 0.04, 0.08, 0.16/)
    real(wp), dimension(ntb_ark), parameter :: &
        ta_Ka = (/0.2, 0.4, 0.6, 0.8/)

    ! For snow moments conversions (from Field et al. 2005)
    real(wp), dimension(10), parameter :: &
        sa = (/ 5.065339, -0.062659, -3.032362, 0.029469, -0.000285, &
        0.31255, 0.000204, 0.003199, 0.0, -0.015952/)
    real(wp), dimension(10), parameter :: &
        sb = (/ 0.476221, -0.015896, 0.165977, 0.007468, -0.000141, &
        0.060366, 0.000079, 0.000594, 0.0, -0.003577/)

    ! Temperatures (5 C interval 0 to -40) used in lookup tables.
    real(wp), dimension(ntb_t), parameter :: &
        Tc = (/-0.01, -5., -10., -15., -20., -25., -30., -35., -40./)

#if defined(ccpp_default)
    ! To permit possible creation of new lookup tables as variables expand/change,
    ! specify a name of external file(s) including version number for pre-computed
    ! Thompson tables.
    character(len=*), parameter :: thomp_table_file = 'thompson_tables_precomp_v2.sl'
    character(len=*), parameter :: qr_acr_qg_file = 'MP_TEMPO_QRacrQG.dat'
    character(len=*), parameter :: qr_acr_qg_hailaware_file = 'MP_TEMPO_HAILAWARE_QRacrQG.dat'
    character(len=*), parameter :: qr_acr_qs_file = 'MP_TEMPO_QRacrQS.dat'
    character(len=*), parameter :: freeze_h2o_file = 'MP_TEMPO_freezeH2O.dat'

    ! Min and max radiative effective radius of cloud water, cloud ice, and snow;
    ! performed by subroutine calc_effectRad. On purpose, these should stay PUBLIC.
    real(wp), parameter :: re_qc_min = 2.50e-6               ! 2.5 microns
    real(wp), parameter :: re_qc_max = 50.0e-6               ! 50 microns
    real(wp), parameter :: re_qi_min = 2.50e-6               ! 2.5 microns
    real(wp), parameter :: re_qi_max = 125.0e-6              ! 125 microns
    real(wp), parameter :: re_qs_min = 5.00e-6               ! 5 microns
    real(wp), parameter :: re_qs_max = 999.0e-6              ! 999 microns (1 mm)

    ! MPI communicator
    ! type(MPI_Comm) :: mpi_communicator

    ! Write tables with master MPI task after computing them in tempo_init
    logical :: thompson_table_writer
#endif

    ! ML data
    integer, parameter :: nc_ml_input = 7
    integer, parameter :: nc_ml_nodes = 24
    integer, parameter :: nc_ml_output = 1

    integer, parameter :: nr_ml_input = 7
    integer, parameter :: nr_ml_nodes = 24
    integer, parameter :: nr_ml_output = 1

    real(wp), dimension(nc_ml_input), parameter :: &
         nc_ml_trans_mean = (/0.000191556196486247, 3.58145042772654e-05, &
         3.12611085359273e-07, 5.74303078738579e-05, 84191.2092225319, &
         279.070551773565, 0.123679354084004/)
    real(wp), dimension(nc_ml_input), parameter :: &
         nc_ml_trans_var = (/5.78143564171777e-08, 3.22834309750552e-08, &
         6.45745893307455e-11, 4.16625579383794e-08, 215694631.771185, &
         94.6576255386858, 0.384841247662964/)

    real(wp), dimension(nc_ml_input * nc_ml_nodes), parameter :: &
         nc_ml_w00 = (/-2.006957, -0.2812008, -0.339073, 1.596426, 2.395225, 1.76315, &
         -0.0626798, 1.267002, -0.02234177, -6.522605e-33, 0.4792154, -0.1253034, &
         -3.217191, -3.092887, -0.0863651, 1.071625, 0.09741028, 0.2255831, &
         -0.6929023, -0.02693799, -3.432344e-33, -0.8791879, -0.9359049, 1.083484, &
         -0.07909214, -0.0122418, -0.02815927, 0.1676407, 0.08252326, 0.6697816, &
         -0.4019359, 0.4687141, 0.001813132, -9.792186e-33, 0.0409322, 0.0113192, &
         -0.01354596, 0.00307771, -0.4635534, 0.03835761, -0.1015553, 0.7316446, &
         -0.05791711, -0.0002690362, -7.920147e-33, -0.1216918, -0.3190572, 0.09809405, &
         -0.16476, -0.03387314, 0.005422261, 0.04043967, 0.03901243, 0.07444729, &
         0.01954299, 0.06918761, 0.04823543, -8.637957e-33, 0.06371575, -0.09250915, &
         -1.109653, -1.373999, -0.2412623, -0.04482195, 0.1584691, 0.06353725, &
         0.0006248798, 0.04593191, -8.878673e-33, -0.4988684, 0.01110262, 0.04623203, &
         0.006581791, 0.03536217, -0.1890567, -0.08839592, 0.1327181, 0.03478973, &
         -0.1565902, -0.100401, -0.1179777, -8.879818e-33, -0.1383738, 0.02847495, &
         -0.005902881, 0.005615512, -0.6308192, -0.02431803, -0.141971, -0.3490018, &
         -0.9850957, -0.1449479, -8.059166e-33, -0.1186465, -1.165381, 0.069015, &
         0.003388841, -0.04041302, 0.1638467, 0.1147008, -0.04833491, -0.07755993, &
         -0.5137688, 0.04546477, 0.04101883, 6.752353e-33, 0.3541977, 0.04880851, &
         -0.00102834, -0.01280629, -0.1116254, -0.02204754, 0.07100908, 0.2354002, &
         0.07129629, 0.2489657, 8.080785e-33, -0.03449865, 0.06037927, -0.02023619, &
         0.7779589, 0.04680278, 0.7492616, 0.6545208, -1.09497, -1.176524, &
         -0.451585, 0.881124, -0.4551499, -8.708624e-33, 1.006558, -0.04979523, &
         -0.0006915367, -0.002993054, 0.01654614, -0.07141764, -0.2216591, 0.8637336, &
         0.8358089, -0.7576646, 8.186339e-33, 0.1161914, 0.7121871, -1.146734, &
         0.03925339, 0.9975697, -0.06953461, -0.07598846, -0.06418022, -0.01897495, &
         -0.03612464, -0.06703389, -0.103049, 9.364858e-33, -0.03949085, 1.005938, &
         0.0001625419, -0.01027657, 0.03823901, 0.3197246, -0.1160718, 0.04449431, &
         0.009831783, 0.6551342, -8.470687e-33, 0.1374123, -0.01782138, 0.01719949/)
    real(wp), dimension(nc_ml_nodes), parameter :: &
         nc_ml_w01 = (/-4.045869, 0.558111, -1.567351, -1.64972, 2.748608, &
         -1.909901, 0.3955558, 1.507247, 0.3599722, -0.0001173223, -0.9398569, &
         -0.6867028, -61.72853, -63.13766, 1.165811, -0.6848684, 0.1931683, &
         1.1208, 2.63087, 0.740169, -21.62499, 1.545568, 3.575141, -1.299604/)
    real(wp), dimension(nc_ml_nodes), parameter :: &
         nc_ml_b00 = (/-0.9842531, 0.3064759, -0.4500185, 1.28336, 1.384105, 0.528031, &
         0.8453538, 1.579872, 2.245679, -0.008679952, 0.4549862, -0.136581, &
         -2.576741, -2.483647, -0.2089484, 0.7607977, 1.847745, 0.7316047, &
         -0.287945, 2.227298, -2.314714, -0.2561245, -0.6993448, -0.1359731/)
    real(wp), dimension(nc_ml_output), parameter :: &
         nc_ml_b01 = (/1.572826/)
    
  ! -------------------------------------------------------------------------------------------------------
  ! -------------------------------------------------------------------------------------------------------
  contains

  subroutine initialize_graupel_vars(hail_flag)
    !! initialize graupel variables based on hail-aware configuration flag

    logical, intent(in) :: hail_flag

    if (hail_flag) then
      dim_nrhg = nrhg
    else
      av_g(idx_bg1) = av_g_old
      bv_g(idx_bg1) = bv_g_old
      dim_nrhg = nrhg1
    endif
  end subroutine initialize_graupel_vars


  subroutine initialize_parameters()
    !! initialize tempo parameters and variables
    
    integer :: m, n

    ! pi could be set by a host model, thus these parameters need to be calculated here
    am_i = pi * rho_i / 6.0_wp
    am_r = pi * rho_w / 6.0_wp
    am_g = [pi * rho_g(1) / 6.0_wp, &
      pi * rho_g(2) / 6.0_wp, &
      pi * rho_g(3) / 6.0_wp, &
      pi * rho_g(4) / 6.0_wp, &
      pi * rho_g(5) / 6.0_wp, &
      pi * rho_g(6) / 6.0_wp, &
      pi * rho_g(7) / 6.0_wp, &
      pi * rho_g(8) / 6.0_wp, &
      pi * rho_g(9) / 6.0_wp]
#if defined(ccpp_default)
    av_i = av_s * d0s ** (bv_s - bv_i)
#endif
    ma_w = m_w / n_avo  
    ar_volume = 4.0_wp / 3.0_wp * pi * (2.5e-6_wp)**3

    lfus = lsub - lvap0
    olfus = 1.0_wp / lfus
    orv = 1.0_wp / rv

    ! Schmidt number to one-third used numerous times
    sc3 = sc**(1.0_wp/3.0_wp)

    ! compute minimum ice diameter from mass and minimum snow/graupel mass from diameter
    d0i = (xm0i/am_i)**(1.0_wp/bm_i)

    ! pre-compute various constants used in the microphysics equations
    xm0s = am_s * d0s**bm_s
    xm0g = am_g(nrhg) * d0g**bm_g
    obmi = 1.0_wp / bm_i
    obmr = 1.0_wp / bm_r
    oams = 1.0_wp / am_s
    obms = 1.0_wp / bm_s
    ocms = oams**obms
    obmg = 1.0_wp / bm_g
    do m = 1, nrhg
      oamg(m) = 1.0_wp / am_g(m)
      ocmg(m) = oamg(m)**obmg
    enddo

    ! gamma functions for cloud water
    do n = 1, 15
      cce(1,n) = n + 1._wp
      cce(2,n) = bm_r + n + 1._wp
      cce(3,n) = bm_r + n + 4._wp
      cce(4,n) = n + bv_c + 1._wp
      cce(5,n) = bm_r + n + bv_c + 1._wp
      ccg(1,n) = gamma(cce(1,n))
      ccg(2,n) = gamma(cce(2,n))
      ccg(3,n) = gamma(cce(3,n))
      ccg(4,n) = gamma(cce(4,n))
      ccg(5,n) = gamma(cce(5,n))
      ocg1(n) = 1.0_wp / ccg(1,n)
      ocg2(n) = 1.0_wp / ccg(2,n)
    enddo

    ! gamma functions for cloud ice
    cie(1) = mu_i + 1._wp
    cie(2) = bm_i + mu_i + 1._wp
    cie(3) = bm_i + mu_i + bv_i + 1._wp
    cie(4) = mu_i + bv_i + 1._wp
    cie(5) = mu_i + 2._wp
    cie(6) = bm_i*0.5_wp + mu_i + bv_i + 1._wp
    cie(7) = bm_i*0.5_wp + mu_i + 1._wp
    cig(1) = gamma(cie(1))
    cig(2) = gamma(cie(2))
    cig(3) = gamma(cie(3))
    cig(4) = gamma(cie(4))
    cig(5) = gamma(cie(5))
    cig(6) = gamma(cie(6))
    cig(7) = gamma(cie(7))
    oig1 = 1.0_wp / cig(1)
    oig2 = 1.0_wp / cig(2)

    ! gamma functions for rain
    cre(1) = bm_r + 1._wp
    cre(2) = mu_r + 1._wp
    cre(3) = bm_r + mu_r + 1._wp
    cre(4) = bm_r*2._wp + mu_r + 1._wp
    cre(5) = mu_r + bv_r + 1._wp
    cre(6) = bm_r + mu_r + bv_r + 1._wp
    cre(7) = bm_r*0.5_wp + mu_r + bv_r + 1._wp
    cre(8) = bm_r + mu_r + bv_r + 3._wp
    cre(9) = mu_r + bv_r + 3._wp
    cre(10) = mu_r + 2._wp
    cre(11) = 0.5_wp*(bv_r + 5._wp + 2._wp*mu_r)
    cre(12) = bm_r*0.5_wp + mu_r + 1._wp
    cre(13) = bm_r*2._wp + mu_r + bv_r + 1._wp

    do n = 1, 13
      crg(n) = gamma(cre(n))
    enddo

    ore1 = 1.0_wp / cre(1)
    org1 = 1.0_wp / crg(1)
    org2 = 1.0_wp / crg(2)
    org3 = 1.0_wp / crg(3)

    ! gamma functions for snow
    cse(1) = bm_s + 1._wp
    cse(2) = bm_s + 2._wp
    cse(3) = bm_s*2._wp
    cse(4) = bm_s + bv_s + 1._wp
    cse(5) = bm_s*2._wp + bv_s + 1._wp
    cse(6) = bm_s*2._wp + 1._wp
    cse(7) = bm_s + mu_s + 1._wp
    cse(8) = bm_s + mu_s + 2._wp
    cse(9) = bm_s + mu_s + 3._wp
    cse(10) = bm_s + mu_s + bv_s + 1._wp
    cse(11) = bm_s*2._wp + mu_s + bv_s + 1._wp
    cse(12) = bm_s*2._wp + mu_s + 1._wp
    cse(13) = bv_s + 2._wp
    cse(14) = bm_s + bv_s
    cse(15) = mu_s + 1._wp
    cse(16) = 1.0_wp + (1.0_wp + bv_s)/2._wp
    cse(17) = bm_s + bv_s + 2._wp

    do n = 1, 17
      csg(n) = gamma(cse(n))
    enddo

    ! gamma functions for graupel
    cge(1,:) = bm_g + 1._wp
    cge(2,:) = mu_g + 1._wp
    cge(3,:) = bm_g + mu_g + 1._wp
    cge(4,:) = bm_g*2. + mu_g + 1._wp
    cge(10,:) = mu_g + 2._wp
    cge(12,:) = bm_g*0.5_wp + mu_g + 1._wp

    do m = 1, nrhg
      cge(5,m) = bm_g*2._wp + mu_g + bv_g(m) + 1._wp
      cge(6,m) = bm_g + mu_g + bv_g(m) + 1._wp
      cge(7,m) = bm_g*0.5_wp + mu_g + bv_g(m) + 1._wp
      cge(8,m) = mu_g + bv_g(m) + 1._wp
      cge(9,m) = mu_g + bv_g(m) + 3._wp
      cge(11,m) = 0.5_wp*(bv_g(m) + 5._wp + 2._wp*mu_g)
    enddo

    do m = 1, nrhg
      do n = 1, 12
        cgg(n,m) = gamma(cge(n,m))
      enddo
    enddo
    oge1 = 1.0_wp / cge(1,1)
    ogg1 = 1.0_wp / cgg(1,1)
    ogg2 = 1.0_wp / cgg(2,1)
    ogg3 = 1.0_wp / cgg(3,1)

    ! rain collecting cloud water and cloud ice
    t1_qr_qc = pi * 0.25_wp * av_r * crg(9)
    t1_qr_qi = pi * 0.25_wp * av_r * crg(9)
    t2_qr_qi = pi * 0.25_wp * am_r*av_r * crg(8)

    ! snow collecting cloud water and cloud ice
    t1_qs_qc = pi * 0.25_wp * av_s
    t1_qs_qi = pi * 0.25_wp * av_s

    ! evaporation of rain; ignore depositional growth of rain.
    t1_qr_ev = 0.78_wp * crg(10)
    t2_qr_ev = 0.308_wp * sc3 * sqrt(av_r) * crg(11)

    ! sublimation/depositional growth of snow
    t1_qs_sd = 0.86_wp
    t2_qs_sd = 0.28_wp * sc3 * sqrt(av_s)

    ! melting of snow
    t1_qs_me = pi * 4._wp *c_sqrd * olfus * 0.86_wp
    t2_qs_me = pi * 4._wp *c_sqrd * olfus * 0.28_wp * sc3 * sqrt(av_s)

    ! sublimation/depositional growth of graupel
    t1_qg_sd = 0.86_wp * cgg(10,1)

    ! melting of graupel
    t1_qg_me = pi * 4._wp * c_cube * olfus * 0.86_wp * cgg(10,1)
  end subroutine initialize_parameters
    

  subroutine initialize_bins_for_tables()
    !! initialize log-spaced bins of hydrometer quantities used for lookup tables

    integer :: n

    ! constants for helping find lookup table indexes.
    nic2 = nint(log10(r_c(1)))
    nii2 = nint(log10(r_i(1)))
    nir2 = nint(log10(r_r(1)))
    nis2 = nint(log10(r_s(1)))
    nig2 = nint(log10(r_g(1)))
    nii3 = nint(log10(nt_i(1)))
    nir3 = nint(log10(n0r_exp(1)))
    nig3 = nint(log10(n0g_exp(1)))
    niin2 = nint(log10(nt_in(1)))

    ! bins of cloud water (from minimum diameter to 100 microns).
    dc(1) = real(d0c, kind=dp)
    dtc(1) = real(d0c, kind=dp)
    do n = 2, nbc
      dc(n) = dc(n-1) + 1.0e-6_dp
      dtc(n) = (dc(n) - dc(n-1))
    enddo

    ! bins of cloud ice (from min diameter up to 2x min snow size).
    call create_bins(numbins=nbi, lowbin=real(d0i, kind=dp), &
      highbin=2.0_dp*d0s, bins=di, deltabins=dti)

    ! bins of rain (from min diameter up to 5 mm).
    call create_bins(numbins=nbr, lowbin=real(d0r, kind=dp), &
      highbin=0.005_dp, bins=dr, deltabins=dtr)

    ! bins of snow (from min diameter up to 2 cm).
    call create_bins(numbins=nbs, lowbin=real(d0s, kind=dp), &
      highbin=0.02_dp, bins=ds, deltabins=dts)

    ! bins of graupel (from min diameter up to 5 cm).
    call create_bins(numbins=nbg, lowbin=real(d0g, kind=dp), &
      highbin=0.05_dp, bins=dg, deltabins=dtg)

    ! bins of cloud droplet number concentration (1 to 3000 per cc).
    call create_bins(numbins=nbc, lowbin=1.0_dp, &
      highbin=3000.0_dp, bins=t_nc)
    t_nc = t_nc * 1.0e6_dp
    nic1 = real(log(t_nc(nbc)/t_nc(1)), kind=dp)
  end subroutine initialize_bins_for_tables


 subroutine create_bins(numbins, lowbin, highbin, bins, deltabins)
    !! calculates log-spaced bins of hydrometer sizes to simplify calculations later
  
    integer, intent(in) :: numbins
    real(dp), intent(in) :: lowbin, highbin

    real(dp), dimension(:), intent(out) :: bins
    real(dp), dimension(:), intent(out), optional :: deltabins
    
    integer :: n
    real(dp), dimension(numbins+1) :: xdx
  
    xdx(1) = lowbin
    xdx(numbins+1) = highbin
    do  n = 2, numbins
      xdx(n) = exp(real(n-1, kind=dp)/real(numbins, kind=dp) * log(xdx(numbins+1)/xdx(1)) + log(xdx(1)))
    enddo

    do n = 1, numbins
      bins(n) = sqrt(xdx(n)*xdx(n+1))
    enddo

    if (present(deltabins)) then
      do n = 1, numbins
        deltabins(n) = xdx(n+1) - xdx(n)
      enddo
    endif
  end subroutine create_bins


  subroutine initialize_arrays_freezewater(table_size)
    !! initialize data arrays for Bigg (1953) freezing of cloud water and rain

    integer, intent(out), optional :: table_size

    ! Cloud water freezing
    if (.not. allocated(tpi_qcfz)) allocate(tpi_qcfz(ntb_c,nbc,ntb_t1,ntb_in), source=0._table_dp)
    if (.not. allocated(tni_qcfz)) allocate(tni_qcfz(ntb_c,nbc,ntb_t1,ntb_in), source=0._table_dp)

    ! Rain freezing
    if (.not. allocated(tpi_qrfz)) allocate(tpi_qrfz(ntb_r,ntb_r1,ntb_t1,ntb_in), source=0._table_dp)
    if (.not. allocated(tpg_qrfz)) allocate(tpg_qrfz(ntb_r,ntb_r1,ntb_t1,ntb_in), source=0._table_dp)
    if (.not. allocated(tni_qrfz)) allocate(tni_qrfz(ntb_r,ntb_r1,ntb_t1,ntb_in), source=0._table_dp)
    if (.not. allocated(tnr_qrfz)) allocate(tnr_qrfz(ntb_r,ntb_r1,ntb_t1,ntb_in), source=0._table_dp)

    ! Table size is precision * entries * dimensions
    if (present(table_size)) then
      table_size = (table_dp * 2 * (ntb_c*nbc*ntb_t1*ntb_in)) + &
        (table_dp * 4 * (ntb_r*ntb_r1*ntb_t1*ntb_in))
    endif 
  end subroutine initialize_arrays_freezewater


  subroutine initialize_arrays_qr_acr_qs(table_size)
    !! initialize data arrays for rain-snow collection

    integer, intent(out), optional :: table_size

    if (.not. allocated(tcs_racs1)) allocate(tcs_racs1(ntb_s,ntb_t,ntb_r1,ntb_r), source=0._table_dp)
    if (.not. allocated(tmr_racs1)) allocate(tmr_racs1(ntb_s,ntb_t,ntb_r1,ntb_r), source=0._table_dp)
    if (.not. allocated(tcs_racs2)) allocate(tcs_racs2(ntb_s,ntb_t,ntb_r1,ntb_r), source=0._table_dp)
    if (.not. allocated(tmr_racs2)) allocate(tmr_racs2(ntb_s,ntb_t,ntb_r1,ntb_r), source=0._table_dp)
    if (.not. allocated(tcr_sacr1)) allocate(tcr_sacr1(ntb_s,ntb_t,ntb_r1,ntb_r), source=0._table_dp)
    if (.not. allocated(tms_sacr1)) allocate(tms_sacr1(ntb_s,ntb_t,ntb_r1,ntb_r), source=0._table_dp)
    if (.not. allocated(tcr_sacr2)) allocate(tcr_sacr2(ntb_s,ntb_t,ntb_r1,ntb_r), source=0._table_dp)
    if (.not. allocated(tms_sacr2)) allocate(tms_sacr2(ntb_s,ntb_t,ntb_r1,ntb_r), source=0._table_dp)
    if (.not. allocated(tnr_racs1)) allocate(tnr_racs1(ntb_s,ntb_t,ntb_r1,ntb_r), source=0._table_dp)
    if (.not. allocated(tnr_racs2)) allocate(tnr_racs2(ntb_s,ntb_t,ntb_r1,ntb_r), source=0._table_dp)
    if (.not. allocated(tnr_sacr1)) allocate(tnr_sacr1(ntb_s,ntb_t,ntb_r1,ntb_r), source=0._table_dp)
    if (.not. allocated(tnr_sacr2)) allocate(tnr_sacr2(ntb_s,ntb_t,ntb_r1,ntb_r), source=0._table_dp)

    ! Table size is precision * entries * dimensions
    if (present(table_size)) then
      table_size = table_dp * 12 * (ntb_s*ntb_t*ntb_r1*ntb_r)
    endif 
  end subroutine initialize_arrays_qr_acr_qs


  subroutine initialize_arrays_qr_acr_qg(table_size)
    !! initialize data arrays for rain-graupel collection

    integer, intent(out), optional :: table_size

     ! Rain-graupel
    if (.not. allocated(tcg_racg)) allocate(tcg_racg(ntb_g1,ntb_g,nrhg,ntb_r1,ntb_r), source=0._table_dp)
    if (.not. allocated(tmr_racg)) allocate(tmr_racg(ntb_g1,ntb_g,nrhg,ntb_r1,ntb_r), source=0._table_dp)
    if (.not. allocated(tcr_gacr)) allocate(tcr_gacr(ntb_g1,ntb_g,nrhg,ntb_r1,ntb_r), source=0._table_dp)
    if (.not. allocated(tnr_racg)) allocate(tnr_racg(ntb_g1,ntb_g,nrhg,ntb_r1,ntb_r), source=0._table_dp)
    if (.not. allocated(tnr_gacr)) allocate(tnr_gacr(ntb_g1,ntb_g,nrhg,ntb_r1,ntb_r), source=0._table_dp)

    ! Table size is precision * entries * dimensions 
    if (present(table_size)) then
      table_size = table_dp * 5 * (ntb_g1*ntb_g*nrhg*ntb_r1*ntb_r)
    endif 
  end subroutine initialize_arrays_qr_acr_qg


  subroutine initialize_arrays_ccn(table_size)
    !! initialize data arrays for ccn lookup table

    integer, intent(out), optional :: table_size

    if (.not. allocated(tnccn_act)) &
      allocate(tnccn_act(ntb_arc,ntb_arw,ntb_art,ntb_arr,ntb_ark), source=0._table_sp)

    ! Table size is precision * entries * dimensions
    if (present(table_size)) then
      table_size = table_sp * 1 * (ntb_arc*ntb_arw*ntb_art*ntb_arr*ntb_ark) + (table_sp + table_sp) * 1
    endif
  end subroutine initialize_arrays_ccn


  subroutine initialize_arrays_drop_evap()
    !! initialize data arrays for drop evaporation data

      if (.not. allocated(tpc_wev)) allocate(tpc_wev(nbc,ntb_c,nbc), source=0._dp)
      if (.not. allocated(tnc_wev)) allocate(tnc_wev(nbc,ntb_c,nbc), source=0._dp)
  end subroutine initialize_arrays_drop_evap

  
  subroutine initialize_array_efsw()
    !! initializes the collision efficiency data array for snow collecting cloud water

    if (.not. allocated(t_efsw)) allocate(t_efsw(nbs,nbc), source=0._dp)
  end subroutine initialize_array_efsw


  subroutine initialize_array_efrw()
    !! initializes the collision efficiency data array for rain collecting cloud water

    if (.not. allocated(t_efrw)) allocate(t_efrw(nbr,nbc), source=0._dp)
  end subroutine initialize_array_efrw


  subroutine initialize_arrays_qi_aut_qs()
    !! initializes data arrays for cloud ice to snow conversion and growth

    if (.not. allocated(tps_iaus)) allocate(tps_iaus(ntb_i,ntb_i1), source=0._dp)
    if (.not. allocated(tni_iaus)) allocate(tni_iaus(ntb_i,ntb_i1), source=0._dp)
    if (.not. allocated(tpi_ide)) allocate(tpi_ide(ntb_i,ntb_i1), source=0._dp)
  end subroutine initialize_arrays_qi_aut_qs

end module module_mp_tempo_params

