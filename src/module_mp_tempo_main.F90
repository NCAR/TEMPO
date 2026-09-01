! ------------------------------------------------------------------
! TEMPO_STRIDE1 compile-time specialization (per-column CPU fast build)
! When -DTEMPO_STRIDE1 is set, the horizontal block iterators its/ite/jts/jte
! (i.e. bi/bj) are hardcoded to a compile-time extent of 1 so the compiler can
! drop the horizontal dimension, remove 3D address arithmetic and specialize the
! vertical (k) loops. The driver forces block_stride=1 in this build, so every
! invocation runs the per-column form regardless of the requested -s stride.
! Without the macro the tokens fall back to the original runtime bounds, so the
! default build is byte-for-byte the original hybrid code.
! ------------------------------------------------------------------
#ifdef TEMPO_STRIDE1
#define TEMPO_ITS 1
#define TEMPO_ITE 1
#define TEMPO_JTS 1
#define TEMPO_JTE 1
#else
#define TEMPO_ITS its
#define TEMPO_ITE ite
#define TEMPO_JTS jts
#define TEMPO_JTE jte
#endif
module module_mp_tempo_main
  !! main tempo microphysics code
  use module_mp_tempo_cfgs, only : ty_tempo_cfgs
  use module_mp_tempo_params, only : wp, sp, dp, &
    min_qv, roverrv, rdry, r1, r2, nt_c_max, nt_c_min, nu_c_scale, nt_c_l, t0, nrhg, rho_g, &
    meters3_to_liters, eps, aero_max, nwfa_default, nifa_default, &
    mu_r, obmr, d0r, d0r_max, org3, fv_r, av_r, low_limit_mass_for_precip, &
    mu_g, ogg3, obmg, idx_bg1, gonv_max, gonv_min, oge1, d0g, &
    av_g_old, bv_g_old, a_coeff, b_coeff, &
    r_uni, ar_volume, &
    ta_na, ntb_arc, ta_ww, ntb_arw, ta_tk, ntb_art, tnccn_act, &
    ntb_s, nis2, ntb_t, &
    nir2, nir3, ntb_r, ntb_r1, org2, org1, bm_r, am_r, crg, cre, &
    nig2, ntb_g, ntb_g1, ogg2, ogg1, bm_g, cgg, am_g, cge, cce, ccg, d0c, ocg1, ocg2, nig3, &
    nbc, ntb_c, r_c, nic2, t_nc, nic1, &
    ntb_i, ntb_i1, nii2, nii3, &
    rho_not0, ntb_in, nt_in, niin2, &
    lsub, orv, pi
  use module_mp_tempo_utils, only : get_nuc, get_constant_cloud_number, snow_moments, calc_rslf, calc_rsif
  use module_mp_nvtx, only : nvtx_range_push, nvtx_range_pop
  use module_mp_tempo_diags, only : reflectivity_10cm, effective_radius, max_hail_diam, &
    freezing_rain
  use module_mp_tempo_aerosols, only : init_ice_friendly_aerosols, init_water_friendly_aerosols, &
    aerosol_collection_efficiency
  use module_mp_tempo_ml, only : tempo_ml_predict_cloud_number
  implicit none
  private

  public :: tempo_main, ty_tempo_main_diags, tempo_tend_init, tempo_tend_finalize

#ifdef FV3
  public :: cloud_check_and_update, ice_check_and_update, snow_check_and_update
#endif
  
#ifdef unit_testing
  public :: get_cloud_table_index, get_snow_table_index, &
    get_temperature_table_index, get_rain_table_index, &
    get_graupel_table_index, get_ice_table_index
#endif
 
  type :: ty_tempo_main_diags
    real(wp), dimension(:,:), allocatable :: rain_precip
    real(wp), dimension(:,:), allocatable :: cloud_precip
    real(wp), dimension(:,:), allocatable :: ice_liquid_equiv_precip
    real(wp), dimension(:,:), allocatable :: snow_liquid_equiv_precip
    real(wp), dimension(:,:), allocatable :: graupel_liquid_equiv_precip
    real(wp), dimension(:,:), allocatable :: frozen_fraction
    real(wp), dimension(:,:), allocatable :: frz_rain_precip
    real(wp), dimension(:,:,:), allocatable :: rain_med_vol_diam
    real(wp), dimension(:,:,:), allocatable :: graupel_med_vol_diam
    real(wp), dimension(:,:,:), allocatable :: refl10cm
    real(wp), dimension(:,:,:), allocatable :: re_cloud
    real(wp), dimension(:,:,:), allocatable :: re_ice
    real(wp), dimension(:,:,:), allocatable :: re_snow
    real(wp), dimension(:,:,:), allocatable :: max_hail_diameter
    real(wp), dimension(:,:,:), allocatable :: cloud_number_mixing_ratio
  end type

  type :: ty_tend
    real(dp), allocatable :: &
      prr_wau(:,:,:), pnr_wau(:,:,:), pnc_wau(:,:,:), prr_rcw(:,:,:), pnc_rcw(:,:,:), pnr_rcr(:,:,:), &
      prs_scw(:,:,:), pnc_scw(:,:,:), png_scw(:,:,:), pbg_scw(:,:,:), prg_gcw(:,:,:), pnc_gcw(:,:,:), pbg_gcw(:,:,:), &
      pri_ihm(:,:,:), pni_ihm(:,:,:), prs_ihm(:,:,:), prg_ihm(:,:,:), prg_scw(:,:,:), &
      prr_rcs(:,:,:), pnr_rcs(:,:,:), prg_rcs(:,:,:), png_rcs(:,:,:), prs_rcs(:,:,:), pbg_rcs(:,:,:), &
      prr_rcg(:,:,:), pnr_rcg(:,:,:), prg_rcg(:,:,:), png_rcg(:,:,:), pbg_rcg(:,:,:), &
      pri_inu(:,:,:), pni_inu(:,:,:), pri_iha(:,:,:), pni_iha(:,:,:), &
      pri_wfz(:,:,:), pni_wfz(:,:,:), &
      prg_rfz(:,:,:), png_rfz(:,:,:), pnr_rfz(:,:,:), pri_rfz(:,:,:), pni_rfz(:,:,:), pbg_rfz(:,:,:), &
      prs_sde(:,:,:), pri_ide(:,:,:), pni_ide(:,:,:), prs_ide(:,:,:), prg_gde(:,:,:), png_gde(:,:,:), &
      pni_iau(:,:,:), prs_iau(:,:,:), &
      prr_sml(:,:,:), prr_gml(:,:,:), pbg_sml(:,:,:), pbg_gml(:,:,:), pnr_sml(:,:,:), pnr_gml(:,:,:), &
      prr_rci(:,:,:), pnr_rci(:,:,:), pri_rci(:,:,:), pni_rci(:,:,:), prg_rci(:,:,:), png_rci(:,:,:), pbg_rci(:,:,:), &
      pni_sci(:,:,:), prs_sci(:,:,:), &
      prw_vcd(:,:,:), pnc_wcd(:,:,:), prv_rev(:,:,:), pnr_rev(:,:,:), &
      pna_rca(:,:,:), pna_sca(:,:,:), pna_gca(:,:,:), pnd_rcd(:,:,:), pnd_scd(:,:,:), pnd_gcd(:,:,:)
  end type


  interface get_cloud_number
    module procedure tempo_ml_predict_cloud_number
    module procedure get_constant_cloud_number
  end interface

  !! Module-level persistent tendency workspace. Allocated once via tempo_tend_init
  !! (called from the host driver before the timestep loop) and freed via tempo_tend_finalize.
  !! tempo_main associates its local `tend` to this module variable.
  type(ty_tend), save :: tend_tile_module
  logical, save :: tend_tile_module_allocated = .false.

  contains

  subroutine ty_tend_allocate_tile(tend, kts, kte, its, ite, jts, jte)
    type(ty_tend), intent(inout) :: tend
    integer, intent(in) :: kts, kte, its, ite, jts, jte
    allocate(tend%prr_wau(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), source=0._dp)
    allocate(tend%pnr_wau(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), source=0._dp)
    allocate(tend%pnc_wau(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), source=0._dp)
    allocate(tend%prr_rcw(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), source=0._dp)
    allocate(tend%pnc_rcw(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), source=0._dp)
    allocate(tend%pnr_rcr(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), source=0._dp)
    allocate(tend%prs_scw(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), source=0._dp)
    allocate(tend%pnc_scw(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), source=0._dp)
    allocate(tend%png_scw(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), source=0._dp)
    allocate(tend%pbg_scw(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), source=0._dp)
    allocate(tend%prg_gcw(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), source=0._dp)
    allocate(tend%pnc_gcw(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), source=0._dp)
    allocate(tend%pbg_gcw(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), source=0._dp)
    allocate(tend%pri_ihm(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), source=0._dp)
    allocate(tend%pni_ihm(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), source=0._dp)
    allocate(tend%prs_ihm(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), source=0._dp)
    allocate(tend%prg_ihm(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), source=0._dp)
    allocate(tend%prg_scw(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), source=0._dp)
    allocate(tend%prr_rcs(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), source=0._dp)
    allocate(tend%pnr_rcs(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), source=0._dp)
    allocate(tend%prg_rcs(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), source=0._dp)
    allocate(tend%png_rcs(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), source=0._dp)
    allocate(tend%prs_rcs(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), source=0._dp)
    allocate(tend%pbg_rcs(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), source=0._dp)
    allocate(tend%prr_rcg(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), source=0._dp)
    allocate(tend%pnr_rcg(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), source=0._dp)
    allocate(tend%prg_rcg(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), source=0._dp)
    allocate(tend%png_rcg(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), source=0._dp)
    allocate(tend%pbg_rcg(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), source=0._dp)
    allocate(tend%pri_inu(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), source=0._dp)
    allocate(tend%pni_inu(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), source=0._dp)
    allocate(tend%pri_iha(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), source=0._dp)
    allocate(tend%pni_iha(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), source=0._dp)
    allocate(tend%pri_wfz(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), source=0._dp)
    allocate(tend%pni_wfz(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), source=0._dp)
    allocate(tend%prg_rfz(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), source=0._dp)
    allocate(tend%png_rfz(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), source=0._dp)
    allocate(tend%pnr_rfz(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), source=0._dp)
    allocate(tend%pri_rfz(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), source=0._dp)
    allocate(tend%pni_rfz(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), source=0._dp)
    allocate(tend%pbg_rfz(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), source=0._dp)
    allocate(tend%prs_sde(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), source=0._dp)
    allocate(tend%pri_ide(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), source=0._dp)
    allocate(tend%pni_ide(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), source=0._dp)
    allocate(tend%prs_ide(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), source=0._dp)
    allocate(tend%prg_gde(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), source=0._dp)
    allocate(tend%png_gde(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), source=0._dp)
    allocate(tend%pni_iau(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), source=0._dp)
    allocate(tend%prs_iau(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), source=0._dp)
    allocate(tend%prr_sml(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), source=0._dp)
    allocate(tend%prr_gml(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), source=0._dp)
    allocate(tend%pbg_sml(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), source=0._dp)
    allocate(tend%pbg_gml(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), source=0._dp)
    allocate(tend%pnr_sml(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), source=0._dp)
    allocate(tend%pnr_gml(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), source=0._dp)
    allocate(tend%prr_rci(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), source=0._dp)
    allocate(tend%pnr_rci(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), source=0._dp)
    allocate(tend%pri_rci(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), source=0._dp)
    allocate(tend%pni_rci(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), source=0._dp)
    allocate(tend%prg_rci(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), source=0._dp)
    allocate(tend%png_rci(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), source=0._dp)
    allocate(tend%pbg_rci(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), source=0._dp)
    allocate(tend%pni_sci(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), source=0._dp)
    allocate(tend%prs_sci(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), source=0._dp)
    allocate(tend%prw_vcd(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), source=0._dp)
    allocate(tend%pnc_wcd(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), source=0._dp)
    allocate(tend%prv_rev(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), source=0._dp)
    allocate(tend%pnr_rev(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), source=0._dp)
    allocate(tend%pna_rca(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), source=0._dp)
    allocate(tend%pna_sca(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), source=0._dp)
    allocate(tend%pna_gca(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), source=0._dp)
    allocate(tend%pnd_rcd(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), source=0._dp)
    allocate(tend%pnd_scd(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), source=0._dp)
    allocate(tend%pnd_gcd(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), source=0._dp)
  end subroutine ty_tend_allocate_tile

  subroutine ty_tend_deallocate_tile(tend)
    type(ty_tend), intent(inout) :: tend
    if (allocated(tend%prr_wau)) deallocate(tend%prr_wau)
    if (allocated(tend%pnr_wau)) deallocate(tend%pnr_wau)
    if (allocated(tend%pnc_wau)) deallocate(tend%pnc_wau)
    if (allocated(tend%prr_rcw)) deallocate(tend%prr_rcw)
    if (allocated(tend%pnc_rcw)) deallocate(tend%pnc_rcw)
    if (allocated(tend%pnr_rcr)) deallocate(tend%pnr_rcr)
    if (allocated(tend%prs_scw)) deallocate(tend%prs_scw)
    if (allocated(tend%pnc_scw)) deallocate(tend%pnc_scw)
    if (allocated(tend%png_scw)) deallocate(tend%png_scw)
    if (allocated(tend%pbg_scw)) deallocate(tend%pbg_scw)
    if (allocated(tend%prg_gcw)) deallocate(tend%prg_gcw)
    if (allocated(tend%pnc_gcw)) deallocate(tend%pnc_gcw)
    if (allocated(tend%pbg_gcw)) deallocate(tend%pbg_gcw)
    if (allocated(tend%pri_ihm)) deallocate(tend%pri_ihm)
    if (allocated(tend%pni_ihm)) deallocate(tend%pni_ihm)
    if (allocated(tend%prs_ihm)) deallocate(tend%prs_ihm)
    if (allocated(tend%prg_ihm)) deallocate(tend%prg_ihm)
    if (allocated(tend%prg_scw)) deallocate(tend%prg_scw)
    if (allocated(tend%prr_rcs)) deallocate(tend%prr_rcs)
    if (allocated(tend%pnr_rcs)) deallocate(tend%pnr_rcs)
    if (allocated(tend%prg_rcs)) deallocate(tend%prg_rcs)
    if (allocated(tend%png_rcs)) deallocate(tend%png_rcs)
    if (allocated(tend%prs_rcs)) deallocate(tend%prs_rcs)
    if (allocated(tend%pbg_rcs)) deallocate(tend%pbg_rcs)
    if (allocated(tend%prr_rcg)) deallocate(tend%prr_rcg)
    if (allocated(tend%pnr_rcg)) deallocate(tend%pnr_rcg)
    if (allocated(tend%prg_rcg)) deallocate(tend%prg_rcg)
    if (allocated(tend%png_rcg)) deallocate(tend%png_rcg)
    if (allocated(tend%pbg_rcg)) deallocate(tend%pbg_rcg)
    if (allocated(tend%pri_inu)) deallocate(tend%pri_inu)
    if (allocated(tend%pni_inu)) deallocate(tend%pni_inu)
    if (allocated(tend%pri_iha)) deallocate(tend%pri_iha)
    if (allocated(tend%pni_iha)) deallocate(tend%pni_iha)
    if (allocated(tend%pri_wfz)) deallocate(tend%pri_wfz)
    if (allocated(tend%pni_wfz)) deallocate(tend%pni_wfz)
    if (allocated(tend%prg_rfz)) deallocate(tend%prg_rfz)
    if (allocated(tend%png_rfz)) deallocate(tend%png_rfz)
    if (allocated(tend%pnr_rfz)) deallocate(tend%pnr_rfz)
    if (allocated(tend%pri_rfz)) deallocate(tend%pri_rfz)
    if (allocated(tend%pni_rfz)) deallocate(tend%pni_rfz)
    if (allocated(tend%pbg_rfz)) deallocate(tend%pbg_rfz)
    if (allocated(tend%prs_sde)) deallocate(tend%prs_sde)
    if (allocated(tend%pri_ide)) deallocate(tend%pri_ide)
    if (allocated(tend%pni_ide)) deallocate(tend%pni_ide)
    if (allocated(tend%prs_ide)) deallocate(tend%prs_ide)
    if (allocated(tend%prg_gde)) deallocate(tend%prg_gde)
    if (allocated(tend%png_gde)) deallocate(tend%png_gde)
    if (allocated(tend%pni_iau)) deallocate(tend%pni_iau)
    if (allocated(tend%prs_iau)) deallocate(tend%prs_iau)
    if (allocated(tend%prr_sml)) deallocate(tend%prr_sml)
    if (allocated(tend%prr_gml)) deallocate(tend%prr_gml)
    if (allocated(tend%pbg_sml)) deallocate(tend%pbg_sml)
    if (allocated(tend%pbg_gml)) deallocate(tend%pbg_gml)
    if (allocated(tend%pnr_sml)) deallocate(tend%pnr_sml)
    if (allocated(tend%pnr_gml)) deallocate(tend%pnr_gml)
    if (allocated(tend%prr_rci)) deallocate(tend%prr_rci)
    if (allocated(tend%pnr_rci)) deallocate(tend%pnr_rci)
    if (allocated(tend%pri_rci)) deallocate(tend%pri_rci)
    if (allocated(tend%pni_rci)) deallocate(tend%pni_rci)
    if (allocated(tend%prg_rci)) deallocate(tend%prg_rci)
    if (allocated(tend%png_rci)) deallocate(tend%png_rci)
    if (allocated(tend%pbg_rci)) deallocate(tend%pbg_rci)
    if (allocated(tend%pni_sci)) deallocate(tend%pni_sci)
    if (allocated(tend%prs_sci)) deallocate(tend%prs_sci)
    if (allocated(tend%prw_vcd)) deallocate(tend%prw_vcd)
    if (allocated(tend%pnc_wcd)) deallocate(tend%pnc_wcd)
    if (allocated(tend%prv_rev)) deallocate(tend%prv_rev)
    if (allocated(tend%pnr_rev)) deallocate(tend%pnr_rev)
    if (allocated(tend%pna_rca)) deallocate(tend%pna_rca)
    if (allocated(tend%pna_sca)) deallocate(tend%pna_sca)
    if (allocated(tend%pna_gca)) deallocate(tend%pna_gca)
    if (allocated(tend%pnd_rcd)) deallocate(tend%pnd_rcd)
    if (allocated(tend%pnd_scd)) deallocate(tend%pnd_scd)
    if (allocated(tend%pnd_gcd)) deallocate(tend%pnd_gcd)
  end subroutine ty_tend_deallocate_tile

  subroutine tempo_tend_init(kts, kte, its, ite, jts, jte)
    !! One-shot allocation of the module-level tendency workspace and matching device-side create.
    !! Call once from the host driver before the timestep loop that invokes tempo_run/tempo_main.
    !! Subsequent calls are no-ops if already allocated.
    integer, intent(in) :: kts, kte, its, ite, jts, jte
    if (tend_tile_module_allocated) return
    call ty_tend_allocate_tile(tend_tile_module, kts, kte, its, ite, jts, jte)
    !$acc enter data copyin(tend_tile_module)
    !$acc enter data create(tend_tile_module%prr_wau, tend_tile_module%pnr_wau, tend_tile_module%pnc_wau, &
    !$acc                   tend_tile_module%prr_rcw, tend_tile_module%pnc_rcw, tend_tile_module%pnr_rcr, &
    !$acc                   tend_tile_module%prs_scw, tend_tile_module%pnc_scw, tend_tile_module%png_scw, &
    !$acc                   tend_tile_module%pbg_scw, tend_tile_module%prg_gcw, tend_tile_module%pnc_gcw, &
    !$acc                   tend_tile_module%pbg_gcw, tend_tile_module%pri_ihm, tend_tile_module%pni_ihm, &
    !$acc                   tend_tile_module%prs_ihm, tend_tile_module%prg_ihm, tend_tile_module%prg_scw, &
    !$acc                   tend_tile_module%prr_rcs, tend_tile_module%pnr_rcs, tend_tile_module%prg_rcs, &
    !$acc                   tend_tile_module%png_rcs, tend_tile_module%prs_rcs, tend_tile_module%pbg_rcs, &
    !$acc                   tend_tile_module%prr_rcg, tend_tile_module%pnr_rcg, tend_tile_module%prg_rcg, &
    !$acc                   tend_tile_module%png_rcg, tend_tile_module%pbg_rcg, &
    !$acc                   tend_tile_module%pri_inu, tend_tile_module%pni_inu, tend_tile_module%pri_iha, &
    !$acc                   tend_tile_module%pni_iha, tend_tile_module%pri_wfz, tend_tile_module%pni_wfz, &
    !$acc                   tend_tile_module%prg_rfz, tend_tile_module%png_rfz, tend_tile_module%pnr_rfz, &
    !$acc                   tend_tile_module%pri_rfz, tend_tile_module%pni_rfz, tend_tile_module%pbg_rfz, &
    !$acc                   tend_tile_module%prs_sde, tend_tile_module%pri_ide, tend_tile_module%pni_ide, &
    !$acc                   tend_tile_module%prs_ide, tend_tile_module%prg_gde, tend_tile_module%png_gde, &
    !$acc                   tend_tile_module%pni_iau, tend_tile_module%prs_iau, &
    !$acc                   tend_tile_module%prr_sml, tend_tile_module%prr_gml, tend_tile_module%pbg_sml, &
    !$acc                   tend_tile_module%pbg_gml, tend_tile_module%pnr_sml, tend_tile_module%pnr_gml, &
    !$acc                   tend_tile_module%prr_rci, tend_tile_module%pnr_rci, tend_tile_module%pri_rci, &
    !$acc                   tend_tile_module%pni_rci, tend_tile_module%prg_rci, tend_tile_module%png_rci, &
    !$acc                   tend_tile_module%pbg_rci, tend_tile_module%pni_sci, tend_tile_module%prs_sci, &
    !$acc                   tend_tile_module%prw_vcd, tend_tile_module%pnc_wcd, tend_tile_module%prv_rev, &
    !$acc                   tend_tile_module%pnr_rev, tend_tile_module%pna_rca, tend_tile_module%pna_sca, &
    !$acc                   tend_tile_module%pna_gca, tend_tile_module%pnd_rcd, tend_tile_module%pnd_scd, &
    !$acc                   tend_tile_module%pnd_gcd)
    tend_tile_module_allocated = .true.
  end subroutine tempo_tend_init

  subroutine tempo_tend_finalize()
    !! Release the device images and host-side allocations of the module-level tendency workspace.
    !! Call once from the host driver after the timestep loop. Safe to call when not allocated (no-op).
    if (.not. tend_tile_module_allocated) return
    !$acc exit data delete(tend_tile_module%prr_wau, tend_tile_module%pnr_wau, tend_tile_module%pnc_wau, &
    !$acc                  tend_tile_module%prr_rcw, tend_tile_module%pnc_rcw, tend_tile_module%pnr_rcr, &
    !$acc                  tend_tile_module%prs_scw, tend_tile_module%pnc_scw, tend_tile_module%png_scw, &
    !$acc                  tend_tile_module%pbg_scw, tend_tile_module%prg_gcw, tend_tile_module%pnc_gcw, &
    !$acc                  tend_tile_module%pbg_gcw, tend_tile_module%pri_ihm, tend_tile_module%pni_ihm, &
    !$acc                  tend_tile_module%prs_ihm, tend_tile_module%prg_ihm, tend_tile_module%prg_scw, &
    !$acc                  tend_tile_module%prr_rcs, tend_tile_module%pnr_rcs, tend_tile_module%prg_rcs, &
    !$acc                  tend_tile_module%png_rcs, tend_tile_module%prs_rcs, tend_tile_module%pbg_rcs, &
    !$acc                  tend_tile_module%prr_rcg, tend_tile_module%pnr_rcg, tend_tile_module%prg_rcg, &
    !$acc                  tend_tile_module%png_rcg, tend_tile_module%pbg_rcg, &
    !$acc                  tend_tile_module%pri_inu, tend_tile_module%pni_inu, tend_tile_module%pri_iha, &
    !$acc                  tend_tile_module%pni_iha, tend_tile_module%pri_wfz, tend_tile_module%pni_wfz, &
    !$acc                  tend_tile_module%prg_rfz, tend_tile_module%png_rfz, tend_tile_module%pnr_rfz, &
    !$acc                  tend_tile_module%pri_rfz, tend_tile_module%pni_rfz, tend_tile_module%pbg_rfz, &
    !$acc                  tend_tile_module%prs_sde, tend_tile_module%pri_ide, tend_tile_module%pni_ide, &
    !$acc                  tend_tile_module%prs_ide, tend_tile_module%prg_gde, tend_tile_module%png_gde, &
    !$acc                  tend_tile_module%pni_iau, tend_tile_module%prs_iau, &
    !$acc                  tend_tile_module%prr_sml, tend_tile_module%prr_gml, tend_tile_module%pbg_sml, &
    !$acc                  tend_tile_module%pbg_gml, tend_tile_module%pnr_sml, tend_tile_module%pnr_gml, &
    !$acc                  tend_tile_module%prr_rci, tend_tile_module%pnr_rci, tend_tile_module%pri_rci, &
    !$acc                  tend_tile_module%pni_rci, tend_tile_module%prg_rci, tend_tile_module%png_rci, &
    !$acc                  tend_tile_module%pbg_rci, tend_tile_module%pni_sci, tend_tile_module%prs_sci, &
    !$acc                  tend_tile_module%prw_vcd, tend_tile_module%pnc_wcd, tend_tile_module%prv_rev, &
    !$acc                  tend_tile_module%pnr_rev, tend_tile_module%pna_rca, tend_tile_module%pna_sca, &
    !$acc                  tend_tile_module%pna_gca, tend_tile_module%pnd_rcd, tend_tile_module%pnd_scd, &
    !$acc                  tend_tile_module%pnd_gcd)
    !$acc exit data delete(tend_tile_module)
    call ty_tend_deallocate_tile(tend_tile_module)
    tend_tile_module_allocated = .false.
  end subroutine tempo_tend_finalize

  subroutine ty_tend_zero_column(tend, its_in, ite_in, jts_in, jte_in)
    type(ty_tend), intent(inout) :: tend
    integer, intent(in), optional :: its_in, ite_in, jts_in, jte_in
    integer :: i, j, k, kts, kte, its, ite, jts, jte

    !! Zero tendency accumulators on device. Components must already be present on the device
    !! (tempo_main does !$acc enter data create(tend%*) immediately after ty_tend_allocate_tile).
    !! The module workspace is allocated for the full tile, but in the hybrid block loop tempo_main
    !! processes only a sub-tile (TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE) per call. Zeroing must be scoped to that sub-tile;
    !! otherwise the full plane is re-zeroed on every block call, giving O(ncols^2) work at stride=1.
    !! Callers pass the active block bounds; absent them we fall back to the full component extent.
    kts = lbound(tend%prr_wau, 1); kte = ubound(tend%prr_wau, 1)
    its = lbound(tend%prr_wau, 2); ite = ubound(tend%prr_wau, 2)
    jts = lbound(tend%prr_wau, 3); jte = ubound(tend%prr_wau, 3)
    if (present(its_in)) its = its_in
    if (present(ite_in)) ite = ite_in
    if (present(jts_in)) jts = jts_in
    if (present(jte_in)) jte = jte_in

    !$acc parallel
    !$acc loop independent collapse(3)
    do j = TEMPO_JTS, TEMPO_JTE
      do i = TEMPO_ITS, TEMPO_ITE
        do k = kts, kte
          tend%prr_wau(k,i,j) = 0._dp
          tend%pnr_wau(k,i,j) = 0._dp
          tend%pnc_wau(k,i,j) = 0._dp
          tend%prr_rcw(k,i,j) = 0._dp
          tend%pnc_rcw(k,i,j) = 0._dp
          tend%pnr_rcr(k,i,j) = 0._dp
          tend%prs_scw(k,i,j) = 0._dp
          tend%pnc_scw(k,i,j) = 0._dp
          tend%png_scw(k,i,j) = 0._dp
          tend%pbg_scw(k,i,j) = 0._dp
          tend%prg_gcw(k,i,j) = 0._dp
          tend%pnc_gcw(k,i,j) = 0._dp
          tend%pbg_gcw(k,i,j) = 0._dp
          tend%pri_ihm(k,i,j) = 0._dp
          tend%pni_ihm(k,i,j) = 0._dp
          tend%prs_ihm(k,i,j) = 0._dp
          tend%prg_ihm(k,i,j) = 0._dp
          tend%prg_scw(k,i,j) = 0._dp
          tend%prr_rcs(k,i,j) = 0._dp
          tend%pnr_rcs(k,i,j) = 0._dp
          tend%prg_rcs(k,i,j) = 0._dp
          tend%png_rcs(k,i,j) = 0._dp
          tend%prs_rcs(k,i,j) = 0._dp
          tend%pbg_rcs(k,i,j) = 0._dp
          tend%prr_rcg(k,i,j) = 0._dp
          tend%pnr_rcg(k,i,j) = 0._dp
          tend%prg_rcg(k,i,j) = 0._dp
          tend%png_rcg(k,i,j) = 0._dp
          tend%pbg_rcg(k,i,j) = 0._dp
          tend%pri_inu(k,i,j) = 0._dp
          tend%pni_inu(k,i,j) = 0._dp
          tend%pri_iha(k,i,j) = 0._dp
          tend%pni_iha(k,i,j) = 0._dp
          tend%pri_wfz(k,i,j) = 0._dp
          tend%pni_wfz(k,i,j) = 0._dp
          tend%prg_rfz(k,i,j) = 0._dp
          tend%png_rfz(k,i,j) = 0._dp
          tend%pnr_rfz(k,i,j) = 0._dp
          tend%pri_rfz(k,i,j) = 0._dp
          tend%pni_rfz(k,i,j) = 0._dp
          tend%pbg_rfz(k,i,j) = 0._dp
          tend%prs_sde(k,i,j) = 0._dp
          tend%pri_ide(k,i,j) = 0._dp
          tend%pni_ide(k,i,j) = 0._dp
          tend%prs_ide(k,i,j) = 0._dp
          tend%prg_gde(k,i,j) = 0._dp
          tend%png_gde(k,i,j) = 0._dp
          tend%pni_iau(k,i,j) = 0._dp
          tend%prs_iau(k,i,j) = 0._dp
          tend%prr_sml(k,i,j) = 0._dp
          tend%prr_gml(k,i,j) = 0._dp
          tend%pbg_sml(k,i,j) = 0._dp
          tend%pbg_gml(k,i,j) = 0._dp
          tend%pnr_sml(k,i,j) = 0._dp
          tend%pnr_gml(k,i,j) = 0._dp
          tend%prr_rci(k,i,j) = 0._dp
          tend%pnr_rci(k,i,j) = 0._dp
          tend%pri_rci(k,i,j) = 0._dp
          tend%pni_rci(k,i,j) = 0._dp
          tend%prg_rci(k,i,j) = 0._dp
          tend%png_rci(k,i,j) = 0._dp
          tend%pbg_rci(k,i,j) = 0._dp
          tend%pni_sci(k,i,j) = 0._dp
          tend%prs_sci(k,i,j) = 0._dp
          tend%prw_vcd(k,i,j) = 0._dp
          tend%pnc_wcd(k,i,j) = 0._dp
          tend%prv_rev(k,i,j) = 0._dp
          tend%pnr_rev(k,i,j) = 0._dp
          tend%pna_rca(k,i,j) = 0._dp
          tend%pna_sca(k,i,j) = 0._dp
          tend%pna_gca(k,i,j) = 0._dp
          tend%pnd_rcd(k,i,j) = 0._dp
          tend%pnd_scd(k,i,j) = 0._dp
          tend%pnd_gcd(k,i,j) = 0._dp
        enddo
      enddo
    enddo
    !$acc end parallel
  end subroutine ty_tend_zero_column

  subroutine tempo_main(tempo_cfgs, &
      qv3d, qc3d, qi3d, qr3d, qs3d, qg3d, qb3d, ni3d, nr3d, nc3d, ng3d, &
      nwfa3d, nifa3d, t3d, p3d, w3d, dz3d, &
      qcfrac3d, qifrac3d, qc_bl3d, qcfrac_bl3d, &
      thten_bl3d, qvten_bl3d, qcten_bl3d, qiten_bl3d, &
      thten_lwrad3d, thten_swrad3d, &
      kts, kte, dt, its, ite, jts, jte, tempo_main_diags)
      
    !! TEMPO microphysics over a horizontal tile (column loop inside).
    type(ty_tempo_cfgs), intent(in) :: tempo_cfgs
    type(ty_tempo_main_diags), intent(inout) :: tempo_main_diags
    integer, intent(in) :: kts, kte, its, ite, jts, jte
    real(wp), intent(in) :: dt
    real(wp), dimension(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(inout), target :: &
      qv3d, qc3d, qi3d, qr3d, qs3d, qg3d, ni3d, nr3d
    real(wp), dimension(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(in), target :: p3d, w3d, dz3d
    real(wp), dimension(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(inout), target :: t3d
    real(wp), dimension(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(inout), optional, target :: &
      qb3d, nc3d, nwfa3d, nifa3d, ng3d, &
      qcfrac3d, qifrac3d, qc_bl3d, qcfrac_bl3d, &
      thten_bl3d, qvten_bl3d, qcten_bl3d, qiten_bl3d, thten_lwrad3d, thten_swrad3d

    !! Module-level tendency workspace is allocated once via tempo_tend_init (called by the host
    !! driver before the timestep loop). Use Fortran associate so the rest of this routine continues
    !! to reference the workspace as `tend` without touching the existing kernel directives.

    real(wp), dimension(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE) :: tten, qvten, qcten, qiten, qrten, qsten, &
      qgten, qbten, niten, nrten !! tendencies (3D over horizontal tile)
    real(wp), dimension(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE) :: ncten, ngten, nwfaten, nifaten !! tendencies (3D tile)

    logical, dimension(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE) :: l_qc, l_qi, l_qr, l_qs, l_qg !! hydrometeor flags (3D tile)
    integer, dimension(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE) :: idx_bg !! graupel density index (3D tile)

    ! thermodynamic variables (3D over horizontal tile)
    real(wp), dimension(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE) :: temp, pres, qv
    real(wp), dimension(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE) :: rho, rhof, rhof2
    real(wp), dimension(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE) :: qvs, qvsi, delqvs
    real(wp), dimension(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE) :: satw, sati, ssatw, ssati
    real(wp), dimension(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE) :: diffu, visco, vsc2, tcond, lvap, ocp, lvt2
 
    real(wp), dimension(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE) :: rc, ri, rr, rs, rg, rb !! local microphysics (3D tile)
    real(wp), dimension(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE) :: ni, nr, nc, ng, nwfa, nifa !! local microphysics (3D tile)

    real(dp), dimension(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE) :: ilamc, ilami, ilamr, ilamg !! inverse lambda (3D tile)
    real(wp), dimension(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE) :: mvd_r, mvd_c, mvd_g !! median volume diameter (3D tile)
    real(dp), dimension(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE) :: smob, smo2, smo1, smo0, smoc, smoe, smof, smog, ns, smoz !! snow moments (3D tile)
    
    real(wp), dimension(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE) :: xrx, xnx !! get_cloud_number work; j1d when no nc3d; j3g6a BL re path when nc3d present
    real(wp), allocatable :: ncsave(:,:,:) !! temporary (3D tile)

    real(wp), dimension(kts:kte+1, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE) :: vtrr, vtnr, vtrs, vtri, vtni, vtrg, vtng, vtrc, vtnc !! fallspeeds (3D tile)
    real(wp), dimension(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE) :: vtboost !! snow fallspeed boost (3D tile)
    integer :: substeps_sedi(TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), ktop_sedi(TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), n !! sedimentation (per horizontal column)
    real(wp) :: semi_sedi_factor(TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE) !! semi-lagrangian sedimentation factor (per column)

    
    ! local variables
    real(wp) :: tempc, tc0, odt
    logical :: tempo_first_main
    logical, save :: first_call_main = .true.
    integer :: i, j, k, nz, substeps_sedi_max
    logical, dimension(TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE) :: column_mp_active
    !! per-column supersaturation flag (tile-sized so split j1/j2 passes do not share one scalar)
    logical, dimension(TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE) :: column_supersaturated
    logical :: mp_active_any_k !! column_mp_active reduction over vertical (OpenACC column_j1f4)
    logical, dimension(TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE) :: qr_col_any !! column_j3a2: column has rain (OpenACC; avoids any() on device)
    logical, dimension(TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE) :: qg_col_any !! column_j3b2: column has graupel (OpenACC; avoids any() on device)
    logical, dimension(TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE) :: qs_col_any !! column_j3c: column has snow (OpenACC; avoids any() on device)
    logical, dimension(TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE) :: qi_col_any !! column_j3c2: column has ice (OpenACC; avoids any() on device)
    logical, dimension(TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE) :: qc_col_any !! column_j3c3: column has cloud (OpenACC; avoids any() on device)

    call nvtx_range_push('tempo_main')

    odt = 1._wp / dt

    nz = kte - kts + 1
    tempo_first_main = first_call_main

    if (.not. tend_tile_module_allocated) then
      !! Defensive fallback: caller forgot to call tempo_tend_init. Do a one-shot init here so the
      !! routine remains usable as a drop-in. Lifetime then matches the program (released via
      !! tempo_tend_finalize from the caller, or implicitly at program exit).
      call tempo_tend_init(kts, kte, its, ite, jts, jte)
    endif

    associate (tend => tend_tile_module)

    ! zero out all mp tendency terms for the tile (once before column loop) -- runs on device
    call nvtx_range_push('ty_tend_zero_column')
    call ty_tend_zero_column(tend, its, ite, jts, jte)
    call nvtx_range_pop()

    if (.not. present(nc3d)) then
      if (.not. allocated(ncsave)) then
        allocate(ncsave(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), source=0._wp)
        !$acc enter data create(ncsave)
      endif
    endif

    !! tempo_main_diags%* are allocated + entered onto the device by tempo_run (driver).
    !! By the time this routine runs, the device images are already present.

    !! Structured data region covering all stack-allocated tile workspace.
    !! Closes immediately before !$acc exit data delete(tend%*) at the bottom of tempo_main so the
    !! device images outlive the last kernel touching them. (Optional argument scalars
    !! qb3d/nc3d/etc. enter the device via tempo_run's outer enter-data block when present.)
    !$acc data create(tten, qvten, qcten, qiten, qrten, qsten, qgten, qbten, niten, nrten, &
    !$acc             ncten, ngten, nwfaten, nifaten, &
    !$acc             l_qc, l_qi, l_qr, l_qs, l_qg, idx_bg, &
    !$acc             temp, pres, qv, &
    !$acc             rho, rhof, rhof2, &
    !$acc             qvs, qvsi, delqvs, &
    !$acc             satw, sati, ssatw, ssati, &
    !$acc             diffu, visco, vsc2, tcond, lvap, ocp, lvt2, &
    !$acc             rc, ri, rr, rs, rg, rb, &
    !$acc             ni, nr, nc, ng, nwfa, nifa, &
    !$acc             ilamc, ilami, ilamr, ilamg, &
    !$acc             mvd_r, mvd_c, mvd_g, &
    !$acc             smob, smo2, smo1, smo0, smoc, smoe, smof, smog, ns, smoz, &
    !$acc             xrx, xnx, &
    !$acc             vtrr, vtnr, vtrs, vtri, vtni, vtrg, vtng, vtrc, vtnc, vtboost, &
    !$acc             substeps_sedi, ktop_sedi, semi_sedi_factor, &
    !$acc             column_mp_active, column_supersaturated, qr_col_any, qg_col_any, qs_col_any, qi_col_any, qc_col_any)

    !! column_j1*, j2*, j3*: multiple full-tile (j,i) passes per stage; order matches former single nest.
    !! tempo_first_main + (jts,its) preserves first_call_main semantics across passes.

    !!! Tile init (column_j1a): five OpenACC loop nests — supersat flags, zero tendencies/snow work,
    !!! zero fallspeed workspace, zero precip diagnostics, thermodynamic init from host 3d state.
    call nvtx_range_push('column_j1a_tile_init')
    !$acc parallel

    !$acc loop gang vector collapse(2)
    do j = TEMPO_JTS, TEMPO_JTE
      do i = TEMPO_ITS, TEMPO_ITE
        column_supersaturated(i,j) = .false.
      enddo
    enddo

    !$acc loop gang vector collapse(3)
    do j = TEMPO_JTS, TEMPO_JTE
      do i = TEMPO_ITS, TEMPO_ITE
        do k = 1, nz
          tten(k,i,j) = 0._wp
          qvten(k,i,j) = 0._wp
          qcten(k,i,j) = 0._wp
          qiten(k,i,j) = 0._wp
          qrten(k,i,j) = 0._wp
          qsten(k,i,j) = 0._wp
          qgten(k,i,j) = 0._wp
          ngten(k,i,j) = 0._wp
          qbten(k,i,j) = 0._wp
          niten(k,i,j) = 0._wp
          nrten(k,i,j) = 0._wp
          ncten(k,i,j) = 0._wp
          nwfaten(k,i,j) = 0._wp
          nifaten(k,i,j) = 0._wp
          smo0(k,i,j) = 0._dp
          smo1(k,i,j) = 0._dp
          smo2(k,i,j) = 0._dp
          smob(k,i,j) = 0._dp
          smoc(k,i,j) = 0._dp
          smoe(k,i,j) = 0._dp
          smof(k,i,j) = 0._dp
          smog(k,i,j) = 0._dp
          smoz(k,i,j) = 0._dp
          ns(k,i,j) = 0._dp
          vtboost(k,i,j) = 1._wp
        enddo
      enddo
    enddo

    !$acc loop gang vector collapse(3)
    do j = TEMPO_JTS, TEMPO_JTE
      do i = TEMPO_ITS, TEMPO_ITE
        do k = 1, nz+1
          vtrr(k,i,j) = 0._wp
          vtnr(k,i,j) = 0._wp
          vtrs(k,i,j) = 0._wp
          vtri(k,i,j) = 0._wp
          vtni(k,i,j) = 0._wp
          vtrc(k,i,j) = 0._wp
          vtnc(k,i,j) = 0._wp
          vtrg(k,i,j) = 0._wp
          vtng(k,i,j) = 0._wp
        enddo
      enddo
    enddo

    !$acc loop gang vector collapse(2)
    do j = TEMPO_JTS, TEMPO_JTE
      do i = TEMPO_ITS, TEMPO_ITE
        tempo_main_diags%rain_precip(i,j) = 0._wp
        tempo_main_diags%cloud_precip(i,j) = 0._wp
        tempo_main_diags%ice_liquid_equiv_precip(i,j) = 0._wp
        tempo_main_diags%snow_liquid_equiv_precip(i,j) = 0._wp
        tempo_main_diags%graupel_liquid_equiv_precip(i,j) = 0._wp
        tempo_main_diags%frz_rain_precip(i,j) = 0._wp
      enddo
    enddo

    !$acc loop gang vector collapse(3)
    do j = TEMPO_JTS, TEMPO_JTE
      do i = TEMPO_ITS, TEMPO_ITE
        do k = 1, nz
          temp(k,i,j) = t3d(k,i,j)
          qv(k,i,j) = max(min_qv, qv3d(k,i,j))
          pres(k,i,j) = p3d(k,i,j)
          rho(k,i,j) = roverrv*pres(k,i,j)/(rdry*temp(k,i,j)*(qv(k,i,j)+roverrv))
        enddo
      enddo
    enddo

    !$acc end parallel
    call nvtx_range_pop()

    call nvtx_range_push('column_j1b_nwfa')
    call init_water_friendly_aerosols(kts, kte, its, ite, jts, jte, tempo_first_main, dz3d, nwfa, nwfa3d=nwfa3d)
    call nvtx_range_pop()

    call nvtx_range_push('column_j1b_nifa')
    call init_ice_friendly_aerosols(kts, kte, its, ite, jts, jte, tempo_first_main, dz3d, nifa, nifa3d=nifa3d)
    call nvtx_range_pop()

    call nvtx_range_push('aerosol_check_and_update')
    call aerosol_check_and_update(dt, kts, kte, its, ite, jts, jte, rho, nwfa, nifa, nwfaten, nifaten, &
      nwfa3d=nwfa3d, nifa3d=nifa3d)
    call nvtx_range_pop()

    call nvtx_range_push('rain_check_and_update')
    call rain_check_and_update(dt, odt, kts, kte, its, ite, jts, jte, rho, l_qr, qr3d, nr3d, rr, nr, qrten, nrten, ilamr, mvd_r)
    call nvtx_range_pop()

    call nvtx_range_push('ice_check_and_update')
    call ice_check_and_update(dt, odt, kts, kte, its, ite, jts, jte, rho, l_qi, qi3d, ni3d, ri, ni, qiten, niten, ilami)
    call nvtx_range_pop()
    call nvtx_range_push('snow_check_and_update')
    call snow_check_and_update(dt, odt, kts, kte, its, ite, jts, jte, rho, l_qs, qs3d, rs, qsten)
    call nvtx_range_pop()

    call nvtx_range_push('snow_moments_column_j1c4')
    !$acc parallel
    !$acc loop gang vector collapse(3) private(tc0)
    column_j1c4: do j = TEMPO_JTS, TEMPO_JTE
      column_i1c4: do i = TEMPO_ITS, TEMPO_ITE
        ! snow moments
        do k = 1, nz
          if (l_qs(k,i,j)) then
            tc0 = min(-0.1, temp(k,i,j)-t0)
            call snow_moments(rs=rs(k,i,j), tc=tc0, &
              smob=smob(k,i,j), smoc=smoc(k,i,j), ns=ns(k,i,j), &
              smo0=smo0(k,i,j), smo1=smo1(k,i,j), smo2=smo2(k,i,j), &
              smoe=smoe(k,i,j), smof=smof(k,i,j), smog=smog(k,i,j))
          endif 
        enddo
      enddo column_i1c4
    enddo column_j1c4
    !$acc end parallel
    call nvtx_range_pop()

    call nvtx_range_push('get_cloud_number_column_j1d')
    !$acc parallel
    !$acc loop gang collapse(2)
    column_j1d: do j = TEMPO_JTS, TEMPO_JTE
      column_i1d: do i = TEMPO_ITS, TEMPO_ITE
        ! set one-moment cloud number concentration
        if (.not. present(nc3d)) then
          if (tempo_cfgs%ml_for_nc_flag) then
            xrx(:,i,j) = qc3d(:,i,j)
            where(xrx(:,i,j) <= 1.e-12_wp) xrx(:,i,j) = 0._wp
            ! ml prediction
            call get_cloud_number(nz, xrx(:,i,j), qr3d(:,i,j), qi3d(:,i,j), qs3d(:,i,j), pres(:,i,j), temp(:,i,j), w3d(:,i,j), xnx(:,i,j))
            nc(:,i,j) = xnx(:,i,j) * rho(:,i,j)
          else
            ! single modment constant value
            call get_cloud_number(nz, nc=nc(:,i,j))
          endif
          ncsave(:,i,j) = nc(:,i,j)
        endif
      enddo column_i1d
    enddo column_j1d
    !$acc end parallel
    call nvtx_range_pop()

    if (present(nc3d)) then
      call nvtx_range_push('cloud_check_and_update')
      call cloud_check_and_update(dt, odt, kts, kte, its, ite, jts, jte, rho, l_qc, qc3d, rc, nc, qcten, ncten, ilamc, mvd_c, &
        nc3d=nc3d, ncsave=ncsave)
      call nvtx_range_pop()
    else
      call nvtx_range_push('cloud_check_and_update')
      call cloud_check_and_update(dt, odt, kts, kte, its, ite, jts, jte, rho, l_qc, qc3d, rc, nc, qcten, ncten, ilamc, mvd_c, &
        ncsave=ncsave)
      call nvtx_range_pop()
    endif

    call nvtx_range_push('graupel_init_column_j1e')
    !$acc parallel
    !$acc loop gang collapse(2)
    column_j1e: do j = TEMPO_JTS, TEMPO_JTE
      column_i1e: do i = TEMPO_ITS, TEMPO_ITE
        ! init ng and qb
        if (tempo_first_main .and. j == jts .and. i == its) then
          if (present(ng3d) .and. present(qb3d)) then
            if (sum(qg3d(:,i,j)) > r1 .and. sum(ng3d(:,i,j)) < eps .and. sum(qb3d(:,i,j)) < eps) then
              call graupel_init(kts, kte, rho(kts:kte,i,j), qg3d(kts:kte,i,j), ng3d(kts:kte,i,j), qb3d(kts:kte,i,j))
            endif 
          endif
        endif 
      enddo column_i1e
    enddo column_j1e
    !$acc end parallel
    call nvtx_range_pop()

    if (present(ng3d) .and. present(qb3d)) then
      call nvtx_range_push('graupel_check_and_update')
      call graupel_check_and_update(dt, odt, kts, kte, its, ite, jts, jte, rho, l_qg, qg3d, rg, ng, rb, idx_bg, qgten, ngten, qbten, &
        ilamg, mvd_g, ng3d=ng3d, qb3d=qb3d)
      call nvtx_range_pop()
    else
      call nvtx_range_push('graupel_check_and_update')
      call graupel_check_and_update(dt, odt, kts, kte, its, ite, jts, jte, rho, l_qg, qg3d, rg, ng, rb, idx_bg, qgten, ngten, qbten, &
        ilamg, mvd_g)
      call nvtx_range_pop()
    endif

    !! column_j1f: re-zero tendencies after initial checks (pure k-inner work)
    !$acc parallel
    !$acc loop gang vector collapse(3)
    do j = TEMPO_JTS, TEMPO_JTE
      do i = TEMPO_ITS, TEMPO_ITE
        do k = 1, nz
          qcten(k,i,j) = 0._wp
          qiten(k,i,j) = 0._wp
          qrten(k,i,j) = 0._wp
          qsten(k,i,j) = 0._wp
          qgten(k,i,j) = 0._wp
          ngten(k,i,j) = 0._wp
          qbten(k,i,j) = 0._wp
          niten(k,i,j) = 0._wp
          nrten(k,i,j) = 0._wp
          ncten(k,i,j) = 0._wp
        enddo
      enddo
    enddo
    !$acc end parallel

    call nvtx_range_push('thermo_vars')
    call thermo_vars(kts, kte, its, ite, jts, jte, qv, temp, pres, rho, rhof, rhof2, qvs, delqvs, qvsi, &
      satw, sati, ssatw, ssati, diffu, visco, vsc2, ocp, lvap, tcond, lvt2, column_supersaturated)
    call nvtx_range_pop()

    column_j1f3: do j = TEMPO_JTS, TEMPO_JTE
      column_i1f3: do i = TEMPO_ITS, TEMPO_ITE
        if (tempo_first_main .and. j == jts .and. i == its) first_call_main = .false.
      enddo column_i1f3
    enddo column_j1f3

    !! column_j1f4: hydrometeor mask OR supersaturation (same as any(l_q*) per category over k)
    !$acc parallel
    !$acc loop gang vector collapse(2) private(mp_active_any_k)
    do j = TEMPO_JTS, TEMPO_JTE
      do i = TEMPO_ITS, TEMPO_ITE
        mp_active_any_k = .false.
        !$acc loop seq
        do k = kts, kte
          mp_active_any_k = mp_active_any_k .or. l_qc(k,i,j) .or. l_qr(k,i,j) .or. l_qi(k,i,j) &
            .or. l_qs(k,i,j) .or. l_qg(k,i,j)
        enddo
        column_mp_active(i,j) = mp_active_any_k .or. column_supersaturated(i,j)
      enddo
    enddo
    !$acc end parallel

    !! column_j2*: full-tile passes over active columns only (same order as former single nest).

    if (.not. tempo_cfgs%turn_off_micro_flag) then
      call nvtx_range_push('warm_rain')
      call warm_rain(kts, kte, its, ite, jts, jte, rhof, l_qc, rc, nc, ilamc, mvd_c, l_qr, rr, nr, mvd_r, tend, odt, &
        column_mp_active=column_mp_active)
      call nvtx_range_pop()
      call nvtx_range_push('rain_snow_rain_graupel')
      call rain_snow_rain_graupel(kts, kte, its, ite, jts, jte, temp, l_qr, rr, nr, ilamr, l_qs, rs, l_qg, rg, ng, ilamg, idx_bg, &
        tend, odt, column_mp_active=column_mp_active)
      call nvtx_range_pop()
      call nvtx_range_push('ice_nucleation')
      call ice_nucleation(kts, kte, its, ite, jts, jte, temp, rho, w3d, qv, qvsi, ssati, ssatw, ni, smo0, rc, nc, rr, nr, ilamr, &
        tend, dt, odt, column_mp_active=column_mp_active, nifa=nifa, nwfa=nwfa)
      call nvtx_range_pop()
      call nvtx_range_push('ice_processes')
      call ice_processes(kts, kte, its, ite, jts, jte, rhof, rhof2, rho, w3d, temp, qv, qvsi, tcond, diffu, vsc2, ssati, l_qi, ri, ni, &
        ilami, l_qs, rs, smoe, smof, smo1, rr, nr, ilamr, mvd_r, l_qg, rg, ng, ilamg, idx_bg, tend, odt, &
        column_mp_active=column_mp_active)
      call nvtx_range_pop()
      call nvtx_range_push('riming')
      call riming(kts, kte, its, ite, jts, jte, temp, rhof, visco, l_qc, rc, nc, ilamc, mvd_c, l_qs, rs, smo0, smob, smoc, smoe, &
        vtboost, l_qg, rg, ng, ilamg, idx_bg, tend, odt, column_mp_active=column_mp_active)
      call nvtx_range_pop()
      call nvtx_range_push('melting')
      call melting(kts, kte, its, ite, jts, jte, rhof2, rho, temp, qvsi, tcond, diffu, vsc2, ssati, delqvs, l_qs, rs, smof, smo0, &
        smo1, l_qg, rg, ng, ilamg, idx_bg, tend, dt, odt, column_mp_active=column_mp_active)
      call nvtx_range_pop()
      call nvtx_range_push('aerosol_scavenging')
      call aerosol_scavenging(kts, kte, its, ite, jts, jte, temp, rho, rhof, visco, nwfa, nifa, l_qr, nr, ilamr, mvd_r, l_qs, rs, &
        smob, smoc, smoe, l_qg, rg, ng, ilamg, idx_bg, tend, odt, column_mp_active=column_mp_active)
      call nvtx_range_pop()
    endif

    call nvtx_range_push('check_over_depletion')
    call check_over_depletion(kts, kte, its, ite, jts, jte, rho, temp, qvsi, qv, l_qc, rc, l_qi, ri, l_qr, rr, l_qs, rs, l_qg, rg, &
      tend, odt, column_mp_active=column_mp_active)
    call nvtx_range_pop()
    call nvtx_range_push('sum_tendencies')
    call sum_tendencies(kts, kte, its, ite, jts, jte, rho, temp, idx_bg, lvap, ocp, tend, tten, qvten, qcten, ncten, qiten, niten, &
      qsten, qrten, nrten, qgten, ngten, qbten, column_mp_active=column_mp_active)
    call nvtx_range_pop()

    !! column_j2b2: update state after tendencies (masked columns)
    !$acc parallel
    !$acc loop gang vector collapse(3) private(tempc)
    do j = TEMPO_JTS, TEMPO_JTE
      do i = TEMPO_ITS, TEMPO_ITE
        do k = 1, nz
          if (column_mp_active(i,j)) then
            temp(k,i,j) = t3d(k,i,j) + tten(k,i,j)*dt
            tempc = temp(k,i,j) - t0
            qv(k,i,j) = max(min_qv, qv3d(k,i,j) + qvten(k,i,j)*dt)
            rho(k,i,j) = roverrv*pres(k,i,j)/(rdry*temp(k,i,j)*(qv(k,i,j)+roverrv))
            nwfaten(k,i,j) = nwfaten(k,i,j) - (tend%pna_rca(k,i,j) + tend%pna_sca(k,i,j) + tend%pna_gca(k,i,j) + &
              tend%pni_iha(k,i,j)) / rho(k,i,j)
            nifaten(k,i,j) = nifaten(k,i,j) - (tend%pnd_rcd(k,i,j) + tend%pnd_scd(k,i,j) + tend%pnd_gcd(k,i,j)) / rho(k,i,j)
          endif
        enddo
      enddo
    enddo
    !$acc end parallel

    call nvtx_range_push('aerosol_check_and_update')
    call aerosol_check_and_update(dt, kts, kte, its, ite, jts, jte, rho, nwfa, nifa, nwfaten, nifaten, &
      nwfa3d=nwfa3d, nifa3d=nifa3d, column_mp_active=column_mp_active)
    call nvtx_range_pop()

    if (present(nc3d)) then
      call nvtx_range_push('cloud_check_and_update')
      call cloud_check_and_update(dt, odt, kts, kte, its, ite, jts, jte, rho, l_qc, qc3d, rc, nc, qcten, ncten, ilamc, mvd_c, &
        nc3d=nc3d, ncsave=ncsave, column_mp_active=column_mp_active, update_qc_nc_state=.false.)
      call nvtx_range_pop()
    else
      call nvtx_range_push('cloud_check_and_update')
      call cloud_check_and_update(dt, odt, kts, kte, its, ite, jts, jte, rho, l_qc, qc3d, rc, nc, qcten, ncten, ilamc, mvd_c, &
        ncsave=ncsave, column_mp_active=column_mp_active, update_qc_nc_state=.false.)
      call nvtx_range_pop()
    endif

    call nvtx_range_push('rain_check_and_update')
    call rain_check_and_update(dt, odt, kts, kte, its, ite, jts, jte, rho, l_qr, qr3d, nr3d, rr, nr, qrten, nrten, ilamr, mvd_r, &
      column_mp_active=column_mp_active, update_qr_nr_state=.false.)
    call nvtx_range_pop()

    call nvtx_range_push('ice_check_and_update')
    call ice_check_and_update(dt, odt, kts, kte, its, ite, jts, jte, rho, l_qi, qi3d, ni3d, ri, ni, qiten, niten, ilami, &
      column_mp_active=column_mp_active, update_qi_ni_state=.false.)
    call nvtx_range_pop()
    call nvtx_range_push('snow_check_and_update')
    call snow_check_and_update(dt, odt, kts, kte, its, ite, jts, jte, rho, l_qs, qs3d, rs, qsten, &
      column_mp_active=column_mp_active, update_qs_state=.false.)
    call nvtx_range_pop()

    call nvtx_range_push('snow_moments_column_j2c5')
    !$acc parallel
    !$acc loop gang collapse(2) 
    column_j2c5: do j = TEMPO_JTS, TEMPO_JTE
      column_i2c5: do i = TEMPO_ITS, TEMPO_ITE
        if (column_mp_active(i,j)) then
          ! snow moments
          !$acc loop vector private(tc0)
          do k = 1, nz
            smo0(k,i,j) = 0._dp
            smo1(k,i,j) = 0._dp
            smo2(k,i,j) = 0._dp
            smob(k,i,j) = 0._dp
            smoc(k,i,j) = 0._dp
            smoe(k,i,j) = 0._dp
            smof(k,i,j) = 0._dp
            smog(k,i,j) = 0._dp
            ns(k,i,j) = 0._dp
            if (l_qs(k,i,j)) then
              tc0 = min(-0.1, temp(k,i,j)-t0)
              call snow_moments(rs=rs(k,i,j), tc=tc0, &
                smob=smob(k,i,j), smoc=smoc(k,i,j), &
                smo2=smo2(k,i,j))
            endif
          enddo
        endif
      enddo column_i2c5
    enddo column_j2c5
    !$acc end parallel
    call nvtx_range_pop()

    if (present(ng3d) .and. present(qb3d)) then
      call nvtx_range_push('graupel_check_and_update')
      call graupel_check_and_update(dt, odt, kts, kte, its, ite, jts, jte, rho, l_qg, qg3d, rg, ng, rb, idx_bg, qgten, ngten, qbten, &
        ilamg, mvd_g, column_mp_active=column_mp_active, update_qg_ng_qb_state=.false., ng3d=ng3d, qb3d=qb3d)
      call nvtx_range_pop()
    else
      call nvtx_range_push('graupel_check_and_update')
      call graupel_check_and_update(dt, odt, kts, kte, its, ite, jts, jte, rho, l_qg, qg3d, rg, ng, rb, idx_bg, qgten, ngten, qbten, &
        ilamg, mvd_g, column_mp_active=column_mp_active, update_qg_ng_qb_state=.false.)
      call nvtx_range_pop()
    endif

    call nvtx_range_push('thermo_vars')
    call thermo_vars(kts, kte, its, ite, jts, jte, qv, temp, pres, rho, rhof, rhof2, qvs, delqvs, qvsi, &
      satw, sati, ssatw, ssati, diffu, visco, vsc2, ocp, lvap, tcond, lvt2, column_supersaturated, &
      column_mp_active=column_mp_active)
    call nvtx_range_pop()

    ! cloud condensation (tile-wide; loops over i,j inside cloud_condensation)
    if (.not. tempo_cfgs%turn_off_micro_flag .and. tempo_cfgs%cloud_condensation_flag) then
      call nvtx_range_push('cloud_condensation')
      call cloud_condensation(kts, kte, its, ite, jts, jte, rho, temp, w3d, &
        ssatw, lvap, tcond, diffu, lvt2, nwfa, qv, qvs, l_qc, rc, nc, &
        tend, dt, odt, column_mp_active=column_mp_active)
      call nvtx_range_pop()
    endif

    !! column_j2e2: accumulate condensation tendencies (masked + flags)
    !$acc parallel
    !$acc loop gang vector collapse(3)
    do j = TEMPO_JTS, TEMPO_JTE
      do i = TEMPO_ITS, TEMPO_ITE
        do k = 1, nz
          if (column_mp_active(i,j) .and. (.not. tempo_cfgs%turn_off_micro_flag) .and. tempo_cfgs%cloud_condensation_flag) then
            qvten(k,i,j) = qvten(k,i,j) - tend%prw_vcd(k,i,j)
            qcten(k,i,j) = qcten(k,i,j) + tend%prw_vcd(k,i,j)
            ncten(k,i,j) = ncten(k,i,j) + tend%pnc_wcd(k,i,j)
            nwfaten(k,i,j) = nwfaten(k,i,j) - tend%pnc_wcd(k,i,j)
            tten(k,i,j) = tten(k,i,j) + lvap(k,i,j)*ocp(k,i,j)*tend%prw_vcd(k,i,j)
          endif
        enddo
      enddo
    enddo
    !$acc end parallel

    if (.not. tempo_cfgs%turn_off_micro_flag .and. tempo_cfgs%cloud_condensation_flag) then
      if (present(nc3d)) then
        call nvtx_range_push('cloud_check_and_update')
        call cloud_check_and_update(dt, odt, kts, kte, its, ite, jts, jte, rho, l_qc, qc3d, rc, nc, qcten, ncten, ilamc, mvd_c, &
          nc3d=nc3d, ncsave=ncsave, column_mp_active=column_mp_active, update_qc_nc_state=.false.)
        call nvtx_range_pop()
      else
        call nvtx_range_push('cloud_check_and_update')
        call cloud_check_and_update(dt, odt, kts, kte, its, ite, jts, jte, rho, l_qc, qc3d, rc, nc, qcten, ncten, ilamc, mvd_c, &
          ncsave=ncsave, column_mp_active=column_mp_active, update_qc_nc_state=.false.)
        call nvtx_range_pop()
      endif
    endif

    !! column_j2e4: apply condensation tendencies to qv, T, rho (masked + flags)
    !$acc parallel
    !$acc loop gang vector collapse(3)
    do j = TEMPO_JTS, TEMPO_JTE
      do i = TEMPO_ITS, TEMPO_ITE
        do k = 1, nz
          if (column_mp_active(i,j) .and. (.not. tempo_cfgs%turn_off_micro_flag) .and. tempo_cfgs%cloud_condensation_flag) then
            qv(k,i,j) = max(min_qv, qv3d(k,i,j) + qvten(k,i,j)*dt)
            temp(k,i,j) = t3d(k,i,j) + tten(k,i,j)*dt
            rho(k,i,j) = roverrv*pres(k,i,j)/(rdry*temp(k,i,j)*(qv(k,i,j)+roverrv))
          endif
        enddo
      enddo
    enddo
    !$acc end parallel

    if (.not. tempo_cfgs%turn_off_micro_flag .and. tempo_cfgs%cloud_condensation_flag) then
      call nvtx_range_push('thermo_vars')
      call thermo_vars(kts, kte, its, ite, jts, jte, qv, temp, pres, rho, rhof, rhof2, qvs, delqvs, qvsi, &
        satw, sati, ssatw, ssati, diffu, visco, vsc2, ocp, lvap, tcond, lvt2, column_supersaturated, &
        column_mp_active=column_mp_active)
      call nvtx_range_pop()
    endif

    if (.not. tempo_cfgs%turn_off_micro_flag) then
      call nvtx_range_push('rain_evaporation')
      call rain_evaporation(kts, kte, its, ite, jts, jte, rho, temp, ssatw, lvap, tcond, diffu, vsc2, rhof2, &
        qv, qvs, l_qr, rr, nr, ilamr, tend, odt, column_mp_active=column_mp_active)
      call nvtx_range_pop()
    endif

    !! column_j2f2: rain evaporation tendency increments (masked + turn_off_micro)
    !$acc parallel
    !$acc loop gang vector collapse(3)
    do j = TEMPO_JTS, TEMPO_JTE
      do i = TEMPO_ITS, TEMPO_ITE
        do k = 1, nz
          if (column_mp_active(i,j) .and. (.not. tempo_cfgs%turn_off_micro_flag)) then
            qrten(k,i,j) = qrten(k,i,j) - tend%prv_rev(k,i,j)
            qvten(k,i,j) = qvten(k,i,j) + tend%prv_rev(k,i,j)
            nrten(k,i,j) = nrten(k,i,j) - tend%pnr_rev(k,i,j)
            nwfaten(k,i,j) = nwfaten(k,i,j) + tend%pnr_rev(k,i,j)
            tten(k,i,j) = tten(k,i,j) - lvap(k,i,j)*ocp(k,i,j)*tend%prv_rev(k,i,j)
          endif
        enddo
      enddo
    enddo
    !$acc end parallel

    if (.not. tempo_cfgs%turn_off_micro_flag) then
      call nvtx_range_push('rain_check_and_update')
      call rain_check_and_update(dt, odt, kts, kte, its, ite, jts, jte, rho, l_qr, qr3d, nr3d, rr, nr, qrten, nrten, ilamr, mvd_r, &
        column_mp_active=column_mp_active, update_qr_nr_state=.false.)
      call nvtx_range_pop()
    endif

    !! column_j2f4: apply rain-evap tendencies to qv, T, rho (masked + turn_off_micro)
    !$acc parallel
    !$acc loop gang vector collapse(3)
    do j = TEMPO_JTS, TEMPO_JTE
      do i = TEMPO_ITS, TEMPO_ITE
        do k = 1, nz
          if (column_mp_active(i,j) .and. (.not. tempo_cfgs%turn_off_micro_flag)) then
            qv(k,i,j) = max(min_qv, qv3d(k,i,j) + qvten(k,i,j)*dt)
            temp(k,i,j) = t3d(k,i,j) + tten(k,i,j)*dt
            rho(k,i,j) = roverrv*pres(k,i,j)/(rdry*temp(k,i,j)*(qv(k,i,j)+roverrv))
          endif
        enddo
      enddo
    enddo
    !$acc end parallel

    if (.not. tempo_cfgs%turn_off_micro_flag) then
      call nvtx_range_push('thermo_vars')
      call thermo_vars(kts, kte, its, ite, jts, jte, qv, temp, pres, rho, rhof, rhof2, qvs, delqvs, qvsi, &
        satw, sati, ssatw, ssati, diffu, visco, vsc2, ocp, lvap, tcond, lvt2, column_supersaturated, &
        column_mp_active=column_mp_active)
      call nvtx_range_pop()
    endif

    !! column_j3*: full-tile passes over active columns (sedimentation, final update, diagnostics).

    !! column_j3a2_init: rain sedimentation control defaults (active columns)
    !$acc parallel
    !$acc loop gang vector collapse(2)
    do j = TEMPO_JTS, TEMPO_JTE
      do i = TEMPO_ITS, TEMPO_ITE
        if (column_mp_active(i,j)) then
          ktop_sedi(i,j) = 1
          substeps_sedi(i,j) = 1
          semi_sedi_factor(i,j) = 10._wp
        endif
      enddo
    enddo
    !$acc end parallel

    !! preliminary rain fallspeed (sets substeps_sedi, ktop_sedi) before semi-Lagrangian / Eulerian sedimentation
    call nvtx_range_push('rain_fallspeed')
    call rain_fallspeed(kts, kte, its, ite, jts, jte, dt, rhof, l_qr, rr, ilamr, dz3d, vtrr, vtnr, substeps_sedi, ktop_sedi, &
      column_mp_active, l_gr_flag=.true.)
    call nvtx_range_pop()

    !! column_j3a2_qr: per-column rain presence before sedimentation (OpenACC; avoids any() on device)
    !$acc parallel
    !$acc loop gang vector collapse(2)
    do j = TEMPO_JTS, TEMPO_JTE
      do i = TEMPO_ITS, TEMPO_ITE
        qr_col_any(i,j) = .false.
        if (column_mp_active(i,j)) then
          !$acc loop seq
          do k = kts, kte
            if (l_qr(k,i,j)) qr_col_any(i,j) = .true.
          enddo
        endif
      enddo
    enddo
    !$acc end parallel

    !! column_j3a2_substeps: semi-Lagrangian rain sedimentation substep count (active columns with rain)
    if (tempo_cfgs%semi_sedi_flag) then
      !$acc parallel
      !$acc loop gang vector collapse(2)
      do j = TEMPO_JTS, TEMPO_JTE
        do i = TEMPO_ITS, TEMPO_ITE
          if (column_mp_active(i,j) .and. qr_col_any(i,j)) then
            substeps_sedi(i,j) = max(int(substeps_sedi(i,j)/semi_sedi_factor(i,j)) + 1, 1)
          endif
        enddo
      enddo
      !$acc end parallel
    endif

    !! column_j3a2_substeps_max: tile max rain sedimentation substeps (active columns with rain)
    substeps_sedi_max = 1
    !$acc parallel
    !$acc loop gang vector collapse(2) reduction(max:substeps_sedi_max)
    do j = TEMPO_JTS, TEMPO_JTE
      do i = TEMPO_ITS, TEMPO_ITE
        if (column_mp_active(i,j) .and. qr_col_any(i,j)) then
          substeps_sedi_max = max(substeps_sedi_max, substeps_sedi(i,j))
        endif
      enddo
    enddo
    !$acc end parallel

    !! column_j3a2: rain sedimentation substeps (OpenACC: n outer; one (j,i) pass per routine)
    call nvtx_range_push('column_j3a2_rain_sedimentation')
    if (tempo_cfgs%semi_sedi_flag) then
      column_j3a2_n: do n = 1, substeps_sedi_max

        !! column_j3a2_rr: semi-Lagrangian sedimentation for rain mass
        call semilagrangian_sedimentation(kts, kte, its, ite, jts, jte, dz3d, rho, rr, qrten, vtrr, substeps_sedi, r1, dt, odt, &
          column_mp_active, qr_col_any, n, precip=tempo_main_diags%rain_precip)

        !! column_j3a2_nr: semi-Lagrangian sedimentation for rain number
        call semilagrangian_sedimentation(kts, kte, its, ite, jts, jte, dz3d, rho, nr, nrten, vtnr, substeps_sedi, r2, dt, odt, &
          column_mp_active, qr_col_any, n)

        !! column_j3a2_vt: zero rain fallspeed workspace between substeps
        !$acc parallel
        !$acc loop gang collapse(2) private(k)
        do j = TEMPO_JTS, TEMPO_JTE
          do i = TEMPO_ITS, TEMPO_ITE
            if (column_mp_active(i,j) .and. qr_col_any(i,j) .and. n <= substeps_sedi(i,j)) then
              do k = kts, kte+1
                vtrr(k,i,j) = 0._wp
                vtnr(k,i,j) = 0._wp
              enddo
            endif
          enddo
        enddo
        !$acc end parallel

        !! column_j3a2_check: rain check and update between substeps
        call rain_check_and_update(dt, odt, kts, kte, its, ite, jts, jte, rho, l_qr, qr3d, nr3d, rr, nr, qrten, nrten, ilamr, mvd_r, &
          column_mp_active=column_mp_active, update_qr_nr_state=.false., substep_mode=.true., col_any=qr_col_any, n=n, &
          substeps_sedi=substeps_sedi)

        !! column_j3a2_fallspeed: rain fallspeed for next substep
        call rain_fallspeed(kts, kte, its, ite, jts, jte, dt, rhof, l_qr, rr, ilamr, dz3d, vtrr, vtnr, substeps_sedi, ktop_sedi, &
          column_mp_active, substep_mode=.true., col_any=qr_col_any, n=n)

      enddo column_j3a2_n
    else
      do n = 1, substeps_sedi_max

        !! column_j3a2_rr: Eulerian sedimentation for rain mass
        call sedimentation(kts, kte, its, ite, jts, jte, dz3d, rho, rr, qrten, vtrr, substeps_sedi, ktop_sedi, r1, dt, &
          column_mp_active, qr_col_any, n, precip=tempo_main_diags%rain_precip)

        !! column_j3a2_nr: Eulerian sedimentation for rain number
        call sedimentation(kts, kte, its, ite, jts, jte, dz3d, rho, nr, nrten, vtnr, substeps_sedi, ktop_sedi, r2, dt, &
          column_mp_active, qr_col_any, n)

        !! column_j3a2_vt: zero rain fallspeed workspace between substeps
        !$acc parallel
        !$acc loop gang collapse(2) private(k)
        do j = TEMPO_JTS, TEMPO_JTE
          do i = TEMPO_ITS, TEMPO_ITE
            if (column_mp_active(i,j) .and. qr_col_any(i,j) .and. n <= substeps_sedi(i,j)) then
              do k = kts, kte+1
                vtrr(k,i,j) = 0._wp
                vtnr(k,i,j) = 0._wp
              enddo
            endif
          enddo
        enddo
        !$acc end parallel

        !! column_j3a2_check: rain check and update between substeps
        call rain_check_and_update(dt, odt, kts, kte, its, ite, jts, jte, rho, l_qr, qr3d, nr3d, rr, nr, qrten, nrten, ilamr, mvd_r, &
          column_mp_active=column_mp_active, update_qr_nr_state=.false., substep_mode=.true., col_any=qr_col_any, n=n, &
          substeps_sedi=substeps_sedi)

        !! column_j3a2_fallspeed: rain fallspeed for next substep
        call rain_fallspeed(kts, kte, its, ite, jts, jte, dt, rhof, l_qr, rr, ilamr, dz3d, vtrr, vtnr, substeps_sedi, ktop_sedi, &
          column_mp_active, substep_mode=.true., col_any=qr_col_any, n=n)

      enddo
    endif
    call nvtx_range_pop()

    !! column_j3b2_init: graupel sedimentation control defaults (active columns)
    !$acc parallel
    !$acc loop gang vector collapse(2)
    do j = TEMPO_JTS, TEMPO_JTE
      do i = TEMPO_ITS, TEMPO_ITE
        if (column_mp_active(i,j)) then
          ktop_sedi(i,j) = 1
          substeps_sedi(i,j) = 1
          semi_sedi_factor(i,j) = 10._wp
        endif
      enddo
    enddo
    !$acc end parallel

    !! preliminary graupel fallspeed (sets substeps_sedi, ktop_sedi) before semi-Lagrangian / Eulerian sedimentation
    if (present(qb3d)) then
      call nvtx_range_push('graupel_fallspeed')
      call graupel_fallspeed(kts, kte, its, ite, jts, jte, dt, rhof, rho, visco, l_qg, rg, rb, idx_bg, ilamg, dz3d, vtrg, vtng, &
        substeps_sedi, ktop_sedi, column_mp_active, l_gg_flag=.true., qb3d=qb3d)
      call nvtx_range_pop()
    else
      call nvtx_range_push('graupel_fallspeed')
      call graupel_fallspeed(kts, kte, its, ite, jts, jte, dt, rhof, rho, visco, l_qg, rg, rb, idx_bg, ilamg, dz3d, vtrg, vtng, &
        substeps_sedi, ktop_sedi, column_mp_active, l_gg_flag=.true.)
      call nvtx_range_pop()
    endif

    !! column_j3b2_qg: per-column graupel presence before sedimentation (OpenACC; avoids any() on device)
    !$acc parallel
    !$acc loop gang vector collapse(2)
    do j = TEMPO_JTS, TEMPO_JTE
      do i = TEMPO_ITS, TEMPO_ITE
        qg_col_any(i,j) = .false.
        if (column_mp_active(i,j)) then
          !$acc loop seq
          do k = kts, kte
            if (l_qg(k,i,j)) qg_col_any(i,j) = .true.
          enddo
        endif
      enddo
    enddo
    !$acc end parallel

    !! column_j3b2_substeps: semi-Lagrangian graupel sedimentation substep count (active columns with graupel)
    if (tempo_cfgs%semi_sedi_flag) then
      !$acc parallel
      !$acc loop gang vector collapse(2)
      do j = TEMPO_JTS, TEMPO_JTE
        do i = TEMPO_ITS, TEMPO_ITE
          if (column_mp_active(i,j) .and. qg_col_any(i,j)) then
            substeps_sedi(i,j) = max(int(substeps_sedi(i,j)/semi_sedi_factor(i,j)) + 1, 1)
          endif
        enddo
      enddo
      !$acc end parallel
    endif

    !! column_j3b2_substeps_max: tile max graupel sedimentation substeps (active columns with graupel)
    substeps_sedi_max = 1
    !$acc parallel
    !$acc loop gang vector collapse(2) reduction(max:substeps_sedi_max)
    do j = TEMPO_JTS, TEMPO_JTE
      do i = TEMPO_ITS, TEMPO_ITE
        if (column_mp_active(i,j) .and. qg_col_any(i,j)) then
          substeps_sedi_max = max(substeps_sedi_max, substeps_sedi(i,j))
        endif
      enddo
    enddo
    !$acc end parallel

    !! column_j3b2: graupel sedimentation substeps (OpenACC: n outer; one (j,i) pass per routine)
    call nvtx_range_push('column_j3b2_graupel_sedimentation')
    if (tempo_cfgs%semi_sedi_flag) then
      column_j3b2_n: do n = 1, substeps_sedi_max

        !! column_j3b2_rg: semi-Lagrangian sedimentation for graupel mass
        call semilagrangian_sedimentation(kts, kte, its, ite, jts, jte, dz3d, rho, rg, qgten, vtrg, substeps_sedi, r1, dt, odt, &
          column_mp_active, qg_col_any, n, precip=tempo_main_diags%graupel_liquid_equiv_precip)

        !! column_j3b2_ng: semi-Lagrangian sedimentation for graupel number
        call semilagrangian_sedimentation(kts, kte, its, ite, jts, jte, dz3d, rho, ng, ngten, vtng, substeps_sedi, r2, dt, odt, &
          column_mp_active, qg_col_any, n)

        !! column_j3b2_rb: semi-Lagrangian sedimentation for graupel volume
        call semilagrangian_sedimentation(kts, kte, its, ite, jts, jte, dz3d, rho, rb, qbten, vtrg, substeps_sedi, &
          meters3_to_liters*r1/rho_g(nrhg), dt, odt, column_mp_active, qg_col_any, n)

        !! column_j3b2_vt: zero graupel fallspeed workspace between substeps
        !$acc parallel
        !$acc loop gang collapse(2) private(k)
        do j = TEMPO_JTS, TEMPO_JTE
          do i = TEMPO_ITS, TEMPO_ITE
            if (column_mp_active(i,j) .and. qg_col_any(i,j) .and. n <= substeps_sedi(i,j)) then
              do k = kts, kte+1
                vtrg(k,i,j) = 0._wp
                vtng(k,i,j) = 0._wp
              enddo
            endif
          enddo
        enddo
        !$acc end parallel

        !! column_j3b2_check: graupel check and update between substeps
        if (present(ng3d) .and. present(qb3d)) then
          call graupel_check_and_update(dt, odt, kts, kte, its, ite, jts, jte, rho, l_qg, qg3d, rg, ng, rb, idx_bg, qgten, ngten, &
            qbten, ilamg, mvd_g, column_mp_active, update_qg_ng_qb_state=.false., ng3d=ng3d, qb3d=qb3d, substep_mode=.true., &
            col_any=qg_col_any, n=n, substeps_sedi=substeps_sedi)
        else
          call graupel_check_and_update(dt, odt, kts, kte, its, ite, jts, jte, rho, l_qg, qg3d, rg, ng, rb, idx_bg, qgten, ngten, &
            qbten, ilamg, mvd_g, column_mp_active, update_qg_ng_qb_state=.false., substep_mode=.true., col_any=qg_col_any, n=n, &
            substeps_sedi=substeps_sedi)
        endif

        !! column_j3b2_fallspeed: graupel fallspeed for next substep
        if (present(qb3d)) then
          call graupel_fallspeed(kts, kte, its, ite, jts, jte, dt, rhof, rho, visco, l_qg, rg, rb, idx_bg, ilamg, dz3d, &
            vtrg, vtng, substeps_sedi, ktop_sedi, column_mp_active, qb3d=qb3d, substep_mode=.true., col_any=qg_col_any, n=n)
        else
          call graupel_fallspeed(kts, kte, its, ite, jts, jte, dt, rhof, rho, visco, l_qg, rg, rb, idx_bg, ilamg, dz3d, &
            vtrg, vtng, substeps_sedi, ktop_sedi, column_mp_active, substep_mode=.true., col_any=qg_col_any, n=n)
        endif

      enddo column_j3b2_n
    else
      do n = 1, substeps_sedi_max

        !! column_j3b2_rg: Eulerian sedimentation for graupel mass
        call sedimentation(kts, kte, its, ite, jts, jte, dz3d, rho, rg, qgten, vtrg, substeps_sedi, ktop_sedi, r1, dt, &
          column_mp_active, qg_col_any, n, precip=tempo_main_diags%graupel_liquid_equiv_precip)

        !! column_j3b2_ng: Eulerian sedimentation for graupel number
        call sedimentation(kts, kte, its, ite, jts, jte, dz3d, rho, ng, ngten, vtng, substeps_sedi, ktop_sedi, r2, dt, &
          column_mp_active, qg_col_any, n)

        !! column_j3b2_rb: Eulerian sedimentation for graupel volume
        call sedimentation(kts, kte, its, ite, jts, jte, dz3d, rho, rb, qbten, vtrg, substeps_sedi, ktop_sedi, &
          meters3_to_liters*r1/rho_g(nrhg), dt, column_mp_active, qg_col_any, n)

        !! column_j3b2_vt: zero graupel fallspeed workspace between substeps
        !$acc parallel
        !$acc loop gang collapse(2) private(k)
        do j = TEMPO_JTS, TEMPO_JTE
          do i = TEMPO_ITS, TEMPO_ITE
            if (column_mp_active(i,j) .and. qg_col_any(i,j) .and. n <= substeps_sedi(i,j)) then
              do k = kts, kte+1
                vtrg(k,i,j) = 0._wp
                vtng(k,i,j) = 0._wp
              enddo
            endif
          enddo
        enddo
        !$acc end parallel

        !! column_j3b2_check: graupel check and update between substeps
        if (present(ng3d) .and. present(qb3d)) then
          call graupel_check_and_update(dt, odt, kts, kte, its, ite, jts, jte, rho, l_qg, qg3d, rg, ng, rb, idx_bg, qgten, ngten, &
            qbten, ilamg, mvd_g, column_mp_active, update_qg_ng_qb_state=.false., ng3d=ng3d, qb3d=qb3d, substep_mode=.true., &
            col_any=qg_col_any, n=n, substeps_sedi=substeps_sedi)
        else
          call graupel_check_and_update(dt, odt, kts, kte, its, ite, jts, jte, rho, l_qg, qg3d, rg, ng, rb, idx_bg, qgten, ngten, &
            qbten, ilamg, mvd_g, column_mp_active, update_qg_ng_qb_state=.false., substep_mode=.true., col_any=qg_col_any, n=n, &
            substeps_sedi=substeps_sedi)
        endif

        !! column_j3b2_fallspeed: graupel fallspeed for next substep
        if (present(qb3d)) then
          call graupel_fallspeed(kts, kte, its, ite, jts, jte, dt, rhof, rho, visco, l_qg, rg, rb, idx_bg, ilamg, dz3d, &
            vtrg, vtng, substeps_sedi, ktop_sedi, column_mp_active, qb3d=qb3d, substep_mode=.true., col_any=qg_col_any, n=n)
        else
          call graupel_fallspeed(kts, kte, its, ite, jts, jte, dt, rhof, rho, visco, l_qg, rg, rb, idx_bg, ilamg, dz3d, &
            vtrg, vtng, substeps_sedi, ktop_sedi, column_mp_active, substep_mode=.true., col_any=qg_col_any, n=n)
        endif

      enddo
    endif
    call nvtx_range_pop()

    !! column_j3c_init: snow sedimentation control defaults (active columns)
    !$acc parallel
    !$acc loop gang vector collapse(2)
    do j = TEMPO_JTS, TEMPO_JTE
      do i = TEMPO_ITS, TEMPO_ITE
        if (column_mp_active(i,j)) then
          ktop_sedi(i,j) = 1
          substeps_sedi(i,j) = 1
        endif
      enddo
    enddo
    !$acc end parallel

    !! preliminary snow fallspeed (sets substeps_sedi, ktop_sedi) before Eulerian sedimentation
    call nvtx_range_push('snow_fallspeed')
    call snow_fallspeed(kts, kte, its, ite, jts, jte, dt, rhof, l_qs, rs, tend%prr_sml, smob, smoc, rr, vtrr, dz3d, vtrs, vtboost, &
      substeps_sedi, ktop_sedi, column_mp_active, l_qs_flag=.true.)
    call nvtx_range_pop()

    !! column_j3c_qs: per-column snow presence before sedimentation (OpenACC; avoids any() on device)
    !$acc parallel
    !$acc loop gang vector collapse(2)
    do j = TEMPO_JTS, TEMPO_JTE
      do i = TEMPO_ITS, TEMPO_ITE
        qs_col_any(i,j) = .false.
        if (column_mp_active(i,j)) then
          !$acc loop seq
          do k = kts, kte
            if (l_qs(k,i,j)) qs_col_any(i,j) = .true.
          enddo
        endif
      enddo
    enddo
    !$acc end parallel

    !! column_j3c_substeps_max: tile max snow sedimentation substeps (active columns with snow)
    substeps_sedi_max = 1
    !$acc parallel
    !$acc loop gang vector collapse(2) reduction(max:substeps_sedi_max)
    do j = TEMPO_JTS, TEMPO_JTE
      do i = TEMPO_ITS, TEMPO_ITE
        if (column_mp_active(i,j) .and. qs_col_any(i,j)) then
          substeps_sedi_max = max(substeps_sedi_max, substeps_sedi(i,j))
        endif
      enddo
    enddo
    !$acc end parallel

    !! column_j3c: snow Eulerian sedimentation substeps (OpenACC: n outer; (j,i) inner)
    call nvtx_range_push('column_j3c_snow_sedimentation')
    column_j3c_n: do n = 1, substeps_sedi_max

      call sedimentation(kts, kte, its, ite, jts, jte, dz3d, rho, rs, qsten, vtrs, substeps_sedi, ktop_sedi, r1, dt, &
        column_mp_active, qs_col_any, n, precip=tempo_main_diags%snow_liquid_equiv_precip)

    enddo column_j3c_n
    call nvtx_range_pop()

    !! column_j3c2_init: ice sedimentation control defaults (active columns)
    !$acc parallel
    !$acc loop gang vector collapse(2)
    do j = TEMPO_JTS, TEMPO_JTE
      do i = TEMPO_ITS, TEMPO_ITE
        if (column_mp_active(i,j)) then
          ktop_sedi(i,j) = 1
          substeps_sedi(i,j) = 1
        endif
      enddo
    enddo
    !$acc end parallel

    !! preliminary ice fallspeed (sets substeps_sedi, ktop_sedi) before Eulerian sedimentation
    call nvtx_range_push('ice_fallspeed')
    call ice_fallspeed(kts, kte, its, ite, jts, jte, dt, rhof, l_qi, ri, ilami, dz3d, vtri, vtni, substeps_sedi, ktop_sedi, &
      column_mp_active, l_qi_flag=.true.)
    call nvtx_range_pop()

    !! column_j3c2_qi: per-column ice presence before sedimentation (OpenACC; avoids any() on device)
    !$acc parallel
    !$acc loop gang vector collapse(2)
    do j = TEMPO_JTS, TEMPO_JTE
      do i = TEMPO_ITS, TEMPO_ITE
        qi_col_any(i,j) = .false.
        if (column_mp_active(i,j)) then
          !$acc loop seq
          do k = kts, kte
            if (l_qi(k,i,j)) qi_col_any(i,j) = .true.
          enddo
        endif
      enddo
    enddo
    !$acc end parallel

    !! column_j3c2: ice Eulerian sedimentation (OpenACC: one (j,i) pass per routine)
    call nvtx_range_push('column_j3c2_ice_sedimentation')

    !! column_j3c2_ri: Eulerian sedimentation for ice mass
    call sedimentation(kts, kte, its, ite, jts, jte, dz3d, rho, ri, qiten, vtri, substeps_sedi, ktop_sedi, r1, dt, &
      column_mp_active, qi_col_any, precip=tempo_main_diags%ice_liquid_equiv_precip)

    !! column_j3c2_ni: Eulerian sedimentation for ice number
    call sedimentation(kts, kte, its, ite, jts, jte, dz3d, rho, ni, niten, vtni, substeps_sedi, ktop_sedi, r2, dt, &
      column_mp_active, qi_col_any)

    call nvtx_range_pop()

    !! column_j3c3_init: cloud sedimentation control defaults (active columns)
    !$acc parallel
    !$acc loop gang vector collapse(2)
    do j = TEMPO_JTS, TEMPO_JTE
      do i = TEMPO_ITS, TEMPO_ITE
        if (column_mp_active(i,j)) then
          ktop_sedi(i,j) = 1
          substeps_sedi(i,j) = 1
        endif
      enddo
    enddo
    !$acc end parallel

    !! preliminary cloud fallspeed (sets ktop_sedi) before Eulerian sedimentation
    call nvtx_range_push('cloud_fallspeed')
    call cloud_fallspeed(kts, kte, its, ite, jts, jte, rhof, w3d, l_qc, rc, nc, ilamc, dz3d, vtrc, vtnc, ktop_sedi, &
      column_mp_active, l_qc_flag=.true.)
    call nvtx_range_pop()

    !! column_j3c3_qc: per-column cloud presence before sedimentation (OpenACC; avoids any() on device)
    !$acc parallel
    !$acc loop gang vector collapse(2)
    do j = TEMPO_JTS, TEMPO_JTE
      do i = TEMPO_ITS, TEMPO_ITE
        qc_col_any(i,j) = .false.
        if (column_mp_active(i,j)) then
          !$acc loop seq
          do k = kts, kte
            if (l_qc(k,i,j)) qc_col_any(i,j) = .true.
          enddo
        endif
      enddo
    enddo
    !$acc end parallel

    !! column_j3c3: cloud Eulerian sedimentation (OpenACC: one (j,i) pass per routine)
    call nvtx_range_push('column_j3c3_cloud_sedimentation')

    !! column_j3c3_rc: Eulerian sedimentation for cloud mass
    call sedimentation(kts, kte, its, ite, jts, jte, dz3d, rho, rc, qcten, vtrc, substeps_sedi, ktop_sedi, r1, dt, &
      column_mp_active, qc_col_any, precip=tempo_main_diags%cloud_precip)

    !! column_j3c3_nc: Eulerian sedimentation for cloud number
    call sedimentation(kts, kte, its, ite, jts, jte, dz3d, rho, nc, ncten, vtnc, substeps_sedi, ktop_sedi, r2, dt, &
      column_mp_active, qc_col_any)

    call nvtx_range_pop()

    !! after sedimentation: freeze cloud water below hgfrz, melt cloud ice above freezing
    if (.not. tempo_cfgs%turn_off_micro_flag) then
      if (present(nc3d)) then
        call nvtx_range_push('freeze_cloud_melt_ice')
        call freeze_cloud_melt_ice(kts, kte, its, ite, jts, jte, temp, rho, ocp, lvap, qi3d, ni3d, qiten, niten, &
          qc3d, qcten, ncten, tten, dt, odt, column_mp_active, nc3d=nc3d)
        call nvtx_range_pop()
      else
        call nvtx_range_push('freeze_cloud_melt_ice')
        call freeze_cloud_melt_ice(kts, kte, its, ite, jts, jte, temp, rho, ocp, lvap, qi3d, ni3d, qiten, niten, &
          qc3d, qcten, ncten, tten, dt, odt, column_mp_active, ncsave=ncsave)
        call nvtx_range_pop()
      endif
    endif

    !! column_j3e: final host-state update (t3d, qv3d, rho, optional aerosol 3d)
    !$acc parallel
    !$acc loop gang vector collapse(3)
    do j = TEMPO_JTS, TEMPO_JTE
      do i = TEMPO_ITS, TEMPO_ITE
        do k = 1, nz
          if (column_mp_active(i,j)) then
            t3d(k,i,j)  = t3d(k,i,j) + tten(k,i,j)*dt
            qv3d(k,i,j) = max(min_qv, (qv3d(k,i,j) + qvten(k,i,j)*dt))
            rho(k,i,j) = roverrv*pres(k,i,j)/(rdry*temp(k,i,j)*(qv(k,i,j)+roverrv))
            if (present(nwfa3d)) then
              nwfa3d(k,i,j) = max(nwfa_default, min(aero_max, (nwfa3d(k,i,j)+nwfaten(k,i,j)*dt)))
            endif
            if (present(nifa3d)) then
              nifa3d(k,i,j) = max(nifa_default, min(aero_max, (nifa3d(k,i,j)+nifaten(k,i,j)*dt)))
            endif
          endif
        enddo
      enddo
    enddo
    !$acc end parallel

    if (present(nc3d)) then
      call nvtx_range_push('cloud_check_and_update')
      call cloud_check_and_update(dt, odt, kts, kte, its, ite, jts, jte, rho, l_qc, qc3d, rc, nc, qcten, ncten, ilamc, mvd_c, &
        nc3d=nc3d, ncsave=ncsave, column_mp_active=column_mp_active)
      call nvtx_range_pop()
    else
      call nvtx_range_push('cloud_check_and_update')
      call cloud_check_and_update(dt, odt, kts, kte, its, ite, jts, jte, rho, l_qc, qc3d, rc, nc, qcten, ncten, ilamc, mvd_c, &
        ncsave=ncsave, column_mp_active=column_mp_active)
      call nvtx_range_pop()
    endif

    call nvtx_range_push('rain_check_and_update')
    call rain_check_and_update(dt, odt, kts, kte, its, ite, jts, jte, rho, l_qr, qr3d, nr3d, rr, nr, qrten, nrten, ilamr, mvd_r, &
      column_mp_active=column_mp_active)
    call nvtx_range_pop()

    call nvtx_range_push('ice_check_and_update')
    call ice_check_and_update(dt, odt, kts, kte, its, ite, jts, jte, rho, l_qi, qi3d, ni3d, ri, ni, qiten, niten, ilami, &
      column_mp_active=column_mp_active)
    call nvtx_range_pop()
    call nvtx_range_push('snow_check_and_update')
    call snow_check_and_update(dt, odt, kts, kte, its, ite, jts, jte, rho, l_qs, qs3d, rs, qsten, column_mp_active=column_mp_active)
    call nvtx_range_pop()

    !! column_j3e6: snow moment workspace after ice/snow checks (gang over columns, vector over k; snow_moments is acc routine)
    call nvtx_range_push('snow_moments_column_j3e6')
    !$acc parallel
    !$acc loop gang collapse(2)
    column_j3e6: do j = TEMPO_JTS, TEMPO_JTE
      column_i3e6: do i = TEMPO_ITS, TEMPO_ITE
        if (column_mp_active(i,j)) then
          !$acc loop vector private(tc0)
          do k = 1, nz
            smo0(k,i,j) = 0._dp
            smo1(k,i,j) = 0._dp
            smo2(k,i,j) = 0._dp
            smob(k,i,j) = 0._dp
            smoc(k,i,j) = 0._dp
            smoe(k,i,j) = 0._dp
            smof(k,i,j) = 0._dp
            smog(k,i,j) = 0._dp
            smoz(k,i,j) = 0._dp
            ns(k,i,j) = 0._dp
            if (l_qs(k,i,j)) then
              tc0 = min(-0.1, temp(k,i,j)-t0)
              call snow_moments(rs=rs(k,i,j), tc=tc0, &
                smob=smob(k,i,j), smoc=smoc(k,i,j), &
                smo2=smo2(k,i,j), smoz=smoz(k,i,j))
            endif
          enddo
        endif
      enddo column_i3e6
    enddo column_j3e6
    !$acc end parallel
    call nvtx_range_pop()

    if (present(ng3d) .and. present(qb3d)) then
      call nvtx_range_push('graupel_check_and_update')
      call graupel_check_and_update(dt, odt, kts, kte, its, ite, jts, jte, rho, l_qg, qg3d, rg, ng, rb, idx_bg, qgten, ngten, qbten, &
        ilamg, mvd_g, column_mp_active=column_mp_active, ng3d=ng3d, qb3d=qb3d)
      call nvtx_range_pop()
    else
      call nvtx_range_push('graupel_check_and_update')
      call graupel_check_and_update(dt, odt, kts, kte, its, ite, jts, jte, rho, l_qg, qg3d, rg, ng, rb, idx_bg, qgten, ngten, qbten, &
        ilamg, mvd_g, column_mp_active=column_mp_active)
      call nvtx_range_pop()
    endif

    !! column_j3f: frozen_fraction diagnostic from precip components
    !$acc parallel
    !$acc loop gang vector collapse(2)
    do j = TEMPO_JTS, TEMPO_JTE
      do i = TEMPO_ITS, TEMPO_ITE
        if (column_mp_active(i,j)) then
          tempo_main_diags%frozen_fraction(i,j) = &
            (tempo_main_diags%ice_liquid_equiv_precip(i,j) + tempo_main_diags%snow_liquid_equiv_precip(i,j) + &
            tempo_main_diags%graupel_liquid_equiv_precip(i,j)) / &
            (tempo_main_diags%ice_liquid_equiv_precip(i,j) + tempo_main_diags%snow_liquid_equiv_precip(i,j) + &
            tempo_main_diags%graupel_liquid_equiv_precip(i,j) + tempo_main_diags%rain_precip(i,j) + r1)
        endif
      enddo
    enddo
    !$acc end parallel

    call nvtx_range_push('freezing_rain')
    call freezing_rain(kts, kte, its, ite, jts, jte, temp, tempo_main_diags%rain_precip, tempo_main_diags%cloud_precip, &
      tempo_main_diags%frz_rain_precip, column_mp_active)
    call nvtx_range_pop()

    !! column_j3g: cloud_number_mixing_ratio diagnostic (full column)
    !$acc parallel
    !$acc loop gang vector collapse(3)
    do j = TEMPO_JTS, TEMPO_JTE
      do i = TEMPO_ITS, TEMPO_ITE
        do k = 1, nz
          if (column_mp_active(i,j) .and. tempo_cfgs%cloud_number_mixing_ratio_flag) then
            tempo_main_diags%cloud_number_mixing_ratio(k, i, j) = nc(k,i,j)*rho(k,i,j)
          endif
        enddo
      enddo
    enddo
    !$acc end parallel

    !! column_j3g2: rain_med_vol_diam diagnostic (full column)
    !$acc parallel
    !$acc loop gang vector collapse(3)
    do j = TEMPO_JTS, TEMPO_JTE
      do i = TEMPO_ITS, TEMPO_ITE
        do k = 1, nz
          if (column_mp_active(i,j) .and. tempo_cfgs%rain_med_vol_diam_flag) then
            tempo_main_diags%rain_med_vol_diam(k, i, j) = mvd_r(k,i,j)
          endif
        enddo
      enddo
    enddo
    !$acc end parallel

    !! column_j3g3: graupel_med_vol_diam diagnostic (full column)
    !$acc parallel
    !$acc loop gang vector collapse(3)
    do j = TEMPO_JTS, TEMPO_JTE
      do i = TEMPO_ITS, TEMPO_ITE
        do k = 1, nz
          if (column_mp_active(i,j) .and. tempo_cfgs%graupel_med_vol_diam_flag) then
            tempo_main_diags%graupel_med_vol_diam(k, i, j) = mvd_g(k,i,j)
          endif
        enddo
      enddo
    enddo
    !$acc end parallel

    if (tempo_cfgs%max_hail_diameter_flag) then
      call nvtx_range_push('max_hail_diam')
      call max_hail_diam(kts, kte, its, ite, jts, jte, rho, rg, ng, ilamg, idx_bg, tempo_main_diags%max_hail_diameter, &
        column_mp_active)
      call nvtx_range_pop()
    endif

    if (tempo_cfgs%refl10cm_flag) then
      call nvtx_range_push('reflectivity_10cm')
      call reflectivity_10cm(kts, kte, its, ite, jts, jte, tempo_cfgs%refl10cm_from_melting_flag, &
        temp, l_qr, rr, nr, ilamr, l_qs, rs, smoc, smob, smoz, l_qg, rg, ng, idx_bg, ilamg, &
        tempo_main_diags%refl10cm, column_mp_active)
      call nvtx_range_pop()
    endif

    !! effective-radius diagnostics: BL / cloud_check prep, then effective_radius (OpenACC: parallel over columns)
    call nvtx_range_push('column_j3g6a_effective_radius_prep')
    !$acc parallel
    !$acc loop gang collapse(2) private(k)
    column_j3g6a: do j = TEMPO_JTS, TEMPO_JTE
      column_i3g6a: do i = TEMPO_ITS, TEMPO_ITE
        if (column_mp_active(i,j)) then
          if ((tempo_cfgs%re_cloud_flag) .and. (tempo_cfgs%re_ice_flag) .and. (tempo_cfgs%re_snow_flag)) then
            if (present(qc_bl3d) .and. present(qcfrac_bl3d) .and. present(nc3d)) then
              xrx(kts:kte,i,j) = qc3d(kts:kte,i,j)
              xnx(kts:kte,i,j) = nc3d(kts:kte,i,j)
              !$acc loop seq
              do k = kts, kte
                if ((xrx(k,i,j) <= r1) .and. &
                  (qc_bl3d(k,i,j) > 1.e-9_wp) .and. (qcfrac_bl3d(k,i,j) > 0._wp)) then
                  xrx(k,i,j) = xrx(k,i,j) + qc_bl3d(k,i,j) / qcfrac_bl3d(k,i,j)
                endif
              enddo
              !$acc loop seq
              do k = kts, kte
                if (xrx(k,i,j) <= 1.e-12_wp) xrx(k,i,j) = 0._wp
              enddo
              call get_cloud_number(nz, xrx(kts:kte,i,j), qr3d(kts:kte,i,j), qi3d(kts:kte,i,j), qs3d(kts:kte,i,j), pres(kts:kte,i,j), &
                temp(kts:kte,i,j), w3d(kts:kte,i,j), xnx(kts:kte,i,j))
              qcten(kts:kte,i,j) = 0._wp
              ncten(kts:kte,i,j) = 0._wp
            endif
          endif
        endif
      enddo column_i3g6a
    enddo column_j3g6a
    !$acc end parallel
    call nvtx_range_pop()

    if ((tempo_cfgs%re_cloud_flag) .and. (tempo_cfgs%re_ice_flag) .and. (tempo_cfgs%re_snow_flag)) then
      if (present(qc_bl3d) .and. present(qcfrac_bl3d) .and. present(nc3d)) then
        call cloud_check_and_update(dt, odt, kts, kte, its, ite, jts, jte, rho, l_qc, xrx, rc, nc, qcten, ncten, ilamc, mvd_c, &
          nc3d=xnx, column_mp_active=column_mp_active, update_qc_nc_state=.false.)
      endif
    endif

    if ((tempo_cfgs%re_cloud_flag) .and. (tempo_cfgs%re_ice_flag) .and. (tempo_cfgs%re_snow_flag)) then
      call nvtx_range_push('effective_radius')
      call effective_radius(kts, kte, its, ite, jts, jte, temp, l_qc, nc, ilamc, l_qi, ilami, l_qs, rs, &
        tempo_main_diags%re_cloud, tempo_main_diags%re_ice, tempo_main_diags%re_snow, column_mp_active)
      call nvtx_range_pop()
    endif

    if (allocated(ncsave)) then
      !$acc exit data delete(ncsave)
      deallocate(ncsave)
    endif

    !! Close structured data region for stack-allocated tile workspace (opened just after
    !! the associate / before column_j1a_tile_init).
    !$acc end data

    end associate

    call nvtx_range_pop()
  end subroutine tempo_main


  subroutine aerosol_check_and_update(dt, kts, kte, its, ite, jts, jte, rho, nwfa, nifa, nwfaten, nifaten, nwfa3d, nifa3d, &
      column_mp_active)
    !! sets aerosol number concentrations and checks bounds over the horizontal tile
    use module_mp_tempo_params, only : nwfa_default, aero_max, nifa_default

    integer, intent(in) :: kts, kte, its, ite, jts, jte
    real(wp), intent(in) :: dt
    real(wp), dimension(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(in) :: rho, nwfaten, nifaten
    real(wp), dimension(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(inout) :: nwfa, nifa
    real(wp), dimension(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(in), optional :: nwfa3d, nifa3d
    logical, dimension(TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(in), optional :: column_mp_active
    integer :: i, j, k
    logical :: use_cmp, use_nwfa3d, use_nifa3d

    use_cmp = present(column_mp_active)
    use_nwfa3d = present(nwfa3d)
    use_nifa3d = present(nifa3d)

    !$acc parallel
    !$acc loop gang vector collapse(3)
    do j = TEMPO_JTS, TEMPO_JTE
      do i = TEMPO_ITS, TEMPO_ITE
        do k = kts, kte
          if (.not. use_cmp .or. column_mp_active(i,j)) then
            if (use_nwfa3d) then
              nwfa(k,i,j) = (nwfa3d(k,i,j)+nwfaten(k,i,j)*dt)*rho(k,i,j)
            endif
            if (use_nifa3d) then
              nifa(k,i,j) = (nifa3d(k,i,j)+nifaten(k,i,j)*dt)*rho(k,i,j)
            endif
            nwfa(k,i,j) = max(nwfa_default*rho(k,i,j), min(aero_max*rho(k,i,j), nwfa(k,i,j)))
            nifa(k,i,j) = max(nifa_default*rho(k,i,j), min(aero_max*rho(k,i,j), nifa(k,i,j)))
          endif
        enddo
      enddo
    enddo
    !$acc end parallel
  end subroutine aerosol_check_and_update


  subroutine cloud_check_and_update(dt, odt, kts, kte, its, ite, jts, jte, rho, l_qc, qc3d, rc, nc, &
      qcten, ncten, ilamc, mvd_c, nc3d, ncsave, column_mp_active, update_qc_nc_state)
    !! Cloud mass/size checks over k and the horizontal tile [TEMPO_ITS:TEMPO_ITE] x [TEMPO_JTS:TEMPO_JTE] (cf. rain_check_and_update).
    !! When update_qc_nc_state is .true. (default), working qc/nc mixing ratios are written back to qc3d / nc3d (two-moment).
    !! When .false., qc3d/nc3d are read but left unchanged (internal work copy). Optional nc3d selects two-moment cloud number.
    use module_mp_tempo_params, only : r1, nt_c_max, nt_c_min, nu_c_scale, &
      am_r, bm_r, cce, ccg, d0c, d0r, ocg1, ocg2, obmr, nt_c_l, d0r, nt_c_l

    integer, intent(in) :: kts, kte, its, ite, jts, jte
    real(wp), intent(in) :: dt, odt
    real(wp), dimension(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(in) :: rho
    logical, dimension(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(inout) :: l_qc
    real(wp), dimension(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(inout) :: qc3d, qcten, ncten, rc, nc
    real(wp), dimension(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(out) :: mvd_c
    real(dp), dimension(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(out) :: ilamc
    real(wp), dimension(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(inout), optional :: nc3d
    real(wp), dimension(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(in), optional :: ncsave
    logical, dimension(TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(in), optional :: column_mp_active
    logical, intent(in), optional :: update_qc_nc_state
    integer :: i, j, k, nu_c
    real(dp) :: lamc, xdc
    logical :: hit_limit
    logical :: write_qc_nc
    logical :: use_cmp, use_nc3d, use_ncsave
    real(wp) :: qcwork(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), ncwork(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE)

    write_qc_nc = .true.
    if (present(update_qc_nc_state)) write_qc_nc = update_qc_nc_state
    use_cmp = present(column_mp_active)
    use_nc3d = present(nc3d)
    use_ncsave = present(ncsave)

    !$acc data create(qcwork, ncwork)
    !$acc parallel
    !$acc loop gang vector collapse(2) private(lamc, xdc, hit_limit, nu_c)
    do j = TEMPO_JTS, TEMPO_JTE
      do i = TEMPO_ITS, TEMPO_ITE
        if (use_cmp) then
          if (.not. column_mp_active(i, j)) cycle
        endif

        qcwork(kts:kte, i, j) = qc3d(kts:kte, i, j)
        if (use_nc3d) then
          ncwork(kts:kte, i, j) = nc3d(kts:kte, i, j)
        endif

        !$acc loop seq
        do k = kts, kte
          hit_limit = .false.
          if (qcwork(k, i, j)+qcten(k, i, j)*dt > r1) then
            l_qc(k, i, j) = .true.
            rc(k, i, j) = (qcwork(k, i, j)+qcten(k, i, j)*dt)*rho(k, i, j)
            qcwork(k, i, j) = qcwork(k, i, j)+qcten(k, i, j)*dt

            if (use_nc3d) then
              nc(k, i, j) = max(nt_c_min, (ncwork(k, i, j)+ncten(k, i, j)*dt)*rho(k, i, j))

              if (nc(k, i, j) <= nt_c_min) then
                hit_limit = .true.
                nc(k, i, j) = nt_c_min
              endif
              if (nc(k, i, j) > nt_c_max) then
                hit_limit = .true.
                nc(k, i, j) = nt_c_max
              endif

              nu_c = get_nuc(nc(k, i, j))
              lamc = (nc(k, i, j)*am_r*ccg(2,nu_c)*ocg1(nu_c)/rc(k, i, j))**obmr
              xdc = (bm_r + nu_c + 1._dp) / lamc
              if (xdc < d0c) then
                lamc = cce(2,nu_c)/d0c
                hit_limit = .true.
              elseif (xDc > d0r*2._dp) then
                lamc = cce(2,nu_c)/(d0r*2._dp)
                hit_limit = .true.
              endif
              nc(k, i, j) = ccg(1,nu_c)*ocg2(nu_c)*rc(k, i, j) / am_r*lamc**bm_r

              if (hit_limit) ncten(k, i, j) = (nc(k, i, j)/rho(k, i, j) - ncwork(k, i, j)) * odt
              ncwork(k, i, j) = max(nt_c_min/rho(k, i, j), &
                min(ccg(1,nu_c)*ocg2(nu_c)*qcwork(k, i, j)/am_r*lamc**bm_r, nt_c_max/rho(k, i, j)))
            else
              if (use_ncsave) then
                nc(k, i, j) = ncsave(k, i, j)
              else
                nc(k, i, j) = nt_c_l
              endif
            endif
            nu_c = get_nuc(nc(k, i, j))
            lamc = (nc(k, i, j)*am_r*ccg(2,nu_c)*ocg1(nu_c)/rc(k, i, j))**obmr
            ilamc(k, i, j) = 1._dp / lamc
            mvd_c(k, i, j) = max(min((3.0_wp + nu_c + 0.672_wp) * ilamc(k, i, j), d0r), d0c)
          else
            l_qc(k, i, j) = .false.
            rc(k, i, j) = r1
            nc(k, i, j) = nt_c_min
            mvd_c(k, i, j) = d0c
            ilamc(k, i, j) = 0._dp
            qcten(k, i, j) = -qcwork(k, i, j) * odt
            qcwork(k, i, j) = 0.0_wp
            if (use_nc3d) then
              ncten(k, i, j) = -ncwork(k, i, j) * odt
              ncwork(k, i, j) = 0.0_wp
            endif
          endif
        enddo

        if (write_qc_nc) then
          qc3d(kts:kte, i, j) = qcwork(kts:kte, i, j)
          if (use_nc3d) nc3d(kts:kte, i, j) = ncwork(kts:kte, i, j)
        endif
      enddo
    enddo
    !$acc end parallel
    !$acc end data
  end subroutine cloud_check_and_update


  subroutine rain_check_and_update(dt, odt, kts, kte, its, ite, jts, jte, rho, l_qr, qr3d, nr3d, rr, nr, qrten, nrten, ilamr, mvd_r, &
      column_mp_active, update_qr_nr_state, substep_mode, col_any, n, substeps_sedi)
    !! Rain mass/size checks over k and the horizontal bounds [TEMPO_ITS:TEMPO_ITE] x [TEMPO_JTS:TEMPO_JTE] (cf. aerosol_check_and_update).
    !! When update_qr_nr_state is .true. (default), working rain mixing ratios / numbers are written back to qr3d/nr3d.
    !! When .false., qr3d/nr3d are read but left unchanged. Substep mode additionally gates work using col_any, n, and substeps_sedi.
    integer, intent(in) :: kts, kte, its, ite, jts, jte
    real(wp), intent(in) :: dt, odt
    real(wp), dimension(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(in) :: rho
    logical, dimension(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(inout) :: l_qr
    real(wp), dimension(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(inout) :: qr3d, nr3d, qrten, nrten
    real(wp), dimension(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(inout) :: rr, nr
    real(dp), dimension(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(inout) :: ilamr
    real(wp), dimension(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(inout) :: mvd_r
    logical, dimension(TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(in), optional :: column_mp_active, col_any
    integer, dimension(TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(in), optional :: substeps_sedi
    logical, intent(in), optional :: update_qr_nr_state, substep_mode
    integer, intent(in), optional :: n
    integer :: i, j, k
    real(dp) :: lamr
    logical :: hit_limit
    logical :: write_qr_nr
    logical :: use_cmp, use_substep
    real(wp) :: qwork(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), nwork(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE)

    write_qr_nr = .true.
    if (present(update_qr_nr_state)) write_qr_nr = update_qr_nr_state
    use_cmp = present(column_mp_active)
    use_substep = .false.
    if (present(substep_mode)) use_substep = substep_mode

    if (use_substep) then
      if (.not. present(column_mp_active) .or. .not. present(col_any) .or. &
          .not. present(n) .or. .not. present(substeps_sedi)) then
        error stop "rain_check_and_update: substep mode requires column_mp_active, col_any, n, and substeps_sedi"
      endif
    endif

    !$acc data create(qwork, nwork)
    !$acc parallel
    !$acc loop gang vector collapse(2) private(lamr, hit_limit)
    do j = TEMPO_JTS, TEMPO_JTE
      do i = TEMPO_ITS, TEMPO_ITE
        if (use_cmp) then
          if (.not. column_mp_active(i, j)) cycle
        endif
        if (use_substep) then
          if (.not. col_any(i, j)) cycle
          if (n > substeps_sedi(i, j)) cycle
        endif

        qwork(kts:kte, i, j) = qr3d(kts:kte, i, j)
        nwork(kts:kte, i, j) = nr3d(kts:kte, i, j)

        !$acc loop seq
        do k = kts, kte
          hit_limit = .false.
          if (qwork(k, i, j)+qrten(k, i, j)*dt > r1) then
            l_qr(k, i, j) = .true.
            rr(k, i, j) = (qwork(k, i, j)+qrten(k, i, j)*dt)*rho(k, i, j)
            qwork(k, i, j) = qwork(k, i, j)+qrten(k, i, j)*dt

            nr(k, i, j) = max(r2, (nwork(k, i, j)+nrten(k, i, j)*dt)*rho(k, i, j))
            if (nr(k, i, j) <= r2) then
              hit_limit = .true.
              mvd_r(k, i, j) = 1.0e-3_wp
              lamr = (3.0_dp + mu_r + 0.672_dp) / mvd_r(k, i, j)
              nr(k, i, j) = crg(2)*org3*rr(k, i, j)*lamr**bm_r / am_r
            endif

            lamr = (am_r*crg(3)*org2*nr(k, i, j)/rr(k, i, j))**obmr
            mvd_r(k, i, j) = (3.0_wp + mu_r + 0.672_wp) / lamr
            if (mvd_r(k, i, j) > d0r_max) then
              hit_limit = .true.
              mvd_r(k, i, j) = d0r_max
              lamr = (3.0_dp + mu_r + 0.672_dp) / mvd_r(k, i, j)
              nr(k, i, j) = crg(2)*org3*rr(k, i, j)*lamr**bm_r / am_r
            elseif (mvd_r(k, i, j) < d0r*0.75_wp) then
              hit_limit = .true.
              mvd_r(k, i, j) = d0r*0.75_wp
              lamr = (3.0_dp + mu_r + 0.672_dp) / mvd_r(k, i, j)
              nr(k, i, j) = crg(2)*org3*rr(k, i, j)*lamr**bm_r / am_r
            endif
            if (hit_limit) nrten(k, i, j) = (nr(k, i, j)/rho(k, i, j) - nwork(k, i, j))*odt
            nwork(k, i, j) = crg(2)*org3*qwork(k, i, j)*lamr**bm_r / am_r
            ilamr(k, i, j) = 1._dp / lamr
          else
            l_qr(k, i, j) = .false.
            rr(k, i, j) = r1
            nr(k, i, j) = r2
            mvd_r(k, i, j) = d0r
            ilamr(k, i, j) = 0._dp
            qrten(k, i, j) = -qwork(k, i, j) * odt
            nrten(k, i, j) = -nwork(k, i, j) * odt
            qwork(k, i, j) = 0.0_wp
            nwork(k, i, j) = 0.0_wp
          endif
        enddo

        if (write_qr_nr) then
          qr3d(kts:kte, i, j) = qwork(kts:kte, i, j)
          nr3d(kts:kte, i, j) = nwork(kts:kte, i, j)
        endif
      enddo
    enddo
    !$acc end parallel
    !$acc end data
  end subroutine rain_check_and_update


  subroutine ice_check_and_update(dt, odt, kts, kte, its, ite, jts, jte, rho, l_qi, qi3d, ni3d, ri, ni, qiten, niten, ilami, &
      column_mp_active, update_qi_ni_state)
    !! Ice mass/size checks over k and the horizontal tile (cf. rain_check_and_update).
    !! When update_qi_ni_state is .true. (default), working values are written back to qi3d/ni3d. When .false., they are unchanged.
    use module_mp_tempo_params, only : max_ni, r1, r2, cie, cig, &
      mu_i, am_i, bm_i, oig1, oig2, obmi , d0s

    integer, intent(in) :: kts, kte, its, ite, jts, jte
    real(wp), intent(in) :: dt, odt
    real(wp), dimension(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(in) :: rho
    logical, dimension(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(inout) :: l_qi
    real(wp), dimension(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(inout) :: qi3d, ni3d, qiten, niten
    real(wp), dimension(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(inout) :: ri, ni
    real(dp), dimension(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(inout) :: ilami
    logical, dimension(TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(in), optional :: column_mp_active
    logical, intent(in), optional :: update_qi_ni_state
    integer :: i, j, k
    real(dp) :: lami, xdi
    logical :: hit_limit
    logical :: write_qi_ni
    logical :: use_cmp
    real(wp) :: qiwork(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), niwork(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE)

    write_qi_ni = .true.
    if (present(update_qi_ni_state)) write_qi_ni = update_qi_ni_state
    use_cmp = present(column_mp_active)

    !$acc data create(qiwork, niwork)
    !$acc parallel
    !$acc loop gang vector collapse(2) private(lami, xdi, hit_limit)
    do j = TEMPO_JTS, TEMPO_JTE
      do i = TEMPO_ITS, TEMPO_ITE
        if (use_cmp) then
          if (.not. column_mp_active(i, j)) cycle
        endif

        qiwork(kts:kte, i, j) = qi3d(kts:kte, i, j)
        niwork(kts:kte, i, j) = ni3d(kts:kte, i, j)

        !$acc loop seq
        do k = kts, kte
          hit_limit = .false.
          if (qiwork(k, i, j)+qiten(k, i, j)*dt > r1) then
            l_qi(k, i, j) = .true.
            ri(k, i, j) = (qiwork(k, i, j)+qiten(k, i, j)*dt)*rho(k, i, j)
            qiwork(k, i, j) = qiwork(k, i, j)+qiten(k, i, j)*dt

            ni(k, i, j) = max(r2, (niwork(k, i, j)+niten(k, i, j)*dt)*rho(k, i, j))

            if (ni(k, i, j) <= r2) then
              hit_limit = .true.
              lami = cie(2)/5.e-6_dp
              ni(k, i, j) = min(max_ni, cig(1)*oig2*ri(k, i, j)/am_i*lami**bm_i)
            endif

            lami = (am_i*cig(2)*oig1*ni(k, i, j)/ri(k, i, j))**obmi
            xdi = (bm_i + mu_i + 1._dp) / lami
            if (xdi < 5.e-6_dp) then
              hit_limit = .true.
              lami = cie(2)/5.e-6_dp
              ni(k, i, j) = min(max_ni, cig(1)*oig2*ri(k, i, j)/am_i*lami**bm_i)
            elseif (xdi > d0s) then
              hit_limit = .true.
              lami = cie(2)/d0s
              ni(k, i, j) = cig(1)*oig2*ri(k, i, j)/am_i*lami**bm_i
            endif

            if (hit_limit) niten(k, i, j) = (ni(k, i, j)/rho(k, i, j) - niwork(k, i, j))*odt
            niwork(k, i, j) = max(r2/rho(k, i, j), &
              min(cig(1)*oig2*qiwork(k, i, j)/am_i*lami**bm_i, max_ni/rho(k, i, j)))
            ilami(k, i, j) = 1._dp / lami
          else
            l_qi(k, i, j) = .false.
            ri(k, i, j) = r1
            ni(k, i, j) = r2
            ilami(k, i, j) = 0._dp
            qiten(k, i, j) = -qiwork(k, i, j) * odt
            niten(k, i, j) = -niwork(k, i, j) * odt
            qiwork(k, i, j) = 0.0_wp
            niwork(k, i, j) = 0.0_wp
          endif
        enddo

        if (write_qi_ni) then
          qi3d(kts:kte, i, j) = qiwork(kts:kte, i, j)
          ni3d(kts:kte, i, j) = niwork(kts:kte, i, j)
        endif
      enddo
    enddo
    !$acc end parallel
    !$acc end data
  end subroutine ice_check_and_update


  subroutine snow_check_and_update(dt, odt, kts, kte, its, ite, jts, jte, rho, l_qs, qs3d, rs, qsten, column_mp_active, update_qs_state)
    !! Snow mass checks over k and the horizontal tile (cf. rain_check_and_update).
    !! When update_qs_state is .true. (default), working values are written back to qs3d. When .false., qs3d is read but unchanged.
    use module_mp_tempo_params, only : max_ni, r1, r2

    integer, intent(in) :: kts, kte, its, ite, jts, jte
    real(wp), intent(in) :: dt, odt
    real(wp), dimension(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(in) :: rho
    logical, dimension(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(inout) :: l_qs
    real(wp), dimension(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(inout) :: qs3d, rs, qsten
    logical, dimension(TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(in), optional :: column_mp_active
    logical, intent(in), optional :: update_qs_state
    integer :: i, j, k
    logical :: write_qs
    logical :: use_cmp
    real(wp) :: qswork(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE)

    write_qs = .true.
    if (present(update_qs_state)) write_qs = update_qs_state
    use_cmp = present(column_mp_active)

    !$acc data create(qswork)
    !$acc parallel
    !$acc loop gang vector collapse(2)
    do j = TEMPO_JTS, TEMPO_JTE
      do i = TEMPO_ITS, TEMPO_ITE
        if (use_cmp) then
          if (.not. column_mp_active(i, j)) cycle
        endif

        qswork(kts:kte, i, j) = qs3d(kts:kte, i, j)

        !$acc loop seq
        do k = kts, kte
          if (qswork(k, i, j)+qsten(k, i, j)*dt > r1) then
            l_qs(k, i, j) = .true.
            rs(k, i, j) = (qswork(k, i, j)+qsten(k, i, j)*dt)*rho(k, i, j)
            qswork(k, i, j) = qswork(k, i, j)+qsten(k, i, j)*dt
          else
            l_qs(k, i, j) = .false.
            rs(k, i, j) = r1
            qsten(k, i, j) = -qswork(k, i, j) * odt
            qswork(k, i, j) = 0.0_wp
          endif
        enddo

        if (write_qs) qs3d(kts:kte, i, j) = qswork(kts:kte, i, j)
      enddo
    enddo
    !$acc end parallel
    !$acc end data
  end subroutine snow_check_and_update


  subroutine graupel_check_and_update(dt, odt, kts, kte, its, ite, jts, jte, rho, l_qg, qg3d, rg, ng, rb, idx_bg, &
      qgten, ngten, qbten, ilamg, mvd_g, column_mp_active, update_qg_ng_qb_state, ng3d, qb3d, substep_mode, col_any, n, &
      substeps_sedi)
    !! Graupel mass/size checks over k and the horizontal bounds [TEMPO_ITS:TEMPO_ITE] x [TEMPO_JTS:TEMPO_JTE] (cf. rain_check_and_update).
    !! When update_qg_ng_qb_state is .true. (default), qg3d/ng3d/qb3d mixing ratios are written back from work copies.
    !! When .false., those prognostics are unchanged (internal work copy only).
    !! Substep mode additionally gates work using col_any, n, and substeps_sedi.
    !! OpenACC: parallel over columns; vertical k is sequential.
    use module_mp_tempo_params, only : r1, r2, nrhg, rho_g, mu_g, &
      am_g, bm_g, ogg3, cgg, ogg2, obmg, d0r, idx_bg1, gonv_max, &
      gonv_min, oge1, ogg1, d0g, meters3_to_liters

    integer, intent(in) :: kts, kte, its, ite, jts, jte
    real(wp), intent(in) :: dt, odt
    real(wp), dimension(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(in) :: rho
    logical, dimension(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(inout) :: l_qg
    real(wp), dimension(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(inout) :: qg3d, qgten, rg, ng, rb, ngten, qbten
    real(wp), dimension(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(inout), optional :: ng3d, qb3d
    real(dp), dimension(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(inout) :: ilamg
    real(wp), dimension(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(inout) :: mvd_g
    integer, dimension(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(inout) :: idx_bg
    logical, dimension(TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(in), optional :: column_mp_active, col_any
    integer, dimension(TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(in), optional :: substeps_sedi
    logical, intent(in), optional :: update_qg_ng_qb_state, substep_mode
    integer, intent(in), optional :: n
    integer :: i, j, k
    real(dp) :: lamg, ygra1, zans1, n0_exp, lam_exp
    logical :: hit_limit
    logical :: write_qg_ng_qb
    logical :: use_cmp, use_ng_qb, use_substep
    real(wp) :: qwork(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE)
    real(wp) :: nwork(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE)
    real(wp) :: qbwork(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE)

    write_qg_ng_qb = .true.
    if (present(update_qg_ng_qb_state)) write_qg_ng_qb = update_qg_ng_qb_state
    use_cmp = present(column_mp_active)
    use_ng_qb = present(ng3d) .and. present(qb3d)
    use_substep = .false.
    if (present(substep_mode)) use_substep = substep_mode

    if (use_substep) then
      if (.not. present(column_mp_active) .or. .not. present(col_any) .or. &
          .not. present(n) .or. .not. present(substeps_sedi)) then
        error stop "graupel_check_and_update: substep mode requires column_mp_active, col_any, n, and substeps_sedi"
      endif
    endif

    !$acc data create(qwork, nwork, qbwork)
    !$acc parallel
    !$acc loop gang vector collapse(2) private(lamg, ygra1, zans1, n0_exp, lam_exp, hit_limit)
    do j = TEMPO_JTS, TEMPO_JTE
      do i = TEMPO_ITS, TEMPO_ITE
        if (use_cmp) then
          if (.not. column_mp_active(i, j)) cycle
        endif
        if (use_substep) then
          if (.not. col_any(i, j)) cycle
          if (n > substeps_sedi(i, j)) cycle
        endif

        qwork(kts:kte, i, j) = qg3d(kts:kte, i, j)
        if (use_ng_qb) then
          nwork(kts:kte, i, j) = ng3d(kts:kte, i, j)
          qbwork(kts:kte, i, j) = qb3d(kts:kte, i, j)
        endif

        !$acc loop seq
        do k = kts, kte
          hit_limit = .false.
          if (qwork(k, i, j)+qgten(k, i, j)*dt > r1) then
            l_qg(k, i, j) = .true.
            rg(k, i, j) = (qwork(k, i, j)+qgten(k, i, j)*dt)*rho(k, i, j)
            qwork(k, i, j) = qwork(k, i, j)+qgten(k, i, j)*dt

            if (use_ng_qb) then
              ng(k, i, j) = max(r2, (nwork(k, i, j)+ngten(k, i, j)*dt)*rho(k, i, j))
              rb(k, i, j) = min(max(rg(k, i, j)*meters3_to_liters/rho_g(nrhg), &
                (qbwork(k, i, j)+qbten(k, i, j)*dt)*rho(k, i, j)), rg(k, i, j)*meters3_to_liters/rho_g(1))
              idx_bg(k, i, j) = max(1, min(nint(10._wp*rg(k, i, j)/rb(k, i, j))+1, nrhg))

              if (ng(k, i, j) <= r2) then
                hit_limit = .true.
                mvd_g(k, i, j) = 1.5e-3_wp
                lamg = (3.0_dp + mu_g + 0.672_dp) / mvd_g(k, i, j)
                ng(k, i, j) = cgg(2,1)*ogg3*rg(k, i, j)*lamg**bm_g / am_g(idx_bg(k, i, j))
              endif

              lamg = (am_g(idx_bg(k, i, j))*cgg(3,1)*ogg2*ng(k, i, j)/rg(k, i, j))**obmg
              mvd_g(k, i, j) = (3.0_wp + mu_g + 0.672_wp) / lamg
              if (mvd_g(k, i, j) > 25.4e-3_wp) then
                hit_limit = .true.
                mvd_g(k, i, j) = 25.4e-3_wp
                lamg = (3.0_dp + mu_g + 0.672_dp) / mvd_g(k, i, j)
                ng(k, i, j) = cgg(2,1)*ogg3*rg(k, i, j)*lamg**bm_g / am_g(idx_bg(k, i, j))
              elseif (mvd_g(k, i, j) < d0r) then
                hit_limit = .true.
                mvd_g(k, i, j) = d0r
                lamg = (3.0_dp + mu_g + 0.672_dp) / mvd_g(k, i, j)
                ng(k, i, j) = cgg(2,1)*ogg3*rg(k, i, j)*lamg**bm_g / am_g(idx_bg(k, i, j))
              endif

              if (hit_limit) ngten(k, i, j) = (ng(k, i, j)/rho(k, i, j) - nwork(k, i, j)) * odt
              nwork(k, i, j) = cgg(2,1)*ogg3*qwork(k, i, j)*lamg**bm_g / am_g(idx_bg(k, i, j))
              qbwork(k, i, j) = min(max(qwork(k, i, j)*meters3_to_liters/rho_g(nrhg), &
                qbwork(k, i, j)+qbten(k, i, j)*dt), meters3_to_liters*qwork(k, i, j)/rho_g(1))
              idx_bg(k, i, j) = max(1, min(nint(10._wp*qwork(k, i, j)/qbwork(k, i, j))+1, nrhg))
            else
              idx_bg(k, i, j) = idx_bg1
              ygra1 = log10(max(1.e-9_dp, real(rg(k, i, j), kind=dp)))
              zans1 = 3.4_dp + 2._dp/7._dp*(ygra1+8._dp)
              n0_exp = max(gonv_min, min(10._dp**(zans1), gonv_max))
              lam_exp = (n0_exp*am_g(idx_bg(k, i, j))*cgg(1,1)/rg(k, i, j))**oge1
              lamg = lam_exp * (cgg(3,1)*ogg2*ogg1)**obmg
              ng(k, i, j) = cgg(2,1)*ogg3*rg(k, i, j)*lamg**bm_g / am_g(idx_bg(k, i, j))
              rb(k, i, j) = meters3_to_liters*rg(k, i, j)/rho_g(idx_bg(k, i, j))
            endif
            ilamg(k, i, j) = 1._dp / lamg
            mvd_g(k, i, j) = (3.0_wp + mu_g + 0.672_wp) * ilamg(k, i, j)
          else
            l_qg(k, i, j) = .false.
            rg(k, i, j) = r1
            ng(k, i, j) = r2
            mvd_g(k, i, j) = d0g
            ilamg(k, i, j) = 0._dp
            idx_bg(k, i, j) = idx_bg1
            rb(k, i, j) = meters3_to_liters*r1/rho_g(idx_bg(k, i, j))
            qgten(k, i, j) = -qwork(k, i, j) * odt
            qwork(k, i, j) = 0.0_wp
            if (use_ng_qb) then
              ngten(k, i, j) = -nwork(k, i, j) * odt
              qbten(k, i, j) = -qbwork(k, i, j) * odt
              nwork(k, i, j) = 0.0_wp
              qbwork(k, i, j) = 0.0_wp
            endif
          endif
        enddo

        if (write_qg_ng_qb) then
          qg3d(kts:kte, i, j) = qwork(kts:kte, i, j)
          if (use_ng_qb) then
            ng3d(kts:kte, i, j) = nwork(kts:kte, i, j)
            qb3d(kts:kte, i, j) = qbwork(kts:kte, i, j)
          endif
        endif
      enddo
    enddo
    !$acc end parallel
    !$acc end data

  end subroutine graupel_check_and_update


  subroutine graupel_init(kts, kte, rho, qg1d, ng1d, qb1d)
  !$acc routine vector
    !! initializes graupel number and volume if both are zero
    !! and hail-aware = true
    use module_mp_tempo_params, only : r1, meters3_to_liters, &
      idx_bg1, am_g, bm_g, mu_g, ogg3, cgg, rho_g, nrhg

    integer, intent(in) :: kts, kte
    real(wp), dimension(kts:kte), intent(in) :: rho
    real(wp), dimension(kts:kte), intent(in) :: qg1d
    real(wp), dimension(kts:kte), intent(inout) :: ng1d, qb1d
    integer :: k, idx
    real(dp) :: lamg
    real(wp) :: mvd_g, rg, ng, rb

    !$acc loop vector private(rg, rb, idx, mvd_g, lamg, ng)
    do k = kts, kte
      if (qg1d(k) > r1) then
        rg = qg1d(k)*rho(k)
        rb = rg*meters3_to_liters/rho_g(idx_bg1)
        idx = max(1, min(nint(10._wp*rg/rb)+1, nrhg))
        mvd_g = 5.e-3_wp
        lamg = (3.0_dp + mu_g + 0.672_dp) / mvd_g
        ng = cgg(2,1)*ogg3*rg*lamg**bm_g / am_g(idx)
        ng1d(k) = ng/rho(k)
        qb1d(k) = rb/rho(k)
      endif
    enddo 
  end subroutine graupel_init


  subroutine thermo_vars(kts, kte, its, ite, jts, jte, qv, temp, pres, rho, rhof, rhof2, qvs, delqvs, qvsi, &
      satw, sati, ssatw, ssati, diffu, visco, vsc2, ocp, lvap, tcond, lvt2, column_supersaturated, column_mp_active)
    !! Thermodynamic auxiliary fields over the horizontal tile (cf. aerosol_check_and_update).
    use module_mp_tempo_params, only : t0, rho_not, eps, cp, lvap0, orv

    integer, intent(in) :: kts, kte, its, ite, jts, jte
    real(wp), dimension(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(in) :: qv, temp, pres, rho
    real(wp), dimension(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(out) :: rhof, rhof2, qvs, &
      delqvs, qvsi, satw, sati, ssatw, ssati, diffu, visco, vsc2, ocp, lvap, tcond, lvt2
    logical, dimension(TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(inout) :: column_supersaturated
    logical, dimension(TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(in), optional :: column_mp_active
    integer :: i, j, k
    real(wp) :: tempc, otemp
    logical :: active_col

    !$acc parallel
    !$acc loop gang collapse(2) private(tempc, otemp, active_col)
    do j = TEMPO_JTS, TEMPO_JTE
      do i = TEMPO_ITS, TEMPO_ITE
        active_col = .true.
        if (present(column_mp_active)) active_col = column_mp_active(i, j)
        if (active_col) then
          column_supersaturated(i, j) = .false.
          !$acc loop vector
          do k = kts, kte
            otemp = 1._wp / temp(k, i, j)
            tempc = temp(k, i, j) - t0
            rhof(k, i, j) = sqrt(rho_not/rho(k, i, j))
            rhof2(k, i, j) = sqrt(rhof(k, i, j))
            qvs(k, i, j) = calc_rslf(pres(k, i, j), temp(k, i, j))
            delqvs(k, i, j) = max(0._wp, calc_rslf(pres(k, i, j), t0)-qv(k, i, j))
            if (tempc <= 0._wp) then
              qvsi(k, i, j) = calc_rsif(pres(k, i, j), temp(k, i, j))
            else
              qvsi(k, i, j) = qvs(k, i, j)
            endif
            satw(k, i, j) = qv(k, i, j)/qvs(k, i, j)
            sati(k, i, j) = qv(k, i, j)/qvsi(k, i, j)
            ssatw(k, i, j) = satw(k, i, j) - 1._wp
            ssati(k, i, j) = sati(k, i, j) - 1._wp
            if (abs(ssatw(k, i, j)) < eps) ssatw(k, i, j) = 0._wp
            if (abs(ssati(k, i, j)) < eps) ssati(k, i, j) = 0._wp
            if (ssati(k, i, j) > 0._wp) column_supersaturated(i, j) = .true.
            diffu(k, i, j) = 2.11e-5_wp*(temp(k, i, j)/t0)**1.94_wp * (101325._wp/pres(k, i, j))
            if (tempc >= 0._wp) then
              visco(k, i, j) = (1.718_wp+0.0049_wp*tempc)*1.0e-5_wp
            else
              visco(k, i, j) = (1.718_wp+0.0049_wp*tempc-1.2e-5_wp*tempc*tempc)*1.0e-5_wp
            endif
            ocp(k, i, j) = 1._wp/(cp*(1._wp+0.887_wp*qv(k, i, j)))
            vsc2(k, i, j) = sqrt(rho(k, i, j)/visco(k, i, j))
            lvap(k, i, j) = lvap0 + (2106.0_wp - 4218.0_wp)*tempc
            tcond(k, i, j) = (5.69_wp + 0.0168_wp*tempc)*1.0e-5_wp * 418.936_wp
            lvt2(k, i, j) = lvap(k, i, j)*lvap(k, i, j)*ocp(k, i, j)*orv*otemp*otemp
          enddo
        endif
      enddo
    enddo
    !$acc end parallel
  end subroutine thermo_vars


  subroutine check_over_depletion(kts, kte, its, ite, jts, jte, rho, temp, qvsi, qv, l_qc, rc, l_qi, ri, &
    l_qr, rr, l_qs, rs, l_qg, rg, tend, odt, column_mp_active)
    !! check to ensure that loss terms don't over-deplete a category (horizontal tile)
    use module_mp_tempo_params, only : eps, rho_i, t0, meters3_to_liters

    integer, intent(in) :: kts, kte, its, ite, jts, jte
    real(wp), intent(in) :: odt
    type(ty_tend), intent(inout) :: tend
    logical, dimension(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(in) :: l_qc, l_qi, l_qr, l_qs, l_qg
    real(wp), dimension(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(in) :: rho, temp, qvsi, qv, rc, ri, rr, rs, rg
    logical, dimension(TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(in), optional :: column_mp_active
    real(wp) :: sump, rate_max, ratio
    integer :: i, j, k
    logical :: active_col

    !$acc parallel
    !$acc loop gang collapse(2) private(active_col)
    do j = TEMPO_JTS, TEMPO_JTE
      do i = TEMPO_ITS, TEMPO_ITE
        active_col = .true.
        if (present(column_mp_active)) active_col = column_mp_active(i, j)
        if (active_col) then

        !$acc loop vector private(sump, rate_max, ratio)
        do k = kts, kte
          sump = tend%pri_inu(k,i,j) + tend%pri_ide(k,i,j) + tend%prs_ide(k,i,j) + &
            tend%prs_sde(k,i,j) + tend%prg_gde(k,i,j) + tend%pri_iha(k,i,j)
          rate_max = (qv(k,i,j)-qvsi(k,i,j))*rho(k,i,j)*odt*0.999_wp
          if ((sump > eps .and. sump > rate_max) .or. &
            (sump < -eps .and. sump < rate_max)) then
            ratio = rate_max/sump
            tend%pri_inu(k,i,j) = tend%pri_inu(k,i,j) * ratio
            tend%pri_ide(k,i,j) = tend%pri_ide(k,i,j) * ratio
            tend%pni_ide(k,i,j) = tend%pni_ide(k,i,j) * ratio
            tend%prs_ide(k,i,j) = tend%prs_ide(k,i,j) * ratio
            tend%prs_sde(k,i,j) = tend%prs_sde(k,i,j) * ratio
            tend%prg_gde(k,i,j) = tend%prg_gde(k,i,j) * ratio
            tend%pri_iha(k,i,j) = tend%pri_iha(k,i,j) * ratio
          endif

          sump = -tend%prr_wau(k,i,j) - tend%pri_wfz(k,i,j) - tend%prr_rcw(k,i,j) - &
            tend%prs_scw(k,i,j) - tend%prg_scw(k,i,j) - tend%prg_gcw(k,i,j)
          rate_max = -rc(k,i,j)*odt
          if (l_qc(k,i,j)) then
            if (sump < rate_max) then
              ratio = rate_max/sump
              tend%prr_wau(k,i,j) = tend%prr_wau(k,i,j) * ratio
              tend%pri_wfz(k,i,j) = tend%pri_wfz(k,i,j) * ratio
              tend%prr_rcw(k,i,j) = tend%prr_rcw(k,i,j) * ratio
              tend%prs_scw(k,i,j) = tend%prs_scw(k,i,j) * ratio
              tend%prg_scw(k,i,j) = tend%prg_scw(k,i,j) * ratio
              tend%prg_gcw(k,i,j) = tend%prg_gcw(k,i,j) * ratio
            endif
          endif

          sump = tend%pri_ide(k,i,j) - tend%prs_iau(k,i,j) - tend%prs_sci(k,i,j) - tend%pri_rci(k,i,j)
          rate_max = -ri(k,i,j)*odt
          if (l_qi(k,i,j)) then
            if (sump < rate_max) then
              ratio = rate_max/sump
              tend%pri_ide(k,i,j) = tend%pri_ide(k,i,j) * ratio
              tend%prs_iau(k,i,j) = tend%prs_iau(k,i,j) * ratio
              tend%prs_sci(k,i,j) = tend%prs_sci(k,i,j) * ratio
              tend%pri_rci(k,i,j) = tend%pri_rci(k,i,j) * ratio
            endif
          endif

          sump = -tend%prg_rfz(k,i,j) - tend%pri_rfz(k,i,j) - tend%prr_rci(k,i,j) + &
            tend%prr_rcs(k,i,j) + tend%prr_rcg(k,i,j)
          rate_max = -rr(k,i,j)*odt
          if (l_qr(k,i,j)) then
            if (sump < rate_max) then
              ratio = rate_max/sump
              tend%prg_rfz(k,i,j) = tend%prg_rfz(k,i,j) * ratio
              tend%pbg_rfz(k,i,j) = tend%pbg_rfz(k,i,j) * ratio
              tend%pri_rfz(k,i,j) = tend%pri_rfz(k,i,j) * ratio
              tend%prr_rci(k,i,j) = tend%prr_rci(k,i,j) * ratio
              tend%prr_rcs(k,i,j) = tend%prr_rcs(k,i,j) * ratio
              tend%prr_rcg(k,i,j) = tend%prr_rcg(k,i,j) * ratio
            endif
          endif

          sump = tend%prs_sde(k,i,j) - tend%prs_ihm(k,i,j) - tend%prr_sml(k,i,j) + &
            tend%prs_rcs(k,i,j)
          rate_max = -rs(k,i,j)*odt
          if (l_qs(k,i,j)) then
            if (sump < rate_max) then
              ratio = rate_max/sump
              tend%prs_sde(k,i,j) = tend%prs_sde(k,i,j) * ratio
              tend%prs_ihm(k,i,j) = tend%prs_ihm(k,i,j) * ratio
              tend%prr_sml(k,i,j) = tend%prr_sml(k,i,j) * ratio
              tend%prs_rcs(k,i,j) = tend%prs_rcs(k,i,j) * ratio
            endif
          endif

          sump = tend%prg_gde(k,i,j) - tend%prg_ihm(k,i,j) - tend%prr_gml(k,i,j) + tend%prg_rcg(k,i,j)
          rate_max = -rg(k,i,j)*odt
          if (l_qg(k,i,j)) then
            if (sump < rate_max) then
              ratio = rate_max/sump
              tend%prg_gde(k,i,j) = tend%prg_gde(k,i,j) * ratio
              tend%prg_ihm(k,i,j) = tend%prg_ihm(k,i,j) * ratio
              tend%prr_gml(k,i,j) = tend%prr_gml(k,i,j) * ratio
              tend%prg_rcg(k,i,j) = tend%prg_rcg(k,i,j) * ratio
              tend%pbg_rcg(k,i,j) = tend%pbg_rcg(k,i,j) * ratio
            endif
          endif

          tend%pri_ihm(k,i,j) = tend%prs_ihm(k,i,j) + tend%prg_ihm(k,i,j)
          ratio = min(abs(tend%prr_rcg(k,i,j)), abs(tend%prg_rcg(k,i,j)))
          tend%prr_rcg(k,i,j) = ratio * sign(1.0_dp, tend%prr_rcg(k,i,j))
          tend%prg_rcg(k,i,j) = -tend%prr_rcg(k,i,j)
          tend%pbg_rcg(k,i,j) = meters3_to_liters*tend%prg_rcg(k,i,j)/rho_i
          if (temp(k,i,j) > t0) then
            ratio = min(abs(tend%prr_rcs(k,i,j)), abs(tend%prs_rcs(k,i,j)))
            tend%prr_rcs(k,i,j) = ratio * sign(1.0_dp, tend%prr_rcs(k,i,j))
            tend%prs_rcs(k,i,j) = -tend%prr_rcs(k,i,j)
          endif
        enddo
        endif
      enddo
    enddo
    !$acc end parallel
  end subroutine check_over_depletion


  subroutine sum_tendencies(kts, kte, its, ite, jts, jte, rho, temp, idx, lvap, ocp, tend, tten, qvten, qcten, &
    ncten, qiten, niten, qsten, qrten, nrten, qgten, ngten, qbten, column_mp_active)
    !! sums tendencies for each hydrometeor category and temperature and moisture (horizontal tile)
    use module_mp_tempo_params, only : lsub, rho_g, t0, lfus, meters3_to_liters

    integer, intent(in) :: kts, kte, its, ite, jts, jte
    type(ty_tend), intent(in) :: tend
    real(wp), dimension(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(in) :: rho, temp, lvap, ocp
    integer, dimension(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(in) :: idx
    real(wp), dimension(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(inout) :: qvten, qcten, ncten, qiten, niten, &
      qsten, qrten, nrten, qgten, ngten, qbten, tten
    logical, dimension(TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(in), optional :: column_mp_active
    real(wp) :: orho, lfus2
    integer :: i, j, k
    logical :: active_col

    !$acc parallel
    !$acc loop gang collapse(2) private(active_col)
    do j = TEMPO_JTS, TEMPO_JTE
      do i = TEMPO_ITS, TEMPO_ITE
        active_col = .true.
        if (present(column_mp_active)) active_col = column_mp_active(i, j)
        if (active_col) then

        !$acc loop vector private(orho, lfus2)
        do k = kts, kte
          orho = 1./rho(k,i,j)
          lfus2 = lsub - lvap(k,i,j)

          qvten(k,i,j) = qvten(k,i,j) + (-tend%pri_inu(k,i,j) - tend%pri_iha(k,i,j) - tend%pri_ide(k,i,j) - &
            tend%prs_ide(k,i,j) - tend%prs_sde(k,i,j) - tend%prg_gde(k,i,j)) * orho

          qcten(k,i,j) = qcten(k,i,j) + (-tend%prr_wau(k,i,j) - tend%pri_wfz(k,i,j) - tend%prr_rcw(k,i,j) - &
            tend%prs_scw(k,i,j) - tend%prg_scw(k,i,j) - tend%prg_gcw(k,i,j)) * orho

          ncten(k,i,j) = ncten(k,i,j) + (-tend%pnc_wau(k,i,j) - tend%pnc_rcw(k,i,j) - tend%pni_wfz(k,i,j) - &
            tend%pnc_scw(k,i,j) - tend%pnc_gcw(k,i,j)) * orho

          qiten(k,i,j) = qiten(k,i,j) + (tend%pri_inu(k,i,j) + tend%pri_iha(k,i,j) + tend%pri_ihm(k,i,j) + &
            tend%pri_wfz(k,i,j) + tend%pri_rfz(k,i,j) + tend%pri_ide(k,i,j) - tend%prs_iau(k,i,j) - &
            tend%prs_sci(k,i,j) - tend%pri_rci(k,i,j)) * orho

          niten(k,i,j) = niten(k,i,j) + (tend%pni_inu(k,i,j) + tend%pni_iha(k,i,j) + tend%pni_ihm(k,i,j) + &
            tend%pni_wfz(k,i,j) + tend%pni_rfz(k,i,j) + tend%pni_ide(k,i,j) - tend%pni_iau(k,i,j) - &
            tend%pni_sci(k,i,j) - tend%pni_rci(k,i,j)) * orho

          qrten(k,i,j) = qrten(k,i,j) + (tend%prr_wau(k,i,j) + tend%prr_rcw(k,i,j) + tend%prr_sml(k,i,j) + &
            tend%prr_gml(k,i,j) + tend%prr_rcs(k,i,j) + tend%prr_rcg(k,i,j) - tend%prg_rfz(k,i,j) - &
            tend%pri_rfz(k,i,j) - tend%prr_rci(k,i,j)) * orho

          nrten(k,i,j) = nrten(k,i,j) + (tend%pnr_wau(k,i,j) + tend%pnr_sml(k,i,j) + tend%pnr_gml(k,i,j) - &
            (tend%pnr_rfz(k,i,j) + tend%pnr_rcr(k,i,j) + tend%pnr_rcg(k,i,j) + tend%pnr_rcs(k,i,j) + &
            tend%pnr_rci(k,i,j) + tend%pni_rfz(k,i,j))) * orho

          qsten(k,i,j) = qsten(k,i,j) + (tend%prs_iau(k,i,j) + tend%prs_sde(k,i,j) + tend%prs_sci(k,i,j) + &
            tend%prs_scw(k,i,j) + tend%prs_rcs(k,i,j) + tend%prs_ide(k,i,j) - tend%prs_ihm(k,i,j) - &
            tend%prr_sml(k,i,j)) * orho

          qgten(k,i,j) = qgten(k,i,j) + (tend%prg_scw(k,i,j) + tend%prg_rfz(k,i,j) + tend%prg_gde(k,i,j) + &
            tend%prg_rcg(k,i,j) + tend%prg_gcw(k,i,j) + tend%prg_rci(k,i,j) + tend%prg_rcs(k,i,j) - &
            tend%prg_ihm(k,i,j) - tend%prr_gml(k,i,j)) * orho

          ngten(k,i,j) = ngten(k,i,j) + (tend%png_scw(k,i,j) + tend%png_rfz(k,i,j) - tend%png_rcg(k,i,j) + &
            tend%png_rci(k,i,j) + tend%png_rcs(k,i,j) + tend%png_gde(k,i,j) - tend%pnr_gml(k,i,j)) * orho

          qbten(k,i,j) = qbten(k,i,j) + (tend%pbg_scw(k,i,j) + tend%pbg_rfz(k,i,j) + tend%pbg_gcw(k,i,j) + &
            tend%pbg_rci(k,i,j) + tend%pbg_rcs(k,i,j) + tend%pbg_rcg(k,i,j) + tend%pbg_sml(k,i,j) - &
            tend%pbg_gml(k,i,j) + meters3_to_liters * (tend%prg_gde(k,i,j) - tend%prg_ihm(k,i,j)) / rho_g(idx(k,i,j))) * orho

          if (temp(k,i,j) < t0) then
            tten(k,i,j) = tten(k,i,j) + &
              (lsub*ocp(k,i,j)*(tend%pri_inu(k,i,j) + tend%pri_ide(k,i,j) + &
              tend%prs_ide(k,i,j) + tend%prs_sde(k,i,j) + tend%prg_gde(k,i,j) + tend%pri_iha(k,i,j)) + &
              lfus2*ocp(k,i,j)*(tend%pri_wfz(k,i,j) + tend%pri_rfz(k,i,j) + tend%prg_rfz(k,i,j) + &
              tend%prs_scw(k,i,j) + tend%prg_scw(k,i,j) + tend%prg_gcw(k,i,j) + tend%prg_rcs(k,i,j) + &
              tend%prs_rcs(k,i,j) + tend%prr_rci(k,i,j) + tend%prg_rcg(k,i,j)))*orho
          else
            tten(k,i,j) = tten(k,i,j) + &
              (lfus*ocp(k,i,j)*(-tend%prr_sml(k,i,j) - tend%prr_gml(k,i,j) - &
              tend%prr_rcg(k,i,j) - tend%prr_rcs(k,i,j)) + &
              lsub*ocp(k,i,j)*(tend%prs_sde(k,i,j) + tend%prg_gde(k,i,j)))*orho
          endif
        enddo
        endif
      enddo
    enddo
    !$acc end parallel
  end subroutine sum_tendencies


  subroutine sedimentation(kts, kte, its, ite, jts, jte, dz3d, rho, xr, xten, vt, substeps_sedi, ktop_sedi, limit, dt, &
      column_mp_active, col_any, n, precip)
    !! Eulerian sedimentation over the horizontal tile (OpenACC: parallel over columns).
    !! When n is present, only columns with n <= substeps_sedi(i,j) are updated (substep loop).

    integer, intent(in) :: kts, kte, its, ite, jts, jte
    real(wp), intent(in) :: dt, limit
    integer, intent(in), optional :: n
    real(wp), dimension(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(in) :: dz3d, rho
    real(wp), dimension(kts:kte+1, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(in) :: vt
    real(wp), dimension(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(inout) :: xr, xten
    integer, dimension(TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(in) :: substeps_sedi, ktop_sedi
    logical, dimension(TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(in) :: column_mp_active, col_any
    real(wp), dimension(:,:), intent(inout), optional :: precip !! assumed-shape so a full-tile diag can be accumulated per sub-tile block (absolute i,j)
    integer :: i, j, k, ktop, steps
    real(wp) :: odz, orho
    real(wp) :: sed_r(kts:kte)
    logical :: use_n

    use_n = present(n)

    !$acc parallel
    !$acc loop gang collapse(2) private(sed_r, ktop, steps, odz, orho, k)
    do j = TEMPO_JTS, TEMPO_JTE
      do i = TEMPO_ITS, TEMPO_ITE
        if (column_mp_active(i,j) .and. col_any(i,j)) then
          if (.not. use_n .or. n <= substeps_sedi(i,j)) then

            steps = substeps_sedi(i,j)
            ktop = ktop_sedi(i,j)

            !$acc loop seq
            do k = kte, kts, -1
              sed_r(k) = vt(k,i,j)*xr(k,i,j)
            enddo
            k = kte
            odz = 1._wp/dz3d(k,i,j)
            orho = 1._wp/rho(k,i,j)
            xten(k,i,j) = xten(k,i,j) - sed_r(k)*odz*(1._wp/real(steps, kind=wp))*orho
            xr(k,i,j) = max(limit, xr(k,i,j) - sed_r(k)*odz*dt*(1._wp/real(steps, kind=wp)))

            !$acc loop seq
            do k = ktop, kts, -1
              odz = 1._wp/dz3d(k,i,j)
              orho = 1._wp/rho(k,i,j)
              xten(k,i,j) = xten(k,i,j) + (sed_r(k+1))*(1._wp/dz3d(k,i,j))*(1._wp/real(steps, kind=wp))/rho(k,i,j)
              xten(k,i,j) = xten(k,i,j) - (sed_r(k))*odz*(1._wp/real(steps, kind=wp))*orho

              xr(k,i,j) = max(limit, xr(k,i,j) + (sed_r(k+1))*(1._wp/dz3d(k,i,j))*dt*(1._wp/real(steps, kind=wp)))
              xr(k,i,j) = max(limit, xr(k,i,j) - (sed_r(k))*odz*dt*(1._wp/real(steps, kind=wp)))
            enddo

            if (present(precip)) then
              if (xr(kts,i,j) > low_limit_mass_for_precip) then
                precip(i,j) = precip(i,j) + sed_r(kts)*dt*(1._wp/real(steps, kind=wp))
              endif
            endif

          endif
        endif
      enddo
    enddo
    !$acc end parallel
  end subroutine sedimentation


  subroutine semilagrangian_sedimentation(kts, kte, its, ite, jts, jte, dz3d, rho, xr, xten, vt, substeps_sedi, limit, dt, odt, &
      column_mp_active, col_any, n, precip)
    !! Semi-Lagrangian sedimentation over the horizontal tile (OpenACC: parallel over columns).
    !! [Juang and Hong (2010)](https://doi.org/10.1175/2009MWR3109.1)

    integer, intent(in) :: kts, kte, its, ite, jts, jte, n
    real(wp), intent(in) :: dt, odt, limit
    real(wp), dimension(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(in) :: dz3d, rho
    real(wp), dimension(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(inout) :: xr, xten
    real(wp), dimension(kts:kte+1, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(in) :: vt
    integer, dimension(TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(in) :: substeps_sedi
    logical, dimension(TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(in) :: column_mp_active, col_any
    real(wp), dimension(:,:), intent(inout), optional :: precip !! assumed-shape so a full-tile diag can be accumulated per sub-tile block (absolute i,j)
    integer :: i, j, k, kk, kb, kt, m, steps
    real(wp) :: fa1, fa2, con1, decfl, dip, slope_dim, tl, th, tl2, th2, qqd, qqh, qql, zsum, qsum, orho, dql, dqh
    real(wp) :: zi(kts:kte+1), wi(kts:kte+1), za(kts:kte+2), dza(kts:kte+1), qa(kts:kte+1), qmi(kts:kte+1), qpi(kts:kte+1), &
      net_flx(kts:kte), precip_flx(kts:kte), rr_save(kts:kte)

    !$acc parallel
    !$acc loop gang collapse(2) private(steps, fa1, fa2, con1, decfl, dip, slope_dim, tl, th, tl2, th2, qqd, qqh, qql, zsum, qsum, &
    !$acc   orho, dql, dqh, zi, wi, za, dza, qa, qmi, qpi, net_flx, precip_flx, rr_save, k, kk, kb, kt, m)
    do j = TEMPO_JTS, TEMPO_JTE
      do i = TEMPO_ITS, TEMPO_ITE
        if (column_mp_active(i,j) .and. col_any(i,j) .and. n <= substeps_sedi(i,j)) then
          steps = substeps_sedi(i,j)

          rr_save = xr(kts:kte,i,j)
          zi = 0._wp
          wi = 0._wp
          za = 0._wp
          dza = 0._wp
          qa = 0._wp
          qmi = 0._wp
          qpi = 0._wp
          net_flx = 0._wp
          precip_flx = 0._wp
          zi(kts) = 0._wp
          do k = kts, kte
            zi(k+1) = zi(k) + dz3d(k,i,j)
          enddo

          fa1 = 9._wp/16._wp
          fa2 = 1._wp/16._wp
          wi(kts) = vt(kts,i,j)
          wi(kts+1) = 0.5_wp*(vt(kts+1,i,j)+vt(kts,i,j))
          do k = kts+2, kte-1
            wi(k) = fa1*(vt(k,i,j)+vt(k-1,i,j))-fa2*(vt(k+1,i,j)+vt(k-2,i,j))
          enddo
          wi(kte) = 0.5_wp*(vt(kte,i,j)+vt(kte-1,i,j))
          wi(kte+1) = vt(kte+1,i,j)

          do k = kts+1, kte
            if(vt(k,i,j) == 0._wp) wi(k) = vt(k-1,i,j)
          enddo

          con1 = 0.05_wp
          do k = kte, kts, -1
            decfl = (wi(k+1)-wi(k))*dt*(1._wp/real(steps, kind=wp))/dz3d(k,i,j)
            if(decfl > con1) then
              wi(k) = wi(k+1) - con1*dz3d(k,i,j)*odt*real(steps, kind=wp)
            endif
          enddo

          do k = kts, kte+1
            za(k) = zi(k) - wi(k)*dt*(1._wp/real(steps, kind=wp))
          enddo
          za(kte+2) = zi(kte+1)
          do k = kts, kte+1
            dza(k) = za(k+1)-za(k)
          enddo

          do k = kts, kte
            qa(k) = xr(k,i,j)*dz3d(k,i,j)/dza(k)
          enddo
          qa(kte+1) = 0._wp

          do k = kts+1, kte
            dip=(qa(k+1)-qa(k))/(dza(k+1)+dza(k))
            slope_dim=(qa(k)-qa(k-1))/(dza(k-1)+dza(k))
            if(dip*slope_dim <= 0._wp) then
              qmi(k)=qa(k)
              qpi(k)=qa(k)
            else
              qpi(k)=qa(k)+0.5_wp*(dip+slope_dim)*dza(k)
              qmi(k)=2._wp*qa(k)-qpi(k)
              if(qpi(k) > 0._wp .or. qmi(k) < 0._wp) then
                qpi(k) = qa(k)
                qmi(k) = qa(k)
              endif
            endif
          enddo
          qpi(kts) = qa(kts)
          qmi(kts) = qa(kts)
          qmi(kte+1) = qa(kte+1)
          qpi(kte+1) = qa(kte+1)

          kb = kts
          kt = kts
          intp : do k = kts, kte
            kb = max(kb-1,kts)
            kt = max(kt-1,kts)
            if(zi(k) >= za(kte+1)) then
              exit intp
            else
              find_kb : do kk = kb, kte
                if(zi(k) <= za(kk+1)) then
                  kb = kk
                  exit find_kb
                endif
              enddo find_kb
              find_kt : do kk = kt, kte+2
                if(zi(k+1) <= za(kk)) then
                  kt = kk
                  exit find_kt
                endif
              enddo find_kt
              kt = kt - 1

              if(kt == kb) then
                tl = (zi(k)-za(kb))/dza(kb)
                th = (zi(k+1)-za(kb))/dza(kb)
                tl2 = tl*tl
                th2 = th*th
                qqd = 0.5_wp*(qpi(kb)-qmi(kb))
                qqh = qqd*th2+qmi(kb)*th
                qql = qqd*tl2+qmi(kb)*tl
                xr(k,i,j) = (qqh-qql)/(th-tl)
              elseif(kt > kb) then
                tl = (zi(k)-za(kb))/dza(kb)
                tl2 = tl*tl
                qqd = 0.5_wp*(qpi(kb)-qmi(kb))
                qql = qqd*tl2+qmi(kb)*tl
                dql = qa(kb)-qql
                zsum = (1._wp-tl)*dza(kb)
                qsum = dql*dza(kb)
                if(kt-kb > 1) then
                  do m = kb+1, kt-1
                    zsum = zsum + dza(m)
                    qsum = qsum + qa(m) * dza(m)
                  enddo
                endif
                th = (zi(k+1)-za(kt))/dza(kt)
                th2 = th*th
                qqd = 0.5_wp*(qpi(kt)-qmi(kt))
                dqh = qqd*th2+qmi(kt)*th
                zsum = zsum + th*dza(kt)
                qsum = qsum + dqh*dza(kt)
                xr(k,i,j) = qsum/zsum
              endif
            endif
            orho = 1._wp / rho(k,i,j)
            xr(k,i,j) = max(xr(k,i,j), limit)
            xten(k,i,j) = xten(k,i,j) + (xr(k,i,j) - rr_save(k)) * orho*odt
          enddo intp

          precip_loop: do k = kts, kte
            if(za(k) < 0._wp .and. za(k+1) <= 0.0_wp) then
              if (present(precip)) precip(i,j) = precip(i,j) + qa(k)*dza(k)
              net_flx(k) = qa(k)*dza(k)
            elseif (za(k) < 0._wp .and. za(k+1) > 0._wp) then
              th = (0._wp-za(k))/dza(k)
              th2 = th*th
              qqd = 0.5_wp*(qpi(k)-qmi(k))
              qqh = qqd*th2+qmi(k)*th
              if (present(precip)) precip(i,j) = precip(i,j) + qqh*dza(k)
              net_flx(k) = qqh*dza(k)
              exit precip_loop
            endif
          enddo precip_loop

          do k = kte, kts, -1
            if(k == kte) then
              precip_flx(k) = net_flx(k)
            else
              precip_flx(k) = precip_flx(k+1) + net_flx(k)
            endif
          enddo
        endif
      enddo
    enddo
    !$acc end parallel
  end subroutine semilagrangian_sedimentation


  subroutine rain_fallspeed(kts, kte, its, ite, jts, jte, dt, rhof, l_qr, rr, ilamr, dz3d, vt, vtn, substeps_sedi, ktop_sedi, &
      column_mp_active, l_gr_flag, substep_mode, col_any, n)
    !! Mass- and number-weighted rain fall speeds over k and the horizontal tile [TEMPO_ITS:TEMPO_ITE] x [TEMPO_JTS:TEMPO_JTE] (cf. graupel_fallspeed).
    !! Full mode sets substeps_sedi/ktop_sedi; substep mode updates only the shared backward-k fall-speed loop.
    !! OpenACC: parallel over columns; vertical k is sequential (backward recurrence on vt/vtn).
    use module_mp_tempo_params, only : crg, av_r, org3, fv_r, cre, r1

    integer, intent(in) :: kts, kte, its, ite, jts, jte
    real(wp), intent(in) :: dt
    real(wp), dimension(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(in) :: rhof, rr
    logical, dimension(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(in) :: l_qr
    real(dp), dimension(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(in) :: ilamr
    real(wp), dimension(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(in) :: dz3d
    real(wp), dimension(kts:kte+1, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(inout) :: vt, vtn
    integer, dimension(TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(inout) :: substeps_sedi, ktop_sedi
    logical, dimension(TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(in), optional :: column_mp_active, col_any
    logical, intent(in), optional :: l_gr_flag, substep_mode !! l_gr_flag: skip column if no rain (any(l_qr))
    integer, intent(in), optional :: n
    integer :: i, j, k
    real(wp) :: dz_by_vt
    real(dp) :: lamr
    logical :: use_cmp, gr_skip_if_empty, col_has_qr, use_substep

    use_cmp = present(column_mp_active)
    gr_skip_if_empty = .false.
    if (present(l_gr_flag)) gr_skip_if_empty = l_gr_flag
    use_substep = .false.
    if (present(substep_mode)) use_substep = substep_mode

    if (use_substep) then
      if (.not. present(column_mp_active) .or. .not. present(col_any) .or. .not. present(n)) then
        error stop "rain_fallspeed: substep mode requires column_mp_active, col_any, and n"
      endif
    endif

    !$acc parallel
    !$acc loop gang vector collapse(2) private(lamr, dz_by_vt, col_has_qr)
    do j = TEMPO_JTS, TEMPO_JTE
      do i = TEMPO_ITS, TEMPO_ITE
        if (use_cmp) then
          if (.not. column_mp_active(i, j)) cycle
        endif
        if (use_substep) then
          if (.not. col_any(i, j)) cycle
          if (n > substeps_sedi(i, j)) cycle
        elseif (gr_skip_if_empty) then
          col_has_qr = .false.
          !$acc loop seq
          do k = kts, kte
            if (l_qr(k, i, j)) then
              col_has_qr = .true.
              exit
            endif
          enddo
          if (.not. col_has_qr) cycle
        endif

        if (.not. use_substep) then
          ktop_sedi(i, j) = 1
          substeps_sedi(i, j) = 1
        endif
        !$acc loop seq
        do k = kte, kts, -1
          if (rr(k, i, j) > r1) then
            lamr = 1._dp / ilamr(k, i, j)
            vt(k, i, j) = rhof(k, i, j)*av_r*crg(6)*org3 * lamr**cre(3) *((lamr+fv_r)**(-cre(6)))
            vtn(k, i, j) = rhof(k, i, j)*av_r*crg(7)/crg(12) * lamr**cre(12)*((lamr+fv_r)**(-cre(7)))
          else
            vt(k, i, j) = vt(k+1, i, j)
            vtn(k, i, j) = vtn(k+1, i, j)
          endif
          if (.not. use_substep) then
            if (max(vt(k, i, j), vtn(k, i, j)) > 1.e-3_wp) then
              ktop_sedi(i, j) = max(ktop_sedi(i, j), k)
              dz_by_vt = dz3d(k, i, j) / (max(vt(k, i, j), vtn(k, i, j)))
              substeps_sedi(i, j) = max(substeps_sedi(i, j), int(dt/dz_by_vt + 1._wp))
            endif
          endif
        enddo
        if (.not. use_substep) then
          if (ktop_sedi(i, j) == kte) ktop_sedi(i, j) = kte - 1
        endif
      enddo
    enddo
    !$acc end parallel
  end subroutine rain_fallspeed


  subroutine graupel_fallspeed(kts, kte, its, ite, jts, jte, dt, rhof, rho, visco, l_qg, rg, rb, idx_bg, ilamg, dz3d, vt, vtn, &
      substeps_sedi, ktop_sedi, column_mp_active, l_gg_flag, qb3d, substep_mode, col_any, n)
    !! Mass- and number-weighted graupel fall speeds over k and the horizontal tile [TEMPO_ITS:TEMPO_ITE] x [TEMPO_JTS:TEMPO_JTE] (cf. cloud_check_and_update).
    !! Full mode sets substeps_sedi/ktop_sedi; substep mode updates only the shared backward-k fall-speed loop.
    !! Optional qb3d selects hail-aware fall-speed coefficients.
    !! OpenACC: parallel over columns; vertical k is sequential (backward recurrence on vt/vtn).
    use module_mp_tempo_params, only : nrhg, rho_g, av_g_old, bv_g_old, &
      cgg, t0, mu_g, ogg2, ogg3, a_coeff, b_coeff, meters3_to_liters, r1

    integer, intent(in) :: kts, kte, its, ite, jts, jte
    real(wp), intent(in) :: dt
    real(wp), dimension(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(in) :: rhof, rho, visco, rg, rb
    logical, dimension(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(in) :: l_qg
    integer, dimension(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(in) :: idx_bg
    real(dp), dimension(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(in) :: ilamg
    real(wp), dimension(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(in) :: dz3d
    real(wp), dimension(kts:kte+1, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(inout) :: vt, vtn
    integer, dimension(TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(inout) :: substeps_sedi, ktop_sedi
    logical, dimension(TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(in), optional :: column_mp_active, col_any
    logical, intent(in), optional :: l_gg_flag, substep_mode !! l_gg_flag: skip column if no graupel (any(l_qg))
    real(wp), dimension(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(in), optional :: qb3d
    integer, intent(in), optional :: n
    integer :: i, j, k
    real(wp) :: dz_by_vt, dens_g, afall, bfall
    logical :: use_cmp, gg_skip_if_empty, col_has_qg, use_qb3d, use_substep

    use_cmp = present(column_mp_active)
    gg_skip_if_empty = .false.
    if (present(l_gg_flag)) gg_skip_if_empty = l_gg_flag
    use_qb3d = present(qb3d)
    use_substep = .false.
    if (present(substep_mode)) use_substep = substep_mode

    if (use_substep) then
      if (.not. present(column_mp_active) .or. .not. present(col_any) .or. .not. present(n)) then
        error stop "graupel_fallspeed: substep mode requires column_mp_active, col_any, and n"
      endif
    endif

    !$acc parallel
    !$acc loop gang vector collapse(2) private(dz_by_vt, dens_g, afall, bfall, col_has_qg)
    do j = TEMPO_JTS, TEMPO_JTE
      do i = TEMPO_ITS, TEMPO_ITE
        if (use_cmp) then
          if (.not. column_mp_active(i, j)) cycle
        endif
        if (use_substep) then
          if (.not. col_any(i, j)) cycle
          if (n > substeps_sedi(i, j)) cycle
        elseif (gg_skip_if_empty) then
          col_has_qg = .false.
          !$acc loop seq
          do k = kts, kte
            if (l_qg(k, i, j)) then
              col_has_qg = .true.
              exit
            endif
          enddo
          if (.not. col_has_qg) cycle
        endif

        if (.not. use_substep) then
          ktop_sedi(i, j) = 1
          substeps_sedi(i, j) = 1
        endif
        !$acc loop seq
        do k = kte, kts, -1
          if (rg(k, i, j) > r1) then
            if (use_qb3d) then
              dens_g = max(rho_g(1), min(meters3_to_liters*rg(k, i, j)/rb(k, i, j), rho_g(nrhg)))
              afall = a_coeff*((4._wp*dens_g*9.8_wp)/(3._wp*rho(k, i, j)))**b_coeff
              afall = afall * visco(k, i, j)**(1._wp-2._wp*b_coeff)
              bfall = 3._wp*b_coeff - 1._wp
            else
              afall = av_g_old
              bfall = bv_g_old
            endif
            vt(k, i, j) = rhof(k, i, j)*afall*cgg(6, idx_bg(k, i, j))*ogg3 * ilamg(k, i, j)**bfall

            if (mu_g == 0) then
              vtn(k, i, j) = rhof(k, i, j)*afall*cgg(7, idx_bg(k, i, j))/cgg(12, idx_bg(k, i, j)) * ilamg(k, i, j)**bfall
            else
              vtn(k, i, j) = rhof(k, i, j)*afall*cgg(8, idx_bg(k, i, j))*ogg2 * ilamg(k, i, j)**bfall
            endif
          else
            vt(k, i, j) = vt(k+1, i, j)
            vtn(k, i, j) = vtn(k+1, i, j)
          endif
          if (.not. use_substep) then
            if (vt(k, i, j) > 1.e-3_wp) then
              ktop_sedi(i, j) = max(ktop_sedi(i, j), k)
              dz_by_vt = dz3d(k, i, j) / vt(k, i, j)
              substeps_sedi(i, j) = max(substeps_sedi(i, j), int(dt/dz_by_vt + 1._wp))
            endif
          endif
        enddo
        if (.not. use_substep) then
          if (ktop_sedi(i, j) == kte) ktop_sedi(i, j) = kte - 1
        endif
      enddo
    enddo
    !$acc end parallel
  end subroutine graupel_fallspeed


  subroutine snow_fallspeed(kts, kte, its, ite, jts, jte, dt, rhof, l_qs, rs, prr_sml, smob, smoc, &
      rr, vtrr, dz3d, vt, vtboost, substeps_sedi, ktop_sedi, column_mp_active, l_qs_flag)
    !! Mass-weighted snow fall speeds over k and the horizontal tile [TEMPO_ITS:TEMPO_ITE] x [TEMPO_JTS:TEMPO_JTE] (cf. rain_fallspeed).
    !! Sets per-column substeps_sedi and ktop_sedi for subsequent sedimentation.
    !! OpenACC: parallel over columns; vertical k is sequential (backward recurrence on vt).
    use module_mp_tempo_params, only : lam0, lam1, fv_s, kap0, kap1, mu_s, &
      cse, csg, av_s, r1

    integer, intent(in) :: kts, kte, its, ite, jts, jte
    real(wp), intent(in) :: dt
    real(wp), dimension(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(in) :: rhof, dz3d, rs, rr, vtboost
    real(wp), dimension(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(in) :: vtrr
    logical, dimension(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(in) :: l_qs
    real(wp), dimension(kts:kte+1, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(inout) :: vt
    real(dp), dimension(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(in) :: smob, smoc
    real(dp), dimension(:,:,:), intent(in) :: prr_sml !! assumed-shape so a full-tile tendency can be read per sub-tile block (absolute i,j)
    integer, dimension(TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(inout) :: substeps_sedi, ktop_sedi
    logical, dimension(TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(in), optional :: column_mp_active
    logical, intent(in), optional :: l_qs_flag !! when .true., skip column if inactive or no snow (any(l_qs))
    integer :: i, j, k
    real(wp) :: dz_by_vt, vts, sr
    real(dp) :: xds, mrat, ils1, ils2, t1_vts, t2_vts, t3_vts, t4_vts
    logical :: use_cmp, qs_skip_if_empty, col_has_qs, compute_sedi

    use_cmp = present(column_mp_active)
    qs_skip_if_empty = .false.
    if (present(l_qs_flag)) qs_skip_if_empty = l_qs_flag

    !$acc parallel
    !$acc loop gang vector collapse(2) private(dz_by_vt, vts, sr, xds, mrat, ils1, ils2, t1_vts, t2_vts, t3_vts, t4_vts, col_has_qs, compute_sedi)
    do j = TEMPO_JTS, TEMPO_JTE
      do i = TEMPO_ITS, TEMPO_ITE
        compute_sedi = .true.
        if (use_cmp) then
          if (.not. column_mp_active(i, j)) then
            ktop_sedi(i, j) = 1
            substeps_sedi(i, j) = 1
            compute_sedi = .false.
          endif
        endif
        if (compute_sedi .and. qs_skip_if_empty) then
          col_has_qs = .false.
          !$acc loop seq
          do k = kts, kte
            if (l_qs(k, i, j)) then
              col_has_qs = .true.
              exit
            endif
          enddo
          if (.not. col_has_qs) then
            ktop_sedi(i, j) = 1
            substeps_sedi(i, j) = 1
            compute_sedi = .false.
          endif
        endif
        if (compute_sedi) then
          ktop_sedi(i, j) = 1
          substeps_sedi(i, j) = 1
          !$acc loop seq
          do k = kte, kts, -1
            if (rs(k, i, j) > r1) then
              xds = smoc(k, i, j) / smob(k, i, j)
              mrat = 1._dp/xds
              ils1 = 1._dp/(mrat*lam0 + fv_s)
              ils2 = 1._dp/(mrat*lam1 + fv_s)
              t1_vts = kap0*csg(4)*ils1**cse(4)
              t2_vts = kap1*mrat**mu_s*csg(10)*ils2**cse(10)
              ils1 = 1._dp/(mrat*lam0)
              ils2 = 1._dp/(mrat*lam1)
              t3_vts = kap0*csg(1)*ils1**cse(1)
              t4_vts = kap1*mrat**mu_s*csg(7)*ils2**cse(7)
              vts = rhof(k, i, j)*av_s * (t1_vts+t2_vts)/(t3_vts+t4_vts)

              if (prr_sml(k, i, j) > 0._dp) then
                sr = rs(k, i, j)/(rs(k, i, j)+rr(k, i, j))
                vt(k, i, j) = vts*sr + (1._wp-sr)*vtrr(k, i, j)
              else
                vt(k, i, j) = vts*vtboost(k, i, j)
              endif
            else
              vt(k, i, j) = vt(k+1, i, j)
            endif
            if (vt(k, i, j) > 1.e-3_wp) then
              ktop_sedi(i, j) = max(ktop_sedi(i, j), k)
              dz_by_vt = dz3d(k, i, j) / vt(k, i, j)
              substeps_sedi(i, j) = max(substeps_sedi(i, j), int(dt/dz_by_vt + 1._wp))
            endif
          enddo
          if (ktop_sedi(i, j) == kte) ktop_sedi(i, j) = kte - 1
        endif
      enddo
    enddo
    !$acc end parallel
  end subroutine snow_fallspeed


  subroutine ice_fallspeed(kts, kte, its, ite, jts, jte, dt, rhof, l_qi, ri, ilami, dz3d, vt, vtn, &
      substeps_sedi, ktop_sedi, column_mp_active, l_qi_flag)
    !! Mass and number weighted ice fall speeds over k and the horizontal tile (cf. snow_fallspeed).
    !! Sets per-column substeps_sedi and ktop_sedi for subsequent sedimentation.
    !! OpenACC: parallel over columns; vertical k is sequential (backward recurrence on vt/vtn).
    use module_mp_tempo_params, only : av_i, cig, oig2, bv_i, r1

    integer, intent(in) :: kts, kte, its, ite, jts, jte
    real(wp), intent(in) :: dt
    real(wp), dimension(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(in) :: rhof, dz3d, ri
    logical, dimension(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(in) :: l_qi
    real(dp), dimension(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(in) :: ilami
    real(wp), dimension(kts:kte+1, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(inout) :: vt, vtn
    integer, dimension(TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(inout) :: substeps_sedi, ktop_sedi
    logical, dimension(TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(in), optional :: column_mp_active
    logical, intent(in), optional :: l_qi_flag !! when .true., skip column if inactive or no ice (any(l_qi))
    integer :: i, j, k
    real(wp) :: dz_by_vt
    logical :: use_cmp, qi_skip_if_empty, col_has_qi, compute_sedi

    use_cmp = present(column_mp_active)
    qi_skip_if_empty = .false.
    if (present(l_qi_flag)) qi_skip_if_empty = l_qi_flag

    !$acc parallel
    !$acc loop gang vector collapse(2) private(dz_by_vt, col_has_qi, compute_sedi)
    do j = TEMPO_JTS, TEMPO_JTE
      do i = TEMPO_ITS, TEMPO_ITE
        compute_sedi = .true.
        if (use_cmp) then
          if (.not. column_mp_active(i, j)) then
            ktop_sedi(i, j) = 1
            substeps_sedi(i, j) = 1
            compute_sedi = .false.
          endif
        endif
        if (compute_sedi .and. qi_skip_if_empty) then
          col_has_qi = .false.
          !$acc loop seq
          do k = kts, kte
            if (l_qi(k, i, j)) then
              col_has_qi = .true.
              exit
            endif
          enddo
          if (.not. col_has_qi) then
            ktop_sedi(i, j) = 1
            substeps_sedi(i, j) = 1
            compute_sedi = .false.
          endif
        endif
        if (compute_sedi) then
          ktop_sedi(i, j) = 1
          substeps_sedi(i, j) = 1
          !$acc loop seq
          do k = kte, kts, -1
            if (ri(k, i, j) > r1) then
              vt(k, i, j) = rhof(k, i, j)*av_i*cig(3)*oig2 * ilami(k, i, j)**bv_i
              vtn(k, i, j) = rhof(k, i, j)*av_i*cig(6)/cig(7) * ilami(k, i, j)**bv_i
            else
              vt(k, i, j) = vt(k+1, i, j)
              vtn(k, i, j) = vtn(k+1, i, j)
            endif
            if (vt(k, i, j) > 1.e-3_wp) then
              ktop_sedi(i, j) = max(ktop_sedi(i, j), k)
              dz_by_vt = dz3d(k, i, j) / vt(k, i, j)
              substeps_sedi(i, j) = max(substeps_sedi(i, j), int(dt/dz_by_vt + 1._wp))
            endif
          enddo
          if (ktop_sedi(i, j) == kte) ktop_sedi(i, j) = kte - 1
        endif
      enddo
    enddo
    !$acc end parallel
  end subroutine ice_fallspeed


  subroutine cloud_fallspeed(kts, kte, its, ite, jts, jte, rhof, w3d, l_qc, rc, nc, ilamc, dz3d, vt, vtn, &
      ktop_sedi, column_mp_active, l_qc_flag)
    !! Mass and number weighted cloud fall speeds over k and the horizontal tile (cf. ice_fallspeed).
    !! Sets per-column ktop_sedi for subsequent sedimentation (substeps unchanged here).
    !! OpenACC: parallel over columns; vertical k loops are sequential (incl. 500 m agl cap without labeled exit).
    use module_mp_tempo_params, only : av_c, ccg, ocg1, ocg2, bv_c, r1, r2

    integer, intent(in) :: kts, kte, its, ite, jts, jte
    real(wp), dimension(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(in) :: rhof, w3d, dz3d, rc, nc
    real(dp), dimension(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(in) :: ilamc
    logical, dimension(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(in) :: l_qc
    real(wp), dimension(kts:kte+1, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(inout) :: vt, vtn
    integer, dimension(TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(inout) :: ktop_sedi
    logical, dimension(TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(in), optional :: column_mp_active
    logical, intent(in), optional :: l_qc_flag !! when .true., skip column if inactive or no cloud (any(l_qc))
    integer :: i, j, k, nu_c
    real(wp) :: hgt
    logical :: use_cmp, qc_skip_if_empty, col_has_qc, compute_sedi, still_within_500m_agl

    use_cmp = present(column_mp_active)
    qc_skip_if_empty = .false.
    if (present(l_qc_flag)) qc_skip_if_empty = l_qc_flag

    !$acc parallel
    !$acc loop gang vector collapse(2) private(hgt, nu_c, col_has_qc, compute_sedi, still_within_500m_agl)
    do j = TEMPO_JTS, TEMPO_JTE
      do i = TEMPO_ITS, TEMPO_ITE
        compute_sedi = .true.
        if (use_cmp) then
          if (.not. column_mp_active(i, j)) then
            ktop_sedi(i, j) = 1
            compute_sedi = .false.
          endif
        endif
        if (compute_sedi .and. qc_skip_if_empty) then
          col_has_qc = .false.
          !$acc loop seq
          do k = kts, kte
            if (l_qc(k, i, j)) then
              col_has_qc = .true.
              exit
            endif
          enddo
          if (.not. col_has_qc) then
            ktop_sedi(i, j) = 1
            compute_sedi = .false.
          endif
        endif
        if (compute_sedi) then
          ktop_sedi(i, j) = 1

          ! clouds/fog settle below 500 m agl
          hgt = 0._wp
          still_within_500m_agl = .true.
          !$acc loop seq
          do k = kts, kte - 1
            if (still_within_500m_agl) then
              if (rc(k, i, j) > r2) ktop_sedi(i, j) = k
              hgt = hgt + dz3d(k, i, j)
              if (hgt > 500._wp) still_within_500m_agl = .false.
            endif
          enddo

          !$acc loop seq
          do k = ktop_sedi(i, j), kts, -1
            if (rc(k, i, j) > r1 .and. w3d(k, i, j) < 0.1_wp) then
              nu_c = get_nuc(nc(k, i, j))
              vt(k, i, j) = rhof(k, i, j)*av_c*ccg(5, nu_c)*ocg2(nu_c) * ilamc(k, i, j)**bv_c
              vtn(k, i, j) = rhof(k, i, j)*av_c*ccg(4, nu_c)*ocg1(nu_c) * ilamc(k, i, j)**bv_c
            endif
          enddo
        endif
      enddo
    enddo
    !$acc end parallel
  end subroutine cloud_fallspeed


  subroutine cloud_condensation(kts, kte, its, ite, jts, jte, rho, temp, w1d, ssatw, lvap, tcond, diffu, lvt2, &
    nwfa, qv, qvs, l_qc, rc, nc, tend, dt, odt, column_mp_active)
    !! cloud condensation and evaporation (horizontal tile)
    !! OpenACC: parallel over columns; vertical k is sequential (early exit from k loop).
    use module_mp_tempo_params, only : eps, r1, t0, orv, pi, rho_w, nbc, &
      tnc_wev, nt_c_min

    integer, intent(in) :: kts, kte, its, ite, jts, jte
    real(wp), intent(in) :: dt, odt
    type(ty_tend), intent(inout) :: tend
    real(wp), dimension(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(in) :: rho, temp, w1d, ssatw, lvap, tcond, diffu, lvt2, &
      nwfa, qv, qvs, rc, nc
    logical, dimension(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(in) :: l_qc
    logical, dimension(TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(in), optional :: column_mp_active
    real(wp) :: clap, fcd, dfcd, xrc, xnc, orho, tempc, otemp, &
      rvs, rvs_p, rvs_pp, gamsc, alphsc, xsat, t1_evap
    real(dp) :: dc_star
    integer :: i, j, k, n, idx_d, idx_n, idx_c
    logical :: active_col

    !$acc parallel
    !$acc loop gang collapse(2) private(active_col)
    do j = TEMPO_JTS, TEMPO_JTE
      do i = TEMPO_ITS, TEMPO_ITE
        active_col = .true.
        if (present(column_mp_active)) active_col = column_mp_active(i, j)
        if (active_col) then

        !$acc loop seq private(clap, fcd, dfcd, xrc, xnc, orho, tempc, otemp, rvs, rvs_p, rvs_pp, gamsc, alphsc, &
        !$acc& xsat, t1_evap, dc_star, idx_d, idx_n, idx_c, n)
        do k = kts, kte
          if (abs(ssatw(k, i, j)) < eps) exit ! RH = 100%

          orho = 1._wp/rho(k, i, j)
          clap = (qv(k, i, j)-qvs(k, i, j))/(1._wp + lvt2(k, i, j)*qvs(k, i, j))
          do n = 1, 3
            fcd = qvs(k, i, j)*exp(lvt2(k, i, j)*clap) - qv(k, i, j) + clap
            dfcd = qvs(k, i, j)*lvt2(k, i, j)*exp(lvt2(k, i, j)*clap) + 1._wp
            clap = clap - fcd/dfcd
          enddo
          xrc = rc(k, i, j) + clap*rho(k, i, j)
          xnc = 0._wp

          if (xrc > r1) then
            ! mass tendency
            tend%prw_vcd(k,i,j) = clap*odt

            if (clap > eps) then ! condensation
              xnc = max(nt_c_min, activate_cloud_number(temp(k, i, j), w1d(k, i, j), nwfa(k, i, j)))
              tend%pnc_wcd(k,i,j) = 0.5_wp*(xnc-nc(k, i, j) + abs(xnc-nc(k, i, j)))*odt*orho
            elseif (l_qc(k, i, j) .and. ssatw(k, i, j) < -1.e-6_wp .and. clap < -eps) then ! evaporation
              tempc = temp(k, i, j) - t0
              otemp = 1._wp/temp(k, i, j)
              rvs = rho(k, i, j)*qvs(k, i, j)
              rvs_p = rvs*otemp*(lvap(k, i, j)*otemp*orv - 1._wp)
              rvs_pp = rvs * (otemp*(lvap(k, i, j)*otemp*oRv - 1._wp) * &
                otemp*(lvap(k, i, j)*otemp*oRv - 1._wp) + &
                (-2._wp*lvap(k, i, j)*otemp*otemp*otemp*oRv) + otemp*otemp)
              gamsc = lvap(k, i, j)*diffu(k, i, j)/tcond(k, i, j) * rvs_p
              alphsc = 0.5_wp*(gamsc/(1._wp+gamsc))*(gamsc/(1._wp+gamsc)) * &
                rvs_pp/rvs_p * rvs/rvs_p
              alphsc = max(1.e-9_wp, alphsc)
              xsat = ssatw(k, i, j)
              if (abs(xsat) < 1.e-9_wp) xsat = 0._wp
              t1_evap = 2._wp*pi*(1.0_wp - alphsc*xsat + 2._wp*alphsc*alphsc*xsat*xsat - &
                5._wp*alphsc*alphsc*alphsc*xsat*xsat*xsat) / (1._wp+gamsc)

              dc_star = sqrt(-2._dp*dt * t1_evap/(2._dp*pi) * &
                4._dp*diffu(k, i, j)*ssatw(k, i, j)*rvs/rho_w)
              idx_d = max(1, min(int(1.e6_wp*dc_star), nbc))
              call get_cloud_table_index(rc(k, i, j), nc(k, i, j), idx_c, idx_n)

              tend%prw_vcd(k,i,j) = max(real(-rc(k, i, j)*0.99_wp*orho*odt, kind=dp), &
                tend%prw_vcd(k,i,j))
              tend%pnc_wcd(k,i,j) = max(real(-nc(k, i, j)*0.99_wp*orho*odt, kind=dp), &
                -tnc_wev(idx_d, idx_c, idx_n)*orho*odt)
            endif
          else
            tend%prw_vcd(k,i,j) = -rc(k, i, j)*orho*odt
            tend%pnc_wcd(k,i,j) = -nc(k, i, j)*orho*odt
          endif
        enddo
        endif
      enddo
    enddo
    !$acc end parallel

  end subroutine cloud_condensation


  subroutine rain_evaporation(kts, kte, its, ite, jts, jte, rho, temp, ssatw, lvap, tcond, diffu, &
      vsc2, rhof2, qv, qvs, l_qr, rr, nr, ilamr, tend, odt, column_mp_active)
    !! rain evaporation that includes reduction in the evaporation rate
    !! in the presence of melting graupel (horizontal tile)
    use module_mp_tempo_params, only : eps, r1, t0, orv, pi, rho_w, &
      org2, cre, t1_qr_ev, t2_qr_ev, fv_r

    integer, intent(in) :: kts, kte, its, ite, jts, jte
    real(wp), intent(in) :: odt
    type(ty_tend), intent(inout) :: tend
    real(wp), dimension(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(in) :: rho, temp, ssatw, lvap, tcond, diffu, vsc2, rhof2, &
      qv, qvs, rr, nr
    logical, dimension(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(in) :: l_qr
    real(dp), dimension(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(in) :: ilamr
    logical, dimension(TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(in), optional :: column_mp_active
    integer :: i, j, k
    real(wp) :: orho, tempc, otemp, &
      rvs, rvs_p, rvs_pp, gamsc, alphsc, xsat, t1_evap, rate_max, eva_factor
    real(dp) :: lamr, n0_r
    logical :: active_col

    !$acc parallel
    !$acc loop gang collapse(2) private(active_col)
    do j = TEMPO_JTS, TEMPO_JTE
      do i = TEMPO_ITS, TEMPO_ITE
        active_col = .true.
        if (present(column_mp_active)) active_col = column_mp_active(i, j)
        if (active_col) then

        !$acc loop vector private(orho, tempc, otemp, rvs, rvs_p, rvs_pp, gamsc, alphsc, xsat, &
        !$acc& t1_evap, rate_max, eva_factor, lamr, n0_r)
        do k = kts, kte
          if (l_qr(k, i, j)) then
            if ((ssatw(k, i, j) < -eps) .and. tend%prw_vcd(k, i, j) <= 0._dp) then
              orho = 1._wp/rho(k, i, j)
              tempc = temp(k, i, j) - t0
              otemp = 1._wp/temp(k, i, j)
              rvs = rho(k, i, j)*qvs(k, i, j)
              rvs_p = rvs*otemp*(lvap(k, i, j)*otemp*orv - 1._wp)
              rvs_pp = rvs * (otemp*(lvap(k, i, j)*otemp*oRv - 1._wp) * &
                otemp*(lvap(k, i, j)*otemp*oRv - 1._wp) + &
                (-2._wp*lvap(k, i, j)*otemp*otemp*otemp*oRv) + otemp*otemp)
              gamsc = lvap(k, i, j)*diffu(k, i, j)/tcond(k, i, j) * rvs_p
              alphsc = 0.5_wp*(gamsc/(1._wp+gamsc))*(gamsc/(1._wp+gamsc)) * &
                rvs_pp/rvs_p * rvs/rvs_p
              alphsc = max(1.e-9_wp, alphsc)
              xsat = min(-1.e-9_wp, ssatw(k, i, j))
              t1_evap = 2._wp*pi*(1.0_wp - alphsc*xsat + 2._wp*alphsc*alphsc*xsat*xsat - &
                5._wp*alphsc*alphsc*alphsc*xsat*xsat*xsat) / (1._wp+gamsc)

              !> @note
              !> rain evaporation rapidly eliminates near zero values when low humidity (<95%)
              !> @endnotes
              if (qv(k, i, j)/qvs(k, i, j) < 0.95_wp .and. rr(k, i, j)*orho <= 1.e-8_wp) then
                tend%prv_rev(k, i, j) = rr(k, i, j)*orho*odt
              else
                lamr = 1._dp/ilamr(k, i, j)
                n0_r = nr(k, i, j)*org2*lamr**cre(2)
                tend%prv_rev(k, i, j) = t1_evap*diffu(k, i, j)*(-ssatw(k, i, j))*n0_r*rvs * &
                  (t1_qr_ev*ilamr(k, i, j)**cre(10) + t2_qr_ev*vsc2(k, i, j)*rhof2(k, i, j)* &
                  ((lamr+0.5*fv_r)**(-cre(11))))
                rate_max = min((rr(k, i, j)*orho*odt), &
                  (qvs(k, i, j)-qv(k, i, j))*odt)
                tend%prv_rev(k, i, j) = min(real(rate_max, kind=dp), tend%prv_rev(k, i, j)*orho)
                if (tend%prr_gml(k, i, j) > 0._dp) then
                  eva_factor = min(1._wp, 0.01_wp+(0.99_wp-0.01_wp)*(tempc/20._wp))
                  tend%prv_rev(k, i, j) = tend%prv_rev(k, i, j)*eva_factor
                endif
              endif
              tend%pnr_rev(k, i, j) = min(real(nr(k, i, j)*0.99*orho*odt, kind=dp),  &
                tend%prv_rev(k, i, j) * nr(k, i, j)/rr(k, i, j))
            endif
          endif
        enddo
        endif
      enddo
    enddo
    !$acc end parallel
  end subroutine rain_evaporation


  subroutine freeze_cloud_melt_ice(kts, kte, its, ite, jts, jte, temp, rho, ocp, lvap, qi3d, ni3d, qiten, niten, &
      qc3d, qcten, ncten, tten, dt, odt, column_mp_active, nc3d, ncsave)
    ! freezes all cloud water and melts all cloud ice instantly given the temperature (horizontal tile)
    !! OpenACC: parallel over columns; vertical k is sequential (no vertical recurrence).
    use module_mp_tempo_params, only : t0, lfus, lsub, hgfrz, nt_c_l

    integer, intent(in) :: kts, kte, its, ite, jts, jte
    real(wp), intent(in) :: dt, odt
    real(wp), dimension(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(in) :: temp, rho, ocp, lvap, qi3d, ni3d, qc3d
    real(wp), dimension(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(inout) :: qiten, niten, qcten, ncten, tten
    real(wp), dimension(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(in), optional :: nc3d, ncsave
    logical, dimension(TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(in), optional :: column_mp_active
    integer :: i, j, k
    real(wp) :: xri, xrc, lfus2, xnc
    logical :: use_cmp, compute_col, use_nc3d, use_ncsave

    use_cmp = present(column_mp_active)
    use_nc3d = present(nc3d)
    use_ncsave = present(ncsave)

    !$acc parallel
    !$acc loop gang vector collapse(2) private(xri, xrc, lfus2, xnc, compute_col)
    do j = TEMPO_JTS, TEMPO_JTE
      do i = TEMPO_ITS, TEMPO_ITE
        compute_col = .true.
        if (use_cmp) then
          if (.not. column_mp_active(i, j)) compute_col = .false.
        endif
        if (compute_col) then
          !$acc loop seq
          do k = kts, kte
            ! instantly melt all cloud ice
            xri = max(0._wp, qi3d(k, i, j)+qiten(k, i, j)*dt)
            if ((temp(k, i, j) > t0) .and. (xri > 0._wp)) then
              qcten(k, i, j) = qcten(k, i, j) + xri*odt
              ncten(k, i, j) = ncten(k, i, j) + ni3d(k, i, j)*odt
              qiten(k, i, j) = qiten(k, i, j) - xri*odt
              niten(k, i, j) = -ni3d(k, i, j)*odt
              tten(k, i, j) = tten(k, i, j) - lfus*ocp(k, i, j)*xri*odt
            endif
            ! instantly freeze all cloud water
            xrc = max(0._wp, qc3d(k, i, j)+qcten(k, i, j)*dt)
            if ((temp(k, i, j) < hgfrz) .and. (xrc > 0._wp)) then
              lfus2 = lsub - lvap(k, i, j)
              if (use_nc3d) then
                xnc = nc3d(k, i, j) + ncten(k, i, j)*dt
              else if (use_ncsave) then
                xnc = ncsave(k, i, j)/rho(k, i, j) + ncten(k, i, j)*dt
              else
                xnc = nt_c_l/rho(k, i, j) + ncten(k, i, j)*dt
              endif
              qiten(k, i, j) = qiten(k, i, j) + xrc*odt
              niten(k, i, j) = niten(k, i, j) + xnc*odt
              qcten(k, i, j) = qcten(k, i, j) - xrc*odt
              ncten(k, i, j) = ncten(k, i, j) - xnc*odt
              tten(k, i, j) = tten(k, i, j) + lfus2*ocp(k, i, j)*xrc*odt
            endif
          enddo
        endif
      enddo
    enddo
    !$acc end parallel
  end subroutine freeze_cloud_melt_ice


  function koop_nucleation(temp, satw, naero, dt) result(nuc)
  !$acc routine seq
    !! aqueous solution freezing of water from
    !! [Koop et al. (2000)](https://doi.org/10.1038/35020537)
    !! newer research suggests that the freezing rate should be lower
    !! than original paper, so J_rate is reduced by two orders of magnitude

    real(wp), intent(in) :: temp, satw, naero, dt
    real(wp) :: xni, mu_diff, a_w_i, delta_aw, log_j_rate, j_rate, prob_h
    real(wp) :: nuc

    xni = 0.0_wp

    mu_diff = 210368._wp + (131.438_wp*temp) - &
      (3.32373e6_wp/temp) - (41729.1_wp*log(temp))
    a_w_i = exp(mu_diff/(r_uni*temp))
    delta_aw = satw - a_w_i

    log_j_rate = -906.7_wp + (8502._wp*delta_aw) - &
      (26924._wp*delta_aw*delta_aw) + (29180._wp*delta_aw*delta_aw*delta_aw)
    log_j_rate = min(20._wp, log_j_rate)
    j_rate = 10._wp**log_j_rate ! cm-3 s-1
    prob_h = min(1._wp-exp(-j_rate*ar_volume*dt), 1._wp)
    if (prob_h > 0._wp) then
      xni = min(prob_h*naero, 1000.e3_wp)
    endif
    nuc = max(0._wp, xni)
  end function koop_nucleation


  function activate_cloud_number(temp, w1d, nwfa, land) result(activ)
  !$acc routine seq
    !! calculations numer of cloud droplets activated

    real(wp), intent(in) :: temp, w1d, nwfa
    integer, intent(in), optional :: land
    real(wp) :: n_local, w_local
    real(wp) :: a, b, c, d, t, u, x1, x2, y1, y2, nx, wy, nuc_frac
    real(wp) :: lower_lim_nuc_frac
    integer :: i, j, k, l, m, n
    real(wp) :: activ

    ! index for number of aerosols
    n_local = nwfa * 1.e-6_wp
    if (n_local >= ta_na(ntb_arc)) then
      n_local = ta_na(ntb_arc) - 1.0_wp
    elseif (n_local <= ta_na(1)) then
      n_local = ta_na(1) + 1.0_wp
    endif
    nindex: do n = 2, ntb_arc
      if (n_local >= ta_na(n-1) .and. n_local < ta_na(n)) exit nindex
    enddo nindex
    i = n
    x1 = log(ta_na(i-1))
    x2 = log(ta_na(i))

    ! index for vertical velocity
    w_local = w1d
    if (w_local >= ta_ww(ntb_arw)) then
        w_local = ta_ww(ntb_arw) - 1.0_wp
    elseif (w_local <= ta_ww(1)) then
        w_local = ta_ww(1) + 0.001_wp
    endif
    windex: do n = 2, ntb_arw
      if (w_local >= ta_ww(n-1) .and. w_local < ta_ww(n)) exit windex
    enddo windex
    j = n
    y1 = log(ta_ww(j-1))
    y2 = log(ta_ww(j))

    k = max(1, min(nint((temp - ta_tk(1))*0.1_wp) + 1, ntb_art))

    ! the next two values are indexes of mean aerosol radius and
    ! hygroscopicity and are currently constant 
    !> @todo
    !> separation tiny size sulfates from larger sea salts
    !> @endtodo
    l = 3
    m = 2

    !> @note
    !> there is a lower limit set for activation over water to improve cloud coverage
    !> @endnote
    lower_lim_nuc_frac = 0.
    if (present(land)) then
      if (land == 1) then ! land
        lower_lim_nuc_frac = 0.
      elseif (land == 0) then ! not land (water/ice)
        lower_lim_nuc_frac = 0.15
      else
        lower_lim_nuc_frac = 0.15 ! catch-all for anything else
      endif
    endif
    
    a = tnccn_act(i-1,j-1,k,l,m)
    b = tnccn_act(i,j-1,k,l,m)
    c = tnccn_act(i,j,k,l,m)
    d = tnccn_act(i-1,j,k,l,m)
    nx = log(n_local)
    wy = log(w_local)
    t = (nx-x1)/(x2-x1)
    u = (wy-y1)/(y2-y1)

    nuc_frac = (1.0_wp-t)*(1.0_wp-u)*a + t*(1.0_wp-u)*b + t*u*c + (1.0_wp-t)*u*d
    nuc_frac = max(nuc_frac, lower_lim_nuc_frac)

    activ = nwfa*nuc_frac
  end function activate_cloud_number


  subroutine warm_rain(kts, kte, its, ite, jts, jte, rhof, l_qc, rc, nc, ilamc, mvd_c, l_qr, rr, nr, mvd_r, tend, odt, column_mp_active)
    !! computes warm-rain process rates -- condensation/evaporation happen later (horizontal tile)
    use module_mp_tempo_params, only : d0r, d0c, r1, nbr, t_efrw, &
      t1_qr_qc, mu_r, am_r, ccg, obmr, ocg2, dr, org2, cre, fv_r, &
      autocon_nr_factor

    integer, intent(in) :: kts, kte, its, ite, jts, jte
    real(wp), intent(in) :: odt
    real(wp), dimension(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(in) :: rhof, mvd_r, mvd_c, rr, nr, rc, nc
    real(dp), dimension(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(in) :: ilamc
    logical, dimension(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(in) :: l_qc, l_qr
    type(ty_tend), intent(inout) :: tend
    logical, dimension(TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(in), optional :: column_mp_active

    real(dp) :: lamr, lamc, n0_r
    real(wp) :: ef_rr, dc_g, dc_b, xdc, zeta1, zeta, taud, tau, ef_rw
    integer :: i, j, k, nu_c, idx
    logical :: active_col

    !$acc parallel
    !$acc loop gang collapse(2) private(active_col)
    do j = TEMPO_JTS, TEMPO_JTE
      do i = TEMPO_ITS, TEMPO_ITE
        active_col = .true.
        if (present(column_mp_active)) active_col = column_mp_active(i, j)
        if (active_col) then
        !> @note
        !> rain self-collection is from
        !> [Seifert and Beheng (2001)](https://doi.org/10.1016/S0169-8095(01)00126-0)
        !> and drop break-up follows
        !> [Verlinde and Cotton (1993)](https://doi.org/10.1175/1520-0493(1993)121<2776:FMOONC>2.0.CO;2)
        !$acc loop vector private(lamr, lamc, n0_r, ef_rr, dc_g, dc_b, xdc, zeta1, zeta, taud, tau, ef_rw, nu_c, idx)
        do k = kts, kte
          if (l_qr(k, i, j)) then
            if (mvd_r(k, i, j) > d0r) then
              ef_rr = max(-0.1_wp, 1.0_wp - exp(2300.0_wp*(mvd_r(k, i, j)-1950.0e-6_wp)))
              tend%pnr_rcr(k, i, j) = ef_rr * 2.0_wp*nr(k, i, j)*rr(k, i, j)
            endif
          endif

          if (l_qc(k, i, j)) then
            if (rc(k, i, j) > 0.01e-3_wp) then
              nu_c = get_nuc(nc(k, i, j))
              lamc = 1._dp / ilamc(k, i, j)
              xdc = max(d0c*1.e6_wp, ((rc(k, i, j)/(am_r*nc(k, i, j)))**obmr) * 1.e6_wp)
              dc_g = ((ccg(3,nu_c)*ocg2(nu_c))**obmr / lamc) * 1.e6_wp
              dc_b = (xdc*xdc*xdc*dc_g*dc_g*dc_g - xdc*xdc*xdc*xdc*xdc*xdc) &
                  **(1._wp/6._wp)
              zeta1 = 0.5_wp*((6.25e-6_wp*xdc*dc_b*dc_b*dc_b - 0.4_wp) &
                  + abs(6.25e-6_wp*xdc*dc_b*dc_b*dc_b - 0.4_wp))
              zeta = 0.027_wp*rc(k, i, j)*zeta1
              taud = 0.5_wp*((0.5_wp*dc_b - 7.5_wp) + abs(0.5_wp*dc_b - 7.5_wp)) + r1
              tau = 3.72_wp/(rc(k, i, j)*taud)
              tend%prr_wau(k, i, j) = zeta/tau
              tend%prr_wau(k, i, j) = min(real(rc(k, i, j)*odt, kind=dp), &
                tend%prr_wau(k, i, j))
              tend%pnr_wau(k, i, j) = tend%prr_wau(k, i, j) / (am_r*nu_c*autocon_nr_factor*d0r*d0r*d0r)
              tend%pnc_wau(k, i, j) = min(real(nc(k, i, j)*odt, kind=dp), &
                tend%prr_wau(k, i, j) / (am_r*mvd_c(k, i, j)*mvd_c(k, i, j)*mvd_c(k, i, j)))
            endif
          endif

          if (l_qr(k, i, j) .and. l_qc(k, i, j)) then
            if (mvd_r(k, i, j) > d0r .and. mvd_c(k, i, j) > d0c) then
              lamr = (3.0_dp + mu_r + 0.672_dp) / mvd_r(k, i, j)
              idx = 1 + int(nbr*log(real(mvd_r(k, i, j)/dr(1), kind=dp)) / &
                log(real(dr(nbr)/dr(1), kind=dp)))
              idx = min(idx, nbr)
              ef_rw = t_efrw(idx, int(mvd_c(k, i, j)*1.e6_wp))
              n0_r = nr(k, i, j)*org2*lamr**cre(2)
              tend%prr_rcw(k, i, j) = rhof(k, i, j)*t1_qr_qc*ef_rw*rc(k, i, j)*n0_r * &
                ((lamr+fv_r)**(-cre(9)))
              tend%prr_rcw(k, i, j) = min(real(rc(k, i, j)*odt, kind=dp), tend%prr_rcw(k, i, j))
              tend%pnc_rcw(k, i, j) = rhof(k, i, j)*t1_qr_qc*ef_rw*nc(k, i, j)*n0_r * &
                ((lamr+fv_r)**(-cre(9)))
              tend%pnc_rcw(k, i, j) = min(real(nc(k, i, j)*odt, kind=dp), tend%pnc_rcw(k, i, j))
            endif
          endif
        enddo
        endif
      enddo
    enddo
    !$acc end parallel
  end subroutine warm_rain


  subroutine riming(kts, kte, its, ite, jts, jte, temp, rhof, visco, l_qc, rc, nc, ilamc, mvd_c, &
    l_qs, rs, smo0, smob, smoc, smoe, vtboost, l_qg, rg, ng, ilamg, idx, tend, odt, column_mp_active)
    !! snow and graupel riming (horizontal tile)
    use module_mp_tempo_params, only : d0c, d0s, nbs, ds, t_efsw, t1_qs_qc, &
      r_g, bm_g, mu_g, av_g, cgg, ogg3, bv_g, rho_w, t0, d0g, pi, cge, ogg2, &
      rime_threshold, rime_conversion, av_s, bv_s, rho_s, xm0i, eps, fv_s, &
      meters3_to_liters

    integer, intent(in) :: kts, kte, its, ite, jts, jte
    real(wp), intent(in) :: odt
    type(ty_tend), intent(inout) :: tend
    real(wp), dimension(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(in) :: rhof, visco, temp, rc, nc, rs, rg, ng
    real(wp), dimension(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(in) :: mvd_c
    real(dp), dimension(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(in) :: smo0, smob, smoc, smoe, ilamg, ilamc
    logical, dimension(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(in) :: l_qc, l_qs, l_qg
    integer, dimension(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(in) :: idx
    real(wp), dimension(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(out) :: vtboost
    logical, dimension(TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(in), optional :: column_mp_active

    real(dp) :: xds, xdg, n0_g, lamc
    real(wp) :: ef_sw, vtg, stoke_g, const_ri, tempc, rime_dens, ef_gw
    real(wp) :: t1_qg_qc, r_frac, g_frac, vts, tf, snow_dens_frac
    integer :: i, j, k, idxs, nu_c
    logical :: active_col

    !$acc parallel
    !$acc loop gang vector collapse(2) private(active_col)
    do j = TEMPO_JTS, TEMPO_JTE
      do i = TEMPO_ITS, TEMPO_ITE
        active_col = .true.
        if (present(column_mp_active)) active_col = column_mp_active(i, j)
        if (active_col) then
        !$acc loop vector private(xds, xdg, n0_g, lamc, ef_sw, vtg, stoke_g, const_ri, tempc, rime_dens, &
        !$acc& ef_gw, t1_qg_qc, r_frac, g_frac, vts, tf, snow_dens_frac, idxs, nu_c)
        do k = kts, kte
      tempc = temp(k, i, j) - t0
      vtboost(k, i, j) = 1._wp
      if (l_qc(k, i, j) .and. l_qs(k, i, j)) then
        nu_c = get_nuc(nc(k, i, j))
        lamc = 1._dp / ilamc(k, i, j)
        xds = smoc(k, i, j) / smob(k, i, j)
        if ((mvd_c(k, i, j) > d0c) .and. (xds > d0s)) then
          !> @note
          !> snow collecting cloud water - assume dc << ds and vtc \(\approx 0\)
          idxs = 1 + int(nbs*log(real(xds/ds(1), kind=dp)) / log(real(ds(nbs)/ds(1), kind=dp)))
          idxs = min(idxs, nbs)
          ef_sw = t_efsw(idxs, int(mvd_c(k, i, j)*1.e6_wp))
          tend%prs_scw(k,i,j) = rhof(k, i, j)*t1_qs_qc*ef_sw*rc(k, i, j)*smoe(k, i, j)
          tend%prs_scw(k,i,j) = min(real(rc(k, i, j)*odt, kind=dp), tend%prs_scw(k,i,j))
          tend%pnc_scw(k,i,j) = rhof(k, i, j)*t1_qs_qc*ef_sw*nc(k, i, j)*smoe(k, i, j)
          tend%pnc_scw(k,i,j) = min(real(nc(k, i, j)*odt, kind=dp), tend%pnc_scw(k,i,j))

          !>
          !> at temperatures below melting, if the riming rate is greater than the depositional
          !> growth rate for snow by a factor rime_threshold, convert a portion of rimed snow
          !> to graupel
          if (temp(k, i, j) < t0) then
            if (tend%prs_scw(k,i,j) > rime_threshold*tend%prs_sde(k,i,j) .and. &
              tend%prs_sde(k,i,j) > eps) then
              r_frac = min(30.0_dp, tend%prs_scw(k,i,j)/tend%prs_sde(k,i,j))
              g_frac = min(rime_conversion, 0.15_wp + (r_frac-2._wp)*.028_wp)
              vtboost(k, i, j) = min(1.5_wp, 1.1_wp + (r_frac-2.)*.016_wp)
              tend%prg_scw(k,i,j) = g_frac*tend%prs_scw(k,i,j)
              tend%png_scw(k,i,j) = tend%prg_scw(k,i,j)*smo0(k, i, j)/rs(k, i, j)
              vts = av_s*xds**bv_s * exp(-fv_s*xds)
              const_ri = -1._wp*(mvd_c(k, i, j)*0.5e6_wp)*vts/min(-0.1_wp,tempc)
              const_ri = max(0.1_wp, min(const_ri, 10._wp))
              rime_dens = (0.051_wp + 0.114_wp*const_ri - 0.0055_wp*const_ri*const_ri)*1000._wp
              if(rime_dens < 150._wp) then
                g_frac = 0._wp
                tend%prg_scw(k,i,j) = 0._dp
                tend%png_scw(k,i,j) = 0._dp
              endif
              snow_dens_frac = min(1._wp, max(0._wp, rs(k, i, j)*odt / &
                (rs(k, i, j)*odt + tend%prg_scw(k,i,j))))
              tend%pbg_scw(k,i,j) = meters3_to_liters*tend%prg_scw(k,i,j) / &
                (rho_s * snow_dens_frac + rime_dens * (1._wp-snow_dens_frac))
              tend%prs_scw(k,i,j) = (1._wp - g_frac)*tend%prs_scw(k,i,j)
            endif
          endif
        endif
      endif

      if (l_qc(k, i, j) .and. l_qg(k, i, j)) then
        !>
        !> graupel collecting cloud water - assume dc << dg and vtc \(\approx 0\)
        if (rg(k, i, j) >= r_g(1) .and. mvd_c(k, i, j) > d0c) then
          xdg = (bm_g + mu_g + 1._wp) * ilamg(k, i, j)
          vtg = rhof(k, i, j)*av_g(idx(k, i, j))*cgg(6,idx(k, i, j))*ogg3 * ilamg(k, i, j)**bv_g(idx(k, i, j))
          stoke_g = mvd_c(k, i, j)*mvd_c(k, i, j)*vtg*rho_w/(9._wp*visco(k, i, j)*xdg)
          !>
          !> rime density formula is from
          !> [Cober and List (1993)](https://doi.org/10.1175/1520-0469(1993)050<1591:MOTHAM>2.0.CO;2)
          const_ri = -1._wp*(mvd_c(k, i, j)*0.5e6_wp)*vtg/min(-0.1_wp, tempc)
          const_ri = max(0.1_wp, min(const_ri, 10._wp))
          rime_dens = (0.051_wp + 0.114_wp*const_ri - 0.0055_wp*const_ri*const_ri)*1000._wp
          if (xdg > d0g) then
            if (stoke_g >= 0.4_wp .and. stoke_g <= 10._wp) then
              ef_gw = 0.55_wp*log10(2.51_wp*stoke_g)
            elseif (stoke_g < 0.4_wp) then
              ef_gw = 0.0_wp
            elseif (stoke_g > 10._wp) then
              ef_gw = 0.77_wp
            endif
            !>
            !> hail size increases below the melting level so the collection efficiency
            !> is reduced (proxy for shedding of collected cloud water)
            if (temp(k, i, j) > t0) ef_gw = ef_gw*0.1_wp
            t1_qg_qc = pi*.25_wp*av_g(idx(k, i, j)) * cgg(9,idx(k, i, j))
            n0_g = ng(k, i, j)*ogg2*(1._wp/ilamg(k, i, j))**cge(2,1)
            tend%prg_gcw(k,i,j) = rhof(k, i, j)*t1_qg_qc*ef_gw*rc(k, i, j)* &
              n0_g*ilamg(k, i, j)**cge(9,idx(k, i, j))
            tend%pnc_gcw(k,i,j) = rhof(k, i, j)*t1_qg_qc*ef_gw*nc(k, i, j)* &
              n0_g*ilamg(k, i, j)**cge(9,idx(k, i, j))
            tend%pnc_gcw(k,i,j) = min(real(nc(k, i, j)*odt, kind=dp), tend%pnc_gcw(k,i,j))
            if (temp(k, i, j) < t0) tend%pbg_gcw(k,i,j) = meters3_to_liters*tend%prg_gcw(k,i,j)/rime_dens

            if (temp(k, i, j) < t0) then
              !>
              !> rime splintering is from
              !> [Hallet and Mossop (1974)](https://doi.org/10.1038/249026a0)
              !> @endnote
              if (tend%prg_gcw(k,i,j) > eps .and. tempc > -8._wp) then
                tf = 0._wp
                if (tempc >= -5._wp .and. tempc < -3._wp) then
                  tf = 0.5_wp*(-3.0_wp - tempc)
                elseif (tempc > -8._wp .and. tempc < -5._wp) then
                  tf = 0.33333333_wp*(8._wp + tempc)
                endif
                tend%pni_ihm(k,i,j) = 3.5e8_wp*tf*tend%prg_gcw(k,i,j)
                tend%pri_ihm(k,i,j) = xm0i*tend%pni_ihm(k,i,j)
                tend%prs_ihm(k,i,j) = tend%prs_scw(k,i,j)/(tend%prs_scw(k,i,j)+tend%prg_gcw(k,i,j)) * &
                  tend%pri_ihm(k,i,j)
                tend%prg_ihm(k,i,j) = tend%prg_gcw(k,i,j)/(tend%prs_scw(k,i,j)+tend%prg_gcw(k,i,j)) * &
                  tend%pri_ihm(k,i,j)
              endif
            endif
          endif
        endif
      endif
        enddo
        endif
      enddo
    enddo
    !$acc end parallel
  end subroutine riming


  subroutine get_snow_table_index(rs, idx_s)
  !$acc routine seq
    !! get snow table index from snow mass

    real(wp), intent(in) :: rs
    integer :: nis, nn, n
    integer, intent(out) :: idx_s

    nis = nint(log10(rs))
    do_loop_rs: do nn = nis-1, nis+1
      n = nn
      if ((rs/10._wp**nn) >= 1._wp .and. (rs/10._wp**nn) < 10._wp) exit do_loop_rs
    enddo do_loop_rs
    idx_s = int(rs/10._wp**n) + 10*(n-nis2) - (n-nis2)
    idx_s = max(1, min(idx_s, ntb_s))
  end subroutine get_snow_table_index


  subroutine get_temperature_table_index(tempk, idx_t)
  !$acc routine seq
    !! get temperature table index
    
    real(wp), intent(in) :: tempk
    real(wp) :: tempc
    integer, intent(out) :: idx_t

    tempc = tempk - t0
    idx_t = int((tempc-2.5_wp)/5._wp) - 1
    idx_t = max(1, -idx_t)
    idx_t = min(idx_t, ntb_t)
  end subroutine get_temperature_table_index


  subroutine get_rain_table_index(rr, ilamr, idx_r, idx_r1)
  !$acc routine seq
    !! get rain table indices from rain mass and lambda

    real(wp), intent(in) :: rr
    real(dp), intent(in) :: ilamr
    real(dp) :: lamr, lam_exp, n0_exp
    integer :: nir, nn, n
    integer, intent(out) :: idx_r, idx_r1

    nir = nint(log10(rr))
    do_loop_rr: do nn = nir-1, nir+1
      n = nn
      if ((rr/10._wp**nn) >= 1._wp .and. (rr/10._wp**nn) < 10._wp) exit do_loop_rr
    enddo do_loop_rr
    idx_r = int(rr/10._wp**n) + 10*(n-nir2) - (n-nir2)
    idx_r = max(1, min(idx_r, ntb_r))

    lamr = 1./ilamr
    lam_exp = lamr * (crg(3)*org2*org1)**bm_r
    n0_exp = org1*rr/am_r * lam_exp**cre(1)
    nir = nint(log10(real(n0_exp, kind=dp)))
    do_loop_nr: do nn = nir-1, nir+1
      n = nn
      if ((n0_exp/10._wp**nn) >= 1._wp .and. (n0_exp/10._wp**nn) < 10._wp) exit do_loop_nr
    enddo do_loop_nr
    idx_r1 = int(n0_exp/10._wp**n) + 10*(n-nir3) - (n-nir3)
    idx_r1 = max(1, min(idx_r1, ntb_r1))
  end subroutine get_rain_table_index


  subroutine get_graupel_table_index(rg, ilamg, idx, idx_g, idx_g1)
  !$acc routine seq
    !! get graupel table indices from graupel mass, lambda, and density index

    real(wp), intent(in) :: rg
    real(dp), intent(in) :: ilamg
    integer, intent(in) :: idx
    real(dp) :: lamg, lam_exp, n0_exp
    integer :: nig, nn, n
    integer, intent(out) :: idx_g, idx_g1

    nig = nint(log10(rg))
    do_loop_rg: do nn = nig-1, nig+1
      n = nn
      if ( (rg/10._wp**nn) >= 1._wp .and. (rg/10._wp**nn).lt.10._wp) exit do_loop_rg
    enddo do_loop_rg
    idx_g = int(rg/10._wp**n) + 10*(n-nig2) - (n-nig2)
    idx_g = max(1, min(idx_g, ntb_g))

    lamg = 1./ilamg
    lam_exp = lamg * (cgg(3,1)*ogg2*ogg1)**bm_g
    n0_exp = ogg1*rg/am_g(idx) * lam_exp**cge(1,1)
    nig = nint(log10(real(n0_exp, kind=dp)))
    do_loop_ng: do nn = nig-1, nig+1
      n = nn
      if ( (n0_exp/10._wp**nn) >= 1._wp .and. (n0_exp/10._wp**nn) < 10._wp) exit do_loop_ng
    enddo do_loop_ng
    idx_g1 = int(n0_exp/10._wp**n) + 10*(n-nig3) - (n-nig3)
    idx_g1 = max(1, min(idx_g1, ntb_g1))
  end subroutine get_graupel_table_index


  subroutine get_cloud_table_index(rc, nc, idx_c, idx_n)
  !$acc routine seq
    !! get cloud table index from mass and number

    real(wp), intent(in) :: rc, nc
    integer, intent(out) :: idx_c, idx_n
    integer :: nic, nn, n

    nic = nint(log10(rc))
    do_loop_rc: do nn = nic-1, nic+1
      n = nn
      if ( (rc/10._wp**nn) >= 1._wp .and. (rc/10._wp**nn) < 10._wp) exit do_loop_rc
    enddo do_loop_rc
    idx_c = int(rc/10._wp**n) + 10*(n-nic2) - (n-nic2)
    idx_c = max(1, min(idx_c, ntb_c))
          
    idx_n = nint(1._wp + real(nbc, kind=wp) * log(real(nc/t_nc(1), kind=dp)) / nic1)
    idx_n = max(1, min(idx_n, nbc))
  end subroutine get_cloud_table_index


  subroutine get_ice_table_index(ri, ni, idx_i, idx_i1)
  !$acc routine seq
    !! get ice table index from mass and number

    real(wp), intent(in) :: ri, ni
    integer, intent(out) :: idx_i, idx_i1
    integer :: nii, nn, n

    nii = nint(log10(ri))
    do_loop_ri: do nn = nii-1, nii+1
      n = nn
      if ( (ri/10._wp**nn) >= 1._wp .and. (ri/10._wp**nn) < 10._wp) exit do_loop_ri
    enddo do_loop_ri
    idx_i = int(ri/10._wp**n) + 10*(n-nii2) - (n-nii2)
    idx_i = max(1, min(idx_i, ntb_i))
  
    nii = nint(log10(ni))
    do_loop_ni: do nn = nii-1, nii+1
      n = nn
      if ( (ni/10._wp**nn) >= 1._wp .and. (ni/10._wp**nn) < 10._wp) exit do_loop_ni
    enddo do_loop_ni
    idx_i1 = int(ni/10._wp**n) + 10*(n-nii3) - (n-nii3)
    idx_i1 = max(1, min(idx_i1, ntb_i1))
  end subroutine get_ice_table_index


  subroutine rain_snow_rain_graupel(kts, kte, its, ite, jts, jte, temp, l_qr, rr, nr, ilamr, l_qs, rs, &
      l_qg, rg, ng, ilamg, idx, tend, odt, column_mp_active)
    !! calculates rain-snow and rain-graupel collection (horizontal tile)
    use module_mp_tempo_params, only : t0, r_r, r_s, r_g, rho_i, rho_g, meters3_to_liters, &
      tmr_racs2, tcr_sacr2, tmr_racs1, tcr_sacr1, tms_sacr1, tcs_racs1, &
      tnr_sacr1, tnr_sacr2, tnr_racs1, tnr_racs2, &
      tcr_gacr, tmr_racg, tcg_racg, tnr_gacr, tnr_racg

    integer, intent(in) :: kts, kte, its, ite, jts, jte
    real(wp), intent(in) :: odt
    type(ty_tend), intent(inout) :: tend
    real(wp), dimension(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(in) :: temp, rr, rs, rg, nr, ng
    real(dp), dimension(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(in) :: ilamr, ilamg
    logical, dimension(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(in) :: l_qr, l_qs, l_qg
    integer, dimension(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(in) :: idx
    logical, dimension(TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(in), optional :: column_mp_active
    integer :: i, j, k, idx_r, idx_r1, idx_s, idx_t, idx_g, idx_g1
    logical :: active_col

    !$acc parallel
    !$acc loop gang collapse(2) private(active_col)
    do j = TEMPO_JTS, TEMPO_JTE
      do i = TEMPO_ITS, TEMPO_ITE
        active_col = .true.
        if (present(column_mp_active)) active_col = column_mp_active(i, j)
        if (active_col) then
        !$acc loop vector private(idx_r, idx_r1, idx_s, idx_t, idx_g, idx_g1)
        do k = kts, kte
          if (l_qr(k, i, j) .and. l_qs(k, i, j)) then
            if (rr(k, i, j) >= r_r(1) .and. rs(k, i, j) >= r_s(1)) then
              call get_temperature_table_index(temp(k, i, j), idx_t)
              call get_rain_table_index(rr(k, i, j), ilamr(k, i, j), idx_r, idx_r1)
              call get_snow_table_index(rs(k, i, j), idx_s)
              if (temp(k, i, j) < t0) then
                tend%prr_rcs(k, i, j) = -(tmr_racs2(idx_s,idx_t,idx_r1,idx_r) &
                  + tcr_sacr2(idx_s,idx_t,idx_r1,idx_r) &
                  + tmr_racs1(idx_s,idx_t,idx_r1,idx_r) &
                  + tcr_sacr1(idx_s,idx_t,idx_r1,idx_r))
                tend%prs_rcs(k, i, j) = tmr_racs2(idx_s,idx_t,idx_r1,idx_r) &
                  + tcr_sacr2(idx_s,idx_t,idx_r1,idx_r) &
                  - tcs_racs1(idx_s,idx_t,idx_r1,idx_r) &
                  - tms_sacr1(idx_s,idx_t,idx_r1,idx_r)
                tend%prg_rcs(k, i, j) = tmr_racs1(idx_s,idx_t,idx_r1,idx_r) &
                  + tcr_sacr1(idx_s,idx_t,idx_r1,idx_r) &
                  + tcs_racs1(idx_s,idx_t,idx_r1,idx_r) &
                  + tms_sacr1(idx_s,idx_t,idx_r1,idx_r)
                tend%prr_rcs(k, i, j) = max(real(-rr(k, i, j)*odt, kind=dp), tend%prr_rcs(k, i, j))
                tend%prs_rcs(k, i, j) = max(real(-rs(k, i, j)*odt, kind=dp), tend%prs_rcs(k, i, j))
                tend%prg_rcs(k, i, j) = min(real((rr(k, i, j)+rs(k, i, j))*odt, kind=dp), &
                  tend%prg_rcs(k, i, j))
                tend%pnr_rcs(k, i, j) = tnr_racs1(idx_s,idx_t,idx_r1,idx_r) &
                  + tnr_racs2(idx_s,idx_t,idx_r1,idx_r) &
                  + tnr_sacr1(idx_s,idx_t,idx_r1,idx_r) &
                  + tnr_sacr2(idx_s,idx_t,idx_r1,idx_r)
                tend%pnr_rcs(k, i, j) = min(real(nr(k, i, j)*odt, kind=dp), tend%pnr_rcs(k, i, j))
                tend%png_rcs(k, i, j) = tend%pnr_rcs(k, i, j)
                tend%pbg_rcs(k, i, j) = meters3_to_liters*tend%prg_rcs(k, i, j)/rho_i
              else
                tend%prs_rcs(k, i, j) = -tcs_racs1(idx_s,idx_t,idx_r1,idx_r) &
                  - tms_sacr1(idx_s,idx_t,idx_r1,idx_r) &
                  + tmr_racs2(idx_s,idx_t,idx_r1,idx_r) &
                  + tcr_sacr2(idx_s,idx_t,idx_r1,idx_r)
                tend%prs_rcs(k, i, j) = max(real(-rs(k, i, j)*odt, kind=dp), tend%prs_rcs(k, i, j))
                tend%prr_rcs(k, i, j) = -tend%prs_rcs(k, i, j)
              endif
            endif
          endif

          if (l_qr(k, i, j) .and. l_qg(k, i, j)) then
            if (rr(k, i, j) >= r_r(1) .and. rg(k, i, j) >= r_g(1)) then
              call get_temperature_table_index(temp(k, i, j), idx_t)
              call get_rain_table_index(rr(k, i, j), ilamr(k, i, j), idx_r, idx_r1)
              call get_graupel_table_index(rg(k, i, j), ilamg(k, i, j), idx(k, i, j), idx_g, idx_g1)
              if (temp(k, i, j) < t0) then
                tend%prg_rcg(k, i, j) = tmr_racg(idx_g1,idx_g,idx(k, i, j),idx_r1,idx_r) &
                  + tcr_gacr(idx_g1,idx_g,idx(k, i, j),idx_r1,idx_r)
                tend%prg_rcg(k, i, j) = min(real(rr(k, i, j)*odt, kind=dp), tend%prg_rcg(k, i, j))
                tend%prr_rcg(k, i, j) = -tend%prg_rcg(k, i, j)
                tend%pnr_rcg(k, i, j) = tnr_racg(idx_g1,idx_g,idx(k, i, j),idx_r1,idx_r) &
                  + tnr_gacr(idx_g1,idx_g,idx(k, i, j),idx_r1,idx_r)
                tend%pnr_rcg(k, i, j) = min(real(nr(k, i, j)*odt, kind=dp), tend%pnr_rcg(k, i, j))
                tend%pbg_rcg(k, i, j) = meters3_to_liters*tend%prg_rcg(k, i, j)/rho_i
              else
                tend%prr_rcg(k, i, j) = tcg_racg(idx_g1,idx_g,idx(k, i, j),idx_r1,idx_r)
                tend%prr_rcg(k, i, j) = min(real(rg(k, i, j)*odt, kind=dp), tend%prr_rcg(k, i, j))
                tend%prg_rcg(k, i, j) = -tend%prr_rcg(k, i, j)
                tend%png_rcg(k, i, j) = tnr_racg(idx_g1,idx_g,idx(k, i, j),idx_r1,idx_r)
                tend%png_rcg(k, i, j) = min(real(ng(k, i, j)*odt, kind=dp), tend%png_rcg(k, i, j))
                tend%pbg_rcg(k, i, j) = meters3_to_liters*tend%prg_rcg(k, i, j)/rho_g(idx(k, i, j))
                !> @note
                !> adds explicit rain drop break-up due to collisions with graupel
                !> at temperatures above melting
                !> @endnote
                tend%pnr_rcg(k, i, j) = -1.5_wp*tnr_gacr(idx_g1,idx_g,idx(k, i, j),idx_r1,idx_r)
              endif
            endif
          endif
        enddo
        endif
      enddo
    enddo
    !$acc end parallel
  end subroutine rain_snow_rain_graupel


  subroutine ice_nucleation(kts, kte, its, ite, jts, jte, temp, rho, w1d, qv, qvsi, ssati, ssatw, &
      ni, smo0, rc, nc, rr, nr, ilamr, tend, dt, odt, column_mp_active, nifa, nwfa)
    !! ice nucleation (horizontal tile)
    use module_mp_tempo_params, only : r1, r_r, r_c, hgfrz, rho_i, xm0i, t0, &
      tpg_qrfz, tpi_qrfz, tni_qrfz, tnr_qrfz, tpi_qcfz, tni_qcfz, &
      demott_nuc_ssati, eps, icenuc_max, tno, ato, max_ni, meters3_to_liters

    integer, intent(in) :: kts, kte, its, ite, jts, jte
    real(wp), intent(in) :: dt, odt
    type(ty_tend), intent(inout) :: tend
    real(wp), dimension(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(in) :: qv, temp, rho, qvsi, rr, nr, rc, nc, w1d, &
      ssati, ssatw, ni
    real(dp), dimension(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(in) :: ilamr, smo0
    logical, dimension(TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(in), optional :: column_mp_active
    real(wp), dimension(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(in), optional :: nifa, nwfa
    real(wp) :: rate_max, tempc, xni, xnc
    integer :: i, j, k, idx_in, idx_r, idx_r1, idx_tc, idx_c, idx_n
    logical :: active_col

    !$acc parallel
    !$acc loop gang vector collapse(2) private(active_col)
    do j = TEMPO_JTS, TEMPO_JTE
      do i = TEMPO_ITS, TEMPO_ITE
        active_col = .true.
        if (present(column_mp_active)) active_col = column_mp_active(i, j)
        if (active_col) then
        !$acc loop vector private(rate_max, tempc, xni, xnc, idx_in, idx_r, idx_r1, idx_tc, idx_c, idx_n)
        do k = kts, kte
          if (temp(k, i, j) < t0) then
            tempc = temp(k, i, j) - t0
            idx_tc = max(1, min(nint(-tempc), 45))
            rate_max = (qv(k, i, j)-qvsi(k, i, j))*rho(k, i, j)*odt*0.999_wp
            if (present(nifa)) then
              xni = demott_nucleation(tempc, rho(k, i, j), nifa(k, i, j))
            else
              xni = 1._wp * 1000._wp ! 1 / Liter
            endif
            call get_in_table_index(xni, idx_in)

            if (rr(k, i, j) > r_r(1)) then
              call get_rain_table_index(rr(k, i, j), ilamr(k, i, j), idx_r, idx_r1)
              tend%prg_rfz(k, i, j) = tpg_qrfz(idx_r,idx_r1,idx_tc,idx_in)*odt
              tend%pri_rfz(k, i, j) = tpi_qrfz(idx_r,idx_r1,idx_tc,idx_in)*odt
              tend%pni_rfz(k, i, j) = tni_qrfz(idx_r,idx_r1,idx_tc,idx_in)*odt
              tend%pnr_rfz(k, i, j) = tnr_qrfz(idx_r,idx_r1,idx_tc,idx_in)*odt
              tend%prg_rfz(k, i, j) = min(real(rr(k, i, j)*odt, kind=dp), tend%prg_rfz(k, i, j))
              tend%pnr_rfz(k, i, j) = min(real(nr(k, i, j)*odt, kind=dp), tend%pnr_rfz(k, i, j))
              tend%png_rfz(k, i, j) = tend%pnr_rfz(k, i, j) * &
                max(min((10._wp**(-0.1_wp*w1d(k, i, j)) + 0.1_wp), 1._wp), 0.1_wp)
            elseif (rr(k, i, j) > r1 .and. temp(k, i, j) < hgfrz) then
              tend%pri_rfz(k, i, j) = rr(k, i, j)*odt
              tend%pni_rfz(k, i, j) = nr(k, i, j)*odt
            endif
            tend%pbg_rfz(k, i, j) = meters3_to_liters*tend%prg_rfz(k, i, j)/rho_i

            if (rc(k, i, j) > r_c(1)) then
              call get_cloud_table_index(rc(k, i, j), nc(k, i, j), idx_c, idx_n)
              tend%pri_wfz(k, i, j) = tpi_qcfz(idx_c,idx_n,idx_tc,idx_in)*odt
              tend%pri_wfz(k, i, j) = min(real(rc(k, i, j)*odt, kind=dp), tend%pri_wfz(k, i, j))
              tend%pni_wfz(k, i, j) = tni_qcfz(idx_c,idx_n,idx_tc,idx_in)*odt
              tend%pni_wfz(k, i, j) = min(real(nc(k, i, j)*odt, kind=dp), &
                tend%pri_wfz(k, i, j)/(2.0_dp*xm0i), tend%pni_wfz(k, i, j))
            elseif (rc(k, i, j) > r1 .and. temp(k, i, j) < hgfrz) then
              tend%pri_wfz(k, i, j) = rc(k, i, j)*odt
              tend%pni_wfz(k, i, j) = nc(k, i, j)*odt
            endif

            if ((ssati(k, i, j) >= demott_nuc_ssati) .or. (ssatw(k, i, j) > eps &
                .and. tempc < -20._wp)) then
              if (present(nifa)) then
                xnc = demott_nucleation(tempc, rho(k, i, j), nifa(k, i, j))
              else
                xnc = min(icenuc_max, tno*exp(ato*(t0-temp(k, i, j))))
              endif
              xni = ni(k, i, j) + (tend%pni_rfz(k, i, j)+tend%pni_wfz(k, i, j))*dt
              tend%pni_inu(k, i, j) = 0.5_wp*(xnc-xni + abs(xnc-xni))*odt
              tend%pri_inu(k, i, j) = min(real(rate_max, kind=dp), xm0i*tend%pni_inu(k, i, j))
              tend%pni_inu(k, i, j) = tend%pri_inu(k, i, j)/xm0i
            endif

            xni = smo0(k, i, j)+ni(k, i, j) + (tend%pni_rfz(k, i, j)+tend%pni_wfz(k, i, j)+tend%pni_inu(k, i, j))*dt
            if (present(nwfa)) then
              if ((xni <= max_ni) .and.(temp(k, i, j) < 238._wp) .and. (ssati(k, i, j) >= 0.4_wp)) then
                xnc = koop_nucleation(temp(k, i, j), ssatw(k, i, j), nwfa(k, i, j), dt)
                tend%pni_iha(k, i, j) = xnc*odt
                tend%pri_iha(k, i, j) = min(real(rate_max, kind=dp), xm0i*0.1_wp*tend%pni_iha(k, i, j))
                tend%pni_iha(k, i, j) = tend%pri_iha(k, i, j)/(xm0i*0.1_wp)
              endif
            endif
          endif
        enddo
        endif
      enddo
    enddo
    !$acc end parallel
  end subroutine ice_nucleation


  function demott_nucleation(tempc, rho, nifa) result(nuc)
  !$acc routine seq
    !! DeMott nucleation

    real(wp), intent(in) :: tempc, rho, nifa
    real(wp) :: xni, nifa_cc
    real(wp) :: nuc

    xni = 0._wp
    nifa_cc = max(0.5_wp, nifa*rho_not0*1.e-6_wp/rho)
    xni = (5.94e-5_wp*(-tempc)**3.33_wp) * (nifa_cc**((-0.0264_wp*(tempc))+0.0033_wp))
    xni = xni*rho/rho_not0 * 1000._wp
    nuc = max(0._wp, xni)
  end function demott_nucleation


  subroutine get_in_table_index(xni, idx_in)
  !$acc routine seq
    !! get ice nuclei table index

    real(wp), intent(in) :: xni
    integer, intent(out) :: idx_in
    integer :: niin, nn, n

    if (xni >  nt_in(1)) then
        niin = nint(log10(xni))
        do_loop_xni: do nn = niin-1, niin+1
          n = nn
          if ( (xni/10._wp**nn) >= 1._wp .and. (xni/10._wp**nn) < 10._wp) exit do_loop_xni
        enddo do_loop_xni
        idx_in = int(xni/10._wp**n) + 10*(n-niin2) - (n-niin2)
        idx_in = max(1, min(idx_in, ntb_in))
    else
        idx_in = 1
    endif
  end subroutine get_in_table_index


  subroutine get_t1_subl(kts, kte, rho, temp, qvsi, tcond, diffu, ssati, t1_subl)
  !$acc routine seq
    !! calculations thermodynamic term used in depositional growth and melting

    integer, intent(in) :: kts, kte
    real(wp), dimension(kts:kte), intent(in) :: rho, temp, qvsi, tcond, diffu, ssati
    real(wp) :: otemp, rvs, rvs_p, rvs_pp, gamsc, alphsc, xsat
    real(wp), dimension(kts:kte), intent(out) :: t1_subl
    integer :: k

    do k = kts, kte
      otemp = 1._wp/temp(k)
      rvs = rho(k)*qvsi(k)
      rvs_p = rvs*otemp*(lsub*otemp*orv - 1._wp)
      rvs_pp = rvs * (otemp*(lsub*otemp*orv - 1._wp) * otemp*(lsub*otemp*orv - 1._wp) + &
        (-2.*lsub*otemp*otemp*otemp*orv) + otemp*otemp)
      gamsc = lsub*diffu(k)/tcond(k) * rvs_p
      alphsc = 0.5_wp*(gamsc/(1._wp+gamsc))*(gamsc/(1._wp+gamsc)) * &
        rvs_pp/rvs_p * rvs/rvs_p
      alphsc = max(1.e-9_wp, alphsc)
      xsat = ssati(k)
      if (abs(xsat) < 1.e-9_wp) xsat = 0._wp
      t1_subl(k) = 4._wp*pi*(1._wp - alphsc*xsat + 2._wp*alphsc*alphsc*xsat*xsat - &
        5._wp*alphsc*alphsc*alphsc*xsat*xsat*xsat) / (1._wp+gamsc)
    end do  
  end subroutine get_t1_subl


  subroutine ice_processes(kts, kte, its, ite, jts, jte, rhof, rhof2, rho, w1d, temp, qv, qvsi, tcond, diffu, &
    vsc2, ssati, l_qi, ri, ni, ilami, l_qs, rs, smoe, smof, smo1, rr, nr, ilamr, &
    mvd_r, l_qg, rg, ng, ilamg, idx, tend, odt, column_mp_active)
    !! ice processes over horizontal tile (cloud ice, snow, graupel)
    use module_mp_tempo_params, only : t0, d0i, bm_i, mu_i, am_i, &
      c_sqrd, c_cube, oig1, cig, d0s, ntb_i, tpi_ide, tps_iaus, tni_iaus, &
      obmi, r_s, ef_si, t1_qs_qi, r_r, org2, cre, t1_qr_qi, t2_qr_qi, &
      fv_r, ef_ri, rho_i, t1_qs_sd, t2_qs_sd, eps, t1_qg_sd, &
      sc3, ogg2, cge, cgg, av_g, rho_w, rho_g

    integer, intent(in) :: kts, kte, its, ite, jts, jte
    real(wp), intent(in) :: odt
    type(ty_tend), intent(inout) :: tend
    logical, dimension(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(in) :: l_qi, l_qs, l_qg
    real(wp), dimension(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(in) :: rhof, rhof2, rho, w1d, ri, ni, rs, rr, nr, &
      temp, qv, qvsi, tcond, diffu, ssati, vsc2, mvd_r, rg, ng
    real(dp), dimension(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(in) :: ilami, smoe, smof, smo1, ilamr, ilamg
    integer, dimension(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(in) :: idx
    logical, dimension(TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(in), optional :: column_mp_active
    real(wp) :: xdi, xmi, oxmi, c_snow, rate_max, otemp, rvs, t2_qg_sd
    real(dp) :: lami, lamr, n0_r, n0_g
    integer :: i, j, k, idx_i, idx_i1
    real(wp), dimension(kts:kte) :: t1_subl
    logical :: active_col

    !$acc parallel
    !$acc loop gang vector collapse(2) private(active_col, t1_subl)
    do j = TEMPO_JTS, TEMPO_JTE
      do i = TEMPO_ITS, TEMPO_ITE
        active_col = .true.
        if (present(column_mp_active)) active_col = column_mp_active(i, j)
        if (active_col) then
        call get_t1_subl(kts, kte, rho(kts:kte, i, j), temp(kts:kte, i, j), qvsi(kts:kte, i, j), tcond(kts:kte, i, j), &
          diffu(kts:kte, i, j), ssati(kts:kte, i, j), t1_subl)

        !$acc loop vector private(xdi, xmi, oxmi, c_snow, rate_max, otemp, rvs, t2_qg_sd, lami, lamr, n0_r, n0_g, idx_i, idx_i1)
        do k = kts, kte
          otemp = 1._wp/temp(k, i, j)
          rvs = rho(k, i, j)*qvsi(k, i, j)
          rate_max = (qv(k, i, j)-qvsi(k, i, j))*rho(k, i, j)*odt*0.999_wp

          if (temp(k, i, j) < t0) then
            if (l_qi(k, i, j)) then
              call get_ice_table_index(ri(k, i, j), ni(k, i, j), idx_i, idx_i1)
              lami = 1._dp/ilami(k, i, j)
              xdi = max(real(d0i, kind=dp), (bm_i + mu_i + 1.) * ilami(k, i, j))
              xmi = am_i*xdi**bm_i
              oxmi = 1._wp/xmi
              tend%pri_ide(k,i,j) = c_cube*t1_subl(k)*diffu(k, i, j)*ssati(k, i, j)*rvs &
                *oig1*cig(5)*ni(k, i, j)*ilami(k, i, j)
              if (tend%pri_ide(k,i,j) < 0._dp) then
                tend%pri_ide(k,i,j) = max(real(-ri(k, i, j)*odt, kind=dp), &
                  tend%pri_ide(k,i,j), real(rate_max, kind=dp))
                tend%pni_ide(k,i,j) = tend%pri_ide(k,i,j)*oxmi
                tend%pni_ide(k,i,j) = max(real(-ni(k, i, j)*odt, kind=dp), tend%pni_ide(k,i,j))
              else
                tend%pri_ide(k,i,j) = min(tend%pri_ide(k,i,j), real(rate_max, kind=dp))
                tend%prs_ide(k,i,j) = (1.0_dp-tpi_ide(idx_i,idx_i1))*tend%pri_ide(k,i,j)
                tend%pri_ide(k,i,j) = tpi_ide(idx_i,idx_i1)*tend%pri_ide(k,i,j)
              endif

              if ((idx_i == ntb_i) .or. (xdi >  5.0_wp*d0s)) then
                tend%prs_iau(k,i,j) = ri(k, i, j)*.99_wp*odt
                tend%pni_iau(k,i,j) = ni(k, i, j)*.95_wp*odt
              elseif (xdi < 0.1_wp*d0s) then
                tend%prs_iau(k,i,j) = 0._dp
                tend%pni_iau(k,i,j) = 0._dp
              else
                tend%prs_iau(k,i,j) = tps_iaus(idx_i,idx_i1)*odt
                tend%prs_iau(k,i,j) = min(real(ri(k, i, j)*.99_wp*odt, kind=dp), tend%prs_iau(k,i,j))
                tend%pni_iau(k,i,j) = tni_iaus(idx_i,idx_i1)*odt
                tend%pni_iau(k,i,j) = min(real(ni(k, i, j)*.95_wp*odt, kind=dp), tend%pni_iau(k,i,j))
              endif

              lami = (am_i*cig(2)*oig1*ni(k, i, j)/ri(k, i, j))**obmi
              xdi = max(real(D0i, kind=dp), (bm_i + mu_i + 1.) * ilami(k, i, j))
              xmi = am_i*xDi**bm_i
              oxmi = 1./xmi
              if (rs(k, i, j) >= r_s(1)) then
                tend%prs_sci(k,i,j) = t1_qs_qi*rhof(k, i, j)*ef_si*ri(k, i, j)*smoe(k, i, j)
                tend%pni_sci(k,i,j) = tend%prs_sci(k,i,j) * oxmi
              endif

              if (rr(k, i, j) >= r_r(1) .and. mvd_r(k, i, j) > 4._wp*xdi) then
                lamr = 1._wp/ilamr(k, i, j)
                n0_r = nr(k, i, j)*org2*lamr**cre(2)
                tend%pri_rci(k,i,j) = rhof(k, i, j)*t1_qr_qi*ef_ri*ri(k, i, j)*n0_r * &
                  ((lamr+fv_r)**(-cre(9)))
                tend%pnr_rci(k,i,j) = rhof(k, i, j)*t1_qr_qi*ef_ri*ni(k, i, j)*n0_r * &
                  ((lamr+fv_r)**(-cre(9)))
                tend%pnr_rci(k,i,j) = min(real(nr(k, i, j)*odt, kind=dp), tend%pnr_rci(k,i,j))
                tend%png_rci(k,i,j) = tend%pnr_rci(k,i,j) * &
                  max(min((10._wp**(-0.1*w1d(k, i, j)) + 0.1_wp), 1._wp), 0.1_wp)
                tend%pni_rci(k,i,j) = tend%pri_rci(k,i,j) * oxmi
                tend%prr_rci(k,i,j) = rhof(k, i, j)*t2_qr_qi*ef_ri*ni(k, i, j)*n0_r * &
                  ((lamr+fv_r)**(-cre(8)))
                tend%prr_rci(k,i,j) = min(real(rr(k, i, j)*odt, kind=dp), tend%prr_rci(k,i,j))
                tend%prg_rci(k,i,j) = tend%pri_rci(k,i,j) + tend%prr_rci(k,i,j)
                tend%pbg_rci(k,i,j) = tend%prg_rci(k,i,j)/rho_i
              endif
            endif

            if (l_qs(k, i, j)) then
              c_snow = c_sqrd + (temp(k, i, j)-t0+1.5_wp)*(c_cube-c_sqrd)/(-30._wp+1.5_wp)
              c_snow = max(c_sqrd, min(c_snow, c_cube))
              tend%prs_sde(k,i,j) = c_snow*t1_subl(k)*diffu(k, i, j)*ssati(k, i, j)*rvs * (t1_qs_sd*smo1(k, i, j) + &
                t2_qs_sd*rhof2(k, i, j)*vsc2(k, i, j)*smof(k, i, j))
              if (tend%prs_sde(k,i,j) < 0._dp) then
                tend%prs_sde(k,i,j) = max(real(-rs(k, i, j)*odt, kind=dp), &
                  tend%prs_sde(k,i,j), real(rate_max, kind=dp))
              else
                tend%prs_sde(k,i,j) = min(tend%prs_sde(k,i,j), real(rate_max, kind=dp))
              endif
            endif
            if (l_qg(k, i, j)) then
              if (ssati(k, i, j) < -eps) then
                n0_g = ng(k, i, j)*ogg2*(1._wp/ilamg(k, i, j))**cge(2,1)
                t2_qg_sd = 0.28_wp*sc3*sqrt(av_g(idx(k, i, j))) * cgg(11,idx(k, i, j))
                tend%prg_gde(k,i,j) = c_cube*t1_subl(k)*diffu(k, i, j)*ssati(k, i, j)*rvs &
                    * n0_g * (t1_qg_sd*ilamg(k, i, j)**cge(10,1) &
                    + t2_qg_sd*vsc2(k, i, j)*rhof2(k, i, j)*ilamg(k, i, j)**cge(11,idx(k, i, j)))
                if (tend%prg_gde(k,i,j) < 0._wp) then
                    tend%prg_gde(k,i,j) = max(real(-rg(k, i, j)*odt, kind=dp), &
                      tend%prg_gde(k,i,j), real(rate_max, kind=dp))
                    tend%png_gde(k,i,j) = tend%prg_gde(k,i,j) * ng(k, i, j)/rg(k, i, j)
                else
                    tend%prg_gde(k,i,j) = min(tend%prg_gde(k,i,j), real(rate_max, kind=dp))
                endif
              endif
            endif
          endif
        enddo
        endif
      enddo
    enddo
    !$acc end parallel
  end subroutine ice_processes


  subroutine melting(kts, kte, its, ite, jts, jte, rhof2, rho, temp, qvsi, tcond, diffu, vsc2, ssati, &
    delqvs, l_qs, rs, smof, smo0, smo1, l_qg, rg, ng, ilamg, idx, tend, dt, odt, column_mp_active)
    !! melting of snow and graupel (horizontal tile)
    use module_mp_tempo_params, only : t0, bm_i, mu_i, pi, c_sqrd, c_cube, d0s, &
      ntb_i, r_s, ef_si, r_r, fv_r, ef_ri, rho_i, t1_qs_sd, t2_qs_sd, eps, &
      t1_qg_sd, sc3, ogg2, cge, cgg, av_g, t1_qs_me, t2_qs_me, lvap0, olfus, &
      t1_qg_me, rho_w, rho_g, meters3_to_liters, timestep_conversion_rime_to_rain

    integer, intent(in) :: kts, kte, its, ite, jts, jte
    real(wp), intent(in) :: dt, odt
    type(ty_tend), intent(inout) :: tend
    logical, dimension(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(in) :: l_qs, l_qg
    real(wp), dimension(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(in) :: rhof2, rho, rs, &
      temp, qvsi, tcond, diffu, ssati, delqvs, vsc2, rg, ng
    real(dp), dimension(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(in) :: smof, smo0, smo1, ilamg
    integer, dimension(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(in) :: idx
    logical, dimension(TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(in), optional :: column_mp_active
    real(wp) :: tempc, otemp, rvs, melt_f, t2_qg_me, t2_qg_sd
    real(dp) :: n0_g, n0_melt, lamg
    integer :: i, j, k
    real(wp), dimension(kts:kte) :: t1_subl
    logical :: active_col

    !$acc parallel
    !$acc loop gang vector collapse(2) private(active_col, t1_subl)
    do j = TEMPO_JTS, TEMPO_JTE
      do i = TEMPO_ITS, TEMPO_ITE
        active_col = .true.
        if (present(column_mp_active)) active_col = column_mp_active(i, j)
        if (active_col) then
        call get_t1_subl(kts, kte, rho(kts:kte, i, j), temp(kts:kte, i, j), qvsi(kts:kte, i, j), tcond(kts:kte, i, j), &
          diffu(kts:kte, i, j), ssati(kts:kte, i, j), t1_subl)

        !$acc loop vector private(tempc, otemp, rvs, melt_f, t2_qg_me, t2_qg_sd, n0_g, n0_melt, lamg)
        do k = kts, kte
          otemp = 1._wp/temp(k, i, j)
          tempc = temp(k, i, j) - t0
          rvs = rho(k, i, j)*qvsi(k, i, j)

          if (temp(k, i, j) > t0) then
            if(l_qs(k, i, j)) then
              tend%prr_sml(k,i,j) = (tempc*tcond(k, i, j)-lvap0*diffu(k, i, j)*delqvs(k, i, j)) * &
                (t1_qs_me*smo1(k, i, j) + t2_qs_me*rhof2(k, i, j)*vsc2(k, i, j)*smof(k, i, j))

              if (tend%prr_sml(k,i,j) > 0._dp) then
                tend%prr_sml(k,i,j) = tend%prr_sml(k,i,j) + 4218._wp*olfus*tempc * &
                  (tend%prr_rcs(k,i,j)+tend%prs_scw(k,i,j))
                tend%prr_sml(k,i,j) = min(real(rs(k, i, j)*odt, kind=dp), &
                  max(0._dp, tend%prr_sml(k,i,j)))
                tend%pnr_sml(k,i,j) = smo0(k, i, j)/rs(k, i, j)*tend%prr_sml(k,i,j) * 10.0_wp**(-0.25_wp*tempc)
                tend%pnr_sml(k,i,j) = min(real(smo0(k, i, j)*odt, kind=dp), tend%pnr_sml(k,i,j))
              else
                tend%prr_sml(k,i,j) = 0._dp
                tend%pnr_sml(k,i,j) = 0._dp
                if (ssati(k, i, j) < 0._wp) then
                  tend%prs_sde(k,i,j) = c_cube*t1_subl(k)*diffu(k, i, j)*ssati(k, i, j)*rvs * &
                    (t1_qs_sd*smo1(k, i, j) + t2_qs_sd*rhof2(k, i, j)*vsc2(k, i, j)*smof(k, i, j))
                  tend%prs_sde(k,i,j) = max(real(-rs(k, i, j)*odt, kind=dp), tend%prs_sde(k,i,j))
                endif
              endif
            endif

            if (l_qg(k, i, j)) then
              n0_g = ng(k, i, j)*ogg2*(1._wp/ilamg(k, i, j))**cge(2,1)
              n0_melt = ng(k, i, j)*ogg2*(1._dp/ilamg(k, i, j))**cge(2,1)
              if ((rg(k, i, j)*ng(k, i, j)) < 1.e-4_wp) then
                lamg = 1./ilamg(k, i, j)
                n0_melt = (1.e-4_wp/rg(k, i, j))*ogg2*lamg**cge(2,1)
              endif
              t2_qg_me = pi*4._wp * c_cube*olfus * &
                0.2_wp*sc3*sqrt(av_g(idx(k, i, j))) * cgg(11,idx(k, i, j))
              tend%prr_gml(k,i,j) = (tempc*tcond(k, i, j)-lvap0*diffu(k, i, j)*delqvs(k, i, j)) * &
                n0_melt*(t1_qg_me*ilamg(k, i, j)**cge(10,1) + &
                t2_qg_me*rhof2(k, i, j)*vsc2(k, i, j)*ilamg(k, i, j)**cge(11,idx(k, i, j)))
              tend%prr_gml(k,i,j) = min(real(rg(k, i, j)*odt, kind=dp), max(0._dp, tend%prr_gml(k,i,j)))
              if (tend%prr_gml(k,i,j) > 0._dp) then
                melt_f = max(0.05_wp, min(tend%prr_gml(k,i,j)*dt/rg(k, i, j),1._wp))
                tend%pbg_gml(k,i,j) = meters3_to_liters*tend%prr_gml(k,i,j) / &
                  max(min(melt_f*rho_g(idx(k, i, j)), rho_w), 50._wp)
                tend%pnr_gml(k,i,j) = tend%prr_gml(k,i,j)*ng(k, i, j)/rg(k, i, j) * 10.0_wp**(-0.33_wp*(temp(k, i, j)-t0))
              else
                tend%prr_gml(k,i,j) = 0._dp
                tend%pnr_gml(k,i,j) = 0._dp
                tend%pbg_gml(k,i,j) = 0._dp
                if (ssati(k, i, j) < 0._wp) then
                  t2_qg_sd = 0.28_wp*Sc3*sqrt(av_g(idx(k, i, j))) * cgg(11,idx(k, i, j))
                  tend%prg_gde(k,i,j) = C_cube*t1_subl(k)*diffu(k, i, j)*ssati(k, i, j)*rvs * n0_g * &
                    (t1_qg_sd*ilamg(k, i, j)**cge(10,1) + &
                    t2_qg_sd*vsc2(k, i, j)*rhof2(k, i, j)*ilamg(k, i, j)**cge(11,idx(k, i, j)))
                  tend%prg_gde(k,i,j) = max(real(-rg(k, i, j)*odt, kind=dp), tend%prg_gde(k,i,j))
                  tend%png_gde(k,i,j) = tend%prg_gde(k,i,j) * ng(k, i, j)/rg(k, i, j)
                endif
              endif
            endif
            if (dt > timestep_conversion_rime_to_rain) then
              tend%prr_rcw(k,i,j) = tend%prr_rcw(k,i,j)+tend%prs_scw(k,i,j)+tend%prg_gcw(k,i,j)
              tend%prs_scw(k,i,j) = 0._dp
              tend%prg_gcw(k,i,j) = 0._dp
            endif
          endif
        enddo
        endif
      enddo
    enddo
    !$acc end parallel
  end subroutine melting


  subroutine aerosol_scavenging(kts, kte, its, ite, jts, jte, temp, rho, rhof, visco, nwfa, nifa, &
    l_qr, nr, ilamr, mvd_r, l_qs, rs, smob, smoc, smoe, &
    l_qg, rg, ng, ilamg, idx, tend, odt, column_mp_active)
    !! scavenging of aerosols by rain, snow, and graupel (horizontal tile)
    use module_mp_tempo_params, only : d0r, t1_qr_qc, fv_r, cre, &
      org2, r_s, t1_qs_qc, r_g, bm_g, mu_g, av_g, cge, cgg, pi, ogg2

    integer, intent(in) :: kts, kte, its, ite, jts, jte
    real(wp), intent(in) :: odt
    real(wp), dimension(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(in) :: temp, rho, rhof, visco, nr, mvd_r, &
      nwfa, nifa, rs, rg, ng
    real(dp), dimension(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(in) :: smob, smoc, smoe, ilamg, ilamr
    integer, dimension(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(in) :: idx
    logical, dimension(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(in) :: l_qr, l_qs, l_qg
    type(ty_tend), intent(inout) :: tend
    logical, dimension(TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(in), optional :: column_mp_active
    real(wp) :: ef_ra, ef_sa, ef_ga, t1_qg_qc
    real(dp) :: n0_r, xds, xdg, n0_g, lamr
    real(wp), parameter :: wf_aerosol_size = 0.04e-6_wp
    real(wp), parameter :: if_aerosol_size = 0.8e-6_wp
    integer :: i, j, k
    logical :: active_col

    !$acc parallel
    !$acc loop gang vector collapse(2) private(active_col)
    do j = TEMPO_JTS, TEMPO_JTE
      do i = TEMPO_ITS, TEMPO_ITE
        active_col = .true.
        if (present(column_mp_active)) active_col = column_mp_active(i, j)
        if (active_col) then
        !$acc loop vector private(ef_ra, ef_sa, ef_ga, t1_qg_qc, n0_r, xds, xdg, n0_g, lamr)
        do k = kts, kte
          if (l_qr(k, i, j) .and. mvd_r(k, i, j).gt. d0r) then
            ef_ra = aerosol_collection_efficiency(real(mvd_r(k, i, j), kind=dp), &
              wf_aerosol_size, visco(k, i, j), rho(k, i, j), temp(k, i, j), 'r')
            lamr = 1._dp/ilamr(k, i, j)
            n0_r = nr(k, i, j)*org2*lamr**cre(2)
            tend%pna_rca(k,i,j) = rhof(k, i, j)*t1_qr_qc*ef_ra*nwfa(k, i, j)*n0_r * &
              ((lamr+fv_r)**(-cre(9)))
            tend%pna_rca(k,i,j) = min(real(nwfa(k, i, j)*odt, kind=dp), &
              tend%pna_rca(k,i,j))
            ef_ra = aerosol_collection_efficiency(real(mvd_r(k, i, j), kind=dp), &
              if_aerosol_size, visco(k, i, j), rho(k, i, j), temp(k, i, j), 'r')
            tend%pnd_rcd(k,i,j) = rhof(k, i, j)*t1_qr_qc*ef_ra*nifa(k, i, j)*n0_r * &
              ((lamr+fv_r)**(-cre(9)))
            tend%pnd_rcd(k,i,j) = min(real(nifa(k, i, j)*odt, kind=dp), &
              tend%pnd_rcd(k,i,j))
          endif

          if (l_qs(k, i, j) .and. rs(k, i, j) > r_s(1)) then
            xds = smoc(k, i, j) / smob(k, i, j)
            ef_sa = aerosol_collection_efficiency(xds,wf_aerosol_size, &
              visco(k, i, j), rho(k, i, j), temp(k, i, j), 's')
            tend%pna_sca(k,i,j) = rhof(k, i, j)*t1_qs_qc*ef_sa*nwfa(k, i, j)*smoe(k, i, j)
            tend%pna_sca(k,i,j) = min(real(nwfa(k, i, j)*odt, kind=dp), &
              tend%pna_sca(k,i,j))
            ef_sa = aerosol_collection_efficiency(xds, if_aerosol_size, &
              visco(k, i, j), rho(k, i, j), temp(k, i, j), 's')
            tend%pnd_scd(k,i,j) = rhof(k, i, j)*t1_qs_qc*ef_sa*nifa(k, i, j)*smoe(k, i, j)
            tend%pnd_scd(k,i,j) = min(real(nifa(k, i, j)*odt, kind=dp), &
              tend%pnd_scd(k,i,j))
          endif

          if (l_qg(k, i, j) .and. rg(k, i, j) > r_g(1)) then
            xdg = (bm_g + mu_g + 1._dp) * ilamg(k, i, j)
            ef_ga = aerosol_collection_efficiency(xdg, wf_aerosol_size, &
              visco(k, i, j), rho(k, i, j), temp(k, i, j), 'g')
            t1_qg_qc = pi*.25_wp*av_g(idx(k, i, j)) * cgg(9,idx(k, i, j))
            n0_g = ng(k, i, j)*ogg2*(1._wp/ilamg(k, i, j))**cge(2,1)
            tend%pna_gca(k,i,j) = rhof(k, i, j)*t1_qg_qc*ef_ga*nwfa(k, i, j)*n0_g * &
              ilamg(k, i, j)**cge(9,idx(k, i, j))
            tend%pna_gca(k,i,j) = min(real(nwfa(k, i, j)*odt, kind=dp), &
              tend%pna_gca(k,i,j))
            ef_ga = aerosol_collection_efficiency(xdg, if_aerosol_size, &
              visco(k, i, j), rho(k, i, j), temp(k, i, j), 'g')
            tend%pnd_gcd(k,i,j) = rhof(k, i, j)*t1_qg_qc*ef_ga*nifa(k, i, j)*n0_g * &
              ilamg(k, i, j)**cge(9,idx(k, i, j))
            tend%pnd_gcd(k,i,j) = min(real(nifa(k, i, j)*odt, kind=dp), &
              tend%pnd_gcd(k,i,j))
          endif
        enddo
        endif
      enddo
    enddo
    !$acc end parallel
  end subroutine aerosol_scavenging

end module module_mp_tempo_main
