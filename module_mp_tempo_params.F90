! Parameters file for TEMPO Microphysics
!=================================================================================================================
module module_mp_tempo_params

#if defined(mpas)
    use mpas_kind_types, only: wp => RKIND, sp => R4KIND, dp => R8KIND
#elif defined(standalone)
    use machine, only: wp => kind_phys, sp => kind_sngl_prec, dp => kind_dbl_prec
#else
    use machine, only: wp => kind_phys, sp => kind_sngl_prec, dp => kind_dbl_prec
#define ccpp_default 1
#endif

#if defined(ccpp_default) && defined(MPI)
    use mpi_f08
#endif

#undef OLD_MPTBLS

    implicit none

    !=================================================================================================================
    ! Parameters needed first by tempo_init()

#if defined(OLD_MPTBLS)
    logical, parameter :: original_thompson = .true.
#else
    logical, parameter :: original_thompson = .false.
#endif

! #########################################################################################
!> \section arg_table_ty_tempo_cfg Argument Table
!! \htmlinclude ty_tempo_cfg.html
!!
! #########################################################################################
    ! Derived data type for configuration flags.
    type ty_tempo_cfg
        logical :: aerosol_aware
        logical :: hail_aware
    end type ty_tempo_cfg

    type(ty_tempo_cfg) configs

    ! Constants that can be defined by the model ===========================
    ! Needed by tempo_init()
    real(wp)            :: PI = 3.1415926536

    ! Enthalpy of sublimation, vaporization, and fusion at 0C.
    real(wp)            :: lsub = 2.834e6
    real(wp)            :: lvap0 = 2.5e6
    real(wp)            :: lfus !(set in mp_tempo_params_init)
    real(wp)            :: olfus !(set in mp_tempo_params_init)

    ! Needed by the driver
    ! Water vapor and air gas constants at constant pressure
    real(wp)            :: Rv = 461.5
    real(wp)            :: oRv !(set in mp_tempo_params_init)
    real(wp)            :: R = 287.04
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
    logical :: sedi_semi = .false.

    ! Hail-aware microphysics options
    logical, parameter :: build_hail_aware_table = .true.
    logical :: using_hail_aware_table = .false.
    integer, parameter :: NRHG = 9
    integer, parameter :: NRHG1 = 1
    integer :: dimNRHG
    integer, parameter :: idx_bg1 = 6 ! => corresponds to graupel density of 500 kg m^-3 from rho_g array

    ! Densities of rain, graupel, and cloud ice.
    real(wp), parameter :: rho_w2 = 1000.0 ! Change to rho_w2 to solve MPAS same var name conflict
    real(wp), parameter :: rho_i = 890.0
    real(wp), dimension(NRHG), parameter :: rho_g = (/50., 100., 200., 300., 400., 500., 600., 700., 800./)

    ! Cloud droplet distribution dispersion parameters
    integer, parameter :: nu_c_max = 15
    integer, parameter :: nu_c_min = 2
    real(wp), parameter :: nu_c_scale = 1000.e6

    ! Generalized gamma distributions for rain, graupel and cloud ice.
    ! N(D) = N_0 * D**mu * exp(-lamda*D);  mu=0 is exponential.
    real(wp), parameter :: mu_r = 0.0
    real(wp), parameter :: mu_s = 0.6357
    real(wp), parameter :: mu_g = 0.0
    real(wp), parameter :: mu_i = 0.0

    ! Mass power law relations:  mass = am*D**bm
    ! Snow from Field et al. (2005), others assume spherical form.
    real(wp)            :: am_r ! = PI * rho_w2 / 6.0 (set in mp_tempo_params_init)
    real(wp), parameter :: bm_r = 3.0
    real(wp), parameter :: am_s = 0.069
    real(wp), parameter :: bm_s = 2.0
    real(wp), dimension (NRHG) :: am_g  !(set in mp_tempo_params_init)
    real(wp), parameter :: bm_g = 3.0
    real(wp)            :: am_i ! = PI * rho_i / 6.0 (set in mp_tempo_params_init)
    real(wp), parameter :: bm_i = 3.0

    ! Fallspeed power laws relations:  v = (av*D**bv)*exp(-fv*D)
    ! Rain from Ferrier (1994), ice, snow, and graupel from
    ! Thompson et al (2008). Coefficient fv is zero for graupel/ice.
    real(wp), parameter :: av_r = 4854.0
    real(wp), parameter :: bv_r = 1.0
    real(wp), parameter :: av_s = 40.0
    real(wp), parameter :: bv_s = 0.55
    real(wp), parameter :: fv_s = 100.0
    real(wp), parameter :: av_g_old = 442.0
    real(wp), parameter :: bv_g_old = 0.89

    ! Variable graupel density av_g and bv_g values
    ! Computed from A. Heymsfield: Best - Reynolds relationship
    real(wp), dimension(NRHG) :: av_g = (/45.9173813, 67.0867386, 98.0158463, 122.353378, &
        143.204224, 161.794724, 178.762115, 194.488785, 209.225876/)
    real(wp), dimension(NRHG) :: bv_g = (/0.640961647, 0.640961647, 0.640961647, 0.640961647, &
        0.640961647, 0.640961647, 0.640961647, 0.640961647, 0.640961647/)

    ! Capacitance of sphere and plates/aggregates: D**3, D**2
    real(wp), parameter :: C_cube = 0.5
    real(wp), parameter :: C_sqrd = 0.15

    ! Schmidt number
    real(wp), parameter :: Sc = 0.632
    real(wp) :: Sc3

    ! Fall speed coefficients for ice and cloud droplets
    real(wp), parameter :: bv_i = 1.0
    real(wp), parameter :: bv_c = 2.0


    ! Ice initiates with this mass (kg), corresponding diameter calc.
    ! Min diameters and mass of cloud, rain, snow, and graupel (m, kg).
    real(wp), parameter :: xm0i = 1.e-12
    real(wp), parameter :: D0c = 1.e-6
    real(wp), parameter :: D0r = 50.e-6
    real(wp), parameter :: D0s = 300.e-6
    real(wp), parameter :: D0g = 350.e-6
    real(wp) :: D0i, xm0s, xm0g

    ! Lookup table dimensions
    integer, parameter :: nbins = 100
    integer, parameter :: nbc = nbins
    integer, parameter :: nbr = nbins
    integer, parameter :: nbs = nbins

    integer, parameter :: nbi = nbins
    integer, parameter :: nbg = nbins
    integer, parameter :: ntb_i = 64
    integer, parameter :: ntb_i1 = 55
    integer, parameter :: ntb_c = 37
    integer, parameter :: ntb_t = 9
    integer, parameter :: ntb_g1 = 37

#if defined(ccpp_default) && defined(OLD_MPTBLS)
    integer, parameter :: ntb_s = 28
    integer, parameter :: ntb_g = 28
#else
    integer, parameter :: ntb_s = 37
    integer, parameter :: ntb_g = 37
#endif
    integer, parameter :: ntb_r = 37
    integer, parameter :: ntb_r1 = 37
    integer, parameter :: ntb_t1 = 45
    integer, parameter :: ntb_IN = 55
    integer, parameter :: ntb_arc = 7
    integer, parameter :: ntb_arw = 9
    integer, parameter :: ntb_art = 7
    integer, parameter :: ntb_arr = 5
    integer, parameter :: ntb_ark = 4

    integer :: nic1, nic2, nii2, nii3, nir2, nir3, nis2, nig2, nig3, niIN2

    !     real(dp), dimension(nbins+1) :: xDx
    real(dp), dimension(nbc) :: Dc, dtc
    real(dp), dimension(nbi) :: Di, dti
    real(dp), dimension(nbr) :: Dr, dtr
    real(dp), dimension(nbs) :: Ds, dts
    real(dp), dimension(nbg) :: Dg, dtg
    real(dp), dimension(nbc) :: t_Nc

    ! Lookup tables for cloud water content (kg/m**3).
    real(wp), dimension(ntb_c), parameter :: &
        r_c = (/1.e-6,2.e-6,3.e-6,4.e-6,5.e-6,6.e-6,7.e-6,8.e-6,9.e-6, &
        1.e-5,2.e-5,3.e-5,4.e-5,5.e-5,6.e-5,7.e-5,8.e-5,9.e-5, &
        1.e-4,2.e-4,3.e-4,4.e-4,5.e-4,6.e-4,7.e-4,8.e-4,9.e-4, &
        1.e-3,2.e-3,3.e-3,4.e-3,5.e-3,6.e-3,7.e-3,8.e-3,9.e-3, &
        1.e-2/)

    ! Lookup tables for cloud ice content (kg/m**3).
    real(wp), dimension(ntb_i), parameter :: &
        r_i = (/1.e-10,2.e-10,3.e-10,4.e-10, &
        5.e-10,6.e-10,7.e-10,8.e-10,9.e-10, &
        1.e-9,2.e-9,3.e-9,4.e-9,5.e-9,6.e-9,7.e-9,8.e-9,9.e-9, &
        1.e-8,2.e-8,3.e-8,4.e-8,5.e-8,6.e-8,7.e-8,8.e-8,9.e-8, &
        1.e-7,2.e-7,3.e-7,4.e-7,5.e-7,6.e-7,7.e-7,8.e-7,9.e-7, &
        1.e-6,2.e-6,3.e-6,4.e-6,5.e-6,6.e-6,7.e-6,8.e-6,9.e-6, &
        1.e-5,2.e-5,3.e-5,4.e-5,5.e-5,6.e-5,7.e-5,8.e-5,9.e-5, &
        1.e-4,2.e-4,3.e-4,4.e-4,5.e-4,6.e-4,7.e-4,8.e-4,9.e-4, &
        1.e-3/)

    ! Lookup tables for rain content (kg/m**3).
    real(wp), dimension(ntb_r), parameter :: &
        r_r = (/1.e-6,2.e-6,3.e-6,4.e-6,5.e-6,6.e-6,7.e-6,8.e-6,9.e-6, &
        1.e-5,2.e-5,3.e-5,4.e-5,5.e-5,6.e-5,7.e-5,8.e-5,9.e-5, &
        1.e-4,2.e-4,3.e-4,4.e-4,5.e-4,6.e-4,7.e-4,8.e-4,9.e-4, &
        1.e-3,2.e-3,3.e-3,4.e-3,5.e-3,6.e-3,7.e-3,8.e-3,9.e-3, &
        1.e-2/)

    ! Lookup tables for graupel y-intercept parameter (/m**4).
    real(wp), dimension(ntb_g1), parameter :: &
        N0g_exp = (/1.e2,2.e2,3.e2,4.e2,5.e2,6.e2,7.e2,8.e2,9.e2, &
        1.e3,2.e3,3.e3,4.e3,5.e3,6.e3,7.e3,8.e3,9.e3, &
        1.e4,2.e4,3.e4,4.e4,5.e4,6.e4,7.e4,8.e4,9.e4, &
        1.e5,2.e5,3.e5,4.e5,5.e5,6.e5,7.e5,8.e5,9.e5, &
        1.e6/)

#if defined(ccpp_default) && defined(OLD_MPTBLS)
    ! Lookup tables for graupel content (kg/m**3).
    real(wp), dimension(ntb_g), parameter :: &
        r_g = (/1.e-5,2.e-5,3.e-5,4.e-5,5.e-5,6.e-5,7.e-5,8.e-5,9.e-5, &
        1.e-4,2.e-4,3.e-4,4.e-4,5.e-4,6.e-4,7.e-4,8.e-4,9.e-4, &
        1.e-3,2.e-3,3.e-3,4.e-3,5.e-3,6.e-3,7.e-3,8.e-3,9.e-3, &
        1.e-2/)

    ! Lookup tables for snow content (kg/m**3).
    real(wp), dimension(ntb_s), parameter :: &
        r_s = (/1.e-5,2.e-5,3.e-5,4.e-5,5.e-5,6.e-5,7.e-5,8.e-5,9.e-5, &
        1.e-4,2.e-4,3.e-4,4.e-4,5.e-4,6.e-4,7.e-4,8.e-4,9.e-4, &
        1.e-3,2.e-3,3.e-3,4.e-3,5.e-3,6.e-3,7.e-3,8.e-3,9.e-3, &
        1.e-2/)
#else
    ! Lookup tables for graupel content (kg/m**3).
    real(wp), dimension(ntb_g), parameter :: &
        r_g = (/1.e-6,2.e-6,3.e-6,4.e-6,5.e-6,6.e-6,7.e-6,8.e-6,9.e-6, &
        1.e-5,2.e-5,3.e-5,4.e-5,5.e-5,6.e-5,7.e-5,8.e-5,9.e-5, &
        1.e-4,2.e-4,3.e-4,4.e-4,5.e-4,6.e-4,7.e-4,8.e-4,9.e-4, &
        1.e-3,2.e-3,3.e-3,4.e-3,5.e-3,6.e-3,7.e-3,8.e-3,9.e-3, &
        1.e-2/)

    ! Lookup tables for snow content (kg/m**3).
    real(wp), dimension(ntb_s), parameter :: &
        r_s = (/1.e-6,2.e-6,3.e-6,4.e-6,5.e-6,6.e-6,7.e-6,8.e-6,9.e-6, &
        1.e-5,2.e-5,3.e-5,4.e-5,5.e-5,6.e-5,7.e-5,8.e-5,9.e-5, &
        1.e-4,2.e-4,3.e-4,4.e-4,5.e-4,6.e-4,7.e-4,8.e-4,9.e-4, &
        1.e-3,2.e-3,3.e-3,4.e-3,5.e-3,6.e-3,7.e-3,8.e-3,9.e-3, &
        1.e-2/)
#endif

    ! Lookup tables for rain y-intercept parameter (/m**4).
    real(wp), dimension(ntb_r1), parameter :: &
        N0r_exp = (/1.e6,2.e6,3.e6,4.e6,5.e6,6.e6,7.e6,8.e6,9.e6, &
        1.e7,2.e7,3.e7,4.e7,5.e7,6.e7,7.e7,8.e7,9.e7, &
        1.e8,2.e8,3.e8,4.e8,5.e8,6.e8,7.e8,8.e8,9.e8, &
        1.e9,2.e9,3.e9,4.e9,5.e9,6.e9,7.e9,8.e9,9.e9, &
        1.e10/)

    ! Lookup tables for ice number concentration (/m**3).
    real(wp), dimension(ntb_i1), parameter :: &
        Nt_i = (/1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0, &
        1.e1,2.e1,3.e1,4.e1,5.e1,6.e1,7.e1,8.e1,9.e1, &
        1.e2,2.e2,3.e2,4.e2,5.e2,6.e2,7.e2,8.e2,9.e2, &
        1.e3,2.e3,3.e3,4.e3,5.e3,6.e3,7.e3,8.e3,9.e3, &
        1.e4,2.e4,3.e4,4.e4,5.e4,6.e4,7.e4,8.e4,9.e4, &
        1.e5,2.e5,3.e5,4.e5,5.e5,6.e5,7.e5,8.e5,9.e5, &
        1.e6/)

    ! Lookup tables for IN concentration (/m**3) from 0.001 to 1000/Liter.
    real(wp), dimension(ntb_IN), parameter :: &
        Nt_IN = (/1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0, &
        1.e1,2.e1,3.e1,4.e1,5.e1,6.e1,7.e1,8.e1,9.e1, &
        1.e2,2.e2,3.e2,4.e2,5.e2,6.e2,7.e2,8.e2,9.e2, &
        1.e3,2.e3,3.e3,4.e3,5.e3,6.e3,7.e3,8.e3,9.e3, &
        1.e4,2.e4,3.e4,4.e4,5.e4,6.e4,7.e4,8.e4,9.e4, &
        1.e5,2.e5,3.e5,4.e5,5.e5,6.e5,7.e5,8.e5,9.e5, &
        1.e6/)

    ! Lookup tables for various accretion/collection terms.
    ! ntb_x refers to the number of elements for rain, snow, graupel,
    ! and temperature array indices.  Variables beginning with t-p/c/m/n
    ! represent lookup tables.  Save compile-time memory by making
    ! allocatable (2009Jun12, J. Michalakes).
    real(dp), allocatable, dimension(:,:,:,:,:) :: tcg_racg, tmr_racg, tcr_gacr, tnr_racg, tnr_gacr
    real(dp), allocatable, dimension(:,:,:,:) :: tcs_racs1, tmr_racs1, tcs_racs2, tmr_racs2, &
        tcr_sacr1, tms_sacr1, tcr_sacr2, tms_sacr2, tnr_racs1, tnr_racs2, tnr_sacr1, tnr_sacr2
    real(dp), allocatable, dimension(:,:,:,:) :: tpi_qcfz, tni_qcfz
    real(dp), allocatable, dimension(:,:,:,:) :: tpi_qrfz, tpg_qrfz, tni_qrfz, tnr_qrfz
    real(dp), allocatable, dimension(:,:) :: tps_iaus, tni_iaus, tpi_ide
    real(dp), allocatable, dimension(:,:) :: t_Efrw, t_Efsw
    real(dp), allocatable, dimension(:,:,:) :: tnr_rev, tpc_wev, tnc_wev
    real(sp), allocatable, dimension(:,:,:,:,:) :: tnccn_act

    ! Variables holding a bunch of exponents and gamma values (cloud water,
    ! cloud ice, rain, snow, then graupel).
    real(wp), dimension(5,15) :: cce, ccg
    real(wp), dimension(15) ::  ocg1, ocg2
    real(wp), dimension(7) :: cie, cig
    real(wp) :: oig1, oig2, obmi
    real(wp), dimension(13) :: cre, crg
    real(wp) :: ore1, org1, org2, org3, obmr
    real(wp) :: oams, obms, ocms
    real(wp), dimension(12,NRHG) :: cge, cgg
    real(wp), dimension(NRHG) :: oamg, ocmg
#if defined(ccpp_default) && defined(OLD_MPTBLS)
    real, dimension(18) :: cse, csg
#else
    real(wp), dimension(17) :: cse, csg
#endif
    real(wp) :: oge1, ogg1, ogg2, ogg3, obmg

    ! Declaration of precomputed constants in various rate eqns.
    real(wp) :: t1_qr_qc, t1_qr_qi, t2_qr_qi, t1_qs_qc, t1_qs_qi
    real(wp) :: t1_qr_ev, t2_qr_ev, t1_qg_qc, t2_qg_sd, t2_qg_me
    real(wp) :: t1_qs_sd, t2_qs_sd, t1_qs_me, t2_qs_me, t1_qg_sd, t1_qg_me

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
    real(wp) ::  mu_c_o, mu_c_l

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
    type(MPI_Comm) :: mpi_communicator

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

  contains
    
    subroutine mp_tempo_params_init()
      ! Any variables defined in the module that are not constants/parameters are set here; if any of the values are overwritten (i.e. by a host model),
      ! then they must be reinitialized once the new values from the host are used
      lfus = lsub - lvap0
      olfus = 1./lfus
      
      oRv = 1.0 / Rv
      
      am_r = PI * rho_w2 / 6.0
      
      am_g = (/PI*rho_g(1)/6.0, &
          PI*rho_g(2)/6.0, &
          PI*rho_g(3)/6.0, &
          PI*rho_g(4)/6.0, &
          PI*rho_g(5)/6.0, &
          PI*rho_g(6)/6.0, &
          PI*rho_g(7)/6.0, &
          PI*rho_g(8)/6.0, &
          PI*rho_g(9)/6.0/)
      
      am_i = PI * rho_i / 6.0
      
      ma_w = M_w / N_avo
      
      ar_volume = 4.0 / 3.0 * PI * (2.5e-6)**3
      
    end subroutine mp_tempo_params_init
    
end module module_mp_tempo_params
