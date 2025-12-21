module module_mp_tempo_utils
  !! utilities for tempo microphysics
  use module_mp_tempo_params !, only : wp, sp, dp
  implicit none
  private
  
  public :: snow_moments, calc_gamma_p, get_nuc, get_cloud_number, &
    calc_rslf, calc_rsif

  contains 

  function get_cloud_number(land) result(num)
    !! returns land-specific value of cloud droplet number concentration
    !! if land == 1, else returns ocean-specific value
    use module_mp_tempo_params, only : nt_c_l, nt_c_o
    
    integer, intent(in), optional :: land
    real(wp) :: num

    num = nt_c_l
    if (present(land)) then
      if (land /= 1) num = nt_c_o
    endif 
  end function


  function calc_gamma_p(a, x) result(gamma_p)
    !! normalized lower gamma function calculated either with a 
    !! series expansion or continued fraction method
    !! input: a = gamma function argument, x = upper limit of integration
    !! output: gamma_p = \(\gamma(a, x) / \Gamma(a)\)
    real(wp), intent(in) :: a, x
    real(wp) :: gamma_p

    if ((x < 0.0_wp) .or. (a <= 0.0_wp)) stop "Invalid arguments for function gamma_p"
    if (x < (a+1.0_wp)) then
      gamma_p = calc_gamma_series(a, x)
    else
      ! gammma_cf computes the upper series
      gamma_p = 1.0_wp - calc_gamma_cf(a, x)
    endif
  end function calc_gamma_p


  function calc_gamma_series(a, x) result(gamma_series)
    !! Solves the normalized lower gamma function: gamma(a,x) / Gamma(a) = x**a * gamma*(a,x)
    ! https://dlmf.nist.gov/8.7 (Equation 8.7.1)
    ! Where gamma*(a,x) = exp(-x) * sum (x**k / Gamma(a+k+1))
    ! Input:
    !   a = gamma function argument, x = upper limit of integration
    ! Output:
    !   Normalized LOWER gamma function: gamma(a, x) / Gamma(a)
    real(wp), intent(in) :: a, x
    integer :: k
    integer, parameter :: it_max = 100
    real(wp), parameter :: smallvalue = 1.e-7_wp
    real(wp) :: ap1, sum_term, sum
    real(wp) :: gamma_series

    if (x <= 0.0_wp) stop "Invalid arguments for function gamma_series"
    ! k = 0 summation term is 1 / Gamma(a+1)
    ap1 = a
    sum_term = 1.0_wp / gamma(ap1+1.0_wp)
    sum = sum_term
    do k = 1, it_max
      ap1 = ap1 + 1.0_wp
      sum_term = sum_term * x / ap1
      sum = sum + sum_term
      if (abs(sum_term) < (abs(sum) * smallvalue)) exit
    enddo
    if (k == it_max) stop "gamma_series solution did not converge"
    gamma_series = sum * x**a * exp(-x)
  end function calc_gamma_series
    

  function calc_gamma_cf(a, x) result(gamma_cf)
  !! Solves the normalized upper gamma function: gamma(a,x) / Gamma(a)
  !! using a continued fractions method (modifed Lentz Algorithm)
  ! http://functions.wolfram.com/06.06.10.0003.01
  ! Input:
  !   a = gamma function argument, x = lower limit of integration
  ! Output:
  !   Normalized UPPER gamma function: gamma(a, x) / Gamma(a)
    real(wp), intent(in) :: a, x
    integer :: k
    integer, parameter :: it_max = 100
    real(wp), parameter :: smallvalue = 1.e-7_wp
    real(wp), parameter :: offset = 1.e-30_wp
    real(wp) :: b, d, h0, c, delta, h, aj
    real(wp) :: gamma_cf
  
    b = 1.0_wp - a + x
    d = 1.0_wp / b
    h0 = offset
    c = b + (1.0_wp/offset)
    delta = c * d
    h = h0 * delta

    do k = 1, it_max
      aj = k * (a-k)
      b = b + 2.0_wp
      d = b + aj*d
      if(abs(d) < offset) d = offset
      c = b + aj/c
      if(abs(c) < offset) c = offset
      d = 1.0_wp / d
      delta = c * d
      h = h * delta
      if (abs(delta-1.0_wp) < smallvalue) exit
    enddo
    if (k == it_max) stop "gamma_cf solution did not converge"
    gamma_cf = exp(-x+a*log(x)) * h / gamma(a)
  end function calc_gamma_cf   


  subroutine snow_moments(rs, tc, smob, smoc, ns, smo0, smo1, smo2, smoe, smof, smog)
    !! Compute snow moments
    ! smo0 = 0th moment
    ! smo1 = 1st moment
    ! smo2 = 2nd moment
    ! ns   = total number concentration (not smo0)
    ! smob = rs*oams (2nd moment when bm_s=2)
    ! smoc = (bm_s+1)th moment
    ! smoe = (bv_s+2)th moment
    ! smof = (1+(bv_s+1)/2)th moment
    ! smog = (bm_s+bv_s+2)th moment
    use module_mp_tempo_params, only : bm_s, sa, sb, &
      oams, cse, csg, lam0, lam1, kap0, kap1, mu_s
    
    real(wp), intent(in) :: rs, tc
    real(dp) :: loga_, a_, b_, smo2_, m0, mrat, slam1, slam2
    real(dp), intent(out) :: smob, smoc
    real(dp), intent(out), optional :: ns, smo0, smo1, smo2, smoe, smof, smog

    ! Second moment and smob
    smob = real(rs*oams, kind=dp)
    if (bm_s > 2.0_wp-1.e-3_wp .and. bm_s < 2.0_wp+1.e-3_wp) then
      loga_ = sa(1) + sa(2)*tc + sa(3)*bm_s &
        + sa(4)*tc*bm_s + sa(5)*tc*tc &
        + sa(6)*bm_s*bm_s + sa(7)*tc*tc*bm_s &
        + sa(8)*tc*bm_s*bm_s + sa(9)*tc*tc*tc &
        + sa(10)*bm_s*bm_s*bm_s
      a_ = 10.0_wp**loga_
      b_ = sb(1) + sb(2)*tc + sb(3)*bm_s &
        + sb(4)*tc*bm_s + sb(5)*tc*tc &
        + sb(6)*bm_s*bm_s + sb(7)*tc*tc*bm_s &
        + sb(8)*tc*bm_s*bm_s + sb(9)*tc*tc*tc &
        + sb(10)*bm_s*bm_s*bm_s
      smo2_ = (smob/a_)**(1._wp/b_)
    else
      smo2_ = smob
    endif
    if (present(smo2)) smo2 = smo2_

    ! bm+1 moment. Useful for diameter calcs.
    loga_ = sa(1) + sa(2)*tc + sa(3)*cse(1) &
      + sa(4)*tc*cse(1) + sa(5)*tc*tc &
      + sa(6)*cse(1)*cse(1) + sa(7)*tc*tc*cse(1) &
      + sa(8)*tc*cse(1)*cse(1) + sa(9)*tc*tc*tc &
      + sa(10)*cse(1)*cse(1)*cse(1)
    a_ = 10.0_dp**loga_
    b_ = sb(1)+sb(2)*tc+sb(3)*cse(1) + sb(4)*tc*cse(1) &
      + sb(5)*tc*tc + sb(6)*cse(1)*cse(1) &
      + sb(7)*tc*tc*cse(1) + sb(8)*tc*cse(1)*cse(1) &
      + sb(9)*tc*tc*tc+sb(10)*cse(1)*cse(1)*cse(1)
    smoc = a_ * smo2_**b_

    ! 0th moment. Represents snow number concentration.
    loga_ = sa(1) + sa(2)*tc + sa(5)*tc*tc + sa(9)*tc*tc*tc
    a_ = 10.0**loga_
    b_ = sb(1) + sb(2)*tc + sb(5)*tc*tc + sb(9)*tc*tc*tc
    if (present(smo0)) smo0 = a_ * smo2_**b_

    ! 1st moment. Useful for depositional growth and melting.
    loga_ = sa(1) + sa(2)*tc + sa(3) &
      + sa(4)*tc + sa(5)*tc*tc &
      + sa(6) + sa(7)*tc*tc &
      + sa(8)*tc + sa(9)*tc*tc*tc &
      + sa(10)
    a_ = 10.0**loga_
    b_ = sb(1)+ sb(2)*tc + sb(3) + sb(4)*tc &
      + sb(5)*tc*tc + sb(6) &
      + sb(7)*tc*tc + sb(8)*tc &
      + sb(9)*tc*tc*tc + sb(10)
    if (present(smo1)) smo1 = a_ * smo2_**b_

    ! snow number concentration (explicit integral, not smo0)
    m0 = smob/smoc
    mrat = smob*m0*m0*m0
    slam1 = m0 * lam0
    slam2 = m0 * lam1
    if (present(ns)) then
      ns = mrat*kap0/slam1 + mrat*kap1*m0**mu_s*csg(15)/slam2**cse(15)
    endif 

    ! bv_s+2 (th) moment. Useful for riming.
    loga_ = sa(1) + sa(2)*tc + sa(3)*cse(13) &
      + sa(4)*tc*cse(13) + sa(5)*tc*tc &
      + sa(6)*cse(13)*cse(13) + sa(7)*tc*tc*cse(13) &
      + sa(8)*tc*cse(13)*cse(13) + sa(9)*tc*tc*tc &
      + sa(10)*cse(13)*cse(13)*cse(13)
    a_ = 10.0**loga_
    b_ = sb(1)+ sb(2)*tc + sb(3)*cse(13) + sb(4)*tc*cse(13) &
      + sb(5)*tc*tc + sb(6)*cse(13)*cse(13) &
      + sb(7)*tc*tc*cse(13) + sb(8)*tc*cse(13)*cse(13) &
      + sb(9)*tc*tc*tc + sb(10)*cse(13)*cse(13)*cse(13)
    if (present(smoe)) smoe = a_ * smo2_**b_

    ! 1+(bv_s+1)/2 (th) moment.  Useful for depositional growth.
    loga_ = sa(1) + sa(2)*tc + sa(3)*cse(16) &
      + sa(4)*tc*cse(16) + sa(5)*tc*tc &
      + sa(6)*cse(16)*cse(16) + sa(7)*tc*tc*cse(16) &
      + sa(8)*tc*cse(16)*cse(16) + sa(9)*tc*tc*tc &
      + sa(10)*cse(16)*cse(16)*cse(16)
    a_ = 10.0**loga_
    b_ = sb(1)+ sb(2)*tc + sb(3)*cse(16) + sb(4)*tc*cse(16) &
      + sb(5)*tc*tc + sb(6)*cse(16)*cse(16) &
      + sb(7)*tc*tc*cse(16) + sb(8)*tc*cse(16)*cse(16) &
      + sb(9)*tc*tc*tc + sb(10)*cse(16)*cse(16)*cse(16)
    if (present(smof)) smof = a_ * smo2_**b_

    ! bm_s + bv_s+2 (th) moment.  Useful for riming into graupel.
    loga_ = sa(1) + sa(2)*tc + sa(3)*cse(17) &
        + sa(4)*tc*cse(17) + sa(5)*tc*tc &
        + sa(6)*cse(17)*cse(17) + sa(7)*tc*tc*cse(17) &
        + sa(8)*tc*cse(17)*cse(17) + sa(9)*tc*tc*tc &
        + sa(10)*cse(17)*cse(17)*cse(17)
    a_ = 10.0**loga_
    b_ = sb(1)+ sb(2)*tc + sb(3)*cse(17) + sb(4)*tc*cse(17) &
        + sb(5)*tc*tc + sb(6)*cse(17)*cse(17) &
        + sb(7)*tc*tc*cse(17) + sb(8)*tc*cse(17)*cse(17) &
        + sb(9)*tc*tc*tc + sb(10)*cse(17)*cse(17)*cse(17)
    if (present(smog)) smog = a_ * smo2_**b_
  end subroutine snow_moments


 
    !=================================================================================================================

    elemental subroutine make_hydrometeor_number_concentrations(qc, qr, qi, nwfa, temp, rhoa, nc, nr, ni)
        implicit none

        real, intent(in) :: qc, qr, qi, nwfa, temp, rhoa
        real, intent(inout) :: nc, nr, ni

    end subroutine make_hydrometeor_number_concentrations

    !=================================================================================================================

    !>\ingroup aathompson
    !!Table of lookup values of radiative effective radius of ice crystals
    !! as a function of Temperature from -94C to 0C.  Taken from WRF RRTMG
    !! radiation code where it is attributed to Jon Egill Kristjansson
    !! and coauthors.
    elemental real function make_IceNumber (Q_ice, temp)

        !IMPLICIT NONE
        REAL, PARAMETER:: Ice_density = 890.0
        REAL, PARAMETER:: PI = 3.1415926536
        real, intent(in):: Q_ice, temp
        integer idx_rei
        real corr, reice, deice
        double precision lambda

        !+---+-----------------------------------------------------------------+
        !..Table of lookup values of radiative effective radius of ice crystals
        !.. as a function of Temperature from -94C to 0C.  Taken from WRF RRTMG
        !.. radiation code where it is attributed to Jon Egill Kristjansson
        !.. and coauthors.
        !+---+-----------------------------------------------------------------+

        !real retab(95)
        !data retab /                                                      &
        !   5.92779, 6.26422, 6.61973, 6.99539, 7.39234,                   &
        !   7.81177, 8.25496, 8.72323, 9.21800, 9.74075, 10.2930,          &
        !   10.8765, 11.4929, 12.1440, 12.8317, 13.5581, 14.2319,          &
        !   15.0351, 15.8799, 16.7674, 17.6986, 18.6744, 19.6955,          &
        !   20.7623, 21.8757, 23.0364, 24.2452, 25.5034, 26.8125,          &
        !   27.7895, 28.6450, 29.4167, 30.1088, 30.7306, 31.2943,          &
        !   31.8151, 32.3077, 32.7870, 33.2657, 33.7540, 34.2601,          &
        !   34.7892, 35.3442, 35.9255, 36.5316, 37.1602, 37.8078,          &
        !   38.4720, 39.1508, 39.8442, 40.5552, 41.2912, 42.0635,          &
        !   42.8876, 43.7863, 44.7853, 45.9170, 47.2165, 48.7221,          &
        !   50.4710, 52.4980, 54.8315, 57.4898, 60.4785, 63.7898,          &
        !   65.5604, 71.2885, 75.4113, 79.7368, 84.2351, 88.8833,          &
        !   93.6658, 98.5739, 103.603, 108.752, 114.025, 119.424,          &
        !   124.954, 130.630, 136.457, 142.446, 148.608, 154.956,          &
        !   161.503, 168.262, 175.248, 182.473, 189.952, 197.699,          &
        !   205.728, 214.055, 222.694, 231.661, 240.971, 250.639/
        real, dimension(95), parameter:: retab = (/                       &
            5.92779, 6.26422, 6.61973, 6.99539, 7.39234,                   &
            7.81177, 8.25496, 8.72323, 9.21800, 9.74075, 10.2930,          &
            10.8765, 11.4929, 12.1440, 12.8317, 13.5581, 14.2319,          &
            15.0351, 15.8799, 16.7674, 17.6986, 18.6744, 19.6955,          &
            20.7623, 21.8757, 23.0364, 24.2452, 25.5034, 26.8125,          &
            27.7895, 28.6450, 29.4167, 30.1088, 30.7306, 31.2943,          &
            31.8151, 32.3077, 32.7870, 33.2657, 33.7540, 34.2601,          &
            34.7892, 35.3442, 35.9255, 36.5316, 37.1602, 37.8078,          &
            38.4720, 39.1508, 39.8442, 40.5552, 41.2912, 42.0635,          &
            42.8876, 43.7863, 44.7853, 45.9170, 47.2165, 48.7221,          &
            50.4710, 52.4980, 54.8315, 57.4898, 60.4785, 63.7898,          &
            65.5604, 71.2885, 75.4113, 79.7368, 84.2351, 88.8833,          &
            93.6658, 98.5739, 103.603, 108.752, 114.025, 119.424,          &
            124.954, 130.630, 136.457, 142.446, 148.608, 154.956,          &
            161.503, 168.262, 175.248, 182.473, 189.952, 197.699,          &
            205.728, 214.055, 222.694, 231.661, 240.971, 250.639 /)

        if (Q_ice == 0) then
            make_IceNumber = 0
            return
        end if

        !+---+-----------------------------------------------------------------+
        !..From the model 3D temperature field, subtract 179K for which
        !.. index value of retab as a start.  Value of corr is for
        !.. interpolating between neighboring values in the table.
        !+---+-----------------------------------------------------------------+

        idx_rei = int(temp-179.)
        idx_rei = min(max(idx_rei,1),94)
        corr = temp - int(temp)
        reice = retab(idx_rei)*(1.-corr) + retab(idx_rei+1)*corr
        deice = 2.*reice * 1.E-6

        !+---+-----------------------------------------------------------------+
        !..Now we have the final radiative effective size of ice (as function
        !.. of temperature only).  This size represents 3rd moment divided by
        !.. second moment of the ice size distribution, so we can compute a
        !.. number concentration from the mean size and mass mixing ratio.
        !.. The mean (radiative effective) diameter is 3./Slope for an inverse
        !.. exponential size distribution.  So, starting with slope, work
        !.. backwords to get number concentration.
        !+---+-----------------------------------------------------------------+

        lambda = 3.0 / deice
        make_IceNumber = Q_ice * lambda*lambda*lambda / (PI*Ice_density)

        !+---+-----------------------------------------------------------------+
        !..Example1: Common ice size coming from Thompson scheme is about 30 microns.
        !.. An example ice mixing ratio could be 0.001 g/kg for a temperature of -50C.
        !.. Remember to convert both into MKS units.  This gives N_ice=357652 per kg.
        !..Example2: Lower in atmosphere at T=-10C matching ~162 microns in retab,
        !.. and assuming we have 0.1 g/kg mixing ratio, then N_ice=28122 per kg,
        !.. which is 28 crystals per liter of air if the air density is 1.0.
        !+---+-----------------------------------------------------------------+

        return
    end function make_IceNumber

    !=================================================================================================================
    elemental real function make_DropletNumber (Q_cloud, qnwfa)

        !IMPLICIT NONE

        real, intent(in):: Q_cloud, qnwfa

        real, parameter:: PI = 3.1415926536
        real, parameter:: am_r = PI*1000./6.
        real, dimension(15), parameter:: g_ratio = (/24,60,120,210,336,   &
        &                504,720,990,1320,1716,2184,2730,3360,4080,4896/)
        double precision:: lambda, qnc
        real:: q_nwfa, x1, xDc
        integer:: nu_c

        if (Q_cloud == 0) then
            make_DropletNumber = 0
            return
        end if

        !+---+

        q_nwfa = MAX(99.E6, MIN(qnwfa,5.E10))
        nu_c = MAX(2, MIN(NINT(2.5E10/q_nwfa), 15))

        x1 = MAX(1., MIN(q_nwfa*1.E-9, 10.)) - 1.
        xDc = (30. - x1*20./9.) * 1.E-6

        lambda = (4.0D0 + nu_c) / xDc
        qnc = Q_cloud / g_ratio(nu_c) * lambda*lambda*lambda / am_r
        make_DropletNumber = SNGL(qnc)

        return
    end function make_DropletNumber

    !=================================================================================================================
    elemental real function make_RainNumber (Q_rain, temp)

        IMPLICIT NONE

        real, intent(in):: Q_rain, temp
        double precision:: lambda, N0, qnr
        real, parameter:: PI = 3.1415926536
        real, parameter:: am_r = PI*1000./6.

        if (Q_rain == 0) then
            make_RainNumber = 0
            return
        end if

        !+---+-----------------------------------------------------------------+
        !.. Not thrilled with it, but set Y-intercept parameter to Marshal-Palmer value
        !.. that basically assumes melting snow becomes typical rain. However, for
        !.. -2C < T < 0C, make linear increase in exponent to attempt to keep
        !.. supercooled collision-coalescence (warm-rain) similar to drizzle rather
        !.. than bigger rain drops.  While this could also exist at T>0C, it is
        !.. more difficult to assume it directly from having mass and not number.
        !+---+-----------------------------------------------------------------+

        N0 = 8.E6

        if (temp .le. 271.15) then
            N0 = 8.E8
        elseif (temp .gt. 271.15 .and. temp.lt.273.15) then
            N0 = 8. * 10**(279.15-temp)
        endif

        lambda = SQRT(SQRT(N0*am_r*6.0/Q_rain))
        qnr = Q_rain / 6.0 * lambda*lambda*lambda / am_r
        make_RainNumber = SNGL(qnr)

        return
    end function make_RainNumber


  function calc_rslf(p, t) result(rslf)
    !! calculates liquid saturation vapor mixing ratio
    real(wp), intent(in) :: p, t
    real(wp) :: esl, x
    real(wp), parameter :: c0 = .611583699E03_wp
    real(wp), parameter :: c1 = .444606896E02_wp
    real(wp), parameter :: c2 = .143177157e01_wp
    real(wp), parameter :: c3 = .264224321e-1_wp
    real(wp), parameter :: c4 = .299291081e-3_wp
    real(wp), parameter :: c5 = .203154182e-5_wp
    real(wp), parameter :: c6 = .702620698e-8_wp
    real(wp), parameter :: c7 = .379534310e-11_wp
    real(wp), parameter :: c8 = -.321582393e-13_wp
    real(wp) :: rslf

    x = max(-80._wp, t-273.16)
    esl = c0+x*(c1+x*(c2+x*(c3+x*(c4+x*(c5+x*(c6+x*(c7+x*c8)))))))
    esl = min(esl, p*0.15) 
    !! @note
    !! Even with P=1050mb and T=55C, the saturation vapor 
    !! pressure only contributes to ~15% of the total pressure
    !! @endnote
    rslf = .622*esl/(p-esl)
  end function calc_rslf


  function calc_rsif(p, t) result(rsif)
    !! calculates liquid saturation vapor mixing ratio
    real(wp), intent(in) :: p, t
    real(wp) :: esi, x
    real(wp), parameter :: c0 = .609868993e03_wp
    real(wp), parameter :: c1 = .499320233e02_wp
    real(wp), parameter :: c2 = .184672631e01_wp
    real(wp), parameter :: c3 = .402737184e-1_wp
    real(wp), parameter :: c4 = .565392987e-3_wp
    real(wp), parameter :: c5 = .521693933e-5_wp
    real(wp), parameter :: c6 = .307839583e-7_wp
    real(wp), parameter :: c7 = .105785160e-9_wp
    real(wp), parameter :: c8 = .161444444e-12_wp
    real(wp) :: rsif

    x = max(-80._wp, t-273.16)
    esi = c0+x*(c1+x*(c2+x*(c3+x*(c4+x*(c5+x*(c6+x*(c7+x*c8)))))))
    esi = min(esi, p*0.15)
    rsif = .622*esi/max(1.e-4_wp,(p-esi))
  end function calc_rsif


  function get_nuc(nc) result(nu_c)
    !! returns nu_c for cloud water (values from 2-15)
    use module_mp_tempo_params, only : nu_c_scale

    real(wp), intent(in) :: nc
    integer :: nu_c
    
    ! write(*,*) 'in nuc', nu_c_scale, nc, nu_c_scale/nc, nint(nu_c_scale/nc), &
    !   nint(nu_c_scale/nc) + 2

    if ((nu_c_scale/nc) >= 12.5_wp) then
      nu_c = 15
    else
      nu_c = min(15, nint(nu_c_scale/nc) + 2)
    endif 
  end function get_nuc

end module module_mp_tempo_utils
