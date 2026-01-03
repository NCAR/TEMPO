module module_mp_tempo_utils
  !! utilities for tempo microphysics
  use module_mp_tempo_params, only : wp, sp, dp
  implicit none
  private
  
  public :: snow_moments, calc_gamma_p, get_nuc, get_constant_cloud_number, &
    calc_rslf, calc_rsif

  contains 

  subroutine get_constant_cloud_number(land, nc)
    !! returns land-specific value of cloud droplet number concentration
    !! when aerosol-aware = false if land = 1, else returns ocean-specific value
    use module_mp_tempo_params, only : nt_c_l, nt_c_o
    
    integer, intent(in), optional :: land
    real(wp), dimension(:), intent(out) :: nc

    nc = nt_c_l
    if (present(land)) then
      if (land /= 1) nc = nt_c_o
    endif 
  end subroutine get_constant_cloud_number


  function calc_gamma_p(a, x) result(gamma_p)
    !! normalized lower gamma function calculated either with a 
    !! series expansion or continued fraction method
    !!
    !! input: a = gamma function argument, x = upper limit of integration
    !!
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
    !! solves the normalized lower gamma function
    !!
    !! \(\gamma(a,x) / \Gamma(a) = x^{a} * \gamma*(a,x)\)
    !!
    !! see [Equation 8.7.1](https://dlmf.nist.gov/8.7)
    !! \(\gamma(a,x) = exp(-x) * \sum_{k=0}^{\infty} \frac{x^{k}}{\Gamma(a+k+1)}\)
    !!
    !! input: a = gamma function argument, x = upper limit of integration
    !!
    !! output: normalized lower gamma function \(\gamma(a, x) / \Gamma(a)\)
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
  !! solves the normalized upper gamma function \(\gamma(a,x) / \Gamma(a)\)
  !! using a continued fractions method
  !! [(modified Lentz Algorithm)](http://functions.wolfram.com/06.06.10.0003.01)
  !!
  !!input: a = gamma function argument, x = lower limit of integration
  !!
  !!output: normalized upper gamma function: \(\gamma(a, x) / \Gamma(a)\)
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


  subroutine snow_moments(rs, tc, smob, smoc, ns, smo0, smo1, smo2, smoe, smof, smog, smoz)
    !! computes snow moments from
    !! [Field et al. (2005)](https://doi.org/10.1256/qj.04.134)
    ! smo0 = 0th moment
    ! smo1 = 1st moment
    ! smo2 = 2nd moment
    ! ns   = total number concentration (not smo0)
    ! smob = rs*oams (2nd moment when bm_s=2)
    ! smoc = (bm_s+1)th moment
    ! smoe = (bv_s+2)th moment
    ! smof = (1+(bv_s+1)/2)th moment
    ! smog = (bm_s+bv_s+2)th moment
    ! somz = (bm**2)th moment for reflectivity
    use module_mp_tempo_params, only : bm_s, sa, sb, &
      oams, cse, csg, lam0, lam1, kap0, kap1, mu_s
    
    real(wp), intent(in) :: rs, tc
    real(dp) :: loga_, a_, b_, smo2_, m0, mrat, slam1, slam2
    real(dp), intent(out) :: smob, smoc
    real(dp), intent(out), optional :: ns, smo0, smo1, smo2, smoe, smof, smog, smoz

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

    !..Calculate bm_s*2 (th) moment.  Useful for reflectivity.
    loga_ = sa(1) + sa(2)*tc + sa(3)*cse(3) &
      + sa(4)*tc*cse(3) + sa(5)*tc*tc &
      + sa(6)*cse(3)*cse(3) + sa(7)*tc*tc*cse(3) &
      + sa(8)*tc*cse(3)*cse(3) + sa(9)*tc*tc*tc &
      + sa(10)*cse(3)*cse(3)*cse(3)
    a_ = 10.0**loga_
    b_ = sb(1)+ sb(2)*tc + sb(3)*cse(3) + sb(4)*tc*cse(3) &
      + sb(5)*tc*tc + sb(6)*cse(3)*cse(3) &
      + sb(7)*tc*tc*cse(3) + sb(8)*tc*cse(3)*cse(3) &
      + sb(9)*tc*tc*tc + sb(10)*cse(3)*cse(3)*cse(3)
    if (present(smoz)) smoz = a_ * smo2**b_

  end subroutine snow_moments


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

    x = max(-80._wp, t-273.16_wp)
    esl = c0+x*(c1+x*(c2+x*(c3+x*(c4+x*(c5+x*(c6+x*(c7+x*c8)))))))
    esl = min(esl, p*0.15) 
    !! @note
    !! even with p = 1050 mb and t = 55 C, the saturation vapor 
    !! pressure only contributes to 15% of the total pressure,
    !! thus the limit on the saturation vapor pressure is set to 15%
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

    x = max(-80._wp, t-273.16_wp)
    esi = c0+x*(c1+x*(c2+x*(c3+x*(c4+x*(c5+x*(c6+x*(c7+x*c8)))))))
    esi = min(esi, p*0.15)
    rsif = .622*esi/max(1.e-4_wp,(p-esi))
  end function calc_rsif


  function get_nuc(nc) result(nu_c)
    !! returns nu_c for cloud water (values from 2-15)
    use module_mp_tempo_params, only : nu_c_scale

    real(wp), intent(in) :: nc
    integer :: nu_c

    if ((nu_c_scale/nc) >= 12.5_wp) then
      nu_c = 15
    else
      nu_c = min(15, nint(nu_c_scale/nc) + 2)
    endif 
  end function get_nuc

end module module_mp_tempo_utils
