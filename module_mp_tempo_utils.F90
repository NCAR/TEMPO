! Utilities for TEMPO Microphysics
!=================================================================================================================
module module_mp_tempo_utils

    use module_mp_tempo_params

#if defined(mpas)
    use mpas_kind_types, only: wp => RKIND, sp => R4KIND, dp => R8KIND
    use mp_radar
#elif defined(standalone)
    use machine, only: wp => kind_phys, sp => kind_sngl_prec, dp => kind_dbl_prec
    use module_mp_radar
#else
    use machine, only: wp => kind_phys, sp => kind_sngl_prec, dp => kind_dbl_prec
    use module_mp_radar
#define ccpp_default 1
#endif

#if defined(ccpp_default) && defined(MPI)
    use mpi_f08
#endif

    implicit none

contains
    !=================================================================================================================
    ! Normalized lower gamma function calculated either with a series expansion or continued fraction method
    ! Input:
    !   a = gamma function argument, x = upper limit of integration
    ! Output:
    !   gamma_p = gamma(a, x) / Gamma(a)

    real function gamma_p(a, x)

        real(wp), intent(in) :: a, x

        if ((x < 0.0) .or. (a <= 0.0)) stop "Invalid arguments for function gamma_p"

        if (x < (a+1.0)) then
            gamma_p = gamma_series(a, x)
        else
            ! gammma_cf computes the upper series
            gamma_p = 1.0 - gamma_cf(a, x)
        endif

    end function gamma_p

    !=================================================================================================================
    ! https://dlmf.nist.gov/8.7 (Equation 8.7.1)
    ! Solves the normalized lower gamma function: gamma(a,x) / Gamma(a) = x**a * gamma*(a,x)
    ! Where gamma*(a,x) = exp(-x) * sum (x**k / Gamma(a+k+1))
    ! Input:
    !   a = gamma function argument, x = upper limit of integration
    ! Output:
    !   Normalized LOWER gamma function: gamma(a, x) / Gamma(a)

    real function gamma_series(a, x)

        real(wp), intent(in) :: a, x
        integer :: k
        integer, parameter :: it_max = 100
        real(wp), parameter :: smallvalue = 1.e-7
        real(wp) :: ap1, sum_term, sum
      
        if (x <= 0.0) stop "Invalid arguments for function gamma_series"

        ! k = 0 summation term is 1 / Gamma(a+1)
        ap1 = a
        sum_term = 1.0 / gamma(ap1+1.0)
        sum = sum_term
        do k = 1, it_max
           ap1 = ap1 + 1.0
           sum_term = sum_term * x / ap1
           sum = sum + sum_term
           if (abs(sum_term) < (abs(sum) * smallvalue)) exit
        enddo
        if (k == it_max) stop "gamma_series solution did not converge"

        gamma_series = sum * x**a * exp(-x)

    end function gamma_series

    !=================================================================================================================
    ! http://functions.wolfram.com/06.06.10.0003.01
    ! Solves the normalized upper gamma function: gamma(a,x) / Gamma(a)
    ! Using a continued fractions method (modifed Lentz Algorithm)
    ! Input:
    !   a = gamma function argument, x = lower limit of integration
    ! Output:
    !   Normalized UPPER gamma function: gamma(a, x) / Gamma(a)

    real function gamma_cf(a, x)

        real(wp), intent(in) :: a, x
        integer :: k
        integer, parameter :: it_max = 100
        real(wp), parameter :: smallvalue = 1.e-7
        real(wp), parameter :: offset = 1.e-30
        real(wp) :: b, d, h0, c, delta, h, aj

        b = 1.0 - a + x
        d = 1.0 / b
        h0 = offset
        c = b + (1.0/offset)
        delta = c * d
        h = h0 * delta

        do k = 1, it_max
            aj = k * (a-k)
            b = b + 2.0
            d = b + aj*d
            if(abs(d) < offset) d = offset
            c = b + aj/c
            if(abs(c) < offset) c = offset
            d = 1.0 / d
            delta = c * d
            h = h * delta
            if (abs(delta-1.0) < smallvalue) exit
        enddo
        if (k == it_max) stop "gamma_cf solution did not converge"

        gamma_cf = exp(-x+a*log(x)) * h / gamma(a)

    end function gamma_cf   

    !=================================================================================================================
    ! Calculates log-spaced bins for hydrometer sizes to simplify calculations later
    ! Input:
    !   numbins, lowbin, highbin
    ! Output:
    !   bins, deltabins

    subroutine create_bins(numbins, lowbin, highbin, bins, deltabins)

        integer, intent(in) :: numbins
        real(dp), intent(in) :: lowbin, highbin

        integer :: n
        real(dp), dimension(numbins+1) :: xDx
        real(dp), dimension(:), intent(out) :: bins
        real(dp), dimension(:), intent(out), optional :: deltabins

        xDx(1) = lowbin
        xDx(numbins+1) = highbin

        do  n = 2, numbins
            xDx(n) = exp(real(n-1, kind=dp)/real(numbins, kind=dp) * log(xDx(numbins+1)/xDx(1)) + log(xDx(1)))
        enddo

        do n = 1, numbins
            bins(n) = sqrt(xDx(n)*xDx(n+1))
        enddo

        if (present(deltabins)) then
            do n = 1, numbins
                deltabins(n) = xDx(n+1) - xDx(n)
            enddo
        endif

    end subroutine create_bins

    !=================================================================================================================
    ! Variable collision efficiency for rain collecting cloud water using method of Beard and Grover, 1974
    ! if a/A less than 0.25; otherwise uses polynomials to get close match of Pruppacher & Klett Fig 14-9.
    ! Output:
    !    t_Efrw

    subroutine table_Efrw

        ! Local variables
        real(dp) :: vtr, stokes, reynolds, Ef_rw
        real(dp) :: p, yc0, F, G, H, z, K0, X
        integer :: i, j

        do j = 1, nbc
            do i = 1, nbr
                Ef_rw = 0.0
                p = Dc(j) / Dr(i)
                if (Dr(i) < 50.e-6 .or. Dc(j) < 3.e-6) then
                    t_Efrw(i,j) = 0.0
                elseif (p > 0.25) then
                    X = Dc(j) * 1.e6_dp
                    if (Dr(i) < 75.e-6) then
                        Ef_rw = 0.026794*X - 0.20604
                    elseif (Dr(i) < 125.e-6) then
                        Ef_rw = -0.00066842*X*X + 0.061542*X - 0.37089
                    elseif (Dr(i) < 175.e-6) then
                        Ef_rw = 4.091e-06*X*X*X*X - 0.00030908*X*X*X + 0.0066237*X*X - 0.0013687*X - 0.073022
                    elseif (Dr(i) < 250.e-6) then
                        Ef_rw = 9.6719e-5*X*X*X - 0.0068901*X*X + 0.17305*X - 0.65988
                    elseif (Dr(i) < 350.e-6) then
                        Ef_rw = 9.0488e-5*X*X*X - 0.006585*X*X + 0.16606*X - 0.56125
                    else
                        Ef_rw = 0.00010721*X*X*X - 0.0072962*X*X + 0.1704*X - 0.46929
                    endif
                else
                    vtr = -0.1021 + 4.932e3*Dr(i) - 0.9551e6*Dr(i)*Dr(i) + 0.07934e9*Dr(i)*Dr(i)*Dr(i) &
                        - 0.002362e12*Dr(i)*Dr(i)*Dr(i)*Dr(i)
                    stokes = Dc(j) * Dc(j) * vtr * rho_w2 / (9.*1.718e-5*Dr(i))
                    reynolds = 9. * stokes / (p*p*rho_w2)

                    F = log(reynolds)
                    G = -0.1007_dp - 0.358_dp*F + 0.0261_dp*F*F
                    K0 = exp(G)
                    z = log(stokes / (K0+1.e-15_dp))
                    H = 0.1465_dp + 1.302_dp*z - 0.607_dp*z*z + 0.293_dp*z*z*z
                    yc0 = 2.0_dp / PI * atan(H)
                    Ef_rw = (yc0+p)*(yc0+p) / ((1.+p)*(1.+p))

                endif
                t_Efrw(i,j) = max(0.0, min(real(Ef_rw, kind=wp), 0.95))
            enddo
        enddo

    end subroutine table_Efrw

    !=================================================================================================================
    ! Variable collision efficiency for snow collecting cloud water using method of Wang and Ji, 2000 except
    ! equate melted snow diameter to their "effective collision cross-section."
    ! Output:
    !    t_Efsw

    subroutine table_Efsw

        ! Local variables
        real(dp) :: Ds_m, vts, vtc, stokes, reynolds, Ef_sw
        real(dp) :: p, yc0, F, G, H, z, K0
        integer :: i, j

        do j = 1, nbc
            vtc = 1.19e4_dp * (1.0e4_dp*Dc(j)*Dc(j)*0.25_dp)
            do i = 1, nbs
                vts = av_s*Ds(i)**bv_s * exp(-fv_s*Ds(i)) - vtc
                Ds_m = (am_s*Ds(i)**bm_s / am_r)**obmr
                p = Dc(j) / Ds_m
                if (p > 0.25 .or. Ds(i) < D0s .or. Dc(j) < 6.e-6 .or. vts < 1.e-3) then
                    t_Efsw(i,j) = 0.0
                else
                    stokes = Dc(j) * Dc(j) * vts * rho_w2 / (9.*1.718e-5*Ds_m)
                    reynolds = 9. * stokes / (p*p*rho_w2)

                    F = log(reynolds)
                    G = -0.1007_dp - 0.358_dp*F + 0.0261_dp*F*F
                    K0 = exp(G)
                    z = log(stokes / (K0+1.e-15_dp))
                    H = 0.1465_dp + 1.302_dp*z - 0.607_dp*z*z + 0.293_dp*z*z*z
                    yc0 = 2.0_dp / PI * atan(H)
                    Ef_sw = (yc0+p)*(yc0+p) / ((1.+p)*(1.+p))

                    t_Efsw(i,j) = max(0.0, min(real(Ef_sw, kind=wp), 0.95))
                endif
            enddo
        enddo

    end subroutine table_Efsw

    !=================================================================================================================
    ! Droplet evaporation
    ! Output:
    !   tpc_wev, tnc_wev

    subroutine table_dropEvap

        ! Local variables
        integer :: i, j, k, n
        real(dp), dimension(nbc) :: N_c, massc
        real(dp) :: summ, summ2, lamc, N0_c
        integer :: nu_c

        do n = 1, nbc
            massc(n) = am_r*Dc(n)**bm_r
        enddo

        do k = 1, nbc
            nu_c = min(nu_c_max, nint(nu_c_scale/t_Nc(k)) + nu_c_min)
            do j = 1, ntb_c
                lamc = (t_Nc(k)*am_r* ccg(2,nu_c)*ocg1(nu_c) / r_c(j))**obmr
                N0_c = t_Nc(k)*ocg1(nu_c) * lamc**cce(1,nu_c)
                do i = 1, nbc
                    N_c(i) = N0_c* Dc(i)**nu_c*exp(-lamc*Dc(i))*dtc(i)
                    summ = 0.
                    summ2 = 0.
                    do n = 1, i
                        summ = summ + massc(n)*N_c(n)
                        summ2 = summ2 + N_c(n)
                    enddo
                    tpc_wev(i,j,k) = summ
                    tnc_wev(i,j,k) = summ2
                enddo
            enddo
        enddo

    end subroutine table_dropEvap

    !=================================================================================================================
    ! Rain collecting graupel (and inverse).  Explicit CE integration.

    subroutine qr_acr_qg(NRHGtable)
        implicit none

        INTEGER, INTENT(IN) ::NRHGtable

        !..Local variables
        INTEGER:: i, j, k, m, n, n2, n3, idx_bg
        INTEGER:: km, km_s, km_e
        DOUBLE PRECISION, DIMENSION(nbg):: N_g
        DOUBLE PRECISION, DIMENSION(nbg,NRHGtable):: vg
        DOUBLE PRECISION, DIMENSION(nbr):: vr, N_r
        DOUBLE PRECISION:: N0_r, N0_g, lam_exp, lamg, lamr
        DOUBLE PRECISION:: massg, massr, dvg, dvr, t1, t2, z1, z2, y1, y2

        !+---+

        do n2 = 1, nbr
            !        vr(n2) = av_r*Dr(n2)**bv_r * DEXP(-fv_r*Dr(n2))
            vr(n2) = -0.1021 + 4.932E3*Dr(n2) - 0.9551E6*Dr(n2)*Dr(n2)     &
                + 0.07934E9*Dr(n2)*Dr(n2)*Dr(n2)                           &
                - 0.002362E12*Dr(n2)*Dr(n2)*Dr(n2)*Dr(n2)
        enddo

        do n3 = 1, NRHGtable
            do n = 1, nbg
                if (NRHGtable == NRHG) idx_bg = n3
                if (NRHGtable == NRHG1) idx_bg = idx_bg1
                vg(n,n3) = av_g(idx_bg)*Dg(n)**bv_g(idx_bg)
            enddo
        enddo

        km_s = 0
        km_e = ntb_r*ntb_r1 - 1

        do km = km_s, km_e
            m = km / ntb_r1 + 1
            k = mod( km , ntb_r1 ) + 1

            lam_exp = (N0r_exp(k)*am_r*crg(1)/r_r(m))**ore1
            lamr = lam_exp * (crg(3)*org2*org1)**obmr
            N0_r = N0r_exp(k)/(crg(2)*lam_exp) * lamr**cre(2)
            do n2 = 1, nbr
                N_r(n2) = N0_r*Dr(n2)**mu_r *DEXP(-lamr*Dr(n2))*dtr(n2)
            enddo

            do n3 = 1, NRHGtable
                if (NRHGtable == NRHG) idx_bg = n3
                if (NRHGtable == NRHG1) idx_bg = idx_bg1

                do j = 1, ntb_g
                    do i = 1, ntb_g1
                        lam_exp = (N0g_exp(i)*am_g(idx_bg)*cgg(1,1)/r_g(j))**oge1
                        lamg = lam_exp * (cgg(3,1)*ogg2*ogg1)**obmg
                        N0_g = N0g_exp(i)/(cgg(2,1)*lam_exp) * lamg**cge(2,1)
                        do n = 1, nbg
                            N_g(n) = N0_g*Dg(n)**mu_g * DEXP(-lamg*Dg(n))*dtg(n)
                        enddo

                        t1 = 0.0d0
                        t2 = 0.0d0
                        z1 = 0.0d0
                        z2 = 0.0d0
                        y1 = 0.0d0
                        y2 = 0.0d0
                        do n2 = 1, nbr
                            massr = am_r * Dr(n2)**bm_r
                            do n = 1, nbg
                                massg = am_g(idx_bg) * Dg(n)**bm_g

                                dvg = 0.5d0*((vr(n2) - vg(n,n3)) + DABS(vr(n2)-vg(n,n3)))
                                dvr = 0.5d0*((vg(n,n3) - vr(n2)) + DABS(vg(n,n3)-vr(n2)))

                                t1 = t1+ PI*.25*Ef_rg*(Dg(n)+Dr(n2))*(Dg(n)+Dr(n2)) &
                                    *dvg*massg * N_g(n)* N_r(n2)
                                z1 = z1+ PI*.25*Ef_rg*(Dg(n)+Dr(n2))*(Dg(n)+Dr(n2)) &
                                    *dvg*massr * N_g(n)* N_r(n2)
                                y1 = y1+ PI*.25*Ef_rg*(Dg(n)+Dr(n2))*(Dg(n)+Dr(n2)) &
                                    *dvg       * N_g(n)* N_r(n2)

                                t2 = t2+ PI*.25*Ef_rg*(Dg(n)+Dr(n2))*(Dg(n)+Dr(n2)) &
                                    *dvr*massr * N_g(n)* N_r(n2)
                                y2 = y2+ PI*.25*Ef_rg*(Dg(n)+Dr(n2))*(Dg(n)+Dr(n2)) &
                                    *dvr       * N_g(n)* N_r(n2)
                                z2 = z2+ PI*.25*Ef_rg*(Dg(n)+Dr(n2))*(Dg(n)+Dr(n2)) &
                                    *dvr*massg * N_g(n)* N_r(n2)
                            enddo
97                          continue
                        enddo
                        tcg_racg(i,j,n3,k,m) = t1
                        tmr_racg(i,j,n3,k,m) = DMIN1(z1, r_r(m)*1.0d0)
                        tcr_gacr(i,j,n3,k,m) = t2
                        tnr_racg(i,j,n3,k,m) = y1
                        tnr_gacr(i,j,n3,k,m) = y2
                    enddo
                enddo
            enddo
        enddo

    end subroutine qr_acr_qg

    !=================================================================================================================
    ! Rain collecting snow (and inverse).  Explicit CE integration.

    subroutine qr_acr_qs

        implicit none

        !..Local variables
        INTEGER:: i, j, k, m, n, n2
        INTEGER:: km, km_s, km_e
        DOUBLE PRECISION, DIMENSION(nbr):: vr, D1, N_r
        DOUBLE PRECISION, DIMENSION(nbs):: vs, N_s
        DOUBLE PRECISION:: loga_, a_, b_, second, M0, M2, M3, Mrat, oM3
        DOUBLE PRECISION:: N0_r, lam_exp, lamr, slam1, slam2
        DOUBLE PRECISION:: dvs, dvr, masss, massr
        DOUBLE PRECISION:: t1, t2, t3, t4, z1, z2, z3, z4
        DOUBLE PRECISION:: y1, y2, y3, y4

        !+---+

        do n2 = 1, nbr
            !        vr(n2) = av_r*Dr(n2)**bv_r * DEXP(-fv_r*Dr(n2))
            vr(n2) = -0.1021 + 4.932E3*Dr(n2) - 0.9551E6*Dr(n2)*Dr(n2)     &
                + 0.07934E9*Dr(n2)*Dr(n2)*Dr(n2)                           &
                - 0.002362E12*Dr(n2)*Dr(n2)*Dr(n2)*Dr(n2)
            D1(n2) = (vr(n2)/av_s)**(1./bv_s)
        enddo
        do n = 1, nbs
            vs(n) = 1.5*av_s*Ds(n)**bv_s * DEXP(-fv_s*Ds(n))
        enddo

        km_s = 0
        km_e = ntb_r*ntb_r1 - 1

        do km = km_s, km_e
            m = km / ntb_r1 + 1
            k = mod( km , ntb_r1 ) + 1

            lam_exp = (N0r_exp(k)*am_r*crg(1)/r_r(m))**ore1
            lamr = lam_exp * (crg(3)*org2*org1)**obmr
            N0_r = N0r_exp(k)/(crg(2)*lam_exp) * lamr**cre(2)
            do n2 = 1, nbr
                N_r(n2) = N0_r*Dr(n2)**mu_r * DEXP(-lamr*Dr(n2))*dtr(n2)
            enddo

            do j = 1, ntb_t
                do i = 1, ntb_s

                    !..From the bm_s moment, compute plus one moment.  If we are not
                    !.. using bm_s=2, then we must transform to the pure 2nd moment
                    !.. (variable called "second") and then to the bm_s+1 moment.

                    M2 = r_s(i)*oams *1.0d0
                    if (bm_s.gt.2.0-1.E-3 .and. bm_s.lt.2.0+1.E-3) then
                        loga_ = sa(1) + sa(2)*Tc(j) + sa(3)*bm_s &
                            + sa(4)*Tc(j)*bm_s + sa(5)*Tc(j)*Tc(j) &
                            + sa(6)*bm_s*bm_s + sa(7)*Tc(j)*Tc(j)*bm_s &
                            + sa(8)*Tc(j)*bm_s*bm_s + sa(9)*Tc(j)*Tc(j)*Tc(j) &
                            + sa(10)*bm_s*bm_s*bm_s
                        a_ = 10.0**loga_
                        b_ = sb(1) + sb(2)*Tc(j) + sb(3)*bm_s &
                            + sb(4)*Tc(j)*bm_s + sb(5)*Tc(j)*Tc(j) &
                            + sb(6)*bm_s*bm_s + sb(7)*Tc(j)*Tc(j)*bm_s &
                            + sb(8)*Tc(j)*bm_s*bm_s + sb(9)*Tc(j)*Tc(j)*Tc(j) &
                            + sb(10)*bm_s*bm_s*bm_s
                        second = (M2/a_)**(1./b_)
                    else
                        second = M2
                    endif

                    loga_ = sa(1) + sa(2)*Tc(j) + sa(3)*cse(1) &
                        + sa(4)*Tc(j)*cse(1) + sa(5)*Tc(j)*Tc(j) &
                        + sa(6)*cse(1)*cse(1) + sa(7)*Tc(j)*Tc(j)*cse(1) &
                        + sa(8)*Tc(j)*cse(1)*cse(1) + sa(9)*Tc(j)*Tc(j)*Tc(j) &
                        + sa(10)*cse(1)*cse(1)*cse(1)
                    a_ = 10.0**loga_
                    b_ = sb(1)+sb(2)*Tc(j)+sb(3)*cse(1) + sb(4)*Tc(j)*cse(1) &
                        + sb(5)*Tc(j)*Tc(j) + sb(6)*cse(1)*cse(1) &
                        + sb(7)*Tc(j)*Tc(j)*cse(1) + sb(8)*Tc(j)*cse(1)*cse(1) &
                        + sb(9)*Tc(j)*Tc(j)*Tc(j)+sb(10)*cse(1)*cse(1)*cse(1)
                    M3 = a_ * second**b_

                    oM3 = 1./M3
                    Mrat = M2*(M2*oM3)*(M2*oM3)*(M2*oM3)
                    M0   = (M2*oM3)**mu_s
                    slam1 = M2 * oM3 * Lam0
                    slam2 = M2 * oM3 * Lam1

                    do n = 1, nbs
                        N_s(n) = Mrat*(Kap0*DEXP(-slam1*Ds(n)) &
                            + Kap1*M0*Ds(n)**mu_s * DEXP(-slam2*Ds(n)))*dts(n)
                    enddo

                    t1 = 0.0d0
                    t2 = 0.0d0
                    t3 = 0.0d0
                    t4 = 0.0d0
                    z1 = 0.0d0
                    z2 = 0.0d0
                    z3 = 0.0d0
                    z4 = 0.0d0
                    y1 = 0.0d0
                    y2 = 0.0d0
                    y3 = 0.0d0
                    y4 = 0.0d0
                    do n2 = 1, nbr
                        massr = am_r * Dr(n2)**bm_r
                        do n = 1, nbs
                            masss = am_s * Ds(n)**bm_s

                            dvs = 0.5d0*((vr(n2) - vs(n)) + DABS(vr(n2)-vs(n)))
                            dvr = 0.5d0*((vs(n) - vr(n2)) + DABS(vs(n)-vr(n2)))
                            if (massr .gt. 1.5*masss) then
                                t1 = t1+ PI*.25*Ef_rs*(Ds(n)+Dr(n2))*(Ds(n)+Dr(n2)) &
                                    *dvs*masss * N_s(n)* N_r(n2)
                                z1 = z1+ PI*.25*Ef_rs*(Ds(n)+Dr(n2))*(Ds(n)+Dr(n2)) &
                                    *dvs*massr * N_s(n)* N_r(n2)
                                y1 = y1+ PI*.25*Ef_rs*(Ds(n)+Dr(n2))*(Ds(n)+Dr(n2)) &
                                    *dvs       * N_s(n)* N_r(n2)
                            else
                                t3 = t3+ PI*.25*Ef_rs*(Ds(n)+Dr(n2))*(Ds(n)+Dr(n2)) &
                                    *dvs*masss * N_s(n)* N_r(n2)
                                z3 = z3+ PI*.25*Ef_rs*(Ds(n)+Dr(n2))*(Ds(n)+Dr(n2)) &
                                    *dvs*massr * N_s(n)* N_r(n2)
                                y3 = y3+ PI*.25*Ef_rs*(Ds(n)+Dr(n2))*(Ds(n)+Dr(n2)) &
                                    *dvs       * N_s(n)* N_r(n2)
                            endif

                            if (massr .gt. 1.5*masss) then
                                t2 = t2+ PI*.25*Ef_rs*(Ds(n)+Dr(n2))*(Ds(n)+Dr(n2)) &
                                    *dvr*massr * N_s(n)* N_r(n2)
                                y2 = y2+ PI*.25*Ef_rs*(Ds(n)+Dr(n2))*(Ds(n)+Dr(n2)) &
                                    *dvr       * N_s(n)* N_r(n2)
                                z2 = z2+ PI*.25*Ef_rs*(Ds(n)+Dr(n2))*(Ds(n)+Dr(n2)) &
                                    *dvr*masss * N_s(n)* N_r(n2)
                            else
                                t4 = t4+ PI*.25*Ef_rs*(Ds(n)+Dr(n2))*(Ds(n)+Dr(n2)) &
                                    *dvr*massr * N_s(n)* N_r(n2)
                                y4 = y4+ PI*.25*Ef_rs*(Ds(n)+Dr(n2))*(Ds(n)+Dr(n2)) &
                                    *dvr       * N_s(n)* N_r(n2)
                                z4 = z4+ PI*.25*Ef_rs*(Ds(n)+Dr(n2))*(Ds(n)+Dr(n2)) &
                                    *dvr*masss * N_s(n)* N_r(n2)
                            endif

                        enddo
                    enddo
                    tcs_racs1(i,j,k,m) = t1
                    tmr_racs1(i,j,k,m) = DMIN1(z1, r_r(m)*1.0d0)
                    tcs_racs2(i,j,k,m) = t3
                    tmr_racs2(i,j,k,m) = z3
                    tcr_sacr1(i,j,k,m) = t2
                    tms_sacr1(i,j,k,m) = z2
                    tcr_sacr2(i,j,k,m) = t4
                    tms_sacr2(i,j,k,m) = z4
                    tnr_racs1(i,j,k,m) = y1
                    tnr_racs2(i,j,k,m) = y3
                    tnr_sacr1(i,j,k,m) = y2
                    tnr_sacr2(i,j,k,m) = y4
                enddo
            enddo
        enddo

    end subroutine qr_acr_qs

    !=================================================================================================================
    ! This is a literal adaptation of Bigg (1954) probability of drops of a particular volume freezing.
    ! Given this probability, simply freeze the proportion of drops summing their masses.

    subroutine freezeH2O

        implicit none

        !..Local variables
        INTEGER:: i, j, k, m, n, n2
        DOUBLE PRECISION :: N_r, N_c
        DOUBLE PRECISION, DIMENSION(nbr):: massr
        DOUBLE PRECISION, DIMENSION(nbc):: massc
        DOUBLE PRECISION:: sum1, sum2, sumn1, sumn2, &
            prob, vol, Texp, orho_w, &
            lam_exp, lamr, N0_r, lamc, N0_c, y
        INTEGER:: nu_c
        REAL:: T_adjust

        orho_w = 1./rho_w2

        do n2 = 1, nbr
            massr(n2) = am_r*Dr(n2)**bm_r
        enddo
        do n = 1, nbc
            massc(n) = am_r*Dc(n)**bm_r
        enddo

        !..Freeze water (smallest drops become cloud ice, otherwise graupel).
        do m = 1, ntb_IN
            T_adjust = MAX(-3.0, MIN(3.0 - ALOG10(Nt_IN(m)), 3.0))
            do k = 1, 45
                !         print*, ' Freezing water for temp = ', -k
                Texp = DEXP( REAL(k,KIND=dp) - T_adjust*1.0D0 ) - 1.0D0
                do j = 1, ntb_r1
                    do i = 1, ntb_r
                        lam_exp = (N0r_exp(j)*am_r*crg(1)/r_r(i))**ore1
                        lamr = lam_exp * (crg(3)*org2*org1)**obmr
                        N0_r = N0r_exp(j)/(crg(2)*lam_exp) * lamr**cre(2)
                        sum1 = 0.0d0
                        sum2 = 0.0d0
                        sumn1 = 0.0d0
                        sumn2 = 0.0d0
                        do n2 = nbr, 1, -1
                            N_r = N0_r*Dr(n2)**mu_r*DEXP(-lamr*Dr(n2))*dtr(n2)
                            vol = massr(n2)*orho_w
                            prob = MAX(0.0D0, 1.0D0 - DEXP(-120.0D0*vol*5.2D-4 * Texp))
                            if (massr(n2) .lt. xm0g) then
                                sumn1 = sumn1 + prob*N_r
                                sum1 = sum1 + prob*N_r*massr(n2)
                            else
                                sumn2 = sumn2 + prob*N_r
                                sum2 = sum2 + prob*N_r*massr(n2)
                            endif
                            if ((sum1+sum2).ge.r_r(i)) EXIT
                        enddo
                        tpi_qrfz(i,j,k,m) = sum1
                        tni_qrfz(i,j,k,m) = sumn1
                        tpg_qrfz(i,j,k,m) = sum2
                        tnr_qrfz(i,j,k,m) = sumn2
                    enddo
                enddo

                do j = 1, nbc
                    nu_c = MIN(15, NINT(1000.E6/t_Nc(j)) + 2)
                    do i = 1, ntb_c
                        lamc = (t_Nc(j)*am_r* ccg(2,nu_c) * ocg1(nu_c) / r_c(i))**obmr
                        N0_c = t_Nc(j)*ocg1(nu_c) * lamc**cce(1,nu_c)
                        sum1 = 0.0d0
                        sumn2 = 0.0d0
                        do n = nbc, 1, -1
                            vol = massc(n)*orho_w
                            prob = MAX(0.0D0, 1.0D0 - DEXP(-120.0D0*vol*5.2D-4 * Texp))
                            N_c = N0_c*Dc(n)**nu_c*EXP(-lamc*Dc(n))*dtc(n)
                            sumn2 = MIN(t_Nc(j), sumn2 + prob*N_c)
                            sum1 = sum1 + prob*N_c*massc(n)
                            if (sum1 .ge. r_c(i)) EXIT
                        enddo
                        tpi_qcfz(i,j,k,m) = sum1
                        tni_qcfz(i,j,k,m) = sumn2
                    enddo
                enddo
            enddo
        enddo

    end subroutine freezeH2O

    !=================================================================================================================
    ! Cloud ice converting to snow since portion greater than min snow size.  Given cloud ice content (kg/m**3),
    ! number concentration (#/m**3) and gamma shape parameter, mu_i, break the distrib into bins and figure out
    ! the mass/number of ice with sizes larger than D0s.  Also, compute incomplete gamma function for the
    ! integration of ice depositional growth from diameter=0 to D0s.  Amount of ice depositional growth is this
    ! portion of distrib while larger diameters contribute to snow growth (as in Harrington et al. 1995).

    subroutine qi_aut_qs

        implicit none

        !..Local variables
        INTEGER:: i, j, n2
        DOUBLE PRECISION, DIMENSION(nbi):: N_i
        DOUBLE PRECISION:: N0_i, lami, Di_mean, t1, t2
        REAL:: xlimit_intg

        do j = 1, ntb_i1
            do i = 1, ntb_i
                lami = (am_i*cig(2)*oig1*Nt_i(j)/r_i(i))**obmi
                Di_mean = (bm_i + mu_i + 1.) / lami
                N0_i = Nt_i(j)*oig1 * lami**cie(1)
                t1 = 0.0d0
                t2 = 0.0d0
                if (SNGL(Di_mean) .gt. 5.*D0s) then
                    t1 = r_i(i)
                    t2 = Nt_i(j)
                    tpi_ide(i,j) = 0.0D0
                elseif (SNGL(Di_mean) .lt. D0i) then
                    t1 = 0.0D0
                    t2 = 0.0D0
                    tpi_ide(i,j) = 1.0D0
                else
                    xlimit_intg = lami*D0s
                    tpi_ide(i,j) = gamma_p(mu_i+2.0, xlimit_intg) * 1.0D0
                    do n2 = 1, nbi
                        N_i(n2) = N0_i*Di(n2)**mu_i * DEXP(-lami*Di(n2))*dti(n2)
                        if (Di(n2).ge.D0s) then
                            t1 = t1 + N_i(n2) * am_i*Di(n2)**bm_i
                            t2 = t2 + N_i(n2)
                        endif
                    enddo
                endif
                tps_iaus(i,j) = t1
                tni_iaus(i,j) = t2
            enddo
        enddo

    end subroutine qi_aut_qs

#if defined (ccpp_default)
    !=================================================================================================================
    !..Fill the table of CCN activation data created from parcel model run
    !.. by Trude Eidhammer with inputs of aerosol number concentration,
    !.. vertical velocity, temperature, lognormal mean aerosol radius, and
    !.. hygroscopicity, kappa.  The data are read from external file and
    !.. contain activated fraction of CCN for given conditions.
    !+---+-----------------------------------------------------------------+
    !+---+-----------------------------------------------------------------+
    !>\ingroup aathompson
    !! Fill the table of CCN activation data created from parcel model run
    !! by Trude Eidhammer with inputs of aerosol number concentration,
    !! vertical velocity, temperature, lognormal mean aerosol radius, and
    !! hygroscopicity, kappa.  The data are read from external file and
    !! contain activated fraction of CCN for given conditions.

    subroutine table_ccnAct(errmess,errflag)

        implicit none

        !..Error handling variables
        CHARACTER(len=*), INTENT(INOUT) :: errmess
        INTEGER,          INTENT(INOUT) :: errflag

        !..Local variables
        INTEGER:: iunit_mp_th1, i
        LOGICAL:: opened

        iunit_mp_th1 = -1
        DO i = 20,99
            INQUIRE ( i , OPENED = opened )
            IF ( .NOT. opened ) THEN
                iunit_mp_th1 = i
                GOTO 2010
            ENDIF
        ENDDO
2010    CONTINUE
        IF ( iunit_mp_th1 < 0 ) THEN
            write(0,*) 'module_mp_tempo: table_ccnAct: '//   &
                'Can not find unused fortran unit to read in lookup table.'
            return
        ENDIF

        !WRITE(*, '(A,I2)') 'module_mp_tempo: opening CCN_ACTIVATE.BIN on unit ',iunit_mp_th1
        OPEN(iunit_mp_th1,FILE='CCN_ACTIVATE.BIN',                      &
            FORM='UNFORMATTED',STATUS='OLD',CONVERT='BIG_ENDIAN',ERR=9009)

        !sms$serial begin
        READ(iunit_mp_th1,ERR=9010) tnccn_act
        !sms$serial end

        RETURN
9009    CONTINUE
        WRITE( errmess , '(A,I2)' ) 'module_mp_tempo: error opening CCN_ACTIVATE.BIN on unit ',iunit_mp_th1
        errflag = 1
        RETURN
9010    CONTINUE
        WRITE( errmess , '(A,I2)' ) 'module_mp_tempo: error reading CCN_ACTIVATE.BIN on unit ',iunit_mp_th1
        errflag = 1
        RETURN

    end subroutine table_ccnAct

    !=================================================================================================================
    ! Rain collecting graupel (and inverse).  Explicit CE integration.
    subroutine qr_acr_qg_par(NRHGtable, qrqg_file)

        implicit none

        INTEGER, INTENT(IN) ::NRHGtable

        !..Local variables
        INTEGER:: i, j, k, m, n, n2, n3, idx_bg
        INTEGER:: km, km_s, km_e
        DOUBLE PRECISION, DIMENSION(nbg):: N_g
        DOUBLE PRECISION, DIMENSION(nbg,NRHGtable):: vg
        DOUBLE PRECISION, DIMENSION(nbr):: vr, N_r
        DOUBLE PRECISION:: N0_r, N0_g, lam_exp, lamg, lamr
        DOUBLE PRECISION:: massg, massr, dvg, dvr, t1, t2, z1, z2, y1, y2
        LOGICAL force_read_thompson, write_thompson_tables
        LOGICAL lexist,lopen
        INTEGER good,ierr
        character(len=*), intent(in) :: qrqg_file

        force_read_thompson = .false.
        write_thompson_tables = .false.
        !+---+

        good = 0
        INQUIRE(FILE=qrqg_file, EXIST=lexist)
#ifdef MPI
        call MPI_BARRIER(mpi_communicator,ierr)
#endif
        IF ( lexist ) THEN
            OPEN(63,file=qrqg_file,form="unformatted",err=1234)
            !sms$serial begin
            READ(63,err=1234) tcg_racg
            READ(63,err=1234) tmr_racg
            READ(63,err=1234) tcr_gacr
!!            READ(63,err=1234) tmg_gacr
            READ(63,err=1234) tnr_racg
            READ(63,err=1234) tnr_gacr
            !sms$serial end
            good = 1
1234        CONTINUE
            IF ( good .NE. 1 ) THEN
                INQUIRE(63,opened=lopen)
                IF (lopen) THEN
                    IF( force_read_thompson ) THEN
                        write(0,*) "Error reading "//qrqg_file//" Aborting because force_read_thompson is .true."
                        return
                    ENDIF
                    CLOSE(63)
                ELSE
                    IF( force_read_thompson ) THEN
                        write(0,*) "Error opening "//qrqg_file//" Aborting because force_read_thompson is .true."
                        return
                    ENDIF
                ENDIF
            ELSE
                INQUIRE(63,opened=lopen)
                IF (lopen) THEN
                    CLOSE(63)
                ENDIF
            ENDIF
        ELSE
            IF( force_read_thompson ) THEN
                write(0,*) "Non-existent "//qrqg_file//" Aborting because force_read_thompson is .true."
                return
            ENDIF
        ENDIF

        IF (.NOT. good .EQ. 1 ) THEN
            if (thompson_table_writer) then
                write_thompson_tables = .true.
                write(0,*) "ThompMP: computing qr_acr_qg"
            endif
            do n2 = 1, nbr
                !        vr(n2) = av_r*Dr(n2)**bv_r * DEXP(-fv_r*Dr(n2))
                vr(n2) = -0.1021 + 4.932E3*Dr(n2) - 0.9551E6*Dr(n2)*Dr(n2)     &
                    + 0.07934E9*Dr(n2)*Dr(n2)*Dr(n2)                          &
                    - 0.002362E12*Dr(n2)*Dr(n2)*Dr(n2)*Dr(n2)
            enddo
            !   do n = 1, nbg
            !    vg(n) = av_g*Dg(n)**bv_g
            !   enddo

            do n3 = 1, NRHGtable
               do n = 1, nbg
                  if (NRHGtable == NRHG) idx_bg = n3
                  if (NRHGtable == NRHG1) idx_bg = idx_bg1
                  vg(n,n3) = av_g(idx_bg)*Dg(n)**bv_g(idx_bg)
               enddo
            enddo

            !..Note values returned from wrf_dm_decomp1d are zero-based, add 1 for
            !.. fortran indices.  J. Michalakes, 2009Oct30.

#if ( defined( DM_PARALLEL ) && ( ! defined( STUBMPI ) ) )
            CALL wrf_dm_decomp1d ( ntb_r*ntb_r1, km_s, km_e )
#else
            km_s = 0
            km_e = ntb_r*ntb_r1 - 1
#endif

            do km = km_s, km_e
                m = km / ntb_r1 + 1
                k = mod( km , ntb_r1 ) + 1

                lam_exp = (N0r_exp(k)*am_r*crg(1)/r_r(m))**ore1
                lamr = lam_exp * (crg(3)*org2*org1)**obmr
                N0_r = N0r_exp(k)/(crg(2)*lam_exp) * lamr**cre(2)
                do n2 = 1, nbr
                    N_r(n2) = N0_r*Dr(n2)**mu_r *DEXP(-lamr*Dr(n2))*dtr(n2)
                enddo

                do n3 = 1, NRHGtable
                   if (NRHGtable == NRHG) idx_bg = n3
                   if (NRHGtable == NRHG1) idx_bg = idx_bg1

                    do j = 1, ntb_g
                        do i = 1, ntb_g1
                            lam_exp = (N0g_exp(i)*am_g(idx_bg)*cgg(1,1)/r_g(j))**oge1
                            lamg = lam_exp * (cgg(3,1)*ogg2*ogg1)**obmg
                            N0_g = N0g_exp(i)/(cgg(2,1)*lam_exp) * lamg**cge(2,1)
                            do n = 1, nbg
                                N_g(n) = N0_g*Dg(n)**mu_g * DEXP(-lamg*Dg(n))*dtg(n)
                            enddo

                            t1 = 0.0d0
                            t2 = 0.0d0
                            z1 = 0.0d0
                            z2 = 0.0d0
                            y1 = 0.0d0
                            y2 = 0.0d0
                            do n2 = 1, nbr
                                massr = am_r * Dr(n2)**bm_r
                                do n = 1, nbg
                                   massg = am_g(idx_bg) * Dg(n)**bm_g

                                    dvg = 0.5d0*((vr(n2) - vg(n,n3)) + DABS(vr(n2)-vg(n,n3)))
                                    dvr = 0.5d0*((vg(n,n3) - vr(n2)) + DABS(vg(n,n3)-vr(n2)))

                                    t1 = t1+ PI*.25*Ef_rg*(Dg(n)+Dr(n2))*(Dg(n)+Dr(n2)) &
                                        *dvg*massg * N_g(n)* N_r(n2)
                                    z1 = z1+ PI*.25*Ef_rg*(Dg(n)+Dr(n2))*(Dg(n)+Dr(n2)) &
                                        *dvg*massr * N_g(n)* N_r(n2)
                                    y1 = y1+ PI*.25*Ef_rg*(Dg(n)+Dr(n2))*(Dg(n)+Dr(n2)) &
                                        *dvg       * N_g(n)* N_r(n2)

                                    t2 = t2+ PI*.25*Ef_rg*(Dg(n)+Dr(n2))*(Dg(n)+Dr(n2)) &
                                        *dvr*massr * N_g(n)* N_r(n2)
                                    y2 = y2+ PI*.25*Ef_rg*(Dg(n)+Dr(n2))*(Dg(n)+Dr(n2)) &
                                        *dvr       * N_g(n)* N_r(n2)
                                    z2 = z2+ PI*.25*Ef_rg*(Dg(n)+Dr(n2))*(Dg(n)+Dr(n2)) &
                                        *dvr*massg * N_g(n)* N_r(n2)
                                enddo
97                              continue
                            enddo
                            tcg_racg(i,j,n3,k,m) = t1
                            tmr_racg(i,j,n3,k,m) = DMIN1(z1, r_r(m)*1.0d0)
                            tcr_gacr(i,j,n3,k,m) = t2
!!                            tmg_gacr(i,j,k,m) = DMIN1(z2, r_g(j)*1.0d0)
                            tnr_racg(i,j,n3,k,m) = y1
                            tnr_gacr(i,j,n3,k,m) = y2
                        enddo
                    enddo
                enddo
            enddo
            IF ( write_thompson_tables ) THEN
                write(0,*) "Writing "//qrqg_file//" in Tempo MP init"
                OPEN(63,file=qrqg_file,form="unformatted",err=9234)
                WRITE(63,err=9234) tcg_racg
                WRITE(63,err=9234) tmr_racg
                WRITE(63,err=9234) tcr_gacr
                WRITE(63,err=9234) tnr_racg
                WRITE(63,err=9234) tnr_gacr
                CLOSE(63)
                RETURN    ! ----- RETURN
9234            CONTINUE
                write(0,*) "Error writing "//qrqg_file
                return
            ENDIF
        ENDIF

    end subroutine qr_acr_qg_par

    !=================================================================================================================
    ! Rain collecting snow (and inverse).  Explicit CE integration.

    subroutine qr_acr_qs_par

        implicit none

        !..Local variables
        INTEGER:: i, j, k, m, n, n2
        INTEGER:: km, km_s, km_e
        DOUBLE PRECISION, DIMENSION(nbr):: vr, D1, N_r
        DOUBLE PRECISION, DIMENSION(nbs):: vs, N_s
        DOUBLE PRECISION:: loga_, a_, b_, second, M0, M2, M3, Mrat, oM3
        DOUBLE PRECISION:: N0_r, lam_exp, lamr, slam1, slam2
        DOUBLE PRECISION:: dvs, dvr, masss, massr
        DOUBLE PRECISION:: t1, t2, t3, t4, z1, z2, z3, z4
        DOUBLE PRECISION:: y1, y2, y3, y4
        LOGICAL force_read_thompson, write_thompson_tables
        LOGICAL lexist,lopen
        INTEGER good,ierr

        !+---+

        force_read_thompson = .false.
        write_thompson_tables = .false.

        good = 0
        INQUIRE(FILE=qr_acr_qs_file, EXIST=lexist)
#ifdef MPI
        call MPI_BARRIER(mpi_communicator,ierr)
#endif
        IF ( lexist ) THEN
            !write(0,*) "ThompMP: read "//qr_acr_qs_file//" instead of computing"
            OPEN(63,file=qr_acr_qs_file,form="unformatted",err=1234)
            !sms$serial begin
            READ(63,err=1234)tcs_racs1
            READ(63,err=1234)tmr_racs1
            READ(63,err=1234)tcs_racs2
            READ(63,err=1234)tmr_racs2
            READ(63,err=1234)tcr_sacr1
            READ(63,err=1234)tms_sacr1
            READ(63,err=1234)tcr_sacr2
            READ(63,err=1234)tms_sacr2
            READ(63,err=1234)tnr_racs1
            READ(63,err=1234)tnr_racs2
            READ(63,err=1234)tnr_sacr1
            READ(63,err=1234)tnr_sacr2
            !sms$serial end
            good = 1
1234        CONTINUE
            IF ( good .NE. 1 ) THEN
                INQUIRE(63,opened=lopen)
                IF (lopen) THEN
                    IF( force_read_thompson ) THEN
                        write(0,*) "Error reading "//qr_acr_qs_file//" Aborting because force_read_thompson is .true."
                        return
                    ENDIF
                    CLOSE(63)
                ELSE
                    IF( force_read_thompson ) THEN
                        write(0,*) "Error opening "//qr_acr_qs_file//" Aborting because force_read_thompson is .true."
                        return
                    ENDIF
                ENDIF
            ELSE
                INQUIRE(63,opened=lopen)
                IF (lopen) THEN
                    CLOSE(63)
                ENDIF
            ENDIF
        ELSE
            IF( force_read_thompson ) THEN
                write(0,*) "Non-existent "//qr_acr_qs_file//" Aborting because force_read_thompson is .true."
                return
            ENDIF
        ENDIF

        IF (.NOT. good .EQ. 1 ) THEN
            if (thompson_table_writer) then
                write_thompson_tables = .true.
                write(0,*) "ThompMP: computing qr_acr_qs"
            endif
            do n2 = 1, nbr
                !        vr(n2) = av_r*Dr(n2)**bv_r * DEXP(-fv_r*Dr(n2))
                vr(n2) = -0.1021 + 4.932E3*Dr(n2) - 0.9551E6*Dr(n2)*Dr(n2)     &
                    + 0.07934E9*Dr(n2)*Dr(n2)*Dr(n2)                          &
                    - 0.002362E12*Dr(n2)*Dr(n2)*Dr(n2)*Dr(n2)
                D1(n2) = (vr(n2)/av_s)**(1./bv_s)
            enddo
            do n = 1, nbs
                vs(n) = 1.5*av_s*Ds(n)**bv_s * DEXP(-fv_s*Ds(n))
            enddo

            !..Note values returned from wrf_dm_decomp1d are zero-based, add 1 for
            !.. fortran indices.  J. Michalakes, 2009Oct30.

#if ( defined( DM_PARALLEL ) && ( ! defined( STUBMPI ) ) )
            CALL wrf_dm_decomp1d ( ntb_r*ntb_r1, km_s, km_e )
#else
            km_s = 0
            km_e = ntb_r*ntb_r1 - 1
#endif

            do km = km_s, km_e
                m = km / ntb_r1 + 1
                k = mod( km , ntb_r1 ) + 1

                lam_exp = (N0r_exp(k)*am_r*crg(1)/r_r(m))**ore1
                lamr = lam_exp * (crg(3)*org2*org1)**obmr
                N0_r = N0r_exp(k)/(crg(2)*lam_exp) * lamr**cre(2)
                do n2 = 1, nbr
                    N_r(n2) = N0_r*Dr(n2)**mu_r * DEXP(-lamr*Dr(n2))*dtr(n2)
                enddo

                do j = 1, ntb_t
                    do i = 1, ntb_s

                        !..From the bm_s moment, compute plus one moment.  If we are not
                        !.. using bm_s=2, then we must transform to the pure 2nd moment
                        !.. (variable called "second") and then to the bm_s+1 moment.

                        M2 = r_s(i)*oams *1.0d0
                        if (bm_s.gt.2.0-1.E-3 .and. bm_s.lt.2.0+1.E-3) then
                            loga_ = sa(1) + sa(2)*Tc(j) + sa(3)*bm_s &
                                + sa(4)*Tc(j)*bm_s + sa(5)*Tc(j)*Tc(j) &
                                + sa(6)*bm_s*bm_s + sa(7)*Tc(j)*Tc(j)*bm_s &
                                + sa(8)*Tc(j)*bm_s*bm_s + sa(9)*Tc(j)*Tc(j)*Tc(j) &
                                + sa(10)*bm_s*bm_s*bm_s
                            a_ = 10.0**loga_
                            b_ = sb(1) + sb(2)*Tc(j) + sb(3)*bm_s &
                                + sb(4)*Tc(j)*bm_s + sb(5)*Tc(j)*Tc(j) &
                                + sb(6)*bm_s*bm_s + sb(7)*Tc(j)*Tc(j)*bm_s &
                                + sb(8)*Tc(j)*bm_s*bm_s + sb(9)*Tc(j)*Tc(j)*Tc(j) &
                                + sb(10)*bm_s*bm_s*bm_s
                            second = (M2/a_)**(1./b_)
                        else
                            second = M2
                        endif

                        loga_ = sa(1) + sa(2)*Tc(j) + sa(3)*cse(1) &
                            + sa(4)*Tc(j)*cse(1) + sa(5)*Tc(j)*Tc(j) &
                            + sa(6)*cse(1)*cse(1) + sa(7)*Tc(j)*Tc(j)*cse(1) &
                            + sa(8)*Tc(j)*cse(1)*cse(1) + sa(9)*Tc(j)*Tc(j)*Tc(j) &
                            + sa(10)*cse(1)*cse(1)*cse(1)
                        a_ = 10.0**loga_
                        b_ = sb(1)+sb(2)*Tc(j)+sb(3)*cse(1) + sb(4)*Tc(j)*cse(1) &
                            + sb(5)*Tc(j)*Tc(j) + sb(6)*cse(1)*cse(1) &
                            + sb(7)*Tc(j)*Tc(j)*cse(1) + sb(8)*Tc(j)*cse(1)*cse(1) &
                            + sb(9)*Tc(j)*Tc(j)*Tc(j)+sb(10)*cse(1)*cse(1)*cse(1)
                        M3 = a_ * second**b_

                        oM3 = 1./M3
                        Mrat = M2*(M2*oM3)*(M2*oM3)*(M2*oM3)
                        M0   = (M2*oM3)**mu_s
                        slam1 = M2 * oM3 * Lam0
                        slam2 = M2 * oM3 * Lam1

                        do n = 1, nbs
                            N_s(n) = Mrat*(Kap0*DEXP(-slam1*Ds(n)) &
                                + Kap1*M0*Ds(n)**mu_s * DEXP(-slam2*Ds(n)))*dts(n)
                        enddo

                        t1 = 0.0d0
                        t2 = 0.0d0
                        t3 = 0.0d0
                        t4 = 0.0d0
                        z1 = 0.0d0
                        z2 = 0.0d0
                        z3 = 0.0d0
                        z4 = 0.0d0
                        y1 = 0.0d0
                        y2 = 0.0d0
                        y3 = 0.0d0
                        y4 = 0.0d0
                        do n2 = 1, nbr
                            massr = am_r * Dr(n2)**bm_r
                            do n = 1, nbs
                                masss = am_s * Ds(n)**bm_s

                                dvs = 0.5d0*((vr(n2) - vs(n)) + DABS(vr(n2)-vs(n)))
                                dvr = 0.5d0*((vs(n) - vr(n2)) + DABS(vs(n)-vr(n2)))

                                if (massr .gt. 1.5*masss) then
                                    t1 = t1+ PI*.25*Ef_rs*(Ds(n)+Dr(n2))*(Ds(n)+Dr(n2)) &
                                        *dvs*masss * N_s(n)* N_r(n2)
                                    z1 = z1+ PI*.25*Ef_rs*(Ds(n)+Dr(n2))*(Ds(n)+Dr(n2)) &
                                        *dvs*massr * N_s(n)* N_r(n2)
                                    y1 = y1+ PI*.25*Ef_rs*(Ds(n)+Dr(n2))*(Ds(n)+Dr(n2)) &
                                        *dvs       * N_s(n)* N_r(n2)
                                else
                                    t3 = t3+ PI*.25*Ef_rs*(Ds(n)+Dr(n2))*(Ds(n)+Dr(n2)) &
                                        *dvs*masss * N_s(n)* N_r(n2)
                                    z3 = z3+ PI*.25*Ef_rs*(Ds(n)+Dr(n2))*(Ds(n)+Dr(n2)) &
                                        *dvs*massr * N_s(n)* N_r(n2)
                                    y3 = y3+ PI*.25*Ef_rs*(Ds(n)+Dr(n2))*(Ds(n)+Dr(n2)) &
                                        *dvs       * N_s(n)* N_r(n2)
                                endif

                                if (massr .gt. 1.5*masss) then
                                    t2 = t2+ PI*.25*Ef_rs*(Ds(n)+Dr(n2))*(Ds(n)+Dr(n2)) &
                                        *dvr*massr * N_s(n)* N_r(n2)
                                    y2 = y2+ PI*.25*Ef_rs*(Ds(n)+Dr(n2))*(Ds(n)+Dr(n2)) &
                                        *dvr       * N_s(n)* N_r(n2)
                                    z2 = z2+ PI*.25*Ef_rs*(Ds(n)+Dr(n2))*(Ds(n)+Dr(n2)) &
                                        *dvr*masss * N_s(n)* N_r(n2)
                                else
                                    t4 = t4+ PI*.25*Ef_rs*(Ds(n)+Dr(n2))*(Ds(n)+Dr(n2)) &
                                        *dvr*massr * N_s(n)* N_r(n2)
                                    y4 = y4+ PI*.25*Ef_rs*(Ds(n)+Dr(n2))*(Ds(n)+Dr(n2)) &
                                        *dvr       * N_s(n)* N_r(n2)
                                    z4 = z4+ PI*.25*Ef_rs*(Ds(n)+Dr(n2))*(Ds(n)+Dr(n2)) &
                                        *dvr*masss * N_s(n)* N_r(n2)
                                endif

                            enddo
                        enddo
                        tcs_racs1(i,j,k,m) = t1
                        tmr_racs1(i,j,k,m) = DMIN1(z1, r_r(m)*1.0d0)
                        tcs_racs2(i,j,k,m) = t3
                        tmr_racs2(i,j,k,m) = z3
                        tcr_sacr1(i,j,k,m) = t2
                        tms_sacr1(i,j,k,m) = z2
                        tcr_sacr2(i,j,k,m) = t4
                        tms_sacr2(i,j,k,m) = z4
                        tnr_racs1(i,j,k,m) = y1
                        tnr_racs2(i,j,k,m) = y3
                        tnr_sacr1(i,j,k,m) = y2
                        tnr_sacr2(i,j,k,m) = y4
                    enddo
                enddo
            enddo

            IF ( write_thompson_tables ) THEN
                write(0,*) "Writing "//qr_acr_qs_file//" in Tempo MP init"
                OPEN(63,file=qr_acr_qs_file,form="unformatted",err=9234)
                WRITE(63,err=9234)tcs_racs1
                WRITE(63,err=9234)tmr_racs1
                WRITE(63,err=9234)tcs_racs2
                WRITE(63,err=9234)tmr_racs2
                WRITE(63,err=9234)tcr_sacr1
                WRITE(63,err=9234)tms_sacr1
                WRITE(63,err=9234)tcr_sacr2
                WRITE(63,err=9234)tms_sacr2
                WRITE(63,err=9234)tnr_racs1
                WRITE(63,err=9234)tnr_racs2
                WRITE(63,err=9234)tnr_sacr1
                WRITE(63,err=9234)tnr_sacr2
                CLOSE(63)
                RETURN    ! ----- RETURN
9234            CONTINUE
                write(0,*) "Error writing "//qr_acr_qs_file
            ENDIF
        ENDIF

    end subroutine qr_acr_qs_par

    !=================================================================================================================
    !! This is a literal adaptation of Bigg (1954) probability of drops of
    !! a particular volume freezing.  Given this probability, simply freeze
    !! the proportion of drops summing their masses.

    subroutine freezeH2O_par(threads)

        implicit none

        !..Interface variables
        INTEGER, INTENT(IN):: threads

        !..Local variables
        INTEGER:: i, j, k, m, n, n2
        DOUBLE PRECISION:: N_r, N_c
        DOUBLE PRECISION, DIMENSION(nbr):: massr
        DOUBLE PRECISION, DIMENSION(nbc):: massc
        DOUBLE PRECISION:: sum1, sum2, sumn1, sumn2, &
            prob, vol, Texp, orho_w, &
            lam_exp, lamr, N0_r, lamc, N0_c, y
        INTEGER:: nu_c
        REAL:: T_adjust
        LOGICAL force_read_thompson, write_thompson_tables
        LOGICAL lexist,lopen
        INTEGER good,ierr

        !+---+
        force_read_thompson = .false.
        write_thompson_tables = .false.

        good = 0
        INQUIRE(FILE=freeze_h2o_file,EXIST=lexist)
#ifdef MPI
        call MPI_BARRIER(mpi_communicator,ierr)
#endif
        IF ( lexist ) THEN
            !write(0,*) "ThompMP: read "//freeze_h2o_file//" instead of computing"
            OPEN(63,file=freeze_h2o_file,form="unformatted",err=1234)
            !sms$serial begin
            READ(63,err=1234)tpi_qrfz
            READ(63,err=1234)tni_qrfz
            READ(63,err=1234)tpg_qrfz
            READ(63,err=1234)tnr_qrfz
            READ(63,err=1234)tpi_qcfz
            READ(63,err=1234)tni_qcfz
            !sms$serial end
            good = 1
1234        CONTINUE
            IF ( good .NE. 1 ) THEN
                INQUIRE(63,opened=lopen)
                IF (lopen) THEN
                    IF( force_read_thompson ) THEN
                        write(0,*) "Error reading "//freeze_h2o_file//" Aborting because force_read_thompson is .true."
                        return
                    ENDIF
                    CLOSE(63)
                ELSE
                    IF( force_read_thompson ) THEN
                        write(0,*) "Error opening "//freeze_h2o_file//" Aborting because force_read_thompson is .true."
                        return
                    ENDIF
                ENDIF
            ELSE
                INQUIRE(63,opened=lopen)
                IF (lopen) THEN
                    CLOSE(63)
                ENDIF
            ENDIF
        ELSE
            IF( force_read_thompson ) THEN
                write(0,*) "Non-existent "//freeze_h2o_file//" Aborting because force_read_thompson is .true."
                return
            ENDIF
        ENDIF

        IF (.NOT. good .EQ. 1 ) THEN
            if (thompson_table_writer) then
                write_thompson_tables = .true.
                write(0,*) "ThompMP: computing freezeH2O"
            endif

            orho_w = 1./rho_w2

            do n2 = 1, nbr
                massr(n2) = am_r*Dr(n2)**bm_r
            enddo
            do n = 1, nbc
                massc(n) = am_r*Dc(n)**bm_r
            enddo

            !..Freeze water (smallest drops become cloud ice, otherwise graupel).
            do m = 1, ntb_IN
                T_adjust = MAX(-3.0, MIN(3.0 - ALOG10(Nt_IN(m)), 3.0))
                do k = 1, 45
                    !         print*, ' Freezing water for temp = ', -k
                    Texp = DEXP( DFLOAT(k) - T_adjust*1.0D0 ) - 1.0D0
                    !$OMP PARALLEL DO SCHEDULE(dynamic) num_threads(threads) &
                    !$OMP PRIVATE(j,i,lam_exp,lamr,N0_r,sum1,sum2,sumn1,sumn2,n2,N_r,vol,prob)
                    do j = 1, ntb_r1
                        do i = 1, ntb_r
                            lam_exp = (N0r_exp(j)*am_r*crg(1)/r_r(i))**ore1
                            lamr = lam_exp * (crg(3)*org2*org1)**obmr
                            N0_r = N0r_exp(j)/(crg(2)*lam_exp) * lamr**cre(2)
                            sum1 = 0.0d0
                            sum2 = 0.0d0
                            sumn1 = 0.0d0
                            sumn2 = 0.0d0
                            do n2 = nbr, 1, -1
                                N_r = N0_r*Dr(n2)**mu_r*DEXP(-lamr*Dr(n2))*dtr(n2)
                                vol = massr(n2)*orho_w
                                prob = MAX(0.0D0, 1.0D0 - DEXP(-120.0D0*vol*5.2D-4 * Texp))
                                if (massr(n2) .lt. xm0g) then
                                    sumn1 = sumn1 + prob*N_r
                                    sum1 = sum1 + prob*N_r*massr(n2)
                                else
                                    sumn2 = sumn2 + prob*N_r
                                    sum2 = sum2 + prob*N_r*massr(n2)
                                endif
                                if ((sum1+sum2).ge.r_r(i)) EXIT
                            enddo
                            tpi_qrfz(i,j,k,m) = sum1
                            tni_qrfz(i,j,k,m) = sumn1
                            tpg_qrfz(i,j,k,m) = sum2
                            tnr_qrfz(i,j,k,m) = sumn2
                        enddo
                    enddo
                    !$OMP END PARALLEL DO

                    !$OMP PARALLEL DO SCHEDULE(dynamic) num_threads(threads) &
                    !$OMP PRIVATE(j,i,nu_c,lamc,N0_c,sum1,sumn2,vol,prob,N_c)
                    do j = 1, nbc
                        nu_c = MIN(15, NINT(1000.E6/t_Nc(j)) + 2)
                        do i = 1, ntb_c
                            lamc = (t_Nc(j)*am_r* ccg(2,nu_c) * ocg1(nu_c) / r_c(i))**obmr
                            N0_c = t_Nc(j)*ocg1(nu_c) * lamc**cce(1,nu_c)
                            sum1 = 0.0d0
                            sumn2 = 0.0d0
                            do n = nbc, 1, -1
                                vol = massc(n)*orho_w
                                prob = MAX(0.0D0, 1.0D0 - DEXP(-120.0D0*vol*5.2D-4 * Texp))
                                N_c = N0_c*Dc(n)**nu_c*EXP(-lamc*Dc(n))*dtc(n)
                                sumn2 = MIN(t_Nc(j), sumn2 + prob*N_c)
                                sum1 = sum1 + prob*N_c*massc(n)
                                if (sum1 .ge. r_c(i)) EXIT
                            enddo
                            tpi_qcfz(i,j,k,m) = sum1
                            tni_qcfz(i,j,k,m) = sumn2
                        enddo
                    enddo
                    !$OMP END PARALLEL DO
                enddo
            enddo

            IF ( write_thompson_tables ) THEN
                write(0,*) "Writing "//freeze_h2o_file//" in Tempo MP init"
                OPEN(63,file=freeze_h2o_file,form="unformatted",err=9234)
                WRITE(63,err=9234)tpi_qrfz
                WRITE(63,err=9234)tni_qrfz
                WRITE(63,err=9234)tpg_qrfz
                WRITE(63,err=9234)tnr_qrfz
                WRITE(63,err=9234)tpi_qcfz
                WRITE(63,err=9234)tni_qcfz
                CLOSE(63)
                RETURN    ! ----- RETURN
9234            CONTINUE
                write(0,*) "Error writing "//freeze_h2o_file
                return
            ENDIF
        ENDIF

    end subroutine freezeH2O_par

#endif

    !=================================================================================================================
    !..Compute _radiation_ effective radii of cloud water, ice, and snow.
    !.. These are entirely consistent with microphysics assumptions, not
    !.. constant or otherwise ad hoc as is internal to most radiation
    !.. schemes.  Since only the smallest snowflakes should impact
    !.. radiation, compute from first portion of complicated Field number
    !.. distribution, not the second part, which is the larger sizes.

    subroutine calc_effectRad (t1d, p1d, qv1d, qc1d, nc1d, qi1d, ni1d, qs1d,   &
    &                re_qc1d, re_qi1d, re_qs1d, kts, kte, lsml, configs)

        IMPLICIT NONE

        !..Sub arguments
        INTEGER, INTENT(IN):: kts, kte
        REAL, DIMENSION(kts:kte), INTENT(IN)::                            &
        &                    t1d, p1d, qv1d, qc1d, nc1d, qi1d, ni1d, qs1d
        REAL, DIMENSION(kts:kte), INTENT(INOUT):: re_qc1d, re_qi1d, re_qs1d
        type(ty_tempo_cfg), intent(in) :: configs
        integer, intent(in), optional :: lsml

        !..Local variables
        INTEGER:: k
        REAL, DIMENSION(kts:kte):: rho, rc, nc, ri, ni, rs
        REAL:: smo2, smob, smoc
        REAL:: tc0, loga_, a_, b_
        DOUBLE PRECISION:: lamc, lami
        LOGICAL:: has_qc, has_qi, has_qs
        INTEGER:: inu_c
        real, dimension(15), parameter:: g_ratio = (/24,60,120,210,336,   &
        &                504,720,990,1320,1716,2184,2730,3360,4080,4896/)

        has_qc = .false.
        has_qi = .false.
        has_qs = .false.

        re_qc1d = 0.
        re_qi1d = 0.
        re_qs1d = 0.

        do k = kts, kte
            rho(k) = 0.622*p1d(k)/(R*t1d(k)*(qv1d(k)+0.622))
            rc(k) = MAX(R1, qc1d(k)*rho(k))
            nc(k) = MAX(2., MIN(nc1d(k)*rho(k), Nt_c_max))
            if (.not. (configs%aerosol_aware .or. merra2_aerosol_aware)) then
               nc(k) = Nt_c
               if (present(lsml)) then
                  if( lsml == 1) then
                     nc(k) = Nt_c_l
                  else
                     nc(k) = Nt_c_o
                  endif
               endif
            endif
            if (rc(k).gt.R1 .and. nc(k).gt.R2) has_qc = .true.
            ri(k) = MAX(R1, qi1d(k)*rho(k))
            ni(k) = MAX(R2, ni1d(k)*rho(k))
            if (ri(k).gt.R1 .and. ni(k).gt.R2) has_qi = .true.
            rs(k) = MAX(R1, qs1d(k)*rho(k))
            if (rs(k).gt.R1) has_qs = .true.
        enddo

        if (has_qc) then
            do k = kts, kte
                if (rc(k).le.R1 .or. nc(k).le.R2) CYCLE
                if (nc(k).lt.100) then
                    inu_c = 15
                elseif (nc(k).gt.1.E10) then
                    inu_c = 2
                else
                    inu_c = MIN(15, NINT(1000.E6/nc(k)) + 2)
                endif
                lamc = (nc(k)*am_r*g_ratio(inu_c)/rc(k))**obmr
                re_qc1d(k) = MAX(2.51E-6, MIN(SNGL(0.5D0 * DBLE(3.+inu_c)/lamc), 50.E-6))
            enddo
        endif

        if (has_qi) then
            do k = kts, kte
                if (ri(k).le.R1 .or. ni(k).le.R2) CYCLE
                lami = (am_i*cig(2)*oig1*ni(k)/ri(k))**obmi
                re_qi1d(k) = MAX(2.51E-6, MIN(SNGL(0.5D0 * DBLE(3.+mu_i)/lami), 125.E-6))
            enddo
        endif

        if (has_qs) then
            do k = kts, kte
                if (rs(k).le.R1) CYCLE
                tc0 = MIN(-0.1, t1d(k)-273.15)
                smob = rs(k)*oams

                !..All other moments based on reference, 2nd moment.  If bm_s.ne.2,
                !.. then we must compute actual 2nd moment and use as reference.
                if (bm_s.gt.(2.0-1.e-3) .and. bm_s.lt.(2.0+1.e-3)) then
                    smo2 = smob
                else
                    loga_ = sa(1) + sa(2)*tc0 + sa(3)*bm_s &
                    &         + sa(4)*tc0*bm_s + sa(5)*tc0*tc0 &
                    &         + sa(6)*bm_s*bm_s + sa(7)*tc0*tc0*bm_s &
                    &         + sa(8)*tc0*bm_s*bm_s + sa(9)*tc0*tc0*tc0 &
                    &         + sa(10)*bm_s*bm_s*bm_s
                    a_ = 10.0**loga_
                    b_ = sb(1) + sb(2)*tc0 + sb(3)*bm_s &
                    &         + sb(4)*tc0*bm_s + sb(5)*tc0*tc0 &
                    &         + sb(6)*bm_s*bm_s + sb(7)*tc0*tc0*bm_s &
                    &         + sb(8)*tc0*bm_s*bm_s + sb(9)*tc0*tc0*tc0 &
                    &         + sb(10)*bm_s*bm_s*bm_s
                    smo2 = (smob/a_)**(1./b_)
                endif
                !..Calculate bm_s+1 (th) moment.  Useful for diameter calcs.
                loga_ = sa(1) + sa(2)*tc0 + sa(3)*cse(1) &
                &         + sa(4)*tc0*cse(1) + sa(5)*tc0*tc0 &
                &         + sa(6)*cse(1)*cse(1) + sa(7)*tc0*tc0*cse(1) &
                &         + sa(8)*tc0*cse(1)*cse(1) + sa(9)*tc0*tc0*tc0 &
                &         + sa(10)*cse(1)*cse(1)*cse(1)
                a_ = 10.0**loga_
                b_ = sb(1)+ sb(2)*tc0 + sb(3)*cse(1) + sb(4)*tc0*cse(1) &
                &        + sb(5)*tc0*tc0 + sb(6)*cse(1)*cse(1) &
                &        + sb(7)*tc0*tc0*cse(1) + sb(8)*tc0*cse(1)*cse(1) &
                &        + sb(9)*tc0*tc0*tc0 + sb(10)*cse(1)*cse(1)*cse(1)
                smoc = a_ * smo2**b_
                re_qs1d(k) = MAX(5.01E-6, MIN(0.5*(smoc/smob), 999.E-6))
            enddo
        endif

    end subroutine calc_effectRad

    !=================================================================================================================
    !..Compute radar reflectivity assuming 10 cm wavelength radar and using
    !.. Rayleigh approximation.  Only complication is melted snow/graupel
    !.. which we treat as water-coated ice spheres and use Uli Blahak's
    !.. library of routines.  The meltwater fraction is simply the amount
    !.. of frozen species remaining from what initially existed at the
    !.. melting level interface.

    subroutine calc_refl10cm (qv1d, qc1d, qr1d, nr1d, qs1d, qg1d, ng1d, qb1d, &
        t1d, p1d, dBZ, kts, kte, ii, jj, configs, rand1, melti, &
        vt_dBZ, first_time_step)

        IMPLICIT NONE

        !..Sub arguments
        INTEGER, INTENT(IN):: kts, kte, ii, jj
        REAL, OPTIONAL, INTENT(IN):: rand1
        REAL, DIMENSION(kts:kte), INTENT(IN)::                            &
            qv1d, qc1d, qr1d, nr1d, qs1d, qg1d, ng1d, qb1d, t1d, p1d
        REAL, DIMENSION(kts:kte), INTENT(INOUT):: dBZ
        REAL, DIMENSION(kts:kte), OPTIONAL, INTENT(INOUT):: vt_dBZ
        LOGICAL, OPTIONAL, INTENT(IN) :: first_time_step

        type(ty_tempo_cfg), intent(in) :: configs

        !..Local variables
        LOGICAL :: do_vt_dBZ
        LOGICAL :: allow_wet_graupel
        LOGICAL :: allow_wet_snow
        REAL, DIMENSION(kts:kte):: temp, pres, qv, rho, rhof
        REAL, DIMENSION(kts:kte):: rc, rr, nr, rs, rg, ng, rb
        INTEGER, DIMENSION(kts:kte):: idx_bg

        DOUBLE PRECISION, DIMENSION(kts:kte):: ilamr, ilamg, N0_r, N0_g
        REAL, DIMENSION(kts:kte):: mvd_r
        REAL, DIMENSION(kts:kte):: smob, smo2, smoc, smoz
        REAL:: oM3, M0, Mrat, slam1, slam2, xDs
        REAL:: ils1, ils2, t1_vts, t2_vts, t3_vts, t4_vts
        REAL:: vtr_dbz_wt, vts_dbz_wt, vtg_dbz_wt

        REAL, DIMENSION(kts:kte):: ze_rain, ze_snow, ze_graupel

        DOUBLE PRECISION:: N0_exp, N0_min, lam_exp, lamr, lamg
        REAL:: a_, b_, loga_, tc0, SR
        DOUBLE PRECISION:: fmelt_s, fmelt_g

        INTEGER:: i, k, k_0, kbot, n
        LOGICAL, OPTIONAL, INTENT(IN):: melti
        LOGICAL, DIMENSION(kts:kte):: L_qr, L_qs, L_qg

        DOUBLE PRECISION:: cback, x, eta, f_d
        REAL:: xslw1, ygra1, zans1

        !+---+
        if (present(vt_dBZ) .and. present(first_time_step)) then
            do_vt_dBZ = .true.
            if (first_time_step) then
                !           no bright banding, to be consistent with hydrometeor retrieval in GSI
                allow_wet_snow = .false.
            else
                allow_wet_snow = .true.
            endif
            allow_wet_graupel = .false.
        else
            do_vt_dBZ = .false.
            allow_wet_snow = .true.
            allow_wet_graupel = .false.
        endif

        do k = kts, kte
            dBZ(k) = -35.0
        enddo

        !+---+-----------------------------------------------------------------+
        !..Put column of data into local arrays.
        !+---+-----------------------------------------------------------------+
        do k = kts, kte
            temp(k) = t1d(k)
            qv(k) = MAX(1.E-10, qv1d(k))
            pres(k) = p1d(k)
            rho(k) = 0.622*pres(k)/(R*temp(k)*(qv(k)+0.622))
            rhof(k) = SQRT(RHO_NOT/rho(k))
            rc(k) = MAX(R1, qc1d(k)*rho(k))
            if (qr1d(k) .gt. R1) then
                rr(k) = qr1d(k)*rho(k)
                nr(k) = MAX(R2, nr1d(k)*rho(k))
                lamr = (am_r*crg(3)*org2*nr(k)/rr(k))**obmr
                ilamr(k) = 1./lamr
                N0_r(k) = nr(k)*org2*lamr**cre(2)
                mvd_r(k) = (3.0 + mu_r + 0.672) * ilamr(k)
                L_qr(k) = .true.
            else
                rr(k) = R1
                nr(k) = R1
                mvd_r(k) = 50.E-6
                L_qr(k) = .false.
            endif
            if (qs1d(k) .gt. R2) then
                rs(k) = qs1d(k)*rho(k)
                L_qs(k) = .true.
            else
                rs(k) = R1
                L_qs(k) = .false.
            endif
            if (qg1d(k) .gt. R2) then
                rg(k) = qg1d(k)*rho(k)
                ng(k) = MAX(R2, ng1d(k)*rho(k))
                rb(k) = MAX(qg1d(k)/rho_g(NRHG), qb1d(k))
                rb(k) = MIN(qg1d(k)/rho_g(1), rb(k))
                idx_bg(k) = MAX(1,MIN(NINT(qg1d(k)/rb(k) *0.01)+1,NRHG))
                if (.not. configs%hail_aware) idx_bg(k) = idx_bg1
                L_qg(k) = .true.
            else
                rg(k) = R1
                ng(k) = R2
                idx_bg(k) = idx_bg1
                L_qg(k) = .false.
            endif
        enddo

        !+---+-----------------------------------------------------------------+
        !..Calculate y-intercept, slope, and useful moments for snow.
        !+---+-----------------------------------------------------------------+
        do k = kts, kte
            smo2(k) = 0.
            smob(k) = 0.
            smoc(k) = 0.
            smoz(k) = 0.
        enddo
        if (ANY(L_qs .eqv. .true.)) then
            do k = kts, kte
                if (.not. L_qs(k)) CYCLE
                tc0 = MIN(-0.1, temp(k)-273.15)
                smob(k) = rs(k)*oams

                !..All other moments based on reference, 2nd moment.  If bm_s.ne.2,
                !.. then we must compute actual 2nd moment and use as reference.
                if (bm_s.gt.(2.0-1.e-3) .and. bm_s.lt.(2.0+1.e-3)) then
                    smo2(k) = smob(k)
                else
                    loga_ = sa(1) + sa(2)*tc0 + sa(3)*bm_s &
                    &         + sa(4)*tc0*bm_s + sa(5)*tc0*tc0 &
                    &         + sa(6)*bm_s*bm_s + sa(7)*tc0*tc0*bm_s &
                    &         + sa(8)*tc0*bm_s*bm_s + sa(9)*tc0*tc0*tc0 &
                    &         + sa(10)*bm_s*bm_s*bm_s
                    a_ = 10.0**loga_
                    b_ = sb(1) + sb(2)*tc0 + sb(3)*bm_s &
                    &         + sb(4)*tc0*bm_s + sb(5)*tc0*tc0 &
                    &         + sb(6)*bm_s*bm_s + sb(7)*tc0*tc0*bm_s &
                    &         + sb(8)*tc0*bm_s*bm_s + sb(9)*tc0*tc0*tc0 &
                    &         + sb(10)*bm_s*bm_s*bm_s
                    smo2(k) = (smob(k)/a_)**(1./b_)
                endif

                !..Calculate bm_s+1 (th) moment.  Useful for diameter calcs.
                loga_ = sa(1) + sa(2)*tc0 + sa(3)*cse(1) &
                &         + sa(4)*tc0*cse(1) + sa(5)*tc0*tc0 &
                &         + sa(6)*cse(1)*cse(1) + sa(7)*tc0*tc0*cse(1) &
                &         + sa(8)*tc0*cse(1)*cse(1) + sa(9)*tc0*tc0*tc0 &
                &         + sa(10)*cse(1)*cse(1)*cse(1)
                a_ = 10.0**loga_
                b_ = sb(1)+ sb(2)*tc0 + sb(3)*cse(1) + sb(4)*tc0*cse(1) &
                &        + sb(5)*tc0*tc0 + sb(6)*cse(1)*cse(1) &
                &        + sb(7)*tc0*tc0*cse(1) + sb(8)*tc0*cse(1)*cse(1) &
                &        + sb(9)*tc0*tc0*tc0 + sb(10)*cse(1)*cse(1)*cse(1)
                smoc(k) = a_ * smo2(k)**b_

                !..Calculate bm_s*2 (th) moment.  Useful for reflectivity.
                loga_ = sa(1) + sa(2)*tc0 + sa(3)*cse(3) &
                &         + sa(4)*tc0*cse(3) + sa(5)*tc0*tc0 &
                &         + sa(6)*cse(3)*cse(3) + sa(7)*tc0*tc0*cse(3) &
                &         + sa(8)*tc0*cse(3)*cse(3) + sa(9)*tc0*tc0*tc0 &
                &         + sa(10)*cse(3)*cse(3)*cse(3)
                a_ = 10.0**loga_
                b_ = sb(1)+ sb(2)*tc0 + sb(3)*cse(3) + sb(4)*tc0*cse(3) &
                &        + sb(5)*tc0*tc0 + sb(6)*cse(3)*cse(3) &
                &        + sb(7)*tc0*tc0*cse(3) + sb(8)*tc0*cse(3)*cse(3) &
                &        + sb(9)*tc0*tc0*tc0 + sb(10)*cse(3)*cse(3)*cse(3)
                smoz(k) = a_ * smo2(k)**b_
            enddo
        endif

        !+---+-----------------------------------------------------------------+
        !..Calculate y-intercept, slope values for graupel.
        !+---+-----------------------------------------------------------------+
        if (ANY(L_qg .eqv. .true.)) then
            do k = kte, kts, -1
                lamg = (am_g(idx_bg(k))*cgg(3,1)*ogg2*ng(k)/rg(k))**obmg
                ilamg(k) = 1./lamg
                N0_g(k) = ng(k)*ogg2*lamg**cge(2,1)
            enddo
        else
            ilamg(:) = 0.
            N0_g(:) = 0.
        endif

        !+---+-----------------------------------------------------------------+
        !..Locate K-level of start of melting (k_0 is level above).
        !+---+-----------------------------------------------------------------+
        k_0 = kts
        if (present(melti)) then
            if ( melti ) then
                K_LOOP:do k = kte-1, kts, -1
                    if ((temp(k).gt.273.15) .and. L_qr(k)                         &
                    &                            .and. (L_qs(k+1).or.L_qg(k+1)) ) then
                        k_0 = MAX(k+1, k_0)
                        EXIT K_LOOP
                    endif
                enddo K_LOOP
            endif
        endif
        !+---+-----------------------------------------------------------------+
        !..Assume Rayleigh approximation at 10 cm wavelength. Rain (all temps)
        !.. and non-water-coated snow and graupel when below freezing are
        !.. simple. Integrations of m(D)*m(D)*N(D)*dD.
        !+---+-----------------------------------------------------------------+

        do k = kts, kte
            ze_rain(k) = 1.e-22
            ze_snow(k) = 1.e-22
            ze_graupel(k) = 1.e-22
            if (L_qr(k)) ze_rain(k) = N0_r(k)*crg(4)*ilamr(k)**cre(4)
            if (L_qs(k)) ze_snow(k) = (0.176/0.93) * (6.0/PI)*(6.0/PI)     &
            &                           * (am_s/900.0)*(am_s/900.0)*smoz(k)
            if (L_qg(k)) ze_graupel(k) = (0.176/0.93) * (6.0/PI)*(6.0/PI)  &
            &               * (am_g(idx_bg(k))/900.0)*(am_g(idx_bg(k))/900.0)  &
            &               * N0_g(k)*cgg(4,1)*ilamg(k)**cge(4,1)
        enddo

        !+---+-----------------------------------------------------------------+
        !..Special case of melting ice (snow/graupel) particles.  Assume the
        !.. ice is surrounded by the liquid water.  Fraction of meltwater is
        !.. extremely simple based on amount found above the melting level.
        !.. Uses code from Uli Blahak (rayleigh_soak_wetgraupel and supporting
        !.. routines).
        !+---+-----------------------------------------------------------------+
        if (present(melti)) then
            if (.not. iiwarm .and. melti .and. k_0.ge.2) then
                do k = k_0-1, kts, -1

                    !..Reflectivity contributed by melting snow
                    if (allow_wet_snow .and. L_qs(k) .and. L_qs(k_0) ) then
                        SR = MAX(0.01, MIN(1.0 - rs(k)/(rs(k) + rr(k)), 0.99))
                        fmelt_s = DBLE(SR*SR)
                        eta = 0.d0
                        oM3 = 1./smoc(k)
                        M0 = (smob(k)*oM3)
                        Mrat = smob(k)*M0*M0*M0
                        slam1 = M0 * Lam0
                        slam2 = M0 * Lam1
                        do n = 1, nrbins
                            x = am_s * xxDs(n)**bm_s
                            call rayleigh_soak_wetgraupel (x, DBLE(ocms), DBLE(obms), &
                            &              fmelt_s, melt_outside_s, m_w_0, m_i_0, lamda_radar, &
                            &              CBACK, mixingrulestring_s, matrixstring_s,          &
                            &              inclusionstring_s, hoststring_s,                    &
                            &              hostmatrixstring_s, hostinclusionstring_s)
                            f_d = Mrat*(Kap0*DEXP(-slam1*xxDs(n))                     &
                            &              + Kap1*(M0*xxDs(n))**mu_s * DEXP(-slam2*xxDs(n)))
                            eta = eta + f_d * CBACK * simpson(n) * xdts(n)
                        enddo
                        ze_snow(k) = SNGL(lamda4 / (pi5 * K_w) * eta)
                    endif

                    !..Reflectivity contributed by melting graupel
                    if (allow_wet_graupel .and. L_qg(k) .and. L_qg(k_0) ) then
                        SR = MAX(0.01, MIN(1.0 - rg(k)/(rg(k) + rr(k)), 0.99))
                        fmelt_g = DBLE(SR*SR)
                        eta = 0.d0
                        lamg = 1./ilamg(k)
                        do n = 1, nrbins
                            x = am_g(idx_bg(k)) * xxDg(n)**bm_g
                            call rayleigh_soak_wetgraupel (x, DBLE(ocmg(idx_bg(k))), DBLE(obmg), &
                            &              fmelt_g, melt_outside_g, m_w_0, m_i_0, lamda_radar, &
                            &              CBACK, mixingrulestring_g, matrixstring_g,          &
                            &              inclusionstring_g, hoststring_g,                    &
                            &              hostmatrixstring_g, hostinclusionstring_g)
                            f_d = N0_g(k)*xxDg(n)**mu_g * DEXP(-lamg*xxDg(n))
                            eta = eta + f_d * CBACK * simpson(n) * xdtg(n)
                        enddo
                        ze_graupel(k) = SNGL(lamda4 / (pi5 * K_w) * eta)
                    endif

                enddo
            endif
        endif

        do k = kte, kts, -1
            dBZ(k) = 10.*log10((ze_rain(k)+ze_snow(k)+ze_graupel(k))*1.d18)
        enddo

        !..Reflectivity-weighted terminal velocity (snow, rain, graupel, mix).
        if (do_vt_dBZ) then
            do k = kte, kts, -1
                vt_dBZ(k) = 1.E-3
                if (rs(k).gt.R2) then
                    Mrat = smob(k) / smoc(k)
                    ils1 = 1./(Mrat*Lam0 + fv_s)
                    ils2 = 1./(Mrat*Lam1 + fv_s)
                    t1_vts = Kap0*csg(5)*ils1**cse(5)
                    t2_vts = Kap1*Mrat**mu_s*csg(11)*ils2**cse(11)
                    ils1 = 1./(Mrat*Lam0)
                    ils2 = 1./(Mrat*Lam1)
                    t3_vts = Kap0*csg(6)*ils1**cse(6)
                    t4_vts = Kap1*Mrat**mu_s*csg(12)*ils2**cse(12)
                    vts_dbz_wt = rhof(k)*av_s * (t1_vts+t2_vts)/(t3_vts+t4_vts)
                    if (temp(k).ge.273.15 .and. temp(k).lt.275.15) then
                        vts_dbz_wt = vts_dbz_wt*1.5
                    elseif (temp(k).ge.275.15) then
                        vts_dbz_wt = vts_dbz_wt*2.0
                    endif
                else
                    vts_dbz_wt = 1.E-3
                endif

                if (rr(k).gt.R1) then
                    lamr = 1./ilamr(k)
                    vtr_dbz_wt = rhof(k)*av_r*crg(13)*(lamr+fv_r)**(-cre(13))      &
                        / (crg(4)*lamr**(-cre(4)))
                else
                    vtr_dbz_wt = 1.E-3
                endif

                if (rg(k).gt.R2) then
                    lamg = 1./ilamg(k)
                    vtg_dbz_wt = rhof(k)*av_g(idx_bg(k))*cgg(5,idx_bg(k))*lamg**(-cge(5,idx_bg(k)))               &
                    &               / (cgg(4,1)*lamg**(-cge(4,1)))
                else
                    vtg_dbz_wt = 1.E-3
                endif

                ! if (rg(k).gt.R2) then
                !     lamg = 1./ilamg(k)
                !     vtg_dbz_wt = rhof(k)*av_g*cgg(5)*lamg**(-cge(5))               &
                !         / (cgg(4)*lamg**(-cge(4)))
                ! else
                !     vtg_dbz_wt = 1.E-3
                ! endif

                vt_dBZ(k) = (vts_dbz_wt*ze_snow(k) + vtr_dbz_wt*ze_rain(k)      &
                    + vtg_dbz_wt*ze_graupel(k))                        &
                    / (ze_rain(k)+ze_snow(k)+ze_graupel(k))
            enddo
        endif

    end subroutine calc_refl10cm

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

!=================================================================================================================
    !+---+-----------------------------------------------------------------+
    ! THIS FUNCTION CALCULATES THE LIQUID SATURATION VAPOR MIXING RATIO AS
    ! A FUNCTION OF TEMPERATURE AND PRESSURE
    !
    REAL FUNCTION RSLF(P,T)

        IMPLICIT NONE
        REAL, INTENT(IN):: P, T
        REAL:: ESL,X
        REAL, PARAMETER:: C0= .611583699E03
        REAL, PARAMETER:: C1= .444606896E02
        REAL, PARAMETER:: C2= .143177157E01
        REAL, PARAMETER:: C3= .264224321E-1
        REAL, PARAMETER:: C4= .299291081E-3
        REAL, PARAMETER:: C5= .203154182E-5
        REAL, PARAMETER:: C6= .702620698E-8
        REAL, PARAMETER:: C7= .379534310E-11
        REAL, PARAMETER:: C8=-.321582393E-13

        X=MAX(-80.,T-273.16)

        !      ESL=612.2*EXP(17.67*X/(T-29.65))
        ESL=C0+X*(C1+X*(C2+X*(C3+X*(C4+X*(C5+X*(C6+X*(C7+X*C8)))))))
        ESL=MIN(ESL, P*0.15)        ! Even with P=1050mb and T=55C, the sat. vap. pres only contributes to ~15% of total pres.
        RSLF=.622*ESL/(P-ESL)

        !    ALTERNATIVE
        !  ; Source: Murphy and Koop, Review of the vapour pressure of ice and
        !             supercooled water for atmospheric applications, Q. J. R.
        !             Meteorol. Soc (2005), 131, pp. 1539-1565.
        !    ESL = EXP(54.842763 - 6763.22 / T - 4.210 * ALOG(T) + 0.000367 * T
        !        + TANH(0.0415 * (T - 218.8)) * (53.878 - 1331.22
        !        / T - 9.44523 * ALOG(T) + 0.014025 * T))

    END FUNCTION RSLF
    !+---+-----------------------------------------------------------------+
    ! THIS FUNCTION CALCULATES THE ICE SATURATION VAPOR MIXING RATIO AS A
    ! FUNCTION OF TEMPERATURE AND PRESSURE
    !
    REAL FUNCTION RSIF(P,T)

        IMPLICIT NONE
        REAL, INTENT(IN):: P, T
        REAL:: ESI,X
        REAL, PARAMETER:: C0= .609868993E03
        REAL, PARAMETER:: C1= .499320233E02
        REAL, PARAMETER:: C2= .184672631E01
        REAL, PARAMETER:: C3= .402737184E-1
        REAL, PARAMETER:: C4= .565392987E-3
        REAL, PARAMETER:: C5= .521693933E-5
        REAL, PARAMETER:: C6= .307839583E-7
        REAL, PARAMETER:: C7= .105785160E-9
        REAL, PARAMETER:: C8= .161444444E-12

        X=MAX(-80.,T-273.16)
        ESI=C0+X*(C1+X*(C2+X*(C3+X*(C4+X*(C5+X*(C6+X*(C7+X*C8)))))))
        ESI=MIN(ESI, P*0.15)
        RSIF=.622*ESI/max(1.e-4,(P-ESI))

        !    ALTERNATIVE
        !  ; Source: Murphy and Koop, Review of the vapour pressure of ice and
        !             supercooled water for atmospheric applications, Q. J. R.
        !             Meteorol. Soc (2005), 131, pp. 1539-1565.
        !     ESI = EXP(9.550426 - 5723.265/T + 3.53068*ALOG(T) - 0.00728332*T)

    END FUNCTION RSIF
    !=================================================================================================================
    
end module module_mp_tempo_utils
