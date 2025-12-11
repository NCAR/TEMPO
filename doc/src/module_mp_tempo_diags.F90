module module_mp_tempo_diags
  !! utilities for tempo microphysics
  use module_mp_tempo_params ! only : wp, sp, dp
  ! use module_mp_tempo_init, only : tempo_init_cfgs
  use module_mp_radar

  implicit none
  private

  public :: calc_effectRad, calc_refl10cm, hail_size_diagnostics

  contains 

   !=================================================================================================================
    !..Compute _radiation_ effective radii of cloud water, ice, and snow.
    !.. These are entirely consistent with microphysics assumptions, not
    !.. constant or otherwise ad hoc as is internal to most radiation
    !.. schemes.  Since only the smallest snowflakes should impact
    !.. radiation, compute from first portion of complicated Field number
    !.. distribution, not the second part, which is the larger sizes.

    subroutine calc_effectRad (t1d, p1d, qv1d, qc1d, nc1d, qi1d, ni1d, qs1d, re_qc1d, re_qi1d, re_qs1d, kts, kte, lsml)

        IMPLICIT NONE

        !..Sub arguments
        INTEGER, INTENT(IN):: kts, kte
        REAL, DIMENSION(kts:kte), INTENT(IN)::                            &
        &                    t1d, p1d, qv1d, qc1d, nc1d, qi1d, ni1d, qs1d
        REAL, DIMENSION(kts:kte), INTENT(INOUT):: re_qc1d, re_qi1d, re_qs1d
        integer, intent(in), optional :: lsml

        !..Local variables
        INTEGER:: k
        REAL, DIMENSION(kts:kte):: rho, rc, nc, ri, ni, rs
        REAL:: smo2, smob, smoc
        REAL:: tc0, loga_, a_, b_
        DOUBLE PRECISION:: lamc, lami
        LOGICAL:: has_qc, has_qi, has_qs
        INTEGER:: inu_c
        real, dimension(15), parameter:: g_ratio = (/24,60,120,210,336, 504,720,990,1320,1716,2184,2730,3360,4080,4896/)

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
            if (.not. (tempo_init_cfgs%aerosolaware_flag .or. merra2_aerosol_aware)) then
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
        t1d, p1d, dBZ, kts, kte, ii, jj, rand1, melti, &
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

      xam_r = am_r
      xbm_r = bm_r
      xmu_r = mu_r
      xam_s = am_s
      xbm_s = bm_s
      xmu_s = mu_s
      xam_g = am_g(idx_bg1)
      xbm_g = bm_g
      xmu_g = mu_g

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
                if (.not. tempo_init_cfgs%hailaware_flag) idx_bg(k) = idx_bg1
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

                    ! ..Reflectivity contributed by melting snow
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
    !..Compute max hail size aloft and at the ground (both 2D fields)

    subroutine hail_size_diagnostics(kts, kte, qg1d, ng1d, qb1d, t1d, p1d, qv1d, qg_max_diam1d)

      implicit none

      integer, intent(in) :: kts, kte
      real(wp), dimension(kts:kte), intent(in) :: qg1d, ng1d, qb1d, t1d, p1d, qv1d
      real(wp), dimension(kts:kte), intent(out) :: qg_max_diam1d

      ! local variables
      real(wp), dimension(kts:kte) :: rho, rg, ng, rb
      integer, dimension(kts:kte) :: idx_bg
      real(dp) :: lamg, N0_g, f_d, sum_nh, sum_t, hail_max
      integer :: k, n
      integer, parameter :: nhbins = 50
      real(dp), dimension(nhbins):: hbins, dhbins
      real(dp), parameter :: lowbin = 500.e-6
      real(dp), parameter :: highbin = 0.075
      real(dp), parameter :: threshold_conc = 0.0005

      ! Binned number distribution method
      call create_bins(numbins=nhbins, lowbin=lowbin, highbin=highbin, bins=hbins, deltabins=dhbins)

      do k = kts, kte
         qg_max_diam1d(k) = 0.
         if(qg1d(k) >= 1.e-6) then
            rho(k) = 0.622*p1d(k)/(R*t1d(k)*(qv1d(k)+0.622))
            rg(k) = qg1d(k)*rho(k)
            ng(k) = max(R2, ng1d(k)*rho(k))
            rb(k) = max(qg1d(k)/rho_g(nrhg), qb1d(k))
            rb(k) = min(qg1d(k)/rho_g(1), rb(k))
            idx_bg(k) = max(1,min(nint(qg1d(k)/rb(k) *0.01)+1,nrhg))
            if (.not. tempo_init_cfgs%hailaware_flag) idx_bg(k) = idx_bg1
            if (rho_g(idx_bg(k)) < 350.) cycle

            lamg = (am_g(idx_bg(k))*cgg(3,1)*ogg2*ng(k)/rg(k))**obmg
            N0_g = ng(k)*ogg2*lamg**cge(2,1)

            sum_nh = 0.
            sum_t = 0.
            do n = nhbins, 1, -1
               f_d = N0_g*hbins(n)**mu_g * exp(-lamg*hbins(n)) * dhbins(n)
               sum_nh = sum_nh + f_d
               if (sum_nh > threshold_conc) exit
               sum_t = sum_nh
            enddo

            if (n >= nhbins) then
               hail_max = hbins(nhbins)
            elseif (hbins(n+1) .gt. 1.e-3) then
               hail_max = hbins(n) - (sum_nh-threshold_conc)/(sum_nh-sum_t) * (hbins(n)-hbins(n+1))
            else
               hail_max = 1.e-4
            endif
            qg_max_diam1d(k) = 1000. * hail_max ! convert to mm
         endif
      enddo

    end subroutine hail_size_diagnostics

end module module_mp_tempo_diags