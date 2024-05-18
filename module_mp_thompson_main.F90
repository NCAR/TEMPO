! 1D Thompson microphysics scheme
!=================================================================================================================
module module_mp_thompson_main

    use module_mp_thompson_params
    use module_mp_thompson_utils, only: rslf, rsif
    
#if defined(mpas)
    use mpas_kind_types, only: wp => RKIND, sp => R4KIND, dp => R8KIND
#else
    use machine, only: wp => kind_phys, sp => kind_sngl_prec, dp => kind_dbl_prec
#endif

    implicit none

contains
    !=================================================================================================================
    ! This subroutine computes the moisture tendencies of water vapor, cloud droplets, rain, cloud ice (pristine),
    ! snow, and graupel. Previously this code was based on Reisner et al (1998), but few of those pieces remain.
    ! A complete description is now found in Thompson et al. (2004, 2008), Thompson and Eidhammer (2014),
    ! and Jensen et al. (2023).

    subroutine mp_thompson_main(qv1d, qc1d, qi1d, qr1d, qs1d, qg1d, qb1d, ni1d, nr1d, nc1d, ng1d, &
        nwfa1d, nifa1d, t1d, p1d, w1d, dzq, pptrain, pptsnow, pptgraul, pptice, &
#if defined(mpas)
        rainprod, evapprod, &
#endif

#if !defined(mpas)
    ! Extended diagnostics, most arrays only
    ! allocated if ext_diag flag is .true.
        ext_diag, sedi_semi, &
        prw_vcdc1, prw_vcde1, &
        tpri_inu1, tpri_ide1_d, tpri_ide1_s, tprs_ide1, &
        tprs_sde1_d, tprs_sde1_s, &
        tprg_gde1_d, tprg_gde1_s, tpri_iha1, tpri_wfz1, &
        tpri_rfz1, tprg_rfz1, tprs_scw1, tprg_scw1,&
        tprg_rcs1, tprs_rcs1, tprr_rci1, &
        tprg_rcg1, tprw_vcd1_c, &
        tprw_vcd1_e, tprr_sml1, tprr_gml1, tprr_rcg1, &
        tprr_rcs1, tprv_rev1, &
        tten1, qvten1, qrten1, qsten1, &
        qgten1, qiten1, niten1, nrten1, ncten1, qcten1, &
#endif
        decfl, pfil1, pfll1, &
        lsml, rand1, rand2, rand3, &
        kts, kte, dt, ii, jj, &
        configs)

#if !defined(mpas) && defined (MPI)
        use mpi_f08
#endif
        implicit none

        ! Subroutine arguments
        integer, intent(in) :: kts, kte, ii, jj
        real(wp), dimension(kts:kte), intent(inout) :: qv1d, qc1d, qi1d, qr1d, qs1d, qg1d, qb1d, &
            ni1d, nr1d, nc1d, ng1d, nwfa1d, nifa1d, t1d
        real(wp), dimension(kts:kte), intent(in) :: p1d, w1d, dzq
        real(wp), intent(inout) :: pptrain, pptsnow, pptgraul, pptice
        real(wp), intent(in) :: dt

        type(config_flags), intent(in) :: configs

        integer, intent(in), optional :: lsml
        real(wp), intent(in), optional :: rand1, rand2, rand3
        real(wp), dimension(kts:kte), intent(out), optional :: pfil1, pfll1
        integer, intent(in), optional :: decfl
#if !defined(mpas)
        ! Extended diagnostics, most arrays only allocated if ext_diag is true
        logical, intent(in) :: ext_diag
        logical, intent(in) :: sedi_semi
        real(wp), dimension(:), intent(out) :: &
            prw_vcdc1, &
            prw_vcde1, tpri_inu1, tpri_ide1_d, &
            tpri_ide1_s, tprs_ide1, &
            tprs_sde1_d, tprs_sde1_s, tprg_gde1_d, &
            tprg_gde1_s, tpri_iha1, tpri_wfz1, &
            tpri_rfz1, tprg_rfz1, tprs_scw1, tprg_scw1, &
            tprg_rcs1, tprs_rcs1, &
            tprr_rci1, tprg_rcg1, &
            tprw_vcd1_c, tprw_vcd1_e, tprr_sml1, &
            tprr_gml1, tprr_rcg1, &
            tprr_rcs1, tprv_rev1, tten1, qvten1, &
            qrten1, qsten1, qgten1, qiten1, niten1, &
            nrten1, ncten1, qcten1
#endif

#if defined(mpas)
        real(wp), dimension(kts:kte), intent(inout) :: rainprod, evapprod
#endif

        !=================================================================================================================
        ! Local variables
        real(wp), dimension(kts:kte) :: tten, qvten, qcten, qiten, qrten, &
            qsten, qgten, qbten, niten, nrten, ncten, ngten, nwfaten, nifaten

        real(dp), dimension(kts:kte) :: prw_vcd

        real(dp), dimension(kts:kte) :: pnc_wcd, pnc_wau, pnc_rcw, pnc_scw, pnc_gcw

        real(dp), dimension(kts:kte) :: pna_rca, pna_sca, pna_gca, pnd_rcd, pnd_scd, pnd_gcd

        real(dp), dimension(kts:kte) :: prr_wau, prr_rcw, prr_rcs, &
            prr_rcg, prr_sml, prr_gml, &
            prr_rci, prv_rev,          &
            pnr_wau, pnr_rcs, pnr_rcg, &
            pnr_rci, pnr_sml, pnr_gml, &
            pnr_rev, pnr_rcr, pnr_rfz

        real(dp), dimension(kts:kte) :: pri_inu, pni_inu, pri_ihm, &
            pni_ihm, pri_wfz, pni_wfz, &
            pri_rfz, pni_rfz, pri_ide, &
            pni_ide, pri_rci, pni_rci, &
            pni_sci, pni_iau, pri_iha, pni_iha

        real(dp), dimension(kts:kte) :: prs_iau, prs_sci, prs_rcs, prs_scw, prs_sde, prs_ihm, prs_ide

        real(dp), dimension(kts:kte) :: prg_scw, prg_rfz, prg_gde, &
            prg_gcw, prg_rci, prg_rcs, prg_rcg, prg_ihm, &
            png_rcs, png_rcg, png_scw, png_gde, &
            pbg_scw, pbg_rfz, pbg_gcw, pbg_rci, pbg_rcs, pbg_rcg, &
            pbg_sml, pbg_gml

        real(dp), parameter :: zerod0 = 0.0d0

        real(wp), dimension(kts:kte) :: pfll, pfil, pdummy
        real(wp) :: dtcfl, rainsfc, graulsfc, orhodt
        integer :: niter
        real(wp), dimension(kts:kte) :: rr_tmp, nr_tmp, rg_tmp

        real(wp), dimension(kts:kte) :: temp, twet, pres, qv
        real(wp), dimension(kts:kte) :: rc, ri, rr, rs, rg, rb
        real(wp), dimension(kts:kte) :: ni, nr, nc, ns, ng, nwfa, nifa
        real(wp), dimension(kts:kte) :: rho, rhof, rhof2
        real(wp), dimension(kts:kte) :: qvs, qvsi, delqvs
        real(wp), dimension(kts:kte) :: satw, sati, ssatw, ssati
        real(wp), dimension(kts:kte) :: diffu, visco, vsc2, &
            tcond, lvap, ocp, lvt2

        real(dp), dimension(kts:kte) :: ilamr, ilamg, n0_r, n0_g
        real(dp) :: n0_melt
        real(wp), dimension(kts:kte) :: mvd_r, mvd_c, mvd_g
        real(wp), dimension(kts:kte) :: smob, smo2, smo1, smo0, &
            smoc, smod, smoe, smof, smog

        real(wp), dimension(kts:kte) :: sed_r, sed_s, sed_g, sed_i, sed_n, sed_c, sed_b

        real(wp) :: rgvm, delta_tp, orho, lfus2
        real(wp), dimension(5):: onstep
        real(dp) :: n0_exp, n0_min, lam_exp, lamc, lamr, lamg
        real(dp) :: lami, ilami, ilamc
        real(wp) :: xdc, dc_b, dc_g, xdi, xdr, xds, xdg, ds_m, dg_m
        real(dp) :: dr_star, dc_star
        real(wp) :: zeta1, zeta, taud, tau
        real(wp) :: stoke_r, stoke_s, stoke_g, stoke_i
        real(wp) :: vti, vtr, vts, vtg, vtc
        real(wp) :: xrho_g, afall, vtg1, vtg2
        real(wp) :: bfall = 3*b_coeff - 1
        real(wp), dimension(kts:kte+1) :: vtik, vtnik, vtrk, vtnrk, vtsk, vtgk, vtngk, vtck, vtnck
        real(wp), dimension(kts:kte) :: vts_boost
        real(wp) :: m0, slam1, slam2
        real(wp) :: mrat, ils1, ils2, t1_vts, t2_vts, t3_vts, t4_vts, c_snow
        real(wp) :: a_, b_, loga_, a1, a2, tf
        real(wp) :: tempc, tc0, r_mvd1, r_mvd2, xkrat
        real(wp) :: dew_t, tlcl, the
        real(wp) :: xnc, xri, xni, xmi, oxmi, xrc, xrr, xnr, xrg, xng, xrb
        real(wp) :: xsat, rate_max, sump, ratio
        real(wp) :: clap, fcd, dfcd
        real(wp) :: otemp, rvs, rvs_p, rvs_pp, gamsc, alphsc, t1_evap, t1_subl
        real(wp) :: r_frac, g_frac, const_Ri, rime_dens
        real(wp) :: Ef_rw, Ef_sw, Ef_gw, Ef_rr
        real(wp) :: Ef_ra, Ef_sa, Ef_ga
        real(wp) :: dtsave, odts, odt, odzq, hgt_agl, SR
        real(wp) :: xslw1, ygra1, zans1, eva_factor
        real(wp) :: melt_f, rand
        integer :: i, k, k2, n, nn, nstep, k_0, kbot, IT, iexfrq, k_melting
        integer, dimension(5) :: ksed1
        integer :: nir, nis, nig, nii, nic, niin
        integer :: idx_tc, idx_t, idx_s, idx_g1, idx_g, idx_r1, idx_r,     &
            idx_i1, idx_i, idx_c, idx, idx_d, idx_n, idx_in
        integer, dimension(kts:kte) :: idx_bg, idx_table

        logical :: melti, no_micro
        logical, dimension(kts:kte) :: l_qc, l_qi, l_qr, l_qs, l_qg
        logical :: debug_flag
        character*256 :: mp_debug
        integer :: nu_c, decfl_

        !=================================================================================================================

        debug_flag = .false.

        no_micro = .true.
        dtsave = dt
        odt = 1./dt
        odts = 1./dtsave
        iexfrq = 1
        rand = 0.0
        decfl_ = 10
        if (present(decfl)) decfl_ = decfl
        
#if !defined(mpas)
        ! Transition value of coefficient matching at crossover from cloud ice to snow
        av_i = av_s * D0s ** (bv_s - bv_i)
#endif

        !=================================================================================================================
        ! Source/sink terms.  First 2 chars: "pr" represents source/sink of
        ! mass while "pn" represents source/sink of number.  Next char is one
        ! of "v" for water vapor, "r" for rain, "i" for cloud ice, "w" for
        ! cloud water, "s" for snow, and "g" for graupel.  Next chars
        ! represent processes: "de" for sublimation/deposition, "ev" for
        ! evaporation, "fz" for freezing, "ml" for melting, "au" for
        ! autoconversion, "nu" for ice nucleation, "hm" for Hallet/Mossop
        ! secondary ice production, and "c" for collection followed by the
        ! character for the species being collected.  ALL of these terms are
        ! positive (except for deposition/sublimation terms which can switch
        ! signs based on super/subsaturation) and are treated as negatives
        ! where necessary in the tendency equations.
        !=================================================================================================================
        ! TODO: Put these in derived data type and add initialization subroutine
        do k = kts, kte
            tten(k) = 0.
            qvten(k) = 0.
            qcten(k) = 0.
            qiten(k) = 0.
            qrten(k) = 0.
            qsten(k) = 0.
            qgten(k) = 0.
            ngten(k) = 0.
            qbten(k) = 0.
            niten(k) = 0.
            nrten(k) = 0.
            ncten(k) = 0.
            nwfaten(k) = 0.
            nifaten(k) = 0.

            prw_vcd(k) = 0.

            pnc_wcd(k) = 0.
            pnc_wau(k) = 0.
            pnc_rcw(k) = 0.
            pnc_scw(k) = 0.
            pnc_gcw(k) = 0.

            prv_rev(k) = 0.
            prr_wau(k) = 0.
            prr_rcw(k) = 0.
            prr_rcs(k) = 0.
            prr_rcg(k) = 0.
            prr_sml(k) = 0.
            prr_gml(k) = 0.
            prr_rci(k) = 0.
            pnr_wau(k) = 0.
            pnr_rcs(k) = 0.
            pnr_rcg(k) = 0.
            pnr_rci(k) = 0.
            pnr_sml(k) = 0.
            pnr_gml(k) = 0.
            pnr_rev(k) = 0.
            pnr_rcr(k) = 0.
            pnr_rfz(k) = 0.

            pri_inu(k) = 0.
            pni_inu(k) = 0.
            pri_ihm(k) = 0.
            pni_ihm(k) = 0.
            pri_wfz(k) = 0.
            pni_wfz(k) = 0.
            pri_rfz(k) = 0.
            pni_rfz(k) = 0.
            pri_ide(k) = 0.
            pni_ide(k) = 0.
            pri_rci(k) = 0.
            pni_rci(k) = 0.
            pni_sci(k) = 0.
            pni_iau(k) = 0.
            pri_iha(k) = 0.
            pni_iha(k) = 0.

            prs_iau(k) = 0.
            prs_sci(k) = 0.
            prs_rcs(k) = 0.
            prs_scw(k) = 0.
            prs_sde(k) = 0.
            prs_ihm(k) = 0.
            prs_ide(k) = 0.

            prg_scw(k) = 0.
            prg_rfz(k) = 0.
            prg_gde(k) = 0.
            prg_gcw(k) = 0.
            prg_rci(k) = 0.
            prg_rcs(k) = 0.
            prg_rcg(k) = 0.
            prg_ihm(k) = 0.
            !   new source/sink terms for 3-moment graupel
            png_scw(k) = 0.
            png_rcs(k) = 0.
            png_rcg(k) = 0.
            png_gde(k) = 0.

            pbg_scw(k) = 0.
            pbg_rfz(k) = 0.
            pbg_gcw(k) = 0.
            pbg_rci(k) = 0.
            pbg_rcs(k) = 0.
            pbg_rcg(k) = 0.
            pbg_sml(k) = 0.
            pbg_gml(k) = 0.

            pna_rca(k) = 0.
            pna_sca(k) = 0.
            pna_gca(k) = 0.

            pnd_rcd(k) = 0.
            pnd_scd(k) = 0.
            pnd_gcd(k) = 0.

            if (present(pfil1)) pfil1(k) = 0.
            if (present(pfll1)) pfll1(k) = 0.
            pfil(k) = 0.
            pfll(k) = 0.
            pdummy(k) = 0.
        enddo
#if defined(mpas)
        do k = kts, kte
            rainprod(k) = 0.
            evapprod(k) = 0.
        enddo
#endif

#if !defined(mpas)
        !Diagnostics
        if (ext_diag) then
            do k = kts, kte
                !vtsk1(k) = 0.
                !txrc1(k) = 0.
                !txri1(k) = 0.
                prw_vcdc1(k) = 0.
                prw_vcde1(k) = 0.
                tpri_inu1(k) = 0.
                tpri_ide1_d(k) = 0.
                tpri_ide1_s(k) = 0.
                tprs_ide1(k) = 0.
                tprs_sde1_d(k) = 0.
                tprs_sde1_s(k) = 0.
                tprg_gde1_d(k) = 0.
                tprg_gde1_s(k) = 0.
                tpri_iha1(k) = 0.
                tpri_wfz1(k) = 0.
                tpri_rfz1(k) = 0.
                tprg_rfz1(k) = 0.
                tprg_scw1(k) = 0.
                tprs_scw1(k) = 0.
                tprg_rcs1(k) = 0.
                tprs_rcs1(k) = 0.
                tprr_rci1(k) = 0.
                tprg_rcg1(k) = 0.
                tprw_vcd1_c(k) = 0.
                tprw_vcd1_e(k) = 0.
                tprr_sml1(k) = 0.
                tprr_gml1(k) = 0.
                tprr_rcg1(k) = 0.
                tprr_rcs1(k) = 0.
                tprv_rev1(k) = 0.
                tten1(k) = 0.
                qvten1(k) = 0.
                qrten1(k) = 0.
                qsten1(k) = 0.
                qgten1(k) = 0.
                qiten1(k) = 0.
                niten1(k) = 0.
                nrten1(k) = 0.
                ncten1(k) = 0.
                qcten1(k) = 0.
            enddo
        endif
#endif
        !..Bug fix (2016Jun15), prevent use of uninitialized value(s) of snow moments.
        do k = kts, kte
            smo0(k) = 0.
            smo1(k) = 0.
            smo2(k) = 0.
            smob(k) = 0.
            smoc(k) = 0.
            smod(k) = 0.
            smoe(k) = 0.
            smof(k) = 0.
            smog(k) = 0.
            ns(k)   = 0.
            mvd_r(k) = 0.
            mvd_c(k) = 0.
        enddo

        !=================================================================================================================
        ! Convert microphysics variables to concentrations (kg / m^3 and # / m^3)
        do k = kts, kte
            temp(k) = t1d(k)
            qv(k) = max(min_qv, qv1d(k))
            pres(k) = p1d(k)
            rho(k) = RoverRv*pres(k)/(r*temp(k)*(qv(k)+RoverRv))

            ! CCPP version has rho(k) multiplier for min and max
            ! nwfa(k) = max(11.1e6, min(9999.e6, nwfa1d(k)*rho(k)))
            ! nifa(k) = max(nain1*0.01, min(9999.e6, nifa1d(k)*rho(k)))
            nwfa(k) = max(nwfa_default*rho(k), min(aero_max*rho(k), nwfa1d(k)*rho(k)))
            nifa(k) = max(nifa_default*rho(k), min(aero_max*rho(k), nifa1d(k)*rho(k)))

            ! From CCPP version
            mvd_r(k) = D0r
            mvd_c(k) = D0c

            if (qc1d(k) .gt. R1) then
                no_micro = .false.
                rc(k) = qc1d(k)*rho(k)
                nc(k) = max(2., min(nc1d(k)*rho(k), nt_c_max))
                l_qc(k) = .true.
                if (nc(k).gt.10000.e6) then
                    nu_c = 2
                elseif (nc(k).lt.100.) then
                    nu_c = 15
                else
                    nu_c = nint(nu_c_scale/nc(k)) + 2
                    rand = 0.0
                    if (present(rand2)) then
                        rand = rand2
                    endif
                    nu_c = max(2, min(nu_c+nint(rand), 15))
                endif
                lamc = (nc(k)*am_r*ccg(2,nu_c)*ocg1(nu_c)/rc(k))**obmr
                xDc = (bm_r + nu_c + 1.) / lamc
                if (xDc.lt. D0c) then
                    lamc = cce(2,nu_c)/D0c
                elseif (xDc.gt. D0r*2.) then
                    lamc = cce(2,nu_c)/(D0r*2.)
                endif
                nc(k) = min(real(nt_c_max, kind=dp), ccg(1,nu_c)*ocg2(nu_c)*rc(k) / am_r*lamc**bm_r)
                ! CCPP version has different values of Nt_c for land/ocean
                if (.not.(configs%aerosol_aware .or. merra2_aerosol_aware)) then
                    nc(k) = Nt_c
                    if (present(lsml)) then
                        if (lsml == 1) then
                            nc(k) = Nt_c_l
                        else
                            nc(k) = Nt_c_o
                        endif
                    endif
                endif
            else
                qc1d(k) = 0.0
                nc1d(k) = 0.0
                rc(k) = R1
                nc(k) = 2.
                L_qc(k) = .false.
            endif

            if (qi1d(k) .gt. R1) then
                no_micro = .false.
                ri(k) = qi1d(k)*rho(k)
                ni(k) = max(r2, ni1d(k)*rho(k))
                if (ni(k).le. r2) then
                    lami = cie(2)/5.e-6
                    ni(k) = min(max_ni, cig(1)*oig2*ri(k)/am_i*lami**bm_i)
                endif
                L_qi(k) = .true.
                lami = (am_i*cig(2)*oig1*ni(k)/ri(k))**obmi
                ilami = 1./lami
                xDi = (bm_i + mu_i + 1.) * ilami
                if (xDi.lt. 5.E-6) then
                    lami = cie(2)/5.E-6
                    ni(k) = min(max_ni, cig(1)*oig2*ri(k)/am_i*lami**bm_i)
                elseif (xDi.gt. 300.E-6) then
                    lami = cie(2)/300.E-6
                    ni(k) = cig(1)*oig2*ri(k)/am_i*lami**bm_i
                endif
            else
                qi1d(k) = 0.0
                ni1d(k) = 0.0
                ri(k) = R1
                ni(k) = R2
                L_qi(k) = .false.
            endif

            if (qr1d(k) .gt. R1) then
                no_micro = .false.
                rr(k) = qr1d(k)*rho(k)
                nr(k) = max(R2, nr1d(k)*rho(k))
                if (nr(k).le. R2) then
                    mvd_r(k) = 1.0E-3
                    lamr = (3.0 + mu_r + 0.672) / mvd_r(k)
                    nr(k) = crg(2)*org3*rr(k)*lamr**bm_r / am_r
                endif
                L_qr(k) = .true.
                lamr = (am_r*crg(3)*org2*nr(k)/rr(k))**obmr
                mvd_r(k) = (3.0 + mu_r + 0.672) / lamr
                if (mvd_r(k) .gt. 2.5E-3) then
                    mvd_r(k) = 2.5E-3
                    lamr = (3.0 + mu_r + 0.672) / mvd_r(k)
                    nr(k) = crg(2)*org3*rr(k)*lamr**bm_r / am_r
                elseif (mvd_r(k) .lt. D0r*0.75) then
                    mvd_r(k) = D0r*0.75
                    lamr = (3.0 + mu_r + 0.672) / mvd_r(k)
                    nr(k) = crg(2)*org3*rr(k)*lamr**bm_r / am_r
                endif
            else
                qr1d(k) = 0.0
                nr1d(k) = 0.0
                rr(k) = R1
                nr(k) = R2
                L_qr(k) = .false.
            endif
            if (qs1d(k) .gt. R1) then
                no_micro = .false.
                rs(k) = qs1d(k)*rho(k)
                L_qs(k) = .true.
            else
                qs1d(k) = 0.0
                rs(k) = R1
                L_qs(k) = .false.
            endif
            if (qg1d(k) .gt. R1) then
                no_micro = .false.
                L_qg(k) = .true.
                rg(k) = qg1d(k)*rho(k)
                ng(k) = max(r2, ng1d(k)*rho(k))
                rb(k) = max(qg1d(k)/rho_g(nrhg), qb1d(k))
                rb(k) = min(qg1d(k)/rho_g(1), rb(k))
                qb1d(k) = rb(k)
                idx_bg(k) = max(1,min(nint(qg1d(k)/rb(k) *0.01)+1,nrhg))
                idx_table(k) = idx_bg(k)
                if (ng(k).le. R2) then
                    mvd_g(k) = 1.5E-3
                    lamg = (3.0 + mu_g + 0.672) / mvd_g(k)
                    ng(k) = cgg(2,1)*ogg3*rg(k)*lamg**bm_g / am_g(idx_bg(k))
                endif
                lamg = (am_g(idx_bg(k))*cgg(3,1)*ogg2*ng(k)/rg(k))**obmg
                mvd_g(k) = (3.0 + mu_g + 0.672) / lamg
                if (mvd_g(k) .gt. 25.4E-3) then
                    mvd_g(k) = 25.4E-3
                    lamg = (3.0 + mu_g + 0.672) / mvd_g(k)
                    ng(k) = cgg(2,1)*ogg3*rg(k)*lamg**bm_g / am_g(idx_bg(k))
                elseif (mvd_g(k) .lt. D0r) then
                    mvd_g(k) = D0r
                    lamg = (3.0 + mu_g + 0.672) / mvd_g(k)
                    ng(k) = cgg(2,1)*ogg3*rg(k)*lamg**bm_g / am_g(idx_bg(k))
                endif
            else
                qg1d(k) = 0.0
                ng1d(k) = 0.0
                qb1d(k) = 0.0
                idx_bg(k) = idx_bg1
                idx_table(k) = idx_bg(k)
                rg(k) = R1
                ng(k) = R2
                rb(k) = R1/rho(k)/rho_g(NRHG)
                L_qg(k) = .false.
            endif
            if (.not. configs%hail_aware) then
                idx_bg(k) = idx_bg1
                idx_table(k) = idx_bg(k)
                ! If dimNRHG = 1, set idx_table(k) = 1,
                ! otherwise idx_bg1
                if(.not. using_hail_aware_table) then
                    idx_table(k) = 1
                endif
            endif
        enddo

        !     if (debug_flag) then
        !      write(mp_debug,*) 'DEBUG-VERBOSE at (i,j) ', ii, ', ', jj
        !      CALL wrf_debug(550, mp_debug)
        !      do k = kts, kte
        !        write(mp_debug, '(a,i3,f8.2,1x,f7.2,1x, 11(1x,e13.6))')        &
        !    &              'VERBOSE: ', k, pres(k)*0.01, temp(k)-273.15, qv(k), rc(k), rr(k), ri(k), rs(k), rg(k), nc(k), nr(k), ni(k), nwfa(k), nifa(k)
        !        CALL wrf_debug(550, mp_debug)
        !      enddo
        !     endif

        !=================================================================================================================
        ! Derive various thermodynamic variables frequently used.
        ! Saturation vapor pressure (mixing ratio) over liquid/ice comes from
        ! Flatau et al. 1992; enthalpy (latent heat) of vaporization from
        ! Bohren & Albrecht 1998; others from Pruppacher & Klett 1978.
        do k = kts, kte
            tempc = temp(k) - 273.15
            rhof(k) = sqrt(rho_not/rho(k))
            rhof2(k) = sqrt(rhof(k))
            qvs(k) = rslf(pres(k), temp(k))
            delqvs(k) = max(0.0, rslf(pres(k), 273.15)-qv(k))
            if (tempc .le. 0.0) then
                qvsi(k) = rsif(pres(k), temp(k))
            else
                qvsi(k) = qvs(k)
            endif
            satw(k) = qv(k)/qvs(k)
            sati(k) = qv(k)/qvsi(k)
            ssatw(k) = satw(k) - 1.
            ssati(k) = sati(k) - 1.
            if (abs(ssatw(k)).lt. eps) ssatw(k) = 0.0
            if (abs(ssati(k)).lt. eps) ssati(k) = 0.0
            if (no_micro .and. ssati(k).gt. 0.0) no_micro = .false.
            diffu(k) = 2.11e-5*(temp(k)/273.15)**1.94 * (101325./pres(k))
            if (tempc .ge. 0.0) then
                visco(k) = (1.718+0.0049*tempc)*1.0e-5
            else
                visco(k) = (1.718+0.0049*tempc-1.2e-5*tempc*tempc)*1.0e-5
            endif
            ocp(k) = 1./(cp2*(1.+0.887*qv(k)))
            vsc2(k) = sqrt(rho(k)/visco(k))
            lvap(k) = lvap0 + (2106.0 - 4218.0)*tempc
            tcond(k) = (5.69 + 0.0168*tempc)*1.0e-5 * 418.936
        enddo

        !=================================================================================================================
        !..If no existing hydrometeor species and no chance to initiate ice or
        !.. condense cloud water, just exit quickly!
        if (no_micro) return

        !..Calculate y-intercept, slope, and useful moments for snow.
        if (.not. iiwarm) then
            do k = kts, kte
                if (.not. L_qs(k)) CYCLE
                tc0 = min(-0.1, temp(k)-273.15)
                smob(k) = rs(k)*oams

                !..All other moments based on reference, 2nd moment.  If bm_s.ne.2,
                !.. then we must compute actual 2nd moment and use as reference.
                if (bm_s.gt.(2.0-1.e-3) .and. bm_s.lt.(2.0+1.e-3)) then
                    smo2(k) = smob(k)
                else
                    loga_ = sa(1) + sa(2)*tc0 + sa(3)*bm_s &
                        + sa(4)*tc0*bm_s + sa(5)*tc0*tc0 &
                        + sa(6)*bm_s*bm_s + sa(7)*tc0*tc0*bm_s &
                        + sa(8)*tc0*bm_s*bm_s + sa(9)*tc0*tc0*tc0 &
                        + sa(10)*bm_s*bm_s*bm_s
                    a_ = 10.0**loga_
                    b_ = sb(1) + sb(2)*tc0 + sb(3)*bm_s &
                        + sb(4)*tc0*bm_s + sb(5)*tc0*tc0 &
                        + sb(6)*bm_s*bm_s + sb(7)*tc0*tc0*bm_s &
                        + sb(8)*tc0*bm_s*bm_s + sb(9)*tc0*tc0*tc0 &
                        + sb(10)*bm_s*bm_s*bm_s
                    smo2(k) = (smob(k)/a_)**(1./b_)
                endif

                !..Calculate 0th moment.  Represents snow number concentration.
                loga_ = sa(1) + sa(2)*tc0 + sa(5)*tc0*tc0 + sa(9)*tc0*tc0*tc0
                a_ = 10.0**loga_
                b_ = sb(1) + sb(2)*tc0 + sb(5)*tc0*tc0 + sb(9)*tc0*tc0*tc0
                smo0(k) = a_ * smo2(k)**b_

                !..Calculate 1st moment.  Useful for depositional growth and melting.
                loga_ = sa(1) + sa(2)*tc0 + sa(3) &
                    + sa(4)*tc0 + sa(5)*tc0*tc0 &
                    + sa(6) + sa(7)*tc0*tc0 &
                    + sa(8)*tc0 + sa(9)*tc0*tc0*tc0 &
                    + sa(10)
                a_ = 10.0**loga_
                b_ = sb(1)+ sb(2)*tc0 + sb(3) + sb(4)*tc0 &
                    + sb(5)*tc0*tc0 + sb(6) &
                    + sb(7)*tc0*tc0 + sb(8)*tc0 &
                    + sb(9)*tc0*tc0*tc0 + sb(10)
                smo1(k) = a_ * smo2(k)**b_

                !..Calculate bm_s+1 (th) moment.  Useful for diameter calcs.
                loga_ = sa(1) + sa(2)*tc0 + sa(3)*cse(1) &
                    + sa(4)*tc0*cse(1) + sa(5)*tc0*tc0 &
                    + sa(6)*cse(1)*cse(1) + sa(7)*tc0*tc0*cse(1) &
                    + sa(8)*tc0*cse(1)*cse(1) + sa(9)*tc0*tc0*tc0 &
                    + sa(10)*cse(1)*cse(1)*cse(1)
                a_ = 10.0**loga_
                b_ = sb(1)+ sb(2)*tc0 + sb(3)*cse(1) + sb(4)*tc0*cse(1) &
                    + sb(5)*tc0*tc0 + sb(6)*cse(1)*cse(1) &
                    + sb(7)*tc0*tc0*cse(1) + sb(8)*tc0*cse(1)*cse(1) &
                    + sb(9)*tc0*tc0*tc0 + sb(10)*cse(1)*cse(1)*cse(1)
                smoc(k) = a_ * smo2(k)**b_
                !..Calculate snow number concentration (explicit integral, not smo0)
                M0 = smob(k)/smoc(k)
                Mrat = smob(k)*M0*M0*M0
                slam1 = M0 * Lam0
                slam2 = M0 * Lam1
                ns(k) = Mrat*Kap0/slam1                                        &
                    + Mrat*Kap1*M0**mu_s*csg(15)/slam2**cse(15)

                !..Calculate bv_s+2 (th) moment.  Useful for riming.
                loga_ = sa(1) + sa(2)*tc0 + sa(3)*cse(13) &
                    + sa(4)*tc0*cse(13) + sa(5)*tc0*tc0 &
                    + sa(6)*cse(13)*cse(13) + sa(7)*tc0*tc0*cse(13) &
                    + sa(8)*tc0*cse(13)*cse(13) + sa(9)*tc0*tc0*tc0 &
                    + sa(10)*cse(13)*cse(13)*cse(13)
                a_ = 10.0**loga_
                b_ = sb(1)+ sb(2)*tc0 + sb(3)*cse(13) + sb(4)*tc0*cse(13) &
                    + sb(5)*tc0*tc0 + sb(6)*cse(13)*cse(13) &
                    + sb(7)*tc0*tc0*cse(13) + sb(8)*tc0*cse(13)*cse(13) &
                    + sb(9)*tc0*tc0*tc0 + sb(10)*cse(13)*cse(13)*cse(13)
                smoe(k) = a_ * smo2(k)**b_

                !..Calculate 1+(bv_s+1)/2 (th) moment.  Useful for depositional growth.
                loga_ = sa(1) + sa(2)*tc0 + sa(3)*cse(16) &
                    + sa(4)*tc0*cse(16) + sa(5)*tc0*tc0 &
                    + sa(6)*cse(16)*cse(16) + sa(7)*tc0*tc0*cse(16) &
                    + sa(8)*tc0*cse(16)*cse(16) + sa(9)*tc0*tc0*tc0 &
                    + sa(10)*cse(16)*cse(16)*cse(16)
                a_ = 10.0**loga_
                b_ = sb(1)+ sb(2)*tc0 + sb(3)*cse(16) + sb(4)*tc0*cse(16) &
                    + sb(5)*tc0*tc0 + sb(6)*cse(16)*cse(16) &
                    + sb(7)*tc0*tc0*cse(16) + sb(8)*tc0*cse(16)*cse(16) &
                    + sb(9)*tc0*tc0*tc0 + sb(10)*cse(16)*cse(16)*cse(16)
                smof(k) = a_ * smo2(k)**b_
                !..Calculate bm_s + bv_s+2 (th) moment.  Useful for riming into graupel.
                loga_ = sa(1) + sa(2)*tc0 + sa(3)*cse(17) &
                    + sa(4)*tc0*cse(17) + sa(5)*tc0*tc0 &
                    + sa(6)*cse(17)*cse(17) + sa(7)*tc0*tc0*cse(17) &
                    + sa(8)*tc0*cse(17)*cse(17) + sa(9)*tc0*tc0*tc0 &
                    + sa(10)*cse(17)*cse(17)*cse(17)
                a_ = 10.0**loga_
                b_ = sb(1)+ sb(2)*tc0 + sb(3)*cse(17) + sb(4)*tc0*cse(17) &
                    + sb(5)*tc0*tc0 + sb(6)*cse(17)*cse(17) &
                    + sb(7)*tc0*tc0*cse(17) + sb(8)*tc0*cse(17)*cse(17) &
                    + sb(9)*tc0*tc0*tc0 + sb(10)*cse(17)*cse(17)*cse(17)
                smog(k) = a_ * smo2(k)**b_

            enddo

            !+---+-----------------------------------------------------------------+
            !..Calculate y-intercept, slope values for graupel.
            !+---+-----------------------------------------------------------------+
            do k = kte, kts, -1
                lamg = (am_g(idx_bg(k))*cgg(3,1)*ogg2*ng(k)/rg(k))**obmg
                ilamg(k) = 1./lamg
                N0_g(k) = ng(k)*ogg2*lamg**cge(2,1)
            enddo
            ! do k = kte, kts, -1
            !     ygra1 = alog10(max(1.e-9, rg(k)))
            !     zans1 = 3.4 + 2./7.*(ygra1+8.) + rand1
            !     N0_exp = 10.**(zans1)
            !     N0_exp = max(dble(gonv_min), min(N0_exp, dble(gonv_max)))
            !     lam_exp = (N0_exp*am_g*cgg(1)/rg(k))**oge1
            !     lamg = lam_exp * (cgg(3)*ogg2*ogg1)**obmg
            !     ilamg(k) = 1./lamg
            !     N0_g(k) = N0_exp/(cgg(2)*lam_exp) * lamg**cge(2)
            !  enddo
        endif

        !+---+-----------------------------------------------------------------+
        !..Calculate y-intercept, slope values for rain.
        !+---+-----------------------------------------------------------------+
        do k = kte, kts, -1
            lamr = (am_r*crg(3)*org2*nr(k)/rr(k))**obmr
            ilamr(k) = 1./lamr
            mvd_r(k) = (3.0 + mu_r + 0.672) / lamr
            N0_r(k) = nr(k)*org2*lamr**cre(2)
        enddo

        !=================================================================================================================
        !..Compute warm-rain process terms (except evap done later).
        !+---+-----------------------------------------------------------------+

        do k = kts, kte

            !..Rain self-collection follows Seifert, 1994 and drop break-up
            !.. follows Verlinde and Cotton, 1993.                                        RAIN2M
            if (L_qr(k) .and. mvd_r(k).gt. D0r) then

                Ef_rr = max(-0.1, 1.0 - exp(2300.0*(mvd_r(k)-1950.0e-6)))
!!!                Ef_rr = 1.0 - exp(2300.0*(mvd_r(k)-1950.0E-6))
                pnr_rcr(k) = Ef_rr * 2.0*nr(k)*rr(k)
            endif

            mvd_c(k) = D0c
            if (l_qc(k)) then
                if (nc(k).gt.10000.e6) then
                    nu_c = 2
                elseif (nc(k).lt.100.) then
                    nu_c = 15
                else
                    nu_c = nint(nu_c_scale/nc(k)) + 2
                    rand = 0.0
                    if (present(rand2)) then
                        rand = rand2
                    endif
                    nu_c = max(2, min(nu_c+nint(rand), 15))
                endif
                xdc = max(D0c*1.e6, ((rc(k)/(am_r*nc(k)))**obmr) * 1.e6)
                lamc = (nc(k)*am_r* ccg(2,nu_c) * ocg1(nu_c) / rc(k))**obmr
                mvd_c(k) = (3.0+nu_c+0.672) / lamc
                mvd_c(k) = max(d0c, min(mvd_c(k), d0r))
            endif

            !..Autoconversion follows Berry & Reinhardt (1974) with characteristic
            !.. diameters correctly computed from gamma distrib of cloud droplets.
            if (rc(k).gt. 0.01e-3) then
                Dc_g = ((ccg(3,nu_c)*ocg2(nu_c))**obmr / lamc) * 1.E6
                Dc_b = (xDc*xDc*xDc*Dc_g*Dc_g*Dc_g - xDc*xDc*xDc*xDc*xDc*xDc) &
                    **(1./6.)
                zeta1 = 0.5*((6.25E-6*xDc*Dc_b*Dc_b*Dc_b - 0.4) &
                    + abs(6.25E-6*xDc*Dc_b*Dc_b*Dc_b - 0.4))
                zeta = 0.027*rc(k)*zeta1
                taud = 0.5*((0.5*Dc_b - 7.5) + abs(0.5*Dc_b - 7.5)) + R1
                tau  = 3.72/(rc(k)*taud)
                prr_wau(k) = zeta/tau
                prr_wau(k) = min(real(rc(k)*odts, kind=dp), prr_wau(k))
                pnr_wau(k) = prr_wau(k) / (am_r*nu_c*10.*D0r*D0r*D0r)           ! RAIN2M
                pnc_wau(k) = min(real(nc(k)*odts, kind=dp), prr_wau(k) / (am_r*mvd_c(k)*mvd_c(k)*mvd_c(k)))  ! Qc2M
            endif

            !>  - Rain collecting cloud water.  In CE, assume Dc<<Dr and vtc=~0.
            if (L_qr(k) .and. mvd_r(k).gt. D0r .and. mvd_c(k).gt. D0c) then
                lamr = 1./ilamr(k)
                idx = 1 + int(nbr*log(real(mvd_r(k)/Dr(1), kind=dp)) / log(real(Dr(nbr)/Dr(1), kind=dp)))
                idx = min(idx, nbr)
                Ef_rw = t_Efrw(idx, int(mvd_c(k)*1.E6))
                prr_rcw(k) = rhof(k)*t1_qr_qc*Ef_rw*rc(k)*N0_r(k) &
                    *((lamr+fv_r)**(-cre(9)))
                prr_rcw(k) = min(real(rc(k)*odts, kind=dp), prr_rcw(k))
                pnc_rcw(k) = rhof(k)*t1_qr_qc*Ef_rw*nc(k)*N0_r(k)             &
                    *((lamr+fv_r)**(-cre(9)))                          ! Qc2M
                pnc_rcw(k) = min(real(nc(k)*odts, kind=dp), pnc_rcw(k))
            endif

            !>  - Rain collecting aerosols, wet scavenging.
            if (L_qr(k) .and. mvd_r(k).gt. D0r) then
                Ef_ra = Eff_aero(mvd_r(k),0.04E-6,visco(k),rho(k),temp(k),'r')
                lamr = 1./ilamr(k)
                pna_rca(k) = rhof(k)*t1_qr_qc*Ef_ra*nwfa(k)*N0_r(k)           &
                    *((lamr+fv_r)**(-cre(9)))
                pna_rca(k) = min(real(nwfa(k)*odts, kind=dp), pna_rca(k))

                Ef_ra = Eff_aero(mvd_r(k),0.8E-6,visco(k),rho(k),temp(k),'r')
                pnd_rcd(k) = rhof(k)*t1_qr_qc*Ef_ra*nifa(k)*N0_r(k)           &
                    *((lamr+fv_r)**(-cre(9)))
                pnd_rcd(k) = min(real(nifa(k)*odts, kind=dp), pnd_rcd(k))
            endif

        enddo

        !=================================================================================================================
        !..Compute all frozen hydrometeor species' process terms.
        !+---+-----------------------------------------------------------------+
        if (.not. iiwarm) then
            do k = kts, kte
                vts_boost(k) = 1.0
                xDs = 0.0
                if (L_qs(k)) xDs = smoc(k) / smob(k)
                !                orho = 1./rho(k)

                ! if (L_qs(k)) then
                !     xDs = smoc(k) / smob(k)
                !     rho_s2 = max(rho_g(1), min(0.13/xDs, rho_i-100.))
                ! else
                !     xDs = 0.
                !     rho_s2 = 100.
                ! endif

                !..Temperature lookup table indexes.
                tempc = temp(k) - 273.15
                idx_tc = max(1, min(nint(-tempc), 45) )
                idx_t = int( (tempc-2.5)/5. ) - 1
                idx_t = max(1, -idx_t)
                idx_t = min(idx_t, ntb_t)
                it = max(1, min(nint(-tempc), 31) )

                !>  - Cloud water lookup table index.
                if (rc(k).gt. r_c(1)) then
                    nic = nint(log10(rc(k)))
                    do_loop_rc: do nn = nic-1, nic+1
                        n = nn
                        if ( (rc(k)/10.**nn).ge.1.0 .and. (rc(k)/10.**nn).lt.10.0 ) exit do_loop_rc
                    enddo do_loop_rc
                    idx_c = int(rc(k)/10.**n) + 10*(n-nic2) - (n-nic2)
                    idx_c = max(1, min(idx_c, ntb_c))
                else
                    idx_c = 1
                endif

                !>  - Cloud droplet number lookup table index.
                idx_n = nint(1.0 + real(nbc, kind=wp) * log(real(nc(k)/t_Nc(1), kind=dp)) / nic1)
                idx_n = max(1, min(idx_n, nbc))

                !>  - Cloud ice lookup table indexes.
                if (ri(k).gt. r_i(1)) then
                    nii = nint(log10(ri(k)))
                    do_loop_ri: do nn = nii-1, nii+1
                        n = nn
                        if ( (ri(k)/10.**nn).ge.1.0 .and. (ri(k)/10.**nn).lt.10.0 ) exit do_loop_ri
                    enddo do_loop_ri
                    idx_i = int(ri(k)/10.**n) + 10*(n-nii2) - (n-nii2)
                    idx_i = max(1, min(idx_i, ntb_i))
                else
                    idx_i = 1
                endif

                if (ni(k).gt. Nt_i(1)) then
                    nii = nint(log10(ni(k)))
                    do_loop_ni: do nn = nii-1, nii+1
                        n = nn
                        if ( (ni(k)/10.**nn).ge.1.0 .and. (ni(k)/10.**nn).lt.10.0 ) exit do_loop_ni
                    enddo do_loop_ni
                    idx_i1 = int(ni(k)/10.**n) + 10*(n-nii3) - (n-nii3)
                    idx_i1 = max(1, min(idx_i1, ntb_i1))
                else
                    idx_i1 = 1
                endif

                !>  - Rain lookup table indexes.
                if (rr(k).gt. r_r(1)) then
                    nir = nint(log10(rr(k)))
                    do_loop_rr: do nn = nir-1, nir+1
                        n = nn
                        if ( (rr(k)/10.**nn).ge.1.0 .and. (rr(k)/10.**nn).lt.10.0 ) exit do_loop_rr
                    enddo do_loop_rr
                    idx_r = int(rr(k)/10.**n) + 10*(n-nir2) - (n-nir2)
                    idx_r = max(1, min(idx_r, ntb_r))

                    lamr = 1./ilamr(k)
                    lam_exp = lamr * (crg(3)*org2*org1)**bm_r
                    N0_exp = org1*rr(k)/am_r * lam_exp**cre(1)
                    nir = nint(log10(real(N0_exp, kind=dp)))
                    do_loop_nr: do nn = nir-1, nir+1
                        n = nn
                        if ( (N0_exp/10.**nn).ge.1.0 .and. (N0_exp/10.**nn).lt.10.0 ) exit do_loop_nr
                    enddo do_loop_nr
                    idx_r1 = int(N0_exp/10.**n) + 10*(n-nir3) - (n-nir3)
                    idx_r1 = max(1, min(idx_r1, ntb_r1))
                else
                    idx_r = 1
                    idx_r1 = ntb_r1
                endif

                !>  - Snow lookup table index.
                if (rs(k).gt. r_s(1)) then
                    nis = nint(log10(rs(k)))
                    do_loop_rs: do nn = nis-1, nis+1
                        n = nn
                        if ( (rs(k)/10.**nn).ge.1.0 .and. (rs(k)/10.**nn).lt.10.0 ) exit do_loop_rs
                    enddo do_loop_rs
                    idx_s = int(rs(k)/10.**n) + 10*(n-nis2) - (n-nis2)
                    idx_s = max(1, min(idx_s, ntb_s))
                else
                    idx_s = 1
                endif

                !>  - Graupel lookup table index.
                if (rg(k).gt. r_g(1)) then
                    nig = nint(log10(rg(k)))
                    do_loop_rg: do nn = nig-1, nig+1
                        n = nn
                        if ( (rg(k)/10.**nn).ge.1.0 .and. (rg(k)/10.**nn).lt.10.0 ) exit do_loop_rg
                    enddo do_loop_rg
                    idx_g = int(rg(k)/10.**n) + 10*(n-nig2) - (n-nig2)
                    idx_g = max(1, min(idx_g, ntb_g))

                    lamg = 1./ilamg(k)
                    lam_exp = lamg * (cgg(3,1)*ogg2*ogg1)**bm_g
                    N0_exp = ogg1*rg(k)/am_g(idx_bg(k)) * lam_exp**cge(1,1)
                    nig = nint(log10(real(N0_exp, kind=dp)))
                    do_loop_ng: do nn = nig-1, nig+1
                        n = nn
                        if ( (N0_exp/10.**nn).ge.1.0 .and. (N0_exp/10.**nn).lt.10.0 ) exit do_loop_ng
                    enddo do_loop_ng
                    idx_g1 = int(N0_exp/10.**n) + 10*(n-nig3) - (n-nig3)
                    idx_g1 = max(1, min(idx_g1, ntb_g1))
                else
                    idx_g = 1
                    idx_g1 = ntb_g1
                endif

                !..Deposition/sublimation prefactor (from Srivastava & Coen 1992).
                otemp = 1./temp(k)
                rvs = rho(k)*qvsi(k)
                rvs_p = rvs*otemp*(lsub*otemp*oRv - 1.)
                rvs_pp = rvs * ( otemp*(lsub*otemp*oRv - 1.) &
                    *otemp*(lsub*otemp*oRv - 1.) &
                    + (-2.*lsub*otemp*otemp*otemp*oRv) &
                    + otemp*otemp)
                gamsc = lsub*diffu(k)/tcond(k) * rvs_p
                alphsc = 0.5*(gamsc/(1.+gamsc))*(gamsc/(1.+gamsc)) &
                    * rvs_pp/rvs_p * rvs/rvs_p
                alphsc = max(1.E-9, alphsc)
                xsat = ssati(k)
                if (abs(xsat).lt. 1.E-9) xsat=0.
                t1_subl = 4.*PI*( 1.0 - alphsc*xsat &
                    + 2.*alphsc*alphsc*xsat*xsat &
                    - 5.*alphsc*alphsc*alphsc*xsat*xsat*xsat ) &
                    / (1.+gamsc)

                !..Snow collecting cloud water.  In CE, assume Dc<<Ds and vtc=~0.
                if (L_qc(k) .and. mvd_c(k).gt. D0c) then
                    ! xDs = 0.0
                    ! if (L_qs(k)) xDs = smoc(k) / smob(k)
                    if (xDs > d0s) then
                        idx = 1 + int(nbs*log(real(xDs/Ds(1), kind=dp)) / log(real(Ds(nbs)/Ds(1), kind=dp)))
                        idx = min(idx, nbs)
                        Ef_sw = t_Efsw(idx, int(mvd_c(k)*1.E6))
                        prs_scw(k) = rhof(k)*t1_qs_qc*Ef_sw*rc(k)*smoe(k)
                        prs_scw(k) = min(real(rc(k)*odts, kind=dp), prs_scw(k))
                        pnc_scw(k) = rhof(k)*t1_qs_qc*Ef_sw*nc(k)*smoe(k)                ! Qc2M
                        pnc_scw(k) = min(real(nc(k)*odts, kind=dp), pnc_scw(k))
                    endif

                    !..Graupel collecting cloud water.  In CE, assume Dc<<Dg and vtc=~0.
                    if (rg(k).ge. r_g(1) .and. mvd_c(k).gt. D0c) then
                        xDg = (bm_g + mu_g + 1.) * ilamg(k)
                        vtg = rhof(k)*av_g(idx_bg(k))*cgg(6,idx_bg(k))*ogg3 * ilamg(k)**bv_g(idx_bg(k))
                        stoke_g = mvd_c(k)*mvd_c(k)*vtg*rho_w2/(9.*visco(k)*xDg)
                        !..Rime density formula of Cober and List (1993) also used by Milbrandt and Morrison (2014).
                        const_Ri = -1.*(mvd_c(k)*0.5E6)*vtg/MIN(-0.1,tempc)
                        const_Ri = MAX(0.1, MIN(const_Ri, 10.))
                        rime_dens = (0.051 + 0.114*const_Ri - 0.0055*const_Ri*const_Ri)*1000.
                        ! CCPP version has check on xDg > D0g
                        if (xDg > D0g) then
                            if (stoke_g.ge.0.4 .and. stoke_g.le.10.) then
                                Ef_gw = 0.55*log10(2.51*stoke_g)
                            elseif (stoke_g.lt.0.4) then
                                Ef_gw = 0.0
                            elseif (stoke_g.gt.10) then
                                Ef_gw = 0.77
                            endif
                            ! Not sure what to do here - hail increases size rapidly here below melting level.
                            if (temp(k).gt.T_0) Ef_gw = Ef_gw*0.1
                            t1_qg_qc = PI*.25*av_g(idx_bg(k)) * cgg(9,idx_bg(k))
                            prg_gcw(k) = rhof(k)*t1_qg_qc*Ef_gw*rc(k)*N0_g(k) &
                                 *ilamg(k)**cge(9,idx_bg(k))
                            pnc_gcw(k) = rhof(k)*t1_qg_qc*Ef_gw*nc(k)*N0_g(k)           &
                                 *ilamg(k)**cge(9,idx_bg(k))                    ! Qc2M
                            pnc_gcw(k) = min(real(nc(k)*odts, kind=dp), pnc_gcw(k))
                            if (temp(k).lt.T_0) pbg_gcw(k) = prg_gcw(k)/rime_dens
                            ! CCPP version has end check on xDg > D0g
                        endif
                    endif
                endif

                !>  - Snow and graupel collecting aerosols, wet scavenging.
                if (rs(k) .gt. r_s(1)) then
                    Ef_sa = Eff_aero(xDs,0.04E-6,visco(k),rho(k),temp(k),'s')
                    pna_sca(k) = rhof(k)*t1_qs_qc*Ef_sa*nwfa(k)*smoe(k)
                    pna_sca(k) = min(real(nwfa(k)*odts, kind=dp), pna_sca(k))

                    Ef_sa = Eff_aero(xDs,0.8E-6,visco(k),rho(k),temp(k),'s')
                    pnd_scd(k) = rhof(k)*t1_qs_qc*Ef_sa*nifa(k)*smoe(k)
                    pnd_scd(k) = min(real(nifa(k)*odts, kind=dp), pnd_scd(k))
                endif
                if (rg(k) .gt. r_g(1)) then
                    xDg = (bm_g + mu_g + 1.) * ilamg(k)
                    Ef_ga = Eff_aero(xDg,0.04E-6,visco(k),rho(k),temp(k),'g')
                    t1_qg_qc = PI*.25*av_g(idx_bg(k)) * cgg(9,idx_bg(k))
                    pna_gca(k) = rhof(k)*t1_qg_qc*Ef_ga*nwfa(k)*N0_g(k)           &
                        *ilamg(k)**cge(9,idx_bg(k))
                    pna_gca(k) = min(real(nwfa(k)*odts, kind=dp), pna_gca(k))

                    Ef_ga = Eff_aero(xDg,0.8E-6,visco(k),rho(k),temp(k),'g')
                    pnd_gcd(k) = rhof(k)*t1_qg_qc*Ef_ga*nifa(k)*N0_g(k)           &
                        *ilamg(k)**cge(9,idx_bg(k))
                    pnd_gcd(k) = min(real(nifa(k)*odts, kind=dp), pnd_gcd(k))
                endif

                !..Rain collecting snow.  Cannot assume Wisner (1972) approximation
                !.. or Mizuno (1990) approach so we solve the CE explicitly and store
                !.. results in lookup table.
                if (rr(k).ge. r_r(1)) then
                    if (rs(k).ge. r_s(1)) then
                        if (temp(k).lt.T_0) then
                            prr_rcs(k) = -(tmr_racs2(idx_s,idx_t,idx_r1,idx_r) &
                                + tcr_sacr2(idx_s,idx_t,idx_r1,idx_r) &
                                + tmr_racs1(idx_s,idx_t,idx_r1,idx_r) &
                                + tcr_sacr1(idx_s,idx_t,idx_r1,idx_r))
                            prs_rcs(k) = tmr_racs2(idx_s,idx_t,idx_r1,idx_r) &
                                + tcr_sacr2(idx_s,idx_t,idx_r1,idx_r) &
                                - tcs_racs1(idx_s,idx_t,idx_r1,idx_r) &
                                - tms_sacr1(idx_s,idx_t,idx_r1,idx_r)
                            prg_rcs(k) = tmr_racs1(idx_s,idx_t,idx_r1,idx_r) &
                                + tcr_sacr1(idx_s,idx_t,idx_r1,idx_r) &
                                + tcs_racs1(idx_s,idx_t,idx_r1,idx_r) &
                                + tms_sacr1(idx_s,idx_t,idx_r1,idx_r)
                            prr_rcs(k) = max(real(-rr(k)*odts, kind=dp), prr_rcs(k))
                            prs_rcs(k) = max(real(-rs(k)*odts, kind=dp), prs_rcs(k))
                            prg_rcs(k) = min(real((rr(k)+rs(k))*odts, kind=dp), prg_rcs(k))
                            pnr_rcs(k) = tnr_racs1(idx_s,idx_t,idx_r1,idx_r)            &   ! RAIN2M
                                + tnr_racs2(idx_s,idx_t,idx_r1,idx_r)          &
                                + tnr_sacr1(idx_s,idx_t,idx_r1,idx_r)          &
                                + tnr_sacr2(idx_s,idx_t,idx_r1,idx_r)
                            pnr_rcs(k) = min(real(nr(k)*odts, kind=dp), pnr_rcs(k))
                            png_rcs(k) = pnr_rcs(k)
                            !-GT        pbg_rcs(k) = prg_rcs(k)/(0.5*(rho_i+rho_s))
                            pbg_rcs(k) = prg_rcs(k)/rho_i
                        else
                            prs_rcs(k) = -tcs_racs1(idx_s,idx_t,idx_r1,idx_r)           &
                                - tms_sacr1(idx_s,idx_t,idx_r1,idx_r)          &
                                + tmr_racs2(idx_s,idx_t,idx_r1,idx_r)          &
                                + tcr_sacr2(idx_s,idx_t,idx_r1,idx_r)
                            prs_rcs(k) = max(real(-rs(k)*odts, kind=dp), prs_rcs(k))
                            prr_rcs(k) = -prs_rcs(k)
                            !pnr_rcs(k) = tnr_racs2(idx_s,idx_t,idx_r1,idx_r)            &   ! RAIN2M
                            !   + tnr_sacr2(idx_s,idx_t,idx_r1,idx_r)
                        endif
                        !pnr_rcs(k) = MIN(DBLE(nr(k)*odts), pnr_rcs(k))
                    endif

                    !..Rain collecting graupel.  Cannot assume Wisner (1972) approximation
                    !.. or Mizuno (1990) approach so we solve the CE explicitly and store
                    !.. results in lookup table.
                    if (rg(k).ge. r_g(1)) then
                        if (temp(k).lt.T_0) then
                            prg_rcg(k) = tmr_racg(idx_g1,idx_g,idx_table(k),idx_r1,idx_r)  &
                                + tcr_gacr(idx_g1,idx_g,idx_table(k),idx_r1,idx_r)
                            prg_rcg(k) = min(real(rr(k)*odts, kind=dp), prg_rcg(k))
                            prr_rcg(k) = -prg_rcg(k)
                            pnr_rcg(k) = tnr_racg(idx_g1,idx_g,idx_table(k),idx_r1,idx_r)  &   ! RAIN2M
                                + tnr_gacr(idx_g1,idx_g,idx_table(k),idx_r1,idx_r)
                            pnr_rcg(k) = min(real(nr(k)*odts, kind=dp), pnr_rcg(k))
                            !-GT        pbg_rcg(k) = prg_rcg(k)/(0.5*(rho_i+rho_g(idx_bg(k))))
                            pbg_rcg(k) = prg_rcg(k)/rho_i
                        else
                            prr_rcg(k) = tcg_racg(idx_g1,idx_g,idx_table(k),idx_r1,idx_r)
                            prr_rcg(k) = min(real(rg(k)*odts, kind=dp), prr_rcg(k))
                            prg_rcg(k) = -prr_rcg(k)
                            png_rcg(k) = tnr_racg(idx_g1,idx_g,idx_table(k),idx_r1,idx_r)
!!!                    + tnr_gacr(idx_g1,idx_g,idx_table(k),idx_r1,idx_r)
                            png_rcg(k) = min(real(ng(k)*odts, kind=dp), png_rcg(k))
                            pbg_rcg(k) = prg_rcg(k)/rho_g(idx_bg(k))
                            !..Put in explicit drop break-up due to collisions.
                            pnr_rcg(k) = -1.5*tnr_gacr(idx_g1,idx_g,idx_table(k),idx_r1,idx_r)  ! RAIN2M
                        endif
                    endif
                endif

                !+---+-----------------------------------------------------------------+
                !..Next IF block handles only those processes below 0C.
                !+---+-----------------------------------------------------------------+

                if (temp(k).lt.T_0) then

                    vts_boost(k) = 1.0
                    rate_max = (qv(k)-qvsi(k))*rho(k)*odts*0.999

                    !+---+---------------- BEGIN NEW ICE NUCLEATION -----------------------+
                    !..Freezing of supercooled water (rain or cloud) is influenced by dust
                    !.. but still using Bigg 1953 with a temperature adjustment of a few
                    !.. degrees depending on dust concentration.  A default value by way
                    !.. of idx_IN is 1.0 per Liter of air is used when dustyIce flag is
                    !.. false.  Next, a combination of deposition/condensation freezing
                    !.. using DeMott et al (2010) dust nucleation when water saturated or
                    !.. Phillips et al (2008) when below water saturation; else, without
                    !.. dustyIce flag, use the previous Cooper (1986) temperature-dependent
                    !.. value.  Lastly, allow homogeneous freezing of deliquesced aerosols
                    !.. following Koop et al. (2001, Nature).
                    !.. Implemented by T. Eidhammer and G. Thompson 2012Dec18
                    !+---+-----------------------------------------------------------------+

                    ! if (dustyIce) then
                    if (dustyIce .and. (configs%aerosol_aware .or. merra2_aerosol_aware)) then
                        xni = iceDeMott(tempc,qvs(k),qvs(k),qvsi(k),rho(k),nifa(k))
                    else
                        xni = 1.0 *1000. ! Default is 1.0 per Liter
                    endif

                    !>  - Ice nuclei lookup table index.
                    if (xni.gt. Nt_IN(1)) then
                        niin = nint(log10(xni))
                        do_loop_xni: do nn = niin-1, niin+1
                            n = nn
                            if ( (xni/10.**nn).ge.1.0 .and. (xni/10.**nn).lt.10.0 ) exit do_loop_xni
                        enddo do_loop_xni
                        idx_IN = int(xni/10.**n) + 10*(n-niin2) - (n-niin2)
                        idx_IN = max(1, min(idx_IN, ntb_IN))
                    else
                        idx_IN = 1
                    endif

                    !..Freezing of water drops into graupel/cloud ice (Bigg 1953).
                    if (rr(k).gt. r_r(1)) then
                        prg_rfz(k) = tpg_qrfz(idx_r,idx_r1,idx_tc,idx_IN)*odts
                        pri_rfz(k) = tpi_qrfz(idx_r,idx_r1,idx_tc,idx_IN)*odts
                        pni_rfz(k) = tni_qrfz(idx_r,idx_r1,idx_tc,idx_IN)*odts
                        pnr_rfz(k) = tnr_qrfz(idx_r,idx_r1,idx_tc,idx_IN)*odts          ! RAIN2M
                        pnr_rfz(k) = min(real(nr(k)*odts, kind=dp), pnr_rfz(k))
                    elseif (rr(k).gt. R1 .and. temp(k).lt.HGFR) then
                        pri_rfz(k) = rr(k)*odts
                        pni_rfz(k) = nr(k)*odts
                    endif
                    pbg_rfz(k) = prg_rfz(k)/rho_i

                    if (rc(k).gt. r_c(1)) then
                        pri_wfz(k) = tpi_qcfz(idx_c,idx_n,idx_tc,idx_IN)*odts
                        pri_wfz(k) = min(real(rc(k)*odts, kind=dp), pri_wfz(k))
                        pni_wfz(k) = tni_qcfz(idx_c,idx_n,idx_tc,idx_IN)*odts
                        pni_wfz(k) = min(real(nc(k)*odts, kind=dp), pri_wfz(k)/(2.0_dp*xm0i), pni_wfz(k))
                    elseif (rc(k).gt. R1 .and. temp(k).lt.HGFR) then
                        pri_wfz(k) = rc(k)*odts
                        pni_wfz(k) = nc(k)*odts
                    endif

                    !..Deposition nucleation of dust/mineral from DeMott et al (2010)
                    !.. we may need to relax the temperature and ssati constraints.
                    if ( (ssati(k).ge. demott_nuc_ssati) .or. (ssatw(k).gt. eps &
                        .and. temp(k).lt.253.15) ) then

                        if (dustyIce) then
                            xnc = iceDeMott(tempc,qv(k),qvs(k),qvsi(k),rho(k),nifa(k))
                            rand = 0.0
                            if (present(rand3)) then
                                rand = rand3
                            endif
                            xnc = xnc*(1.0 + 50.*rand)
                        else
                            xnc = min(icenuc_max, tno*exp(ato*(t_0-temp(k))))
                        endif
                        xni = ni(k) + (pni_rfz(k)+pni_wfz(k))*dtsave
                        pni_inu(k) = 0.5*(xnc-xni + abs(xnc-xni))*odts
                        pri_inu(k) = min(real(rate_max, kind=dp), xm0i*pni_inu(k))
                        pni_inu(k) = pri_inu(k)/xm0i
                    endif

                    !..Freezing of aqueous aerosols based on Koop et al (2001, Nature)
                    xni = ns(k)+ni(k) + (pni_rfz(k)+pni_wfz(k)+pni_inu(k))*dtsave
                    if ((configs%aerosol_aware .or. merra2_aerosol_aware) .and. homogIce .and. &
                        (xni.le.max_ni) .and.(temp(k).lt.238.).and.(ssati(k).ge.0.4)) then

                        xnc = iceKoop(temp(k),qv(k),qvs(k),nwfa(k), dtsave)
                        pni_iha(k) = xnc*odts
                        pri_iha(k) = min(real(rate_max, kind=dp), xm0i*0.1*pni_iha(k))
                        pni_iha(k) = pri_iha(k)/(xm0i*0.1)
                    endif
                    !+---+------------------ END NEW ICE NUCLEATION -----------------------+


                    !..Deposition/sublimation of cloud ice (Srivastava & Coen 1992).
                    if (L_qi(k)) then
                        lami = (am_i*cig(2)*oig1*ni(k)/ri(k))**obmi
                        ilami = 1./lami
                        xDi = max(real(D0i, kind=dp), (bm_i + mu_i + 1.) * ilami)
                        xmi = am_i*xDi**bm_i
                        oxmi = 1./xmi
                        pri_ide(k) = C_cube*t1_subl*diffu(k)*ssati(k)*rvs &
                            *oig1*cig(5)*ni(k)*ilami

                        if (pri_ide(k) .lt. 0.0) then
                            pri_ide(k) = max(real(-ri(k)*odts, kind=dp), pri_ide(k), real(rate_max, kind=dp))
                            pni_ide(k) = pri_ide(k)*oxmi
                            pni_ide(k) = max(real(-ni(k)*odts, kind=dp), pni_ide(k))
                        else
                            pri_ide(k) = min(pri_ide(k), real(rate_max, kind=dp))
                            prs_ide(k) = (1.0_dp-tpi_ide(idx_i,idx_i1))*pri_ide(k)
                            pri_ide(k) = tpi_ide(idx_i,idx_i1)*pri_ide(k)
                        endif

                        !..Some cloud ice needs to move into the snow category.  Use lookup
                        !.. table that resulted from explicit bin representation of distrib.
                        if ( (idx_i.eq. ntb_i) .or. (xDi.gt. 5.0*D0s) ) then
                            prs_iau(k) = ri(k)*.99*odts
                            pni_iau(k) = ni(k)*.95*odts
                        elseif (xDi.lt. 0.1*D0s) then
                            prs_iau(k) = 0.
                            pni_iau(k) = 0.
                        else
                            prs_iau(k) = tps_iaus(idx_i,idx_i1)*odts
                            prs_iau(k) = min(real(ri(k)*.99*odts, kind=dp), prs_iau(k))
                            pni_iau(k) = tni_iaus(idx_i,idx_i1)*odts
                            pni_iau(k) = min(real(ni(k)*.95*odts, kind=dp), pni_iau(k))
                        endif
                    endif

                    !..Deposition/sublimation of snow/graupel follows Srivastava & Coen
                    !.. (1992).
                    if (l_qs(k)) then
                        c_snow = c_sqrd + (tempc+1.5)*(c_cube-c_sqrd)/(-30.+1.5)
                        c_snow = max(c_sqrd, min(c_snow, c_cube))
                        prs_sde(k) = c_snow*t1_subl*diffu(k)*ssati(k)*rvs &
                            * (t1_qs_sd*smo1(k) &
                            + t2_qs_sd*rhof2(k)*vsc2(k)*smof(k))
                        if (prs_sde(k).lt. 0.) then
                            prs_sde(k) = max(real(-rs(k)*odts, kind=dp), prs_sde(k), real(rate_max, kind=dp))
                        else
                            prs_sde(k) = min(prs_sde(k), real(rate_max, kind=dp))
                        endif
                    endif

                    if (l_qg(k) .and. ssati(k).lt. -eps) then
                        t2_qg_sd = 0.28*sc3*sqrt(av_g(idx_bg(k))) * cgg(11,idx_bg(k))
                        prg_gde(k) = c_cube*t1_subl*diffu(k)*ssati(k)*rvs &
                            * n0_g(k) * (t1_qg_sd*ilamg(k)**cge(10,1) &
                            + t2_qg_sd*vsc2(k)*rhof2(k)*ilamg(k)**cge(11,idx_bg(k)))
                        if (prg_gde(k).lt. 0.) then
                            prg_gde(k) = max(real(-rg(k)*odts, kind=dp), prg_gde(k), real(rate_max, kind=dp))
                            png_gde(k) = prg_gde(k) * ng(k)/rg(k)
                        else
                            prg_gde(k) = min(prg_gde(k), real(rate_max, kind=dp))
                        endif
                    endif

                    !..Snow collecting cloud ice.  In CE, assume Di<<Ds and vti=~0.
                    if (L_qi(k)) then
                        lami = (am_i*cig(2)*oig1*ni(k)/ri(k))**obmi
                        ilami = 1./lami
                        xDi = max(real(D0i, kind=dp), (bm_i + mu_i + 1.) * ilami)
                        xmi = am_i*xDi**bm_i
                        oxmi = 1./xmi
                        if (rs(k).ge. r_s(1)) then
                            prs_sci(k) = t1_qs_qi*rhof(k)*Ef_si*ri(k)*smoe(k)
                            pni_sci(k) = prs_sci(k) * oxmi
                        endif

                        !..Rain collecting cloud ice.  In CE, assume Di<<Dr and vti=~0.
                        if (rr(k).ge. r_r(1) .and. mvd_r(k).gt. 4.*xDi) then
                            lamr = 1./ilamr(k)
                            pri_rci(k) = rhof(k)*t1_qr_qi*Ef_ri*ri(k)*N0_r(k) &
                                *((lamr+fv_r)**(-cre(9)))
                            pnr_rci(k) = rhof(k)*t1_qr_qi*Ef_ri*ni(k)*N0_r(k)           &   ! RAIN2M
                                *((lamr+fv_r)**(-cre(9)))
                            pnr_rci(k) = min(real(nr(k)*odts, kind=dp), pnr_rci(k))
                            pni_rci(k) = pri_rci(k) * oxmi
                            prr_rci(k) = rhof(k)*t2_qr_qi*Ef_ri*ni(k)*N0_r(k) &
                                *((lamr+fv_r)**(-cre(8)))
                            prr_rci(k) = min(real(rr(k)*odts, kind=dp), prr_rci(k))
                            prg_rci(k) = pri_rci(k) + prr_rci(k)
                            pbg_rci(k) = prg_rci(k)/rho_i
                        endif
                    endif

                    !..Ice multiplication from rime-splinters (Hallet & Mossop 1974).
                    if (prg_gcw(k).gt. eps .and. tempc.gt.-8.0) then
                        tf = 0.
                        if (tempc.ge.-5.0 .and. tempc.lt.-3.0) then
                            tf = 0.5*(-3.0 - tempc)
                        elseif (tempc.gt.-8.0 .and. tempc.lt.-5.0) then
                            tf = 0.33333333*(8.0 + tempc)
                        endif
                        pni_ihm(k) = 3.5E8*tf*prg_gcw(k)
                        pri_ihm(k) = xm0i*pni_ihm(k)
                        prs_ihm(k) = prs_scw(k)/(prs_scw(k)+prg_gcw(k)) &
                            * pri_ihm(k)
                        prg_ihm(k) = prg_gcw(k)/(prs_scw(k)+prg_gcw(k)) &
                            * pri_ihm(k)
                    endif

                    !..A portion of rimed snow converts to graupel but some remains snow.
                    !.. Interp from 15 to 95% as riming factor increases from 2.0 to 30.0
                    !.. 0.028 came from (.95-.15)/(30.-2.).  This remains ad-hoc and should
                    !.. be revisited.
                    if (prs_scw(k).gt.rime_threshold*prs_sde(k) .and. prs_sde(k).gt.eps) then
                        r_frac = min(30.0_dp, prs_scw(k)/prs_sde(k))
                        g_frac = min(rime_conversion, 0.15 + (r_frac-2.)*.028)
                        vts_boost(k) = min(1.5, 1.1 + (r_frac-2.)*.016)
                        prg_scw(k) = g_frac*prs_scw(k)
                        png_scw(k) = prg_scw(k)*smo0(k)/rs(k)
                        !..gt      png_scw(k) = prg_scw(k)*ns(k)/rs(k)
                        vts = av_s*xds**bv_s * exp(-fv_s*xds)
                        const_ri = -1.*(mvd_c(k)*0.5e6)*vts/min(-0.1,tempc)
                        const_ri = max(0.1, min(const_ri, 10.))
                        rime_dens = (0.051 + 0.114*const_Ri - 0.0055*const_Ri*const_Ri)*1000.
                        if(rime_dens .lt. 150.) then                                  ! Idea of A. Jensen
                            g_frac = 0.
                            prg_scw(k)=0.
                            png_scw(k)=0.
                        endif
                        pbg_scw(k) = prg_scw(k)/(0.5*(rime_dens+rho_s2))
                        prs_scw(k) = (1. - g_frac)*prs_scw(k)
                    endif
                else

                    !..Melt snow and graupel and enhance from collisions with liquid.
                    !.. We also need to sublimate snow and graupel if subsaturated.
                    if (l_qs(k)) then
                        prr_sml(k) = (tempc*tcond(k)-lvap0*diffu(k)*delQvs(k))       &
                            * (t1_qs_me*smo1(k) + t2_qs_me*rhof2(k)*vsc2(k)*smof(k))
                        if (prr_sml(k) .gt. 0.) then
                            prr_sml(k) = prr_sml(k) + 4218.*olfus*(temp(k)-T_0)       &
                                * (prr_rcs(k)+prs_scw(k))
                        endif
                        prr_sml(k) = min(real(rs(k)*odts, kind=dp), max(zerod0, prr_sml(k)))

                        if (prr_sml(k).gt.0.0) then
                            pnr_sml(k) = smo0(k)/rs(k)*prr_sml(k) * 10.0**(-0.25*(temp(k)-T_0))
                        elseif (ssati(k).lt. 0.) then
                            prs_sde(k) = C_sqrd*t1_subl*diffu(k)*ssati(k)*rvs &
                                * (t1_qs_sd*smo1(k) &
                                + t2_qs_sd*rhof2(k)*vsc2(k)*smof(k))
                            prs_sde(k) = max(real(-rs(k)*odts, kind=dp), prs_sde(k))
                        endif
                    endif

                    if (l_qg(k)) then
                        n0_melt = n0_g(k)
                        if ((rg(k)*ng(k)) .lt. 1.e-4) then
                            lamg = 1./ilamg(k)
                            n0_melt = (1.e-4/rg(k))*ogg2*lamg**cge(2,1)
                        endif
                        t2_qg_me = pi*4.*c_cube*olfus * 0.28*sc3*sqrt(av_g(idx_bg(k))) * cgg(11,idx_bg(k))
                        prr_gml(k) = (tempc*tcond(k)-lvap0*diffu(k)*delQvs(k))       &
                            * N0_melt*(t1_qg_me*ilamg(k)**cge(10,1)           &
                            + t2_qg_me*rhof2(k)*vsc2(k)*ilamg(k)**cge(11,idx_bg(k)))
                        !          if (prr_gml(k) .gt. 0.) then
                        !             prr_gml(k) = prr_gml(k) + 4218.*olfus*(twet(k)-T_0)       &
                        !                                     * (prr_rcg(k)+prg_gcw(k))
                        !          endif
                        prr_gml(k) = min(real(rg(k)*odts, kind=dp), max(0.D0, prr_gml(k)))
                        !           pnr_gml(k) = N0_g(k)*cgg(2)*ilamg(k)**cge(2) / rg(k)         &   ! RAIN2M
                        !                      * prr_gml(k) * 10.0**(-0.5*tempc)

                        if (prr_gml(k) .gt. 0.0) then
                            melt_f = max(0.05, min(prr_gml(k)*dt/rg(k),1.0))
                            !..1000 is density water, 50 is lower limit (max ice density is 800)
                            pbg_gml(k) = prr_gml(k) / max(min(melt_f*rho_g(idx_bg(k)),1000.),50.)
                            !-GT        pnr_gml(k) = prr_gml(k)*ng(k)/rg(k)
                            pnr_gml(k) = prr_gml(k)*ng(k)/rg(k) * 10.0**(-0.33*(temp(k)-T_0))
                        elseif (ssati(k).lt. 0.) then
                            t2_qg_sd = 0.28*Sc3*sqrt(av_g(idx_bg(k))) * cgg(11,idx_bg(k))
                            prg_gde(k) = C_cube*t1_subl*diffu(k)*ssati(k)*rvs &
                                * N0_g(k) * (t1_qg_sd*ilamg(k)**cge(10,1) &
                                + t2_qg_sd*vsc2(k)*rhof2(k)*ilamg(k)**cge(11,idx_bg(k)))
                            prg_gde(k) = max(real(-rg(k)*odts, kind=dp), prg_gde(k))
                            png_gde(k) = prg_gde(k) * ng(k)/rg(k)
                        endif
                    endif

                    !.. This change will be required if users run adaptive time step that
                    !.. results in delta-t that is generally too long to allow cloud water
                    !.. collection by snow/graupel above melting temperature.
                    !.. Credit to Bjorn-Egil Nygaard for discovering.
                    if (dt .gt. 120.) then
                        prr_rcw(k)=prr_rcw(k)+prs_scw(k)+prg_gcw(k)
                        prs_scw(k)=0.
                        prg_gcw(k)=0.
                    endif

                endif
                if (.not. configs%hail_aware) idx_bg(k) = idx_bg1
            enddo
        endif

        !=================================================================================================================
        !..Ensure we do not deplete more hydrometeor species than exists.
        !+---+-----------------------------------------------------------------+
        do k = kts, kte

            !..If ice supersaturated, ensure sum of depos growth terms does not
            !.. deplete more vapor than possibly exists.  If subsaturated, limit
            !.. sum of sublimation terms such that vapor does not reproduce ice
            !.. supersat again.
            sump = pri_inu(k) + pri_ide(k) + prs_ide(k) &
                + prs_sde(k) + prg_gde(k) + pri_iha(k)
            rate_max = (qv(k)-qvsi(k))*rho(k)*odts*0.999
            if ( (sump.gt. eps .and. sump.gt. rate_max) .or. &
                (sump.lt. -eps .and. sump.lt. rate_max) ) then
                ratio = rate_max/sump
                pri_inu(k) = pri_inu(k) * ratio
                pri_ide(k) = pri_ide(k) * ratio
                pni_ide(k) = pni_ide(k) * ratio
                prs_ide(k) = prs_ide(k) * ratio
                prs_sde(k) = prs_sde(k) * ratio
                prg_gde(k) = prg_gde(k) * ratio
                pri_iha(k) = pri_iha(k) * ratio
            endif

            !..Cloud water conservation.
            sump = -prr_wau(k) - pri_wfz(k) - prr_rcw(k) &
                - prs_scw(k) - prg_scw(k) - prg_gcw(k)
            rate_max = -rc(k)*odts
            if (sump.lt. rate_max .and. L_qc(k)) then
                ratio = rate_max/sump
                prr_wau(k) = prr_wau(k) * ratio
                pri_wfz(k) = pri_wfz(k) * ratio
                prr_rcw(k) = prr_rcw(k) * ratio
                prs_scw(k) = prs_scw(k) * ratio
                prg_scw(k) = prg_scw(k) * ratio
                prg_gcw(k) = prg_gcw(k) * ratio
            endif

            !..Cloud ice conservation.
            sump = pri_ide(k) - prs_iau(k) - prs_sci(k) &
                - pri_rci(k)
            rate_max = -ri(k)*odts
            if (sump.lt. rate_max .and. L_qi(k)) then
                ratio = rate_max/sump
                pri_ide(k) = pri_ide(k) * ratio
                prs_iau(k) = prs_iau(k) * ratio
                prs_sci(k) = prs_sci(k) * ratio
                pri_rci(k) = pri_rci(k) * ratio
            endif

            !..Rain conservation.
            sump = -prg_rfz(k) - pri_rfz(k) - prr_rci(k) &
                + prr_rcs(k) + prr_rcg(k)
            rate_max = -rr(k)*odts
            if (sump.lt. rate_max .and. L_qr(k)) then
                ratio = rate_max/sump
                prg_rfz(k) = prg_rfz(k) * ratio
                pri_rfz(k) = pri_rfz(k) * ratio
                prr_rci(k) = prr_rci(k) * ratio
                prr_rcs(k) = prr_rcs(k) * ratio
                prr_rcg(k) = prr_rcg(k) * ratio
            endif

            !..Snow conservation.
            sump = prs_sde(k) - prs_ihm(k) - prr_sml(k) &
                + prs_rcs(k)
            rate_max = -rs(k)*odts
            if (sump.lt. rate_max .and. L_qs(k)) then
                ratio = rate_max/sump
                prs_sde(k) = prs_sde(k) * ratio
                prs_ihm(k) = prs_ihm(k) * ratio
                prr_sml(k) = prr_sml(k) * ratio
                prs_rcs(k) = prs_rcs(k) * ratio
            endif

            !..Graupel conservation.
            sump = prg_gde(k) - prg_ihm(k) - prr_gml(k) &
                + prg_rcg(k)
            rate_max = -rg(k)*odts
            if (sump.lt. rate_max .and. L_qg(k)) then
                ratio = rate_max/sump
                prg_gde(k) = prg_gde(k) * ratio
                prg_ihm(k) = prg_ihm(k) * ratio
                prr_gml(k) = prr_gml(k) * ratio
                prg_rcg(k) = prg_rcg(k) * ratio
            endif

            !..Re-enforce proper mass conservation for subsequent elements in case
            !.. any of the above terms were altered.  Thanks P. Blossey. 2009Sep28
            pri_ihm(k) = prs_ihm(k) + prg_ihm(k)
            ratio = min(abs(prr_rcg(k)), abs(prg_rcg(k)) )
            prr_rcg(k) = ratio * sign(1.0, sngl(prr_rcg(k)))
            prg_rcg(k) = -prr_rcg(k)
            if (temp(k).gt.t_0) then
                ratio = min(abs(prr_rcs(k)), abs(prs_rcs(k)) )
                prr_rcs(k) = ratio * sign(1.0, sngl(prr_rcs(k)))
                prs_rcs(k) = -prr_rcs(k)
            endif

        enddo

        !=================================================================================================================
        !..Calculate tendencies of all species but constrain the number of ice
        !.. to reasonable values.
        !+---+-----------------------------------------------------------------+
        do k = kts, kte
            orho = 1./rho(k)
            lfus2 = lsub - lvap(k)

            !..Aerosol number tendency
            if (configs%aerosol_aware) then
                nwfaten(k) = nwfaten(k) - (pna_rca(k) + pna_sca(k)          &
                    + pna_gca(k) + pni_iha(k)) * orho
                nifaten(k) = nifaten(k) - (pnd_rcd(k) + pnd_scd(k)          &
                    + pnd_gcd(k)) * orho
                if (dustyIce) then
                    nifaten(k) = nifaten(k) - pni_inu(k)*orho
                else
                    nifaten(k) = 0.
                endif
            endif

            !..Water vapor tendency
            qvten(k) = qvten(k) + (-pri_inu(k) - pri_iha(k) - pri_ide(k)   &
                - prs_ide(k) - prs_sde(k) - prg_gde(k)) &
                * orho

            !..Cloud water tendency
            qcten(k) = qcten(k) + (-prr_wau(k) - pri_wfz(k) &
                - prr_rcw(k) - prs_scw(k) - prg_scw(k) &
                - prg_gcw(k)) &
                * orho

            !..Cloud water number tendency
            ncten(k) = ncten(k) + (-pnc_wau(k) - pnc_rcw(k) &
                - pni_wfz(k) - pnc_scw(k) - pnc_gcw(k)) &
                * orho

            !..Cloud water mass/number balance; keep mass-wt mean size between
            !.. 1 and 50 microns.  Also no more than Nt_c_max drops total.
            xrc=max(r1, (qc1d(k) + qcten(k)*dtsave)*rho(k))
            xnc=max(2., (nc1d(k) + ncten(k)*dtsave)*rho(k))
            if (xrc .gt. r1) then
                if (xnc.gt.10000.e6) then
                    nu_c = 2
                elseif (xnc.lt.100.) then
                    nu_c = 15
                else
                    nu_c = nint(nu_c_scale/xnc) + 2
                    rand = 0.0
                    if (present(rand2)) then
                        rand = rand2
                    endif
                    nu_c = max(2, min(nu_c+nint(rand), 15))
                endif
                lamc = (xnc*am_r*ccg(2,nu_c)*ocg1(nu_c)/rc(k))**obmr
                xDc = (bm_r + nu_c + 1.) / lamc
                if (xDc.lt. D0c) then
                    lamc = cce(2,nu_c)/D0c
                    xnc = ccg(1,nu_c)*ocg2(nu_c)*xrc/am_r*lamc**bm_r
                    ncten(k) = (xnc-nc1d(k)*rho(k))*odts*orho
                elseif (xDc.gt. D0r*2.) then
                    lamc = cce(2,nu_c)/(D0r*2.)
                    xnc = ccg(1,nu_c)*ocg2(nu_c)*xrc/am_r*lamc**bm_r
                    ncten(k) = (xnc-nc1d(k)*rho(k))*odts*orho
                endif
            else
                ncten(k) = -nc1d(k)*odts
            endif
            xnc=max(0., (nc1d(k) + ncten(k)*dtsave)*rho(k))
            if (xnc.gt.Nt_c_max) &
                ncten(k) = (Nt_c_max-nc1d(k)*rho(k))*odts*orho

            !..Cloud ice mixing ratio tendency
            qiten(k) = qiten(k) + (pri_inu(k) + pri_iha(k) + pri_ihm(k)    &
                + pri_wfz(k) + pri_rfz(k) + pri_ide(k) &
                - prs_iau(k) - prs_sci(k) - pri_rci(k)) &
                * orho

            !..Cloud ice number tendency.
            niten(k) = niten(k) + (pni_inu(k) + pni_iha(k) + pni_ihm(k)    &
                + pni_wfz(k) + pni_rfz(k) + pni_ide(k) &
                - pni_iau(k) - pni_sci(k) - pni_rci(k)) &
                * orho

            !..Cloud ice mass/number balance; keep mass-wt mean size between
            !.. 5 and 300 microns.  Also no more than 500 xtals per liter.
            xri=max(r1,(qi1d(k) + qiten(k)*dtsave)*rho(k))
            xni=max(r2,(ni1d(k) + niten(k)*dtsave)*rho(k))
            if (xri.gt. r1) then
                lami = (am_i*cig(2)*oig1*xni/xri)**obmi
                ilami = 1./lami
                xdi = (bm_i + mu_i + 1.) * ilami
                if (xdi.lt. 5.e-6) then
                    lami = cie(2)/5.e-6
                    xni = min(max_ni, cig(1)*oig2*xri/am_i*lami**bm_i)
                    niten(k) = (xni-ni1d(k)*rho(k))*odts*orho
                elseif (xdi.gt. 300.e-6) then
                    lami = cie(2)/300.e-6
                    xni = cig(1)*oig2*xri/am_i*lami**bm_i
                    niten(k) = (xni-ni1d(k)*rho(k))*odts*orho
                endif
            else
                niten(k) = -ni1d(k)*odts
            endif
            xni=max(0.,(ni1d(k) + niten(k)*dtsave)*rho(k))
            if (xni.gt.max_ni) &
                niten(k) = (max_ni-ni1d(k)*rho(k))*odts*orho

            !..Rain tendency
            qrten(k) = qrten(k) + (prr_wau(k) + prr_rcw(k) &
                + prr_sml(k) + prr_gml(k) + prr_rcs(k) &
                + prr_rcg(k) - prg_rfz(k) &
                - pri_rfz(k) - prr_rci(k)) &
                * orho

            !..Rain number tendency
            nrten(k) = nrten(k) + (pnr_wau(k) + pnr_sml(k) + pnr_gml(k) &
                - (pnr_rfz(k) + pnr_rcr(k) + pnr_rcg(k) &
                + pnr_rcs(k) + pnr_rci(k) + pni_rfz(k)) ) &
                * orho

            !..Rain mass/number balance; keep median volume diameter between
            !.. 37 microns (D0r*0.75) and 2.5 mm.
            xrr=max(R1,(qr1d(k) + qrten(k)*dtsave)*rho(k))
            xnr=max(R2,(nr1d(k) + nrten(k)*dtsave)*rho(k))
            if (xrr.gt. R1) then
                lamr = (am_r*crg(3)*org2*xnr/xrr)**obmr
                mvd_r(k) = (3.0 + mu_r + 0.672) / lamr
                if (mvd_r(k) .gt. 2.5E-3) then
                    mvd_r(k) = 2.5E-3
                    lamr = (3.0 + mu_r + 0.672) / mvd_r(k)
                    xnr = crg(2)*org3*xrr*lamr**bm_r / am_r
                    nrten(k) = (xnr-nr1d(k)*rho(k))*odts*orho
                elseif (mvd_r(k) .lt. D0r*0.75) then
                    mvd_r(k) = D0r*0.75
                    lamr = (3.0 + mu_r + 0.672) / mvd_r(k)
                    xnr = crg(2)*org3*xrr*lamr**bm_r / am_r
                    nrten(k) = (xnr-nr1d(k)*rho(k))*odts*orho
                endif
            else
                qrten(k) = -qr1d(k)*odts
                nrten(k) = -nr1d(k)*odts
            endif

            !..Snow tendency
            qsten(k) = qsten(k) + (prs_iau(k) + prs_sde(k) + prs_sci(k) &
                + prs_scw(k) + prs_rcs(k) + prs_ide(k) &
                - prs_ihm(k) - prr_sml(k)) &
                * orho

            !..Graupel tendency
            qgten(k) = qgten(k) + (prg_scw(k) + prg_rfz(k) &
                + prg_gde(k) + prg_rcg(k) + prg_gcw(k) &
                + prg_rci(k) + prg_rcs(k) - prg_ihm(k) &
                - prr_gml(k)) &
                * orho

            !..Graupel number tendency
            ngten(k) = ngten(k) + (png_scw(k) + pnr_rfz(k) - png_rcg(k) &
                + pnr_rci(k) + png_rcs(k) + png_gde(k) &
                - pnr_gml(k)) * orho

            !..Graupel volume mixing ratio tendency
            qbten(k) = qbten(k) + (pbg_scw(k) + pbg_rfz(k) &
                + pbg_gcw(k) + pbg_rci(k) + pbg_rcs(k) &
                + pbg_rcg(k) + pbg_sml(k) - pbg_gml(k) &
                + (prg_gde(k) - prg_ihm(k)) /rho_g(idx_bg(k)) ) &
                * orho

            !..Graupel mass/number balance; keep its median volume diameter between
            !.. 3.0 times minimum size (D0g) and 25 mm.
            xrg=max(r1,(qg1d(k) + qgten(k)*dtsave)*rho(k))
            xng=max(r2,(ng1d(k) + ngten(k)*dtsave)*rho(k))
            xrb=max(xrg/rho(k)/rho_g(nrhg),(qb1d(k) + qbten(k)*dtsave))
            xrb=min(xrg/rho(k)/rho_g(1), xrb)

            if (xrg .gt. R1) then
                lamg = (am_g(idx_bg(k))*cgg(3,1)*ogg2*xng/xrg)**obmg
                mvd_g(k) = (3.0 + mu_g + 0.672) / lamg

                if (mvd_g(k) .gt. 25.4E-3) then
                    mvd_g(k) = 25.4E-3
                    lamg = (3.0 + mu_g + 0.672) / mvd_g(k)
                    xng = cgg(2,1)*ogg3*xrg*lamg**bm_g / am_g(idx_bg(k))
                    ngten(k) = (xng-ng1d(k)*rho(k))*odts*orho
                elseif (mvd_g(k) .lt. D0r) then
                    mvd_g(k) = D0r
                    lamg = (3.0 + mu_g + 0.672) / mvd_g(k)
                    xng = cgg(2,1)*ogg3*xrg*lamg**bm_g / am_g(idx_bg(k))
                    ngten(k) = (xng-ng1d(k)*rho(k))*odts*orho
                endif

            else
                qgten(k) = -qg1d(k)*odts
                ngten(k) = -ng1d(k)*odts
                qbten(k) = -qb1d(k)*odts
            endif

            !..Temperature tendency
            if (temp(k).lt.T_0) then
                tten(k) = tten(k) &
                    + ( lsub*ocp(k)*(pri_inu(k) + pri_ide(k) &
                    + prs_ide(k) + prs_sde(k) &
                    + prg_gde(k) + pri_iha(k)) &
                    + lfus2*ocp(k)*(pri_wfz(k) + pri_rfz(k) &
                    + prg_rfz(k) + prs_scw(k) &
                    + prg_scw(k) + prg_gcw(k) &
                    + prg_rcs(k) + prs_rcs(k) &
                    + prr_rci(k) + prg_rcg(k)) &
                    )*orho * (1-IFDRY)
            else
                tten(k) = tten(k) &
                    + ( lfus*ocp(k)*(-prr_sml(k) - prr_gml(k) &
                    - prr_rcg(k) - prr_rcs(k)) &
                    + lsub*ocp(k)*(prs_sde(k) + prg_gde(k)) &
                    )*orho * (1-IFDRY)
            endif

        enddo

        !=================================================================================================================
        !..Update variables for TAU+1 before condensation & sedimention.
        !+---+-----------------------------------------------------------------+
        do k = kts, kte
            temp(k) = t1d(k) + dt*tten(k)
            otemp = 1./temp(k)
            tempc = temp(k) - 273.15
            qv(k) = max(min_qv, qv1d(k) + dt*qvten(k))
            rho(k) = roverrv*pres(k)/(r*temp(k)*(qv(k)+roverrv))
            rhof(k) = sqrt(rho_not/rho(k))
            rhof2(k) = sqrt(rhof(k))
            qvs(k) = rslf(pres(k), temp(k))
            ssatw(k) = qv(k)/qvs(k) - 1.
            if (abs(ssatw(k)).lt. eps) ssatw(k) = 0.0
            diffu(k) = 2.11e-5*(temp(k)/273.15)**1.94 * (101325./pres(k))
            if (tempc .ge. 0.0) then
                visco(k) = (1.718+0.0049*tempc)*1.0e-5
            else
                visco(k) = (1.718+0.0049*tempc-1.2e-5*tempc*tempc)*1.0e-5
            endif
            vsc2(k) = sqrt(rho(k)/visco(k))
            lvap(k) = lvap0 + (2106.0 - 4218.0)*tempc
            tcond(k) = (5.69 + 0.0168*tempc)*1.0e-5 * 418.936
            ocp(k) = 1./(cp2*(1.+0.887*qv(k)))
            lvt2(k)=lvap(k)*lvap(k)*ocp(k)*orv*otemp*otemp

            if (configs%aerosol_aware) then
                nwfa(k) = MAX(nwfa_default*rho(k), (nwfa1d(k) + nwfaten(k)*DT)*rho(k))
            endif
        enddo

        do k = kts, kte
            if ((qc1d(k) + qcten(k)*dt) .gt. r1) then
                rc(k) = (qc1d(k) + qcten(k)*dt)*rho(k)
                nc(k) = max(2., min((nc1d(k)+ncten(k)*dt)*rho(k), nt_c_max))
                if (.not.(configs%aerosol_aware .or. merra2_aerosol_aware)) then
                    nc(k) = Nt_c
                    if (present(lsml)) then
                        if(lsml == 1) then
                            nc(k) = Nt_c_l
                        else
                            nc(k) = Nt_c_o
                        endif
                    endif
                endif
                L_qc(k) = .true.
            else
                rc(k) = R1
                nc(k) = 2.
                L_qc(k) = .false.
            endif

            if ((qi1d(k) + qiten(k)*DT) .gt. R1) then
                ri(k) = (qi1d(k) + qiten(k)*DT)*rho(k)
                ni(k) = max(R2, (ni1d(k) + niten(k)*DT)*rho(k))
                L_qi(k) = .true.
            else
                ri(k) = R1
                ni(k) = R2
                L_qi(k) = .false.
            endif

            if ((qr1d(k) + qrten(k)*DT) .gt. R1) then
                rr(k) = (qr1d(k) + qrten(k)*DT)*rho(k)
                nr(k) = max(R2, (nr1d(k) + nrten(k)*DT)*rho(k))
                L_qr(k) = .true.
                lamr = (am_r*crg(3)*org2*nr(k)/rr(k))**obmr
                mvd_r(k) = (3.0 + mu_r + 0.672) / lamr
                if (mvd_r(k) .gt. 2.5E-3) then
                    mvd_r(k) = 2.5E-3
                    lamr = (3.0 + mu_r + 0.672) / mvd_r(k)
                    nr(k) = crg(2)*org3*rr(k)*lamr**bm_r / am_r
                elseif (mvd_r(k) .lt. D0r*0.75) then
                    mvd_r(k) = D0r*0.75
                    lamr = (3.0 + mu_r + 0.672) / mvd_r(k)
                    nr(k) = crg(2)*org3*rr(k)*lamr**bm_r / am_r
                endif
            else
                rr(k) = R1
                nr(k) = R2
                L_qr(k) = .false.
            endif

            if ((qs1d(k) + qsten(k)*DT) .gt. R1) then
                rs(k) = (qs1d(k) + qsten(k)*DT)*rho(k)
                L_qs(k) = .true.
            else
                rs(k) = R1
                L_qs(k) = .false.
            endif
        enddo

        if (configs%hail_aware) then
            do k = kts, kte
                if ((qg1d(k) + qgten(k)*dt) .gt. r1) then
                    l_qg(k) = .true.
                    rg(k) = (qg1d(k) + qgten(k)*dt)*rho(k)
                    ng(k) = max(r2, (ng1d(k) + ngten(k)*dt)*rho(k))
                    rb(k) = max(rg(k)/rho(k)/rho_g(nrhg), qb1d(k) + qbten(k)*dt)
                    rb(k) = min(rg(k)/rho(k)/rho_g(1), rb(k))
                    idx_bg(k) = max(1,min(nint(rg(k)/rho(k)/rb(k) *0.01)+1,nrhg))
                else
                    rg(k) = r1
                    ng(k) = r2
                    rb(k) = r1/rho(k)/rho_g(nrhg)
                    idx_bg(k) = idx_bg1
                    l_qg(k) = .false.
                endif
            enddo
        else
            do k = kte, kts, -1
                idx_bg(k) = idx_bg1
            enddo
            do k = kte, kts, -1
                if ((qg1d(k) + qgten(k)*dt) .gt. r1) then
                    rg(k) = (qg1d(k) + qgten(k)*dt)*rho(k)
                    ygra1 = log10(max(1.e-9, rg(k)))
                    zans1 = 3.0 + 2./7.*(ygra1+8.)
                    zans1 = max(2., min(zans1, 6.))
                    n0_exp = 10.**(zans1)
                    lam_exp = (n0_exp*am_g(idx_bg(k))*cgg(1,1)/rg(k))**oge1
                    lamg = lam_exp * (cgg(3,1)*ogg2*ogg1)**obmg
                    ng(k) = cgg(2,1)*ogg3*rg(k)*lamg**bm_g / am_g(idx_bg(k))
                    rb(k) = rg(k)/rho(k)/rho_g(idx_bg(k))
                else
                    rg(k) = R1
                    ng(k) = R2
                    rb(k) = R1/rho(k)/rho_g(NRHG)
                    L_qg(k) = .false.
                endif
            enddo
        endif

        !=================================================================================================================
        !..With tendency-updated mixing ratios, recalculate snow moments and
        !.. intercepts/slopes of graupel and rain.
        !+---+-----------------------------------------------------------------+
        if (.not. iiwarm) then
            do k = kts, kte
                smo2(k) = 0.
                smob(k) = 0.
                smoc(k) = 0.
                smod(k) = 0.
            enddo
            do k = kts, kte
                if (.not. L_qs(k)) CYCLE
                tc0 = min(-0.1, temp(k)-273.15)
                smob(k) = rs(k)*oams

                !..All other moments based on reference, 2nd moment.  If bm_s.ne.2,
                !.. then we must compute actual 2nd moment and use as reference.
                if (bm_s.gt.(2.0-1.e-3) .and. bm_s.lt.(2.0+1.e-3)) then
                    smo2(k) = smob(k)
                else
                    loga_ = sa(1) + sa(2)*tc0 + sa(3)*bm_s &
                        + sa(4)*tc0*bm_s + sa(5)*tc0*tc0 &
                        + sa(6)*bm_s*bm_s + sa(7)*tc0*tc0*bm_s &
                        + sa(8)*tc0*bm_s*bm_s + sa(9)*tc0*tc0*tc0 &
                        + sa(10)*bm_s*bm_s*bm_s
                    a_ = 10.0**loga_
                    b_ = sb(1) + sb(2)*tc0 + sb(3)*bm_s &
                        + sb(4)*tc0*bm_s + sb(5)*tc0*tc0 &
                        + sb(6)*bm_s*bm_s + sb(7)*tc0*tc0*bm_s &
                        + sb(8)*tc0*bm_s*bm_s + sb(9)*tc0*tc0*tc0 &
                        + sb(10)*bm_s*bm_s*bm_s
                    smo2(k) = (smob(k)/a_)**(1./b_)
                endif

                !..Calculate bm_s+1 (th) moment.  Useful for diameter calcs.
                loga_ = sa(1) + sa(2)*tc0 + sa(3)*cse(1) &
                    + sa(4)*tc0*cse(1) + sa(5)*tc0*tc0 &
                    + sa(6)*cse(1)*cse(1) + sa(7)*tc0*tc0*cse(1) &
                    + sa(8)*tc0*cse(1)*cse(1) + sa(9)*tc0*tc0*tc0 &
                    + sa(10)*cse(1)*cse(1)*cse(1)
                a_ = 10.0**loga_
                b_ = sb(1)+ sb(2)*tc0 + sb(3)*cse(1) + sb(4)*tc0*cse(1) &
                    + sb(5)*tc0*tc0 + sb(6)*cse(1)*cse(1) &
                    + sb(7)*tc0*tc0*cse(1) + sb(8)*tc0*cse(1)*cse(1) &
                    + sb(9)*tc0*tc0*tc0 + sb(10)*cse(1)*cse(1)*cse(1)
                smoc(k) = a_ * smo2(k)**b_

                !..Calculate bm_s+bv_s (th) moment.  Useful for sedimentation.
                loga_ = sa(1) + sa(2)*tc0 + sa(3)*cse(14) &
                    + sa(4)*tc0*cse(14) + sa(5)*tc0*tc0 &
                    + sa(6)*cse(14)*cse(14) + sa(7)*tc0*tc0*cse(14) &
                    + sa(8)*tc0*cse(14)*cse(14) + sa(9)*tc0*tc0*tc0 &
                    + sa(10)*cse(14)*cse(14)*cse(14)
                a_ = 10.0**loga_
                b_ = sb(1)+ sb(2)*tc0 + sb(3)*cse(14) + sb(4)*tc0*cse(14) &
                    + sb(5)*tc0*tc0 + sb(6)*cse(14)*cse(14) &
                    + sb(7)*tc0*tc0*cse(14) + sb(8)*tc0*cse(14)*cse(14) &
                    + sb(9)*tc0*tc0*tc0 + sb(10)*cse(14)*cse(14)*cse(14)
                smod(k) = a_ * smo2(k)**b_
            enddo

            !+---+-----------------------------------------------------------------+
            !..Calculate y-intercept, slope values for graupel.
            !+---+-----------------------------------------------------------------+
            do k = kte, kts, -1
                lamg = (am_g(idx_bg(k))*cgg(3,1)*ogg2*ng(k)/rg(k))**obmg
                ilamg(k) = 1./lamg
                N0_g(k) = ng(k)*ogg2*lamg**cge(2,1)
            enddo

        endif

        !+---+-----------------------------------------------------------------+
        !..Calculate y-intercept, slope values for rain.
        !+---+-----------------------------------------------------------------+
        do k = kte, kts, -1
            lamr = (am_r*crg(3)*org2*nr(k)/rr(k))**obmr
            ilamr(k) = 1./lamr
            mvd_r(k) = (3.0 + mu_r + 0.672) / lamr
            N0_r(k) = nr(k)*org2*lamr**cre(2)
        enddo

        !=================================================================================================================
        !..Cloud water condensation and evaporation.  Nucleate cloud droplets
        !.. using explicit CCN aerosols with hygroscopicity like sulfates using
        !.. parcel model lookup table results provided by T. Eidhammer.  Evap
        !.. drops using calculation of max drop size capable of evaporating in
        !.. single timestep and explicit number of drops smaller than Dc_star
        !.. from lookup table.
        !+---+-----------------------------------------------------------------+
        do k = kts, kte
            orho = 1./rho(k)
            if ((ssatw(k).gt. eps) .or. (ssatw(k).lt. -eps .and. L_qc(k))) then
                clap = (qv(k)-qvs(k))/(1. + lvt2(k)*qvs(k))
                do n = 1, 3
                    fcd = qvs(k)* exp(lvt2(k)*clap) - qv(k) + clap
                    dfcd = qvs(k)*lvt2(k)* exp(lvt2(k)*clap) + 1.
                    clap = clap - fcd/dfcd
                enddo
                xrc = rc(k) + clap*rho(k)
                xnc = 0.
                if (xrc > R1) then
                    prw_vcd(k) = clap*odt
                    !+---+-----------------------------------------------------------------+ !  DROPLET NUCLEATION
                    if (clap .gt. eps) then
                        if (configs%aerosol_aware .or. merra2_aerosol_aware) then
                            rand = 0.0
                            if (present(rand3)) then
                                rand = rand3
                            endif
                            xnc = max(2., activ_ncloud(temp(k), w1d(k)+rand, nwfa(k), lsml))
                        else
                            xnc = Nt_c
                            if (present(lsml)) then
                                if(lsml == 1) then
                                    xnc = Nt_c_l
                                else
                                    xnc = Nt_c_o
                                endif
                            endif
                        endif
                        pnc_wcd(k) = 0.5*(xnc-nc(k) + abs(xnc-nc(k)))*odts*orho

                        ! Be careful here: initial cloud evaporation can increase aerosols
                        if (configs%aerosol_aware) nwfaten(k) = nwfaten(k) - pnc_wcd(k)

                        !+---+-----------------------------------------------------------------+ !  EVAPORATION
                    elseif (clap .lt. -eps .AND. ssatw(k).lt.-1.e-6 .and. configs%aerosol_aware) then
                        tempc = temp(k) - 273.15
                        otemp = 1./temp(k)
                        rvs = rho(k)*qvs(k)
                        rvs_p = rvs*otemp*(lvap(k)*otemp*oRv - 1.)
                        rvs_pp = rvs * ( otemp*(lvap(k)*otemp*oRv - 1.) &
                            *otemp*(lvap(k)*otemp*oRv - 1.) &
                            + (-2.*lvap(k)*otemp*otemp*otemp*oRv) &
                            + otemp*otemp)
                        gamsc = lvap(k)*diffu(k)/tcond(k) * rvs_p
                        alphsc = 0.5*(gamsc/(1.+gamsc))*(gamsc/(1.+gamsc)) &
                            * rvs_pp/rvs_p * rvs/rvs_p
                        alphsc = max(1.e-9, alphsc)
                        xsat = ssatw(k)
                        if (abs(xsat).lt. 1.e-9) xsat=0.
                        t1_evap = 2.*pi*( 1.0 - alphsc*xsat  &
                            + 2.*alphsc*alphsc*xsat*xsat  &
                            - 5.*alphsc*alphsc*alphsc*xsat*xsat*xsat ) &
                            / (1.+gamsc)

                        dc_star = dsqrt(-2.d0*dt * t1_evap/(2.*pi) &
                            * 4.*diffu(k)*ssatw(k)*rvs/rho_w2)
                        idx_d = max(1, min(int(1.e6*dc_star), nbc))

                        idx_n = nint(1.0 + float(nbc) * dlog(nc(k)/t_nc(1)) / nic1)
                        idx_n = max(1, min(idx_n, nbc))

                        !>  - Cloud water lookup table index.
                        if (rc(k).gt. r_c(1)) then
                            nic = nint(log10(rc(k)))
                            do_loop_rc_cond: do nn = nic-1, nic+1
                                n = nn
                                if ( (rc(k)/10.**nn).ge.1.0 .and. (rc(k)/10.**nn).lt.10.0 ) exit do_loop_rc_cond
                            enddo do_loop_rc_cond
                            idx_c = int(rc(k)/10.**n) + 10*(n-nic2) - (n-nic2)
                            idx_c = max(1, min(idx_c, ntb_c))
                        else
                            idx_c = 1
                        endif

                        !prw_vcd(k) = MAX(DBLE(-rc(k)*orho*odt),                     &
                        !           -tpc_wev(idx_d, idx_c, idx_n)*orho*odt)
                        prw_vcd(k) = max(real(-rc(k)*0.99*orho*odt, kind=dp), prw_vcd(k))
                        pnc_wcd(k) = max(real(-nc(k)*0.99*orho*odt, &
                            kind=dp), real(-tnc_wev(idx_d, idx_c, idx_n)*orho*odt, kind=dp))
                        ! Be careful here: initial cloud evaporation can increase aerosols
                        if (configs%aerosol_aware) nwfaten(k) = nwfaten(k) - pnc_wcd(k)
                    endif
                else
                    prw_vcd(k) = -rc(k)*orho*odt
                    pnc_wcd(k) = -nc(k)*orho*odt
                    ! Be careful here: initial cloud evaporation can increase aerosols
                    ! if (configs%aerosol_aware) &
                    ! nwfaten(k) = nwfaten(k) - pnc_wcd(k)
                endif

                !+---+-----------------------------------------------------------------+

                qvten(k) = qvten(k) - prw_vcd(k)
                qcten(k) = qcten(k) + prw_vcd(k)
                ncten(k) = ncten(k) + pnc_wcd(k)

                tten(k) = tten(k) + lvap(k)*ocp(k)*prw_vcd(k)*(1-IFDRY)
                rc(k) = max(R1, (qc1d(k) + dt*qcten(k))*rho(k))
                if (rc(k).eq.R1) l_qc(k) = .false.
                nc(k) = max(2., min((nc1d(k)+ncten(k)*dt)*rho(k), nt_c_max))
                if (.not.(configs%aerosol_aware .or. merra2_aerosol_aware)) then
                    nc(k) = Nt_c
                    if (present(lsml)) then
                        if(lsml == 1) then
                            nc(k) = Nt_c_l
                        else
                            nc(k) = Nt_c_o
                        endif
                    endif
                endif                   
                qv(k) = max(min_qv, qv1d(k) + DT*qvten(k))
                temp(k) = t1d(k) + DT*tten(k)
                rho(k) = RoverRv*pres(k)/(R*temp(k)*(qv(k)+RoverRv))
                qvs(k) = rslf(pres(k), temp(k))
                ssatw(k) = qv(k)/qvs(k) - 1.
            endif
        enddo

        !=================================================================================================================
        !.. If still subsaturated, allow rain to evaporate, following
        !.. Srivastava & Coen (1992).
        !+---+-----------------------------------------------------------------+
        do k = kts, kte
            if ( (ssatw(k).lt. -eps) .and. L_qr(k) &
                .and. (.not.(prw_vcd(k).gt. 0.)) ) then
                tempc = temp(k) - 273.15
                otemp = 1./temp(k)
                orho = 1./rho(k)
                rhof(k) = sqrt(rho_not*orho)
                rhof2(k) = sqrt(rhof(k))
                diffu(k) = 2.11e-5*(temp(k)/273.15)**1.94 * (101325./pres(k))
                if (tempc .ge. 0.0) then
                    visco(k) = (1.718+0.0049*tempc)*1.0e-5
                else
                    visco(k) = (1.718+0.0049*tempc-1.2e-5*tempc*tempc)*1.0e-5
                endif
                vsc2(k) = sqrt(rho(k)/visco(k))
                lvap(k) = lvap0 + (2106.0 - 4218.0)*tempc
                tcond(k) = (5.69 + 0.0168*tempc)*1.0E-5 * 418.936
                ocp(k) = 1./(Cp2*(1.+0.887*qv(k)))

                rvs = rho(k)*qvs(k)
                rvs_p = rvs*otemp*(lvap(k)*otemp*oRv - 1.)
                rvs_pp = rvs * ( otemp*(lvap(k)*otemp*oRv - 1.) &
                    *otemp*(lvap(k)*otemp*oRv - 1.) &
                    + (-2.*lvap(k)*otemp*otemp*otemp*oRv) &
                    + otemp*otemp)
                gamsc = lvap(k)*diffu(k)/tcond(k) * rvs_p
                alphsc = 0.5*(gamsc/(1.+gamsc))*(gamsc/(1.+gamsc)) &
                    * rvs_pp/rvs_p * rvs/rvs_p
                alphsc = max(1.E-9, alphsc)
                xsat   = min(-1.E-9, ssatw(k))
                t1_evap = 2.*PI*( 1.0 - alphsc*xsat  &
                    + 2.*alphsc*alphsc*xsat*xsat  &
                    - 5.*alphsc*alphsc*alphsc*xsat*xsat*xsat ) &
                    / (1.+gamsc)

                lamr = 1./ilamr(k)
                !..Rapidly eliminate near zero values when low humidity (<95%)
                if (qv(k)/qvs(k) .lt. 0.95 .AND. rr(k)*orho.le.1.E-8) then
                    prv_rev(k) = rr(k)*orho*odts
                else
                    prv_rev(k) = t1_evap*diffu(k)*(-ssatw(k))*N0_r(k)*rvs &
                        * (t1_qr_ev*ilamr(k)**cre(10) &
                        + t2_qr_ev*vsc2(k)*rhof2(k)*((lamr+0.5*fv_r)**(-cre(11))))
                    rate_max = min((rr(k)*orho*odts), (qvs(k)-qv(k))*odts)
                    prv_rev(k) = min(real(rate_max, kind=dp), prv_rev(k)*orho)

                    !..TEST: G. Thompson  10 May 2013
                    !..Reduce the rain evaporation in same places as melting graupel occurs.
                    !..Rationale: falling and simultaneous melting graupel in subsaturated
                    !..regions will not melt as fast because particle temperature stays
                    !..at 0C.  Also not much shedding of the water from the graupel so
                    !..likely that the water-coated graupel evaporating much slower than
                    !..if the water was immediately shed off.
                    IF (prr_gml(k).gt.0.0) THEN
                        eva_factor = min(1.0, 0.01+(0.99-0.01)*(tempc/20.0))
                        prv_rev(k) = prv_rev(k)*eva_factor
                    ENDIF
                endif

                pnr_rev(k) = min(real(nr(k)*0.99*orho*odts, kind=dp),                  &   ! RAIN2M
                    prv_rev(k) * nr(k)/rr(k))

                qrten(k) = qrten(k) - prv_rev(k)
                qvten(k) = qvten(k) + prv_rev(k)
                nrten(k) = nrten(k) - pnr_rev(k)
                nwfaten(k) = nwfaten(k) + pnr_rev(k)
                tten(k) = tten(k) - lvap(k)*ocp(k)*prv_rev(k)*(1-IFDRY)

                rr(k) = max(r1, (qr1d(k) + dt*qrten(k))*rho(k))
                qv(k) = max(1.e-10, qv1d(k) + dt*qvten(k))
                nr(k) = max(r2, (nr1d(k) + dt*nrten(k))*rho(k))
                temp(k) = t1d(k) + DT*tten(k)
                rho(k) = RoverRv*pres(k)/(R*temp(k)*(qv(k)+RoverRv))
            endif
        enddo
#if defined(mpas)
        do k = kts, kte
            evapprod(k) = prv_rev(k) - (min(zeroD0,prs_sde(k)) + &
                min(zeroD0,prg_gde(k)))
            rainprod(k) = prr_wau(k) + prr_rcw(k) + prs_scw(k) + &
                prg_scw(k) + prs_iau(k) + &
                prg_gcw(k) + prs_sci(k) + &
                pri_rci(k)
        enddo
#endif

        !=================================================================================================================
        !..Find max terminal fallspeed (distribution mass-weighted mean
        !.. velocity) and use it to determine if we need to split the timestep
        !.. (var nstep>1).  Either way, only bother to do sedimentation below
        !.. 1st level that contains any sedimenting particles (k=ksed1 on down).
        !.. New in v3.0+ is computing separate for rain, ice, snow, and
        !.. graupel species thus making code faster with credit to J. Schmidt.
        !+---+-----------------------------------------------------------------+
        nstep = 0
        onstep(:) = 1.0
        ksed1(:) = 1
        do k = kte+1, kts, -1
            vtrk(k) = 0.
            vtnrk(k) = 0.
            vtik(k) = 0.
            vtnik(k) = 0.
            vtsk(k) = 0.
            vtgk(k) = 0.
            vtngk(k) = 0.
            vtck(k) = 0.
            vtnck(k) = 0.
        enddo
        if (any(l_qr .eqv. .true.)) then
            do k = kte, kts, -1
                vtr = 0.
                rhof(k) = sqrt(rho_not/rho(k))

                if (rr(k).gt. R1) then
                    lamr = (am_r*crg(3)*org2*nr(k)/rr(k))**obmr
                    vtr = rhof(k)*av_r*crg(6)*org3 * lamr**cre(3)                 &
                        *((lamr+fv_r)**(-cre(6)))
                    vtrk(k) = vtr
                    ! First below is technically correct:
                    !         vtr = rhof(k)*av_r*crg(5)*org2 * lamr**cre(2)                 &
                    !                     *((lamr+fv_r)**(-cre(5)))
                    ! Test: make number fall faster (but still slower than mass)
                    ! Goal: less prominent size sorting
                    vtr = rhof(k)*av_r*crg(7)/crg(12) * lamr**cre(12)             &
                        *((lamr+fv_r)**(-cre(7)))
                    vtnrk(k) = vtr
                else
                    vtrk(k) = vtrk(k+1)
                    vtnrk(k) = vtnrk(k+1)
                endif

                if (max(vtrk(k),vtnrk(k)) .gt. 1.e-3) then
                    ksed1(1) = max(ksed1(1), k)
                    delta_tp = dzq(k)/(max(vtrk(k),vtnrk(k)))
                    nstep = max(nstep, int(dt/delta_tp + 1.))
                endif
            enddo
            if (ksed1(1) .eq. kte) ksed1(1) = kte-1
            if (nstep .gt. 0) onstep(1) = 1./real(nstep)
        endif

    !+---+-----------------------------------------------------------------+

        if (any(l_qc .eqv. .true.)) then
            hgt_agl = 0.
            do_loop_hgt_agl : do k = kts, kte-1
                if (rc(k) .gt. R2) ksed1(5) = k
                hgt_agl = hgt_agl + dzq(k)
                if (hgt_agl .gt. 500.0) exit do_loop_hgt_agl
            enddo do_loop_hgt_agl

            do k = ksed1(5), kts, -1
                vtc = 0.
                if (rc(k) .gt. R1 .and. w1d(k) .lt. 1.E-1) then
                    if (nc(k).gt.10000.e6) then
                        nu_c = 2
                    elseif (nc(k).lt.100.) then
                        nu_c = 15
                    else
                        nu_c = nint(nu_c_scale/nc(k)) + 2
                        rand = 0.0
                        if (present(rand2)) then
                            rand = rand2
                        endif
                        nu_c = max(2, min(nu_c+nint(rand), 15))
                    endif
                    lamc = (nc(k)*am_r*ccg(2,nu_c)*ocg1(nu_c)/rc(k))**obmr
                    ilamc = 1./lamc
                    vtc = rhof(k)*av_c*ccg(5,nu_c)*ocg2(nu_c) * ilamc**bv_c
                    vtck(k) = vtc
                    vtc = rhof(k)*av_c*ccg(4,nu_c)*ocg1(nu_c) * ilamc**bv_c
                    vtnck(k) = vtc
                endif
            enddo
        endif

    !+---+-----------------------------------------------------------------+

        if (.not. iiwarm) then
            if (any(l_qi .eqv. .true.)) then

                nstep = 0
                do k = kte, kts, -1
                    vti = 0.

                    if (ri(k).gt. R1) then
                        lami = (am_i*cig(2)*oig1*ni(k)/ri(k))**obmi
                        ilami = 1./lami
                        vti = rhof(k)*av_i*cig(3)*oig2 * ilami**bv_i
                        vtik(k) = vti
                        ! First below is technically correct:
                        !          vti = rhof(k)*av_i*cig(4)*oig1 * ilami**bv_i
                        ! Goal: less prominent size sorting
                        vti = rhof(k)*av_i*cig(6)/cig(7) * ilami**bv_i
                        vtnik(k) = vti
                    else
                        vtik(k) = vtik(k+1)
                        vtnik(k) = vtnik(k+1)
                    endif

                    if (vtik(k) .gt. 1.e-3) then
                        ksed1(2) = max(ksed1(2), k)
                        delta_tp = dzq(k)/vtik(k)
                        nstep = max(nstep, int(dt/delta_tp + 1.))
                    endif
                enddo
                if (ksed1(2) .eq. kte) ksed1(2) = kte-1
                if (nstep .gt. 0) onstep(2) = 1./real(nstep)
            endif

        !+---+-----------------------------------------------------------------+

            if (any(l_qs .eqv. .true.)) then

                nstep = 0
                do k = kte, kts, -1
                    vts = 0.

                    if (rs(k).gt. R1) then
                        xDs = smoc(k) / smob(k)
                        Mrat = 1./xDs
                        ils1 = 1./(Mrat*Lam0 + fv_s)
                        ils2 = 1./(Mrat*Lam1 + fv_s)
                        t1_vts = Kap0*csg(4)*ils1**cse(4)
                        t2_vts = Kap1*Mrat**mu_s*csg(10)*ils2**cse(10)
                        ils1 = 1./(Mrat*Lam0)
                        ils2 = 1./(Mrat*Lam1)
                        t3_vts = Kap0*csg(1)*ils1**cse(1)
                        t4_vts = Kap1*Mrat**mu_s*csg(7)*ils2**cse(7)
                        vts = rhof(k)*av_s * (t1_vts+t2_vts)/(t3_vts+t4_vts)
                        if (prr_sml(k) .gt. 0.0) then
                            SR = rs(k)/(rs(k)+rr(k))
                            vtsk(k) = vts*SR + (1.-SR)*vtrk(k)
                        else
                            vtsk(k) = vts*vts_boost(k)
                        endif
                    else
                        vtsk(k) = vtsk(k+1)
                    endif

                    if (vtsk(k) .gt. 1.e-3) then
                        ksed1(3) = max(ksed1(3), k)
                        delta_tp = dzq(k)/vtsk(k)
                        nstep = max(nstep, int(dt/delta_tp + 1.))
                    endif
                enddo
                if (ksed1(3) .eq. kte) ksed1(3) = kte-1
                if (nstep .gt. 0) onstep(3) = 1./real(nstep)
            endif

        !+---+-----------------------------------------------------------------+

            if (ANY(L_qg .eqv. .true.)) then
                nstep = 0
                do k = kte, kts, -1
                    vtg = 0.

                    if (rg(k).gt. R1) then
                        if (configs%hail_aware) then
                            xrho_g = MAX(rho_g(1),MIN(rg(k)/rho(k)/rb(k),rho_g(NRHG)))
                            afall = a_coeff*((4.0*xrho_g*9.8)/(3.0*rho(k)))**b_coeff
                            afall = afall * visco(k)**(1.0-2.0*b_coeff)
                        else
                            afall = av_g_old
                            bfall = bv_g_old
                        endif
                        vtg = rhof(k)*afall*cgg(6,idx_bg(k))*ogg3 * ilamg(k)**bfall
                        vtgk(k) = vtg
                        ! Goal: less prominent size sorting
                        !    the ELSE section below is technically (mathematically) correct:
                        if (mu_g .eq. 0) then
                            vtg = rhof(k)*afall*cgg(7,idx_bg(k))/cgg(12,idx_bg(k)) * ilamg(k)**bfall
                        else
                            vtg = rhof(k)*afall*cgg(8,idx_bg(k))*ogg2 * ilamg(k)**bfall
                        endif
                        vtngk(k) = vtg
#if !defined(mpas)
                        if (temp(k).gt. T_0) then
                            vtgk(k) = MAX(vtg, vtrk(k))
                        else
                            vtgk(k) = vtg
                        endif
#endif                         
                    else
                        vtgk(k) = vtgk(k+1)
                        vtngk(k) = vtngk(k+1)
                    endif

                    if (vtgk(k) .gt. 1.e-3) then
                        ksed1(4) = max(ksed1(4), k)
                        delta_tp = dzq(k)/vtgk(k)
                        nstep = max(nstep, int(dt/delta_tp + 1.))
                    endif
                enddo
                if (ksed1(4) .eq. kte) ksed1(4) = kte-1
                if (nstep .gt. 0) onstep(4) = 1./real(nstep)
            endif
        endif

    !=================================================================================================================
    !..Sedimentation of mixing ratio is the integral of v(D)*m(D)*N(D)*dD,
    !.. whereas neglect m(D) term for number concentration.  Therefore,
    !.. cloud ice has proper differential sedimentation.
    !.. New in v3.0+ is computing separate for rain, ice, snow, and
    !.. graupel species thus making code faster with credit to J. Schmidt.
    !.. Bug fix, 2013Nov01 to tendencies using rho(k+1) correction thanks to
    !.. Eric Skyllingstad.
    !+---+-----------------------------------------------------------------+

        if (any(l_qr .eqv. .true.)) then
            nstep = nint(1./onstep(1))
            if(.not. sedi_semi) then
                do n = 1, nstep
                    do k = kte, kts, -1
                        sed_r(k) = vtrk(k)*rr(k)
                        sed_n(k) = vtnrk(k)*nr(k)
                    enddo
                    k = kte
                    odzq = 1./dzq(k)
                    orho = 1./rho(k)
                    qrten(k) = qrten(k) - sed_r(k)*odzq*onstep(1)*orho
                    nrten(k) = nrten(k) - sed_n(k)*odzq*onstep(1)*orho
                    rr(k) = max(r1, rr(k) - sed_r(k)*odzq*dt*onstep(1))
                    nr(k) = max(r2, nr(k) - sed_n(k)*odzq*dt*onstep(1))
#if !defined(mpas)
                    pfll1(k) = pfll1(k) + sed_r(k)*DT*onstep(1)
#endif
                    do k = ksed1(1), kts, -1
                        odzq = 1./dzq(k)
                        orho = 1./rho(k)
                        qrten(k) = qrten(k) + (sed_r(k+1)-sed_r(k)) &
                            *odzq*onstep(1)*orho
                        nrten(k) = nrten(k) + (sed_n(k+1)-sed_n(k)) &
                            *odzq*onstep(1)*orho
                        rr(k) = max(r1, rr(k) + (sed_r(k+1)-sed_r(k)) &
                            *odzq*dt*onstep(1))
                        nr(k) = max(r2, nr(k) + (sed_n(k+1)-sed_n(k)) &
                            *odzq*DT*onstep(1))
#if !defined(mpas)
                        pfll1(k) = pfll1(k) + sed_r(k)*DT*onstep(1)
#endif
                    enddo

                    if (rr(kts).gt.R1*1000.) &
                        pptrain = pptrain + sed_r(kts)*DT*onstep(1)
                enddo
            else !if(.not. sedi_semi)
                niter = 1
                dtcfl = dt
                niter = int(nstep/max(decfl_,1)) + 1
                dtcfl = dt/niter
                do n = 1, niter
                    rr_tmp(:) = rr(:)
                    nr_tmp(:) = nr(:)
                    call semi_lagrange_sedim(kte,dzq,vtrk,rr,rainsfc,pfll,dtcfl,R1)
                    call semi_lagrange_sedim(kte,dzq,vtnrk,nr,vtr,pdummy,dtcfl,R2)
                    do k = kts, kte
                        orhodt = 1./(rho(k)*dt)
                        qrten(k) = qrten(k) + (rr(k) - rr_tmp(k)) * orhodt
                        nrten(k) = nrten(k) + (nr(k) - nr_tmp(k)) * orhodt
                        pfll1(k) = pfll1(k) + pfll(k)
                    enddo
                    pptrain = pptrain + rainsfc

                    do k = kte+1, kts, -1
                        vtrk(k) = 0.
                        vtnrk(k) = 0.
                    enddo
                    do k = kte, kts, -1
                        vtr = 0.
                        if (rr(k).gt. R1) then
                            lamr = (am_r*crg(3)*org2*nr(k)/rr(k))**obmr
                            vtr = rhof(k)*av_r*crg(6)*org3 * lamr**cre(3)           &
                                *((lamr+fv_r)**(-cre(6)))
                            vtrk(k) = vtr
                            ! First below is technically correct:
                            !         vtr = rhof(k)*av_r*crg(5)*org2 * lamr**cre(2)                &
                            !                     *((lamr+fv_r)**(-cre(5)))
                            ! Test: make number fall faster (but still slower than mass)
                            ! Goal: less prominent size sorting
                            vtr = rhof(k)*av_r*crg(7)/crg(12) * lamr**cre(12)       &
                                *((lamr+fv_r)**(-cre(7)))
                            vtnrk(k) = vtr
                        endif
                    enddo
                enddo
            endif! if(.not. sedi_semi)
        endif
    !+---+-----------------------------------------------------------------+

        if (any(l_qc .eqv. .true.)) then

            do k = kte, kts, -1
                sed_c(k) = vtck(k)*rc(k)
                sed_n(k) = vtnck(k)*nc(k)
            enddo
            do k = ksed1(5), kts, -1
                odzq = 1./dzq(k)
                orho = 1./rho(k)
                qcten(k) = qcten(k) + (sed_c(k+1)-sed_c(k)) *odzq*orho
                ncten(k) = ncten(k) + (sed_n(k+1)-sed_n(k)) *odzq*orho
                rc(k) = max(r1, rc(k) + (sed_c(k+1)-sed_c(k)) *odzq*dt)
                nc(k) = max(10., nc(k) + (sed_n(k+1)-sed_n(k)) *odzq*dt)
            enddo
        endif

        !+---+-----------------------------------------------------------------+

        if (any(l_qi .eqv. .true.)) then

            nstep = nint(1./onstep(2))
            do n = 1, nstep
                do k = kte, kts, -1
                    sed_i(k) = vtik(k)*ri(k)
                    sed_n(k) = vtnik(k)*ni(k)
                enddo
                k = kte
                odzq = 1./dzq(k)
                orho = 1./rho(k)
                qiten(k) = qiten(k) - sed_i(k)*odzq*onstep(2)*orho
                niten(k) = niten(k) - sed_n(k)*odzq*onstep(2)*orho
                ri(k) = max(r1, ri(k) - sed_i(k)*odzq*dt*onstep(2))
                ni(k) = max(r2, ni(k) - sed_n(k)*odzq*dt*onstep(2))
#if !defined(mpas)
                pfil1(k) = pfil1(k) + sed_i(k)*DT*onstep(2)
#endif
                do k = ksed1(2), kts, -1
                    odzq = 1./dzq(k)
                    orho = 1./rho(k)
                    qiten(k) = qiten(k) + (sed_i(k+1)-sed_i(k)) &
                        *odzq*onstep(2)*orho
                    niten(k) = niten(k) + (sed_n(k+1)-sed_n(k)) &
                        *odzq*onstep(2)*orho
                    ri(k) = max(r1, ri(k) + (sed_i(k+1)-sed_i(k)) &
                        *odzq*dt*onstep(2))
                    ni(k) = max(r2, ni(k) + (sed_n(k+1)-sed_n(k)) &
                        *odzq*DT*onstep(2))
#if !defined(mpas)
                    pfil1(k) = pfil1(k) + sed_i(k)*DT*onstep(2)
#endif
                enddo

                if (ri(kts).gt.R1*1000.) &
                    pptice = pptice + sed_i(kts)*DT*onstep(2)
            enddo
        endif

        !+---+-----------------------------------------------------------------+

        if (any(l_qs .eqv. .true.)) then

            nstep = nint(1./onstep(3))
            do n = 1, nstep
                do k = kte, kts, -1
                    sed_s(k) = vtsk(k)*rs(k)
                enddo
                k = kte
                odzq = 1./dzq(k)
                orho = 1./rho(k)
                qsten(k) = qsten(k) - sed_s(k)*odzq*onstep(3)*orho
                rs(k) = max(r1, rs(k) - sed_s(k)*odzq*dt*onstep(3))
#if !defined(mpas)
                pfil1(k) = pfil1(k) + sed_s(k)*DT*onstep(3)
#endif
                do k = ksed1(3), kts, -1
                    odzq = 1./dzq(k)
                    orho = 1./rho(k)
                    qsten(k) = qsten(k) + (sed_s(k+1)-sed_s(k)) &
                        *odzq*onstep(3)*orho
                    rs(k) = max(r1, rs(k) + (sed_s(k+1)-sed_s(k)) &
                        *odzq*DT*onstep(3))
#if !defined(mpas)
                    pfil1(k) = pfil1(k) + sed_s(k)*DT*onstep(3)
#endif
                enddo

                if (rs(kts).gt.R1*1000.) &
                    pptsnow = pptsnow + sed_s(kts)*DT*onstep(3)
            enddo
        endif

        !+---+-----------------------------------------------------------------+

        if (any(l_qg .eqv. .true.)) then
            nstep = nint(1./onstep(4))
            if(.not. sedi_semi) then

                do n = 1, nstep
                    do k = kte, kts, -1
                        sed_g(k) = vtgk(k)*rg(k)
                        sed_n(k) = vtngk(k)*ng(k)
                        sed_b(k) = vtgk(k)*rb(k)
                    enddo
                    k = kte
                    odzq = 1./dzq(k)
                    orho = 1./rho(k)
                    qgten(k) = qgten(k) - sed_g(k)*odzq*onstep(4)*orho
                    ngten(k) = ngten(k) - sed_n(k)*odzq*onstep(4)*orho
                    qbten(k) = qbten(k) - sed_b(k)*odzq*onstep(4)
                    rg(k) = max(r1, rg(k) - sed_g(k)*odzq*dt*onstep(4))
                    ng(k) = max(r2, ng(k) - sed_n(k)*odzq*dt*onstep(4))
                    rb(k) = max(r1/rho(k)/rho_g(nrhg), rb(k) - sed_b(k)*odzq*dt*onstep(4))
#if !defined(mpas)
                    pfil1(k) = pfil1(k) + sed_g(k)*DT*onstep(4)
#endif
                    do k = ksed1(4), kts, -1
                        odzq = 1./dzq(k)
                        orho = 1./rho(k)
                        qgten(k) = qgten(k) + (sed_g(k+1)-sed_g(k)) &
                            *odzq*onstep(4)*orho
                        ngten(k) = ngten(k) + (sed_n(k+1)-sed_n(k)) &
                            *odzq*onstep(4)*orho
                        qbten(k) = qbten(k) + (sed_b(k+1)-sed_b(k)) &
                            *odzq*onstep(4)
                        rg(k) = max(r1, rg(k) + (sed_g(k+1)-sed_g(k)) &
                            *odzq*dt*onstep(4))
                        ng(k) = max(r2, ng(k) + (sed_n(k+1)-sed_n(k)) &
                            *odzq*dt*onstep(4))
                        rb(k) = max(rg(k)/rho(k)/rho_g(nrhg), rb(k) + (sed_b(k+1)-sed_b(k))  &
                            *odzq*DT*onstep(4))
#if !defined(mpas)
                        pfil1(k) = pfil1(k) + sed_g(k)*DT*onstep(4)
#endif
                    enddo

                    if (rg(kts).gt.R1*1000.) &
                        pptgraul = pptgraul + sed_g(kts)*DT*onstep(4)
                enddo
            else ! if(.not. sedi_semi) then
                niter = 1
                dtcfl = dt
                niter = int(nstep/max(decfl_,1)) + 1
                dtcfl = dt/niter

                do n = 1, niter
                    rg_tmp(:) = rg(:)
                    call semi_lagrange_sedim(kte,dzq,vtgk,rg,graulsfc,pfil,dtcfl,R1)
                    do k = kts, kte
                        orhodt = 1./(rho(k)*dt)
                        qgten(k) = qgten(k) + (rg(k) - rg_tmp(k))*orhodt
                        pfil1(k) = pfil1(k) + pfil(k)
                    enddo
                    pptgraul = pptgraul + graulsfc
                    do k = kte+1, kts, -1
                        vtgk(k) = 0.
                    enddo
                    do k = kte, kts, -1
                        vtg = 0.
                        if (rg(k).gt. R1) then
                            ygra1 = alog10(max(1.E-9, rg(k)))
                            rand = 0.0
                            if (present(rand1)) then
                                rand = rand1
                            endif
                            lamg = (am_g(idx_bg(k))*cgg(3,1)*ogg2*ng(k)/rg(k))**obmg

                            ! zans1 = 3.4 + 2./7.*(ygra1+8.) + rand
                            ! ilamg(k) = 1./lamg
                            ! N0_g(k) = ng(k)*ogg2*lamg**cge(2,1)
                
                            ! N0_exp = 10.**(zans1)
                            ! N0_exp = MAX(DBLE(gonv_min), MIN(N0_exp, DBLE(gonv_max)))
                            ! lam_exp = (N0_exp*am_g*cgg(1)/rg(k))**oge1
                            ! lamg = lam_exp * (cgg(3)*ogg2*ogg1)**obmg

                            vtg = rhof(k)*afall*cgg(6,idx_bg(k))*ogg3 * ilamg(k)**bfall
                            ! vtg = rhof(k)*av_g*cgg(6)*ogg3 * (1./lamg)**bv_g
                            if (temp(k).gt. T_0) then
                                vtgk(k) = MAX(vtg, vtrk(k))
                            else
                                vtgk(k) = vtg
                            endif
                        endif
                    enddo
                enddo
            endif ! if(.not. sedi_semi) then
        endif

        !+---+-----------------------------------------------------------------+
        !.. Instantly melt any cloud ice into cloud water if above 0C and
        !.. instantly freeze any cloud water found below HGFR.
        !+---+-----------------------------------------------------------------+
        if (.not. iiwarm) then
            do k = kts, kte
                xri = max(0.0, qi1d(k) + qiten(k)*DT)
                if ( (temp(k).gt. T_0) .and. (xri.gt. 0.0) ) then
                    qcten(k) = qcten(k) + xri*odt
                    ncten(k) = ncten(k) + ni1d(k)*odt
                    qiten(k) = qiten(k) - xri*odt
                    niten(k) = -ni1d(k)*odt
                    tten(k) = tten(k) - lfus*ocp(k)*xri*odt*(1-IFDRY)
                endif

                xrc = max(0.0, qc1d(k) + qcten(k)*DT)
                if ( (temp(k).lt. HGFR) .and. (xrc.gt. 0.0) ) then
                    lfus2 = lsub - lvap(k)
                    xnc = nc1d(k) + ncten(k)*DT
                    qiten(k) = qiten(k) + xrc*odt
                    niten(k) = niten(k) + xnc*odt
                    qcten(k) = qcten(k) - xrc*odt
                    ncten(k) = ncten(k) - xnc*odt
                    tten(k) = tten(k) + lfus2*ocp(k)*xrc*odt*(1-IFDRY)
                endif
            enddo
        endif

        !=================================================================================================================
        !.. All tendencies computed, apply and pass back final values to parent.
        !+---+-----------------------------------------------------------------+
        do k = kts, kte
            t1d(k)  = t1d(k) + tten(k)*DT
            qv1d(k) = max(min_qv, qv1d(k) + qvten(k)*dt)
            qc1d(k) = qc1d(k) + qcten(k)*dt
            nc1d(k) = max(2./rho(k), min(nc1d(k) + ncten(k)*dt, nt_c_max))
            if (configs%aerosol_aware) then
                nwfa1d(k) = max(nwfa_default, min(aero_max, (nwfa1d(k)+nwfaten(k)*dt)))
                nifa1d(k) = max(nifa_default, min(aero_max, (nifa1d(k)+nifaten(k)*dt)))
            endif
            if (qc1d(k) .le. R1) then
                qc1d(k) = 0.0
                nc1d(k) = 0.0
            else
                if (nc1d(k)*rho(k).gt.10000.e6) then
                    nu_c = 2
                elseif (nc1d(k)*rho(k).lt.100.) then
                    nu_c = 15
                else
                    nu_c = nint(nu_c_scale/(nc1d(k)*rho(k))) + 2
                    rand = 0.0
                    if (present(rand2)) then
                        rand = rand2
                    endif
                    nu_c = max(2, min(nu_c+nint(rand), 15))
                endif

                lamc = (am_r*ccg(2,nu_c)*ocg1(nu_c)*nc1d(k)/qc1d(k))**obmr
                xDc = (bm_r + nu_c + 1.) / lamc
                if (xDc.lt. D0c) then
                    lamc = cce(2,nu_c)/D0c
                elseif (xDc.gt. D0r*2.) then
                    lamc = cce(2,nu_c)/(D0r*2.)
                endif
                nc1d(k) = min(ccg(1,nu_c)*ocg2(nu_c)*qc1d(k)/am_r*lamc**bm_r,&
                    real(Nt_c_max, kind=dp)/rho(k))
            endif

            qi1d(k) = qi1d(k) + qiten(k)*DT
            ni1d(k) = max(R2/rho(k), ni1d(k) + niten(k)*DT)
            if (qi1d(k) .le. R1) then
                qi1d(k) = 0.0
                ni1d(k) = 0.0
            else
                lami = (am_i*cig(2)*oig1*ni1d(k)/qi1d(k))**obmi
                ilami = 1./lami
                xDi = (bm_i + mu_i + 1.) * ilami
                if (xDi.lt. 5.E-6) then
                    lami = cie(2)/5.E-6
                elseif (xDi.gt. 300.E-6) then
                    lami = cie(2)/300.E-6
                endif
                ni1d(k) = min(cig(1)*oig2*qi1d(k)/am_i*lami**bm_i, max_ni/rho(k))
            endif
            qr1d(k) = qr1d(k) + qrten(k)*DT
            nr1d(k) = max(R2/rho(k), nr1d(k) + nrten(k)*DT)
            if (qr1d(k) .le. R1) then
                qr1d(k) = 0.0
                nr1d(k) = 0.0
            else
                lamr = (am_r*crg(3)*org2*nr1d(k)/qr1d(k))**obmr
                mvd_r(k) = (3.0 + mu_r + 0.672) / lamr
                if (mvd_r(k) .gt. 2.5E-3) then
                    mvd_r(k) = 2.5E-3
                elseif (mvd_r(k) .lt. D0r*0.75) then
                    mvd_r(k) = D0r*0.75
                endif
                lamr = (3.0 + mu_r + 0.672) / mvd_r(k)
                nr1d(k) = crg(2)*org3*qr1d(k)*lamr**bm_r / am_r
            endif
            qs1d(k) = qs1d(k) + qsten(k)*DT
            if (qs1d(k) .le. R1) qs1d(k) = 0.0
            qg1d(k) = qg1d(k) + qgten(k)*DT
            ng1d(k) = MAX(R2/rho(k), ng1d(k) + ngten(k)*DT)
            if (qg1d(k) .le. R1) then
                qg1d(k) = 0.0
                ng1d(k) = 0.0
                qb1d(k) = 0.0
            else
                qb1d(k) = max(qg1d(k)/rho_g(nrhg), qb1d(k) + qbten(k)*dt)
                qb1d(k) = min(qg1d(k)/rho_g(1), qb1d(k))
                idx_bg(k) = max(1,min(nint(qg1d(k)/qb1d(k) *0.01)+1,nrhg))
                if (.not. configs%hail_aware) idx_bg(k) = idx_bg1
                lamg = (am_g(idx_bg(k))*cgg(3,1)*ogg2*ng1d(k)/qg1d(k))**obmg
                mvd_g(k) = (3.0 + mu_g + 0.672) / lamg
                if (mvd_g(k) .gt. 25.4E-3) then
                    mvd_g(k) = 25.4E-3
                elseif (mvd_g(k) .lt. D0r) then
                    mvd_g(k) = D0r
                endif
                lamg = (3.0 + mu_g + 0.672) / mvd_g(k)
                ng1d(k) = cgg(2,1)*ogg3*qg1d(k)*lamg**bm_g / am_g(idx_bg(k))
            endif

        enddo

#if !defined(mpas)
        ! Diagnostics
        calculate_extended_diagnostics: if (ext_diag) then
            do k = kts, kte
                if(prw_vcd(k).gt.0)then
                    prw_vcdc1(k) = prw_vcd(k)*dt
                elseif(prw_vcd(k).lt.0)then
                    prw_vcde1(k) = -1*prw_vcd(k)*dt
                endif
                !heating/cooling diagnostics
                tpri_inu1(k) = pri_inu(k)*lsub*ocp(k)*orho * (1-IFDRY)*DT

                if(pri_ide(k).gt.0)then
                    tpri_ide1_d(k) = pri_ide(k)*lsub*ocp(k)*orho * (1-IFDRY)*DT
                else
                    tpri_ide1_s(k) = -pri_ide(k)*lsub*ocp(k)*orho * (1-IFDRY)*DT
                endif

                if(temp(k).lt.T_0)then
                    tprs_ide1(k) = prs_ide(k)*lsub*ocp(k)*orho * (1-IFDRY)*DT
                endif

                if(prs_sde(k).gt.0)then
                    tprs_sde1_d(k) = prs_sde(k)*lsub*ocp(k)*orho * (1-IFDRY)*DT
                else
                    tprs_sde1_s(k) = -prs_sde(k)*lsub*ocp(k)*orho * (1-IFDRY)*DT
                endif

                if(prg_gde(k).gt.0)then
                    tprg_gde1_d(k) = prg_gde(k)*lsub*ocp(k)*orho * (1-IFDRY)*DT
                else
                    tprg_gde1_s(k) = -prg_gde(k)*lsub*ocp(k)*orho * (1-IFDRY)*DT
                endif

                tpri_iha1(k) = pri_iha(k)*lsub*ocp(k)*orho * (1-IFDRY)*DT
                tpri_wfz1(k) = pri_wfz(k)*lfus2*ocp(k)*orho * (1-IFDRY)*DT
                tpri_rfz1(k) = pri_rfz(k)*lfus2*ocp(k)*orho * (1-IFDRY)*DT
                tprg_rfz1(k) = prg_rfz(k)*lfus2*ocp(k)*orho * (1-IFDRY)*DT
                tprs_scw1(k) = prs_scw(k)*lfus2*ocp(k)*orho * (1-IFDRY)*DT
                tprg_scw1(k) = prg_scw(k)*lfus2*ocp(k)*orho * (1-IFDRY)*DT
                tprg_rcs1(k) = prg_rcs(k)*lfus2*ocp(k)*orho * (1-IFDRY)*DT

                if(temp(k).lt.T_0)then
                    tprs_rcs1(k) = prs_rcs(k)*lfus2*ocp(k)*orho * (1-IFDRY)*DT
                endif

                tprr_rci1(k) = prr_rci(k)*lfus2*ocp(k)*orho * (1-IFDRY)*DT

                if(temp(k).lt.T_0)then
                    tprg_rcg1(k) = prg_rcg(k)*lfus2*ocp(k)*orho * (1-IFDRY)*DT
                endif

                if(prw_vcd(k).gt.0)then
                    tprw_vcd1_c(k) = lvap(k)*ocp(k)*prw_vcd(k)*(1-IFDRY)*DT
                else
                    tprw_vcd1_e(k) = -lvap(k)*ocp(k)*prw_vcd(k)*(1-IFDRY)*DT
                endif

                ! cooling terms
                tprr_sml1(k) = prr_sml(k)*lfus*ocp(k)*orho * (1-IFDRY)*DT
                tprr_gml1(k) = prr_gml(k)*lfus*ocp(k)*orho * (1-IFDRY)*DT

                if(temp(k).ge.T_0)then
                    tprr_rcg1(k) = -prr_rcg(k)*lfus*ocp(k)*orho * (1-IFDRY)*DT
                endif

                if(temp(k).ge.T_0)then
                    tprr_rcs1(k) = -prr_rcs(k)*lfus*ocp(k)*orho * (1-IFDRY)*DT
                endif

                tprv_rev1(k) = lvap(k)*ocp(k)*prv_rev(k)*(1-IFDRY)*DT
                tten1(k) = tten(k)*DT
                qvten1(k) = qvten(k)*DT
                qiten1(k) = qiten(k)*DT
                qrten1(k) = qrten(k)*DT
                qsten1(k) = qsten(k)*DT
                qgten1(k) = qgten(k)*DT
                niten1(k) = niten(k)*DT
                nrten1(k) = nrten(k)*DT
                ncten1(k) = ncten(k)*DT
                qcten1(k) = qcten(k)*DT
            enddo
        endif calculate_extended_diagnostics
#endif

    end subroutine mp_thompson_main
    !=================================================================================================================
    !..Function to compute collision efficiency of collector species (rain,
    !.. snow, graupel) of aerosols.  Follows Wang et al, 2010, ACP, which
    !.. follows Slinn (1983).
    !+---+-----------------------------------------------------------------+

    real function Eff_aero(D, Da, visc,rhoa,Temp,species)

        implicit none
        real(wp) :: D, Da, visc, rhoa, Temp
        character(LEN=1) :: species
        real(wp) :: aval, Cc, diff, Re, Sc, St, St2, vt, Eff
        real(wp), parameter :: boltzman = 1.3806503E-23
        real(wp), parameter :: meanPath = 0.0256E-6

        vt = 1.
        if (species .eq. 'r') then
            vt = -0.1021 + 4.932E3*D - 0.9551E6*D*D                        &
                + 0.07934E9*D*D*D - 0.002362E12*D*D*D*D
        elseif (species .eq. 's') then
            vt = av_s*D**bv_s
        elseif (species .eq. 'g') then
            vt = av_g(idx_bg1)*D**bv_g(idx_bg1)
        endif

        Cc    = 1. + 2.*meanPath/Da *(1.257+0.4*exp(-0.55*Da/meanPath))
        diff  = boltzman*Temp*Cc/(3.*PI*visc*Da)

        Re    = 0.5*rhoa*D*vt/visc
        Sc    = visc/(rhoa*diff)

        St    = Da*Da*vt*1000./(9.*visc*D)
        aval  = 1.+LOG(1.+Re)
        St2   = (1.2 + 1./12.*aval)/(1.+aval)

        eff = 4./(re*sc) * (1. + 0.4*sqrt(re)*sc**0.3333                  &
            + 0.16*sqrt(re)*sqrt(sc))                  &
            + 4.*da/d * (0.02 + da/d*(1.+2.*sqrt(re)))

        if (St.gt.St2) Eff = Eff  + ( (St-St2)/(St-St2+0.666667))**1.5
        eff_aero = max(1.e-5, min(eff, 1.0))

    end function Eff_aero

    !=================================================================================================================
    !..Retrieve fraction of CCN that gets activated given the model temp,
    !.. vertical velocity, and available CCN concentration.  The lookup
    !.. table (read from external file) has CCN concentration varying the
    !.. quickest, then updraft, then temperature, then mean aerosol radius,
    !.. and finally hygroscopicity, kappa.
    !.. TO_DO ITEM:  For radiation cooling producing fog, in which case the
    !.. updraft velocity could easily be negative, we could use the temp
    !.. and its tendency to diagnose a pretend postive updraft velocity.
    !+---+-----------------------------------------------------------------+
    real function activ_ncloud(Tt, Ww, NCCN,lsm_in)

        implicit none
        REAL, INTENT(IN):: Tt, Ww, NCCN
        INTEGER, INTENT(IN), optional:: lsm_in

        REAL:: n_local, w_local
        INTEGER:: i, j, k, l, m, n
        REAL:: A, B, C, D, t, u, x1, x2, y1, y2, nx, wy, fraction
        REAL:: lower_lim_nuc_frac


        !     ta_Na = (/10.0, 31.6, 100.0, 316.0, 1000.0, 3160.0, 10000.0/)  ntb_arc
        !     ta_Ww = (/0.01, 0.0316, 0.1, 0.316, 1.0, 3.16, 10.0, 31.6, 100.0/)  ntb_arw
        !     ta_Tk = (/243.15, 253.15, 263.15, 273.15, 283.15, 293.15, 303.15/)  ntb_art
        !     ta_Ra = (/0.01, 0.02, 0.04, 0.08, 0.16/)  ntb_arr
        !     ta_Ka = (/0.2, 0.4, 0.6, 0.8/)  ntb_ark

        n_local = NCCN * 1.E-6
        w_local = Ww

        if (n_local .ge. ta_Na(ntb_arc)) then
            n_local = ta_Na(ntb_arc) - 1.0
        elseif (n_local .le. ta_Na(1)) then
            n_local = ta_Na(1) + 1.0
        endif
        do n = 2, ntb_arc
            if (n_local.ge.ta_Na(n-1) .and. n_local.lt.ta_Na(n)) goto 8003
        enddo
8003    continue
        i = n
        x1 = LOG(ta_Na(i-1))
        x2 = LOG(ta_Na(i))

        if (w_local .ge. ta_Ww(ntb_arw)) then
            w_local = ta_Ww(ntb_arw) - 1.0
        elseif (w_local .le. ta_Ww(1)) then
            w_local = ta_Ww(1) + 0.001
        endif
        do n = 2, ntb_arw
            if (w_local.ge.ta_Ww(n-1) .and. w_local.lt.ta_Ww(n)) goto 8005
        enddo
8005    continue
        j = n
        y1 = LOG(ta_Ww(j-1))
        y2 = LOG(ta_Ww(j))

        k = MAX(1, MIN( NINT( (Tt - ta_Tk(1))*0.1) + 1, ntb_art))

        !..The next two values are indexes of mean aerosol radius and
        !.. hygroscopicity.  Currently these are constant but a future version
        !.. should implement other variables to allow more freedom such as
        !.. at least simple separation of tiny size sulfates from larger
        !.. sea salts.
        l = 3
        m = 2

        lower_lim_nuc_frac = 0.
        if (present(lsm_in)) then
            if (lsm_in .eq. 1) then       ! land
                lower_lim_nuc_frac = 0.
            else if (lsm_in .eq. 0) then  ! water
                lower_lim_nuc_frac = 0.15
            else
                lower_lim_nuc_frac = 0.15  ! catch-all for anything else
            endif
        endif
        
        A = tnccn_act(i-1,j-1,k,l,m)
        B = tnccn_act(i,j-1,k,l,m)
        C = tnccn_act(i,j,k,l,m)
        D = tnccn_act(i-1,j,k,l,m)
        nx = LOG(n_local)
        wy = LOG(w_local)

        t = (nx-x1)/(x2-x1)
        u = (wy-y1)/(y2-y1)

        !     t = (n_local-ta(Na(i-1))/(ta_Na(i)-ta_Na(i-1))
        !     u = (w_local-ta_Ww(j-1))/(ta_Ww(j)-ta_Ww(j-1))

        fraction = (1.0-t)*(1.0-u)*A + t*(1.0-u)*B + t*u*C + (1.0-t)*u*D
        fraction = MAX(fraction, lower_lim_nuc_frac)


        !     if (NCCN*fraction .gt. 0.75*Nt_c_max) then
        !        write(*,*) ' DEBUG-GT ', n_local, w_local, Tt, i, j, k
        !     endif

        activ_ncloud = NCCN*fraction

    end function activ_ncloud

    !=================================================================================================================

    !+---+-----------------------------------------------------------------+
    real function iceDeMott(tempc, qv, qvs, qvsi, rho, nifa)
        implicit none

        REAL, INTENT(IN):: tempc, qv, qvs, qvsi, rho, nifa

        !..Local vars
        REAL:: satw, sati, siw, p_x, si0x, dtt, dsi, dsw, dab, fc, hx
        REAL:: ntilde, n_in, nmax, nhat, mux, xni, nifa_cc
        REAL, PARAMETER:: p_c1    = 1000.
        REAL, PARAMETER:: p_rho_c = 0.76
        REAL, PARAMETER:: p_alpha = 1.0
        REAL, PARAMETER:: p_gam   = 2.
        REAL, PARAMETER:: delT    = 5.
        REAL, PARAMETER:: T0x     = -40.
        REAL, PARAMETER:: Sw0x    = 0.97
        REAL, PARAMETER:: delSi   = 0.1
        REAL, PARAMETER:: hdm     = 0.15
        REAL, PARAMETER:: p_psi   = 0.058707*p_gam/p_rho_c
        REAL, PARAMETER:: aap     = 1.
        REAL, PARAMETER:: bbp     = 0.
        REAL, PARAMETER:: y1p     = -35.
        REAL, PARAMETER:: y2p     = -25.
        REAL, PARAMETER:: rho_not0 = 101325./(287.05*273.15)

        !+---+

        xni = 0.0
        !     satw = qv/qvs
        !     sati = qv/qvsi
        !     siw = qvs/qvsi
        !     p_x = -1.0261+(3.1656e-3*tempc)+(5.3938e-4*(tempc*tempc))         &
        !                +  (8.2584e-6*(tempc*tempc*tempc))
        !     si0x = 1.+(10.**p_x)
        !     if (sati.ge.si0x .and. satw.lt.0.985) then
        !        dtt = delta_p (tempc, T0x, T0x+delT, 1., hdm)
        !        dsi = delta_p (sati, Si0x, Si0x+delSi, 0., 1.)
        !        dsw = delta_p (satw, Sw0x, 1., 0., 1.)
        !        fc = dtt*dsi*0.5
        !        hx = min(fc+((1.-fc)*dsw), 1.)
        !        ntilde = p_c1*p_gam*((exp(12.96*(sati-1.1)))**0.3) / p_rho_c
        !        if (tempc .le. y1p) then
        !           n_in = ntilde
        !        elseif (tempc .ge. y2p) then
        !           n_in = p_psi*p_c1*exp(12.96*(sati-1.)-0.639)
        !        else
        !           if (tempc .le. -30.) then
        !              nmax = p_c1*p_gam*(exp(12.96*(siw-1.1)))**0.3/p_rho_c
        !           else
        !              nmax = p_psi*p_c1*exp(12.96*(siw-1.)-0.639)
        !           endif
        !           ntilde = MIN(ntilde, nmax)
        !           nhat = MIN(p_psi*p_c1*exp(12.96*(sati-1.)-0.639), nmax)
        !           dab = delta_p (tempc, y1p, y2p, aap, bbp)
        !           n_in = MIN(nhat*(ntilde/nhat)**dab, nmax)
        !        endif
        !        mux = hx*p_alpha*n_in*rho
        !        xni = mux*((6700.*nifa)-200.)/((6700.*5.E5)-200.)
        !     elseif (satw.ge.0.985 .and. tempc.gt.HGFR-273.15) then
        nifa_cc = MAX(0.5, nifa*RHO_NOT0*1.E-6/rho)
        !        xni  = 3.*nifa_cc**(1.25)*exp((0.46*(-tempc))-11.6)              !  [DeMott, 2015]
        xni = (5.94e-5*(-tempc)**3.33)                                 & !  [DeMott, 2010]
            * (nifa_cc**((-0.0264*(tempc))+0.0033))
        xni = xni*rho/RHO_NOT0 * 1000.
        !     endif

        iceDeMott = MAX(0., xni)

    end FUNCTION iceDeMott

    !+---+-----------------------------------------------------------------+
    !..Newer research since Koop et al (2001) suggests that the freezing
    !.. rate should be lower than original paper, so J_rate is reduced
    !.. by two orders of magnitude.

    real function iceKoop(temp, qv, qvs, naero, dt)
        implicit none

        REAL, INTENT(IN):: temp, qv, qvs, naero, DT
        REAL:: mu_diff, a_w_i, delta_aw, log_J_rate, J_rate, prob_h, satw
        REAL:: xni

        xni = 0.0
        satw = qv/qvs
        mu_diff    = 210368.0 + (131.438*temp) - (3.32373E6/temp)         &
        &           - (41729.1*alog(temp))
        a_w_i      = exp(mu_diff/(R_uni*temp))
        delta_aw   = satw - a_w_i
        log_J_rate = -906.7 + (8502.0*delta_aw)                           &
        &           - (26924.0*delta_aw*delta_aw)                          &
        &           + (29180.0*delta_aw*delta_aw*delta_aw)
        log_J_rate = MIN(20.0, log_J_rate)
        J_rate     = 10.**log_J_rate                                       ! cm-3 s-1
        prob_h     = MIN(1.-exp(-J_rate*ar_volume*DT), 1.)
        if (prob_h .gt. 0.) then
            xni     = MIN(prob_h*naero, 1000.E3)
        endif

        iceKoop = MAX(0.0, xni)

    end FUNCTION iceKoop

    !+---+-----------------------------------------------------------------+
    !.. Helper routine for Phillips et al (2008) ice nucleation.  Trude

    REAL FUNCTION delta_p (yy, y1, y2, aa, bb)
        IMPLICIT NONE

        REAL, INTENT(IN):: yy, y1, y2, aa, bb
        REAL:: dab, A, B, a0, a1, a2, a3

        A   = 6.*(aa-bb)/((y2-y1)*(y2-y1)*(y2-y1))
        B   = aa+(A*y1*y1*y1/6.)-(A*y1*y1*y2*0.5)
        a0  = B
        a1  = A*y1*y2
        a2  = -A*(y1+y2)*0.5
        a3  = A/3.

        if (yy.le.y1) then
            dab = aa
        else if (yy.ge.y2) then
            dab = bb
        else
            dab = a0+(a1*yy)+(a2*yy*yy)+(a3*yy*yy*yy)
        endif

        if (dab.lt.aa) then
            dab = aa
        endif
        if (dab.gt.bb) then
            dab = bb
        endif
        delta_p = dab

    END FUNCTION delta_p

    !+---+-----------------------------------------------------------------+
!-------------------------------------------------------------------
      SUBROUTINE semi_lagrange_sedim(km,dzl,wwl,rql,precip,pfsan,dt,R1)
!-------------------------------------------------------------------
!
! This routine is a semi-Lagrangain forward advection for hydrometeors
! with mass conservation and positive definite advection
! 2nd order interpolation with monotonic piecewise parabolic method is used.
! This routine is under assumption of decfl < 1 for semi_Lagrangian
!
! dzl    depth of model layer in meter
! wwl    terminal velocity at model layer m/s
! rql    dry air density*mixing ratio
! precip precipitation at surface 
! dt     time step
!
! author: hann-ming henry juang <henry.juang@noaa.gov>
!         implemented by song-you hong
! reference: Juang, H.-M., and S.-Y. Hong, 2010: Forward semi-Lagrangian advection
!         with mass conservation and positive definiteness for falling
!         hydrometeors. *Mon.  Wea. Rev.*, *138*, 1778-1791
!
      implicit none

      integer, intent(in) :: km
      real, intent(in) ::  dt, R1
      real, intent(in) :: dzl(km),wwl(km)
      real, intent(out) :: precip
      real, intent(inout) :: rql(km)
      real, intent(out)  :: pfsan(km)
      integer  k,m,kk,kb,kt
      real  tl,tl2,qql,dql,qqd
      real  th,th2,qqh,dqh
      real  zsum,qsum,dim,dip,con1,fa1,fa2
      real  allold, decfl
      real  dz(km), ww(km), qq(km)
      real  wi(km+1), zi(km+1), za(km+2)
      real  qn(km)
      real  dza(km+1), qa(km+1), qmi(km+1), qpi(km+1)
      real  net_flx(km)
!
      precip = 0.0
      qa(:) = 0.0
      qq(:) = 0.0
      dz(:) = dzl(:)
      ww(:) = wwl(:)
      do k = 1,km
        if(rql(k).gt.R1) then 
          qq(k) = rql(k) 
        else 
          ww(k) = 0.0 
        endif
        pfsan(k) = 0.0
        net_flx(k) = 0.0
      enddo
! skip for no precipitation for all layers
      allold = 0.0
      do k=1,km
        allold = allold + qq(k)
      enddo
      if(allold.le.0.0) then
         return 
      endif
!
! compute interface values
      zi(1)=0.0
      do k=1,km
        zi(k+1) = zi(k)+dz(k)
      enddo
!     n=1
! plm is 2nd order, we can use 2nd order wi or 3rd order wi
! 2nd order interpolation to get wi
      wi(1) = ww(1)
      wi(km+1) = ww(km)
      do k=2,km
        wi(k) = (ww(k)*dz(k-1)+ww(k-1)*dz(k))/(dz(k-1)+dz(k))
      enddo
! 3rd order interpolation to get wi
      fa1 = 9./16.
      fa2 = 1./16.
      wi(1) = ww(1)
      wi(2) = 0.5*(ww(2)+ww(1))
      do k=3,km-1
        wi(k) = fa1*(ww(k)+ww(k-1))-fa2*(ww(k+1)+ww(k-2))
      enddo
      wi(km) = 0.5*(ww(km)+ww(km-1))
      wi(km+1) = ww(km)

! terminate of top of raingroup
      do k=2,km
        if( ww(k).eq.0.0 ) wi(k)=ww(k-1)
      enddo

! diffusivity of wi
      con1 = 0.05
      do k=km,1,-1
        decfl = (wi(k+1)-wi(k))*dt/dz(k)
        if( decfl .gt. con1 ) then
          wi(k) = wi(k+1) - con1*dz(k)/dt
        endif
      enddo
! compute arrival point
      do k=1,km+1
        za(k) = zi(k) - wi(k)*dt
      enddo
      za(km+2) = zi(km+1)

      do k=1,km+1
        dza(k) = za(k+1)-za(k)
      enddo

! computer deformation at arrival point
      do k=1,km
        qa(k) = qq(k)*dz(k)/dza(k)
      enddo
      qa(km+1) = 0.0

! estimate values at arrival cell interface with monotone
      do k=2,km
        dip=(qa(k+1)-qa(k))/(dza(k+1)+dza(k))
        dim=(qa(k)-qa(k-1))/(dza(k-1)+dza(k))
        if( dip*dim.le.0.0 ) then
          qmi(k)=qa(k)
          qpi(k)=qa(k)
        else
          qpi(k)=qa(k)+0.5*(dip+dim)*dza(k)
          qmi(k)=2.0*qa(k)-qpi(k)
          if( qpi(k).lt.0.0 .or. qmi(k).lt.0.0 ) then
            qpi(k) = qa(k)
            qmi(k) = qa(k)
          endif
        endif
      enddo
      qpi(1)=qa(1)
      qmi(1)=qa(1)
      qmi(km+1)=qa(km+1)
      qpi(km+1)=qa(km+1)

! interpolation to regular point
      qn = 0.0
      kb=1
      kt=1
      intp : do k=1,km
             kb=max(kb-1,1)
             kt=max(kt-1,1)
! find kb and kt
             if( zi(k).ge.za(km+1) ) then
               exit intp
             else
               find_kb : do kk=kb,km
                         if( zi(k).le.za(kk+1) ) then
                           kb = kk
                           exit find_kb
                         else
                           cycle find_kb
                         endif
               enddo find_kb
               find_kt : do kk=kt,km+2
                         if( zi(k+1).le.za(kk) ) then
                           kt = kk
                           exit find_kt
                         else
                           cycle find_kt
                         endif
               enddo find_kt
               kt = kt - 1
! compute q with piecewise constant method
               if( kt.eq.kb ) then
                 tl=(zi(k)-za(kb))/dza(kb)
                 th=(zi(k+1)-za(kb))/dza(kb)
                 tl2=tl*tl
                 th2=th*th
                 qqd=0.5*(qpi(kb)-qmi(kb))
                 qqh=qqd*th2+qmi(kb)*th
                 qql=qqd*tl2+qmi(kb)*tl
                 qn(k) = (qqh-qql)/(th-tl)
               else if( kt.gt.kb ) then
                 tl=(zi(k)-za(kb))/dza(kb)
                 tl2=tl*tl
                 qqd=0.5*(qpi(kb)-qmi(kb))
                 qql=qqd*tl2+qmi(kb)*tl
                 dql = qa(kb)-qql
                 zsum  = (1.-tl)*dza(kb)
                 qsum  = dql*dza(kb)
                 if( kt-kb.gt.1 ) then
                 do m=kb+1,kt-1
                   zsum = zsum + dza(m)
                   qsum = qsum + qa(m) * dza(m)
                 enddo
                 endif
                 th=(zi(k+1)-za(kt))/dza(kt)
                 th2=th*th
                 qqd=0.5*(qpi(kt)-qmi(kt))
                 dqh=qqd*th2+qmi(kt)*th
                 zsum  = zsum + th*dza(kt)
                 qsum  = qsum + dqh*dza(kt)
                 qn(k) = qsum/zsum
               endif
               cycle intp
             endif

       enddo intp

! rain out
      sum_precip: do k=1,km
                    if( za(k).lt.0.0 .and. za(k+1).le.0.0 ) then
                      precip = precip + qa(k)*dza(k)
                      net_flx(k) =  qa(k)*dza(k)
                      cycle sum_precip
                    else if ( za(k).lt.0.0 .and. za(k+1).gt.0.0 ) then
                      th = (0.0-za(k))/dza(k)
                      th2 = th*th
                      qqd = 0.5*(qpi(k)-qmi(k))
                      qqh = qqd*th2+qmi(k)*th
                      precip = precip + qqh*dza(k)
                      net_flx(k) = qqh*dza(k)
                      exit sum_precip
                    endif
                    exit sum_precip
      enddo sum_precip

! calculating precipitation fluxes
      do k=km,1,-1
         if(k == km) then
           pfsan(k) = net_flx(k)
         else
           pfsan(k) = pfsan(k+1) + net_flx(k)
         end if
      enddo
!
! replace the new values
      rql(:) = max(qn(:),R1)

  END SUBROUTINE semi_lagrange_sedim

  !------------------
  
end module module_mp_thompson_main
