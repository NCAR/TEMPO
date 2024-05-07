! 3D Thompson-Eidhammer Microphysics Driver for MPAS
!=================================================================================================================
module module_mp_thompson

    use mpas_kind_types, only: wp => RKIND, sp => R4KIND, dp => R8KIND

    use module_mp_thompson_params
    use module_mp_thompson_main
    use module_mp_thompson_utils
    use mpas_atmphys_functions, only: gammp, wgamma, rslf, rsif
    use mpas_atmphys_utilities
    use mpas_io_units, only : mpas_new_unit, mpas_release_unit
    use mp_radar, only : radar_init

    type(config_flags) configs

contains
!=================================================================================================================
! This subroutine handles initialzation of the microphysics scheme including building of lookup tables,
! allocating arrays for the microphysics scheme, and defining gamma function variables. 

! Input:
!   l_mp_tables = .false. to build lookup tables. If l_mp_tables = .true., lookup tables are not built.

    ! AAJ No support yet for hail_aware in microphysics driver
    subroutine thompson_init(l_mp_tables, hail_aware_flag, aerosol_aware_flag)
    
        implicit none

! Input arguments:
        logical, intent(in) :: l_mp_tables, hail_aware_flag
        logical, optional, intent(in) :: aerosol_aware_flag

        integer, parameter :: open_OK = 0
        integer, parameter :: num_records = 5
        integer :: qr_acr_qg_filesize, qr_acr_qg_check, qr_acr_qg_dim1size, qr_acr_qg_dim9size
        logical :: qr_acr_qg_exists
        integer :: i, j, k, l, m, n
        integer :: istat
        logical :: micro_init
        integer :: mp_unit
        character(len=132) :: message

        ! If lookup tables are already built
        if (l_mp_tables) then
           configs%hail_aware = hail_aware_flag
           write(message, '(L1)') configs%hail_aware
           call physics_message('--- thompson_init() called with hail_aware_flag = ' // trim(message))

           if (present(aerosol_aware_flag)) then
              configs%aerosol_aware = aerosol_aware_flag
              write(message, '(L1)') configs%aerosol_aware
              call physics_message('--- thompson_init() called with aerosol_aware_flag = ' // trim(message))
           endif
        endif
        
! Allocate space for lookup tables (J. Michalakes 2009Jun08).
        if (hail_aware_flag) then
            dimNRHG = NRHG
        else
            av_g(idx_bg1) = av_g_old
            bv_g(idx_bg1) = bv_g_old
            dimNRHG = NRHG1
        endif

        micro_init = .false.

!=================================================================================================================
! Check the qr_acr_qg lookup table to make sure it is compatible with runtime options

        ! If lookup tables are already built
        if (l_mp_tables) then

           inquire(file='MP_THOMPSON_QRacrQG_DATA.DBL', exist=qr_acr_qg_exists)
           if (qr_acr_qg_exists) then ! Check again that file exists

! Add check on qr_ac_qg filesize to determine if table includes hail-awareness (dimNRHG=9)
              qr_acr_qg_check = dp * num_records * (dimNRHG * ntb_g1 * ntb_g * ntb_r1 * ntb_r + 1)
              qr_acr_qg_dim1size = dp * num_records * (NRHG1 * ntb_g1 * ntb_g * ntb_r1 * ntb_r + 1)
              qr_acr_qg_dim9size = dp * num_records * (NRHG * ntb_g1 * ntb_g * ntb_r1 * ntb_r + 1)

              inquire(file='MP_THOMPSON_QRacrQG_DATA.DBL', size=qr_acr_qg_filesize)

              if (qr_acr_qg_filesize == qr_acr_qg_dim1size) then
                 using_hail_aware_table = .false.
                 call physics_message('--- thompson_init() ' // &
                      'Lookup table for qr_acr_qg is not hail aware.')
                 if (hail_aware_flag) then
                    call physics_error_fatal('--- thompson_init() Cannot use hail-aware microphysics ' // &
                         'with non hail-aware qr_acr_qg lookup table. ' // &
                         'Please rebuild table with parameter build_hail_aware_table set to true.')
                 endif
              elseif (qr_acr_qg_filesize == qr_acr_qg_dim9size) then
                 using_hail_aware_table = .true.
                 call physics_message('--- thompson_init() ' // &
                      'Lookup table for qr_acr_qg is hail aware.')
              else
                 using_hail_aware_table = .false.
                 if (hail_aware_flag) using_hail_aware_table = .true.
                 call physics_message('--- thompson_init() ' // &
                      'Could not determine if lookup table for qr_acr_qg is hail aware based on file size.')
              endif
           endif
        endif
!=================================================================================================================
! Allocate space for lookup tables (J. Michalakes 2009Jun08).
        if (.not. allocated(tcg_racg)) then
            allocate(tcg_racg(ntb_g1,ntb_g,dimNRHG,ntb_r1,ntb_r))
            micro_init = .true.
        endif

        if (.not. allocated(tmr_racg)) allocate(tmr_racg(ntb_g1,ntb_g,dimNRHG,ntb_r1,ntb_r))
        if (.not. allocated(tcr_gacr)) allocate(tcr_gacr(ntb_g1,ntb_g,dimNRHG,ntb_r1,ntb_r))
        if (.not. allocated(tnr_racg)) allocate(tnr_racg(ntb_g1,ntb_g,dimNRHG,ntb_r1,ntb_r))
        if (.not. allocated(tnr_gacr)) allocate(tnr_gacr(ntb_g1,ntb_g,dimNRHG,ntb_r1,ntb_r))

        if (.not. allocated(tcs_racs1)) allocate(tcs_racs1(ntb_s,ntb_t,ntb_r1,ntb_r))
        if (.not. allocated(tmr_racs1)) allocate(tmr_racs1(ntb_s,ntb_t,ntb_r1,ntb_r))
        if (.not. allocated(tcs_racs2)) allocate(tcs_racs2(ntb_s,ntb_t,ntb_r1,ntb_r))
        if (.not. allocated(tmr_racs2)) allocate(tmr_racs2(ntb_s,ntb_t,ntb_r1,ntb_r))
        if (.not. allocated(tcr_sacr1)) allocate(tcr_sacr1(ntb_s,ntb_t,ntb_r1,ntb_r))
        if (.not. allocated(tms_sacr1)) allocate(tms_sacr1(ntb_s,ntb_t,ntb_r1,ntb_r))
        if (.not. allocated(tcr_sacr2)) allocate(tcr_sacr2(ntb_s,ntb_t,ntb_r1,ntb_r))
        if (.not. allocated(tms_sacr2)) allocate(tms_sacr2(ntb_s,ntb_t,ntb_r1,ntb_r))
        if (.not. allocated(tnr_racs1)) allocate(tnr_racs1(ntb_s,ntb_t,ntb_r1,ntb_r))
        if (.not. allocated(tnr_racs2)) allocate(tnr_racs2(ntb_s,ntb_t,ntb_r1,ntb_r))
        if (.not. allocated(tnr_sacr1)) allocate(tnr_sacr1(ntb_s,ntb_t,ntb_r1,ntb_r))
        if (.not. allocated(tnr_sacr2)) allocate(tnr_sacr2(ntb_s,ntb_t,ntb_r1,ntb_r))

        if (.not. allocated(tpi_qcfz)) allocate(tpi_qcfz(ntb_c,nbc,45,ntb_in))
        if (.not. allocated(tni_qcfz)) allocate(tni_qcfz(ntb_c,nbc,45,ntb_in))

        if (.not. allocated(tpi_qrfz)) allocate(tpi_qrfz(ntb_r,ntb_r1,45,ntb_in))
        if (.not. allocated(tpg_qrfz)) allocate(tpg_qrfz(ntb_r,ntb_r1,45,ntb_in))
        if (.not. allocated(tni_qrfz)) allocate(tni_qrfz(ntb_r,ntb_r1,45,ntb_in))
        if (.not. allocated(tnr_qrfz)) allocate(tnr_qrfz(ntb_r,ntb_r1,45,ntb_in))

        if (.not. allocated(tps_iaus)) allocate(tps_iaus(ntb_i,ntb_i1))
        if (.not. allocated(tni_iaus)) allocate(tni_iaus(ntb_i,ntb_i1))
        if (.not. allocated(tpi_ide)) allocate(tpi_ide(ntb_i,ntb_i1))

        if (.not. allocated(t_efrw)) allocate(t_efrw(nbr,nbc))
        if (.not. allocated(t_efsw)) allocate(t_efsw(nbs,nbc))

        if (.not. allocated(tnr_rev)) allocate(tnr_rev(nbr, ntb_r1, ntb_r))
        if (.not. allocated(tpc_wev)) allocate(tpc_wev(nbc,ntb_c,nbc))
        if (.not. allocated(tnc_wev)) allocate(tnc_wev(nbc,ntb_c,nbc))

        if (.not. allocated(tnccn_act)) allocate(tnccn_act(ntb_arc,ntb_arw,ntb_art,ntb_arr,ntb_ark))

!=================================================================================================================
        if (micro_init) then

! Schmidt number to one-third used numerous times.
            Sc3 = Sc**(1./3.)

! Compute min ice diameter from mass and minimum snow/graupel mass from diameter
            D0i = (xm0i/am_i)**(1.0/bm_i)
            xm0s = am_s * D0s**bm_s
            xm0g = am_g(NRHG) * D0g**bm_g

! These constants various exponents and gamma() assoc with cloud, rain, snow, and graupel.
            do n = 1, 15
                cce(1,n) = n + 1.
                cce(2,n) = bm_r + n + 1.
                cce(3,n) = bm_r + n + 4.
                cce(4,n) = n + bv_c + 1.
                cce(5,n) = bm_r + n + bv_c + 1.
                ccg(1,n) = WGAMMA(cce(1,n))
                ccg(2,n) = WGAMMA(cce(2,n))
                ccg(3,n) = WGAMMA(cce(3,n))
                ccg(4,n) = WGAMMA(cce(4,n))
                ccg(5,n) = WGAMMA(cce(5,n))
                ocg1(n) = 1./ccg(1,n)
                ocg2(n) = 1./ccg(2,n)
            enddo

            cie(1) = mu_i + 1.
            cie(2) = bm_i + mu_i + 1.
            cie(3) = bm_i + mu_i + bv_i + 1.
            cie(4) = mu_i + bv_i + 1.
            cie(5) = mu_i + 2.
            cie(6) = bm_i*0.5 + mu_i + bv_i + 1.
            cie(7) = bm_i*0.5 + mu_i + 1.
            cig(1) = WGAMMA(cie(1))
            cig(2) = WGAMMA(cie(2))
            cig(3) = WGAMMA(cie(3))
            cig(4) = WGAMMA(cie(4))
            cig(5) = WGAMMA(cie(5))
            cig(6) = WGAMMA(cie(6))
            cig(7) = WGAMMA(cie(7))
            oig1 = 1./cig(1)
            oig2 = 1./cig(2)
            obmi = 1./bm_i

            cre(1) = bm_r + 1.
            cre(2) = mu_r + 1.
            cre(3) = bm_r + mu_r + 1.
            cre(4) = bm_r*2. + mu_r + 1.
            cre(5) = mu_r + bv_r + 1.
            cre(6) = bm_r + mu_r + bv_r + 1.
            cre(7) = bm_r*0.5 + mu_r + bv_r + 1.
            cre(8) = bm_r + mu_r + bv_r + 3.
            cre(9) = mu_r + bv_r + 3.
            cre(10) = mu_r + 2.
            cre(11) = 0.5*(bv_r + 5. + 2.*mu_r)
            cre(12) = bm_r*0.5 + mu_r + 1.
            cre(13) = bm_r*2. + mu_r + bv_r + 1.

            do n = 1, 13
                crg(n) = WGAMMA(cre(n))
            enddo

            obmr = 1./bm_r
            ore1 = 1./cre(1)
            org1 = 1./crg(1)
            org2 = 1./crg(2)
            org3 = 1./crg(3)

            cse(1) = bm_s + 1.
            cse(2) = bm_s + 2.
            cse(3) = bm_s*2.
            cse(4) = bm_s + bv_s + 1.
            cse(5) = bm_s*2. + bv_s + 1.
            cse(6) = bm_s*2. + 1.
            cse(7) = bm_s + mu_s + 1.
            cse(8) = bm_s + mu_s + 2.
            cse(9) = bm_s + mu_s + 3.
            cse(10) = bm_s + mu_s + bv_s + 1.
            cse(11) = bm_s*2. + mu_s + bv_s + 1.
            cse(12) = bm_s*2. + mu_s + 1.
            cse(13) = bv_s + 2.
            cse(14) = bm_s + bv_s
            cse(15) = mu_s + 1.
            cse(16) = 1.0 + (1.0 + bv_s)/2.
            cse(17) = bm_s + bv_s + 2.

            do n = 1, 17
                csg(n) = WGAMMA(cse(n))
            enddo

            oams = 1./am_s
            obms = 1./bm_s
            ocms = oams**obms

            cge(1,:) = bm_g + 1.
            cge(2,:) = mu_g + 1.
            cge(3,:) = bm_g + mu_g + 1.
            cge(4,:) = bm_g*2. + mu_g + 1.
            cge(10,:) = mu_g + 2.
            cge(12,:) = bm_g*0.5 + mu_g + 1.

            do m = 1, NRHG
                cge(5,m) = bm_g*2. + mu_g + bv_g(m) + 1.
                cge(6,m) = bm_g + mu_g + bv_g(m) + 1.
                cge(7,m) = bm_g*0.5 + mu_g + bv_g(m) + 1.
                cge(8,m) = mu_g + bv_g(m) + 1.      ! not used
                cge(9,m) = mu_g + bv_g(m) + 3.
                cge(11,m) = 0.5*(bv_g(m) + 5. + 2.*mu_g)
            enddo

            do m = 1, NRHG
                do n = 1, 12
                    cgg(n,m) = WGAMMA(cge(n,m))
                enddo
            enddo

            oamg = 1./am_g
            obmg = 1./bm_g

            do m = 1, NRHG
                oamg(m) = 1./am_g(m)
                ocmg(m) = oamg(m)**obmg
            enddo

            oge1 = 1./cge(1,1)
            ogg1 = 1./cgg(1,1)
            ogg2 = 1./cgg(2,1)
            ogg3 = 1./cgg(3,1)

!=================================================================================================================
! Simplify various rate eqns the best we can now.

! Rain collecting cloud water and cloud ice
            t1_qr_qc = PI * 0.25 * av_r * crg(9)
            t1_qr_qi = PI * 0.25 * av_r * crg(9)
            t2_qr_qi = PI * 0.25 * am_r*av_r * crg(8)

! Graupel collecting cloud water
!     t1_qg_qc = PI*.25*av_g * cgg(9)

! Snow collecting cloud water
            t1_qs_qc = PI * 0.25 * av_s

! Snow collecting cloud ice
            t1_qs_qi = PI * 0.25 * av_s

! Evaporation of rain; ignore depositional growth of rain.
            t1_qr_ev = 0.78 * crg(10)
            t2_qr_ev = 0.308*Sc3*SQRT(av_r) * crg(11)

!..Sublimation/depositional growth of snow
            t1_qs_sd = 0.86
            t2_qs_sd = 0.28*Sc3*SQRT(av_s)

!..Melting of snow
            t1_qs_me = PI*4.*C_sqrd*olfus * 0.86
            t2_qs_me = PI*4.*C_sqrd*olfus * 0.28*Sc3*SQRT(av_s)

!..Sublimation/depositional growth of graupel
            t1_qg_sd = 0.86 * cgg(10,1)
!     t2_qg_sd = 0.28*Sc3*SQRT(av_g) * cgg(11)

!..Melting of graupel
            t1_qg_me = PI*4.*C_cube*olfus * 0.86 * cgg(10,1)
!     t2_qg_me = PI*4.*C_cube*olfus * 0.28*Sc3*SQRT(av_g) * cgg(11)

!..Constants for helping find lookup table indexes.
            nic2 = nint(alog10(r_c(1)))
            nii2 = nint(alog10(r_i(1)))
            nii3 = nint(alog10(Nt_i(1)))
            nir2 = nint(alog10(r_r(1)))
            nir3 = nint(alog10(N0r_exp(1)))
            nis2 = nint(alog10(r_s(1)))
            nig2 = nint(alog10(r_g(1)))
            nig3 = nint(alog10(N0g_exp(1)))
            niIN2 = nint(alog10(Nt_IN(1)))

!..Create bins of cloud water (from min diameter up to 100 microns).
            Dc(1) = D0c*1.0d0
            dtc(1) = D0c*1.0d0
            do n = 2, nbc
                Dc(n) = Dc(n-1) + 1.0D-6
                dtc(n) = (Dc(n) - Dc(n-1))
            enddo

!..Create bins of cloud ice (from min diameter up to 5x min snow size).
            xDx(1) = D0i*1.0d0
            xDx(nbi+1) = 2.0d0*D0s
            do n = 2, nbi
                xDx(n) = DEXP(REAL(n-1,KIND=dp)/REAL(nbi,KIND=dp) & 
                        *DLOG(xDx(nbi+1)/xDx(1)) +DLOG(xDx(1)))
            enddo
            do n = 1, nbi
                Di(n) = DSQRT(xDx(n)*xDx(n+1))
                dti(n) = xDx(n+1) - xDx(n)
            enddo

!..Create bins of rain (from min diameter up to 5 mm).
            xDx(1) = D0r*1.0d0
            xDx(nbr+1) = 0.005d0
            do n = 2, nbr
                xDx(n) = DEXP(REAL(n-1,KIND=dp)/REAL(nbr,KIND=dp) &
                    *DLOG(xDx(nbr+1)/xDx(1)) +DLOG(xDx(1)))
            enddo
            do n = 1, nbr
                Dr(n) = DSQRT(xDx(n)*xDx(n+1))
                dtr(n) = xDx(n+1) - xDx(n)
            enddo

!..Create bins of snow (from min diameter up to 2 cm).
            xDx(1) = D0s*1.0d0
            xDx(nbs+1) = 0.02d0
            do n = 2, nbs
                xDx(n) = DEXP(REAL(n-1,KIND=dp)/REAL(nbs,KIND=dp) &
                    *DLOG(xDx(nbs+1)/xDx(1)) +DLOG(xDx(1)))
            enddo
            do n = 1, nbs
                Ds(n) = DSQRT(xDx(n)*xDx(n+1))
                dts(n) = xDx(n+1) - xDx(n)
            enddo

!..Create bins of graupel (from min diameter up to 5 cm).
            xDx(1) = D0g*1.0d0
            xDx(nbg+1) = 0.05d0
            do n = 2, nbg
                xDx(n) = DEXP(REAL(n-1,KIND=dp)/REAL(nbg,KIND=dp) &
                    *DLOG(xDx(nbg+1)/xDx(1)) +DLOG(xDx(1)))
            enddo
            do n = 1, nbg
                Dg(n) = DSQRT(xDx(n)*xDx(n+1))
                dtg(n) = xDx(n+1) - xDx(n)
            enddo

!..Create bins of cloud droplet number concentration (1 to 3000 per cc).
            xDx(1) = 1.0d0
            xDx(nbc+1) = 3000.0d0
            do n = 2, nbc
                xDx(n) = DEXP(REAL(n-1,KIND=dp)/REAL(nbc,KIND=dp) &
                    *DLOG(xDx(nbc+1)/xDx(1)) +DLOG(xDx(1)))
            enddo
            do n = 1, nbc
                t_Nc(n) = DSQRT(xDx(n)*xDx(n+1)) * 1.D6
            enddo
            nic1 = DLOG(t_Nc(nbc)/t_Nc(1))

!+---+-----------------------------------------------------------------+
!..Create lookup tables for most costly calculations.
!+---+-----------------------------------------------------------------+

            do m = 1, ntb_r
                do k = 1, ntb_r1
                    do n = 1, dimNRHG
                        do j = 1, ntb_g
                            do i = 1, ntb_g1
                                tcg_racg(i,j,n,k,m) = 0.0d0
                                tmr_racg(i,j,n,k,m) = 0.0d0
                                tcr_gacr(i,j,n,k,m) = 0.0d0
                                tnr_racg(i,j,n,k,m) = 0.0d0
                                tnr_gacr(i,j,n,k,m) = 0.0d0
                            enddo
                        enddo
                    enddo
                enddo
            enddo

            do m = 1, ntb_r
                do k = 1, ntb_r1
                    do j = 1, ntb_t
                        do i = 1, ntb_s
                            tcs_racs1(i,j,k,m) = 0.0d0
                            tmr_racs1(i,j,k,m) = 0.0d0
                            tcs_racs2(i,j,k,m) = 0.0d0
                            tmr_racs2(i,j,k,m) = 0.0d0
                            tcr_sacr1(i,j,k,m) = 0.0d0
                            tms_sacr1(i,j,k,m) = 0.0d0
                            tcr_sacr2(i,j,k,m) = 0.0d0
                            tms_sacr2(i,j,k,m) = 0.0d0
                            tnr_racs1(i,j,k,m) = 0.0d0
                            tnr_racs2(i,j,k,m) = 0.0d0
                            tnr_sacr1(i,j,k,m) = 0.0d0
                            tnr_sacr2(i,j,k,m) = 0.0d0
                        enddo
                    enddo
                enddo
            enddo

            do m = 1, ntb_IN
                do k = 1, 45
                    do j = 1, ntb_r1
                        do i = 1, ntb_r
                            tpi_qrfz(i,j,k,m) = 0.0d0
                            tni_qrfz(i,j,k,m) = 0.0d0
                            tpg_qrfz(i,j,k,m) = 0.0d0
                            tnr_qrfz(i,j,k,m) = 0.0d0
                        enddo
                    enddo
                    do j = 1, nbc
                        do i = 1, ntb_c
                            tpi_qcfz(i,j,k,m) = 0.0d0
                            tni_qcfz(i,j,k,m) = 0.0d0
                        enddo
                    enddo
                enddo
            enddo

            do j = 1, ntb_i1
                do i = 1, ntb_i
                    tps_iaus(i,j) = 0.0d0
                    tni_iaus(i,j) = 0.0d0
                    tpi_ide(i,j) = 0.0d0
                enddo
            enddo

            do j = 1, nbc
                do i = 1, nbr
                    t_Efrw(i,j) = 0.0
                enddo
                do i = 1, nbs
                    t_Efsw(i,j) = 0.0
                enddo
            enddo

            do k = 1, ntb_r
                do j = 1, ntb_r1
                    do i = 1, nbr
                        tnr_rev(i,j,k) = 0.0d0
                    enddo
                enddo
            enddo

            do k = 1, nbc
                do j = 1, ntb_c
                    do i = 1, nbc
                        tpc_wev(i,j,k) = 0.0d0
                        tnc_wev(i,j,k) = 0.0d0
                    enddo
                enddo
            enddo

            do m = 1, ntb_ark
                do l = 1, ntb_arr
                    do k = 1, ntb_art
                        do j = 1, ntb_arw
                            do i = 1, ntb_arc
                                tnccn_act(i,j,k,l,m) = 1.0
                            enddo
                        enddo
                    enddo
                enddo
            enddo

!..Check that the look-up tables are available.
            if(.not. l_mp_tables) return

!..Read a static file containing CCN activation of aerosols. The
!.. data were created from a parcel model by Feingold & Heymsfield with
!.. further changes by Eidhammer and Kriedenweis.
            
! Add this to the build table functionality for MPAS
            ! if (is_aerosol_aware) then
            !     call table_ccnAct
            !  endif

!..Collision efficiency between rain/snow and cloud water.
!     call physics_message('--- creating qc collision eff tables')
            call table_Efrw
            call table_Efsw

!..Drop evaporation.
!     call physics_message('--- creating rain evap table')
            call table_dropEvap

!..Rain collecting graupel & graupel collecting rain.
#if defined(mpas)
            call mpas_new_unit(mp_unit, unformatted = .true.)
#else
            mp_unit = 11
#endif
            open(unit=mp_unit,file='CCN_ACTIVATE.BIN',form='UNFORMATTED',status='OLD',action='READ', &
            iostat = istat)
            if(istat /= open_OK) &
            call physics_error_fatal('subroutine thompson_init: ' // &
            'failure opening CCN_ACTIVATE.BIN')
            read(mp_unit) tnccn_act
            close(unit=mp_unit)

            open(unit=mp_unit,file='MP_THOMPSON_QRacrQG_DATA.DBL',form='UNFORMATTED',status='OLD',action='READ', &
                iostat = istat)
            if(istat /= open_OK) &
                call physics_error_fatal('subroutine thompson_init: ' // &
                'failure opening MP_THOMPSON_QRacrQG.DBL')
            read(mp_unit) tcg_racg
            read(mp_unit) tmr_racg
            read(mp_unit) tcr_gacr
            read(mp_unit) tnr_racg
            read(mp_unit) tnr_gacr
            close(unit=mp_unit)
!     write(0,*) '--- end read MP_THOMPSON_QRacrQG.DBL'
!     write(0,*) 'max tcg_racg =',maxval(tcg_racg)
!     write(0,*) 'min tcg_racg =',minval(tcg_racg)
!     write(0,*) 'max tmr_racg =',maxval(tmr_racg)
!     write(0,*) 'min tmr_racg =',minval(tmr_racg)
!     write(0,*) 'max tcr_gacr =',maxval(tcr_gacr)
!     write(0,*) 'min tcr_gacr =',minval(tcr_gacr)
!     write(0,*) 'max tmg_gacr =',maxval(tmg_gacr)
!     write(0,*) 'min tmg_gacr =',minval(tmg_gacr)
!     write(0,*) 'max tnr_racg =',maxval(tnr_racg)
!     write(0,*) 'min tnr_racg =',minval(tnr_racg)
!     write(0,*) 'max tnr_gacr =',maxval(tnr_gacr)
!     write(0,*) 'min tnr_gacr =',minval(tnr_gacr)

!..Rain collecting snow & snow collecting rain.
            open(unit=mp_unit,file='MP_THOMPSON_QRacrQS_DATA.DBL',form='UNFORMATTED',status='OLD',action='READ', &
                iostat=istat)
            if(istat /= open_OK) &
                call physics_error_fatal('subroutine thompson_init: ' // &
                'failure opening MP_THOMPSON_QRacrQS.DBL')
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
!     write(0,*) '--- end read MP_THOMPSON_QRacrQS.DBL'
!     write(0,*) 'max tcs_racs1 =',maxval(tcs_racs1)
!     write(0,*) 'min tcs_racs1 =',minval(tcs_racs1)
!     write(0,*) 'max tmr_racs1 =',maxval(tmr_racs1)
!     write(0,*) 'min tmr_racs1 =',minval(tmr_racs1)
!     write(0,*) 'max tcs_racs2 =',maxval(tcs_racs2)
!     write(0,*) 'min tcs_racs2 =',minval(tcs_racs2)
!     write(0,*) 'max tmr_racs2 =',maxval(tmr_racs2)
!     write(0,*) 'min tmr_racs2 =',minval(tmr_racs2)
!     write(0,*) 'max tcr_sacr1 =',maxval(tcr_sacr1)
!     write(0,*) 'min tcr_sacr1 =',minval(tcr_sacr1)
!     write(0,*) 'max tms_sacr1 =',maxval(tms_sacr1)
!     write(0,*) 'min tms_sacr1 =',minval(tms_sacr1)
!     write(0,*) 'max tcr_sacr2 =',maxval(tcr_sacr2)
!     write(0,*) 'min tcr_sacr2 =',minval(tcr_sacr2)
!     write(0,*) 'max tms_sacr2 =',maxval(tms_sacr2)
!     write(0,*) 'min tms_sacr2 =',minval(tms_sacr2)
!     write(0,*) 'max tnr_racs1 =',maxval(tnr_racs1)
!     write(0,*) 'min tnr_racs1 =',minval(tnr_racs1)
!     write(0,*) 'max tnr_racs2 =',maxval(tnr_racs2)
!     write(0,*) 'min tnr_racs2 =',minval(tnr_racs2)
!     write(0,*) 'max tnr_sacr1 =',maxval(tnr_sacr1)
!     write(0,*) 'min tnr_sacr1 =',minval(tnr_sacr1)
!     write(0,*) 'max tnr_sacr2 =',maxval(tnr_sacr2)
!     write(0,*) 'min tnr_sacr2 =',minval(tnr_sacr2)

!..Cloud water and rain freezing (Bigg, 1953).
            open(unit=mp_unit,file='MP_THOMPSON_freezeH2O_DATA.DBL',form='UNFORMATTED',status='OLD',action='READ', &
                iostat=istat)
            if(istat /= open_OK) &
                call physics_error_fatal('subroutine thompson_init: ' // &
                'failure opening MP_THOMPSON_freezeH2O.DBL')
            read(mp_unit) tpi_qrfz
            read(mp_unit) tni_qrfz
            read(mp_unit) tpg_qrfz
            read(mp_unit) tnr_qrfz
            read(mp_unit) tpi_qcfz
            read(mp_unit) tni_qcfz
            close(unit=mp_unit)
!     write(0,*) '--- end read MP_THOMPSON_freezeH2O.DBL:'
!     write(0,*) 'max tpi_qrfz =',maxval(tpi_qrfz)
!     write(0,*) 'min tpi_qrfz =',minval(tpi_qrfz)
!     write(0,*) 'max tni_qrfz =',maxval(tni_qrfz)
!     write(0,*) 'min tni_qrfz =',minval(tni_qrfz)
!     write(0,*) 'max tpg_qrfz =',maxval(tpg_qrfz)
!     write(0,*) 'min tpg_qrfz =',minval(tpg_qrfz)
!     write(0,*) 'max tnr_qrfz =',maxval(tnr_qrfz)
!     write(0,*) 'min tnr_qrfz =',minval(tnr_qrfz)
!     write(0,*) 'max tpi_qcfz =',maxval(tpi_qcfz)
!     write(0,*) 'min tpi_qcfz =',minval(tpi_qcfz)
!     write(0,*) 'max tni_qcfz =',maxval(tni_qcfz)
!     write(0,*) 'min tni_qcfz =',minval(tni_qcfz)

!..Conversion of some ice mass into snow category.
            open(unit=mp_unit,file='MP_THOMPSON_QIautQS_DATA.DBL',form='UNFORMATTED',status='OLD',action='READ', &
                iostat=istat)
            if(istat /= open_OK) &
                call physics_error_fatal('subroutine thompson_init: ' // &
                'failure opening MP_THOMPSON_QIautQS.DBL')
            read(mp_unit) tpi_ide
            read(mp_unit) tps_iaus
            read(mp_unit) tni_iaus
            close(unit=mp_unit)
#if defined(mpas)
            call mpas_release_unit(mp_unit)
#endif
!     write(0,*) '--- end read MP_THOMPSON_QIautQS.DBL '
!     write(0,*) 'max tps_iaus =',maxval(tps_iaus)
!     write(0,*) 'min tps_iaus =',minval(tps_iaus)
!     write(0,*) 'max tni_iaus =',maxval(tni_iaus)
!     write(0,*) 'min tni_iaus =',minval(tni_iaus)

!..Initialize various constants for computing radar reflectivity.
            xam_r = am_r
            xbm_r = bm_r
            xmu_r = mu_r
            xam_s = am_s
            xbm_s = bm_s
            xmu_s = mu_s
            xam_g = am_g(idx_bg1)
            xbm_g = bm_g
            xmu_g = mu_g
            call radar_init

        endif

    end subroutine thompson_init

!=================================================================================================================
! This is a wrapper routine designed to transfer values from 3D to 1D.
! Required microphysics variables are qv, qc, qr, nr, qi, ni, qs, qg
! Optional microphysics variables are aerosol aware (nc, nwfa, nifa, nwfa2d, nifa2d), and hail aware (ng, qg)

    subroutine thompson_3d_to_1d_driver(qv, qc, qr, qi, qs, qg, qb, ni, nr, nc, ng, &
        nwfa, nifa, nwfa2d, nifa2d, th, pii, p, w, dz, dt_in, itimestep, &
        rainnc, rainncv, snownc, snowncv, graupelnc, graupelncv, sr, &
        refl_10cm, diagflag, do_radar_ref, re_cloud, re_ice, re_snow, &
        has_reqc, has_reqi, has_reqs, ntc, muc, rainprod, evapprod, &
        ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)

        implicit none

! Subroutine (3D) arguments
        integer, intent(in) :: ids,ide, jds,jde, kds,kde, ims,ime, jms,jme, kms,kme, its,ite, jts,jte, kts,kte
        real, dimension(ims:ime, kms:kme, jms:jme), intent(inout) :: qv, qc, qr, qi, qs, qg, ni, nr, th
        real, dimension(ims:ime, kms:kme, jms:jme), intent(inout) :: re_cloud, re_ice, re_snow
        integer, intent(in) :: has_reqc, has_reqi, has_reqs
        real, dimension(ims:ime, kms:kme, jms:jme), intent(in) :: pii, p, w, dz
        real, dimension(ims:ime, jms:jme), intent(inout) :: rainnc, rainncv, sr
        real, dimension(ims:ime, jms:jme), intent(in) :: ntc, muc
        real, dimension(ims:ime, kms:kme, jms:jme), intent(inout) :: rainprod, evapprod
        real, dimension(ims:ime, kms:kme, jms:jme), optional, intent(inout) :: nc, nwfa, nifa, qb, ng
        real, dimension(ims:ime, jms:jme), optional, intent(in) :: nwfa2d, nifa2d
        real, dimension(ims:ime, kms:kme, jms:jme), intent(inout), optional :: refl_10cm
        real, dimension(ims:ime, jms:jme), optional, intent(inout) :: snownc, snowncv, graupelnc, graupelncv
        real, intent(in) :: dt_in
        integer, intent(in) :: itimestep

! Local (1d) variables
        real, dimension(kts:kte) :: qv1d, qc1d, qi1d, qr1d, qs1d, qg1d, qb1d, ni1d, nr1d, nc1d, ng1d, &
            nwfa1d, nifa1d, t1d, p1d, w1d, dz1d, rho, dbz
        real, dimension(kts:kte) :: re_qc1d, re_qi1d, re_qs1d
        real, dimension(kts:kte):: rainprod1d, evapprod1d
        real, dimension(its:ite, jts:jte) :: pcp_ra, pcp_sn, pcp_gr, pcp_ic
        real :: dt, pptrain, pptsnow, pptgraul, pptice
        real :: qc_max, qr_max, qs_max, qi_max, qg_max, ni_max, nr_max
        real :: nwfa1
        real :: ygra1, zans1
        real :: graupel_vol
        double precision :: lamg, lam_exp, lamr, n0_min, n0_exp
        integer :: i, j, k
        integer :: imax_qc, imax_qr, imax_qi, imax_qs, imax_qg, imax_ni, imax_nr
        integer :: jmax_qc, jmax_qr, jmax_qi, jmax_qs, jmax_qg, jmax_ni, jmax_nr
        integer :: kmax_qc, kmax_qr, kmax_qi, kmax_qs, kmax_qg, kmax_ni, kmax_nr
        integer :: i_start, j_start, i_end, j_end
        logical, optional, intent(in) :: diagflag
        integer, optional, intent(in) :: do_radar_ref
        character*256:: mp_debug

!=================================================================================================================        
        i_start = its
        j_start = jts
        i_end = min(ite, ide-1)
        j_end = min(jte, jde-1)
        dt = dt_in

        qc_max = 0.0
        qr_max = 0.0
        qs_max = 0.0
        qi_max = 0.0
        qg_max = 0.0
        ni_max = 0.0
        nr_max = 0.0
        imax_qc = 0
        imax_qr = 0
        imax_qi = 0
        imax_qs = 0
        imax_qg = 0
        imax_ni = 0
        imax_nr = 0
        jmax_qc = 0
        jmax_qr = 0
        jmax_qi = 0
        jmax_qs = 0
        jmax_qg = 0
        jmax_ni = 0
        jmax_nr = 0
        kmax_qc = 0
        kmax_qr = 0
        kmax_qi = 0
        kmax_qs = 0
        kmax_qg = 0
        kmax_ni = 0
        kmax_nr = 0
        
        do i = 1, 256
            mp_debug(i:i) = char(0)
        enddo

!     if (.NOT. is_aerosol_aware .AND. PRESENT(nc) .AND. PRESENT(nwfa)  &
!               .AND. PRESENT(nifa) .AND. PRESENT(nwfa2d)) then
!        write(mp_debug,*) 'WARNING, nc-nwfa-nifa-nwfa2d present but is_aerosol_aware is FALSE'
!        CALL wrf_debug(0, mp_debug)
!     endif

!=================================================================================================================
        j_loop:  do j = j_start, j_end
            i_loop:  do i = i_start, i_end

                pptrain = 0.0
                pptsnow = 0.0
                pptgraul = 0.0
                pptice = 0.0
                rainncv(i,j) = 0.0
                if (present(snowncv)) then
                    snowncv(i,j) = 0.0
                endif
                if (present(graupelncv)) then
                    graupelncv(i,j) = 0.0
                endif
                sr(i,j) = 0.0

                ! ntc and muc are defined in mpas submodule
                ! ntc is different for land versus ocean
                nt_c = ntc(i,j)
                mu_c = muc(i,j)
                
                do k = kts, kte
                    t1d(k) = th(i,k,j) * pii(i,k,j)
                    p1d(k) = p(i,k,j)
                    w1d(k) = w(i,k,j)
                    dz1d(k) = dz(i,k,j)
                    qv1d(k) = qv(i,k,j)
                    qc1d(k) = qc(i,k,j)
                    qi1d(k) = qi(i,k,j)
                    qr1d(k) = qr(i,k,j)
                    qs1d(k) = qs(i,k,j)
                    qg1d(k) = qg(i,k,j)
                    ni1d(k) = ni(i,k,j)
                    nr1d(k) = nr(i,k,j)
                    rho(k) = 0.622 * p1d(k) / (R * t1d(k) * (qv1d(k)+0.622))

! nwfa, nifa, and nc are optional aerosol-aware variables
                    if (present(nwfa)) then
                       nwfa1d(k) = nwfa(i,k,j)
                    else
                       nwfa1d(k) = 11.1e6 / rho(k)
                       configs%aerosol_aware = .false.
                    endif

                    if (present(nifa)) then
                       nifa1d(k) = nifa(i,k,j)
                    else
                       nifa1d(k) = nain1*0.01 / rho(k)
                       configs%aerosol_aware = .false.
                    endif

                    if (present(nc)) then
                       nc1d(k) = nc(i,k,j)
                    else
                       nc1d(k) = nt_c / rho(k)
                       configs%aerosol_aware = .false.
                    endif
                 enddo

! ng and qb are optional hail-aware variables
                if ((present(ng)) .and. (present(qb))) then
                    configs%hail_aware = .true.
                    do k = kts, kte
                        ng1d(k) = ng(i,k,j)
                        qb1d(k) = qb(i,k,j)
                    enddo
                else
                    do k = kte, kts, -1
! This is the one-moment graupel formulation
                        if (qg1d(k) > R1) then
                            ygra1 = alog10(max(1.e-9, qg1d(k)*rho(k)))
                            zans1 = 3.0 + 2.0/7.0*(ygra1+8.0)
                            zans1 = max(2.0, min(zans1, 6.0))
                            n0_exp = 10.0**(zans1)
                            lam_exp = (n0_exp*am_g(idx_bg1)*cgg(1,1) / (rho(k)*qg1d(k)))**oge1
                            lamg = lam_exp * (cgg(3,1)*ogg2*ogg1)**obmg
                            ng1d(k) = cgg(2,1) * ogg3*rho(k) * qg1d(k) * lamg**bm_g / am_g(idx_bg1)
                            ng1d(k) = max(R2, (ng1d(k)/rho(k)))
                            qb1d(k) = qg1d(k) / rho_g(idx_bg1)
                        else
                            ng1d(k) = 0
                            qb1d(k) = 0
                        endif
                    enddo
                endif

!=================================================================================================================
! Main call to the 1D microphysics
                call mp_thompson_main(qv1d, qc1d, qi1d, qr1d, qs1d, qg1d, qb1d, ni1d, nr1d, nc1d, ng1d, &
                    nwfa1d, nifa1d, t1d, p1d, w1d, dz1d, pptrain, pptsnow, pptgraul, pptice, &
                    rainprod1d, evapprod1d, kts, kte, dt, i, j, configs)
                    
!=================================================================================================================
! Compute diagnostics and return output to 3D

                pcp_ra(i,j) = pptrain
                pcp_sn(i,j) = pptsnow
                pcp_gr(i,j) = pptgraul
                pcp_ic(i,j) = pptice
                rainncv(i,j) = pptrain + pptsnow + pptgraul + pptice
                rainnc(i,j) = rainnc(i,j) + pptrain + pptsnow + pptgraul + pptice
                if (present(snowncv) .and. present(snownc)) then
                    snowncv(i,j) = pptsnow + pptice
                    snownc(i,j) = snownc(i,j) + pptsnow + pptice
                endif
                if (present(graupelncv) .and. present(graupelnc)) then
                    graupelncv(i,j) = pptgraul
                    graupelnc(i,j) = graupelnc(i,j) + pptgraul
                endif
                sr(i,j) = (pptsnow + pptgraul + pptice) / (rainncv(i,j)+R1)

                if ((present(ng)) .and. (present(qb))) then
                    do k = kts, kte
                        ng(i,k,j) = ng1d(k)
                        qb(i,k,j) = qb1d(k)
                    enddo
                endif

                do k = kts, kte
                    if (present(nc)) nc(i,k,j) = nc1d(k)
                    if (present(nwfa)) nwfa(i,k,j) = nwfa1d(k)
                    if (present(nifa)) nifa(i,k,j) = nifa1d(k)
                    qv(i,k,j) = qv1d(k)
                    qc(i,k,j) = qc1d(k)
                    qi(i,k,j) = qi1d(k)
                    qr(i,k,j) = qr1d(k)
                    qs(i,k,j) = qs1d(k)
                    qg(i,k,j) = qg1d(k)
                    ni(i,k,j) = ni1d(k)
                    nr(i,k,j) = nr1d(k)
                    th(i,k,j) = t1d(k) / pii(i,k,j)
                    rainprod(i,k,j) = rainprod1d(k)
                    evapprod(i,k,j) = evapprod1d(k)

                    if (qc1d(k) .gt. qc_max) then
                        imax_qc = i
                        jmax_qc = j
                        kmax_qc = k
                        qc_max = qc1d(k)
                    elseif (qc1d(k) .lt. 0.0) then
                        write(mp_debug,*) 'WARNING, negative qc ', qc1d(k),        &
                            ' at i,j,k=', i,j,k
!            CALL wrf_debug(150, mp_debug)
                    endif
                    if (qr1d(k) .gt. qr_max) then
                        imax_qr = i
                        jmax_qr = j
                        kmax_qr = k
                        qr_max = qr1d(k)
                    elseif (qr1d(k) .lt. 0.0) then
                        write(mp_debug,*) 'WARNING, negative qr ', qr1d(k),        &
                            ' at i,j,k=', i,j,k
!            CALL wrf_debug(150, mp_debug)
                    endif
                    if (nr1d(k) .gt. nr_max) then
                        imax_nr = i
                        jmax_nr = j
                        kmax_nr = k
                        nr_max = nr1d(k)
                    elseif (nr1d(k) .lt. 0.0) then
                        write(mp_debug,*) 'WARNING, negative nr ', nr1d(k),        &
                            ' at i,j,k=', i,j,k
!            CALL wrf_debug(150, mp_debug)
                    endif
                    if (qs1d(k) .gt. qs_max) then
                        imax_qs = i
                        jmax_qs = j
                        kmax_qs = k
                        qs_max = qs1d(k)
                    elseif (qs1d(k) .lt. 0.0) then
                        write(mp_debug,*) 'WARNING, negative qs ', qs1d(k),        &
                            ' at i,j,k=', i,j,k
!            CALL wrf_debug(150, mp_debug)
                    endif
                    if (qi1d(k) .gt. qi_max) then
                        imax_qi = i
                        jmax_qi = j
                        kmax_qi = k
                        qi_max = qi1d(k)
                    elseif (qi1d(k) .lt. 0.0) then
                        write(mp_debug,*) 'WARNING, negative qi ', qi1d(k),        &
                            ' at i,j,k=', i,j,k
!            CALL wrf_debug(150, mp_debug)
                    endif
                    if (qg1d(k) .gt. qg_max) then
                        imax_qg = i
                        jmax_qg = j
                        kmax_qg = k
                        qg_max = qg1d(k)
                    elseif (qg1d(k) .lt. 0.0) then
                        write(mp_debug,*) 'WARNING, negative qg ', qg1d(k),        &
                            ' at i,j,k=', i,j,k
!            CALL wrf_debug(150, mp_debug)
                    endif
                    if (ni1d(k) .gt. ni_max) then
                        imax_ni = i
                        jmax_ni = j
                        kmax_ni = k
                        ni_max = ni1d(k)
                    elseif (ni1d(k) .lt. 0.0) then
                        write(mp_debug,*) 'WARNING, negative ni ', ni1d(k),        &
                            ' at i,j,k=', i,j,k
!            CALL wrf_debug(150, mp_debug)
                    endif
                    if (qv1d(k) .lt. 0.0) then
                        write(mp_debug,*) 'WARNING, negative qv ', qv1d(k),        &
                            ' at i,j,k=', i,j,k
!            CALL wrf_debug(150, mp_debug)
                        if (k.lt.kte-2 .and. k.gt.kts+1) then
                            write(mp_debug,*) '   below and above are: ', qv(i,k-1,j), qv(i,k+1,j)
!               CALL wrf_debug(150, mp_debug)
                            qv(i,k,j) = MAX(1.E-7, 0.5*(qv(i,k-1,j) + qv(i,k+1,j)))
                        else
                            qv(i,k,j) = 1.E-7
                        endif
                    endif
                enddo

!=================================================================================================================
! Reflectivity
                call calc_refl10cm (qv1d, qc1d, qr1d, nr1d, qs1d, qg1d, ng1d, qb1d, &
                     t1d, p1d, dBZ, kts, kte, i, j, configs)
                do k = kts, kte
                   refl_10cm(i,k,j) = max(-35., dBZ(k))
                enddo

! Cloud, ice, and snow effective radius
                if (has_reqc.ne.0 .and. has_reqi.ne.0 .and. has_reqs.ne.0) then
                    do k = kts, kte
                        re_qc1d(k) = 2.49E-6
                        re_qi1d(k) = 4.99E-6
                        re_qs1d(k) = 9.99E-6
                    enddo
                    call calc_effectRad (t1d, p1d, qv1d, qc1d, nc1d, qi1d, ni1d, qs1d,  &
                        re_qc1d, re_qi1d, re_qs1d, kts, kte, configs)
                    do k = kts, kte
                        re_cloud(i,k,j) = MAX(2.49E-6, MIN(re_qc1d(k), 50.E-6))
                        re_ice(i,k,j)   = MAX(4.99E-6, MIN(re_qi1d(k), 125.E-6))
                        re_snow(i,k,j)  = MAX(9.99E-6, MIN(re_qs1d(k), 999.E-6))
                    enddo
                endif

            enddo i_loop
        enddo j_loop

! DEBUG - GT
!        write(mp_debug,'(a,7(a,e13.6,1x,a,i3,a,i3,a,i3,a,1x))') 'MP-GT:', &
!            'qc: ', qc_max, '(', imax_qc, ',', jmax_qc, ',', kmax_qc, ')', &
!            'qr: ', qr_max, '(', imax_qr, ',', jmax_qr, ',', kmax_qr, ')', &
!            'qi: ', qi_max, '(', imax_qi, ',', jmax_qi, ',', kmax_qi, ')', &
!            'qs: ', qs_max, '(', imax_qs, ',', jmax_qs, ',', kmax_qs, ')', &
!            'qg: ', qg_max, '(', imax_qg, ',', jmax_qg, ',', kmax_qg, ')', &
!            'ni: ', ni_max, '(', imax_ni, ',', jmax_ni, ',', kmax_ni, ')', &
!            'nr: ', nr_max, '(', imax_nr, ',', jmax_nr, ',', kmax_nr, ')'
!     CALL wrf_debug(150, mp_debug)
! END DEBUG - GT

        do i = 1, 256
            mp_debug(i:i) = char(0)
        enddo

    end subroutine thompson_3d_to_1d_driver
!=================================================================================================================

end module module_mp_thompson
