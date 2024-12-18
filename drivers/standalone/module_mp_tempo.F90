! 3D TEMPO Driver for Standalone Testing
!=================================================================================================================
module module_mp_tempo

    use machine, only: wp => kind_phys, sp => kind_sngl_prec, dp => kind_dbl_prec
    use module_mp_tempo_params
    use module_mp_tempo_utils, only : create_bins, table_Efrw, table_Efsw, table_dropEvap, &
         calc_refl10cm, calc_effectRad
    use module_mp_tempo_main, only : mp_tempo_main
    ! use mp_radar

    implicit none

contains
    !=================================================================================================================
    ! This subroutine handles initialzation of the microphysics scheme including building of lookup tables,
    ! allocating arrays for the microphysics scheme, and defining gamma function variables.

    ! Input:
    !   l_mp_tables = .false. to build lookup tables. If l_mp_tables = .true., lookup tables are not built.

    ! AAJ No support yet for hail_aware in microphysics driver
    subroutine tempo_init(l_mp_tables, hail_aware_flag, aerosol_aware_flag)

        ! Input arguments:
        logical, intent(in) :: l_mp_tables, hail_aware_flag
        logical, intent(in), optional :: aerosol_aware_flag

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
            write(*,*) '--- tempo_init() called with hail_aware_flag = ' // trim(message)

            if (present(aerosol_aware_flag)) then
                configs%aerosol_aware = aerosol_aware_flag
                write(message, '(L1)') configs%aerosol_aware
                write(*,*) '--- tempo_init() called with aerosol_aware_flag = ' // trim(message)
            endif
        endif

        if (configs%hail_aware) then
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

            inquire(file='MP_TEMPO_QRacrQG_DATA.DBL', exist=qr_acr_qg_exists)
            if (qr_acr_qg_exists) then ! Check again that file exists

                ! Add check on qr_ac_qg filesize to determine if table includes hail-awareness (dimNRHG=9)
                qr_acr_qg_check = dp * num_records * (dimNRHG * ntb_g1 * ntb_g * ntb_r1 * ntb_r + 1)
                qr_acr_qg_dim1size = dp * num_records * (NRHG1 * ntb_g1 * ntb_g * ntb_r1 * ntb_r + 1)
                qr_acr_qg_dim9size = dp * num_records * (NRHG * ntb_g1 * ntb_g * ntb_r1 * ntb_r + 1)
                
                inquire(file='MP_TEMPO_QRacrQG_DATA.DBL', size=qr_acr_qg_filesize)

                if (qr_acr_qg_filesize == qr_acr_qg_dim1size) then
                    using_hail_aware_table = .false.
                    write(*,*) '--- tempo_init() ' // &
                        'Lookup table for qr_acr_qg is not hail aware.'
                    dimNRHG = NRHG1
                    if (hail_aware_flag) then
                        write(*,*) '--- tempo_init() Cannot use hail-aware microphysics ' // &
                            'with non hail-aware qr_acr_qg lookup table. ' // &
                            'Please rebuild table with parameter build_hail_aware_table set to true.'
                    endif
                elseif (qr_acr_qg_filesize == qr_acr_qg_dim9size) then
                    using_hail_aware_table = .true.
                    write(*,*) '--- tempo_init() ' // &
                        'Lookup table for qr_acr_qg is hail aware.'
                    dimNRHG = NRHG
                else
                    using_hail_aware_table = .false.
                    if (hail_aware_flag) using_hail_aware_table = .true.
                    write(*,*) '--- tempo_init() ' // &
                        'Could not determine if lookup table for qr_acr_qg is hail aware based on file size.'
                endif
            endif
        endif

        !=================================================================================================================
        ! Allocate space for lookup tables (J. Michalakes 2009Jun08).
        if (.not. allocated(tcg_racg)) then
            allocate(tcg_racg(ntb_g1,ntb_g,dimNRHG,ntb_r1,ntb_r))
            micro_init = .true.
        endif

        ! Rain-graupel (including array above tcg_racg)
        if (.not. allocated(tmr_racg)) allocate(tmr_racg(ntb_g1,ntb_g,dimNRHG,ntb_r1,ntb_r))
        if (.not. allocated(tcr_gacr)) allocate(tcr_gacr(ntb_g1,ntb_g,dimNRHG,ntb_r1,ntb_r))
        if (.not. allocated(tnr_racg)) allocate(tnr_racg(ntb_g1,ntb_g,dimNRHG,ntb_r1,ntb_r))
        if (.not. allocated(tnr_gacr)) allocate(tnr_gacr(ntb_g1,ntb_g,dimNRHG,ntb_r1,ntb_r))

        ! Rain-snow
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

        ! Cloud water freezing
        if (.not. allocated(tpi_qcfz)) allocate(tpi_qcfz(ntb_c,nbc,ntb_t1,ntb_IN))
        if (.not. allocated(tni_qcfz)) allocate(tni_qcfz(ntb_c,nbc,ntb_t1,ntb_IN))

        ! Rain freezing
        if (.not. allocated(tpi_qrfz)) allocate(tpi_qrfz(ntb_r,ntb_r1,ntb_t1,ntb_IN))
        if (.not. allocated(tpg_qrfz)) allocate(tpg_qrfz(ntb_r,ntb_r1,ntb_t1,ntb_IN))
        if (.not. allocated(tni_qrfz)) allocate(tni_qrfz(ntb_r,ntb_r1,ntb_t1,ntb_IN))
        if (.not. allocated(tnr_qrfz)) allocate(tnr_qrfz(ntb_r,ntb_r1,ntb_t1,ntb_IN))

        ! Ice growth and conversion to snow
        if (.not. allocated(tps_iaus)) allocate(tps_iaus(ntb_i,ntb_i1))
        if (.not. allocated(tni_iaus)) allocate(tni_iaus(ntb_i,ntb_i1))
        if (.not. allocated(tpi_ide)) allocate(tpi_ide(ntb_i,ntb_i1))

        ! Collision efficiencies
        if (.not. allocated(t_efrw)) allocate(t_efrw(nbr,nbc))
        if (.not. allocated(t_efsw)) allocate(t_efsw(nbs,nbc))

        ! Cloud water
        if (.not. allocated(tnr_rev)) allocate(tnr_rev(nbr,ntb_r1,ntb_r))
        if (.not. allocated(tpc_wev)) allocate(tpc_wev(nbc,ntb_c,nbc))
        if (.not. allocated(tnc_wev)) allocate(tnc_wev(nbc,ntb_c,nbc))

        ! CCN
        if (.not. allocated(tnccn_act)) allocate(tnccn_act(ntb_arc,ntb_arw,ntb_art,ntb_arr,ntb_ark))
        !=================================================================================================================
        if (micro_init) then

            ! Schmidt number to one-third used numerous times.
            Sc3 = Sc**(1./3.)

            ! Compute minimum ice diameter from mass and minimum snow/graupel mass from diameter
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
                ccg(1,n) = gamma(cce(1,n))
                ccg(2,n) = gamma(cce(2,n))
                ccg(3,n) = gamma(cce(3,n))
                ccg(4,n) = gamma(cce(4,n))
                ccg(5,n) = gamma(cce(5,n))
                ocg1(n) = 1.0 / ccg(1,n)
                ocg2(n) = 1.0 / ccg(2,n)
            enddo

            cie(1) = mu_i + 1.
            cie(2) = bm_i + mu_i + 1.
            cie(3) = bm_i + mu_i + bv_i + 1.
            cie(4) = mu_i + bv_i + 1.
            cie(5) = mu_i + 2.
            cie(6) = bm_i*0.5 + mu_i + bv_i + 1.
            cie(7) = bm_i*0.5 + mu_i + 1.
            cig(1) = gamma(cie(1))
            cig(2) = gamma(cie(2))
            cig(3) = gamma(cie(3))
            cig(4) = gamma(cie(4))
            cig(5) = gamma(cie(5))
            cig(6) = gamma(cie(6))
            cig(7) = gamma(cie(7))
            oig1 = 1.0 / cig(1)
            oig2 = 1.0 / cig(2)
            obmi = 1.0 / bm_i

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
                crg(n) = gamma(cre(n))
            enddo

            obmr = 1.0 / bm_r
            ore1 = 1.0 / cre(1)
            org1 = 1.0 / crg(1)
            org2 = 1.0 / crg(2)
            org3 = 1.0 / crg(3)

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
                csg(n) = gamma(cse(n))
            enddo

            oams = 1.0 / am_s
            obms = 1.0 / bm_s
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
                    cgg(n,m) = gamma(cge(n,m))
                enddo
            enddo

            oamg = 1.0 / am_g
            obmg = 1.0 / bm_g

            do m = 1, NRHG
                oamg(m) = 1.0 / am_g(m)
                ocmg(m) = oamg(m)**obmg
            enddo

            oge1 = 1.0 / cge(1,1)
            ogg1 = 1.0 / cgg(1,1)
            ogg2 = 1.0 / cgg(2,1)
            ogg3 = 1.0 / cgg(3,1)

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
            t2_qr_ev = 0.308 * Sc3 * SQRT(av_r) * crg(11)

            ! Sublimation/depositional growth of snow
            t1_qs_sd = 0.86
            t2_qs_sd = 0.28 * Sc3 * SQRT(av_s)

            ! Melting of snow
            t1_qs_me = PI * 4. *C_sqrd * olfus * 0.86
            t2_qs_me = PI * 4. *C_sqrd * olfus * 0.28 * Sc3 * SQRT(av_s)

            ! Sublimation/depositional growth of graupel
            t1_qg_sd = 0.86 * cgg(10,1)
            !     t2_qg_sd = 0.28*Sc3*SQRT(av_g) * cgg(11)

            ! Melting of graupel
            t1_qg_me = PI * 4. * C_cube * olfus * 0.86 * cgg(10,1)
            !     t2_qg_me = PI*4.*C_cube*olfus * 0.28*Sc3*SQRT(av_g) * cgg(11)


            ! Constants for helping find lookup table indexes.
            nic2 = nint(log10(r_c(1)))
            nii2 = nint(log10(r_i(1)))
            nii3 = nint(log10(Nt_i(1)))
            nir2 = nint(log10(r_r(1)))
            nir3 = nint(log10(N0r_exp(1)))
            nis2 = nint(log10(r_s(1)))
            nig2 = nint(log10(r_g(1)))
            nig3 = nint(log10(N0g_exp(1)))
            niIN2 = nint(log10(Nt_IN(1)))

            ! Create bins of cloud water (from minimum diameter to 100 microns).
            Dc(1) = D0c*1.0_dp
            dtc(1) = D0c*1.0_dp
            do n = 2, nbc
                Dc(n) = Dc(n-1) + 1.0e-6_dp
                dtc(n) = (Dc(n) - Dc(n-1))
            enddo

            ! Create bins of cloud ice (from min diameter up to 2x min snow size).
            call create_bins(numbins=nbi, lowbin=D0i*1.0_dp, highbin=D0s*2.0_dp, &
                bins=Di, deltabins=dti)

            ! Create bins of rain (from min diameter up to 5 mm).
            call create_bins(numbins=nbr, lowbin=D0r*1.0_dp, highbin=0.005_dp, &
                bins=Dr, deltabins=dtr)

            ! Create bins of snow (from min diameter up to 2 cm).
            call create_bins(numbins=nbs, lowbin=D0s*1.0_dp, highbin=0.02_dp, &
                bins=Ds, deltabins=dts)

            ! Create bins of graupel (from min diameter up to 5 cm).
            call create_bins(numbins=nbg, lowbin=D0g*1.0_dp, highbin=0.05_dp, &
                bins=Dg, deltabins=dtg)

            ! Create bins of cloud droplet number concentration (1 to 3000 per cc).
            call create_bins(numbins=nbc, lowbin=1.0_dp, highbin=3000.0_dp, &
                bins=t_Nc)
            t_Nc = t_Nc * 1.0e6_dp
            nic1 = log(t_Nc(nbc)/t_Nc(1))

            !=================================================================================================================
            ! Create lookup tables for most costly calculations.

            do m = 1, ntb_r
                do k = 1, ntb_r1
                    do n = 1, dimNRHG
                        do j = 1, ntb_g
                            do i = 1, ntb_g1
                                tcg_racg(i,j,n,k,m) = 0.0_dp
                                tmr_racg(i,j,n,k,m) = 0.0_dp
                                tcr_gacr(i,j,n,k,m) = 0.0_dp
                                tnr_racg(i,j,n,k,m) = 0.0_dp
                                tnr_gacr(i,j,n,k,m) = 0.0_dp
                            enddo
                        enddo
                    enddo
                enddo
            enddo

            do m = 1, ntb_r
                do k = 1, ntb_r1
                    do j = 1, ntb_t
                        do i = 1, ntb_s
                            tcs_racs1(i,j,k,m) = 0.0_dp
                            tmr_racs1(i,j,k,m) = 0.0_dp
                            tcs_racs2(i,j,k,m) = 0.0_dp
                            tmr_racs2(i,j,k,m) = 0.0_dp
                            tcr_sacr1(i,j,k,m) = 0.0_dp
                            tms_sacr1(i,j,k,m) = 0.0_dp
                            tcr_sacr2(i,j,k,m) = 0.0_dp
                            tms_sacr2(i,j,k,m) = 0.0_dp
                            tnr_racs1(i,j,k,m) = 0.0_dp
                            tnr_racs2(i,j,k,m) = 0.0_dp
                            tnr_sacr1(i,j,k,m) = 0.0_dp
                            tnr_sacr2(i,j,k,m) = 0.0_dp
                        enddo
                    enddo
                enddo
            enddo

            do m = 1, ntb_IN
                do k = 1, ntb_t1
                    do j = 1, ntb_r1
                        do i = 1, ntb_r
                            tpi_qrfz(i,j,k,m) = 0.0_dp
                            tni_qrfz(i,j,k,m) = 0.0_dp
                            tpg_qrfz(i,j,k,m) = 0.0_dp
                            tnr_qrfz(i,j,k,m) = 0.0_dp
                        enddo
                    enddo
                    do j = 1, nbc
                        do i = 1, ntb_c
                            tpi_qcfz(i,j,k,m) = 0.0_dp
                            tni_qcfz(i,j,k,m) = 0.0_dp
                        enddo
                    enddo
                enddo
            enddo

            do j = 1, ntb_i1
                do i = 1, ntb_i
                    tps_iaus(i,j) = 0.0_dp
                    tni_iaus(i,j) = 0.0_dp
                    tpi_ide(i,j) = 0.0_dp
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
                        tnr_rev(i,j,k) = 0.0_dp
                    enddo
                enddo
            enddo

            do k = 1, nbc
                do j = 1, ntb_c
                    do i = 1, nbc
                        tpc_wev(i,j,k) = 0.0_dp
                        tnc_wev(i,j,k) = 0.0_dp
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

            !=================================================================================================================
            ! Check that the look-up tables are available.
            if (.not. l_mp_tables) return

            ! Collision efficiency between rain/snow and cloud water.
            call table_Efrw ! => fills t_Efrw
            call table_Efsw ! => fills t_Efsw

            ! Drop evaporation
            call table_dropEvap ! => fills tpc_wev and tnc_wev

            ! Read a static file containing CCN activation of aerosols. The data were created from a parcel model
            ! by Feingold & Heymsfield with further changes by Eidhammer and Kriedenweis.
            mp_unit = 11
            
            open(unit=mp_unit,file='CCN_ACTIVATE.BIN',form='unformatted',status='old',action='read',iostat=istat, &
                 convert='big_endian')
            if (istat /= open_OK) then
                write(*,*) '--- tempo_init() failure opening CCN_ACTIVATE.BIN'
            endif
            read(mp_unit) tnccn_act
            close(unit=mp_unit)

            ! Rain collecting graupel & graupel collecting rain.

            open(unit=mp_unit,file='MP_TEMPO_QRacrQG_DATA.DBL',form='unformatted',status='old',action='read', &
                iostat=istat,                  convert='big_endian')
            if (istat /= open_OK) then
                write(*,*) '--- tempo_init() failure opening MP_TEMPO_QRacrQG.DBL'
            endif
            read(mp_unit) tcg_racg
            read(mp_unit) tmr_racg
            read(mp_unit) tcr_gacr
            read(mp_unit) tnr_racg
            read(mp_unit) tnr_gacr
            close(unit=mp_unit)

            ! Rain collecting snow & snow collecting rain.
            open(unit=mp_unit,file='MP_TEMPO_QRacrQS_DATA.DBL',form='unformatted',status='old',action='read', &
                iostat=istat,                  convert='big_endian')
            if (istat /= open_OK) then
                write(*,*) '--- tempo_init() failure opening MP_TEMPO_QRacrQS.DBL'
            endif
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

            ! Cloud water and rain freezing (Bigg, 1953).
            open(unit=mp_unit,file='MP_TEMPO_freezeH2O_DATA.DBL',form='unformatted',status='old',action='read', &
                iostat=istat,                  convert='big_endian')
            if (istat /= open_OK) then
                write(*,*) '--- tempo_init() failure opening MP_TEMPO_freezeH2O.DBL'
            endif
            read(mp_unit) tpi_qrfz
            read(mp_unit) tni_qrfz
            read(mp_unit) tpg_qrfz
            read(mp_unit) tnr_qrfz
            read(mp_unit) tpi_qcfz
            read(mp_unit) tni_qcfz
            close(unit=mp_unit)

            ! Conversion of some ice mass into snow category.
            open(unit=mp_unit,file='MP_TEMPO_QIautQS_DATA.DBL',form='unformatted',status='old',action='read', &
                iostat=istat,                  convert='big_endian')
            if (istat /= open_OK) then
                write(*,*) '--- tempo_init() failure opening MP_TEMPO_QIautQS.DBL'
            endif
            read(mp_unit) tpi_ide
            read(mp_unit) tps_iaus
            read(mp_unit) tni_iaus
            close(unit=mp_unit)

            ! Initialize various constants for computing radar reflectivity.
            ! xam_r = am_r
            ! xbm_r = bm_r
            ! xmu_r = mu_r
            ! xam_s = am_s
            ! xbm_s = bm_s
            ! xmu_s = mu_s
            ! xam_g = am_g(idx_bg1)
            ! xbm_g = bm_g
            ! xmu_g = mu_g
            ! call radar_init

        endif ! micro_init

    end subroutine tempo_init

    !=================================================================================================================
    ! This is a wrapper routine designed to transfer values from 3D to 1D.
    ! Required microphysics variables are qv, qc, qr, nr, qi, ni, qs, qg
    ! Optional microphysics variables are aerosol aware (nc, nwfa, nifa, nwfa2d, nifa2d), and hail aware (ng, qg)

    subroutine tempo_3d_to_1d_driver(qv, qc, qr, qi, qs, qg, qb, ni, nr, nc, ng, &
        nwfa, nifa, nwfa2d, nifa2d, th, pii, p, w, dz, dt_in, itimestep, &
        rainnc, rainncv, snownc, snowncv, graupelnc, graupelncv, sr, &
        refl_10cm, diagflag, do_radar_ref, re_cloud, re_ice, re_snow, &
        has_reqc, has_reqi, has_reqs, ntc, muc, &
        ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)

        ! Subroutine (3D) arguments
        integer, intent(in) :: ids,ide, jds,jde, kds,kde, ims,ime, jms,jme, kms,kme, its,ite, jts,jte, kts,kte
        real, dimension(ims:ime, kms:kme, jms:jme), intent(inout) :: qv, qc, qr, qi, qs, qg, ni, nr, th
        real, dimension(ims:ime, kms:kme, jms:jme), intent(inout) :: re_cloud, re_ice, re_snow
        integer, intent(in) :: has_reqc, has_reqi, has_reqs
        real, dimension(ims:ime, kms:kme, jms:jme), intent(in) :: pii, p, w, dz
        real, dimension(ims:ime, jms:jme), intent(inout) :: rainnc, rainncv, sr
        real, dimension(ims:ime, jms:jme), intent(in), optional :: ntc, muc
        real, dimension(ims:ime, kms:kme, jms:jme), intent(inout), optional :: nc, nwfa, nifa, qb, ng
        real, dimension(ims:ime, jms:jme), intent(in), optional :: nwfa2d, nifa2d
        real, dimension(ims:ime, kms:kme, jms:jme), intent(inout), optional :: refl_10cm
        real, dimension(ims:ime, jms:jme), intent(inout), optional :: snownc, snowncv, graupelnc, graupelncv
        real, intent(in) :: dt_in
        integer, intent(in) :: itimestep

        ! Local (1d) variables
        real, dimension(kts:kte) :: qv1d, qc1d, qi1d, qr1d, qs1d, qg1d, qb1d, ni1d, nr1d, nc1d, ng1d, &
            nwfa1d, nifa1d, t1d, p1d, w1d, dz1d, rho, dbz
        real, dimension(kts:kte) :: re_qc1d, re_qi1d, re_qs1d
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
        character(len=132) :: message

        !=================================================================================================================
        i_start = its
        j_start = jts
        i_end = max(ite, ide-1)
        j_end = max(jte, jde-1)
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

                ! ntc and muc are defined in mpas submodule based on landmask
                if (present(ntc)) then
                    Nt_c = ntc(i,j)
                    mu_c = muc(i,j)
                else
                    Nt_c = Nt_c_o
                    mu_c = 4
                endif

                !=================================================================================================================
                ! Begin k loop
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
                    rho(k) = RoverRv * p1d(k) / (R * t1d(k) * (qv1d(k)+RoverRv))

                    ! nwfa, nifa, and nc are optional aerosol-aware variables
                    if (present(nwfa)) then
                        nwfa1d(k) = nwfa(i,k,j)
                    else
                        nwfa1d(k) = nwfa_default / rho(k)
                        configs%aerosol_aware = .false.
                    endif

                    if (present(nifa)) then
                        nifa1d(k) = nifa(i,k,j)
                    else
                        nifa1d(k) = nifa_default / rho(k)
                        configs%aerosol_aware = .false.
                    endif

                    if (present(nc)) then
                        nc1d(k) = nc(i,k,j)
                    else
                        nc1d(k) = Nt_c / rho(k)
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
                            ygra1 = log10(max(1.e-9, qg1d(k)*rho(k)))
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
                
                if (itimestep == 1) then
                   write(*,*) '--- tempo_3d_to_1d_driver() configuration...'
                   write(message, '(L1)') configs%hail_aware
                   write(*,*) '       hail_aware_flag = ' // trim(message)
                   write(message, '(L1)') configs%aerosol_aware
                   write(*,*) '       aerosol_aware_flag = ' // trim(message)
                   write(*,*) 'calling mp_tempo_main() at itimestep = 1'
                endif

                !=================================================================================================================
                ! Main call to the 1D microphysics
                call mp_tempo_main(qv1d=qv1d, qc1d=qc1d, qi1d=qi1d, qr1d=qr1d, qs1d=qs1d, qg1d=qg1d, qb1d=qb1d, &
                           ni1d=ni1d, nr1d=nr1d, nc1d=nc1d, ng1d=ng1d, nwfa1d=nwfa1d, nifa1d=nifa1d, t1d=t1d, p1d=p1d, &
                           w1d=w1d, dzq=dz1d, pptrain=pptrain, pptsnow=pptsnow, pptgraul=pptgraul, pptice=pptice, &
                           kts=kts, kte=kte, dt=dt, ii=i, jj=j, configs=configs)

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
                sr(i,j) = (pptsnow + pptgraul + pptice) / (rainncv(i,j) + R1)

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
                enddo

                !=================================================================================================================
                ! Reflectivity
                call calc_refl10cm (qv1d, qc1d, qr1d, nr1d, qs1d, qg1d, ng1d, qb1d, &
                    t1d, p1d, dBZ, kts, kte, i, j, configs)
                do k = kts, kte
                    refl_10cm(i,k,j) = max(-35.0_wp, dBZ(k))
                enddo

                ! Cloud, ice, and snow effective radius
                if (has_reqc /= 0 .and. has_reqi /= 0 .and. has_reqs /= 0) then
                    do k = kts, kte
                        re_qc1d(k) = 2.49e-6
                        re_qi1d(k) = 4.99e-6
                        re_qs1d(k) = 9.99e-6
                    enddo
                    call calc_effectRad (t1d=t1d, p1d=p1d, qv1d=qv1d, qc1d=qc1d, nc1d=nc1d, qi1d=qi1d, ni1d=ni1d, qs1d=qs1d,  &
                        re_qc1d=re_qc1d, re_qi1d=re_qi1d, re_qs1d=re_qs1d, kts=kts, kte=kte, configs=configs)
                    do k = kts, kte
                        re_cloud(i,k,j) = max(2.49e-6, min(re_qc1d(k), 50.e-6))
                        re_ice(i,k,j)   = max(4.99e-6, min(re_qi1d(k), 125.e-6))
                        re_snow(i,k,j)  = max(9.99e-6, min(re_qs1d(k), 999.e-6))
                    enddo
                endif

            enddo i_loop
        enddo j_loop

    end subroutine tempo_3d_to_1d_driver
    !=================================================================================================================

end module module_mp_tempo
