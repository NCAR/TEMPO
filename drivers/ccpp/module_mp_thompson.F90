! 3D Thompson-Eidhammer Microphysics Driver for MPAS
!=================================================================================================================
module module_mp_thompson

    use module_mp_thompson_params
    use module_mp_thompson_utils, only : wgamma, create_bins, table_Efrw, table_Efsw, table_dropEvap,
    table_ccnAct, qi_aut_qs, qr_acr_qg_par, qr_acr_qs_par, freezeH2O_par
    use module_mp_thompson_main, only : mp_thompson_main
    use mp_radar, only : radar_init

    implicit none

    type(config_flags) configs

contains
    !=================================================================================================================
    ! This subroutine handles initialzation of the microphysics scheme including building of lookup tables,
    ! allocating arrays for the microphysics scheme, and defining gamma function variables.
    SUBROUTINE thompson_init(is_aerosol_aware_in,       &
        merra2_aerosol_aware_in,   &
        mpicomm, mpirank, mpiroot, &
        threads, errmsg, errflg)

        IMPLICIT NONE

        LOGICAL, INTENT(IN) :: is_aerosol_aware_in
        LOGICAL, INTENT(IN) :: merra2_aerosol_aware_in
        TYPE(MPI_Comm), INTENT(IN) :: mpicomm
        INTEGER, INTENT(IN) :: mpirank, mpiroot
        INTEGER, INTENT(IN) :: threads
        CHARACTER(len=*), INTENT(INOUT) :: errmsg
        INTEGER,          INTENT(INOUT) :: errflg

        INTEGER:: i, j, k, l, m, n
        LOGICAL:: micro_init
        real :: stime, etime
        LOGICAL, PARAMETER :: precomputed_tables = .FALSE.

        ! Set module variable is_aerosol_aware/merra2_aerosol_aware
        is_aerosol_aware = is_aerosol_aware_in
        merra2_aerosol_aware = merra2_aerosol_aware_in
        if (is_aerosol_aware .and. merra2_aerosol_aware) then
            errmsg = 'Logic error in thompson_init: only one of the two options can be true, ' // &
                'not both: is_aerosol_aware or merra2_aerosol_aware'
            errflg = 1
            return
        end if
        if (mpirank==mpiroot) then
            if (is_aerosol_aware) then
                write (*,'(a)') 'Using aerosol-aware version of Thompson microphysics'
            else if(merra2_aerosol_aware) then
                write (*,'(a)') 'Using merra2 aerosol-aware version of Thompson microphysics'
            else
                write (*,'(a)') 'Using non-aerosol-aware version of Thompson microphysics'
            end if
        end if

        micro_init = .FALSE.

        av_g(idx_bg1) = av_g_old
        bv_g(idx_bg1) = bv_g_old

        !> - Allocate space for lookup tables (J. Michalakes 2009Jun08).

        if (.NOT. ALLOCATED(tcg_racg) ) then
            ALLOCATE(tcg_racg(ntb_g1,ntb_g,1,ntb_r1,ntb_r))
            micro_init = .TRUE.
        endif

        if (.NOT. ALLOCATED(tmr_racg)) ALLOCATE(tmr_racg(ntb_g1,ntb_g,1,ntb_r1,ntb_r))
        if (.NOT. ALLOCATED(tcr_gacr)) ALLOCATE(tcr_gacr(ntb_g1,ntb_g,1,ntb_r1,ntb_r))
        if (.NOT. ALLOCATED(tnr_racg)) ALLOCATE(tnr_racg(ntb_g1,ntb_g,1,ntb_r1,ntb_r))
        if (.NOT. ALLOCATED(tnr_gacr)) ALLOCATE(tnr_gacr(ntb_g1,ntb_g,1,ntb_r1,ntb_r))

        if (.NOT. ALLOCATED(tcs_racs1)) ALLOCATE(tcs_racs1(ntb_s,ntb_t,ntb_r1,ntb_r))
        if (.NOT. ALLOCATED(tmr_racs1)) ALLOCATE(tmr_racs1(ntb_s,ntb_t,ntb_r1,ntb_r))
        if (.NOT. ALLOCATED(tcs_racs2)) ALLOCATE(tcs_racs2(ntb_s,ntb_t,ntb_r1,ntb_r))
        if (.NOT. ALLOCATED(tmr_racs2)) ALLOCATE(tmr_racs2(ntb_s,ntb_t,ntb_r1,ntb_r))
        if (.NOT. ALLOCATED(tcr_sacr1)) ALLOCATE(tcr_sacr1(ntb_s,ntb_t,ntb_r1,ntb_r))
        if (.NOT. ALLOCATED(tms_sacr1)) ALLOCATE(tms_sacr1(ntb_s,ntb_t,ntb_r1,ntb_r))
        if (.NOT. ALLOCATED(tcr_sacr2)) ALLOCATE(tcr_sacr2(ntb_s,ntb_t,ntb_r1,ntb_r))
        if (.NOT. ALLOCATED(tms_sacr2)) ALLOCATE(tms_sacr2(ntb_s,ntb_t,ntb_r1,ntb_r))
        if (.NOT. ALLOCATED(tnr_racs1)) ALLOCATE(tnr_racs1(ntb_s,ntb_t,ntb_r1,ntb_r))
        if (.NOT. ALLOCATED(tnr_racs2)) ALLOCATE(tnr_racs2(ntb_s,ntb_t,ntb_r1,ntb_r))
        if (.NOT. ALLOCATED(tnr_sacr1)) ALLOCATE(tnr_sacr1(ntb_s,ntb_t,ntb_r1,ntb_r))
        if (.NOT. ALLOCATED(tnr_sacr2)) ALLOCATE(tnr_sacr2(ntb_s,ntb_t,ntb_r1,ntb_r))

        if (.NOT. ALLOCATED(tpi_qcfz)) ALLOCATE(tpi_qcfz(ntb_c,nbc,45,ntb_IN))
        if (.NOT. ALLOCATED(tni_qcfz)) ALLOCATE(tni_qcfz(ntb_c,nbc,45,ntb_IN))

        if (.NOT. ALLOCATED(tpi_qrfz)) ALLOCATE(tpi_qrfz(ntb_r,ntb_r1,45,ntb_IN))
        if (.NOT. ALLOCATED(tpg_qrfz)) ALLOCATE(tpg_qrfz(ntb_r,ntb_r1,45,ntb_IN))
        if (.NOT. ALLOCATED(tni_qrfz)) ALLOCATE(tni_qrfz(ntb_r,ntb_r1,45,ntb_IN))
        if (.NOT. ALLOCATED(tnr_qrfz)) ALLOCATE(tnr_qrfz(ntb_r,ntb_r1,45,ntb_IN))

        if (.NOT. ALLOCATED(tps_iaus)) ALLOCATE(tps_iaus(ntb_i,ntb_i1))
        if (.NOT. ALLOCATED(tni_iaus)) ALLOCATE(tni_iaus(ntb_i,ntb_i1))
        if (.NOT. ALLOCATED(tpi_ide)) ALLOCATE(tpi_ide(ntb_i,ntb_i1))

        if (.NOT. ALLOCATED(t_Efrw)) ALLOCATE(t_Efrw(nbr,nbc))
        if (.NOT. ALLOCATED(t_Efsw)) ALLOCATE(t_Efsw(nbs,nbc))

        if (.NOT. ALLOCATED(tnr_rev)) ALLOCATE(tnr_rev(nbr, ntb_r1, ntb_r))
        if (.NOT. ALLOCATED(tpc_wev)) ALLOCATE(tpc_wev(nbc,ntb_c,nbc))
        if (.NOT. ALLOCATED(tnc_wev)) ALLOCATE(tnc_wev(nbc,ntb_c,nbc))

        if (.NOT. ALLOCATED(tnccn_act))                                   &
            ALLOCATE(tnccn_act(ntb_arc,ntb_arw,ntb_art,ntb_arr,ntb_ark))

        if_micro_init: if (micro_init) then

            !> - From Martin et al. (1994), assign gamma shape parameter mu for cloud
            !! drops according to general dispersion characteristics (disp=~0.25
            !! for maritime and 0.45 for continental)
            !.. disp=SQRT((mu+2)/(mu+1) - 1) so mu varies from 15 for Maritime
            !.. to 2 for really dirty air.  This not used in 2-moment cloud water
            !.. scheme and nu_c used instead and varies from 2 to 15 (integer-only).
            mu_c_l = MIN(15., (1000.E6/Nt_c_l + 2.))
            mu_c_o = MIN(15., (1000.E6/Nt_c_o + 2.))

            !> - Compute Schmidt number to one-third used numerous times
            Sc3 = Sc**(1./3.)

            !> - Compute minimum ice diam from mass, min snow/graupel mass from diam
            D0i = (xm0i/am_i)**(1./bm_i)
            xm0s = am_s * D0s**bm_s
            xm0g = am_g * D0g**bm_g

            !> - Compute constants various exponents and gamma() associated with cloud,
            !! rain, snow, and graupel
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
            cse(17) = cse(16) + mu_s + 1.
            cse(18) = bv_s + mu_s + 3.
            do n = 1, 18
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
                cge(5,m) = bm_g*2. + mu_g + bv_g(idx_bg1) + 1.
                cge(6,m) = bm_g + mu_g + bv_g(idx_bg1) + 1.
                cge(7,m) = bm_g*0.5 + mu_g + bv_g(idx_bg1) + 1.
                cge(8,m) = mu_g + bv_g(idx_bg1) + 1.      ! not used
                cge(9,m) = mu_g + bv_g(idx_bg1) + 3.
                cge(11,m) = 0.5*(bv_g(idx_bg1) + 5. + 2.*mu_g)
            enddo

            do m = 1, NRHG
                do n = 1, 12
                    cgg(n,m) = wgamma(cge(n,m))
                enddo
            enddo

            do m = 1, NRHG
                oamg(m) = 1.0 / am_g(idx_bg1)
                ocmg(m) = oamg(m)**obmg
            enddo

            ! cge(1) = bm_g + 1.
            ! cge(2) = mu_g + 1.
            ! cge(3) = bm_g + mu_g + 1.
            ! cge(4) = bm_g*2. + mu_g + 1.
            ! cge(5) = bm_g*2. + mu_g + bv_g + 1.
            ! cge(6) = bm_g + mu_g + bv_g + 1.
            ! cge(7) = bm_g + mu_g + bv_g + 2.
            ! cge(8) = bm_g + mu_g + bv_g + 3.
            ! cge(9) = mu_g + bv_g + 3.
            ! cge(10) = mu_g + 2.
            ! cge(11) = 0.5*(bv_g + 5. + 2.*mu_g)
            ! cge(12) = 0.5*(bv_g + 5.) + mu_g
            ! do n = 1, 12
            !    cgg(n) = WGAMMA(cge(n))
            ! enddo

            oamg = 1./am_g
            obmg = 1./bm_g
            ocmg = oamg**obmg
            oge1 = 1./cge(1)
            ogg1 = 1./cgg(1)
            ogg2 = 1./cgg(2)
            ogg3 = 1./cgg(3)

            !+---+-----------------------------------------------------------------+
            !> - Simplify various rate equations
            !+---+-----------------------------------------------------------------+

            !>  - Compute rain collecting cloud water and cloud ice
            t1_qr_qc = PI*.25*av_r * crg(9)
            t1_qr_qi = PI*.25*av_r * crg(9)
            t2_qr_qi = PI*.25*am_r*av_r * crg(8)

            !>  - Compute graupel collecting cloud water
            t1_qg_qc = PI*.25*av_g * cgg(9)

            !>  - Compute snow collecting cloud water
            t1_qs_qc = PI*.25*av_s

            !>  - Compute snow collecting cloud ice
            t1_qs_qi = PI*.25*av_s

            !>  - Compute evaporation of rain; ignore depositional growth of rain
            t1_qr_ev = 0.78 * crg(10)
            t2_qr_ev = 0.308*Sc3*SQRT(av_r) * crg(11)

            !>  - Compute sublimation/depositional growth of snow
            t1_qs_sd = 0.86
            t2_qs_sd = 0.28*Sc3*SQRT(av_s)

            !>  - Compute melting of snow
            t1_qs_me = PI*4.*C_sqrd*olfus * 0.86
            t2_qs_me = PI*4.*C_sqrd*olfus * 0.28*Sc3*SQRT(av_s)

            !>  - Compute sublimation/depositional growth of graupel
            t1_qg_sd = 0.86 * cgg(10)
            t2_qg_sd = 0.28*Sc3*SQRT(av_g) * cgg(11)

            !>  - Compute melting of graupel
            t1_qg_me = PI*4.*C_cube*olfus * 0.86 * cgg(10)
            t2_qg_me = PI*4.*C_cube*olfus * 0.28*Sc3*SQRT(av_g) * cgg(11)

            !>  - Compute constants for helping find lookup table indexes
            nic2 = NINT(ALOG10(r_c(1)))
            nii2 = NINT(ALOG10(r_i(1)))
            nii3 = NINT(ALOG10(Nt_i(1)))
            nir2 = NINT(ALOG10(r_r(1)))
            nir3 = NINT(ALOG10(N0r_exp(1)))
            nis2 = NINT(ALOG10(r_s(1)))
            nig2 = NINT(ALOG10(r_g(1)))
            nig3 = NINT(ALOG10(N0g_exp(1)))
            niIN2 = NINT(ALOG10(Nt_IN(1)))

            !>  - Create bins of cloud water (from min diameter up to 100 microns)
            Dc(1) = D0c*1.0d0
            dtc(1) = D0c*1.0d0
            do n = 2, nbc
                Dc(n) = Dc(n-1) + 1.0D-6
                dtc(n) = (Dc(n) - Dc(n-1))
            enddo

            !>  - Create bins of cloud ice (from min diameter up to 5x min snow size)
            call create_bins(numbins=nbi, lowbin=D0i*1.0_dp, highbin=D0s*2.0_dp, &
                bins=Di, deltabins=dti)

            !>  - Create bins of rain (from min diameter up to 5 mm)
            call create_bins(numbins=nbr, lowbin=D0r*1.0_dp, highbin=0.005_dp, &
                bins=Dr, deltabins=dtr)

            !>  - Create bins of snow (from min diameter up to 2 cm)
            call create_bins(numbins=nbs, lowbin=D0s*1.0_dp, highbin=0.02_dp, &
                bins=Ds, deltabins=dts)

            !>  - Create bins of graupel (from min diameter up to 5 cm)
            call create_bins(numbins=nbg, lowbin=D0g*1.0_dp, highbin=0.05_dp, &
                bins=Dg, deltabins=dtg)

            !>  - Create bins of cloud droplet number concentration (1 to 3000 per cc)
            call create_bins(numbins=nbc, lowbin=1.0_dp, highbin=3000.0_dp, &
                bins=t_Nc)
            t_Nc = t_Nc * 1.0e6_dp
            nic1 = DLOG(t_Nc(nbc)/t_Nc(1))

            !+---+-----------------------------------------------------------------+
            !> - Create lookup tables for most costly calculations
            !+---+-----------------------------------------------------------------+

            ! Assign mpicomm to module variable
            mpi_communicator = mpicomm

            ! Standard tables are only written by master MPI task;
            ! (physics init cannot be called by multiple threads,
            !  hence no need to test for a specific thread number)
            if (mpirank==mpiroot) then
                thompson_table_writer = .true.
            else
                thompson_table_writer = .false.
            end if

            precomputed_tables_1: if (.not.precomputed_tables) then

                call cpu_time(stime)

                do m = 1, ntb_r
                    do k = 1, ntb_r1
                        do n = 1, 1
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

                if (mpirank==mpiroot) write (*,*)'creating microphysics lookup tables ... '
                if (mpirank==mpiroot) write (*,'(a, f5.2, a, f5.2, a, f5.2, a, f5.2)') &
                    ' using: mu_c_o=',mu_c_o,' mu_i=',mu_i,' mu_r=',mu_r,' mu_g=',mu_g

                !>  - Call table_ccnact() to read a static file containing CCN activation of aerosols. The
                !! data were created from a parcel model by Feingold & Heymsfield with
                !! further changes by Eidhammer and Kriedenweis
                if (mpirank==mpiroot) write(*,*) '  calling table_ccnAct routine'
                call table_ccnAct(errmsg,errflg)
                if (.not. errflg==0) return

                !>  - Call table_efrw() and table_efsw() to creat collision efficiency table
                !! between rain/snow and cloud water
                if (mpirank==mpiroot) write(*,*) '  creating qc collision eff tables'
                call table_Efrw
                call table_Efsw

                !>  - Call table_dropevap() to creat rain drop evaporation table
                if (mpirank==mpiroot) write(*,*) '  creating rain evap table'
                call table_dropEvap

                !>  - Call qi_aut_qs() to create conversion of some ice mass into snow category
                if (mpirank==mpiroot) write(*,*) '  creating ice converting to snow table'
                call qi_aut_qs

                call cpu_time(etime)
                if (mpirank==mpiroot) print '("Calculating Thompson tables part 1 took ",f10.3," seconds.")', etime-stime

            end if precomputed_tables_1

            !>  - Call radar_init() to initialize various constants for computing radar reflectivity
            call cpu_time(stime)
            xam_r = am_r
            xbm_r = bm_r
            xmu_r = mu_r
            xam_s = am_s
            xbm_s = bm_s
            xmu_s = mu_s
            xam_g = am_g
            xbm_g = bm_g
            xmu_g = mu_g
            call radar_init
            call cpu_time(etime)
            if (mpirank==mpiroot) print '("Calling radar_init took ",f10.3," seconds.")', etime-stime


            if_not_iiwarm: if (.not. iiwarm) then

                precomputed_tables_2: if (.not.precomputed_tables) then

                    call cpu_time(stime)

                    !>  - Call qr_acr_qg() to create rain collecting graupel & graupel collecting rain table
                    if (mpirank==mpiroot) write(*,*) '  creating rain collecting graupel table'
                    call cpu_time(stime)
                    call qr_acr_qg_par
                    call cpu_time(etime)
                    if (mpirank==mpiroot) print '("Computing rain collecting graupel table took ",f10.3," seconds.")', etime-stime

                    !>  - Call qr_acr_qs() to create rain collecting snow & snow collecting rain table
                    if (mpirank==mpiroot) write (*,*) '  creating rain collecting snow table'
                    call cpu_time(stime)
                    call qr_acr_qs_par
                    call cpu_time(etime)
                    if (mpirank==mpiroot) print '("Computing rain collecting snow table took ",f10.3," seconds.")', etime-stime

                    !>  - Call freezeh2o() to create cloud water and rain freezing (Bigg, 1953) table
                    if (mpirank==mpiroot) write(*,*) '  creating freezing of water drops table'
                    call cpu_time(stime)
                    call freezeH2O_par(threads)
                    call cpu_time(etime)
                    if (mpirank==mpiroot) print '("Computing freezing of water drops table took ",f10.3," seconds.")', etime-stime

                    call cpu_time(etime)
                    if (mpirank==mpiroot) print '("Calculating Thompson tables part 2 took ",f10.3," seconds.")', etime-stime

                end if precomputed_tables_2

            endif if_not_iiwarm

            if (mpirank==mpiroot) write(*,*) ' ... DONE microphysical lookup tables'

        endif if_micro_init

    end SUBROUTINE thompson_init

    !>\ingroup aathompson
    !!This is a wrapper routine designed to transfer values from 3D to 1D.
    !!\section gen_mpgtdriver Thompson mp_gt_driver General Algorithm
    !> @{
    SUBROUTINE mp_gt_driver(qv, qc, qr, qi, qs, qg, ni, nr, nc,     &
        nwfa, nifa, nwfa2d, nifa2d,             &
        tt, th, pii,                            &
        p, w, dz, dt_in, dt_inner,              &
        sedi_semi, decfl, lsm,                  &
        RAINNC, RAINNCV,                        &
        SNOWNC, SNOWNCV,                        &
        ICENC, ICENCV,                          &
        GRAUPELNC, GRAUPELNCV, SR,              &
#if ( WRF_CHEM == 1 )
        rainprod, evapprod,                     &
#endif
        refl_10cm, diagflag, do_radar_ref,      &
        max_hail_diam_sfc,                      &
        vt_dbz_wt, first_time_step,             &
        re_cloud, re_ice, re_snow,              &
        has_reqc, has_reqi, has_reqs,           &
        aero_ind_fdb, rand_perturb_on,          &
        kme_stoch,                              &
        rand_pert, spp_prt_list, spp_var_list,  &
        spp_stddev_cutoff, n_var_spp,           &
        ids,ide, jds,jde, kds,kde,              &  ! domain dims
        ims,ime, jms,jme, kms,kme,              &  ! memory dims
        its,ite, jts,jte, kts,kte,              &  ! tile dims
        fullradar_diag, istep, nsteps,          &
        errmsg, errflg,                         &
    ! Extended diagnostics, array pointers
    ! only associated if ext_diag flag is .true.
        ext_diag,                               &
    !vts1, txri, txrc,                       &
        prw_vcdc,                               &
        prw_vcde, tpri_inu, tpri_ide_d,         &
        tpri_ide_s, tprs_ide, tprs_sde_d,       &
        tprs_sde_s, tprg_gde_d,                 &
        tprg_gde_s, tpri_iha, tpri_wfz,         &
        tpri_rfz, tprg_rfz, tprs_scw, tprg_scw, &
        tprg_rcs, tprs_rcs,                     &
        tprr_rci, tprg_rcg,                     &
        tprw_vcd_c, tprw_vcd_e, tprr_sml,       &
        tprr_gml, tprr_rcg,                     &
        tprr_rcs, tprv_rev, tten3, qvten3,      &
        qrten3, qsten3, qgten3, qiten3, niten3, &
        nrten3, ncten3, qcten3,                 &
        pfils, pflls)

        implicit none

        !..Subroutine arguments
        INTEGER, INTENT(IN):: ids,ide, jds,jde, kds,kde, &
            ims,ime, jms,jme, kms,kme, &
            its,ite, jts,jte, kts,kte
        REAL, DIMENSION(ims:ime, kms:kme, jms:jme), INTENT(INOUT):: &
            qv, qc, qr, qi, qs, qg, ni, nr
        REAL, DIMENSION(ims:ime, kms:kme, jms:jme), OPTIONAL, INTENT(INOUT):: &
            tt, th
        REAL, DIMENSION(ims:ime, kms:kme, jms:jme), OPTIONAL, INTENT(IN):: &
            pii
        REAL, DIMENSION(ims:ime, kms:kme, jms:jme), OPTIONAL, INTENT(INOUT):: &
            nc, nwfa, nifa
        REAL, DIMENSION(ims:ime, jms:jme), OPTIONAL, INTENT(IN):: nwfa2d, nifa2d
        INTEGER, DIMENSION(ims:ime, jms:jme), INTENT(IN):: lsm
        REAL, DIMENSION(ims:ime, kms:kme, jms:jme), OPTIONAL, INTENT(INOUT):: &
            re_cloud, re_ice, re_snow
        REAL, DIMENSION(ims:ime, kms:kme, jms:jme), INTENT(INOUT):: pfils, pflls
        INTEGER, INTENT(IN) :: rand_perturb_on, kme_stoch, n_var_spp
        REAL, DIMENSION(:,:), INTENT(IN) :: rand_pert
        REAL, DIMENSION(:), INTENT(IN) :: spp_prt_list, spp_stddev_cutoff
        CHARACTER(len=10), DIMENSION(:), INTENT(IN) :: spp_var_list
        INTEGER, INTENT(IN):: has_reqc, has_reqi, has_reqs
#if ( WRF_CHEM == 1 )
        REAL, DIMENSION(ims:ime, kms:kme, jms:jme), INTENT(INOUT):: &
            rainprod, evapprod
#endif
        REAL, DIMENSION(ims:ime, kms:kme, jms:jme), INTENT(IN):: &
            p, w, dz
        REAL, DIMENSION(ims:ime, jms:jme), INTENT(INOUT):: &
            RAINNC, RAINNCV, SR
        REAL, DIMENSION(ims:ime, jms:jme), OPTIONAL, INTENT(INOUT)::      &
            SNOWNC, SNOWNCV,                              &
            ICENC, ICENCV,                                &
            GRAUPELNC, GRAUPELNCV
        REAL, DIMENSION(ims:ime, kms:kme, jms:jme), INTENT(INOUT)::       &
            refl_10cm
        REAL, DIMENSION(ims:ime, jms:jme), INTENT(INOUT)::       &
            max_hail_diam_sfc
        REAL, DIMENSION(ims:ime, kms:kme, jms:jme), OPTIONAL, INTENT(INOUT):: &
            vt_dbz_wt
        LOGICAL, INTENT(IN) :: first_time_step
        REAL, INTENT(IN):: dt_in, dt_inner
        LOGICAL, INTENT(IN) :: sedi_semi
        INTEGER, INTENT(IN) :: decfl
        ! To support subcycling: current step and maximum number of steps
        INTEGER, INTENT (IN) :: istep, nsteps
        LOGICAL, INTENT (IN) :: fullradar_diag
        ! Extended diagnostics, array pointers only associated if ext_diag flag is .true.
        LOGICAL, INTENT (IN) :: ext_diag
        LOGICAL, OPTIONAL, INTENT(IN):: aero_ind_fdb
        REAL, DIMENSION(:,:,:), INTENT(INOUT)::                     &
        !vts1, txri, txrc,                       &
            prw_vcdc,                               &
            prw_vcde, tpri_inu, tpri_ide_d,         &
            tpri_ide_s, tprs_ide,                   &
            tprs_sde_d, tprs_sde_s, tprg_gde_d,     &
            tprg_gde_s, tpri_iha, tpri_wfz,         &
            tpri_rfz, tprg_rfz, tprs_scw, tprg_scw, &
            tprg_rcs, tprs_rcs,                     &
            tprr_rci, tprg_rcg,                     &
            tprw_vcd_c, tprw_vcd_e, tprr_sml,       &
            tprr_gml, tprr_rcg,                     &
            tprr_rcs, tprv_rev, tten3, qvten3,      &
            qrten3, qsten3, qgten3, qiten3, niten3, &
            nrten3, ncten3, qcten3

        !..Local variables
        REAL, DIMENSION(kts:kte):: &
            qv1d, qc1d, qi1d, qr1d, qs1d, qg1d, ni1d,     &
            nr1d, nc1d, nwfa1d, nifa1d,                   &
            t1d, p1d, w1d, dz1d, rho, dBZ, pfil1, pfll1
        !..Extended diagnostics, single column arrays
        REAL, DIMENSION(:), ALLOCATABLE::                              &
        !vtsk1, txri1, txrc1,                       &
            prw_vcdc1,                                 &
            prw_vcde1, tpri_inu1, tpri_ide1_d,         &
            tpri_ide1_s, tprs_ide1,                    &
            tprs_sde1_d, tprs_sde1_s, tprg_gde1_d,     &
            tprg_gde1_s, tpri_iha1, tpri_wfz1,         &
            tpri_rfz1, tprg_rfz1, tprs_scw1, tprg_scw1,&
            tprg_rcs1, tprs_rcs1,                      &
            tprr_rci1, tprg_rcg1,                      &
            tprw_vcd1_c, tprw_vcd1_e, tprr_sml1,       &
            tprr_gml1, tprr_rcg1,                      &
            tprr_rcs1, tprv_rev1,  tten1, qvten1,      &
            qrten1, qsten1, qgten1, qiten1, niten1,    &
            nrten1, ncten1, qcten1

        REAL, DIMENSION(kts:kte):: re_qc1d, re_qi1d, re_qs1d
#if ( WRF_CHEM == 1 )
        REAL, DIMENSION(kts:kte):: &
            rainprod1d, evapprod1d
#endif
        REAL, DIMENSION(its:ite, jts:jte):: pcp_ra, pcp_sn, pcp_gr, pcp_ic
        REAL:: dt, pptrain, pptsnow, pptgraul, pptice
        REAL:: qc_max, qr_max, qs_max, qi_max, qg_max, ni_max, nr_max
        INTEGER:: lsml
        REAL:: rand1, rand2, rand3, rand_pert_max
        INTEGER:: i, j, k, m
        INTEGER:: imax_qc,imax_qr,imax_qi,imax_qs,imax_qg,imax_ni,imax_nr
        INTEGER:: jmax_qc,jmax_qr,jmax_qi,jmax_qs,jmax_qg,jmax_ni,jmax_nr
        INTEGER:: kmax_qc,kmax_qr,kmax_qi,kmax_qs,kmax_qg,kmax_ni,kmax_nr
        INTEGER:: i_start, j_start, i_end, j_end
        LOGICAL, OPTIONAL, INTENT(IN) :: diagflag
        INTEGER, OPTIONAL, INTENT(IN) :: do_radar_ref
        logical :: melti = .false.
        INTEGER :: ndt, it

        ! CCPP error handling
        character(len=*), optional, intent(  out) :: errmsg
        integer,          optional, intent(  out) :: errflg

        ! CCPP
        if (present(errmsg)) errmsg = ''
        if (present(errflg)) errflg = 0

        ! No need to test for every subcycling step
        test_only_once: if (first_time_step .and. istep==1) then
            ! Activate this code when removing the guard above

            if ( (present(tt) .and. (present(th) .or. present(pii))) .or. &
                (.not.present(tt) .and. .not.(present(th) .and. present(pii))) ) then
                if (present(errmsg) .and. present(errflg)) then
                    write(errmsg, '(a)') 'Logic error in mp_gt_driver: provide either tt or th+pii'
                    errflg = 1
                    return
                else
                    write(*,'(a)') 'Logic error in mp_gt_driver: provide either tt or th+pii'
                    stop
                end if
            end if

            if (is_aerosol_aware .and. (.not.present(nc)     .or. &
                .not.present(nwfa)   .or. &
                .not.present(nifa)   .or. &
                .not.present(nwfa2d) .or. &
                .not.present(nifa2d)      )) then
                if (present(errmsg) .and. present(errflg)) then
                    write(errmsg, '(*(a))') 'Logic error in mp_gt_driver: provide nc, nwfa, nifa, nwfa2d', &
                        ' and nifa2d for aerosol-aware version of Thompson microphysics'
                    errflg = 1
                    return
                else
                    write(*, '(*(a))') 'Logic error in mp_gt_driver: provide nc, nwfa, nifa, nwfa2d', &
                        ' and nifa2d for aerosol-aware version of Thompson microphysics'
                    stop
                end if
            else if (merra2_aerosol_aware .and. (.not.present(nc)   .or. &
                .not.present(nwfa) .or. &
                .not.present(nifa)      )) then
                if (present(errmsg) .and. present(errflg)) then
                    write(errmsg, '(*(a))') 'Logic error in mp_gt_driver: provide nc, nwfa, and nifa', &
                        ' for merra2 aerosol-aware version of Thompson microphysics'
                    errflg = 1
                    return
                else
                    write(*, '(*(a))') 'Logic error in mp_gt_driver: provide nc, nwfa, and nifa', &
                        ' for merra2 aerosol-aware version of Thompson microphysics'
                    stop
                end if
            else if (.not.is_aerosol_aware .and. .not.merra2_aerosol_aware .and. &
                (present(nwfa) .or. present(nifa) .or. present(nwfa2d) .or. present(nifa2d))) then
                write(*,*) 'WARNING, nc/nwfa/nifa/nwfa2d/nifa2d present but is_aerosol_aware/merra2_aerosol_aware are FALSE'
            end if
        end if test_only_once

        ! These must be alwyas allocated
        !allocate (vtsk1(kts:kte))
        !allocate (txri1(kts:kte))
        !allocate (txrc1(kts:kte))
        allocate_extended_diagnostics: if (ext_diag) then
            allocate (prw_vcdc1(kts:kte))
            allocate (prw_vcde1(kts:kte))
            allocate (tpri_inu1(kts:kte))
            allocate (tpri_ide1_d(kts:kte))
            allocate (tpri_ide1_s(kts:kte))
            allocate (tprs_ide1(kts:kte))
            allocate (tprs_sde1_d(kts:kte))
            allocate (tprs_sde1_s(kts:kte))
            allocate (tprg_gde1_d(kts:kte))
            allocate (tprg_gde1_s(kts:kte))
            allocate (tpri_iha1(kts:kte))
            allocate (tpri_wfz1(kts:kte))
            allocate (tpri_rfz1(kts:kte))
            allocate (tprg_rfz1(kts:kte))
            allocate (tprs_scw1(kts:kte))
            allocate (tprg_scw1(kts:kte))
            allocate (tprg_rcs1(kts:kte))
            allocate (tprs_rcs1(kts:kte))
            allocate (tprr_rci1(kts:kte))
            allocate (tprg_rcg1(kts:kte))
            allocate (tprw_vcd1_c(kts:kte))
            allocate (tprw_vcd1_e(kts:kte))
            allocate (tprr_sml1(kts:kte))
            allocate (tprr_gml1(kts:kte))
            allocate (tprr_rcg1(kts:kte))
            allocate (tprr_rcs1(kts:kte))
            allocate (tprv_rev1(kts:kte))
            allocate (tten1(kts:kte))
            allocate (qvten1(kts:kte))
            allocate (qrten1(kts:kte))
            allocate (qsten1(kts:kte))
            allocate (qgten1(kts:kte))
            allocate (qiten1(kts:kte))
            allocate (niten1(kts:kte))
            allocate (nrten1(kts:kte))
            allocate (ncten1(kts:kte))
            allocate (qcten1(kts:kte))
        end if allocate_extended_diagnostics

        !+---+
        i_start = its
        j_start = jts
        i_end   = ite
        j_end   = jte

        !..For idealized testing by developer.
        !     if ( (ide-ids+1).gt.4 .and. (jde-jds+1).lt.4 .and.                &
        !          ids.eq.its.and.ide.eq.ite.and.jds.eq.jts.and.jde.eq.jte) then
        !        i_start = its + 2
        !        i_end   = ite
        !        j_start = jts
        !        j_end   = jte
        !     endif

        !     dt = dt_in
        RAINNC(:,:) = 0.0
        SNOWNC(:,:) = 0.0
        ICENC(:,:) = 0.0
        GRAUPELNC(:,:) = 0.0
        pcp_ra(:,:) = 0.0
        pcp_sn(:,:) = 0.0
        pcp_gr(:,:) = 0.0
        pcp_ic(:,:) = 0.0
        pfils(:,:,:) = 0.0
        pflls(:,:,:) = 0.0
        rand_pert_max = 0.0
        ndt = max(nint(dt_in/dt_inner),1)
        dt = dt_in/ndt
        if(dt_in .le. dt_inner) dt= dt_in

        !Get the Thompson MP SPP magnitude and standard deviation cutoff,
        !then compute rand_pert_max

        if (rand_perturb_on .ne. 0) then
            do k =1,n_var_spp
                select case (spp_var_list(k))
                  case('mp')
                    rand_pert_max = spp_prt_list(k)*spp_stddev_cutoff(k)
                end select
            enddo
        endif

        do it = 1, ndt

            qc_max = 0.
            qr_max = 0.
            qs_max = 0.
            qi_max = 0.
            qg_max = 0
            ni_max = 0.
            nr_max = 0.
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

            j_loop:  do j = j_start, j_end
                i_loop:  do i = i_start, i_end

                    !+---+-----------------------------------------------------------------+
                    !..Introduce stochastic parameter perturbations by creating as many scalar rand1, rand2, ...
                    !.. variables as needed to perturb different pieces of microphysics. gthompsn  21Mar2018
                    ! Setting spp_mp_opt to 1 gives graupel Y-intercept pertubations (2^0)
                    !                   2 gives cloud water distribution gamma shape parameter perturbations (2^1)
                    !                   4 gives CCN & IN activation perturbations (2^2)
                    !                   3 gives both 1+2
                    !                   5 gives both 1+4
                    !                   6 gives both 2+4
                    !                   7 gives all 1+2+4
                    ! For now (22Mar2018), standard deviation should be up to 0.75 and cut-off at 3.0
                    ! stddev in order to constrain the various perturbations from being too extreme.
                    !+---+-----------------------------------------------------------------+
                    rand1 = 0.0
                    rand2 = 0.0
                    rand3 = 0.0
                    if (rand_perturb_on .ne. 0) then
                        if (MOD(rand_perturb_on,2) .ne. 0) rand1 = rand_pert(i,1)
                        m = RSHIFT(ABS(rand_perturb_on),1)
                        if (MOD(m,2) .ne. 0) rand2 = rand_pert(i,1)*2.
                        m = RSHIFT(ABS(rand_perturb_on),2)
                        if (MOD(m,2) .ne. 0) rand3 = 0.25*(rand_pert(i,1)+rand_pert_max)
                        m = RSHIFT(ABS(rand_perturb_on),3)
                    endif
                    !+---+-----------------------------------------------------------------+

                    pptrain = 0.
                    pptsnow = 0.
                    pptgraul = 0.
                    pptice = 0.
                    RAINNCV(i,j) = 0.
                    IF ( PRESENT (snowncv) ) THEN
                        SNOWNCV(i,j) = 0.
                    ENDIF
                    IF ( PRESENT (icencv) ) THEN
                        ICENCV(i,j) = 0.
                    ENDIF
                    IF ( PRESENT (graupelncv) ) THEN
                        GRAUPELNCV(i,j) = 0.
                    ENDIF
                    SR(i,j) = 0.

                    do k = kts, kte
                        if (present(tt)) then
                            t1d(k) = tt(i,k,j)
                        else
                            t1d(k) = th(i,k,j)*pii(i,k,j)
                        end if
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
                        rho(k) = 0.622*p1d(k)/(R*t1d(k)*(qv1d(k)+0.622))

                        ! These arrays are always allocated and must be initialized
                        !vtsk1(k) = 0.
                        !txrc1(k) = 0.
                        !txri1(k) = 0.
                        initialize_extended_diagnostics: if (ext_diag) then
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
                            tprs_scw1(k) = 0.
                            tprg_scw1(k) = 0.
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
                        endif initialize_extended_diagnostics
                    enddo
                    lsml = lsm(i,j)
                    if (configs%aerosol_aware .or. merra2_aerosol_aware) then
                        do k = kts, kte
                            nc1d(k) = nc(i,k,j)
                            nwfa1d(k) = nwfa(i,k,j)
                            nifa1d(k) = nifa(i,k,j)
                        enddo
                    else
                        do k = kts, kte
                            if(lsml == 1) then
                                nc1d(k) = Nt_c_l/rho(k)
                            else
                                nc1d(k) = Nt_c_o/rho(k)
                            endif
                            nwfa1d(k) = 11.1E6
                            nifa1d(k) = naIN1*0.01
                        enddo
                    endif

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

                    !> - Call mp_thompson()
                    call mp_thompson_main(qv1d, qc1d, qi1d, qr1d, qs1d, qg1d, qb1d, ni1d,     &
                        nr1d, nc1d, ng1d, nwfa1d, nifa1d, t1d, p1d, w1d, dz1d,  &
                        pptrain, pptsnow, pptgraul, pptice, &
#if ( WRF_CHEM == 1 )
                        rainprod1d, evapprod1d, &
#endif
                        rand1, rand2, rand3, &
                        ext_diag,                    &
                        sedi_semi, decfl,                                &
                        prw_vcdc1, prw_vcde1,                            &
                        tpri_inu1, tpri_ide1_d, tpri_ide1_s, tprs_ide1,  &
                        tprs_sde1_d, tprs_sde1_s,                        &
                        tprg_gde1_d, tprg_gde1_s, tpri_iha1, tpri_wfz1,  &
                        tpri_rfz1, tprg_rfz1, tprs_scw1, tprg_scw1,      &
                        tprg_rcs1, tprs_rcs1, tprr_rci1,                 &
                        tprg_rcg1, tprw_vcd1_c,                          &
                        tprw_vcd1_e, tprr_sml1, tprr_gml1, tprr_rcg1,    &
                        tprr_rcs1, tprv_rev1,                            &
                        tten1, qvten1, qrten1, qsten1,                   &
                        qgten1, qiten1, niten1, nrten1, ncten1, qcten1,  &
                        pfil1, pfll1, lsml, &
                        kts, kte, dt, i, j, configs)

                    pcp_ra(i,j) = pcp_ra(i,j) + pptrain
                    pcp_sn(i,j) = pcp_sn(i,j) + pptsnow
                    pcp_gr(i,j) = pcp_gr(i,j) + pptgraul
                    pcp_ic(i,j) = pcp_ic(i,j) + pptice
                    RAINNCV(i,j) = pptrain + pptsnow + pptgraul + pptice
                    RAINNC(i,j) = RAINNC(i,j) + pptrain + pptsnow + pptgraul + pptice
                    IF ( PRESENT(snowncv) .AND. PRESENT(snownc) ) THEN
                        ! Add ice to snow if separate ice not present
                        IF ( .NOT.PRESENT(icencv) .OR. .NOT.PRESENT(icenc) ) THEN
                            SNOWNCV(i,j) = pptsnow + pptice
                            SNOWNC(i,j) = SNOWNC(i,j) + pptsnow + pptice
                        ELSE
                            SNOWNCV(i,j) = pptsnow
                            SNOWNC(i,j) = SNOWNC(i,j) + pptsnow
                        ENDIF
                    ENDIF
                    ! Use separate ice if present (as in FV3)
                    IF ( PRESENT(icencv) .AND. PRESENT(icenc) ) THEN
                        ICENCV(i,j) = pptice
                        ICENC(i,j) = ICENC(i,j) + pptice
                    ENDIF
                    IF ( PRESENT(graupelncv) .AND. PRESENT(graupelnc) ) THEN
                        GRAUPELNCV(i,j) = pptgraul
                        GRAUPELNC(i,j) = GRAUPELNC(i,j) + pptgraul
                    ENDIF
                    SR(i,j) = (pptsnow + pptgraul + pptice)/(RAINNCV(i,j)+1.e-12)



                    !..Reset lowest model level to initial state aerosols (fake sfc source).
                    !.. Changed 13 May 2013 to fake emissions in which nwfa2d is aerosol
                    !.. number tendency (number per kg per second).
                    if (configs%aerosol_aware) then
                        if ( PRESENT (aero_ind_fdb) ) then
                            if ( .not. aero_ind_fdb) then
                                nwfa1d(kts) = nwfa1d(kts) + nwfa2d(i,j)*dt
                                nifa1d(kts) = nifa1d(kts) + nifa2d(i,j)*dt
                            endif
                        else
                            nwfa1d(kts) = nwfa1d(kts) + nwfa2d(i,j)*dt
                            nifa1d(kts) = nifa1d(kts) + nifa2d(i,j)*dt
                        end if

                        do k = kts, kte
                            nc(i,k,j) = nc1d(k)
                            nwfa(i,k,j) = nwfa1d(k)
                            nifa(i,k,j) = nifa1d(k)
                        enddo
                    endif

                    if (merra2_aerosol_aware) then
                        do k = kts, kte
                            nc(i,k,j) = nc1d(k)
                            nwfa(i,k,j) = nwfa1d(k)
                            nifa(i,k,j) = nifa1d(k)
                        enddo
                    endif

                    if ((present(ng)) .and. (present(qb))) then
                        do k = kts, kte
                            ng(i,k,j) = ng1d(k)
                            qb(i,k,j) = qb1d(k)
                        enddo
                    endif

                    do k = kts, kte
                        qv(i,k,j) = qv1d(k)
                        qc(i,k,j) = qc1d(k)
                        qi(i,k,j) = qi1d(k)
                        qr(i,k,j) = qr1d(k)
                        qs(i,k,j) = qs1d(k)
                        qg(i,k,j) = qg1d(k)
                        ni(i,k,j) = ni1d(k)
                        nr(i,k,j) = nr1d(k)
                        pfils(i,k,j) = pfils(i,k,j) + pfil1(k)
                        pflls(i,k,j) = pflls(i,k,j) + pfll1(k)
                        if (present(tt)) then
                            tt(i,k,j) = t1d(k)
                        else
                            th(i,k,j) = t1d(k)/pii(i,k,j)
                        end if
#if ( WRF_CHEM == 1 )
                        rainprod(i,k,j) = rainprod1d(k)
                        evapprod(i,k,j) = evapprod1d(k)
#endif
                        if (qc1d(k) .gt. qc_max) then
                            imax_qc = i
                            jmax_qc = j
                            kmax_qc = k
                            qc_max = qc1d(k)
                        elseif (qc1d(k) .lt. 0.0) then
                            write(*,'(a,e16.7,a,3i8)') 'WARNING, negative qc ', qc1d(k),        &
                                ' at i,j,k=', i,j,k
                        endif
                        if (qr1d(k) .gt. qr_max) then
                            imax_qr = i
                            jmax_qr = j
                            kmax_qr = k
                            qr_max = qr1d(k)
                        elseif (qr1d(k) .lt. 0.0) then
                            write(*,'(a,e16.7,a,3i8)') 'WARNING, negative qr ', qr1d(k),        &
                                ' at i,j,k=', i,j,k
                        endif
                        if (nr1d(k) .gt. nr_max) then
                            imax_nr = i
                            jmax_nr = j
                            kmax_nr = k
                            nr_max = nr1d(k)
                        elseif (nr1d(k) .lt. 0.0) then
                            write(*,'(a,e16.7,a,3i8)') 'WARNING, negative nr ', nr1d(k),        &
                                ' at i,j,k=', i,j,k
                        endif
                        if (qs1d(k) .gt. qs_max) then
                            imax_qs = i
                            jmax_qs = j
                            kmax_qs = k
                            qs_max = qs1d(k)
                        elseif (qs1d(k) .lt. 0.0) then
                            write(*,'(a,e16.7,a,3i8)') 'WARNING, negative qs ', qs1d(k),        &
                                ' at i,j,k=', i,j,k
                        endif
                        if (qi1d(k) .gt. qi_max) then
                            imax_qi = i
                            jmax_qi = j
                            kmax_qi = k
                            qi_max = qi1d(k)
                        elseif (qi1d(k) .lt. 0.0) then
                            write(*,'(a,e16.7,a,3i8)') 'WARNING, negative qi ', qi1d(k),        &
                                ' at i,j,k=', i,j,k
                        endif
                        if (qg1d(k) .gt. qg_max) then
                            imax_qg = i
                            jmax_qg = j
                            kmax_qg = k
                            qg_max = qg1d(k)
                        elseif (qg1d(k) .lt. 0.0) then
                            write(*,'(a,e16.7,a,3i8)') 'WARNING, negative qg ', qg1d(k),        &
                                ' at i,j,k=', i,j,k
                        endif
                        if (ni1d(k) .gt. ni_max) then
                            imax_ni = i
                            jmax_ni = j
                            kmax_ni = k
                            ni_max = ni1d(k)
                        elseif (ni1d(k) .lt. 0.0) then
                            write(*,'(a,e16.7,a,3i8)') 'WARNING, negative ni ', ni1d(k),        &
                                ' at i,j,k=', i,j,k
                        endif
                        if (qv1d(k) .lt. 0.0) then
                            write(*,'(a,e16.7,a,3i8)') 'WARNING, negative qv ', qv1d(k),        &
                                ' at i,j,k=', i,j,k
                            if (k.lt.kte-2 .and. k.gt.kts+1) then
                                write(*,*) '   below and above are: ', qv(i,k-1,j), qv(i,k+1,j)
                                qv(i,k,j) = MAX(1.E-7, 0.5*(qv(i,k-1,j) + qv(i,k+1,j)))
                            else
                                qv(i,k,j) = 1.E-7
                            endif
                        endif
                    enddo

                    assign_extended_diagnostics: if (ext_diag) then
                        do k=kts,kte
                            !vts1(i,k,j)       = vtsk1(k)
                            !txri(i,k,j)       = txri(i,k,j)       + txri1(k)
                            !txrc(i,k,j)       = txrc(i,k,j)       + txrc1(k)
                            prw_vcdc(i,k,j)   = prw_vcdc(i,k,j)   + prw_vcdc1(k)
                            prw_vcde(i,k,j)   = prw_vcde(i,k,j)   + prw_vcde1(k)
                            tpri_inu(i,k,j)   = tpri_inu(i,k,j)   + tpri_inu1(k)
                            tpri_ide_d(i,k,j) = tpri_ide_d(i,k,j) + tpri_ide1_d(k)
                            tpri_ide_s(i,k,j) = tpri_ide_s(i,k,j) + tpri_ide1_s(k)
                            tprs_ide(i,k,j)   = tprs_ide(i,k,j)   + tprs_ide1(k)
                            tprs_sde_s(i,k,j) = tprs_sde_s(i,k,j) + tprs_sde1_s(k)
                            tprs_sde_d(i,k,j) = tprs_sde_d(i,k,j) + tprs_sde1_d(k)
                            tprg_gde_d(i,k,j) = tprg_gde_d(i,k,j) + tprg_gde1_d(k)
                            tprg_gde_s(i,k,j) = tprg_gde_s(i,k,j) + tprg_gde1_s(k)
                            tpri_iha(i,k,j)   = tpri_iha(i,k,j)   + tpri_iha1(k)
                            tpri_wfz(i,k,j)   = tpri_wfz(i,k,j)   + tpri_wfz1(k)
                            tpri_rfz(i,k,j)   = tpri_rfz(i,k,j)   + tpri_rfz1(k)
                            tprg_rfz(i,k,j)   = tprg_rfz(i,k,j)   + tprg_rfz1(k)
                            tprs_scw(i,k,j)   = tprs_scw(i,k,j)   + tprs_scw1(k)
                            tprg_scw(i,k,j)   = tprg_scw(i,k,j)   + tprg_scw1(k)
                            tprg_rcs(i,k,j)   = tprg_rcs(i,k,j)   + tprg_rcs1(k)
                            tprs_rcs(i,k,j)   = tprs_rcs(i,k,j)   + tprs_rcs1(k)
                            tprr_rci(i,k,j)   = tprr_rci(i,k,j)   + tprr_rci1(k)
                            tprg_rcg(i,k,j)   = tprg_rcg(i,k,j)   + tprg_rcg1(k)
                            tprw_vcd_c(i,k,j) = tprw_vcd_c(i,k,j) + tprw_vcd1_c(k)
                            tprw_vcd_e(i,k,j) = tprw_vcd_e(i,k,j) + tprw_vcd1_e(k)
                            tprr_sml(i,k,j)   = tprr_sml(i,k,j)   + tprr_sml1(k)
                            tprr_gml(i,k,j)   = tprr_gml(i,k,j)   + tprr_gml1(k)
                            tprr_rcg(i,k,j)   = tprr_rcg(i,k,j)   + tprr_rcg1(k)
                            tprr_rcs(i,k,j)   = tprr_rcs(i,k,j)   + tprr_rcs1(k)
                            tprv_rev(i,k,j)   = tprv_rev(i,k,j)   + tprv_rev1(k)
                            tten3(i,k,j)      = tten3(i,k,j)      + tten1(k)
                            qvten3(i,k,j)     = qvten3(i,k,j)     + qvten1(k)
                            qrten3(i,k,j)     = qrten3(i,k,j)     + qrten1(k)
                            qsten3(i,k,j)     = qsten3(i,k,j)     + qsten1(k)
                            qgten3(i,k,j)     = qgten3(i,k,j)     + qgten1(k)
                            qiten3(i,k,j)     = qiten3(i,k,j)     + qiten1(k)
                            niten3(i,k,j)     = niten3(i,k,j)     + niten1(k)
                            nrten3(i,k,j)     = nrten3(i,k,j)     + nrten1(k)
                            ncten3(i,k,j)     = ncten3(i,k,j)     + ncten1(k)
                            qcten3(i,k,j)     = qcten3(i,k,j)     + qcten1(k)

                        enddo
                    endif assign_extended_diagnostics

                    if (ndt>1 .and. it==ndt) then

                        SR(i,j) = (pcp_sn(i,j) + pcp_gr(i,j) + pcp_ic(i,j))/(RAINNC(i,j)+1.e-12)
                        RAINNCV(i,j) = RAINNC(i,j)
                        IF ( PRESENT (snowncv) ) THEN
                            SNOWNCV(i,j) = SNOWNC(i,j)
                        ENDIF
                        IF ( PRESENT (icencv) ) THEN
                            ICENCV(i,j) = ICENC(i,j)
                        ENDIF
                        IF ( PRESENT (graupelncv) ) THEN
                            GRAUPELNCV(i,j) = GRAUPELNC(i,j)
                        ENDIF
                    endif

                    ! Diagnostic calculations only for last step
                    ! if Thompson MP is called multiple times
                    last_step_only: IF ((ndt>1 .and. it==ndt) .or. &
                        (nsteps>1 .and. istep==nsteps) .or. &
                        (nsteps==1 .and. ndt==1)) THEN

                        max_hail_diam_sfc(i,j) = hail_mass_99th_percentile(kts, kte, qg1d, t1d, p1d, qv1d)

                        !> - Call calc_refl10cm()

                        diagflag_present: IF ( PRESENT (diagflag) ) THEN
                            if (diagflag .and. do_radar_ref == 1) then
                                !
                                ! Only set melti to true at the output times
                                if (fullradar_diag) then
                                    melti=.true.
                                else
                                    melti=.false.
                                endif
                                !
                                if (present(vt_dbz_wt)) then
                                    call calc_refl10cm (qv1d, qc1d, qr1d, nr1d, qs1d, qg1d,   &
                                        t1d, p1d, dBZ, rand1, kts, kte, i, j, &
                                        melti, vt_dbz_wt(i,:,j),              &
                                        first_time_step)
                                else
                                    call calc_refl10cm (qv1d, qc1d, qr1d, nr1d, qs1d, qg1d,   &
                                        t1d, p1d, dBZ, rand1, kts, kte, i, j, &
                                        melti)
                                end if
                                do k = kts, kte
                                    refl_10cm(i,k,j) = MAX(-35., dBZ(k))
                                enddo
                            endif
                        ENDIF diagflag_present

                        IF (has_reqc.ne.0 .and. has_reqi.ne.0 .and. has_reqs.ne.0) THEN
                            do k = kts, kte
                                re_qc1d(k) = re_qc_min
                                re_qi1d(k) = re_qi_min
                                re_qs1d(k) = re_qs_min
                            enddo
                            !> - Call calc_effectrad()
                            call calc_effectRad (t1d, p1d, qv1d, qc1d, nc1d, qi1d, ni1d, qs1d,  &
                                re_qc1d, re_qi1d, re_qs1d, lsml, kts, kte)
                            do k = kts, kte
                                re_cloud(i,k,j) = MAX(re_qc_min, MIN(re_qc1d(k), re_qc_max))
                                re_ice(i,k,j)   = MAX(re_qi_min, MIN(re_qi1d(k), re_qi_max))
                                re_snow(i,k,j)  = MAX(re_qs_min, MIN(re_qs1d(k), re_qs_max))
                            enddo
                        ENDIF
                    ENDIF last_step_only

                enddo i_loop
            enddo j_loop

            ! DEBUG - GT
            !      write(*,'(a,7(a,e13.6,1x,a,i3,a,i3,a,i3,a,1x))') 'MP-GT:', &
            !         'qc: ', qc_max, '(', imax_qc, ',', jmax_qc, ',', kmax_qc, ')', &
            !         'qr: ', qr_max, '(', imax_qr, ',', jmax_qr, ',', kmax_qr, ')', &
            !         'qi: ', qi_max, '(', imax_qi, ',', jmax_qi, ',', kmax_qi, ')', &
            !         'qs: ', qs_max, '(', imax_qs, ',', jmax_qs, ',', kmax_qs, ')', &
            !         'qg: ', qg_max, '(', imax_qg, ',', jmax_qg, ',', kmax_qg, ')', &
            !         'ni: ', ni_max, '(', imax_ni, ',', jmax_ni, ',', kmax_ni, ')', &
            !         'nr: ', nr_max, '(', imax_nr, ',', jmax_nr, ',', kmax_nr, ')'
            ! END DEBUG - GT
        enddo ! end of nt loop

        do j = j_start, j_end
            do k = kts, kte
                do i = i_start, i_end
                    pfils(i,k,j) = pfils(i,k,j)/dt_in
                    pflls(i,k,j) = pflls(i,k,j)/dt_in
                enddo
            enddo
        enddo

        ! These are always allocated
        !deallocate (vtsk1)
        !deallocate (txri1)
        !deallocate (txrc1)
        deallocate_extended_diagnostics: if (ext_diag) then
            deallocate (prw_vcdc1)
            deallocate (prw_vcde1)
            deallocate (tpri_inu1)
            deallocate (tpri_ide1_d)
            deallocate (tpri_ide1_s)
            deallocate (tprs_ide1)
            deallocate (tprs_sde1_d)
            deallocate (tprs_sde1_s)
            deallocate (tprg_gde1_d)
            deallocate (tprg_gde1_s)
            deallocate (tpri_iha1)
            deallocate (tpri_wfz1)
            deallocate (tpri_rfz1)
            deallocate (tprg_rfz1)
            deallocate (tprs_scw1)
            deallocate (tprg_scw1)
            deallocate (tprg_rcs1)
            deallocate (tprs_rcs1)
            deallocate (tprr_rci1)
            deallocate (tprg_rcg1)
            deallocate (tprw_vcd1_c)
            deallocate (tprw_vcd1_e)
            deallocate (tprr_sml1)
            deallocate (tprr_gml1)
            deallocate (tprr_rcg1)
            deallocate (tprr_rcs1)
            deallocate (tprv_rev1)
            deallocate (tten1)
            deallocate (qvten1)
            deallocate (qrten1)
            deallocate (qsten1)
            deallocate (qgten1)
            deallocate (qiten1)
            deallocate (niten1)
            deallocate (nrten1)
            deallocate (ncten1)
            deallocate (qcten1)
        end if deallocate_extended_diagnostics

    END SUBROUTINE mp_gt_driver
    !> @}

    !>\ingroup aathompson
    SUBROUTINE thompson_finalize()

        IMPLICIT NONE

        if (ALLOCATED(tcg_racg)) DEALLOCATE(tcg_racg)
        if (ALLOCATED(tmr_racg)) DEALLOCATE(tmr_racg)
        if (ALLOCATED(tcr_gacr)) DEALLOCATE(tcr_gacr)
        if (ALLOCATED(tnr_racg)) DEALLOCATE(tnr_racg)
        if (ALLOCATED(tnr_gacr)) DEALLOCATE(tnr_gacr)

        if (ALLOCATED(tcs_racs1)) DEALLOCATE(tcs_racs1)
        if (ALLOCATED(tmr_racs1)) DEALLOCATE(tmr_racs1)
        if (ALLOCATED(tcs_racs2)) DEALLOCATE(tcs_racs2)
        if (ALLOCATED(tmr_racs2)) DEALLOCATE(tmr_racs2)
        if (ALLOCATED(tcr_sacr1)) DEALLOCATE(tcr_sacr1)
        if (ALLOCATED(tms_sacr1)) DEALLOCATE(tms_sacr1)
        if (ALLOCATED(tcr_sacr2)) DEALLOCATE(tcr_sacr2)
        if (ALLOCATED(tms_sacr2)) DEALLOCATE(tms_sacr2)
        if (ALLOCATED(tnr_racs1)) DEALLOCATE(tnr_racs1)
        if (ALLOCATED(tnr_racs2)) DEALLOCATE(tnr_racs2)
        if (ALLOCATED(tnr_sacr1)) DEALLOCATE(tnr_sacr1)
        if (ALLOCATED(tnr_sacr2)) DEALLOCATE(tnr_sacr2)

        if (ALLOCATED(tpi_qcfz)) DEALLOCATE(tpi_qcfz)
        if (ALLOCATED(tni_qcfz)) DEALLOCATE(tni_qcfz)

        if (ALLOCATED(tpi_qrfz)) DEALLOCATE(tpi_qrfz)
        if (ALLOCATED(tpg_qrfz)) DEALLOCATE(tpg_qrfz)
        if (ALLOCATED(tni_qrfz)) DEALLOCATE(tni_qrfz)
        if (ALLOCATED(tnr_qrfz)) DEALLOCATE(tnr_qrfz)

        if (ALLOCATED(tps_iaus)) DEALLOCATE(tps_iaus)
        if (ALLOCATED(tni_iaus)) DEALLOCATE(tni_iaus)
        if (ALLOCATED(tpi_ide))  DEALLOCATE(tpi_ide)

        if (ALLOCATED(t_Efrw)) DEALLOCATE(t_Efrw)
        if (ALLOCATED(t_Efsw)) DEALLOCATE(t_Efsw)

        if (ALLOCATED(tnr_rev)) DEALLOCATE(tnr_rev)
        if (ALLOCATED(tpc_wev)) DEALLOCATE(tpc_wev)
        if (ALLOCATED(tnc_wev)) DEALLOCATE(tnc_wev)

        if (ALLOCATED(tnccn_act)) DEALLOCATE(tnccn_act)

    END SUBROUTINE thompson_finalize

end module module_mp_thompson
 !+---+-----------------------------------------------------------------+
 !ctrlL
 !+---+-----------------------------------------------------------------+
 !+---+-----------------------------------------------------------------+

