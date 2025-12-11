! 3D TEMPO Driver for MPAS
!=================================================================================================================
module module_mp_tempo_driver

    ! use mpas_kind_types, only: wp => RKIND, sp => R4KIND, dp => R8KIND
    use module_mp_tempo_params
    ! use module_mp_tempo_utils
    use module_mp_tempo_diags, only : calc_effectRad, calc_refl10cm, hail_size_diagnostics
    use module_mp_tempo_main, only : mp_tempo_main
    use module_mp_tempo_ml, only : predict_number_sub
    ! use module_mp_tempo_init, only : tempo_init_cfgs

    !use mpas_atmphys_utilities, only : physics_message, physics_error_fatal
    !use mpas_io_units, only : mpas_new_unit, mpas_release_unit
    !use mp_radar

    implicit none

contains


    !=================================================================================================================
    ! This is a wrapper routine designed to transfer values from 3D to 1D.
    ! Required microphysics variables are qv, qc, qr, nr, qi, ni, qs, qg
    ! Optional microphysics variables are aerosol aware (nc, nwfa, nifa, nwfa2d, nifa2d), and hail aware (ng, qg)

    subroutine tempo_3d_to_1d_driver(qv, qc, qr, qi, qs, qg, qb, ni, nr, nc, ng, &
        nwfa, nifa, nwfa2d, nifa2d, th, pii, p, w, dz, dt_in, itimestep, &
        rainnc, rainncv, snownc, snowncv, graupelnc, graupelncv, sr, frainnc, &
        refl_10cm, diagflag, do_radar_ref, re_cloud, re_ice, re_snow, qcbl, cldfrac, &
        has_reqc, has_reqi, has_reqs, ntc, muc, rainprod, evapprod, &
        max_hail_diameter_column, max_hail_diameter_sfc, &
        ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)

        ! Subroutine (3D) arguments
        integer, intent(in) :: ids,ide, jds,jde, kds,kde, ims,ime, jms,jme, kms,kme, its,ite, jts,jte, kts,kte
        real, dimension(ims:ime, kms:kme, jms:jme), intent(inout) :: qv, qc, qr, qi, qs, qg, ni, nr, th
        real, dimension(ims:ime, kms:kme, jms:jme), intent(inout) :: re_cloud, re_ice, re_snow
        integer, intent(in) :: has_reqc, has_reqi, has_reqs
        real, dimension(ims:ime, kms:kme, jms:jme), intent(in) :: pii, p, w, dz
        real, dimension(ims:ime, jms:jme), intent(inout) :: rainnc, rainncv, sr
        real, optional, dimension(ims:ime,jms:jme), intent(inout) :: frainnc, max_hail_diameter_column, max_hail_diameter_sfc
        real, dimension(ims:ime, kms:kme, jms:jme), intent(inout) :: rainprod, evapprod
        real, dimension(ims:ime, jms:jme), intent(in), optional :: ntc, muc
        real, dimension(ims:ime, kms:kme, jms:jme), intent(inout), optional :: nc, nwfa, nifa, qb, ng
        real, dimension(ims:ime, jms:jme), intent(in), optional :: nwfa2d, nifa2d
        real, dimension(ims:ime, kms:kme, jms:jme), intent(inout), optional :: refl_10cm
        real, dimension(ims:ime, kms:kme, jms:jme), intent(in), optional :: qcbl, cldfrac
        real, dimension(ims:ime, jms:jme), intent(inout), optional :: snownc, snowncv, graupelnc, graupelncv
        real, intent(in) :: dt_in
        integer, intent(in) :: itimestep

        ! Local (1d) variables
        real, dimension(kts:kte) :: qv1d, qc1d, qi1d, qr1d, qs1d, qg1d, qb1d, ni1d, nr1d, nc1d, ng1d, &
            nwfa1d, nifa1d, t1d, p1d, w1d, dz1d, rho, dbz, qcbl1d, cldfrac1d, qg_max_diam1d
        real, dimension(kts:kte) :: re_qc1d, re_qi1d, re_qs1d
        real, dimension(kts:kte):: rainprod1d, evapprod1d
        double precision, dimension(kts:kte) :: ncbl1d
        real, dimension(its:ite, jts:jte) :: pcp_ra, pcp_sn, pcp_gr, pcp_ic, frain
        real :: dt, pptrain, pptsnow, pptgraul, pptice
        real :: qc_max, qr_max, qs_max, qi_max, qg_max, ni_max, nr_max
        real :: nwfa1
        real :: ygra1, zans1
        real :: graupel_vol
        real :: tmprc, tmpnc, xDc
        integer :: nu_c
        logical, dimension(kts:kte) :: sgs_clouds
        double precision :: lamg, lam_exp, lamr, n0_min, n0_exp, lamc
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

                    sgs_clouds(k) = .false.
                    if (present(qcbl) .and. present(cldfrac)) then
                       qcbl1d(k) = qcbl(i,k,j)
                       cldfrac1d(k) = cldfrac(i,k,j)
                       ncbl1d(k) = 0.
                    endif

                    ! nwfa, nifa, and nc are optional aerosol-aware variables
                    if (present(nwfa)) then
                        if (present(nwfa2d)) then
                           if (k == kts) then
                              nwfa(i,k,j) = nwfa(i,k,j) + nwfa2d(i,j) * dt
                           endif
                        endif
                        nwfa(i,k,j) = max(nwfa_default, min(aero_max, nwfa(i,k,j)))
                        nwfa1d(k) = nwfa(i,k,j)
                    else
                        nwfa1d(k) = nwfa_default / rho(k)
                        !configs%aerosol_aware = .false.
                    endif

                    if (present(nifa)) then
                        nifa1d(k) = nifa(i,k,j)
                    else
                        nifa1d(k) = nifa_default / rho(k)
                        !configs%aerosol_aware = .false.
                    endif

                    if (present(nc)) then
                        nc1d(k) = nc(i,k,j)
                    else
                        nc1d(k) = Nt_c / rho(k)
                        !configs%aerosol_aware = .false.
                    endif
                enddo

                ! ng and qb are optional hail-aware variables
                if ((present(ng)) .and. (present(qb))) then
                    !configs%hail_aware = .true.
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
                            ng1d(k) = 0.
                            qb1d(k) = 0.
                        endif
                    enddo
                endif
                
                !if (itimestep == 1) then
                !   call physics_message('--- tempo_3d_to_1d_driver() configuration...')
                !   write(message, '(L1)') configs%hail_aware
                !   call physics_message('       hail_aware_flag = ' // trim(message))
                !   write(message, '(L1)') configs%aerosol_aware
                !   call physics_message('       aerosol_aware_flag = ' // trim(message))
                !   call physics_message('calling mp_tempo_main() at itimestep = 1')
                !endif

                !=================================================================================================================
                ! Main call to the 1D microphysics
                call mp_tempo_main(qv1d=qv1d, qc1d=qc1d, qi1d=qi1d, qr1d=qr1d, qs1d=qs1d, qg1d=qg1d, qb1d=qb1d, &
                           ni1d=ni1d, nr1d=nr1d, nc1d=nc1d, ng1d=ng1d, nwfa1d=nwfa1d, nifa1d=nifa1d, t1d=t1d, p1d=p1d, &
                           w1d=w1d, dzq=dz1d, pptrain=pptrain, pptsnow=pptsnow, pptgraul=pptgraul, pptice=pptice, &
                           kts=kts, kte=kte, dt=dt, ii=i, jj=j) !, configs=configs)

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
                if (present(frainnc)) then
                   frain(i,j) = 0.
                   if(t1d(1) <= 273.) then
                      frain(i,j) = pcp_ra(i,j)
                   endif
                   frainnc(i,j) = frainnc(i,j) + frain(i,j)
                endif

                sr(i,j) = (pptsnow + pptgraul + pptice) / (rainncv(i,j) + R1)

                ! ng and qb are optional hail-aware variables
                if ((present(ng)) .and. (present(qb))) then
                    do k = kts, kte
                        ng(i,k,j) = ng1d(k)
                        qb(i,k,j) = qb1d(k)
                    enddo
                else
                    do k = kte, kts, -1
                        ! This is the one-moment graupel formulation
                        if (qg1d(k) > R1) then
                            rho(k) = RoverRv * p1d(k) / (R * t1d(k) * (qv1d(k)+RoverRv))
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
                            ng1d(k) = 0.
                            qb1d(k) = 0.
                        endif
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

                    if (present(qcbl) .and. present(cldfrac)) then
                       if ((qc1d(k) <= R1) .and. (qcbl1d(k) > 1.e-9) .and. (cldfrac1d(k) > 0.)) then
                          qc1d(k) = qc1d(k) + qcbl1d(k)/cldfrac1d(k) ! Uses in-cloud PBL mass
                          sgs_clouds(k) = .true.
                       else
                          sgs_clouds(k) = .false.
                       endif
                    else
                       sgs_clouds(k) = .false.
                    endif
                 enddo

                 if (any(sgs_clouds)) then
                    ! return array of ncbl1d
                    call predict_number_sub(kts, kte, qc1d, qr1d, qi1d, qs1d, p1d, t1d, w1d, &
                         ncbl1d, predict_nc=.true.)
                    do k = kts, kte
                       if (sgs_clouds(k)) then
                          nc1d(k) = nc1d(k) + real(ncbl1d(k))
                          rho(k) = RoverRv * p1d(k) / (R * t1d(k) * (qv1d(k)+RoverRv))
                          tmprc = qc1d(k)*rho(k)
                          tmpnc = max(2., min(nc1d(k)*rho(k), nt_c_max))
                          if (tmpnc.gt.10000.e6) then
                             nu_c = 2
                          elseif (tmpnc.lt.100.) then
                             nu_c = 15
                          else
                             nu_c = nint(nu_c_scale/tmpnc) + 2
                             nu_c = max(2, min(nu_c, 15))
                          endif
                          lamc = (tmpnc*am_r*ccg(2,nu_c)*ocg1(nu_c)/tmprc)**obmr
                          xDc = (bm_r + nu_c + 1.) / lamc
                          if (xDc .lt. D0c) then
                             lamc = cce(2,nu_c)/D0c
                          elseif (xDc.gt. D0r*2.) then
                             lamc = cce(2,nu_c)/(D0r*2.)
                          endif
                          tmpnc = min(real(nt_c_max, kind=dp), ccg(1,nu_c)*ocg2(nu_c)*tmprc / am_r*lamc**bm_r)
                          nc1d(k) = tmpnc/rho(k) ! Update nc1d for calc_effectRad
                       endif
                    enddo
                 endif
                    ! if (present(qcbl) .and. present(cldfrac)) then
                    !    if ((qc1d(k) <= R1) .and. (qcbl(i,k,j) > 1.e-9) .and. (cldfrac(i,k,j) > 0.)) then
                    !       qc1d(k) = qc1d(k) + qcbl(i,k,j)/cldfrac(i,k,j) ! Uses in-cloud PBL mass
                    !       ! ML prediction of number concentration (Don't add in qibl for now)
                    !       ! COULD BE DLOW FROM CALLING CALLING SUBROUTINE FOR EACH i,k,j
                    !       ncbl(k) = predict_number(qc1d(k), qr1d(k), qi1d(k), qs1d(k), &
                    !            p1d(k), t1d(k), w1d(k), predict_nc=.true.)
                    !       nc1d(k) = nc1d(k) + ncbl(k)
                    !       rho(k) = RoverRv * p1d(k) / (R * t1d(k) * (qv1d(k)+RoverRv))
                    !       tmprc = qc1d(k)*rho(k)
                    !       tmpnc = max(2., min(nc1d(k)*rho(k), nt_c_max))
                    !       if (tmpnc.gt.10000.e6) then
                    !          nu_c = 2
                    !       elseif (tmpnc.lt.100.) then
                    !          nu_c = 15
                    !       else
                    !          nu_c = nint(nu_c_scale/tmpnc) + 2
                    !          nu_c = max(2, min(nu_c, 15))
                    !       endif
                    !       lamc = (tmpnc*am_r*ccg(2,nu_c)*ocg1(nu_c)/tmprc)**obmr
                    !       xDc = (bm_r + nu_c + 1.) / lamc
                    !       if (xDc .lt. D0c) then
                    !          lamc = cce(2,nu_c)/D0c
                    !       elseif (xDc.gt. D0r*2.) then
                    !          lamc = cce(2,nu_c)/(D0r*2.)
                    !       endif
                    !       tmpnc = min(real(nt_c_max, kind=dp), ccg(1,nu_c)*ocg2(nu_c)*tmprc / am_r*lamc**bm_r)
                    !       nc1d(k) = tmpnc/rho(k) ! Update nc1d for calc_effectRad
                    !    endif
                    ! endif
!!! K LOOP                enddo

                !=================================================================================================================
                ! Reflectivity
                call calc_refl10cm (qv1d=qv1d, qc1d=qc1d, qr1d=qr1d, nr1d=nr1d, qs1d=qs1d, qg1d=qg1d, ng1d=ng1d, qb1d=qb1d, &
                    t1d=t1d, p1d=p1d, dBZ=dBZ, kts=kts, kte=kte, ii=i, jj=j) !, configs=configs)
                do k = kts, kte
                    refl_10cm(i,k,j) = max(-35.0_wp, dBZ(k))
                enddo

                if ((present(max_hail_diameter_sfc)) .and. (present(max_hail_diameter_column))) then
                   ! Maximium hail size
                   call hail_size_diagnostics(kts=kts, kte=kte, qg1d=qg1d, ng1d=ng1d, qb1d=qb1d, t1d=t1d, p1d=p1d, qv1d=qv1d, &
                        qg_max_diam1d=qg_max_diam1d) !, configs=configs)

                   max_hail_diameter_sfc(i,j) = max(0.0_wp, qg_max_diam1d(kts))
                   max_hail_diameter_column(i,j) = max(0.0_wp, maxval(qg_max_diam1d))
                endif

                ! Cloud, ice, and snow effective radius
                if (has_reqc /= 0 .and. has_reqi /= 0 .and. has_reqs /= 0) then
                    do k = kts, kte
                        re_qc1d(k) = 2.49e-6
                        re_qi1d(k) = 4.99e-6
                        re_qs1d(k) = 9.99e-6
                    enddo
                    call calc_effectRad (t1d=t1d, p1d=p1d, qv1d=qv1d, qc1d=qc1d, nc1d=nc1d, qi1d=qi1d, &
                         ni1d=ni1d, qs1d=qs1d, re_qc1d=re_qc1d, re_qi1d=re_qi1d, re_qs1d=re_qs1d, &
                         kts=kts, kte=kte) !, configs=configs)
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

end module module_mp_tempo_driver
