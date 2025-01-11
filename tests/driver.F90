! gfortran -c machine.F90 ../module_mp_tempo_params.F90 ../module_mp_tempo_utils.F90 ../drivers/standalone/module_mp_tempo.F90 driver.F90
! gfortran machine.o module_mp_tempo_params.o module_mp_tempo_utils.o module_mp_tempo.o driver.o
! ./a.out                                                                                                                                                       
!=================================================================================================================                                              

program driver
  use machine, only: wp => kind_phys, sp => kind_sngl_prec, dp => kind_dbl_prec
  use module_mp_tempo_params
  use module_mp_tempo, only : tempo_init, tempo_3d_to_1d_driver

  implicit none

  ! call build_lookup_tables                                                                                                                                    
  call run_profile()

contains
!=================================================================================================================                                              

  subroutine build_lookup_tables


  end subroutine build_lookup_tables
!=================================================================================================================                                              

subroutine run_profile

    logical :: l_mp_tables
    logical :: hail_aware_flag
    logical :: aerosol_aware_flag

    integer, parameter :: ids=1, ide=1, jds=1, jde=1, kds=1, kde=50
    integer, parameter :: ims=1, ime=1, jms=1, jme=1, kms=1, kme=50
    integer, parameter :: its=1, ite=1, jts=1, jte=1, kts=1, kte=50
    real(wp), dimension(its:ite,kts:kte,jts:jte) :: qv, qc, qr, qi, qs, qg, ni, nr, nc, nwfa, nifa
    real(wp), dimension(its:ite,kts:kte,jts:jte) :: th, pii, p, w, dz, refl_10cm
    real(wp), dimension(its:ite,kts:kte,jts:jte) :: re_cloud, re_ice, re_snow
    real(wp), dimension(its:ite,jts:jte) :: nwfa2d, nifa2d, rainnc, rainncv, sr, ntc, muc
    real, parameter :: dt_in = 1
    integer, parameter :: has_reqc = 1
    integer, parameter :: has_reqi = 1
    integer, parameter :: has_reqs = 1

    integer :: k, t, i_start, j_start, i_end, j_end, i, j, itimestep

    l_mp_tables = .true.
    hail_aware_flag = .false.
    aerosol_aware_flag = .true.
    write(*,*) '--- calling tempo_init()'
    call tempo_init(l_mp_tables, hail_aware_flag, aerosol_aware_flag)

    i_start = its
    j_start = jts
    i_end = ite
    j_end = jte
    j_loop:  do j = j_start, j_end
      i_loop:  do i = i_start, i_end
        do k = kts, kte
            qv(i,k,j) = 5.0e-3
            qc(i,k,j) = 1.0e-3
            qr(i,k,j) = 1.0e-3
            qi(i,k,j) = 1.0e-3
            qs(i,k,j) = 1.0e-3
            qg(i,k,j) = 1.0e-3
            ni(i,k,j) = 1.0e3
            nr(i,k,j) = 1.0e3
            nc(i,k,j) = 1.0e6
            nwfa(i,k,j) = 0.0
            nifa(i,k,j) = 0.0
            th(i,k,j) = 280.0
            p(i,k,j) = 80000.
            pii(i,k,j) = (p(i,k,j)/100000.0)**0.286
            w(i,k,j) = 0.0
            dz(i,k,j) = 50.0
            refl_10cm(i,j,j) = -35.0
            re_cloud(i,k,j) = 2.5
            re_ice(i,k,j) = 2.5
            re_snow(i,k,j) = 2.5
        enddo
        nwfa2d(i,j) = 0.0
        nwfa2d(i,j) = 0.0
        rainnc(i,j) = 0.0
        rainncv(i,j) = 0.0
        sr(i,j) = 0.0
        ntc(i,j) = 300.e6
        muc(i,j) = 4.0
      enddo i_loop
    enddo j_loop

    ! Time integration                                                                                                                                          
    do t = 1, 1800
      itimestep = t
      write(*,*) 'Timestep: ', t

      call tempo_3d_to_1d_driver(qv=qv, qc=qc, qr=qr, qi=qi, qs=qs, qg=qg, ni=ni, nr=nr, nc=nc, &
          nwfa=nwfa, nifa=nifa, nwfa2d=nwfa2d, nifa2d=nwfa2d, th=th, pii=pii, p=p, w=w, dz=dz, dt_in=dt_in, &
          itimestep=itimestep, rainnc=rainnc, rainncv=rainncv, sr=sr, refl_10cm=refl_10cm, &
          re_cloud=re_cloud, re_ice=re_ice, re_snow=re_snow, has_reqc=has_reqc, has_reqi=has_reqi, has_reqs=has_reqs, &
          ntc=ntc, muc=muc, &
          ids=ids, ide=ide, jds=jds, jde=jde, kds=kds, kde=kde, &
          ims=ims, ime=ime, jms=jms, jme=jme, kms=kms, kme=kme, its=its, ite=ite, jts=jts, jte=jte, kts=kts, kte=kte)
    enddo

  end subroutine run_profile

!=================================================================================================================                                              
  
end program driver
