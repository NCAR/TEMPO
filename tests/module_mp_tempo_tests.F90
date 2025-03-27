! Module for TEMPO Microphysics tests
!=================================================================================================================
module module_mp_tempo_tests

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

    implicit none
  
    integer, parameter :: ids=1, ide=1, jds=1, jde=1, kds=1
    integer, parameter :: ims=1, ime=1, jms=1, jme=1, kms=1
    integer, parameter :: its=1, ite=1, jts=1, jte=1, kts=1
    real, parameter :: dt_in = 1
    integer, parameter :: integration_length_sec = 3600
    integer :: kde, kme, kte
    logical :: l_mp_tables, hail_aware_flag, aerosol_aware_flag
    real(wp), allocatable :: qv(:,:,:), qc(:,:,:), qr(:,:,:), qi(:,:,:), qs(:,:,:), qg(:,:,:), &
         ni(:,:,:), nr(:,:,:), nc(:,:,:), nwfa(:,:,:), nifa(:,:,:), ng(:,:,:), qb(:,:,:)
    real(wp), allocatable :: th(:,:,:), pii(:,:,:), p(:,:,:), w(:,:,:), dz(:,:,:), refl_10cm(:,:,:)
    real(wp), allocatable :: re_cloud(:,:,:), re_ice(:,:,:), re_snow(:,:,:)
    real(wp), dimension(its:ite,jts:jte) :: nwfa2d, nifa2d, rainnc, rainncv, sr
    
  contains

    !=================================================================================================================    

    subroutine init_tempo_flags_for_test_all_true()
      
      l_mp_tables = .true.
      aerosol_aware_flag = .true.
      hail_aware_flag = .true.
      
    end subroutine init_tempo_flags_for_test_all_true

    !=================================================================================================================
    
    subroutine init_mpas_59lev_convective_data_for_test(nlev)

      integer, intent(in) :: nlev
      real(wp), dimension(nlev) :: klevs, qv_in, qc_in, qr_in, qi_in, qs_in, qg_in, ni_in, nr_in, nc_in, &
           nwfa_in, nifa_in, theta_in, ng_in, volg_in, pressure_in, w_in, dz_in
      integer :: k

      kde = nlev
      kme = nlev
      kte = nlev
      
      ! opening data file
      open (2, file = './data/mpas_59lev_test.txt', status = 'old')
      read(2,*) ! header
      
      do k = 1, nlev
         read(2,*) klevs(k), qv_in(k), qc_in(k), qr_in(k), qi_in(k), qs_in(k), &
              qg_in(k), ni_in(k), nr_in(k), nc_in(k), nwfa_in(k), nifa_in(k), &
              theta_in(k), ng_in(k), volg_in(k), pressure_in(k), w_in(k), dz_in(k)
      end do
      close(2)

      if(.not. allocated(qv)) allocate(qv(ite-its+1, kte-kts+1, jte-jts+1))
      if(.not. allocated(qc)) allocate(qc(ite-its+1, kte-kts+1, jte-jts+1))
      if(.not. allocated(qr)) allocate(qr(ite-its+1, kte-kts+1, jte-jts+1))
      if(.not. allocated(qi)) allocate(qi(ite-its+1, kte-kts+1, jte-jts+1))
      if(.not. allocated(qs)) allocate(qs(ite-its+1, kte-kts+1, jte-jts+1))
      if(.not. allocated(qg)) allocate(qg(ite-its+1, kte-kts+1, jte-jts+1))
      if(.not. allocated(ni)) allocate(ni(ite-its+1, kte-kts+1, jte-jts+1))
      if(.not. allocated(nr)) allocate(nr(ite-its+1, kte-kts+1, jte-jts+1))
      if(.not. allocated(nc)) allocate(nc(ite-its+1, kte-kts+1, jte-jts+1))
      if(.not. allocated(nwfa)) allocate(nwfa(ite-its+1, kte-kts+1, jte-jts+1))
      if(.not. allocated(nifa)) allocate(nifa(ite-its+1, kte-kts+1, jte-jts+1))
      if(.not. allocated(ng)) allocate(ng(ite-its+1, kte-kts+1, jte-jts+1))
      if(.not. allocated(qb)) allocate(qb(ite-its+1, kte-kts+1, jte-jts+1))
      if(.not. allocated(th)) allocate(th(ite-its+1, kte-kts+1, jte-jts+1))
      if(.not. allocated(p)) allocate(p(ite-its+1, kte-kts+1, jte-jts+1))
      if(.not. allocated(pii)) allocate(pii(ite-its+1, kte-kts+1, jte-jts+1))      
      if(.not. allocated(w)) allocate(w(ite-its+1, kte-kts+1, jte-jts+1))
      if(.not. allocated(dz)) allocate(dz(ite-its+1, kte-kts+1, jte-jts+1))
      if(.not. allocated(refl_10cm)) allocate(refl_10cm(ite-its+1, kte-kts+1, jte-jts+1))
      if(.not. allocated(re_cloud)) allocate(re_cloud(ite-its+1, kte-kts+1, jte-jts+1))
      if(.not. allocated(re_ice)) allocate(re_ice(ite-its+1, kte-kts+1, jte-jts+1))
      if(.not. allocated(re_snow)) allocate(re_snow(ite-its+1, kte-kts+1, jte-jts+1))                        
      
      qv = reshape(qv_in, (/ite-its+1, kte-kts+1, jte-jts+1/))
      qc = reshape(qc_in, (/ite-its+1, kte-kts+1, jte-jts+1/))
      qr = reshape(qr_in, (/ite-its+1, kte-kts+1, jte-jts+1/))
      qi = reshape(qi_in, (/ite-its+1, kte-kts+1, jte-jts+1/))
      qs = reshape(qs_in, (/ite-its+1, kte-kts+1, jte-jts+1/))
      qg = reshape(qg_in, (/ite-its+1, kte-kts+1, jte-jts+1/))
      ni = reshape(ni_in, (/ite-its+1, kte-kts+1, jte-jts+1/))
      nr = reshape(nr_in, (/ite-its+1, kte-kts+1, jte-jts+1/))
      nc = reshape(nc_in, (/ite-its+1, kte-kts+1, jte-jts+1/))
      nwfa = reshape(nwfa_in, (/ite-its+1, kte-kts+1, jte-jts+1/))
      nifa = reshape(nifa_in, (/ite-its+1, kte-kts+1, jte-jts+1/))
      ng = reshape(ng_in, (/ite-its+1, kte-kts+1, jte-jts+1/))
      qb = reshape(volg_in, (/ite-its+1, kte-kts+1, jte-jts+1/))                  
      th = reshape(theta_in, (/ite-its+1, kte-kts+1, jte-jts+1/))
      p = reshape(pressure_in, (/ite-its+1, kte-kts+1, jte-jts+1/))
      w = reshape(w_in, (/ite-its+1, kte-kts+1, jte-jts+1/))
      dz = reshape(dz_in, (/ite-its+1, kte-kts+1, jte-jts+1/))

      pii = (p/100000.0)**0.286
      refl_10cm = 0.
      nwfa2d = 0.
      nifa2d = 0.
      rainnc = 0.
      rainncv = 0.
      sr = 0.
      re_cloud = 0.
      re_ice = 0.
      re_snow = 0.
      
    end subroutine init_mpas_59lev_convective_data_for_test
    
    !=================================================================================================================

    subroutine mpas_test()
      
      use module_mp_tempo_params
      use module_mp_tempo, only : tempo_init, tempo_3d_to_1d_driver
     
      ! local variables
      integer :: t, itimestep, i, j, k
      
      ! Initialize input data for tests
      call init_mpas_59lev_convective_data_for_test(nlev=59)
      
      ! Initialize TEMPO flags for tests
      call init_tempo_flags_for_test_all_true()
      
      ! Initialize TEMPO
      write(*,*) '--- calling tempo_init()'
      call tempo_init(l_mp_tables, hail_aware_flag, aerosol_aware_flag)

      ! Time integration
      do t = 1, integration_length_sec
         itimestep = t

         if (t == 1) then
            write(*,*) '--- calling tempo_3d_to_1d_driver()'
         endif

         if (t == integration_length_sec) then
            write(*,*) 'Final timestep for TEMPO microphysics: ', t
         endif
         
         call tempo_3d_to_1d_driver(qv=qv, qc=qc, qr=qr, qi=qi, qs=qs, qg=qg, ni=ni, nr=nr, nc=nc, ng=ng, qb=qb, &
              nwfa=nwfa, nifa=nifa, nwfa2d=nwfa2d, nifa2d=nwfa2d, th=th, pii=pii, p=p, w=w, dz=dz, dt_in=dt_in, &
              itimestep=itimestep, rainnc=rainnc, rainncv=rainncv, sr=sr, &
              refl_10cm=refl_10cm, re_cloud=re_cloud, re_ice=re_ice, re_snow=re_snow, &
              has_reqc=0, has_reqi=0, has_reqs=0, &
              ids=ids, ide=ide, jds=jds, jde=jde, kds=kds, kde=kde, &
              ims=ims, ime=ime, jms=jms, jme=jme, kms=kms, kme=kme, &
              its=its, ite=ite, jts=jts, jte=jte, kts=kts, kte=kte)
      enddo

      ! output data to file 
      open(1, file = 'mpas_59lev_test_results.txt')
      write(1,*) 'klev qv qc qr qi qs qg ni nr nc nwfa nifa theta ng volg pressure w dz'
      do j = jts, jte
         do i = its, ite
            ! k-loop
            do k = kts, kte
               write(1,*) k, qv(i,k,j), qc(i,k,j), qr(i,k,j), qi(i,k,j), qs(i,k,j), qg(i,k,j), ni(i,k,j), nr(i,k,j), &
                    nc(i,k,j), nwfa(i,k,j), nifa(i,k,j), th(i,k,j), ng(i,k,j), qb(i,k,j), p(i,k,j), w(i,k,j), dz(i,k,j)
            enddo
         enddo
      enddo
      
      close(1)

    end subroutine mpas_test
    !=================================================================================================================
    
  end module module_mp_tempo_tests
