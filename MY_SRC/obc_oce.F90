MODULE obc_oce
   !!==============================================================================
   !!                       ***  MODULE obc_oce   ***
   !! Open Boundary Cond. :   define related variables
   !!==============================================================================
   !! history :  OPA  ! 1991-01 (CLIPPER)  Original code 
   !!   NEMO     1.0  ! 2002-02 (C. Talandier)  modules, F90
   !!----------------------------------------------------------------------
#if defined key_obc
   !!----------------------------------------------------------------------
   !!   'key_obc' :                                 Open Boundary Condition
   !!----------------------------------------------------------------------
   USE par_oce         ! ocean parameters
   USE obc_par         ! open boundary condition parameters
   USE par_trc         ! TOP parameters

   IMPLICIT NONE
   PUBLIC
   
   PUBLIC   obc_oce_alloc   ! called by obcini.F90 module
   PUBLIC   obc_oce_trc_alloc   ! called by obcini.F90 module


   !!----------------------------------------------------------------------
   !! open boundary variables
   !!----------------------------------------------------------------------
   !
   !                                            !!* Namelist namobc: open boundary condition *
   INTEGER           ::   nn_nbobc    = 1        !: number of open boundaries ( 1=< nbobc =< 4 ) 
   INTEGER           ::   nn_obcdta   = 0        !:  = 0 use the initial state as obc data
   !                                             !   = 1 read obc data in obcxxx.dta files
#if defined key_obc_trc
!  INTEGER ::
   INTEGER, DIMENSION(jptra) ::              & !
      nn_obcdta_trc = 0           !:  = 0 use the initial state as obc data
      !                           !   = 1 read obc data in obcxxx.dta files
   CHARACTER(len=20) ::    & !
      cn_obcdta_trc = 'annual'    !: set to annual  if obc datafile hold 1 year of data
      !                           !  set to monthly if obc datafile hold 1 month of data
!   LOGICAL ::              & !
   LOGICAL, DIMENSION(jptra) ::    & !
      ln_obc_ini = .true.        ! use the initial state as a climatology
   CHARACTER(len=12), PUBLIC, DIMENSION(jptra) ::   nam_obc_trc     !: tracer name for OBC
#endif
   CHARACTER(len=20) ::   cn_obcdta   = 'annual' !: set to annual  if obc datafile hold 1 year of data
   !                                             !  set to monthly if obc datafile hold 1 month of data
   LOGICAL           ::   ln_obc_clim = .true.   !:  obc data files are climatological
   LOGICAL           ::   ln_obc_fla  = .false.  !:  Flather open boundary condition not used
   LOGICAL           ::   ln_vol_cst  = .true.   !:  Conservation of the whole volume
#if defined key_obc_trc
   LOGICAL, PUBLIC, DIMENSION(jptra) ::   ln_obc_clim_trc = .true.  !:  obc data files are climatological
!  LOGICAL   ::   ln_obc_clim_trc = .true.
   LOGICAL                           ::   ln_obc_fla_trc  = .false. !:  Flather open boundary condition not used
   LOGICAL                           ::   ln_vol_cst_trc  = .true.  !:  Conservation of the whole volume
#endif
   REAL(wp)          ::   rn_dpein    =  1.      !: damping time scale for inflow at East open boundary
   REAL(wp)          ::   rn_dpwin    =  1.      !:    "                      "   at West open boundary
   REAL(wp)          ::   rn_dpsin    =  1.      !:    "                      "   at South open boundary
   REAL(wp)          ::   rn_dpnin    =  1.      !:    "                      "   at North open boundary
   REAL(wp)          ::   rn_dpeob    = 15.      !: damping time scale for the climatology at East open boundary
   REAL(wp)          ::   rn_dpwob    = 15.      !:    "                           "       at West open boundary
   REAL(wp)          ::   rn_dpsob    = 15.      !:    "                           "       at South open boundary
   REAL(wp)          ::   rn_dpnob    = 15.      !:    "                           "       at North open boundary
   REAL(wp)          ::   rn_volemp   =  1.      !: = 0 the total volume will have the variability of the 
   !                                             !      surface Flux E-P else (volemp = 1) the volume will be constant
   !                                             !  = 1 the volume will be constant during all the integration.
#if defined key_obc_trc
   REAL(wp) ::             & !!: open boundary namelist (namobc)
      rn_dpein_trc =  1.  ,      &  !: damping time scale for inflow at East open boundary
      rn_dpwin_trc =  1.  ,      &  !:    "                      "   at West open boundary
      rn_dpsin_trc =  1.  ,      &  !:    "                      "   at South open boundary
      rn_dpnin_trc =  1.  ,      &  !:    "                      "   at North open boundary
      rn_dpeob_trc = 15.  ,      &  !: damping time scale for the climatology at East open boundary
      rn_dpwob_trc = 15.  ,      &  !:    "                           "       at West open boundary
      rn_dpsob_trc = 15.  ,      &  !:    "                           "       at South open boundary
      rn_dpnob_trc = 15.            !:    "                           "       at North open boundary
#endif
   !                                  !!! OLD non-DOCTOR name of namelist variables
   INTEGER  ::   nbobc                 !: number of open boundaries ( 1=< nbobc =< 4 ) 
   INTEGER  ::   nobc_dta              !:  = 0 use the initial state as obc data
   REAL(wp) ::   rdpein                !: damping time scale for inflow at East open boundary
   REAL(wp) ::   rdpwin                !:    "                      "   at West open boundary
   REAL(wp) ::   rdpsin                !:    "                      "   at South open boundary
   REAL(wp) ::   rdpnin                !:    "                      "   at North open boundary
   REAL(wp) ::   rdpeob                !: damping time scale for the climatology at East open boundary
   REAL(wp) ::   rdpwob                !:    "                           "       at West open boundary
   REAL(wp) ::   rdpsob                !:    "                           "       at South open boundary
   REAL(wp) ::   rdpnob                !:    "                           "       at North open boundary
   REAL(wp) ::   volemp                !: = 0 the total volume will have the variability of the 
#if defined key_obc_trc
!   INTEGER  ::   nobc_dta_trc              !:  = 0 use the initial state as obc data
   INTEGER, DIMENSION(jptra) :: nobc_dta_trc
   REAL(wp) ::   rdpein_trc                !: damping time scale for inflow at East open boundary
   REAL(wp) ::   rdpwin_trc                !:    "                      "   at West open boundary
   REAL(wp) ::   rdpsin_trc                !:    "                      "   at South open boundary
   REAL(wp) ::   rdpnin_trc                !:    "                      "   at North open boundary
   REAL(wp) ::   rdpeob_trc                !: damping time scale for the climatology at East open boundary
   REAL(wp) ::   rdpwob_trc                !:    "                           "       at West open boundary
   REAL(wp) ::   rdpsob_trc                !:    "                           "       at South open boundary
   REAL(wp) ::   rdpnob_trc                !:    "                           "       at North open boundary
   REAL(wp) ::   volemp_trc                !: = 0 the total volume will have the variability of the 
   !                                   !      surface Flux E-P else (volemp = 1) the volume will be constant
   !                                   !  = 1 the volume will be constant during all the integration.
   CHARACTER ( len=20 ) :: cffile_trc
#endif
   CHARACTER(len=20) :: cffile


   !!General variables for open boundaries:
   !!--------------------------------------
   LOGICAL ::   lfbceast, lfbcwest      !: logical flag for a fixed East and West open boundaries       
   LOGICAL ::   lfbcnorth, lfbcsouth    !: logical flag for a fixed North and South open boundaries       
   !                                    !  These logical flags are set to 'true' if damping time 
   !                                    !  scale are set to 0 in the namelist, for both inflow and outflow).
#if defined key_obc_trc
   LOGICAL ::   lfbceast_trc, lfbcwest_trc      !: logical flag for a fixed East and West open boundaries       
   LOGICAL ::   lfbcnorth_trc, lfbcsouth_trc    !: logical flag for a fixed North and South open boundaries       
   !                                    !  These logical flags are set to 'true' if damping time 
   !                                    !  scale are set to 0 in the namelist, for both inflow and outflow).
#endif

 ! PERIANT CONFIG
   REAL(wp), PUBLIC ::   obcsurf_indian   !: lateral surface of open boundaries Indian
   REAL(wp), PUBLIC ::   obcsurf_pacif    !: lateral surface of open boundaries Pacific
   REAL(wp), PUBLIC ::   obcsurf_atl      !: lateral surface of open boundaries Atlantic

   REAL(wp), PUBLIC ::   obcsurftot       !: Total lateral surface of open boundaries
   
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) :: &  !:
      obctmsk,            &  !: mask array identical to tmask, execpt along OBC where it is set to 0
      !                      !  it used to calculate the cumulate flux E-P in the obcvol.F90 routine
      obcumask, obcvmask     !: u-, v- Force filtering mask for the open
      !                      !  boundary condition on grad D

   !!--------------------
   !! East open boundary:
   !!--------------------
   INTEGER ::   nie0  , nie1      !: do loop index in mpp case for jpieob
   INTEGER ::   nie0p1, nie1p1    !: do loop index in mpp case for jpieob+1
   INTEGER ::   nie0m1, nie1m1    !: do loop index in mpp case for jpieob-1
   INTEGER ::   nje0  , nje1      !: do loop index in mpp case for jpjed, jpjef
   INTEGER ::   nje0p1, nje1m1    !: do loop index in mpp case for jpjedp1,jpjefm1
   INTEGER ::   nje1m2, nje0m1    !: do loop index in mpp case for jpjefm1-1,jpjed

   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:) ::   &  !:
      sshfoe,           & !: now climatology of the east boundary sea surface height
      ubtfoe,vbtfoe       !: now climatology of the east boundary barotropic transport
     
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:) ::   &  !:
      ufoe, vfoe,       & !: now climatology of the east boundary velocities 
      tfoe, sfoe,       & !: now climatology of the east boundary temperature and salinity
      uclie               !: baroclinic componant of the zonal velocity after radiation 
      !                   ! in the obcdyn.F90 routine

#if defined key_obc_trc
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   &  !:
      trcfoe               !: now climatology of the east boundary passive tracers
#endif
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:) ::   sshfoe_b   !: east boundary ssh correction averaged over the barotropic loop
      !                                            !  (if Flather's algoritm applied at open boundary)

   !!-------------------------------
   !! Arrays for radiative East OBC: 
   !!-------------------------------
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:,:) ::   uebnd, vebnd      !: baroclinic u & v component of the velocity over 3 rows 
      !                                                    !  and 3 time step (now, before, and before before)
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:,:) ::   tebnd, sebnd      !: East boundary temperature and salinity over 2 rows 
      !                                                    !  and 2 time step (now and before)
   REAL(wp), ALLOCATABLE, SAVE,     DIMENSION(:,:) ::   u_cxebnd, v_cxebnd    !: Zonal component of the phase speed ratio computed with 
      !                                                    !  radiation of u and v velocity (respectively) at the 
      !                                                    !  east open boundary (u_cxebnd = cx rdt )
#if defined key_obc_trc
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:,:,:) ::   &  !:
      trcebnd                             !: East boundary temperature and salinity over 2 rows 
                                          ! and 2 time step (now and before)
#endif
   REAL(wp), ALLOCATABLE, SAVE,     DIMENSION(:,:) ::   uemsk, vemsk, temsk   !: 2D mask for the East OB

   ! Note that those arrays are optimized for mpp case 
   ! (hence the dimension jpj is the size of one processor subdomain)

   !!--------------------
   !! West open boundary
   !!--------------------
   INTEGER ::   niw0  , niw1       !: do loop index in mpp case for jpiwob
   INTEGER ::   niw0p1, niw1p1     !: do loop index in mpp case for jpiwob+1
   INTEGER ::   njw0  , njw1       !: do loop index in mpp case for jpjwd, jpjwf
   INTEGER ::   njw0p1, njw1m1     !: do loop index in mpp case for jpjwdp1,jpjwfm1
   INTEGER ::   njw1m2, njw0m1     !: do loop index in mpp case for jpjwfm2,jpjwd

   REAL(wp), ALLOCATABLE, SAVE,   DIMENSION(:) ::   &  !:
      sshfow,           & !: now climatology of the west boundary sea surface height
      ubtfow,vbtfow       !: now climatology of the west boundary barotropic transport

   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:) ::   &  !:
      ufow, vfow,       & !: now climatology of the west velocities 
      tfow, sfow,       & !: now climatology of the west temperature and salinity
      ucliw               !: baroclinic componant of the zonal velocity after the radiation 
      !                   !  in the obcdyn.F90 routine

#if defined key_obc_trc
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   &  !:
      trcfow               !: now climatology of the west boundary passive tracers
#endif
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:) ::   sshfow_b    !: west boundary ssh correction averaged over the barotropic loop
      !                                          !  (if Flather's algoritm applied at open boundary)

   !!-------------------------------
   !! Arrays for radiative West OBC
   !!-------------------------------
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:,:) ::   uwbnd, vwbnd     !: baroclinic u & v components of the velocity over 3 rows 
      !                                                   !  and 3 time step (now, before, and before before)
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:,:) ::   twbnd, swbnd     !: west boundary temperature and salinity over 2 rows and 
      !                                                   !  2 time step (now and before)
#if defined key_obc_trc
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:,:,:) ::   trcwbnd    !: West boundary passive tracers over 2 rows 
      !                                                    !  and 2 time step (now and before)
#endif
   REAL(wp), ALLOCATABLE, SAVE,     DIMENSION(:,:) ::   u_cxwbnd, v_cxwbnd   !: Zonal component of the phase speed ratio computed with 
      !                                                   !  radiation of zonal and meridional velocity (respectively) 
      !                                                   !  at the west open boundary (u_cxwbnd = cx rdt )
   REAL(wp), ALLOCATABLE, SAVE,     DIMENSION(:,:) ::   uwmsk, vwmsk, twmsk  !: 2D mask for the West OB

   ! Note that those arrays are optimized for mpp case 
   ! (hence the dimension jpj is the size of one processor subdomain)

   !!---------------------
   !! North open boundary
   !!---------------------
   INTEGER ::   nin0  , nin1       !: do loop index in mpp case for jpind, jpinf
   INTEGER ::   nin0p1, nin1m1     !: do loop index in mpp case for jpindp1, jpinfm1
   INTEGER ::   nin1m2, nin0m1     !: do loop index in mpp case for jpinfm1-1,jpind
   INTEGER ::   njn0  , njn1       !: do loop index in mpp case for jpnob
   INTEGER ::   njn0p1, njn1p1     !: do loop index in mpp case for jpnob+1
   INTEGER ::   njn0m1, njn1m1     !: do loop index in mpp case for jpnob-1

   REAL(wp), ALLOCATABLE, SAVE,   DIMENSION(:) ::   &  !:
      sshfon,           & !: now climatology of the north boundary sea surface height
      ubtfon,vbtfon       !: now climatology of the north boundary barotropic transport

   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:) ::   &    !:
      ufon, vfon,       & !: now climatology of the north boundary velocities
      tfon, sfon,       & !: now climatology of the north boundary temperature and salinity
      vclin               !: baroclinic componant of the meridian velocity after the radiation
      !                   !  in yhe obcdyn.F90 routine
   ! PERIANT control at north boundary
   REAL (wp ) :: &        !: Observed barotropic transport 
      trp_indian_obs ,  & !: Indian ocean
      trp_pacif_obs  ,  & !: Pacific ocean
      trp_atl_obs         !: Atlantic ocean

#if defined key_obc_trc
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   &  !:
      trcfon               !: now climatology of the north boundary passive tracers
#endif
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:) ::   sshfon_b      !: north boundary ssh correction averaged over the barotropic loop
      !                                            !  (if Flather's algoritm applied at open boundary)

   !!--------------------------------
   !! Arrays for radiative North OBC
   !!--------------------------------
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:,:) ::   unbnd, vnbnd      !: baroclinic u & v components of the velocity over 3
      !                                                    !  rows and 3 time step (now, before, and before before)
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:,:) ::   tnbnd, snbnd      !: north boundary temperature and salinity over
      !                                                    !  2 rows and 2 time step (now and before)
#if defined key_obc_trc
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:,:,:) ::   trcnbnd    !: North boundary passive tracers over 2 rows 
      !                                                    !  and 2 time step (now and before)
#endif
   REAL(wp), ALLOCATABLE, SAVE,     DIMENSION(:,:) ::   u_cynbnd, v_cynbnd    !: Meridional component of the phase speed ratio compu-
      !                                                    !  ted with radiation of zonal and meridional velocity 
      !                                                    !  (respectively) at the north OB (u_cynbnd = cx rdt )
   REAL(wp), ALLOCATABLE, SAVE,     DIMENSION(:,:) ::   unmsk, vnmsk, tnmsk   !: 2D mask for the North OB

   ! Note that those arrays are optimized for mpp case 
   ! (hence the dimension jpj is the size of one processor subdomain)
   
   !!---------------------
   !! South open boundary
   !!---------------------
   INTEGER ::   nis0  , nis1       !: do loop index in mpp case for jpisd, jpisf
   INTEGER ::   nis0p1, nis1m1     !: do loop index in mpp case for jpisdp1, jpisfm1
   INTEGER ::   nis1m2, nis0m1     !: do loop index in mpp case for jpisfm1-1,jpisd
   INTEGER ::   njs0  , njs1       !: do loop index in mpp case for jpsob
   INTEGER ::   njs0p1, njs1p1     !: do loop index in mpp case for jpsob+1

   REAL(wp), ALLOCATABLE, SAVE,   DIMENSION(:) ::    &   !:
      sshfos,           & !: now climatology of the south boundary sea surface height
      ubtfos,vbtfos       !: now climatology of the south boundary barotropic transport

   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:) ::    &   !:
      ufos, vfos,       & !: now climatology of the south boundary velocities 
      tfos, sfos,       & !: now climatology of the south boundary temperature and salinity
      vclis               !: baroclinic componant of the meridian velocity after the radiation 
      !                   !  in the obcdyn.F90 routine

#if defined key_obc_trc
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   &  !:
      trcfos               !: now climatology of the south boundary passive tracers
#endif

   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:) ::   sshfos_b     !: south boundary ssh correction averaged over the barotropic loop
      !                                           !  (if Flather's algoritm applied at open boundary)

   !!--------------------------------
   !! Arrays for radiative South OBC   (computed by the forward time step in dynspg)
   !!--------------------------------
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:,:) ::   usbnd, vsbnd     !: baroclinic u & v components of the velocity over 3 
      !                                                   !  rows and 3 time step (now, before, and before before)
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:,:) ::   tsbnd, ssbnd     !: south boundary temperature and salinity over
      !                                                   !  2 rows and 2 time step (now and before)
#if defined key_obc_trc
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:,:,:) ::   trcsbnd    !: South boundary passive tracers over 2 rows 
      !                                                    !  and 2 time step (now and before)
#endif
   REAL(wp), ALLOCATABLE, SAVE,     DIMENSION(:,:) ::   u_cysbnd, v_cysbnd   !: Meridional component of the phase speed ratio
      !                                                   !  computed with radiation of zonal and meridional velocity 
      !                                                   !  (repsectively) at the south OB (u_cynbnd = cx rdt )
   REAL(wp), ALLOCATABLE, SAVE,     DIMENSION(:,:) ::   usmsk, vsmsk, tsmsk  !: 2D mask for the South OB

   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: obc_oce.F90 2715 2011-03-30 15:58:35Z rblod $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   INTEGER FUNCTION obc_oce_alloc()
      !!----------------------------------------------------------------------
      !!               ***  FUNCTION obc_oce_alloc  ***
      !!----------------------------------------------------------------------

      ALLOCATE(                                                               &
              !! East open boundary
              obctmsk(jpi,jpj), obcumask(jpi,jpj), obcvmask(jpi,jpj),        &
              sshfoe(jpjed:jpjef), ubtfoe(jpjed:jpjef), vbtfoe(jpjed:jpjef), &
              ufoe(jpj,jpk), vfoe(jpj,jpk), tfoe(jpj,jpk), sfoe(jpj,jpk),    &
              uclie(jpj,jpk), sshfoe_b(jpjed:jpjef,jpj),                     &
              !! Arrays for radiative East OBC
              uebnd(jpj,jpk,3,3), vebnd(jpj,jpk,3,3) ,                       &
              tebnd(jpj,jpk,2,2), sebnd(jpj,jpk,2,2),                        &
              u_cxebnd(jpj,jpk), v_cxebnd(jpj,jpk),                          &
              uemsk(jpj,jpk), vemsk(jpj,jpk), temsk(jpj,jpk),                &
              !! West open boundary
              sshfow(jpjwd:jpjwf), ubtfow(jpjwd:jpjwf), vbtfow(jpjwd:jpjwf), &
              ufow(jpj,jpk), vfow(jpj,jpk), tfow(jpj,jpk),                   &
              sfow(jpj,jpk), ucliw(jpj,jpk), sshfow_b(jpjwd:jpjwf,jpj),      &
              !! Arrays for radiative West OBC
              uwbnd(jpj,jpk,3,3), vwbnd(jpj,jpk,3,3),                        &
              twbnd(jpj,jpk,2,2), swbnd(jpj,jpk,2,2),                        &
              u_cxwbnd(jpj,jpk), v_cxwbnd(jpj,jpk),                          &
              uwmsk(jpj,jpk), vwmsk(jpj,jpk), twmsk(jpj,jpk),                &
              !! North open boundary
              sshfon(jpind:jpinf), ubtfon(jpind:jpinf), vbtfon(jpind:jpinf), &
              ufon(jpi,jpk), vfon(jpi,jpk), tfon(jpi,jpk),                   &
              sfon(jpi,jpk), vclin(jpi,jpk), sshfon_b(jpind:jpinf,jpj),      &
              !! Arrays for radiative North OBC
              unbnd(jpi,jpk,3,3), vnbnd(jpi,jpk,3,3),                        &
              tnbnd(jpi,jpk,2,2), snbnd(jpi,jpk,2,2),                        &
              u_cynbnd(jpi,jpk), v_cynbnd(jpi,jpk),                          &
              unmsk(jpi,jpk), vnmsk(jpi,jpk), tnmsk (jpi,jpk),               &
              !! South open boundary
              sshfos(jpisd:jpisf), ubtfos(jpisd:jpisf), vbtfos(jpisd:jpisf), &
              ufos(jpi,jpk), vfos(jpi,jpk), tfos(jpi,jpk),                   &
              sfos(jpi,jpk), vclis(jpi,jpk),                                 &
              sshfos_b(jpisd:jpisf,jpj),                                     &
              !! Arrays for radiative South OBC 
              usbnd(jpi,jpk,3,3), vsbnd(jpi,jpk,3,3),                        &
              tsbnd(jpi,jpk,2,2), ssbnd(jpi,jpk,2,2),                        &
              u_cysbnd(jpi,jpk), v_cysbnd(jpi,jpk),                          &
              usmsk(jpi,jpk), vsmsk(jpi,jpk), tsmsk(jpi,jpk),                &
              !!
              STAT=obc_oce_alloc )
      !
   END FUNCTION obc_oce_alloc

   INTEGER FUNCTION obc_oce_trc_alloc()
      !!----------------------------------------------------------------------
      !!               ***  FUNCTION obc_oce_trc_alloc  ***
      !!----------------------------------------------------------------------
      ALLOCATE(                                                               &
              !! East open boundary
              trcfoe(jpj,jpk,jptra),                                         &
              !! Arrays for radiative East OBC
              trcebnd(jpj,jpk,2,2,jptra),                                     &
              !! West open boundary
              trcfow(jpj,jpk,jptra),                                         &
              !! Arrays for radiative West OBC
              trcwbnd(jpj,jpk,2,2,jptra),                                     &
              !! North open boundary
              trcfon(jpi,jpk,jptra),                                         &
              !! Arrays for radiative North OBC
              trcnbnd(jpi,jpk,2,2,jptra),                                     &
              !! South open boundary
              trcfos(jpi,jpk,jptra),                                         &
              !! Arrays for radiative South OBC 
              trcsbnd(jpi,jpk,2,2,jptra),                                     &
              !!
              STAT=obc_oce_trc_alloc )
   END FUNCTION obc_oce_trc_alloc
   
#else
   !!----------------------------------------------------------------------
   !!   Default option         Empty module                          No OBC
   !!----------------------------------------------------------------------
#endif

   !!======================================================================
END MODULE obc_oce
