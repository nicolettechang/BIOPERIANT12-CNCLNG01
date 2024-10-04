MODULE obcdta_trc
  !!==================================================================================
  !!                            ***  MODULE obcdta_trc  ***
  !! Open boundary data : read the data for the open boundaries.
  !! $hist: C. Dufour 2009-10-07 obcdta_trc inspired from obcdta
  !!==================================================================================
#if defined key_obc_trc
  !!----------------------------------------------------------------------------------
  !!   'key_obc_trc'         :                                Open Boundary Conditions
  !!----------------------------------------------------------------------------------
  !!   obc_dta_trc           : read passive tracers data along each open boundary
  !!----------------------------------------------------------------------------------
  !! * Modules used
   USE oce             ! ocean dynamics and tracers 
   USE trc             ! passive tracer
   USE dom_oce         ! ocean space and time domain
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE phycst          ! physical constants
   USE obc_par         ! ocean open boundary conditions
   USE obc_oce         ! ocean open boundary conditions
   USE in_out_manager  ! I/O logical units
   USE lib_mpp         ! distributed memory computing
   USE dynspg_oce      ! ocean: surface pressure gradient
   USE ioipsl          ! now only for  ymds2ju function 
   USE iom             ! 

   IMPLICIT NONE
   PRIVATE

   PUBLIC   obc_dta_trc         ! routine  called by step.F90
   PUBLIC   obc_dta_trc_alloc   ! function called by obcini.F90

   !$AGRIF_DO_NOT_TREAT
   REAL(wp),  DIMENSION(2)              ::   zjcnes_obc   ! 
   REAL(wp),  DIMENSION(:), ALLOCATABLE :: ztcobc_tmp
   REAL(wp),  DIMENSION(75,jptra) :: ztcobc=0
   !$AGRIF_END_DO_NOT_TREAT
   REAL(wp) :: rdt_obc
   REAL(wp) :: zjcnes
   INTEGER :: imm0, iyy0, idd0, iyy, imm, idd
   INTEGER :: nt_a=2, nt_b=1, itobc, ndate0_cnes, nday_year0
   INTEGER ::  itobce, itobcw, itobcs, itobcn, itobc_b  ! number of time steps in OBC files

  INTEGER                   :: ntobc     !:  where we are in the obc file
  INTEGER, DIMENSION(jptra) :: ntobc_b   !:  first record used
  INTEGER, DIMENSION(jptra) :: ntobc_a   !:  second record used

   CHARACTER (len=40) ::   cl_obc_eTR  ! name of data files
   CHARACTER (len=40) ::   cl_obc_wTR  !   -       -
   CHARACTER (len=40) ::   cl_obc_nTR  !   -       -
   CHARACTER (len=40) ::   cl_obc_sTR  !   -       -

   ! bt arrays for interpolating time dependent data on the boundaries
   INTEGER ::   nt_m=0, ntobc_m
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:)   :: trcbtedta    ! East
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:)   :: trcbtwdta    ! West
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:)   :: trcbtndta    ! North
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:)   :: trcbtsdta    ! South
   ! arrays used for interpolating time dependent data on the boundaries
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:,:) :: trcedta    ! East
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:,:) :: trcwdta    ! West
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:,:) :: trcndta    ! North
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:,:) :: trcsdta    ! South

   ! Masks set to .TRUE. after successful allocation below
   LOGICAL , ALLOCATABLE, SAVE, DIMENSION(:,:)   :: ltemsk  ! boolean msks
   LOGICAL , ALLOCATABLE, SAVE, DIMENSION(:,:)   :: ltwmsk  ! used for outliers
   LOGICAL , ALLOCATABLE, SAVE, DIMENSION(:,:)   :: ltnmsk  ! checks
   LOGICAL , ALLOCATABLE, SAVE, DIMENSION(:,:)   :: ltsmsk

   !! * Substitutions
#  include "obc_vectopt_loop_substitute.h90"
#  include "domzgr_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: obcdta.F90 3116 2011-11-15 20:55:40Z cetlod $
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS
   INTEGER FUNCTION obc_dta_trc_alloc()
      !!-------------------------------------------------------------------
      !!                     ***  ROUTINE obc_dta_trc_alloc  ***
      !!-------------------------------------------------------------------
      INTEGER :: ierr(2)
      !!-------------------------------------------------------------------
# if defined key_dynspg_ts
      ALLOCATE(   &     ! time-splitting : 0:jptobc
         ! bt arrays for interpolating time dependent data on the boundaries
         &      trcbtedta  (jpj,0:jptobc,jptra)  ,    &
         &      trcbtwdta  (jpj,0:jptobc,jptra)  ,    &
         &      trcbtndta  (jpi,0:jptobc,jptra)  ,    &
         &      trcbtsdta  (jpi,0:jptobc,jptra)  ,    &
         ! arrays used for interpolating time dependent data on the boundaries
         &      trcedta(jpj,jpk,0:jptobc,jptra) ,      &
         &      trcwdta(jpj,jpk,0:jptobc,jptra) ,      &
         &      trcndta(jpi,jpk,0:jptobc,jptra) ,      &
         &      trcsdta(jpi,jpk,0:jptobc,jptra) ,  STAT=ierr(1) )
# else
      ALLOCATE(   &     ! no time splitting : 1:jptobc
         ! bt arrays for interpolating time dependent data on the boundaries
         &      trcbtedta  (jpj,jptobc,jptra)   ,     &
         &      trcbtwdta  (jpj,jptobc,jptra)   ,     &
         &      trcbtndta  (jpi,jptobc,jptra)   ,     &
         &      trcbtsdta  (jpi,jptobc,jptra)   ,     &
         ! arrays used for interpolating time dependent data on the boundaries
         &      trcedta(jpj,jpk,jptobc,jptra) ,      &
         &      trcwdta(jpj,jpk,jptobc,jptra) ,      &
         &      trcndta(jpi,jpk,jptobc,jptra) ,      &
         &      trcsdta(jpi,jpk,jptobc,jptra) ,  STAT=ierr(1) )
# endif

      ALLOCATE( ltemsk(jpj,jpk)  ,     &
         &      ltwmsk(jpj,jpk)  ,     &
         &      ltnmsk(jpi,jpk)  ,     &
         &      ltsmsk(jpi,jpk)  , STAT=ierr(2) )

      obc_dta_trc_alloc = MAXVAL( ierr )
      IF( lk_mpp )   CALL mpp_sum( obc_dta_trc_alloc )

      IF( obc_dta_trc_alloc == 0 )  THEN         ! Initialise mask values following successful allocation
         !      east            !          west            !          north           !          south           !
         ltemsk(:,:) = .TRUE.   ;   ltwmsk(:,:) = .TRUE.   ;   ltnmsk(:,:) = .TRUE.   ;   ltsmsk(:,:) = .TRUE.
      END IF
      !
   END FUNCTION obc_dta_trc_alloc

  SUBROUTINE obc_dta_trc( kt )
    !!----------------------------------------------------------------------------------
    !!                      ***  SUBROUTINE obc_dta_trc  ***
    !!                    
    !! ** Purpose :   Call obc_dta_trc_single which computes OBC for each tracer 
    !!                since PISCES tracers have different conditions at boundary
    !!                ex: CFC : obc file with 75 time steps (output from global)
    !!                    DIC,Alk, etc : monthly climatology
    !!                    CaCO3, POC, etc: constantes
    !!
    !! History :
    !!        !  10-07 (C. Dufour, R. Dussin) adaptation to PISCES tracers
    !!-----------------------------------------------------------------------------------
    !! * Arguments
    INTEGER, INTENT( in ) ::   kt                  ! ocean time-step index
    !! * Local declarations  
    INTEGER ::   jt , ji , jj                      ! passive tracer index
    INTEGER, DIMENSION(jptra) :: zitobc_trc  ! contains value of itobc for each tracer

    IF(lwp) WRITE(numout,*)
    IF(lwp) WRITE(numout,*)  'obc_dta_trc : find boundary data for passive tracers'
    IF(lwp) WRITE(numout,*)  '~~~~~~~'

    IF (kt == nit000 ) zitobc_trc(:)=0

    DO jt = 1,jptra
       CALL obc_dta_trc_single(kt, jt, zitobc_trc)
       IF (kt == nit000 )  zitobc_trc(jt) = itobc
    ENDDO
  END SUBROUTINE obc_dta_trc

  SUBROUTINE obc_dta_trc_single(kt, jc, zitobc_trc)
    !!-----------------------------------------------------------------------------------
    !!                      ***  SUBROUTINE obc_dta_trc  ***
    !!                    
    !! ** Purpose :   Find the climatological  boundary arrays for the specified date, 
    !!                The boundary arrays are netcdf files. Three possible cases: 
    !!                - one time frame only in the file (time dimension = 1).
    !!                in that case the boundary data does not change in time.
    !!                - many time frames. In that case,  if we have 12 frames
    !!                we assume monthly fields. 
    !!                Else, we assume that time_counter is in seconds 
    !!                since the beginning of either the current year or a reference
    !!                year given in the namelist.
    !!                (no check is done so far but one would have to check the "unit"
    !!                 attribute of variable time_counter).
    !!
    !!
    !! History :
    !!        !  98-05 (J.M. Molines) Original code
    !!   8.5  !  02-10 (C. Talandier, A-M. Treguier) Free surface, F90
    !!
    !!   9.0  !  04-06 (F. Durand, A-M. Treguier) Netcdf BC files on input
    !!        !  2007-2008 (C. Langlais, P. Mathiot, J.M. Molines) high frequency boundaries data
    !!        !  09-11 (C. Dufour) adapted for passive tracers
    !!        !  10-07 (C. Dufour) adaptation to PISCES tracers
    !!-----------------------------------------------------------------------------------
    !! * Arguments
    INTEGER, INTENT( in ) ::   kt , jc          ! ocean time-step index
    !! * Local declarations
    INTEGER, DIMENSION(jptra) :: zitobc_trc  ! contains value of itobc for each tracer
    INTEGER ::   it  ! dummy loop indices
!    INTEGER ::  ikprint                                ! frequency for printouts.
    INTEGER, SAVE :: immfile, iyyfile                     !
    INTEGER :: nt              !  record indices (incrementation)
!    INTEGER :: istop           ! local error check

    REAL(wp) ::   zxy, znum, zden ! time interpolation weight

    LOGICAL :: lldebug = .true.
    ! variables for the julian day calculation
!    INTEGER :: iyear, imonth, iday
    REAL(wp) :: zsec

    ! IOM STUFF
!    INTEGER ::  idvar, id_e, id_n, id_x       ! file identifiers
!    INTEGER, DIMENSION(1)  :: itmp
!    CHARACTER(LEN=25) :: cl_vname

    !!---------------------------------------------------------------------------

    ! 0.  initialisation :
    ! --------------------
    IF( lldebug ) THEN
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*)  'obc_dta_trc_single : find boundary data for each tracer'
      IF(lwp) WRITE(numout,*)  '~~~~~~~'
    ENDIF

    itobc = zitobc_trc(jc)

    IF ( kt == nit000  )  CALL obc_dta_ini_trc_single ( kt, jc )
    IF ( nobc_dta_trc(jc) == 0 ) THEN ! already done in obc_dta_ini
       RETURN
    ENDIF
    IF ( itobc == 1    )   THEN       ! case of only one time frame in file done in obc_dta_ini
       RETURN
    ENDIF

    ! in the following code, we assume that obc data are read from files, with more than 1 time frame in it

    iyyfile=iyy ; immfile = 00  ! set component of the current file name
    IF ( cffile_trc /= 'annual') immfile = imm   ! 
    IF ( ln_obc_clim_trc(jc)   ) iyyfile = 0000  ! assume that climatological files are labeled y0000

    ! 1. Synchronize time of run with time of data files
    !---------------------------------------------------
    ! nday_year is the day number in the current year ( 1 for 01/01 )
    zsec=MOD( (kt-nit000)*rdt - (nday_year - nday_year0 )*rday, rday ) ! number of seconds in the current day
    IF (ln_obc_clim_trc(jc))  THEN
      zjcnes = nday_year - 1  + zsec/rday
    ELSE
      ! CD: debug -> recompute zjcnes0 before iterating
      zjcnes = ndate0_cnes + (nday_year - nday_year0 ) + zsec/rday
      zjcnes = zjcnes + rdt/rday
    ENDIF

    ! look for 'before' record number in the current file
    ntobc = nrecbef_trc_single (jc)  ! this function return the record number for 'before', relative to zjcnes

    IF (MOD(kt-1,10)==0) THEN
       IF (lwp) WRITE(numout,*) 'kt= ',kt,' zjcnes =', zjcnes,' ndastp =',ndastp, 'mm =',imm
    END IF

    ! 2. read a new data if necessary 
    !--------------------------------
    IF ( ntobc /= ntobc_b(jc) ) THEN
    ! we need to read the 'after' record
    ! swap working index:
    nt=nt_b ; nt_b=nt_a ; nt_a=nt
    ntobc_b(jc) = ntobc
    ! new record number :
    ntobc_a(jc) = ntobc_a(jc) + 1

    ! all tricky things related to record number, changing files etc... are managed by obc_read

    CALL obc_read_trc_single (kt, nt_a, ntobc_a(jc), iyyfile, immfile, jc )

    ! update zjcnes_obc
    zjcnes_obc(nt_b)= ztcobc(ntobc_b(jc), jc)
    zjcnes_obc(nt_a)= ztcobc(ntobc_a(jc), jc)
    ENDIF

    ! 3.   interpolation at each time step
    ! ------------------------------------
    IF ( ln_obc_clim_trc(jc)) THEN
      znum= MOD(zjcnes           - zjcnes_obc(nt_b), REAL(nyear_len(1),wp) ) ; IF ( znum < 0 ) znum = znum + REAL(nyear_len(1),wp)
      zden= MOD(zjcnes_obc(nt_a) - zjcnes_obc(nt_b), REAL(nyear_len(1),wp) ) ; IF ( zden < 0 ) zden = zden + REAL(nyear_len(1),wp)
    ELSE
      znum= zjcnes           - zjcnes_obc(nt_b)
      zden= zjcnes_obc(nt_a) - zjcnes_obc(nt_b)
    ENDIF 
    zxy = znum / zden

    IF( lp_obc_east_trc ) THEN
       !  fills  trfoe
!       DO jn = 1,jptra
         trcfoe(:,:,jc) = zxy * trcedta (:,:,nt_a,jc) + (1. - zxy)*trcedta(:,:,nt_b,jc)
!       END DO
    ENDIF

    IF( lp_obc_west_trc) THEN
       !  fills trfow
         trcfow(:,:,jc) = zxy * trcwdta (:,:,nt_a,jc) + (1. - zxy)*trcwdta(:,:,nt_b,jc)
    ENDIF

    IF( lp_obc_north_trc) THEN
       !  fills trfon
         trcfon(:,:,jc) = zxy * trcndta (:,:,nt_a,jc) + (1. - zxy)*trcndta(:,:,nt_b,jc)
    ENDIF

    IF( lp_obc_south_trc) THEN
       !  fills tfos
         trcfos(:,:,jc) = zxy * trcsdta (:,:,nt_a,jc) + (1. - zxy)*trcsdta(:,:,nt_b,jc)
    ENDIF

  END SUBROUTINE obc_dta_trc_single

  SUBROUTINE obc_dta_ini_trc_single (kt, jn)
    !!-----------------------------------------------------------------------------
    !!                       ***  SUBROUTINE obc_dta_ini_trc  ***
    !!
    !! ** Purpose :
    !!      When obc_dta_trc first call, realize some data initialization
    !!
    !! ** Method :
    !!
    !! History :
    !!   9.0  ! 07-10 (J.M. Molines )
    !!        ! 10-07 (C. Dufour) adaptation to PISCES tracers
    !!----------------------------------------------------------------------------
    !! * Argument
    INTEGER, INTENT(in)  :: kt      ! ocean time-step index
    INTEGER  :: jn      ! ocean time-step index

    !! * Local declarations
    INTEGER ::   ji,jj,jk,it   ! dummy loop indices

    REAL(wp) ::   zxy                                    ! time interpolation weight


    INTEGER ::   inum               ! temporary logical unit
    INTEGER ::   ibloc, nreclo, jrec
    INTEGER ::   jfoe, jfow, ifon, ifos
    INTEGER, SAVE :: immfile, iyyfile                     !
    INTEGER :: nt              !  record indices (incrementation)
    LOGICAL :: lldebug = .true.
    ! variables for the julian day calculation
    INTEGER :: iyear, imonth, iday
    REAL(wp) :: zsec , zjulian, zjuliancnes


    IF( lldebug ) THEN
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*)  'obc_dta_ini_trc_single : initialization of boundary data for each tracer'
      IF(lwp) WRITE(numout,*)  '~~~~~~~'
      IF (lwp) THEN
         IF ( nobc_dta_trc(jn) == 0 ) THEN
            WRITE(numout,*)  '          OBC data taken from initial conditions.'
         ELSE
            WRITE(numout,*)  '          OBC data taken from netcdf files.'
         ENDIF
      ENDIF
    ENDIF
    nday_year0 = nday_year  ! to remember the day when kt=nit000

    trcedta(:,:,:,jn) = 0.e0  ! East
    trcwdta(:,:,:,jn) = 0.e0  ! West
    trcndta(:,:,:,jn) = 0.e0  ! North
    trcsdta(:,:,:,jn) = 0.e0  ! South

    trcfoe(:,:,jn) = 0.e0     ! East
    trcfow(:,:,jn) = 0.e0     ! West
    trcfon(:,:,jn) = 0.e0     ! North
    trcfos(:,:,jn) = 0.e0     ! South
    ibloc  = 4096
    nreclo = ibloc*( ( jpk *jpbyt -1)/ibloc + 1)

    IF (nobc_dta_trc(jn) == 0 ) THEN   ! boundary data are the initial data of this run (set only at nit000)

       IF (lp_obc_east_trc) THEN  ! East
           DO ji = nie0 , nie1
             trcfoe(nje0:nje1,:,jn) = temsk(nje0:nje1,:) * trn (ji+1 , nje0:nje1 , :,jn) * tmask(ji+1,nje0:nje1 , :)
           END DO

         ! save trcfoe into a climatological type file if ln_obc_ini is true
         IF (ln_obc_ini(jn)) THEN
           IF(lwp) WRITE(numout,*)  ' ln_obc_ini = TRUE -> initial state is saved like a climatology for ', ctrcnm(jn)
           CALL ctl_opn( inum, 'obc_east_'//TRIM(ctrcnm(jn))//'_y0000m00', 'UNKNOWN', 'UNFORMATTED', 'DIRECT',   &
           &         nreclo, numout, lwp, 0 )
           ! Write header
           ! ----------------
           WRITE (inum,REC=1) nreclo,jpjglo,jpk,jpjed,jpjef
           IF( jpieob /= 0 ) THEN
             IF( nje0+njmpp-1  == jpjed .AND. nie1 >= nie0 ) THEN
               ! ... dump of jpjed if it is on this proc.
               jrec = 2
               jfoe = jpjed - njmpp + 1
               WRITE(inum,REC=jrec)( trcfoe(jfoe,jk,jn),jk=1,jpk)
             ENDIF
             DO ji = nie0, nie1
               DO jj = nje0, nje1
                 ! ... only interested processors go through the following lines
                 !           jfoe = jj + njmpp -1
                 jfoe = jj
                 jrec = 2 + jj + njmpp -1 -jpjed
                 WRITE (inum,REC=jrec)( trcfoe(jfoe,jk,jn),jk=1,jpk)
               END DO
             END DO
           END IF
         ENDIF
       ENDIF

       IF (lp_obc_west_trc) THEN  ! West
           DO ji = niw0 , niw1
             trcfow(njw0:njw1,:,jn) = twmsk(njw0:njw1,:) * trn (ji , njw0:njw1 , :,jn) * tmask(ji , njw0:njw1 , :)
           END DO

         ! save trcfow into a climatological type file if ln_obc_ini is true
         IF (ln_obc_ini(jn)) THEN
           IF(lwp) WRITE(numout,*)  ' ln_obc_ini = TRUE -> initial state is saved like a climatology for ', ctrcnm(jn)
           CALL ctl_opn( inum, 'obc_west_'//TRIM(ctrcnm(jn))//'_y0000m00', 'UNKNOWN', 'UNFORMATTED', 'DIRECT',   &
           &         nreclo, numout, lwp, 0 )
           ! Write header
           ! ----------------
!           WRITE (inum,REC=1) nreclo,jpjglo,jpk,jptra,jpjwd,jpjwf
           WRITE (inum,REC=1) nreclo,jpjglo,jpk,jpjwd,jpjwf
           IF( jpiwob /= 0 ) THEN
             IF( njw0+njmpp+1 == jpjwd .AND. niw1 >= niw0 ) THEN
               ! ... dump of jpjed if it is on this proc.
               jrec = 2
               jfow = jpjwd -njmpp + 1
!               WRITE(inum,REC=jrec)(( trcfow(jfow,jk,jn),jk=1,jpk),jn=1,jptra)
               WRITE(inum,REC=jrec)( trcfow(jfow,jk,jn),jk=1,jpk)
             ENDIF
             DO ji = niw0, niw1
                DO jj = njw0, njw1
                  ! ... only interested processors go through the following lines
                  !           jfow = jj + njmpp -1
                  jfow = jj
                  jrec = 2 + jj + njmpp -1 -jpjwd
!                  WRITE (inum,REC=jrec)(( trcfow(jfow,jk,jn),jk=1,jpk),jn=1,jptra)
                  WRITE (inum,REC=jrec)( trcfow(jfow,jk,jn),jk=1,jpk)
                END DO
             END DO
           END IF
         ENDIF
       ENDIF

       IF (lp_obc_north_trc) THEN ! North
           DO jj = njn0 , njn1
             trcfon(nin0:nin1,:,jn) = tnmsk(nin0:nin1,:) * trn (nin0:nin1 , jj+1 , :,jn) * tmask(nin0:nin1 , jj+1 , :)
           END DO
         ! save trcfon into a climatological type file if ln_obc_ini is true
         IF (ln_obc_ini(jn)) THEN
           IF(lwp) WRITE(numout,*)  ' ln_obc_ini = TRUE -> initial state is saved like a climatology for ', ctrcnm(jn)
           CALL ctl_opn( inum, 'obc_north_'//TRIM(ctrcnm(jn))//'_y0000m00', 'UNKNOWN', 'UNFORMATTED', 'DIRECT',   &
           &         nreclo, numout, lwp, 0 )
           ! Write header
           ! ----------------
!           WRITE (inum,REC=1) nreclo,jpiglo,jpk,jptra,jpind,jpinf
           WRITE (inum,REC=1) nreclo,jpiglo,jpk,jpind,jpinf
           IF( jpjnob /= 0 ) THEN
             IF( nin0+nimpp-1 == jpind .AND. njn1 >= njn0 ) THEN
                ! ... dump of jpind if it is on this proc.
                jrec = 2
                ifon = jpind -nimpp +1
!                WRITE(inum,REC=jrec)(( trcfon(ifon,jk,jn),jk=1,jpk),jn=1,jptra)
                WRITE(inum,REC=jrec)( trcfon(ifon,jk,jn),jk=1,jpk)
             ENDIF
             DO jj = njn0, njn1
                DO ji = nin0, nin1
                   ! ... only interested processors go through the following lines
                   !          ifon = ji + nimpp -1   
                   ifon = ji
                   jrec = 2 + ji + nimpp -1 -jpind
!                   WRITE (inum,REC=jrec)(( trcfon(ifon,jk,jn),jk=1,jpk),jn=1,jptra)
                   WRITE (inum,REC=jrec)( trcfon(ifon,jk,jn),jk=1,jpk)
                END DO
             END DO
           END IF
         ENDIF
       ENDIF

       IF (lp_obc_south_trc) THEN ! South
           DO jj = njs0 , njs1
             trcfos(nis0:nis1,:,jn) = tsmsk(nis0:nis1,:) * trn (nis0:nis1 , jj , :,jn) * tmask(nis0:nis1 , jj , :)
           END DO
         ! save trcfos into a climatological type file if ln_obc_ini is true
         IF (ln_obc_ini(jn)) THEN
           IF(lwp) WRITE(numout,*)  ' ln_obc_ini = TRUE -> initial state is saved like a climatology for  ', ctrcnm(jn)
           CALL ctl_opn( inum, 'obc_south_'//TRIM(ctrcnm(jn))//'_y0000m00', 'UNKNOWN', 'UNFORMATTED', 'DIRECT',   &
           &         nreclo, numout, lwp, 0 )
           ! Write header
           ! ----------------
!           WRITE (inum,REC=1) nreclo,jpiglo,jpk,jptra,jpisd,jpisf
           WRITE (inum,REC=1) nreclo,jpiglo,jpk,jpisd,jpisf
           IF( jpjsob /= 0 ) THEN
              IF( nis0+nimpp-1 == jpisd .AND. njs1 >= njs0 ) THEN
                 ! ... dump of jpisd if it is on this proc.
                 jrec = 2
                 ifos = jpisd -nimpp +1
!                 WRITE(inum,REC=jrec)(( trcfos(ifos,jk,jn),jk=1,jpk),jn=1,jptra)
                 WRITE(inum,REC=jrec)( trcfos(ifos,jk,jn),jk=1,jpk)
              ENDIF
              DO jj = njs0, njs1
                DO ji = nis0, nis1
                  ! ... only interested processors go through the following lines
                  !          ifos = ji + nimpp -1   
                  ifos = ji
                  jrec = 2 + ji + nimpp -1 -jpisd
!                  WRITE (inum,REC=jrec)(( trcfos(ifos,jk,jn),jk=1,jpk),jn=1,jptra)
                  WRITE (inum,REC=jrec)( trcfos(ifos,jk,jn),jk=1,jpk)
                END DO
              END DO
           END IF
         ENDIF
       ENDIF
       RETURN  ! exit the routine all is done
    ENDIF  ! nobc_dta_trc = 0 

!!!! In the following OBC data are read from files.
    ! all logical-mask are initialzed to true when declared
    WHERE ( temsk == 0 ) ltemsk=.FALSE.
    WHERE ( twmsk == 0 ) ltwmsk=.FALSE.
    WHERE ( tnmsk == 0 ) ltnmsk=.FALSE.
    WHERE ( tsmsk == 0 ) ltsmsk=.FALSE.


    iyear=1950;  imonth=01; iday=01;  zsec=0.
    ! zjuliancnes : julian day corresonding  to  01/01/1950
    CALL ymds2ju(iyear, imonth, iday,zsec , zjuliancnes)

    !current year and curent month 
    iyy=INT(ndastp/10000) ; imm=INT((ndastp -iyy*10000)/100) ; idd=(ndastp-iyy*10000-imm*100)
    IF (iyy <  1900)  iyy = iyy+1900  ! always assume that years are on 4 digits.
    CALL ymds2ju(iyy, imm, idd ,zsec , zjulian)
    ndate0_cnes = zjulian - zjuliancnes   ! jcnes day when call to obc_dta_ini

    iyyfile=iyy ; immfile=0  ! set component of the current file name
    IF ( cffile_trc /= 'annual') immfile=imm
    IF ( ln_obc_clim_trc(jn)) iyyfile = 0  ! assume that climatological files are labeled y0000

    CALL obc_dta_chktime_trc_single ( iyyfile, immfile, jn )

    IF ( itobc == 1 ) THEN
       ! in this case we will provide boundary data only once.
       nt_a=1 ; ntobc_a(jn)=1
       CALL obc_read_trc_single (nit000, nt_a, ntobc_a(jn), iyyfile, immfile, jn)
       IF( lp_obc_east_trc ) THEN
          !  fills trcfoe
          trcfoe(:,:,jn) =  trcedta (:,:,1,jn)
       ENDIF

       IF( lp_obc_west_trc) THEN
          !  fills trcfow
          trcfow(:,:,jn) =  trcwdta (:,:,1,jn)
       ENDIF

       IF( lp_obc_north_trc) THEN
          !  fills trcfon
          trcfon(:,:,jn) =  trcndta (:,:,1,jn)
       ENDIF

       IF( lp_obc_south_trc) THEN
          !  fills trcfos
          trcfos(:,:,jn) =  trcsdta (:,:,1,jn)
       ENDIF
       RETURN  ! we go out of obc_dta_ini_trc -------------------------------------->>>>>
    ENDIF

    ! nday_year is the day number in the current year ( 1 for 01/01 )
    ! we suppose that we always start from the begining of a day
    !    zsec=MOD( (kt-nit000)*rdt - (nday_year - nday_year0 )*rday, rday ) ! number of seconds in the current day
    zsec=0.e0  ! here, kt=nit000, nday_year = ndat_year0 

    IF (ln_obc_clim_trc(jn))  THEN
      zjcnes = nday_year - 1  + zsec/rday  ! for clim file time is in days in a year
    ELSE
      zjcnes = ndate0_cnes + (nday_year - nday_year0 ) + zsec/rday
    ENDIF

    ! look for 'before' record number in the current file
    ntobc = nrecbef_trc_single (jn)

    IF (lwp) WRITE(numout,*) 'obc files frequency :',cffile_trc
    IF (lwp) WRITE(numout,*) ' zjcnes0 =',zjcnes,' ndastp0 =',ndastp
    IF (lwp) WRITE(numout,*) ' annee0 ',iyy,' month0 ', imm,' day0 ', idd
    IF (lwp) WRITE(numout,*) 'first file open :',cl_obc_nTR

    ! record initialisation
    !--------------------
    nt_b = 1 ; nt_a = 2

    ntobc_a(jn) = ntobc + 1
    ntobc_b(jn) = ntobc

    CALL obc_read_trc_single (kt, nt_b, ntobc_b(jn), iyyfile, immfile, jn)  ! read 'before' fields
    CALL obc_read_trc_single (kt, nt_a, ntobc_a(jn), iyyfile, immfile, jn)  ! read 'after' fields

    zjcnes_obc(nt_b)= ztcobc(ntobc_b(jn),jn)
    zjcnes_obc(nt_a)= ztcobc(ntobc_a(jn),jn)
    ! 
  END SUBROUTINE obc_dta_ini_trc_single


  SUBROUTINE obc_dta_chktime_trc_single (kyyfile, kmmfile, jtrc)
   !
   ! check the number of time steps in the files and read ztcobc 
   !
   ! * Arguments
   INTEGER, INTENT(in) :: kyyfile, kmmfile, jtrc
   ! * local variables
   INTEGER :: istop       ! error control
   INTEGER :: ji          ! dummy loop index

    INTEGER ::  idvar, id_e, id_w, id_n, id_s, id_x       ! file identifiers
    INTEGER, DIMENSION(1)  :: itmp
    CHARACTER(LEN=25) :: cl_vname

    ntobc_a(jtrc) = 0; itobce =0 ; itobcw = 0; itobcn = 0; itobcs = 0
    ! build file name
    WRITE(cl_obc_eTR ,'("obc_east_'//TRIM(ctrcnm(jtrc))//'_y",i4.4,"m",i2.2,".nc")'  ) kyyfile,kmmfile
    WRITE(cl_obc_wTR ,'("obc_west_'//TRIM(ctrcnm(jtrc))//'_y",i4.4,"m",i2.2,".nc")'  ) kyyfile,kmmfile
    WRITE(cl_obc_nTR ,'("obc_north_'//TRIM(ctrcnm(jtrc))//'_y",i4.4,"m",i2.2,".nc")' ) kyyfile,kmmfile
    WRITE(cl_obc_sTR ,'("obc_south_'//TRIM(ctrcnm(jtrc))//'_y",i4.4,"m",i2.2,".nc")' ) kyyfile,kmmfile

    cl_vname = 'time_counter'
    IF ( lp_obc_east_trc ) THEN
       CALL iom_open ( cl_obc_eTR , id_e )
       idvar = iom_varid( id_e, cl_vname, kdimsz = itmp ); itobce=itmp(1)
    ENDIF
    IF ( lp_obc_west_trc ) THEN
       CALL iom_open ( cl_obc_wTR , id_w )
       idvar = iom_varid( id_w, cl_vname, kdimsz = itmp ) ; itobcw=itmp(1)
    ENDIF
    IF ( lp_obc_north_trc ) THEN
       CALL iom_open ( cl_obc_nTR , id_n )
       idvar = iom_varid( id_n, cl_vname, kdimsz = itmp ) ; itobcn=itmp(1)
    ENDIF
    IF ( lp_obc_south_trc ) THEN
       CALL iom_open ( cl_obc_sTR , id_s )
       idvar = iom_varid( id_s, cl_vname, kdimsz = itmp ) ; itobcs=itmp(1)
    ENDIF

    itobc = MAX( itobce, itobcw, itobcn, itobcs )
    istop = 0
    IF ( lp_obc_east_trc  .AND. itobce /= itobc ) istop = istop+1
    IF ( lp_obc_west_trc  .AND. itobcw /= itobc ) istop = istop+1
    IF ( lp_obc_north_trc .AND. itobcn /= itobc ) istop = istop+1
    IF ( lp_obc_south_trc .AND. itobcs /= itobc ) istop = istop+1
    nstop = nstop + istop

    IF ( istop /=  0 )  THEN
       WRITE(ctmp1,*) ' east, west, north, south: ', itobce, itobcw, itobcn, itobcs
       CALL ctl_stop( 'obcdta : all files must have the same number of time steps', ctmp1 )
    ENDIF

    IF ( itobc == 1 ) THEN
       IF (lwp) THEN
          WRITE(numout,*) ' obcdta found one time step only in the OBC files'
          IF (ln_obc_clim_trc(jtrc)) THEN
             ! OK no problem
          ELSE
             ln_obc_clim_trc(jtrc)=.true.
             WRITE(numout,*) ' we force ln_obc_clim_trc to T for ', ctrcnm(jtrc)
          ENDIF
       ENDIF
! CD: debug in ocean.output
       IF ( lp_obc_east_trc )  CALL iom_close (id_e)
       IF ( lp_obc_west_trc )  CALL iom_close (id_w)
       IF ( lp_obc_north_trc ) CALL iom_close (id_n)
       IF ( lp_obc_south_trc ) CALL iom_close (id_s)
!
    ELSE
       IF ( ALLOCATED(ztcobc_tmp) ) DEALLOCATE ( ztcobc_tmp )
       ALLOCATE (ztcobc_tmp(itobc))
       DO ji=1,1   ! use a dummy loop to read ztcobc only once
          IF ( lp_obc_east_trc ) THEN
             CALL iom_gettime ( id_e, ztcobc_tmp, cl_vname ) ; CALL iom_close (id_e) ; EXIT
          ENDIF
          IF ( lp_obc_west_trc ) THEN
             CALL iom_gettime ( id_w, ztcobc_tmp, cl_vname ) ; CALL iom_close (id_w) ; EXIT
          ENDIF
          IF ( lp_obc_north_trc ) THEN
             CALL iom_gettime ( id_n, ztcobc_tmp, cl_vname ) ; CALL iom_close (id_n) ; EXIT
          ENDIF
          IF ( lp_obc_south_trc ) THEN
             CALL iom_gettime ( id_s, ztcobc_tmp, cl_vname ) ; CALL iom_close (id_s) ; EXIT
          ENDIF
       END DO
       ztcobc(1:itobc,jtrc)=ztcobc_tmp(:)
       rdt_obc = ztcobc(2,jtrc)-ztcobc(1,jtrc)  !  just an information, not used for any computation
       IF (lwp) WRITE(numout,*) ' For tracer', ctrcnm(jtrc)
       IF (lwp) WRITE(numout,*) ' obcdta_trc found', itobc,' time steps in the OBC files'
       IF (lwp) WRITE(numout,*) ' time step of obc data :', rdt_obc,' days'
       IF (lwp) WRITE(numout,*) ' ztcobc(:,jtrc) :', ztcobc(:,jtrc)
     ENDIF
     zjcnes = zjcnes - rdt/rday  ! trick : zcnes is always incremented by rdt/rday in obc_dta!
  END SUBROUTINE obc_dta_chktime_trc_single

  !!==============================================================================
  SUBROUTINE obc_read_trc_single (kt, nt_x, ntobc_x, iyy, imm, jtra)
     !!-------------------------------------------------------------------------
     !!                      ***  ROUTINE obc_read_trc  ***
     !!
     !! ** Purpose :  Read the boundary data in files identified by iyy and imm
     !!               According to the validated open boundaries, return the 
     !!               following arrays :
     !!                trcedta : East OBC passive tracers
     !!                trcwdta : West OBC passive tracers
     !!                trcndta : North OBC passive tracers
     !!                trcsdta : South OBC passive tracers
     !!
     !! ** Method  :  These fields are read in the record ntobc_x of the files.
     !!               The number of records is already known. If  ntobc_x is greater
     !!               than the number of record, this routine will look for next file,
     !!               updating the indices (case of inter-annual obcs) or loop at the
     !!               begining in case of climatological file (ln_obc_clim_trc = true ).
     !! -------------------------------------------------------------------------
     !! History:     !  2005  ( P. Mathiot, C. Langlais ) Original code
     !!              !  2008  ( J,M, Molines ) Use IOM and cleaning
     !!              !  2009  (C. Dufour) adapted for passive tracers 
     !!              !  2010  (C. Dufour) adaptation to PISCES tracers
     !!--------------------------------------------------------------------------

    ! * Arguments
    INTEGER, INTENT( in ) :: kt, nt_x, jtra
    INTEGER, INTENT( inout ) :: ntobc_x , iyy, imm      ! yes ! inout !

    ! * Local variables
    CHARACTER (len=40) :: &    ! file names
         cl_obc_eTR , cl_obc_wTR , cl_obc_nTR , cl_obc_sTR

    INTEGER :: ikprint
!    REAL(wp) :: zmin, zmax   ! control of boundary values

    !IOM stuff
    INTEGER :: id_e, id_w, id_n, id_s
    INTEGER, DIMENSION(2) :: istart, icount
    LOGICAL :: lldebug = .true.

    !--------------------------------------------------------------------------
    IF ( ntobc_x > itobc ) THEN
      IF (ln_obc_clim_trc(jtra)) THEN  ! just loop on the same file
        ntobc_x = 1
      ELSE
        ! need to change file : it is always for an 'after' data
        IF ( cffile_trc == 'annual' ) THEN ! go to next year file
          iyy = iyy + 1
        ELSE IF ( cffile_trc =='monthly' ) THEN  ! go to next month file
          imm = imm + 1
          IF ( imm == 13 ) THEN
            imm = 1 ; iyy = iyy + 1
          ENDIF
        ELSE
         ctmp1='obcread_trc : this type of obc file is not supported :( '
         ctmp2=TRIM(cffile_trc)
         CALL ctl_stop (ctmp1, ctmp2)
         ! cffile_trc should be either annual or monthly ...
        ENDIF
       ! as the file is changed, need to update itobc etc ...
        CALL obc_dta_chktime_trc_single (iyy,imm,jtra)
        ntobc_x = nrecbef_trc_single(jtra) + 1 ! remember : this case occur for an after data
      ENDIF
    ENDIF

    IF ( lp_obc_east_trc ) THEN
    IF( lldebug ) THEN
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*)  'cl_obc_eTR : ',cl_obc_eTR 
      IF(lwp) WRITE(numout,*)  'iyy = ',iyy
      IF(lwp) WRITE(numout,*)  'imm = ',imm
    ENDIF
       ! ... Read datafile and set passive tracers values 
       ! ... initialise the tredta arrays
       WRITE(cl_obc_eTR ,'("obc_east_'//TRIM(ctrcnm(jtra))//'_y"  ,i4.4,"m",i2.2,".nc")' ) iyy,imm
       ! JMM this may change depending on the obc data format ...
       istart(:)=(/nje0+njmpp-1,1/) ; icount(:)=(/nje1-nje0 +1,jpk/)
       IF (lwp) WRITE(numout,*) 'read data in :', TRIM(cl_obc_eTR)
       ! Several files for tracers?
       IF (nje1 >= nje0 ) THEN
          CALL iom_open ( cl_obc_eTR , id_e )
          CALL iom_get ( id_e, jpdom_unknown, TRIM(ctrcnm(jtra)), trcedta(nje0:nje1,:,nt_x,jtra), &
                 &               ktime=ntobc_x , kstart=istart, kcount= icount )
          CALL iom_close (id_e)

          ! mask the boundary values
            trcedta(:,:,nt_x,jtra) = trcedta(:,:,nt_x,jtra)*temsk(:,:)

          ! check any outliers 
          !! Depends on tracers... Important for PISCES ???
!          zmin=MINVAL( tredta(:,:,nt_x,jn), mask=ltemsk ) ; zmax=MAXVAL(tredta(:,:,nt_x,jn), mask=ltemsk)
!          IF (  zmin < -10. .OR. zmax > 40)   THEN
!             CALL ctl_stop('Error in tedta',' routine obcdta')
!          ENDIF

          !               Usually printout is done only once at kt = nit000, unless nprint (namelist) > 1      
          IF ( lwp .AND.  ( kt == nit000 .OR. nprint /= 0 )  ) THEN
             WRITE(numout,*)
             WRITE(numout,*) ' Read East OBC data records ', ntobc_x
             ikprint = jpj/20 +1
             WRITE(numout,*) ' Passive tracer  record 1 - printout every 3 level'
             CALL prihre( trcedta(:,:,nt_x,:), jpj, jpk, 1, jpj, ikprint, jpk, 1, -3, 1., numout )
             WRITE(numout,*)
          ENDIF
       ENDIF
    ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IF ( lp_obc_west_trc ) THEN
       ! ... Read datafile and set passive tracers values
       ! ... initialise the trwdta arrays
       WRITE(cl_obc_wTR ,'("obc_west_'//TRIM(ctrcnm(jtra))//'_y"  ,i4.4,"m",i2.2,".nc")' ) iyy,imm
       istart(:)=(/njw0+njmpp-1,1/) ; icount(:)=(/njw1-njw0 +1,jpk/)
       IF (lwp) WRITE(numout,*) 'read data in :', TRIM(cl_obc_wTR)
       !! Several files???
       IF ( njw1 >= njw0 ) THEN
          CALL iom_open ( cl_obc_wTR , id_w )
          CALL iom_get ( id_w, jpdom_unknown, TRIM(ctrcnm(jtra)), trcwdta(njw0:njw1,:,nt_x,jtra), &
                 &               ktime=ntobc_x , kstart=istart, kcount= icount )
          CALL iom_close ( id_w )

          ! mask the boundary values
            trcwdta(:,:,nt_x,jtra) = trcwdta(:,:,nt_x,jtra)*twmsk(:,:)

          ! check any outliers
!          zmin=MINVAL( twdta(:,:,nt_x), mask=ltwmsk ) ; zmax=MAXVAL(twdta(:,:,nt_x), mask=ltwmsk)
!          IF (  zmin < -10. .OR. zmax > 40)   THEN
!             CALL ctl_stop('Error in twdta',' routine obcdta')
!          ENDIF


          IF ( lwp .AND.  ( kt == nit000 .OR. nprint /= 0 )  ) THEN
             WRITE(numout,*)
             WRITE(numout,*) ' Read West OBC data records ', ntobc_x
             ikprint = jpj/20 +1
             WRITE(numout,*) ' Passive tracers  record 1 - printout every 3 level'
             CALL prihre( trcwdta(:,:,nt_x,:), jpj, jpk, 1, jpj, ikprint, jpk, 1, -3, 1., numout )
             WRITE(numout,*)
          ENDIF
       END IF
    ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IF( lp_obc_north_trc) THEN
       WRITE(cl_obc_nTR ,'("obc_north_'//TRIM(ctrcnm(jtra))//'_y" ,i4.4,"m",i2.2,".nc")' ) iyy,imm
       istart(:)=(/nin0+nimpp-1,1/) ; icount(:)=(/nin1-nin0 +1,jpk/)
       IF (lwp) WRITE(numout,*) 'read data in :', TRIM(cl_obc_nTR)
       IF ( nin1 >= nin0 ) THEN
          CALL iom_open ( cl_obc_nTR , id_n )
          CALL iom_get ( id_n, jpdom_unknown, TRIM(ctrcnm(jtra)), trcndta(nin0:nin1,:,nt_x,jtra), &
                 &               ktime=ntobc_x , kstart=istart, kcount= icount )
          CALL iom_close (id_n)

          ! mask the boundary values
            trcndta(:,:,nt_x,jtra) = trcndta(:,:,nt_x,jtra)*tnmsk(:,:)

          ! check any outliers
!          zmin=MINVAL( tndta(:,:,nt_x), mask=ltnmsk ) ; zmax=MAXVAL(tndta(:,:,nt_x), mask=ltnmsk)
!          IF (  zmin < -10. .OR. zmax > 40)   THEN
!             CALL ctl_stop('Error in tndta',' routine obcdta')
!          ENDIF

          IF ( lwp .AND.  ( kt == nit000 .OR. nprint /= 0 )  ) THEN
             WRITE(numout,*)
             WRITE(numout,*) ' Read North OBC data records ', ntobc_x
             ikprint = jpi/20 +1
             WRITE(numout,*) ' Passive tracers  record 1 - printout every 3 level'
             CALL prihre( trcndta(:,:,nt_x,:), jpi, jpk, 1, jpi, ikprint, jpk, 1, -3, 1., numout )
          ENDIF
       ENDIF
    ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IF( lp_obc_south_trc) THEN
       WRITE(cl_obc_sTR ,'("obc_south_'//TRIM(ctrcnm(jtra))//'_y" ,i4.4,"m",i2.2,".nc")' ) iyy,imm
       istart(:)=(/nis0+nimpp-1,1/) ; icount(:)=(/nis1-nis0 +1,jpk/)
       IF (lwp) WRITE(numout,*) 'read data in :', TRIM(cl_obc_sTR)
       IF ( nis1 >= nis0 ) THEN
          CALL iom_open ( cl_obc_sTR , id_s )
          CALL iom_get ( id_s, jpdom_unknown, TRIM(ctrcnm(jtra)), trcsdta(nis0:nis1,:,nt_x,jtra), &
                 &               ktime=ntobc_x , kstart=istart, kcount= icount )
          CALL iom_close (id_s)

          ! mask the boundary values
            trcsdta(:,:,nt_x,jtra) = trcsdta(:,:,nt_x,jtra)*tsmsk(:,:)

          ! check any outliers
!          zmin=MINVAL( trsdta(:,:,nt_x,jn), mask=ltsmsk ) ; zmax=MAXVAL(trsdta(:,:,nt_x,jn), mask=ltsmsk)
!          IF (  zmin < -10. .OR. zmax > 40)   THEN
!             CALL ctl_stop('Error in trsdta',' routine obcdta')
!          ENDIF

          IF ( lwp .AND.  ( kt == nit000 .OR. nprint /= 0 )  ) THEN
             WRITE(numout,*)
             WRITE(numout,*) ' Read South OBC data records ', ntobc_x
             ikprint = jpi/20 +1
             WRITE(numout,*) ' Passive tracers  record 1 - printout every 3 level'
             CALL prihre( trcsdta(:,:,nt_x,:), jpi, jpk, 1, jpi, ikprint, jpk, 1, -3, 1., numout )
             WRITE(numout,*)
          ENDIF
       ENDIF
    ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  END SUBROUTINE obc_read_trc_single
  INTEGER FUNCTION nrecbef_trc_single(jtr)
      !!-----------------------------------------------------------------------
      !!                     ***    FUNCTION nrecbef_trc   ***
      !!
      !!  Purpose : - provide the before record number in files, with respect to zjcnes
      !!
      !!    History : 2008-04 : ( J.M. Molines ) Original code
      !!-----------------------------------------------------------------------

      INTEGER, INTENT( in ) ::   jtr
      INTEGER :: it , idum

    idum = itobc
    DO it =1, itobc
       IF ( ztcobc(it,jtr) > zjcnes ) THEN ;  idum = it - 1 ; EXIT ;  ENDIF
    ENDDO
    ! idum can be 0 (climato, before first record)
    IF ( idum == 0 ) THEN
       IF ( ln_obc_clim_trc(jtr) ) THEN
         idum = itobc
       ELSE
         ctmp1='obc_dta_trc: find ntobc == 0 for  non climatological file '
         ctmp2='consider adding a first record in your data file '
         CALL ctl_stop(ctmp1, ctmp2)
       ENDIF
    ENDIF
    ! idum can be itobc ( zjcnes > ztcobc (itobc) )
    !  This is not a problem ...
    nrecbef_trc_single = idum

  END FUNCTION nrecbef_trc_single

#else
  !!------------------------------------------------------------------------------
  !!   default option:           Dummy module          NO Open Boundary Conditions
  !!------------------------------------------------------------------------------
CONTAINS
  SUBROUTINE obc_dta_trc( kt )             ! Dummy routine
    INTEGER, INTENT (in) :: kt
    WRITE(*,*) 'obc_dta_trc: You should not have seen this print! error?', kt
  END SUBROUTINE obc_dta_trc
#endif
END MODULE obcdta_trc

