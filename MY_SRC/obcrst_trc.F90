MODULE obcrst_trc
#if defined key_obc_trc
   !!=================================================================================
   !!                       ***  MODULE  obcrst_trc  ***
   !! Ocean dynamic :  Input/Output files for restart on OBC passive tracers
   !!=================================================================================

   !!---------------------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers variables
   USE dom_oce         ! ocean space and time domain variables
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE phycst          ! physical constants
   USE obc_oce         ! ocean open boundary conditions
   USE lib_mpp         ! for mppobc
   USE in_out_manager  ! I/O manager

   IMPLICIT NONE
   PRIVATE

   !! * Accessibility
   PUBLIC obc_rst_read_trc       ! routine called by obc_ini
   PUBLIC obc_rst_write_trc      ! routine called by step

   !!---------------------------------------------------------------------------------
   !!   OPA 9.0 , LOCEAN-IPSL (2005) 
   !! $Id: obcrst.F90 1152 2008-06-26 14:11:13Z rblod $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!---------------------------------------------------------------------------------

CONTAINS

   SUBROUTINE obc_rst_write_trc ( kt )
      !!--------------------------------------------------------------------------------
      !!                  ***  SUBROUTINE obc_rst_write_trc  ***
      !!                
      !! ** Purpose :   Write open boundary restart fields in restart.obc.output file 
      !!
      !! ** Method  :   restart.obc.output file: Direct access non formatted file.
      !!      Each nstock time step , save fields which are necessary for restart.
      !!      - This routine is called if at least the key_obc is defined. It is called
      !!        at the same time step than rstwrite.
      !!      - First record holds OBC parameters nbobc,jpieob,jpiwob,jpjnob,jpjsob and 
      !!        the OBC layout jpjed, jpjef ... for checking purposes.
      !!      - Following records hold the boundary arrays, in the order east west north
      !!        south, if they exist.
      !!      - The writing is realised by vertical slab across the boundary, for bsf, u,
      !!        v, t, and s boundary arrays. Each record hold a vertical slab.
      !!      - For mpp, this allows each processor to write only the correct informations
      !!        it hold. If a processor has no valid informations on boundary, it just 
      !!        skip the writing part (automatically).
      !!      - Special care is taken for dumping the starting point of a boundary (jpjed,
      !!        jpjwd, jpind, jpisd) because of the general definition of nje0 njw0,nin0,
      !!        nis0. This is done to avoid records to be written by 2 adjacent processors.
      !!
      !!  History :
      !!         ! 97-11 (J.M. Molines) Original code
      !!         ! 98-11 (J.M. Molines) Bug fix for adjacent processors
      !!   8.5   ! 02-10 (C. Talandier, A-M. Treguier) F90
      !!         ! 03-06 (J.M. Molines) Bug fix for adjacent processors
      !!   9.0   ! 04-02 (G. Madec)  suppression of numwob, use inum
      !!         ! 09-11 (C. Dufour) adapted for passive tracers
      !!-----------------------------------------------------------------------------------
      !! * Arguments
      INTEGER, INTENT( in ) ::   kt   ! ocean time-step index

      !! * Local declarations
      INTEGER ::   ji, jj, jk, jn
      INTEGER ::   inum               ! temporary logical unit
      INTEGER ::   ibloc, nreclo, jrec, jt, jb 
      INTEGER ::   jfoe, jfow, ifon, ifos
      INTEGER ::   ino0, it0
      !!-----------------------------------------------------------------------------

      ! 1. Output of restart fields (inum)
      ! ------------------------------------
 
      IF( ( mod(kt,nstock) == 0 ) .OR. ( kt == nitend ) ) THEN

         ! 1.0 Initializations
         ! -------------------
         IF(lwp) THEN
              WRITE(numout,*) ' '
              WRITE(numout,*) 'obcrst_trc: OBC output for restart with obc_rst_write_trc routine'
              WRITE(numout,*) '~~~~~~'
              WRITE(numout,*) '        output done in restart.obc.output file at it= ', kt, ' date= ', ndastp
         END IF

         ino0 = no
         it0  = kt
         ibloc  = 4096*4
         nreclo = ibloc*( ( ( 4*jptra *jpk )*jpbyt -1)/ibloc + 1)
         IF(lwp) WRITE(numout,*) '             '
         IF(lwp) WRITE(numout,*) '        OBC restart file opened with nreclo = ',nreclo

         ! 1.1 Open file
         ! -------------

         CALL ctl_opn( inum, 'restart.obc.trc.output', 'UNKNOWN', 'UNFORMATTED', 'DIRECT',   &
            &         nreclo, numout, lwp, 1 )
 
         ! 1.2 Write header
         ! ----------------
         WRITE (inum,REC=1) ino0,it0,nbobc,jpieob,jpiwob,jpjnob,jpjsob,     &
                              jpjed,jpjef,jpjwd,jpjwf,jpind,jpinf,jpisd,jpisf

         ! 1.3 Write east boundary array if any.
         ! -------------------------------------
         IF( lp_obc_east_trc ) THEN
            IF( lfbceast_trc ) THEN
               IF(lwp) THEN
                  WRITE(numout,*) ' '
                  WRITE(numout,*) '        No restart file for the fixed east OBC for tracers'
               END IF
            ELSE
               IF( jpieob /= 0 ) THEN
                  IF( nje0+njmpp-1  == jpjed .AND. nie1 >= nie0 ) THEN
            ! ... dump of jpjed if it is on this proc.
                     jrec = 2
                     jfoe = jpjed - njmpp + 1
                     PRINT *,'Narea =',narea,' write jrec =2 east'
                     WRITE(inum,REC=jrec)                                    &
                        (((( trcebnd(jfoe,jk,jb,jt,jn),jk=1,jpk),jb=1,2),jt=1,2),jn=1,jptra)
                  ENDIF
                  DO ji = nie0, nie1
                     DO jj = nje0, nje1
            ! ... only interested processors go through the following lines
            !           jfoe = jj + njmpp -1
                        jfoe = jj 
                        jrec = 2 + jj + njmpp -1 -jpjed
                        WRITE (inum,REC=jrec)                                   &
                          (((( trcebnd(jfoe,jk,jb,jt,jn),jk=1,jpk),jb=1,2),jt=1,2),jn=1,jptra) 
                     END DO
                  END DO
               END IF
            END IF
         END IF
 
         ! 1.4 Write west boundary arrays if any
         ! -------------------------------------
         IF( lp_obc_west_trc ) THEN
            IF( lfbcwest_trc ) THEN
               IF(lwp) THEN
                  WRITE(numout,*) ' '
                  WRITE(numout,*) '        No restart file for the fixed west OBC for tracers'
               END IF
            ELSE
               IF( jpiwob /= 0 ) THEN
                  IF( njw0+njmpp+1 == jpjwd .AND. niw1 >= niw0 ) THEN
            ! ... dump of jpjwd if it is on this proc.
                     jrec = 3 + jpjef - jpjed
            !        jfow = jpjwd
                     jfow = jpjwd -njmpp + 1
                     PRINT *,'Narea =',narea,' write jrec =',jrec,' west'
                     WRITE (inum,REC=jrec)                                   &
                       (((( trcwbnd(jfow,jk,jb,jt,jn),jk=1,jpk),jb=1,2),jt=1,2),jn=1,jptra)
                  END IF
                  DO ji = niw0, niw1
                     DO jj = njw0, njw1
            ! ... only interested processors go through the following lines
            !           jfow = jj + njmpp -1
                        jfow = jj 
                        jrec = 3 + jpjef -jpjed + jj + njmpp -1 -jpjwd
                        WRITE (inum,REC=jrec)                                   &
                          (((( trcwbnd(jfow,jk,jb,jt,jn),jk=1,jpk),jb=1,2),jt=1,2),jn=1,jptra)
                     END DO
                  END DO
               END IF
            END IF
         END IF
 
         ! 1.5 Write north boundary arrays if any
         ! --------------------------------------
         IF( lp_obc_north_trc ) THEN
            IF( lfbcnorth_trc ) THEN
               IF(lwp) THEN
                  WRITE(numout,*) ' '
                  WRITE(numout,*) '        No restart file for the fixed north OBC for tracers'
               END IF
            ELSE
               IF( jpjnob /= 0) THEN
                  IF( nin0+nimpp-1 == jpind .AND. njn1 >= njn0 ) THEN
            ! ... dump of jpind if it is on this proc.
                     jrec = 4 + jpjef -jpjed + jpjwf -jpjwd
            !        ifon = jpind
                     ifon = jpind -nimpp +1
                     WRITE (inum,REC=jrec)                                   &
                       (((( trcnbnd(ifon,jk,jb,jt,jn),jk=1,jpk),jb=1,2),jt=1,2),jn=1,jptra)
                  END IF
                  DO jj = njn0, njn1
                     DO ji = nin0, nin1
            ! ... only interested processors go through the following lines
            !           ifon = ji + nimpp -1
                        ifon = ji 
                        jrec = 4 + jpjef -jpjed + jpjwf -jpjwd +ji + nimpp -1  -jpind
                        WRITE (inum,REC=jrec)                                   &
                          (((( trcnbnd(ifon,jk,jb,jt,jn),jk=1,jpk),jb=1,2),jt=1,2),jn=1,jptra)
                     END DO
                  END DO
               END IF
            END IF
         END IF
 
         ! 1.6 Write south boundary arrays if any
         ! --------------------------------------
         IF( lp_obc_south_trc ) THEN
            IF( lfbcsouth_trc ) THEN
               IF(lwp) THEN
                  WRITE(numout,*) ' '
                  WRITE(numout,*) '        No restart file for the fixed south OBC for tracers'
                  WRITE(numout,*) ' '
               END IF
            ELSE
               IF( jpjsob /= 0 ) THEN
                  IF( nis0+nimpp-1 == jpisd .AND. njs1 >= njs0 ) THEN
            ! ... dump of jpisd if it is on this proc.
                     jrec = 5 + jpjef -jpjed + jpjwf -jpjwd +jpinf -jpind
            !        ifos = jpisd
                     ifos = jpisd -nimpp + 1
                     WRITE (inum,REC=jrec)                                   &
                       (((( trcsbnd(ifos,jk,jb,jt,jn),jk=1,jpk),jb=1,2),jt=1,2),jn=1,jptra)
                  END IF
                  DO jj = njs0, njs1
                     DO ji = nis0, nis1
            ! ... only interested processors go through the following lines
            !           ifos = ji + nimpp -1
                        ifos = ji 
                        jrec = 5 + jpjef -jpjed + jpjwf -jpjwd +jpinf -jpind + &
                              ji + nimpp -1 -jpisd
                        WRITE (inum,REC=jrec) &
                          (((( trcsbnd(ifos,jk,jb,jt,jn),jk=1,jpk),jb=1,2),jt=1,2),jn=1,jptra)
                     END DO
                  END DO
               END IF
            END IF
         END IF
      CLOSE(inum)
      END IF

   END SUBROUTINE obc_rst_write_trc


   SUBROUTINE obc_rst_read_trc
      !!----------------------------------------------------------------------------
      !!                   ***  SUBROUTINE obc_rst_read_trc  ***
      !!                   
      !! ** Purpose :   Read files for restart at open boundaries
      !!
      !! ** Method  :   Read the previous boundary arrays on unit inum
      !!      The first record indicates previous characterics
      !!
      !! History :
      !!        ! 97-11 (J.M. Molines) Original code
      !!   8.5  ! 02-10 (C. Talandier, A-M. Treguier) F90
      !!        ! 09-11 (C. Dufour) adapted for passive tracers
      !!----------------------------------------------------------------------------
      !! * Local declarations
      INTEGER ::   inum = 11            ! temporary logical unit
      INTEGER ::   ji,jj,jk,jn,ios
      INTEGER ::   ino0,it0,nbobc0,jpieob0,jpiwob0,jpjnob0,jpjsob0
      INTEGER ::   ino1,it1,nbobc1,jpieob1,jpiwob1,jpjnob1,jpjsob1
      INTEGER ::   ied0,ief0,iwd0,iwf0,ind0,inf0,isd0,isf0
      INTEGER ::   ied1,ief1,iwd1,iwf1,ind1,inf1,isd1,isf1
      INTEGER ::   ibloc, nreclo, jrec, jt, jb
      INTEGER ::   jfoe, jfow, ifon, ifos
      !!-----------------------------------------------------------------------------

      ! 0. Initialisations
      ! ------------------
 
      ino0    = no
      it0     = nit000
      nbobc0  = nbobc
      jpieob0 = jpieob
      jpiwob0 = jpiwob
      jpjnob0 = jpjnob
      jpjsob0 = jpjsob
 
      ied0   = jpjed
      ief0   = jpjef
      iwd0   = jpjwd
      iwf0   = jpjwf
      ind0   = jpind
      inf0   = jpinf
      isd0   = jpisd
      isf0   = jpisf
 
      ibloc  = 4096*4
      nreclo = ibloc *( ( ( 4 * jptra *jpk )*jpbyt -1)/ibloc + 1)
 
      IF(lwp) THEN
         WRITE(numout,*) 'obcrst_trc: beginning of restart with obc_rst_read_trc routine'
         WRITE(numout,*) '~~~~~~'
         WRITE(numout,*) ' '
         WRITE(numout,*) '        The present run :'
         WRITE(numout,*) '        number job is  : ',no 
         WRITE(numout,*) '        with the time nit000 : ',nit000
         WRITE(numout,*) '        OBC restart file opened with nreclo = ',nreclo 
      END IF
 
      ! 0.1 Open files
      ! ---------------
      CALL ctl_opn( inum, 'restart.obc.trc', 'UNKNOWN', 'UNFORMATTED', 'DIRECT',   &
         &         nreclo, numout, lwp, 1 )

      ! 1. Read
      ! -------
 
      ! 1.1 First record
      ! -----------------
      READ(inum,REC=1) ino1,it1,nbobc1,jpieob1,jpiwob1,jpjnob1,     &
                         jpjsob1,ied1,ief1,iwd1,iwf1,ind1,inf1,isd1,isf1
 
      IF(lwp) THEN
         WRITE(numout,*) ' '
         WRITE(numout,*) '        READ inum with number job : ',ino1,' with the time it: ',it1
         WRITE(numout,*) ' '
      END IF
 
      ! 1.2 Control of date
      ! --------------------
      IF( ( it0-it1 ) /= 1 .AND. abs(nrstdt) == 1 ) THEN
          CALL ctl_stop( '        ===>>>> : problem with nit000 for the restart',   &
               &         '        ==============',   &
               &         '        we stop in obc_rst_read_trc routine. Verify the file or rerun with the value',   &
               &         '        0 for the control of time parameter nrstdt' )
             
      END IF
 
      ! 1.3 Control of number of open boundaries
      ! ----------------------------------------
      IF( nbobc1 /= nbobc0 ) THEN
         IF(lwp) THEN
            WRITE(numout,*) '        ===> W A R N I N G: The number of OBC have changed:'
            WRITE(numout,*) '        Last run : ',nbobc0,' obcs'
            WRITE(numout,*) '        This run : ',nbobc1,' obcs'
         END IF
      END IF
 
      ! 1.4 Control of which boundary is open
      ! -------------------------------------
      IF( lp_obc_east_trc .AND. ( jpieob1 /= 0 ) ) THEN
         IF(lwp) THEN
            WRITE(numout,*) '         '
            WRITE(numout,*) '        East open boundary'
            IF( jpieob0 /= jpieob1 ) CALL ctl_stop( '         ==>>>> : Problem in obc_rst_read_trc, jpieob have changed' )
         END IF
      END IF
 
      IF( lp_obc_west_trc .AND. ( jpiwob1 /= 0 ) ) THEN
         IF(lwp) THEN
            WRITE(numout,*) '         '
            WRITE(numout,*) '        West open boundary'
            IF( jpiwob0 /= jpiwob1 ) CALL ctl_stop( '        ==>>>> : Problem in obc_rst_read_trc, jpiwob has changed' )
         END IF
      END IF
 
      IF( lp_obc_north_trc .AND. ( jpjnob1 /= 0 ) ) THEN
         IF(lwp) THEN
            WRITE(numout,*) '         '
            WRITE(numout,*) '        North open boundary'
            IF( jpjnob0 /= jpjnob1 ) CALL ctl_stop( '        ==>>>> : Problem in obc_rst_read_trc, jpjnob has changed' )
         END IF
      END IF
 
      IF( lp_obc_south_trc .AND. ( jpjsob1 /= 0 ) ) THEN
         IF(lwp) THEN
            WRITE(numout,*) '         '
            WRITE(numout,*) '        South open boundary'
            IF( jpjsob0 /= jpjsob1) CALL ctl_stop( '        ==>>>> : Problem in obc_rst_read_trc, jpjsob has changed' )
         END IF
      END IF
 
 
      ! 1.5 Control of the limit of the boundaries
      ! ------------------------------------------
      IF( lp_obc_east_trc .AND. ( jpieob1 /= 0 ) ) THEN
         IF( ied1 /= ied0 ) CALL ctl_stop( '        ==>>>> : Problem in obc_rst_read_trc, jpjed has changed' )
         IF( ief1 /= ief0 ) CALL ctl_stop( '        ==>>>> : Problem in obc_rst_read_trc, jpjef has changed' )
      END IF

      IF( lp_obc_west_trc .AND. ( jpiwob1 /= 0 ) ) THEN
         IF( iwd1 /= iwd0 ) CALL ctl_stop( '        ==>>>> : Problem in obc_rst_read_trc, jpjwd has changed' )
         IF( iwf1 /= iwf0 ) CALL ctl_stop( '        ==>>>> : Problem in obc_rst_read_trc, jpjwf has changed' )
      END IF
 
      IF( lp_obc_north_trc .AND. ( jpjnob1 /= 0 ) ) THEN
         IF( ind1 /= ind0 ) CALL ctl_stop( '        ==>>>> : Problem in obc_rst_read_trc, jpind has changed' )
         IF( inf1 /= inf0 ) CALL ctl_stop( '        ==>>>> : Problem in obc_rst_read_trc, jpinf has changed' )
      END IF
 
      IF( lp_obc_south_trc .AND. ( jpjsob1 /= 0 ) ) THEN
         IF( isd1 /= isd0 ) CALL ctl_stop( '        ==>>>> : Problem in obc_rst_read_trc, jpisd has changed' )
         IF( isf1 /= isf0 ) CALL ctl_stop( '        ==>>>> : Problem in obc_rst_read_trc, jpisf has changed' )
      END IF
 
 
      ! 2. Now read the boundary arrays
      ! -------------------------------
 
      ! 2.1 Read east boundary array if any.
      ! ------------------------------------
      IF( lp_obc_east_trc ) THEN
         IF( jpieob1 /= 0) THEN
            IF( nje0+njmpp-1 == jpjed .AND. nie1 >= nie0 ) THEN
      ! ... read of jpjed if it is on this proc.
               jrec = 2
      !        jfoe = jpjed
               jfoe = jpjed -njmpp + 1
               READ (inum,REC=jrec)                                   &
                 (((( trcebnd(jfoe,jk,jb,jt,jn),jk=1,jpk),jb=1,2),jt=1,2),jn=1,jptra)
            END IF
            DO ji = nie0, nie1
               DO jj = nje0, nje1
      ! ... only interested processors go through the following lines
      !           jfoe = jj + njmpp -1
                  jfoe = jj 
                  jrec = 2 + jj + njmpp -1 -jpjed
                  READ (inum,REC=jrec)                                   &
                    (((( trcebnd(jfoe,jk,jb,jt,jn),jk=1,jpk),jb=1,2),jt=1,2),jn=1,jptra)
               END DO
            END DO

         ELSE

            !  lp_obc_east was not TRUE previously
         END IF

      END IF
 
      ! 2.2 Read west boundary arrays if any.
      ! -------------------------------------
      IF( lp_obc_west_trc ) THEN
         IF( jpiwob1 /= 0) THEN
            IF( njw0+njmpp-1 == jpjwd .AND. niw1 >= niw0 ) THEN
      ! ... read of jpjwd if it is on this proc.
               jrec = 3 + jpjef - jpjed
      !        jfow = jpjwd
               jfow = jpjwd -njmpp + 1
               READ (inum,REC=jrec)                                   &
                 (((( trcwbnd(jfow,jk,jb,jt,jn),jk=1,jpk),jb=1,2),jt=1,2),jn=1,jptra)
            END IF
            DO ji = niw0, niw1
               DO jj = njw0, njw1
      ! ... only interested processors go through the following lines
      !           jfow = jj + njmpp -1
                  jfow = jj 
                  jrec = 3 + jpjef -jpjed + jj + njmpp -1 -jpjwd
                  READ (inum,REC=jrec)                                   &
                    (((( trcwbnd(jfow,jk,jb,jt,jn),jk=1,jpk),jb=1,2),jt=1,2),jn=1,jptra)
               END DO
            END DO

         ELSE

            !  lp_obc_west was not TRUE previously
         END IF

      END IF
 
      ! 2.3 Read north boundary arrays if any.
      ! --------------------------------------
      IF( lp_obc_north_trc ) THEN
         IF( jpjnob1 /= 0) THEN
            IF( nin0+nimpp-1 == jpind .AND. njn1 >= njn0 ) THEN
      ! ... read of jpind if it is on this proc.
               jrec = 4 + jpjef -jpjed + jpjwf -jpjwd
      !        ifon = jpind
               ifon = jpind -nimpp +1
               READ (inum,REC=jrec)                                   &
                 (((( trcnbnd(ifon,jk,jb,jt,jn),jk=1,jpk),jb=1,2),jt=1,2),jn=1,jptra)
            END IF
            DO jj = njn0, njn1
               DO ji = nin0, nin1
      ! ... only interested processors go through the following lines
      !           ifon = ji + nimpp -1
                  ifon = ji 
                  jrec = 4 + jpjef -jpjed + jpjwf -jpjwd +ji + nimpp -1  -jpind
                  READ (inum,REC=jrec)                                   & 
                    (((( trcnbnd(ifon,jk,jb,jt,jn),jk=1,jpk),jb=1,2),jt=1,2),jn=1,jptra)
               END DO
            END DO

         ELSE

           !  lp_obc_north was not TRUE previously
         END IF

      END IF
 
      ! 2.4 Read south boundary arrays if any.
      ! -------------------------------------
      IF( lp_obc_south_trc ) THEN
         IF( jpjsob1 /= 0) THEN
            IF( nis0+nimpp-1 == jpisd .AND. njs1 >= njs0 ) THEN
      ! ... read of jpisd if it is on this proc.
               jrec = 5 + jpjef -jpjed + jpjwf -jpjwd +jpinf -jpind
      !        ifos = jpisd
               ifos = jpisd -nimpp + 1
               READ (inum,REC=jrec)                                   &
                 (((( trcsbnd(ifos,jk,jb,jt,jn),jk=1,jpk),jb=1,2),jt=1,2),jn=1,jptra)
            END IF
            DO jj = njs0, njs1
               DO ji = nis0, nis1
      ! ... only interested processors go through the following lines
      !           ifos = ji + nimpp -1
                  ifos = ji 
                  jrec = 5 + jpjef -jpjed + jpjwf -jpjwd +jpinf -jpind +  &
                        ji + nimpp -1 -jpisd
                  READ (inum,REC=jrec)                                   & 
                    (((( trcsbnd(ifos,jk,jb,jt,jn),jk=1,jpk),jb=1,2),jt=1,2),jn=1,jptra)
               END DO
            END DO
         ELSE
            !  lp_obc_south was not TRUE previously
         END IF

      END IF
      CLOSE(inum)

      IF( lk_mpp ) THEN
         IF( lp_obc_east_trc ) THEN
            CALL mppobc(trcebnd,jpjed,jpjef,jpieob+1,jpk*2*2*jptra,2,jpj, numout )
         ENDIF
         IF( lp_obc_west_trc ) THEN
            CALL mppobc(trcwbnd,jpjwd,jpjwf,jpiwob,jpk*2*2*jptra,2,jpj, numout )
         ENDIF
         IF( lp_obc_north_trc ) THEN 
            CALL mppobc(trcnbnd,jpind,jpinf,jpjnob+1,jpk*2*2*jptra,1,jpi, numout )
         ENDIF
         IF( lp_obc_south_trc ) THEN
            CALL mppobc(trcsbnd,jpisd,jpisf,jpjsob,jpk*2*2*jptra,1,jpi, numout )
         ENDIF
      ENDIF
 
   END SUBROUTINE obc_rst_read_trc
#else
   !!=================================================================================
   !!                       ***  MODULE  obcrst_trc  ***
   !! Ocean dynamic :  Input/Output files for restart on OBC passive tracers
   !!=================================================================================
CONTAINS
   SUBROUTINE obc_rst_write_trc( kt )           !  No Open boundary ==> empty routine
      INTEGER,INTENT(in) :: kt
      WRITE(*,*) 'obc_rst_write_trc: You should not have seen this print! error?', kt
   END SUBROUTINE obc_rst_write_trc
   SUBROUTINE obc_rst_read_trc                 !  No Open boundary ==> empty routine
   END SUBROUTINE obc_rst_read_trc
#endif

   !!=================================================================================
END MODULE obcrst_trc
