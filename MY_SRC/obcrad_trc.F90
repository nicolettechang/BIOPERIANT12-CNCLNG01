MODULE obcrad_trc 
   !!=================================================================================
   !!                       ***  MODULE  obcrad_trc  ***
   !! Ocean dynamic :   Phase velocities for each open boundary
   !!=================================================================================
#if defined key_obc_trc
   !!---------------------------------------------------------------------------------
   !!   obc_rad_trc        : call the subroutine for each open boundary
   !!   obc_rad_east_trc   : compute the east phase velocities
   !!   obc_rad_west_trc   : compute the west phase velocities
   !!   obc_rad_north_trc  : compute the north phase velocities
   !!   obc_rad_south_trc  : compute the south phase velocities
   !!---------------------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers variables
   USE trc             ! passive tracers
   USE dom_oce         ! ocean space and time domain variables
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE phycst          ! physical constants
   USE obc_oce         ! ocean open boundary conditions
   USE lib_mpp         ! for mppobc
   USE in_out_manager  ! I/O units

   IMPLICIT NONE
   PRIVATE

   !! * Accessibility
   PUBLIC obc_rad_trc    ! routine called by step.F90

   !! * Module variables
   INTEGER ::   ji, jj, jk, jn     ! dummy loop indices

   INTEGER ::      & ! ... boundary space indices 
      nib   = 1,   & ! nib   = boundary point
      nibm  = 2,   & ! nibm  = 1st interior point
      nibm2 = 3,   & ! nibm2 = 2nd interior point
                     ! ... boundary time indices 
      nit   = 1,   & ! nit    = now
      nitm  = 2,   & ! nitm   = before
      nitm2 = 3      ! nitm2  = before-before

   !! * Substitutions
#  include "obc_vectopt_loop_substitute.h90"
   !!---------------------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: obcrad.F90 2977 2011-10-22 13:46:41Z cetlod $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!---------------------------------------------------------------------------------

CONTAINS

   SUBROUTINE obc_rad_trc ( kt )
      !!------------------------------------------------------------------------------
      !!                     SUBROUTINE obc_rad_trc
      !!                    ************************
      !! ** Purpose :
      !!      Perform swap of arrays to calculate radiative phase speeds at the open 
      !!      boundaries and calculate those phase speeds if the open boundaries are 
      !!      not fixed. In case of fixed open boundaries does nothing.
      !!
      !!     The logical variable lp_obc_east, and/or lp_obc_west, and/or lp_obc_north,
      !!     and/or lp_obc_south allow the user to determine which boundary is an
      !!     open one (must be done in the param_obc.h90 file).
      !! 
      !! ** Reference : 
      !!     Marchesiello P., 1995, these de l'universite J. Fourier, Grenoble, France.
      !!
      !!  History :
      !!    9.0  !  09-11  (C. Dufour) adapted for passive tracers
      !!------------------------------------------------------------------------------
      INTEGER, INTENT( in ) ::   kt
      !!----------------------------------------------------------------------

      IF( lp_obc_east_trc  .AND. .NOT.lfbceast_trc  )   CALL obc_rad_east_trc ( kt )   ! East open boundary

      IF( lp_obc_west_trc  .AND. .NOT.lfbcwest_trc  )   CALL obc_rad_west_trc ( kt )   ! West open boundary

      IF( lp_obc_north_trc .AND. .NOT.lfbcnorth_trc )   CALL obc_rad_north_trc( kt )   ! North open boundary

      IF( lp_obc_south_trc .AND. .NOT.lfbcsouth_trc )   CALL obc_rad_south_trc( kt )   ! South open boundary

   END SUBROUTINE obc_rad_trc


   SUBROUTINE obc_rad_east_trc ( kt )
      !!------------------------------------------------------------------------------
      !!                     ***  SUBROUTINE obc_rad_east_trc  ***
      !!                   
      !! ** Purpose :
      !!      Perform swap of arrays to calculate radiative phase speeds at the open 
      !!      east boundary and calculate those phase speeds if this OBC is not fixed.
      !!      In case of fixed OBC, this subrountine is not called.
      !!
      !!  History :
      !!         ! 95-03 (J.-M. Molines) Original from SPEM
      !!         ! 97-07 (G. Madec, J.-M. Molines) additions
      !!         ! 97-12 (M. Imbard) Mpp adaptation
      !!         ! 00-06 (J.-M. Molines) 
      !!    8.5  ! 02-10 (C. Talandier, A-M. Treguier) Free surface, F90
      !!------------------------------------------------------------------------------
      !! * Arguments
      INTEGER, INTENT( in ) ::   kt

      !! * Local declarations
      INTEGER  ::   ij
      !!------------------------------------------------------------------------------

      ! 1. Swap arrays before calculating radiative velocities
      ! ------------------------------------------------------

         ! Passive tracers
         ! ----------------------------

      IF( kt > nit000 ) THEN

         ! ... advance in time (time filter, array swap)
         DO jn = 1, jptra
           DO jk = 1, jpkm1
            DO jj = 1, jpj
         ! ... fields nitm <== nit  plus time filter at the boundary
               trcebnd(jj,jk,nib,nitm,jn) = trcebnd(jj,jk,nib,nit,jn)*temsk(jj,jk)
            END DO
           END DO
         END DO

         DO jn = 1, jptra
          DO ji = fs_nie0+1, fs_nie1+1 ! Vector opt.
            DO jk = 1, jpkm1
               DO jj = 1, jpj
                  trcebnd(jj,jk,nibm,nitm,jn) = trcebnd(jj,jk,nibm,nit,jn)*temsk(jj,jk)
         ! ... fields nit <== now (kt+1)
                  trcebnd(jj,jk,nib  ,nit,jn) = trn(ji  ,jj,jk,jn)*temsk(jj,jk)
                  trcebnd(jj,jk,nibm ,nit,jn) = trn(ji-1,jj,jk,jn)*temsk(jj,jk)
               END DO
            END DO
           END DO
         END DO
         IF( lk_mpp )   CALL mppobc(trcebnd,jpjed,jpjef,jpieob+1,jpk*2*2*jptra,2,jpj, numout )

         ! ... extremeties nie0, nie1
         ij = jpjed +1 - njmpp
         IF( ij >= 2 .AND. ij < jpjm1 ) THEN
           DO jn = 1, jptra
            DO jk = 1,jpkm1
               trcebnd(ij,jk,nibm,nitm,jn) = trcebnd(ij+1 ,jk,nibm,nitm,jn)
            END DO
           END DO
         END IF
         ij = jpjef +1 - njmpp
         IF( ij >= 2 .AND. ij < jpjm1 ) THEN
           DO jn = 1, jptra
            DO jk = 1,jpkm1
               trcebnd(ij,jk,nibm,nitm,jn) = trcebnd(ij-1 ,jk,nibm,nitm,jn)
            END DO
           END DO
         END IF

      END IF     ! End of array swap


   END SUBROUTINE obc_rad_east_trc


   SUBROUTINE obc_rad_west_trc ( kt )
      !!------------------------------------------------------------------------------
      !!                  ***  SUBROUTINE obc_rad_west_trc  ***
      !!                    
      !! ** Purpose :
      !!      Perform swap of arrays to calculate radiative phase speeds at the open 
      !!      west boundary and calculate those phase speeds if this OBC is not fixed.
      !!      In case of fixed OBC, this subrountine is not called.
      !!
      !!  History :
      !!         ! 95-03 (J.-M. Molines) Original from SPEM
      !!         ! 97-07 (G. Madec, J.-M. Molines) additions
      !!         ! 97-12 (M. Imbard) Mpp adaptation
      !!         ! 00-06 (J.-M. Molines) 
      !!    8.5  ! 02-10 (C. Talandier, A-M. Treguier) Free surface, F90
      !!------------------------------------------------------------------------------
      !! * Arguments
      INTEGER, INTENT( in ) ::   kt

      !! * Local declarations
      INTEGER ::   ij
      !!------------------------------------------------------------------------------

      ! 1. Swap arrays before calculating radiative velocities
      ! ------------------------------------------------------

         ! Passive tracers
         ! ----------------------------
 
      IF( kt > nit000 ) THEN
         ! ... advance in time (time filter, array swap)
         DO jn = 1, jptra
          DO jk = 1, jpkm1
            DO jj = 1, jpj
         ! ... fields nitm <== nit  plus time filter at the boundary
               trcwbnd(jj,jk,nib,nitm,jn) = trcwbnd(jj,jk,nib,nit,jn)*twmsk(jj,jk)
            END DO
          END DO
         END DO 

         DO jn = 1, jptra
          DO ji = fs_niw0, fs_niw1 ! Vector opt.
            DO jk = 1, jpkm1
               DO jj = 1, jpj
                  trcwbnd(jj,jk,nibm ,nitm,jn) = trcwbnd(jj,jk,nibm ,nit,jn)*twmsk(jj,jk)
         ! ... fields nit <== now (kt+1)
                  trcwbnd(jj,jk,nib  ,nit,jn) = trn(ji   ,jj,jk,jn)*twmsk(jj,jk)
                  trcwbnd(jj,jk,nibm ,nit,jn) = trn(ji+1 ,jj,jk,jn)*twmsk(jj,jk)
               END DO
            END DO
          END DO
         END DO
         IF( lk_mpp )   CALL mppobc(trcwbnd,jpjwd,jpjwf,jpiwob,jpk*2*2*jptra,2,jpj, numout )

         ! ... extremeties niw0, niw1
         ij = jpjwd +1 - njmpp
         IF( ij >= 2 .AND. ij < jpjm1 ) THEN
           DO jn = 1, jptra
            DO jk = 1,jpkm1
               trcwbnd(ij,jk,nibm,nitm,jn) = trcwbnd(ij+1 ,jk,nibm,nitm,jn)
            END DO
           END DO
         END IF
         ij = jpjwf +1 - njmpp
         IF( ij >= 2 .AND. ij < jpjm1 ) THEN
           DO jn = 1, jptra
            DO jk = 1,jpkm1
               trcwbnd(ij,jk,nibm,nitm,jn) = trcwbnd(ij-1 ,jk,nibm,nitm,jn)
            END DO
           END DO
         END IF
 
      END IF     ! End of array swap


   END SUBROUTINE obc_rad_west_trc


   SUBROUTINE obc_rad_north_trc ( kt )
      !!------------------------------------------------------------------------------
      !!                  ***  SUBROUTINE obc_rad_north_trc  ***
      !!           
      !! ** Purpose :
      !!      Perform swap of arrays to calculate radiative phase speeds at the open 
      !!      north boundary and calculate those phase speeds if this OBC is not fixed.
      !!      In case of fixed OBC, this subrountine is not called.
      !!
      !!  History :
      !!         ! 95-03 (J.-M. Molines) Original from SPEM
      !!         ! 97-07 (G. Madec, J.-M. Molines) additions
      !!         ! 97-12 (M. Imbard) Mpp adaptation
      !!         ! 00-06 (J.-M. Molines) 
      !!    8.5  ! 02-10 (C. Talandier, A-M. Treguier) Free surface, F90
      !!------------------------------------------------------------------------------
      !! * Arguments
      INTEGER, INTENT( in ) ::   kt

      !! * Local declarations
      INTEGER  ::   ii
      !!------------------------------------------------------------------------------

      ! 1. Swap arrays before calculating radiative velocities
      ! ------------------------------------------------------

         ! Passive tracers
         ! ----------------------------

      IF( kt > nit000 ) THEN
         ! ... advance in time (time filter, array swap)
         DO jn = 1, jptra
          DO jk = 1, jpkm1
            DO ji = 1, jpi
         ! ... fields nitm <== nit  plus time filter at the boundary
               trcnbnd(ji,jk,nib ,nitm,jn) = trcnbnd(ji,jk,nib,nit,jn)*tnmsk(ji,jk)
            END DO
          END DO
         END DO

         DO jn = 1, jptra
          DO jj = fs_njn0+1, fs_njn1+1  ! Vector opt.
            DO jk = 1, jpkm1
               DO ji = 1, jpi
                  trcnbnd(ji,jk,nibm ,nitm,jn) = trcnbnd(ji,jk,nibm ,nit,jn)*tnmsk(ji,jk)
         ! ... fields nit <== now (kt+1)
                  trcnbnd(ji,jk,nib  ,nit,jn) = trn(ji,jj,  jk,jn)*tnmsk(ji,jk)
                  trcnbnd(ji,jk,nibm ,nit,jn) = trn(ji,jj-1,jk,jn)*tnmsk(ji,jk)
               END DO
            END DO
          END DO
         END DO
         IF( lk_mpp )   CALL mppobc(trcnbnd,jpind,jpinf,jpjnob+1,jpk*2*2*jptra,1,jpi, numout )

         ! ... extremeties  njn0,njn1
         ii = jpind + 1 - nimpp
         IF( ii >= 2 .AND. ii < jpim1 ) THEN
           DO jn = 1, jptra
            DO jk = 1, jpkm1
               trcnbnd(ii,jk,nibm,nitm,jn) = trcnbnd(ii+1,jk,nibm,nitm,jn)
            END DO
           END DO
         END IF
         ii = jpinf + 1 - nimpp
         IF( ii >= 2 .AND. ii < jpim1 ) THEN
           DO jn = 1, jptra
            DO jk = 1, jpkm1
               trcnbnd(ii,jk,nibm,nitm,jn) = trcnbnd(ii-1,jk,nibm,nitm,jn)
            END DO
           END DO
         END IF

      END IF     ! End of array swap


   END SUBROUTINE obc_rad_north_trc


   SUBROUTINE obc_rad_south_trc ( kt )
      !!------------------------------------------------------------------------------
      !!                  ***  SUBROUTINE obc_rad_south_trc  ***
      !!           
      !! ** Purpose :
      !!      Perform swap of arrays to calculate radiative phase speeds at the open 
      !!      south boundary and calculate those phase speeds if this OBC is not fixed.
      !!      In case of fixed OBC, this subrountine is not called.
      !!
      !!  History :
      !!         ! 95-03 (J.-M. Molines) Original from SPEM
      !!         ! 97-07 (G. Madec, J.-M. Molines) additions
      !!         ! 97-12 (M. Imbard) Mpp adaptation
      !!         ! 00-06 (J.-M. Molines) 
      !!    8.5  ! 02-10 (C. Talandier, A-M. Treguier) Free surface, F90
      !!------------------------------------------------------------------------------
      !! * Arguments
      INTEGER, INTENT( in ) ::   kt

      !! * Local declarations
      INTEGER ::   ii
      REAL(wp) ::   z05cx, zdt, z4nor2, z2dx, z2dy
      REAL(wp) ::   zvcb, zvcbm, zvcbm2
      !!------------------------------------------------------------------------------

      ! 1. Swap arrays before calculating radiative velocities
      ! ------------------------------------------------------

         ! Passive tracers
         ! ----------------------------

      IF( kt > nit000 ) THEN
         ! ... advance in time (time filter, array swap)
         DO jn = 1, jptra
          DO jk = 1, jpkm1
            DO ji = 1, jpi
         ! ... fields nitm <== nit  plus time filter at the boundary
               trcsbnd(ji,jk,nib,nitm,jn) = trcsbnd(ji,jk,nib,nit,jn)*tsmsk(ji,jk)
            END DO
          END DO
         END DO

         DO jn = 1, jptra
          DO jj = fs_njs0, fs_njs1  ! Vector opt.
            DO jk = 1, jpkm1
               DO ji = 1, jpi
                  trcsbnd(ji,jk,nibm ,nitm,jn) = trcsbnd(ji,jk,nibm ,nit,jn)*tsmsk(ji,jk)
         ! ... fields nit <== now (kt+1)
                  trcsbnd(ji,jk,nib  ,nit,jn) = trn(ji,jj   ,jk,jn)*tsmsk(ji,jk)
                  trcsbnd(ji,jk,nibm ,nit,jn) = trn(ji,jj+1 ,jk,jn)*tsmsk(ji,jk)
               END DO
            END DO
          END DO
         END DO
         IF( lk_mpp )   CALL mppobc(trcsbnd,jpisd,jpisf,jpjsob,jpk*2*2*jptra,1,jpi, numout )

         ! ... extremeties  njs0,njs1
         ii = jpisd + 1 - nimpp
         IF( ii >= 2 .AND. ii < jpim1 ) THEN
           DO jn = 1, jptra
            DO jk = 1, jpkm1
               trcsbnd(ii,jk,nibm,nitm,jn) = trcsbnd(ii+1,jk,nibm,nitm,jn)
            END DO
           END DO
         END IF
         ii = jpisf + 1 - nimpp
         IF( ii >= 2 .AND. ii < jpim1 ) THEN
           DO jn = 1, jptra
            DO jk = 1, jpkm1
               trcsbnd(ii,jk,nibm,nitm,jn) = trcsbnd(ii-1,jk,nibm,nitm,jn)
            END DO
           END DO
         END IF

      END IF     ! End of array swap

 
   END SUBROUTINE obc_rad_south_trc

#else
   !!=================================================================================
   !!                       ***  MODULE  obcrad_trc  ***
   !! Ocean dynamic :   Phase velocities for each open boundary
   !!=================================================================================
CONTAINS
   SUBROUTINE obc_rad_trc( kt )            ! No open boundaries ==> empty routine
      INTEGER, INTENT(in) :: kt
      WRITE(*,*) 'obc_rad_trc: You should not have seen this print! error?', kt
   END SUBROUTINE obc_rad_trc
#endif

END MODULE obcrad_trc
