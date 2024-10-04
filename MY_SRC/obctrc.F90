MODULE obctrc
   !!=====================================================================================
   !!                       ***  MODULE  obctrc  ***
   !! Ocean tracers:   Radiation of passive tracers on each open boundary
   !! Hist: obctrc.F90 2009-10-05 16:05 inspired from obctra.F90 
   !!=====================================================================================
#if defined key_obc_trc
   !!-------------------------------------------------------------------------------------
   !!   'key_obc_trc'      :                                      Open Boundary Conditions
   !!-------------------------------------------------------------------------------------
   !!   obc_trc        : call the subroutine for each open boundary
   !!   obc_trc_east   : radiation of the east open boundary passive tracers
   !!   obc_trc_west   : radiation of the west open boundary passive tracers
   !!   obc_trc_north  : radiation of the north open boundary passive tracers
   !!   obc_trc_south  : radiation of the south open boundary passive tracers
   !!-------------------------------------------------------------------------------------
   !! * Modules used
   USE oce             ! ocean dynamics and tracers variables
   USE dom_oce         ! ocean space and time domain variables 
   USE phycst          ! physical constants
   USE obc_oce         ! ocean open boundary conditions
   USE lib_mpp         ! ???
   USE lbclnk          ! ???
   USE in_out_manager  ! I/O manager
   USE oce_trc
!  USE sms
   USE trc
   USE trcdta          ! read climatology

   IMPLICIT NONE
   PRIVATE

   !! * Accessibility
   PUBLIC obc_trc     ! routine called in trcnxt.F90 

   !! * Module variables
   INTEGER ::      & ! ... boundary space indices 
      nib   = 1,   & ! nib   = boundary point
      nibm  = 2,   & ! nibm  = 1st interior point
      nibm2 = 3,   & ! nibm2 = 2nd interior point
                     ! ... boundary time indices 
      nit   = 1,   & ! nit    = now
      nitm  = 2,   & ! nitm   = before
      nitm2 = 3      ! nitm2  = before-before

   INTEGER ::  jn

   REAL(wp) ::     &
      rtaue  , rtauw  , rtaun  , rtaus  ,  &  ! Boundary restoring coefficient
      rtauein, rtauwin, rtaunin, rtausin      ! Boundary restoring coefficient for inflow 

   !! * Substitutions
#  include "obc_vectopt_loop_substitute.h90"
   !!---------------------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: obctra.F90 2977 2011-10-22 13:46:41Z cetlod $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!---------------------------------------------------------------------------------

CONTAINS

   SUBROUTINE obc_trc( kt )
      !!-------------------------------------------------------------------------------
      !!                 ***  SUBROUTINE obc_trc  ***
      !!                    
      !! ** Purpose :   Compute passive tracer fields along the open boundaries.
      !!      This routine is called by the trcnxt.F routine and updates tra
      !!      which are the actual passive tracers fields.
      !!        The logical variable lp_obc_east, and/or lp_obc_west, and/or lp_obc_north,
      !!      and/or lp_obc_south allow the user to determine which boundary is an
      !!      open one (must be done in the param_obc.h90 file).
      !!
      !! Reference : 
      !!   Marchesiello P., 1995, these de l'universite J. Fourier, Grenoble, France.
      !!
      !!  History :
      !!        !  09-11 (C. Dufour) wrote it from obctra.F90
      !!----------------------------------------------------------------------
      !! * Arguments
      INTEGER, INTENT( in ) ::   kt
      !!----------------------------------------------------------------------

      ! 0. Local constant initialization

      IF( kt == nit000 .OR. ln_rstart) THEN
         ! ... Boundary restoring coefficient
         rtaue = 2. * rdt / rdpeob_trc
         rtauw = 2. * rdt / rdpwob_trc
         rtaun = 2. * rdt / rdpnob_trc
         rtaus = 2. * rdt / rdpsob_trc
         ! ... Boundary restoring coefficient for inflow ( all boundaries)
         rtauein = 2. * rdt / rdpein_trc
         rtauwin = 2. * rdt / rdpwin_trc
         rtaunin = 2. * rdt / rdpnin_trc
         rtausin = 2. * rdt / rdpsin_trc
      END IF

      IF( lp_obc_east_trc  )   CALL obc_trc_east ( kt )    ! East open boundary

      IF( lp_obc_west_trc  )   CALL obc_trc_west ( kt )    ! West open boundary

      IF( lp_obc_north_trc )   CALL obc_trc_north( kt )    ! North open boundary

      IF( lp_obc_south_trc )   CALL obc_trc_south( kt )    ! South open boundary


      IF( lk_mpp ) THEN                  !!bug ???
       DO jn=1, jptra
         IF( kt >= nit000+3 .AND. ln_rsttr ) THEN
            CALL lbc_lnk( trb(:,:,:,jn), 'T', 1. )
!            CALL lbc_lnk( sb, 'T', 1. )
         END IF
         CALL lbc_lnk( tra(:,:,:,jn), 'T', 1. )
!         CALL lbc_lnk( sa, 'T', 1. )
       ENDDO
      ENDIF

   END SUBROUTINE obc_trc

   SUBROUTINE obc_trc_east ( kt )
      !!------------------------------------------------------------------------------
      !!                ***  SUBROUTINE obc_trc_east  ***
      !!                  
      !! ** Purpose :
      !!      Apply the radiation algorithm on east OBC tracers ta, sa using the 
      !!      phase velocities calculated in obc_rad_east subroutine in obcrad.F90 module
      !!      If the logical lfbceast is .TRUE., there is no radiation but only fixed OBC
      !!
      !!  History :
      !!         ! 95-03 (J.-M. Molines) Original from SPEM
      !!         ! 97-07 (G. Madec, J.-M. Molines) additions
      !!         ! 97-12 (M. Imbard) Mpp adaptation
      !!         ! 00-06 (J.-M. Molines) 
      !!    8.5  ! 02-10 (C. Talandier, A-M. Treguier) F90
      !!------------------------------------------------------------------------------
      !! * Arguments
      INTEGER, INTENT( in ) ::   kt

      !! * Local declaration
      INTEGER ::   ji, jj, jk      ! dummy loop indices
      REAL(wp) ::   z05cx, ztau, zin
      !!------------------------------------------------------------------------------

      ! 1. First three time steps and more if lfbceast is .TRUE.
      !    In that case open boundary conditions are FIXED.
      ! --------------------------------------------------------

      IF( ( kt < nit000+3 .AND. .NOT.ln_rstart ) .OR. lfbceast_trc ) THEN
       DO jn = 1, jptra
         DO ji = fs_nie0+1, fs_nie1+1 ! Vector opt.
            DO jk = 1, jpkm1
               DO jj = 1, jpj
                  tra(ji,jj,jk,jn) = tra(ji,jj,jk,jn) * (1. - temsk(jj,jk)) + &
                                 trcfoe(jj,jk,jn)*temsk(jj,jk)
               END DO
            END DO
         END DO
        END DO

      ELSE

      ! 2. Beyond the fourth time step if lfbceast is .FALSE.
      ! -----------------------------------------------------

         ! Temperature and salinity radiation
         ! ----------------------------------
         !
         !            nibm2      nibm      nib
         !              |   nibm  |   nib///|///
         !              |    |    |    |////|///
         !  jj   line --v----f----v----f----v---
         !              |    |    |    |////|///
         !                   |         |///   //
         !  jj   line   T    u    T    u/// T //
         !                   |         |///   //
         !              |    |    |    |////|///
         !  jj-1 line --v----f----v----f----v---
         !              |    |    |    |////|///
         !                jpieob-1    jpieob / ///
         !              |         |         |
         !           jpieob-1    jpieob     jpieob+1
         !
         ! ... radiative conditions + relaxation toward a climatology
         !     the phase velocity is taken as the phase velocity of the tangen-
         !     tial velocity (here vn), which have been saved in (u_cxebnd,v_cxebnd)
         ! ... (jpjedp1, jpjefm1), jpieob+1
         DO ji = fs_nie0+1, fs_nie1+1 ! Vector opt.
            DO jk = 1, jpkm1
               DO jj = 2, jpjm1
         ! ... i-phase speed ratio (from averaged of v_cxebnd)
                  z05cx = ( 0.5 * ( v_cxebnd(jj,jk) + v_cxebnd(jj-1,jk) ) ) / e1t(ji-1,jj)
                  z05cx = min( z05cx, 1. )
         ! ... z05cx=< 0, inflow  zin=0, ztau=1    
         !           > 0, outflow zin=1, ztau=rtaue
                  zin = sign( 1., z05cx )
                  zin = 0.5*( zin + abs(zin) )
         ! ... for inflow rtauein is used for relaxation coefficient else rtaue
                  ztau = (1.-zin ) * rtauein  + zin * rtaue
                  z05cx = z05cx * zin
         ! ... update ( ta, sa ) with radiative or climatological (t, s)
                 DO jn = 1, jptra
                  tra(ji,jj,jk,jn) = tra(ji,jj,jk,jn) * (1. - temsk(jj,jk)) +           &
                                     temsk(jj,jk) * ( ( 1. - z05cx - ztau )         &
                                     * trcebnd(jj,jk,nib ,nitm,jn) + 2.*z05cx              &
                                     * trcebnd(jj,jk,nibm,nit,jn ) + ztau * trcfoe (jj,jk,jn) ) &
                                     / (1. + z05cx)
                 END DO
               END DO
            END DO
         END DO

      END IF

   END SUBROUTINE obc_trc_east


   SUBROUTINE obc_trc_west ( kt )
      !!------------------------------------------------------------------------------
      !!                 ***  SUBROUTINE obc_trc_west  ***
      !!           
      !! ** Purpose :
      !!      Apply the radiation algorithm on west OBC tracers ta, sa using the 
      !!      phase velocities calculated in obc_rad_west subroutine in obcrad.F90 module
      !!      If the logical lfbcwest is .TRUE., there is no radiation but only fixed OBC
      !!
      !!  History :
      !!         ! 95-03 (J.-M. Molines) Original from SPEM
      !!         ! 97-07 (G. Madec, J.-M. Molines) additions
      !!         ! 97-12 (M. Imbard) Mpp adaptation
      !!         ! 00-06 (J.-M. Molines) 
      !!    8.5  ! 02-10 (C. Talandier, A-M. Treguier) F90
      !!------------------------------------------------------------------------------
      !! * Arguments
      INTEGER, INTENT( in ) ::   kt

      !! * Local declaration
      INTEGER ::   ji, jj, jk      ! dummy loop indices
      REAL(wp) ::   z05cx, ztau, zin
      !!------------------------------------------------------------------------------

      ! 1. First three time steps and more if lfbcwest is .TRUE.
      !    In that case open boundary conditions are FIXED.
      ! --------------------------------------------------------

      IF( ( kt < nit000+3 .AND. .NOT.ln_rstart ) .OR. lfbcwest_trc ) THEN

       DO jn = 1, jptra
         DO ji = fs_niw0, fs_niw1 ! Vector opt.
            DO jk = 1, jpkm1
               DO jj = 1, jpj
                  tra(ji,jj,jk,jn) = tra(ji,jj,jk,jn) * (1. - twmsk(jj,jk)) + &
                                     trcfow(jj,jk,jn)*twmsk(jj,jk)
               END DO
            END DO
         END DO
       END DO

      ELSE

      ! 2. Beyond the fourth time step if lfbcwest_trc is .FALSE.
      ! -----------------------------------------------------
          
         ! Temperature and salinity radiation
         ! ----------------------------------
         !
         !          nib       nibm     nibm2
         !     nib///|   nibm  |  nibm2  |
         !   ///|////|    |    |    |    |   
         !   ---v----f----v----f----v----f-- jj   line
         !   ///|////|    |    |    |    |   
         !   //   ///|         |         |   
         !   // T ///u    T    u    T    u   jj   line
         !   //   ///|         |         |   
         !   ///|////|    |    |    |    |   
         !   ---v----f----v----f----v----f-- jj-1 line
         !   ///|////|    |    |    |    |   
         !         jpiwob    jpiwob+1    jpiwob+2
         !      |         |         |        
         !    jpiwob    jpiwob+1   jpiwob+2
         !
         ! ... radiative conditions + relaxation toward a climatology
         ! ... the phase velocity is taken as the phase velocity of the tangen-
         ! ... tial velocity (here vn), which have been saved in (v_cxwbnd)
         DO ji = fs_niw0, fs_niw1 ! Vector opt.
            DO jk = 1, jpkm1
               DO jj = 2, jpjm1
         ! ... i-phase speed ratio (from averaged of v_cxwbnd)
                  z05cx = (  0.5 * ( v_cxwbnd(jj,jk) + v_cxwbnd(jj-1,jk) ) ) / e1t(ji+1,jj)
                  z05cx = max( z05cx, -1. )
         ! ... z05cx > 0, inflow  zin=0, ztau=1    
         !           < 0, outflow zin=1, ztau=rtauw
                  zin = sign( 1., -1.* z05cx )
                  zin = 0.5*( zin + abs(zin) )
                  ztau = (1.-zin )*rtauwin + zin * rtauw
                  z05cx = z05cx * zin
         ! ... update (ta,sa) with radiative or climatological (t, s)
                  DO jn = 1, jptra
                   tra(ji,jj,jk,jn) = tra(ji,jj,jk,jn) * (1. - twmsk(jj,jk)) +           &
                                      twmsk(jj,jk) * ( ( 1. + z05cx - ztau )         &
                                      * trcwbnd(jj,jk,nib ,nitm,jn) - 2.*z05cx              &
                                      * trcwbnd(jj,jk,nibm,nit,jn ) + ztau * trcfow (jj,jk,jn) ) &
                                      / (1. - z05cx)
                  END DO
               END DO
            END DO
         END DO

      END IF

   END SUBROUTINE obc_trc_west


   SUBROUTINE obc_trc_north ( kt )
      !!------------------------------------------------------------------------------
      !!                 ***  SUBROUTINE obc_trc_north  ***
      !!
      !! ** Purpose :
      !!      Apply the radiation algorithm on north OBC passive tracers tra using the 
      !!      phase velocities calculated in obc_rad_north subroutine in obcrad.F90 module
      !!      If the logical lfbcnorth is .TRUE., there is no radiation but only fixed OBC
      !!
      !!  History :
      !!         ! 95-03 (J.-M. Molines) Original from SPEM
      !!         ! 97-07 (G. Madec, J.-M. Molines) additions
      !!         ! 97-12 (M. Imbard) Mpp adaptation
      !!         ! 00-06 (J.-M. Molines) 
      !!    8.5  ! 02-10 (C. Talandier, A-M. Treguier) F90
      !!         ! 09-11 (C. Dufour) adaptation for passive tracers
      !!------------------------------------------------------------------------------
      !! * Arguments
      INTEGER, INTENT( in ) ::   kt

      !! * Local declaration
      INTEGER ::   ji, jj, jk      ! dummy loop indices
      REAL(wp) ::   z05cx, ztau, zin
      !!------------------------------------------------------------------------------

      ! 1. First three time steps and more if lfbcnorth_trc is .TRUE.
      !    In that case open boundary conditions are FIXED.
      ! --------------------------------------------------------

      IF( ( kt < nit000+3 .AND. .NOT.ln_rstart ) .OR. lfbcnorth_trc ) THEN

        DO jn = 1, jptra
         DO jj = fs_njn0+1, fs_njn1+1  ! Vector opt.
            DO jk = 1, jpkm1
               DO ji = 1, jpi
                  tra(ji,jj,jk,jn)= tra(ji,jj,jk,jn) * (1.-tnmsk(ji,jk)) + &
                                tnmsk(ji,jk) * trcfon(ji,jk,jn)
               END DO
            END DO
         END DO
        END DO 

      ELSE

      ! 2. Beyond the fourth time step if lfbcnorth_trc is .FALSE.
      ! -------------------------------------------------------
          
         ! Passive tracers radiation
         ! ----------------------------------
         !
         !           ji-1   ji   ji   ji +1
         !             |
         !    nib //// u // T // u // T //   jpjnob + 1
         !        /////|//////////////////
         !    nib  ----f----v----f----v---   jpjnob
         !             |         |       
         !      nibm-- u -- T -- u -- T --   jpjnob
         !             |         |            
         !   nibm  ----f----v----f----v---  jpjnob-1
         !             |         |      
         !     nibm2-- u -- T -- T -- T --  jpjnob-1
         !             |         |    
         !   nibm2 ----f----v----f----v---  jpjnob-2
         !             |         |
         !
         ! ... radiative conditions + relaxation toward a climatology
         ! ... the phase velocity is taken as the normal phase velocity of the tangen-
         ! ... tial velocity (here un), which has been saved in (u_cynbnd)
         ! ... jpjnob+1,(jpindp1, jpinfm1)
         DO jj = fs_njn0+1, fs_njn1+1 ! Vector opt.
            DO jk = 1, jpkm1
               DO ji = 2, jpim1
         ! ... j-phase speed ratio (from averaged of vtnbnd)
         !        (bounded by 1)
                  z05cx = ( 0.5 * ( u_cynbnd(ji,jk) + u_cynbnd(ji-1,jk) ) ) / e2t(ji,jj-1)
                  z05cx = min( z05cx, 1. )
         ! ... z05cx=< 0, inflow  zin=0, ztau=1    
         !           > 0, outflow zin=1, ztau=rtaun
                  zin = sign( 1., z05cx )
                  zin = 0.5*( zin + abs(zin) )
         ! ... for inflow rtaunin is used for relaxation coefficient else rtaun
                  ztau = (1.-zin ) * rtaunin + zin * rtaun
                  z05cx = z05cx * zin
         ! ... update (tra) with radiative or climatological (passive tracers)
                  DO jn = 1, jptra
                   tra(ji,jj,jk,jn) = tra(ji,jj,jk,jn) * (1.-tnmsk(ji,jk)) +             &
                                      tnmsk(ji,jk) * ( ( 1. - z05cx - ztau )             &
                                      * trcnbnd(ji,jk,nib ,nitm,jn) + 2.*z05cx                   &
                                      * trcnbnd(ji,jk,nibm,nit,jn ) + ztau * trcfon (ji,jk,jn) ) &
                                      / (1. + z05cx)
                  END DO
               END DO
            END DO
         END DO

      END IF

   END SUBROUTINE obc_trc_north


   SUBROUTINE obc_trc_south ( kt )
      !!------------------------------------------------------------------------------
      !!                ***  SUBROUTINE obc_trc_south  ***
      !!     
      !! ** Purpose :
      !!      Apply the radiation algorithm on south OBC tracers ta, sa using the 
      !!      phase velocities calculated in obc_rad_south subroutine in obcrad.F90 module
      !!      If the logical lfbcsouth is .TRUE., there is no radiation but only fixed OBC
      !!
      !!  History :
      !!         ! 95-03 (J.-M. Molines) Original from SPEM
      !!         ! 97-07 (G. Madec, J.-M. Molines) additions
      !!         ! 97-12 (M. Imbard) Mpp adaptation
      !!         ! 00-06 (J.-M. Molines) 
      !!    8.5  ! 02-10 (C. Talandier, A-M Treguier) F90
      !!------------------------------------------------------------------------------
      !! * Arguments
      INTEGER, INTENT( in ) ::   kt

      !! * Local declaration
      INTEGER ::   ji, jj, jk      ! dummy loop indices
      REAL(wp) ::   z05cx, ztau, zin
      !!------------------------------------------------------------------------------

      ! 1. First three time steps and more if lfbcsouth is .TRUE.
      !    In that case open boundary conditions are FIXED.
      ! --------------------------------------------------------

      IF( ( kt < nit000+3 .AND. .NOT.ln_rstart ) .OR. lfbcsouth_trc ) THEN

        DO jn = 1, jptra
         DO jj = fs_njs0, fs_njs1  ! Vector opt.
            DO jk = 1, jpkm1
               DO ji = 1, jpi
                  tra(ji,jj,jk,jn)= tra(ji,jj,jk,jn) * (1.-tsmsk(ji,jk)) + &
                                    tsmsk(ji,jk) * trcfos(ji,jk,jn)
               END DO
            END DO
         END DO
       END DO

      ELSE

      ! 2. Beyond the fourth time step if lfbcsouth is .FALSE.
      ! -------------------------------------------------------
          
         ! Temperature and salinity radiation
         ! ----------------------------------
         !
         !           ji-1   ji   ji   ji +1
         !             |         |
         !   nibm2 ----f----v----f----v---   jpjsob+2
         !             |         |       
         !   nibm2 --  u -- T -- u -- T --   jpjsob+2
         !             |         |            
         !   nibm  ----f----v----f----v---   jpjsob+1
         !             |         |      
         !    nibm --  u -- T -- T -- T --   jpjsob+1
         !             |         |    
         !   nib  -----f----v----f----v---   jpjsob
         !       //////|/////////|//////// 
         !    nib //// u // T // u // T //   jpjsob 
         !
         !... radiative conditions + relaxation toward a climatology
         !... the phase velocity is taken as the phase velocity of the tangen-
         !... tial velocity (here un), which has been saved in (u_cysbnd)
         !... jpjsob,(jpisdp1, jpisfm1)
         DO jj = fs_njs0, fs_njs1  ! Vector opt.
            DO jk = 1, jpkm1
               DO ji = 2, jpim1
         !... j-phase speed ratio (from averaged of u_cysbnd)
         !       (bounded by 1)
                  z05cx = ( 0.5 * ( u_cysbnd(ji,jk) + u_cysbnd(ji-1,jk) ) ) / e2t(ji,jj+1)
                  z05cx = max( z05cx, -1. )
         !... z05cx > 0, inflow  zin=0, ztau=1
         !          < 0, outflow zin=1, ztau=rtaus
                  zin = sign( 1., -1.* z05cx )
                  zin = 0.5*( zin + abs(zin) )
                  ztau = (1.-zin ) * rtausin + zin * rtaus
                  z05cx = z05cx * zin
         !... update (ta,sa) with radiative or climatological (t, s)
                  DO jn = 1, jptra
                   tra(ji,jj,jk,jn) = tra(ji,jj,jk,jn) * (1.-tsmsk(ji,jk)) +             &
                                      tsmsk(ji,jk) * ( ( 1. + z05cx - ztau )         &
                                      * trcsbnd(ji,jk,nib ,nitm,jn) - 2.*z05cx              &
                                      * trcsbnd(ji,jk,nibm,nit,jn ) + ztau * trcfos (ji,jk,jn) ) &
                                      / (1. - z05cx)
                  END DO
               END DO
            END DO
         END DO

      END IF

   END SUBROUTINE obc_trc_south

#else
   !!---------------------------------------------------------------------------------
   !!   Default option                                                    Empty module
   !!---------------------------------------------------------------------------------
CONTAINS
   SUBROUTINE obc_trc      ! Empty routine
   END SUBROUTINE obc_trc
#endif

   !!=================================================================================
END MODULE obctrc
