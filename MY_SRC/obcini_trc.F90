 MODULE obcini_trc
   !!==========================================================================================
   !!                       ***  MODULE  obcini_trc  ***
   !! OBC initial state :  Open boundary initial state for passive tracers
   !!==========================================================================================
   !!  History :       !  2009-11  (C. Dufour)  adapted from obcini.F90
   !!                  !  2010-03  (C. Dufour)  adapted to NEMOv3.2
   !!                  !  2011-05  (A. Albert)  adapted to NEMOv3.4
   !!------------------------------------------------------------------------------------------
#if defined key_obc_trc
   !!------------------------------------------------------------------------------------------
   !!   'key_obc_trc'                                  Open Boundary Conditions
   !!------------------------------------------------------------------------------------------
   !!   obc_init_trc       : initialization for the open boundary condition for passive tracers
   !!------------------------------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers variables
   USE dom_oce         ! ocean space and time domain variables
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE phycst          ! physical constants
   USE obc_oce         ! open boundary condition: ocean
   USE obcdta          ! open boundary condition: data
   USE obcdta_trc          ! open boundary condition: data
   USE in_out_manager  ! I/O units
   USE lib_mpp         ! MPP library
   USE dynspg_oce      ! flag lk_dynspg_flt
   USE obcini          ! initialization for dynamics tracers
   USE trc             ! TOP

   IMPLICIT NONE
   PRIVATE

   PUBLIC   obc_init_trc   ! routine called by opa.F90

   !! * Substitutions
#  include "obc_vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: obcini.F90 2715 2011-03-30 15:58:35Z rblod $ 
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS
   
   SUBROUTINE obc_init_trc
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE obc_init_trc  ***
      !!         
      !! ** Purpose :   Initialization of the dynamics and tracer fields at 
      !!              the open boundaries.
      !!
      !! ** Method  :   initialization of open boundary variables
      !!      (u, v) over 3 time step and 3 rows
      !!      (t, s) over 2 time step and 2 rows
      !!      if ln_rstart = .FALSE. : no restart, fields set to zero
      !!      if ln_rstart = .TRUE.  : restart, fields are read in a file 
      !!      if rdpxxx = 0 then lfbc is set true for this boundary.
      !!
      !! ** Input   :   restart.obc file, restart file for open boundaries 
      !!----------------------------------------------------------------------
      USE obcrst_trc,   ONLY :   obc_rst_read_trc   ! Make obc_rst_read routine available
      !!
      INTEGER  ::   ji, jj, istop , inumfbc, jn
      INTEGER, DIMENSION(4) ::   icorner
      REAL(wp), DIMENSION(2) ::   ztestmask

      ! Definition of a tracer as a structure
      TYPE PTRACER
         CHARACTER(len = 20)  :: namt   !: name of tracer
         LOGICAL              :: clim   !: climatological file or not
         INTEGER              :: nobc   !: read in a file or not
         LOGICAL              :: bini   !: obc from the very initial state
      END TYPE PTRACER

      TYPE(PTRACER) , DIMENSION(jptra) :: obc_tracer
      !!
      NAMELIST/namobc_trc/ rn_dpein_trc, rn_dpwin_trc, rn_dpnin_trc, rn_dpsin_trc,       &
         &                 rn_dpeob_trc, rn_dpwob_trc, rn_dpnob_trc, rn_dpsob_trc,       &
         &                 cn_obcdta_trc,ln_vol_cst_trc, ln_obc_fla_trc, obc_tracer
!         &                 nn_obcdta_trc, cn_obcdta_trc,                                 &
!         &                 ln_obc_clim_trc, ln_vol_cst_trc, ln_obc_fla_trc, ln_obc_ini
      !!----------------------------------------------------------------------

      ! CD add CALL ctl_opn to debug
      CALL ctl_opn( numnat, 'namelist_top', 'OLD', 'FORMATTED', 'SEQUENTIAL', 1, numout, .FALSE. )
!      REWIND( numnat )              ! Namelist namobc : open boundaries
!      READ  ( numnat, namobc_trc )
  
      ! Namelist default values
      DO jn = 1, jptra
        nn_obcdta_trc(jn) = 1
        ln_obc_clim_trc(jn) = .TRUE.
        ln_obc_ini(jn) = .FALSE.
      ENDDO

      REWIND( numnat )              ! Namelist namobc : open boundaries
      READ  ( numnat, namobc_trc )

      ! Upload the parameter values from the namelist
      DO jn = 1, jptra
         nam_obc_trc(jn)     = obc_tracer(jn)%namt
         ln_obc_clim_trc(jn) = obc_tracer(jn)%clim
         nn_obcdta_trc(jn)   = obc_tracer(jn)%nobc
         ln_obc_ini(jn)      = obc_tracer(jn)%bini
      END DO

      IF(lwp) THEN                   ! control print
          WRITE(numout,*)
          WRITE(numout,*) ' Namelist obc: namobc_trc'
          DO jn = 1, jptra
            WRITE(numout,*) '   short name of tracer            : ', TRIM(nam_obc_trc(jn))
            WRITE(numout,*) '   climatological file             : ', ln_obc_clim_trc(jn)
            WRITE(numout,*) '   read a file or not              : ', nn_obcdta_trc(jn)
            WRITE(numout,*) '   obc from very initial state     : ', ln_obc_ini(jn)
          ENDDO
      ENDIF

      CLOSE ( numnat)

      ! convert DOCTOR namelist name into the OLD names
      nobc_dta_trc = nn_obcdta_trc
      cffile_trc   = cn_obcdta_trc
      rdpein_trc   = rn_dpein_trc
      rdpwin_trc   = rn_dpwin_trc
      rdpsin_trc   = rn_dpsin_trc
      rdpnin_trc   = rn_dpnin_trc
      rdpeob_trc   = rn_dpeob_trc
      rdpwob_trc   = rn_dpwob_trc
      rdpsob_trc   = rn_dpsob_trc
      rdpnob_trc   = rn_dpnob_trc
      volemp   = rn_volemp
      !                              ! allocate obc arrays
      IF( obc_oce_trc_alloc() /= 0 )   CALL ctl_stop( 'STOP', 'obc_init_trc : unable to allocate obc_oce arrays' )
      IF( obc_dta_trc_alloc() /= 0 )   CALL ctl_stop( 'STOP', 'obc_init_trc : unable to allocate obc_dta arrays' )


      ! By security we set rdpxin and rdpxob respectively to 1. and 15. if the corresponding OBC is not activated
      IF( .NOT.lp_obc_east_trc  ) THEN   ;   rdpein_trc = 1.   ;   rdpeob_trc = 15.   ;   END IF
      IF( .NOT.lp_obc_west_trc  ) THEN   ;   rdpwin_trc = 1.   ;   rdpwob_trc = 15.   ;   END IF
      IF( .NOT.lp_obc_north_trc ) THEN   ;   rdpnin_trc = 1.   ;   rdpnob_trc = 15.   ;   END IF
      IF( .NOT.lp_obc_south_trc ) THEN   ;   rdpsin_trc = 1.   ;   rdpsob_trc = 15.   ;   END IF

      ! number of open boudaries and open boundary indicators
      nbobc = 0
      IF( lp_obc_east_trc  )   nbobc = nbobc + 1
      IF( lp_obc_west_trc  )   nbobc = nbobc + 1
      IF( lp_obc_north_trc )   nbobc = nbobc + 1
      IF( lp_obc_south_trc )   nbobc = nbobc + 1

      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'obc_init_trc : initialization of open boundaries for passive tracers'
      IF(lwp) WRITE(numout,*) '~~~~~~~~'
      IF(lwp) WRITE(numout,*) '   Number of open boundaries    nbobc = ', nbobc
      IF(lwp) WRITE(numout,*)

      ! control prints
      IF(lwp) WRITE(numout,*) '   Namelist namobc_trc'
      IF(lwp) WRITE(numout,*) '      name of obc passive tracers                    name obc tracer = ', nam_obc_trc
      IF(lwp) WRITE(numout,*) '      climatology (true) or not                      ln_obc_clim_trc = ', ln_obc_clim_trc
      IF(lwp) WRITE(numout,*) '      data in file (=1) or initial state used (=0)   nn_obcdta_trc   = ', nn_obcdta_trc
      IF(lwp) WRITE(numout,*) '      climatology (true) or not                      ln_obc_ini      = ', ln_obc_ini
      IF(lwp) WRITE(numout,*) '      vol_cst (true) or not:                         ln_vol_cst_trc  = ', ln_vol_cst_trc
      IF(lwp) WRITE(numout,*) ' '
      IF(lwp) WRITE(numout,*) '   WARNING                                                  '
      IF(lwp) WRITE(numout,*) '      Flather"s algorithm is applied with explicit free surface scheme                 '
      IF(lwp) WRITE(numout,*) '      or with free surface time-splitting scheme                                       '
      IF(lwp) WRITE(numout,*) '      Nor radiation neither relaxation is allowed with explicit free surface scheme:   '
      IF(lwp) WRITE(numout,*) '      Radiation and/or relaxation is allowed with free surface time-splitting scheme '
      IF(lwp) WRITE(numout,*) '      depending of the choice of rdpXin = rdpXob  = 0. for open boundaries             '
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) '      For the filtered free surface case,                                              '
      IF(lwp) WRITE(numout,*) '      radiation, relaxation or presciption of data can be applied                      '
      IF(lwp) WRITE(numout,*)

      IF( lwp.AND.lp_obc_east_trc ) THEN
         WRITE(numout,*) '      East open boundary for passive tracers :'
         WRITE(numout,*) '         i index                    jpieob   = ', jpieob
         WRITE(numout,*) '         damping time scale (days)  rn_dpeob_trc = ', rn_dpeob_trc
         WRITE(numout,*) '         damping time scale (days)  rn_dpein_trc = ', rn_dpein_trc
      ENDIF

      IF( lwp.AND.lp_obc_west_trc ) THEN
         WRITE(numout,*) '      West open boundary for passive tracers :'
         WRITE(numout,*) '         i index                    jpiwob   = ', jpiwob
         WRITE(numout,*) '         damping time scale (days)  rn_dpwob_trc = ', rn_dpwob_trc
         WRITE(numout,*) '         damping time scale (days)  rn_dpwin_trc = ', rn_dpwin_trc
      ENDIF

      IF( lwp.AND.lp_obc_north_trc ) THEN
         WRITE(numout,*) '      North open boundary for passive tracers :'
         WRITE(numout,*) '         j index                    jpjnob   = ', jpjnob
         WRITE(numout,*) '         damping time scale (days)  rn_dpnob_trc = ', rn_dpnob_trc
         WRITE(numout,*) '         damping time scale (days)  rn_dpnin_trc = ', rn_dpnin_trc
      ENDIF

      IF( lwp.AND.lp_obc_south_trc ) THEN
         WRITE(numout,*) '      South open boundary for passive tracers:'
         WRITE(numout,*) '         j index                    jpjsob   = ', jpjsob
         WRITE(numout,*) '         damping time scale (days)  rn_dpsob_trc = ', rn_dpsob_trc
         WRITE(numout,*) '         damping time scale (days)  rn_dpsin_trc = ', rn_dpsin_trc
         WRITE(numout,*)
      ENDIF

      IF( nbobc >= 2 .AND. jperio /= 0 )   &
         &   CALL ctl_stop( ' Cyclic or symmetric, and open boundary condition are not compatible' )

      ! 1. Initialisation of constants 
      ! ------------------------------
      ! ...                          convert rdp$ob in seconds
      ! Fixed Bdy flag              inbound                outbound
      lfbceast_trc  = .FALSE.   ;   rdpein_trc = rdpein_trc * rday    ;   rdpeob_trc = rdpeob_trc * rday
      lfbcwest_trc  = .FALSE.   ;   rdpwin_trc = rdpwin_trc * rday    ;   rdpwob_trc = rdpwob_trc * rday
      lfbcnorth_trc = .FALSE.   ;   rdpnin_trc = rdpnin_trc * rday    ;   rdpnob_trc = rdpnob_trc * rday
      lfbcsouth_trc = .FALSE.   ;   rdpsin_trc = rdpsin_trc * rday    ;   rdpsob_trc = rdpsob_trc * rday
      inumfbc = 0
      ! ... look for Fixed Boundaries (rdp = 0 )
      ! ... When specified, lbcxxx flags are set to TRUE and rdpxxx are set to
      ! ...  a small arbitrary value, (to avoid division by zero further on). 
      ! ...  rdpxxx is not used anymore.
      IF( lp_obc_east_trc )  THEN
         IF( (rdpein_trc+rdpeob_trc) == 0 )  THEN
            lfbceast_trc = .TRUE.   ;   rdpein_trc = 1e-3   ;   rdpeob_trc = 1e-3
            inumfbc = inumfbc+1
         ELSEIF ( (rdpein_trc*rdpeob_trc) == 0 )  THEN
            CALL ctl_stop( 'obc_init_trc : rn_dpein_trc & rn_dpeob_trc must be both zero or non zero' )
         END IF
      END IF

      IF( lp_obc_west_trc )  THEN
         IF( (rdpwin_trc + rdpwob_trc) == 0 )  THEN
            lfbcwest_trc = .TRUE.     ;     rdpwin_trc = 1e-3     ;     rdpwob_trc = 1e-3
            inumfbc = inumfbc+1
         ELSEIF ( (rdpwin_trc*rdpwob_trc) == 0 )  THEN
            CALL ctl_stop( 'obc_init_trc : rn_dpwin_trc & rn_dpwob_trc must be both zero or non zero' )
         END IF
      END IF
      IF( lp_obc_north_trc )  THEN
         IF( (rdpnin_trc + rdpnob_trc) == 0 )  THEN
            lfbcnorth_trc = .TRUE.     ;     rdpnin_trc = 1e-3     ;     rdpnob_trc = 1e-3
            inumfbc = inumfbc+1
         ELSEIF ( (rdpnin_trc*rdpnob_trc) == 0 )  THEN
            CALL ctl_stop( 'obc_init_trc : rn_dpnin_trc & rn_dpnob_trc must be both zero or non zero' )
         END IF
      END IF
      IF( lp_obc_south_trc )  THEN
         IF( (rdpsin_trc + rdpsob_trc) == 0 )  THEN
            lfbcsouth_trc = .TRUE.   ;   rdpsin_trc = 1e-3   ;   rdpsob_trc = 1e-3
            inumfbc = inumfbc+1
         ELSEIF ( (rdpsin_trc*rdpsob_trc) == 0 )  THEN
            CALL ctl_stop( 'obc_init_trc : rn_dpsin_trc & rn_dpsob_trc must be both zero or non zero' )
         END IF
      END IF

      ! 2.  Clever mpp indices for loops on the open boundaries. 
      !     The loops will be performed only on the processors 
      !     that contain a given open boundary.
      ! --------------------------------------------------------

      IF( lp_obc_east_trc ) THEN
         ! ...   mpp initialization
         nie0   = max( 1, min(jpieob   - nimpp+1, jpi     ) )
         nie1   = max( 0, min(jpieob   - nimpp+1, jpi - 1 ) )
         nie0p1 = max( 1, min(jpieob+1 - nimpp+1, jpi     ) )
         nie1p1 = max( 0, min(jpieob+1 - nimpp+1, jpi - 1 ) )
         nie0m1 = max( 1, min(jpieob-1 - nimpp+1, jpi     ) )
         nie1m1 = max( 0, min(jpieob-1 - nimpp+1, jpi - 1 ) )
         nje0   = max( 2, min(jpjed    - njmpp+1, jpj     ) )
         nje1   = max( 0, min(jpjef    - njmpp+1, jpj - 1 ) )
         nje0p1 = max( 1, min(jpjedp1  - njmpp+1, jpj     ) )
         nje0m1 = max( 1, min(jpjed    - njmpp+1, jpj     ) )
         nje1m1 = max( 0, min(jpjefm1  - njmpp+1, jpj - 1 ) )
         nje1m2 = max( 0, min(jpjefm1-1- njmpp+1, jpj - 1 ) )
         IF(lwp) THEN
            IF( lfbceast_trc ) THEN
               WRITE(numout,*)'     '
               WRITE(numout,*)'         Specified East Open Boundary'
            ELSE
               WRITE(numout,*)'     '
               WRITE(numout,*)'         Radiative East Open Boundary'
            END IF
         END IF
      END IF

      IF( lp_obc_west_trc ) THEN
         ! ...   mpp initialization
         niw0   = max( 1, min(jpiwob   - nimpp+1, jpi     ) )
         niw1   = max( 0, min(jpiwob   - nimpp+1, jpi - 1 ) )
         niw0p1 = max( 1, min(jpiwob+1 - nimpp+1, jpi     ) )
         niw1p1 = max( 0, min(jpiwob+1 - nimpp+1, jpi - 1 ) )
         njw0   = max( 2, min(jpjwd    - njmpp+1, jpj     ) )
         njw1   = max( 0, min(jpjwf    - njmpp+1, jpj - 1 ) )
         njw0p1 = max( 1, min(jpjwdp1  - njmpp+1, jpj     ) )
         njw0m1 = max( 1, min(jpjwd    - njmpp+1, jpj     ) )
         njw1m1 = max( 0, min(jpjwfm1  - njmpp+1, jpj - 1 ) )
         njw1m2 = max( 0, min(jpjwfm1-1- njmpp+1, jpj - 1 ) )
         IF(lwp) THEN
            IF( lfbcwest_trc ) THEN
               WRITE(numout,*)'     '
               WRITE(numout,*)'         Specified West Open Boundary'
            ELSE
               WRITE(numout,*)'     '
               WRITE(numout,*)'         Radiative West Open Boundary'
            END IF
         END IF
      END IF
 
      IF( lp_obc_north_trc ) THEN
         ! ...   mpp initialization
         nin0   = max( 2, min(jpind    - nimpp+1, jpi     ) )
         nin1   = max( 0, min(jpinf    - nimpp+1, jpi - 1 ) )
         nin0p1 = max( 1, min(jpindp1  - nimpp+1, jpi     ) )
         nin0m1 = max( 1, min(jpind    - nimpp+1, jpi     ) )
         nin1m1 = max( 0, min(jpinfm1  - nimpp+1, jpi - 1 ) )
         nin1m2 = max( 0, min(jpinfm1-1- nimpp+1, jpi - 1 ) )
         njn0   = max( 1, min(jpjnob   - njmpp+1, jpj     ) )
         njn1   = max( 0, min(jpjnob   - njmpp+1, jpj - 1 ) )
         njn0p1 = max( 1, min(jpjnob+1 - njmpp+1, jpj     ) )
         njn1p1 = max( 0, min(jpjnob+1 - njmpp+1, jpj - 1 ) )
         njn0m1 = max( 1, min(jpjnob-1 - njmpp+1, jpj     ) )
         njn1m1 = max( 0, min(jpjnob-1 - njmpp+1, jpj - 1 ) )
         IF(lwp) THEN
            IF( lfbcnorth_trc ) THEN
               WRITE(numout,*)'     '
               WRITE(numout,*)'         Specified North Open Boundary'
            ELSE
               WRITE(numout,*)'     '
               WRITE(numout,*)'         Radiative North Open Boundary'
            END IF
         END IF
      END IF

      IF( lp_obc_south_trc ) THEN
         ! ...   mpp initialization
         nis0   = max( 2, min(jpisd    - nimpp+1, jpi     ) )
         nis1   = max( 0, min(jpisf    - nimpp+1, jpi - 1 ) )
         nis0p1 = max( 1, min(jpisdp1  - nimpp+1, jpi     ) )
         nis0m1 = max( 1, min(jpisd    - nimpp+1, jpi     ) )
         nis1m1 = max( 0, min(jpisfm1  - nimpp+1, jpi - 1 ) )
         nis1m2 = max( 0, min(jpisfm1-1- nimpp+1, jpi - 1 ) )
         njs0   = max( 1, min(jpjsob   - njmpp+1, jpj     ) )
         njs1   = max( 0, min(jpjsob   - njmpp+1, jpj - 1 ) )
         njs0p1 = max( 1, min(jpjsob+1 - njmpp+1, jpj     ) )
         njs1p1 = max( 0, min(jpjsob+1 - njmpp+1, jpj - 1 ) )
         IF(lwp) THEN
            IF( lfbcsouth_trc ) THEN
               WRITE(numout,*)'     '
               WRITE(numout,*)'         Specified South Open Boundary'
            ELSE
               WRITE(numout,*)'     '
               WRITE(numout,*)'         Radiative South Open Boundary'
            END IF
         END IF
      END IF

      ! 3. mask correction for OBCs
      ! ---------------------------
      !
      ! Already in OBC dynamics!!!!!
      !
      ! 6. Initialization of open boundary variables (u, v, t, s)
      ! --------------------------------------------------------------
      !   only if at least one boundary is  radiative 
      IF ( inumfbc < nbobc .AND.  ln_rstart ) THEN
         !  Restart from restart.obc
         CALL obc_rst_read_trc
      ELSE

          ! ... Initialization to zero of radiation arrays
          !     for passive tracers.
          
          trcebnd(:,:,:,:,:)   = 0.e0   ;   trcnbnd(:,:,:,:,:)   = 0.e0
          trcwbnd(:,:,:,:,:)   = 0.e0   ;   trcsbnd(:,:,:,:,:)   = 0.e0 

      END IF

      ! 7. Control print
      ! -----------------------------------------------------------------
      ! Already done in OBC dynamics!!!!   

   END SUBROUTINE obc_init_trc

#else
   !!---------------------------------------------------------------------------------
   !!   Dummy module                                                NO open boundaries
   !!---------------------------------------------------------------------------------
CONTAINS
   SUBROUTINE obc_init_trc     ! Dummy routine
   END SUBROUTINE obc_init_trc
#endif

   !!=================================================================================
END MODULE obcini_trc
