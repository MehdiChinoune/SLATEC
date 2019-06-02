!** DPOPT
SUBROUTINE DPOPT(Prgopt,Mrelas,Nvars,Info,Csc,Ibasis,Ropt,Intopt,Lopt)
  !>
  !  Subsidiary to DSPLP
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      DOUBLE PRECISION (SPOPT-S, DPOPT-D)
  !***
  ! **Author:**  (UNKNOWN)
  !***
  ! **Description:**
  !
  !     THE EDITING REQUIRED TO CONVERT THIS SUBROUTINE FROM SINGLE TO
  !     DOUBLE PRECISION INVOLVES THE FOLLOWING CHARACTER STRING CHANGES.
  !
  !     USE AN EDITING COMMAND (CHANGE) /STRING-1/(TO)STRING-2/.
  !     /REAL (12 BLANKS)/DOUBLE PRECISION/,/R1MACH/D1MACH/,/E0/D0/
  !
  !     REVISED 821122-1045
  !     REVISED YYMMDD-HHMM
  !
  !     THIS SUBROUTINE PROCESSES THE OPTION VECTOR, PRGOPT(*),
  !     AND VALIDATES ANY MODIFIED DATA.
  !
  !***
  ! **See also:**  DSPLP
  !***
  ! **Routines called:**  D1MACH, XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   811215  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890605  Removed unreferenced labels.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900328  Added TYPE section.  (WRB)
  !   900510  Fixed an error message.  (RWC)
  USE service, ONLY : XERMSG, D1MACH
  INTEGER :: Info, Mrelas, Nvars
  INTEGER :: Ibasis(Nvars+Mrelas), Intopt(08)
  REAL(8) :: Csc(Nvars), Prgopt(:), Ropt(07)
  LOGICAL :: Lopt(8)
  INTEGER :: i, iadbig, ictmax, ictopt, idg, iopt, ipagef, isave, itbrc, itest, &
    j, key, kprint, last, lds, lprg, mxitlp, n20043, n20053, n20096, nerr, next, npp
  REAL(8) :: abig, asmall, costsc, eps, one, tolls, tune, zero, tolabs
  LOGICAL :: contin, usrbas, sizeup, savedt, colscp, cstscp, minprb, stpedg
  !
  !* FIRST EXECUTABLE STATEMENT  DPOPT
  iopt = 1
  zero = 0.D0
  one = 1.D0
  asmall = 0.D0
  abig = 0.D0
  costsc = 0.D0
  !
  !
  !     PROCEDURE (INITIALIZE PARAMETERS AND PROCESS USER OPTIONS)
  contin = .FALSE.
  usrbas = .FALSE.
  sizeup = .FALSE.
  savedt = .FALSE.
  colscp = .FALSE.
  cstscp = .FALSE.
  minprb = .TRUE.
  stpedg = .TRUE.
  !
  !     GET THE MACHINE REL. FLOATING POINT ACCURACY VALUE FROM THE
  !     LIBRARY SUBPROGRAM, D1MACH( ).
  eps = D1MACH(4)
  tolls = D1MACH(4)
  tune = one
  tolabs = zero
  !
  !     DEFINE NOMINAL FILE NUMBERS FOR MATRIX PAGES AND DATA SAVING.
  ipagef = 1
  isave = 2
  itbrc = 10
  mxitlp = 3*(Nvars+Mrelas)
  kprint = 0
  idg = -4
  npp = Nvars
  lprg = 0
  !
  last = 1
  iadbig = 10000
  ictmax = 1000
  ictopt = 0
  DO
    next = INT( Prgopt(last) )
    IF ( next<=0.OR.next>iadbig ) THEN
      !
      !     THE CHECKS FOR SMALL OR LARGE VALUES OF NEXT ARE TO PREVENT
      !     WORKING WITH UNDEFINED DATA.
      nerr = 14
      CALL XERMSG('DPOPT',&
        'IN DSPLP, THE USER OPTION ARRAY HAS UNDEFINED DATA.',nerr,iopt)
      Info = -nerr
      RETURN
    ELSEIF ( next==1 ) THEN
      !
      !     PROCEDURE (VALIDATE OPTIONALLY MODIFIED DATA)
      !
      !     IF USER HAS DEFINED THE BASIS, CHECK FOR VALIDITY OF INDICES.
      IF ( usrbas ) THEN
        i = 1
        n20096 = Mrelas
        DO WHILE ( (n20096-i)>=0 )
          itest = Ibasis(i)
          IF ( itest>0.AND.itest<=(Nvars+Mrelas) ) THEN
            i = i + 1
          ELSE
            nerr = 16
            CALL XERMSG('DPOPT',&
              'IN DSPLP, AN INDEX OF USER-SUPPLIED BASIS IS OUT OF RANGE.'&
              ,nerr,iopt)
            Info = -nerr
            RETURN
          END IF
        END DO
      END IF
      !
      !     IF USER HAS PROVIDED SIZE PARAMETERS, MAKE SURE THEY ARE ORDERED
      !     AND POSITIVE.
      IF ( sizeup ) THEN
        IF ( asmall<=zero.OR.abig<asmall ) THEN
          nerr = 17
          CALL XERMSG('DPOPT',&
            'IN DSPLP, SIZE PARAMETERS FOR MATRIX MUST BE SMALLEST AND LARGEST MAGNITUDES OF NONZERO ENTRIES.',nerr,iopt)
          Info = -nerr
          RETURN
        END IF
      END IF
      !
      !     THE NUMBER OF ITERATIONS OF REV. SIMPLEX STEPS MUST BE POSITIVE.
      IF ( mxitlp>0 ) THEN
        !
        !     CHECK THAT SAVE AND PAGE FILE NUMBERS ARE DEFINED AND NOT EQUAL.
        IF ( isave<=0.OR.ipagef<=0.OR.(isave==ipagef) ) EXIT
        !
        Lopt(1) = contin
        Lopt(2) = usrbas
        Lopt(3) = sizeup
        Lopt(4) = savedt
        Lopt(5) = colscp
        Lopt(6) = cstscp
        Lopt(7) = minprb
        Lopt(8) = stpedg
        !
        Intopt(1) = idg
        Intopt(2) = ipagef
        Intopt(3) = isave
        Intopt(4) = mxitlp
        Intopt(5) = kprint
        Intopt(6) = itbrc
        Intopt(7) = npp
        Intopt(8) = lprg
        !
        Ropt(1) = eps
        Ropt(2) = asmall
        Ropt(3) = abig
        Ropt(4) = costsc
        Ropt(5) = tolls
        Ropt(6) = tune
        Ropt(7) = tolabs
        RETURN
      ELSE
        nerr = 18
        CALL XERMSG('DPOPT',&
          'IN DSPLP, THE NUMBER OF REVISED SIMPLEX STEPS BETWEEN CHECK-POINTS MUST BE POSITIVE.',nerr,iopt)
        Info = -nerr
        RETURN
      END IF
    ELSEIF ( ictopt<=ictmax ) THEN
      key = INT( Prgopt(last+1) )
      !
      !     IF KEY = 50, THIS IS TO BE A MAXIMIZATION PROBLEM
      !     INSTEAD OF A MINIMIZATION PROBLEM.
      IF ( key==50 ) THEN
        minprb = Prgopt(last+2)==zero
        lds = 3
        !
        !     IF KEY = 51, THE LEVEL OF OUTPUT IS BEING MODIFIED.
        !     KPRINT = 0, NO OUTPUT
        !            = 1, SUMMARY OUTPUT
        !            = 2, LOTS OF OUTPUT
        !            = 3, EVEN MORE OUTPUT
      ELSEIF ( key==51 ) THEN
        kprint = INT( Prgopt(last+2) )
        lds = 3
        !
        !     IF KEY = 52, REDEFINE THE FORMAT AND PRECISION USED
        !     IN THE OUTPUT.
      ELSEIF ( key==52 ) THEN
        IF ( Prgopt(last+2)/=zero ) idg = INT( Prgopt(last+3) )
        lds = 4
        !
        !     IF KEY = 53, THE ALLOTTED SPACE FOR THE SPARSE MATRIX
        !     STORAGE AND/OR SPARSE EQUATION SOLVING HAS BEEN CHANGED.
        !    (PROCESSED IN DSPLP(). THIS IS TO COMPUTE THE LENGTH OF PRGOPT(*).)
      ELSEIF ( key==53 ) THEN
        lds = 5
        !
        !     IF KEY = 54, REDEFINE THE FILE NUMBER WHERE THE PAGES
        !     FOR THE SPARSE MATRIX ARE STORED.
      ELSEIF ( key==54 ) THEN
        IF ( Prgopt(last+2)/=zero ) ipagef = INT( Prgopt(last+3) )
        lds = 4
        !
        !     IF KEY = 55,  A CONTINUATION FOR A PROBLEM MAY BE REQUESTED.
      ELSEIF ( key==55 ) THEN
        contin = Prgopt(last+2)/=zero
        lds = 3
        !
        !     IF KEY = 56, REDEFINE THE FILE NUMBER WHERE THE SAVED DATA
        !     WILL BE STORED.
      ELSEIF ( key==56 ) THEN
        IF ( Prgopt(last+2)/=zero ) isave = INT( Prgopt(last+3) )
        lds = 4
        !
        !     IF KEY = 57, SAVE DATA (ON EXTERNAL FILE)  AT MXITLP ITERATIONS OR
        !     THE OPTIMUM, WHICHEVER COMES FIRST.
      ELSEIF ( key==57 ) THEN
        savedt = Prgopt(last+2)/=zero
        lds = 3
        !
        !     IF KEY = 58,  SEE IF PROBLEM IS TO RUN ONLY A GIVEN
        !     NUMBER OF ITERATIONS.
      ELSEIF ( key==58 ) THEN
        IF ( Prgopt(last+2)/=zero ) mxitlp = INT( Prgopt(last+3) )
        lds = 4
        !
        !     IF KEY = 59,  SEE IF USER PROVIDES THE BASIS INDICES.
      ELSEIF ( key==59 ) THEN
        usrbas = Prgopt(last+2)/=zero
        IF ( .NOT.(usrbas) ) THEN
          lds = Mrelas + 3
        ELSE
          i = 1
          n20043 = Mrelas
          DO WHILE ( (n20043-i)>=0 )
            Ibasis(i) = INT( Prgopt(last+2+i) )
            i = i + 1
          END DO
          lds = Mrelas + 3
        END IF
        !
        !     IF KEY = 60,  SEE IF USER HAS PROVIDED SCALING OF COLUMNS.
      ELSEIF ( key==60 ) THEN
        colscp = Prgopt(last+2)/=zero
        IF ( .NOT.(colscp) ) THEN
          lds = Nvars + 3
        ELSE
          j = 1
          n20053 = Nvars
          DO WHILE ( (n20053-j)>=0 )
            Csc(j) = ABS(Prgopt(last+2+j))
            j = j + 1
          END DO
          lds = Nvars + 3
        END IF
        !
        !     IF KEY = 61,  SEE IF USER HAS PROVIDED SCALING OF COSTS.
      ELSEIF ( key==61 ) THEN
        cstscp = Prgopt(last+2)/=zero
        IF ( cstscp ) costsc = Prgopt(last+3)
        lds = 4
        !
        !     IF KEY = 62,  SEE IF SIZE PARAMETERS ARE PROVIDED WITH THE DATA.
        !     THESE WILL BE CHECKED AGAINST THE MATRIX ELEMENT SIZES LATER.
      ELSEIF ( key==62 ) THEN
        sizeup = Prgopt(last+2)/=zero
        IF ( sizeup ) THEN
          asmall = Prgopt(last+3)
          abig = Prgopt(last+4)
        END IF
        lds = 5
        !
        !     IF KEY = 63, SEE IF TOLERANCE FOR LINEAR SYSTEM RESIDUAL ERROR IS
        !     PROVIDED.
      ELSEIF ( key==63 ) THEN
        IF ( Prgopt(last+2)/=zero ) tolls = MAX(eps,Prgopt(last+3))
        lds = 4
        !
        !     IF KEY = 64,  SEE IF MINIMUM REDUCED COST OR STEEPEST EDGE
        !     DESCENT IS TO BE USED FOR SELECTING VARIABLES TO ENTER BASIS.
      ELSEIF ( key==64 ) THEN
        stpedg = Prgopt(last+2)==zero
        lds = 3
        !
        !     IF KEY = 65, SET THE NUMBER OF ITERATIONS BETWEEN RECALCULATING
        !     THE ERROR IN THE PRIMAL SOLUTION.
      ELSEIF ( key==65 ) THEN
        IF ( Prgopt(last+2)/=zero ) itbrc = INT( MAX(one,Prgopt(last+3)) )
        lds = 4
        !
        !     IF KEY = 66, SET THE NUMBER OF NEGATIVE REDUCED COSTS TO BE FOUND
        !     IN THE PARTIAL PRICING STRATEGY.
      ELSEIF ( key==66 ) THEN
        IF ( Prgopt(last+2)/=zero ) THEN
          npp = INT( MAX(Prgopt(last+3),one) )
          npp = MIN(npp,Nvars)
        END IF
        lds = 4
        !     IF KEY = 67, CHANGE THE TUNING PARAMETER TO APPLY TO THE ERROR
        !     ESTIMATES FOR THE PRIMAL AND DUAL SYSTEMS.
      ELSEIF ( key==67 ) THEN
        IF ( Prgopt(last+2)/=zero ) tune = ABS(Prgopt(last+3))
        lds = 4
      ELSEIF ( key==68 ) THEN
        lds = 6
        !
        !     RESET THE ABSOLUTE TOLERANCE TO BE USED ON THE FEASIBILITY
        !     DECISION PROVIDED THE RELATIVE ERROR TEST FAILED.
      ELSEIF ( key==69 ) THEN
        IF ( Prgopt(last+2)/=zero ) tolabs = Prgopt(last+3)
        lds = 4
      END IF
      !
      ictopt = ictopt + 1
      last = next
      lprg = lprg + lds
    ELSE
      nerr = 15
      CALL XERMSG('DPOPT',&
        'IN DSPLP, OPTION ARRAY PROCESSING IS CYCLING.',nerr,iopt)
      Info = -nerr
      RETURN
    END IF
  END DO
  nerr = 19
  CALL XERMSG('DPOPT',&
    'IN DSPLP, FILE NUMBERS FOR SAVED DATA AND MATRIX PAGES MUST BE POSITIVE AND NOT EQUAL.',nerr,iopt)
  Info = -nerr
  RETURN
END SUBROUTINE DPOPT
