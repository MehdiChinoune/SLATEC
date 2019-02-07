!*==EVCHCK.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK EVCHCK
SUBROUTINE EVCHCK(Lout,Kprint,Npts,Xev,Fev,Dev,Fev2,Fail)
  IMPLICIT NONE
  !*--EVCHCK5
  !***BEGIN PROLOGUE  EVCHCK
  !***SUBSIDIARY
  !***PURPOSE  Test evaluation accuracy of CHFDV and CHFEV for PCHQK1.
  !***LIBRARY   SLATEC (PCHIP)
  !***TYPE      SINGLE PRECISION (EVCHCK-S, DEVCHK-D)
  !***KEYWORDS  PCHIP EVALUATOR QUICK CHECK
  !***AUTHOR  Fritsch, F. N., (LLNL)
  !***DESCRIPTION
  !
  ! -------- CODE TO TEST EVALUATION ACCURACY OF CHFDV AND CHFEV --------
  !
  !     USING FUNCTION AND DERIVATIVE VALUES FROM A CUBIC (COMPUTED IN
  !     DOUBLE PRECISION) AT NINT DIFFERENT (X1,X2) PAIRS:
  !     1. CHECKS THAT CHFDV AND CHFEV BOTH REPRODUCE ENDPOINT VALUES.
  !     2. EVALUATES AT NPTS POINTS, 10 OF WHICH ARE OUTSIDE THE INTERVAL
  !        AND:
  !        A. CHECKS ACCURACY OF CHFDV FUNCTION AND DERIVATIVE VALUES
  !           AGAINST EXACT VALUES.
  !        B. CHECKS THAT RETURNED VALUES OF NEXT SUM TO 10.
  !        C. CHECKS THAT FUNCTION VALUES FROM CHFEV AGREE WITH THOSE
  !           FROM CHFDV.
  !
  !
  !     FORTRAN INTRINSICS USED:  ABS, MAX, MIN.
  !     FORTRAN LIBRARY ROUTINES USED:  SQRT, (READ), (WRITE).
  !     SLATEC LIBRARY ROUTINES USED:  CHFDV, CHFEV, R1MACH, RAND.
  !     OTHER ROUTINES USED:  FDTRUE.
  !
  !***ROUTINES CALLED  CHFDV, CHFEV, FDTRUE, R1MACH, RAND
  !***REVISION HISTORY  (YYMMDD)
  !   820601  DATE WRITTEN
  !   820624  CONVERTED TO QUICK CHECK FOR SLATEC LIBRARY.
  !   820630  1. MODIFIED DEFINITIONS OF RELATIVE ERROR AND TEST
  !             TOLERANCES.
  !           2. VARIOUS IMPROVEMENTS TO OUTPUT FORMATS.
  !   820716  1. SET MACHEP VIA A CALL TO R1MACH.
  !           2. CHANGED FROM FORTLIB'S RANF TO SLATEC'S RAND.
  !   890629  1. Appended E0 to real constants to reduce S.P./D.P.
  !             differences.
  !           2. Other minor cosmetic changes.
  !   890831  Modified array declarations.  (WRB)
  !   890911  Removed unnecessary intrinsics.  (WRB)
  !   891004  Cosmetic changes to prologue.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  Revised prologue and improved some output formats.  (FNF)
  !           Also moved formats to end to be consistent with other PCHIP
  !           quick checks.
  !   900316  Additional minor cosmetic changes.  (FNF)
  !   900321  Made miscellaneous cosmetic changes.  (FNF)
  !   901130  Added 1P's to formats and revised some to reduce maximum
  !           line length.  (FNF)
  !   910708  Minor modifications in use of KPRINT.  (WRB)
  !   910801  Added EXTERNAL statement for RAND due to problem on IBM
  !           RS 6000.  (WRB)
  !***END PROLOGUE  EVCHCK
  !
  !  Declare arguments.
  !
  INTEGER Lout , Kprint , Npts
  REAL Xev(*) , Fev(*) , Dev(*) , Fev2(*)
  LOGICAL Fail
  !
  !  DECLARATIONS.
  !
  INTEGER i , ierr , iint , next(2) , next2(2) , nint
  REAL aed , aed2 , aedmax , aedmin , aef , aef2 , aefmax , aefmin ,&
    check(2) , checkf(2) , checkd(2) , d1 , d2 , dermax , dtrue , dx ,&
    eps1 , eps2 , f1 , f2 , fact , fermax , floord , floorf , four ,&
    ftrue , left(3) , machep , one , red , red2 , redmax , redmin , ref ,&
    ref2 , refmax , refmin , right(3) , small , ten , tol1 , tol2 , x1 ,&
    x2 , xadmax , xadmin , xafmax , xafmin , xrdmax , xrdmin , xrfmax ,&
    xrfmin , zero
  LOGICAL failoc , failnx
  !
  REAL R1MACH
  !       The following should stay REAL (no D.P. equivalent).
  REAL RAND
  EXTERNAL RAND
  !
  !  INITIALIZE.
  !
  DATA zero/0.E0/ , one/1.E0/ , four/4.E0/ , ten/10.E0/
  DATA small/1.0E-10/
  DATA nint/3/
  DATA left/ - 1.5E0 , 2.0E-10 , 1.0E0/
  DATA right/2.5E0 , 3.0E-10 , 1.0E+8/
  !
  !***FIRST EXECUTABLE STATEMENT  EVCHCK
  machep = R1MACH(4)
  eps1 = four*machep
  eps2 = ten*machep
  !
  Fail = .FALSE.
  !
  IF ( Kprint>=2 ) WRITE (Lout,99001)
  99001 FORMAT (//10X,'EVCHCK RESULTS'/10X,'--------------')
  !
  !  CYCLE OVER INTERVALS.
  !
  DO iint = 1 , nint
    x1 = left(iint)
    x2 = right(iint)
    !
    fact = MAX(SQRT(x2-x1),one)
    tol1 = eps1*fact
    tol2 = eps2*fact
    !
    !  COMPUTE AND PRINT ENDPOINT VALUES.
    !
    CALL FDTRUE(x1,f1,d1)
    CALL FDTRUE(x2,f2,d2)
    !
    IF ( Kprint>=3 ) THEN
      IF ( iint==1 ) WRITE (Lout,99002)
      !
      !  FORMATS.
      !
      99002     FORMAT (/10X,'CHFDV ACCURACY TEST')
      WRITE (Lout,'(/)')
      WRITE (Lout,99017) 'X1' , x1 , 'X2' , x2
      WRITE (Lout,99017) 'F1' , f1 , 'F2' , f2
      WRITE (Lout,99017) 'D1' , d1 , 'D2' , d2
    ENDIF
    !
    IF ( Kprint>=2 ) WRITE (Lout,99003) x1 , x2
    99003   FORMAT (/10X,'INTERVAL = (',1P,E12.5,',',E12.5,' ):')
    !
    !  COMPUTE FLOORS FOR RELATIVE ERRORS.
    !
    floorf = MAX(MIN(ABS(f1),ABS(f2)),small)
    floord = MAX(MIN(ABS(d1),ABS(d2)),small)
    !
    !  CHECK REPRODUCTION OF ENDPOINT VALUES.
    !
    Xev(1) = x1
    Xev(2) = x2
    !     -----------------------------------------------------------
    CALL CHFDV(x1,x2,f1,f2,d1,d2,2,Xev,checkf,checkd,next,ierr)
    !     -----------------------------------------------------------
    aef = checkf(1) - f1
    ref = RERR(aef,f1,floorf)
    aef2 = checkf(2) - f2
    ref2 = RERR(aef2,f2,floorf)
    aed = checkd(1) - d1
    red = RERR(aed,d1,floord)
    aed2 = checkd(2) - d2
    red2 = RERR(aed2,d2,floord)
    !
    failoc = MAX(ABS(ref),ABS(ref2),ABS(red),ABS(red2))>tol1
    Fail = Fail .OR. failoc
    !
    IF ( Kprint>=3 ) THEN
      WRITE (Lout,99004) next , aef , aef2 , aed , aed2
      99004     FORMAT (/' ERRORS AT ENDPOINTS:',40X,'(NEXT =',2I3,')'//1P,4X,'F1:',&
        E13.5,4X,'F2:',E13.5,4X,'D1:',E13.5,4X,'D2:',E13.5)
      WRITE (Lout,99005) ref , ref2 , red , red2
      99005     FORMAT (1P,4(7X,E13.5))
    ENDIF
    !
    IF ( failoc.AND.(Kprint>=2) ) WRITE (Lout,99006)
    99006   FORMAT (/' ***** CHFDV FAILED TO REPRODUCE ENDPOINT VALUES.')
    !
    !  CHFEV SHOULD AGREE EXACTLY WITH CHFDV.
    !                     -------
    !     --------------------------------------------------------------
    CALL CHFEV(x1,x2,f1,f2,d1,d2,2,Xev,check,next,ierr)
    !     --------------------------------------------------------------
    failoc = (check(1)/=checkf(1)) .OR. (check(2)/=checkf(2))
    Fail = Fail .OR. failoc
    !
    IF ( failoc.AND.(Kprint>=2) ) WRITE (Lout,99007)
    99007   FORMAT (/' ***** CHFEV DOES NOT AGREE WITH CHFDV AT ENDPOINTS.')
    !
    !  EVALUATE AT NPTS 'UNIFORMLY RANDOM' POINTS IN (X1,X2).
    !     THIS VERSION EXTENDS EVALUATION DOMAIN BY ADDING 4 SUBINTERVALS
    !     TO LEFT AND 6 TO RIGHT OF [X1,X2].
    !
    dx = (x2-x1)/(Npts-10)
    DO i = 1 , Npts
      Xev(i) = (x1+(i-5)*dx) + dx*RAND(zero)
    ENDDO
    !     --------------------------------------------------------
    CALL CHFDV(x1,x2,f1,f2,d1,d2,Npts,Xev,Fev,Dev,next,ierr)
    !     --------------------------------------------------------
    IF ( ierr/=0 ) THEN
      failoc = .TRUE.
      IF ( Kprint>=2 ) WRITE (Lout,99008) ierr
      99008     FORMAT (/' ***** ERROR ***** CHFDV RETURNED IERR =',I5)
    ELSE
      !
      !     CUMULATE LARGEST AND SMALLEST ERRORS FOR SUMMARY.
      !
      DO i = 1 , Npts
        CALL FDTRUE(Xev(i),ftrue,dtrue)
        aef = Fev(i) - ftrue
        ref = RERR(aef,ftrue,floorf)
        aed = Dev(i) - dtrue
        red = RERR(aed,dtrue,floord)
        !
        IF ( i==1 ) THEN
          !            INITIALIZE.
          aefmin = aef
          aefmax = aef
          aedmin = aed
          aedmax = aed
          refmin = ref
          refmax = ref
          redmin = red
          redmax = red
          xafmin = Xev(1)
          xafmax = Xev(1)
          xadmin = Xev(1)
          xadmax = Xev(1)
          xrfmin = Xev(1)
          xrfmax = Xev(1)
          xrdmin = Xev(1)
          xrdmax = Xev(1)
        ELSE
          !            SELECT.
          IF ( aef<aefmin ) THEN
            aefmin = aef
            xafmin = Xev(i)
          ELSEIF ( aef>aefmax ) THEN
            aefmax = aef
            xafmax = Xev(i)
          ENDIF
          IF ( aed<aedmin ) THEN
            aedmin = aed
            xadmin = Xev(i)
          ELSEIF ( aed>aedmax ) THEN
            aedmax = aed
            xadmax = Xev(i)
          ENDIF
          IF ( ref<refmin ) THEN
            refmin = ref
            xrfmin = Xev(i)
          ELSEIF ( ref>refmax ) THEN
            refmax = ref
            xrfmax = Xev(i)
          ENDIF
          IF ( red<redmin ) THEN
            redmin = red
            xrdmin = Xev(i)
          ELSEIF ( red>redmax ) THEN
            redmax = red
            xrdmax = Xev(i)
          ENDIF
        ENDIF
      ENDDO
      !
      fermax = MAX(ABS(refmax),ABS(refmin))
      dermax = MAX(ABS(redmax),ABS(redmin))
      !
      failnx = (next(1)+next(2))/=10
      failoc = failnx .OR. (MAX(fermax,dermax)>tol2)
    ENDIF
    Fail = Fail .OR. failoc
    !
    !  PRINT SUMMARY.
    !
    IF ( Kprint>=3 ) THEN
      WRITE (Lout,99009) Npts - 10 , next
      99009     FORMAT (/' ERRORS AT ',I5,' INTERIOR POINTS + 10 OUTSIDE:',15X,&
        '(NEXT =',2I3,')'//30X,'FUNCTION',17X,'DERIVATIVE'/15X,&
        2(11X,'ABS',9X,'REL'))
      !
      WRITE (Lout,99018) 'MIN' , aefmin , refmin , aedmin , redmin
      WRITE (Lout,99019) xafmin , xrfmin , xadmin , xrdmin
      WRITE (Lout,99018) 'MAX' , aefmax , refmax , aedmax , redmax
      WRITE (Lout,99019) xafmax , xrfmax , xadmax , xrdmax
    ENDIF
    !
    IF ( Kprint>=2 ) THEN
      IF ( failoc ) THEN
        IF ( fermax>tol2 ) WRITE (Lout,99020) 'F' , fermax , tol2
        IF ( dermax>tol2 ) WRITE (Lout,99020) 'D' , dermax , tol2
        IF ( failnx ) WRITE (Lout,99010) next
        99010       FORMAT (/' ***** REPORTED NEXT =',2I5,'   RATHER THAN    4    6')
      ELSE
        WRITE (Lout,99011)
        99011       FORMAT (/' CHFDV RESULTS OK.')
      ENDIF
    ENDIF
    !
    !  CHECK THAT CHFEV AGREES WITH CHFDV.
    !
    !     -----------------------------------------------------------------
    CALL CHFEV(x1,x2,f1,f2,d1,d2,Npts,Xev,Fev2,next2,ierr)
    !     -----------------------------------------------------------------
    IF ( ierr/=0 ) THEN
      failoc = .TRUE.
      IF ( Kprint>=2 ) WRITE (Lout,99012) ierr
      99012     FORMAT (/' ***** ERROR ***** CHFEV RETURNED IERR =',I5)
    ELSE
      aefmax = ABS(Fev2(1)-Fev(1))
      xafmax = Xev(1)
      DO i = 2 , Npts
        aef = ABS(Fev2(i)-Fev(i))
        IF ( aef>aefmax ) THEN
          aefmax = aef
          xafmax = Xev(i)
        ENDIF
      ENDDO
      failnx = (next2(1)/=next(1)) .OR. (next2(2)/=next(2))
      failoc = failnx .OR. (aefmax/=zero)
      IF ( Kprint>=2 ) THEN
        IF ( failoc ) THEN
          WRITE (Lout,99013)
          99013         FORMAT (/' ***** CHFEV DID NOT AGREE WITH CHFDV:')
          IF ( aefmax/=zero ) WRITE (Lout,99014) aefmax , xafmax
          99014         FORMAT (7X,'MAXIMUM DIFFERENCE ',1P,E12.5,'; OCCURRED AT X =',&
            E12.5)
          IF ( failnx ) WRITE (Lout,99015) next2 , next
          99015         FORMAT (7X,'REPORTED NEXT =',2I3,'   RATHER THAN ',2I3)
        ELSE
          WRITE (Lout,99016)
          99016         FORMAT (/' CHFEV AGREES WITH CHFDV.')
        ENDIF
      ENDIF
    ENDIF
    !
    Fail = Fail .OR. failoc
    !
    !  GO BACK FOR ANOTHER INTERVAL.
    !
  ENDDO
  !
  RETURN
  99017 FORMAT (10X,A2,' =',1P,E18.10,5X,A2,' =',E18.10)
  99018 FORMAT (/5X,A3,'IMUM ERROR:  ',1P,2E12.4,2X,2E12.4)
  99019 FORMAT (5X,'LOCATED AT X =  ',1P,2E12.4,2X,2E12.4)
  99020 FORMAT (/' ***** MAXIMUM RELATIVE ERROR IN ',A1,' =',1P,E12.5,','/17X,&
    'EXCEEDS TOLERANCE =',E12.5)
  !------------- LAST LINE OF EVCHCK FOLLOWS -----------------------------
CONTAINS

  !
  !  DEFINE RELATIVE ERROR WITH FLOOR.
  !
  REAL FUNCTION RERR(err,value,floor)
    REAL, INTENT(IN) :: err , value , floor
    RERR = err/MAX(ABS(value),floor)
  END FUNCTION RERR
END SUBROUTINE EVCHCK
