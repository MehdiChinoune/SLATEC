MODULE TEST33_MOD
  USE service, ONLY : SP, DP
  IMPLICIT NONE

CONTAINS
  !** DFDTRU
  SUBROUTINE DFDTRU(X,F,D)
    !> Compute exact function values for DEVCHK.
    !***
    ! **Library:**   SLATEC (PCHIP)
    !***
    ! **Type:**      DOUBLE PRECISION (FDTRUE-S, DFDTRU-D)
    !***
    ! **Keywords:**  PCHIP EVALUATOR QUICK CHECK
    !***
    ! **Author:**  Fritsch, F. N., (LLNL)
    !***
    ! **Description:**
    !
    !        COMPUTE EXACT FUNCTION VALUES IN DOUBLE PRECISION.
    !
    !                   F(X) = X*(X+1)*(X-2)
    !
    !***
    ! **Routines called:**  (NONE)

    !* REVISION HISTORY  (YYMMDD)
    !   820601  DATE WRITTEN
    !   890618  REVISION DATE from Version 3.2
    !   890706  Cosmetic changes to prologue.  (WRB)
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    !   900315  Revised prologue.  (FNF)
    !   900316  Deleted variables ONE and TWO.  (FNF)
    !   900321  Changed name of d.p. version from DFTRUE to DFDTRU.

    REAL(DP) :: X, F, D
    REAL(DP) :: fact1, fact2, xx
    !
    !* FIRST EXECUTABLE STATEMENT  DFDTRU
    xx = X
    fact1 = xx + 1
    fact2 = xx - 2
    F = xx*fact1*fact2
    D = fact1*fact2 + xx*(fact1+fact2)
    !
    !------------- LAST LINE OF DFDTRU FOLLOWS -----------------------------
  END SUBROUTINE DFDTRU
  !** DEVCHK
  SUBROUTINE DEVCHK(Lout,Kprint,Npts,Xev,Fev,Dev,Fev2,Fail)
    !> Test evaluation accuracy of DCHFDV and DCHFEV for DPCHQ1.
    !***
    ! **Library:**   SLATEC (PCHIP)
    !***
    ! **Type:**      DOUBLE PRECISION (EVCHCK-S, DEVCHK-D)
    !***
    ! **Keywords:**  PCHIP EVALUATOR QUICK CHECK
    !***
    ! **Author:**  Fritsch, F. N., (LLNL)
    !***
    ! **Description:**
    !
    ! -------- CODE TO TEST EVALUATION ACCURACY OF DCHFDV AND DCHFEV -------
    !
    !     USING FUNCTION AND DERIVATIVE VALUES FROM A CUBIC (COMPUTED IN
    !     DOUBLE PRECISION) AT NINT DIFFERENT (X1,X2) PAIRS:
    !     1. CHECKS THAT DCHFDV AND DCHFEV BOTH REPRODUCE ENDPOINT VALUES.
    !     2. EVALUATES AT NPTS POINTS, 10 OF WHICH ARE OUTSIDE THE INTERVAL
    !        AND:
    !        A. CHECKS ACCURACY OF DCHFDV FUNCTION AND DERIVATIVE VALUES
    !           AGAINST EXACT VALUES.
    !        B. CHECKS THAT RETURNED VALUES OF NEXT SUM TO 10.
    !        C. CHECKS THAT FUNCTION VALUES FROM DCHFEV AGREE WITH THOSE
    !           FROM DCHFDV.
    !
    !
    !     FORTRAN INTRINSICS USED:  ABS, MAX, MIN.
    !     FORTRAN LIBRARY ROUTINES USED:  SQRT, (READ), (WRITE).
    !     SLATEC LIBRARY ROUTINES USED:  DCHFDV, DCHFEV, D1MACH, RAND.
    !     OTHER ROUTINES USED:  DFDTRU.
    !
    !***
    ! **Routines called:**  D1MACH, DCHFDV, DCHFEV, DFDTRU, RAND

    !* REVISION HISTORY  (YYMMDD)
    !   820601  DATE WRITTEN
    !   820624  CONVERTED TO QUICK CHECK FOR SLATEC LIBRARY.
    !   820630  1. MODIFIED DEFINITIONS OF RELATIVE ERROR AND TEST
    !             TOLERANCES.
    !           2. VARIOUS IMPROVEMENTS TO OUTPUT FORMATS.
    !   820716  1. SET MACHEP VIA A CALL TO D1MACH.
    !           2. CHANGED FROM FORTLIB'S RANF TO SLATEC'S RAND.
    !   890628  1. Removed unnecessary IMPLICIT declaration.
    !           2. Removed unnecessary variable NEV.
    !           3. Other changes to reduce S.P./D.P. differences.
    !   890629  Added RERR to DOUBLE PRECISION declaration.
    !   890706  Cosmetic changes to prologue.  (WRB)
    !   890831  Modified array declarations.  (WRB)
    !   890911  Removed unnecessary intrinsics.  (WRB)
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    !   900315  Revised prologue and improved some output formats.  (FNF)
    !           Also moved formats to end to be consistent with other PCHIP
    !           quick checks.
    !   900316  Additional minor cosmetic changes.  (FNF)
    !   900321  Changed name of DFTRUE to DFDTRU and made additional minor
    !           cosmetic changes.  (FNF)
    !   901130  Added 1P's to formats and revised some to reduce maximum
    !           line length.  (FNF)
    !   910708  Minor modifications in use of KPRINT.  (WRB)
    !   910801  Added EXTERNAL statement for RAND due to problem on IBM
    !           RS 6000.  (WRB)
    !   910819  Changed argument to RAND function from a D.P. zero to a
    !           S.P. zero.  (WRB)
    USE slatec, ONLY : D1MACH, DCHFDV, DCHFEV, RAND
    !
    !  Declare arguments.
    !
    INTEGER :: Lout, Kprint, Npts
    REAL(DP) :: Xev(*), Fev(*), Dev(*), Fev2(*)
    LOGICAL :: Fail
    !
    !  DECLARATIONS.
    !
    INTEGER :: i, ierr, iintt, next(2), next2(2)
    REAL(DP) :: aed, aed2, aedmax, aedmin, aef, aef2, aefmax, aefmin, check(2), &
      checkf(2), checkd(2), d1, d2, dermax, dtrue, dx, eps1, eps2, f1, f2, fact, &
      fermax, floord, floorf, ftrue, machep, red, red2, redmax, redmin, &
      ref, ref2, refmax, refmin, tol1, tol2, x1, x2, xadmax, xadmin, &
      xafmax, xafmin, xrdmax, xrdmin, xrfmax, xrfmin
    LOGICAL :: failoc, failnx
    !
    !  INITIALIZE.
    !
    REAL(DP), PARAMETER :: zero = 0._DP, one = 1._DP, four = 4._DP, ten = 10._DP
    REAL(DP), PARAMETER :: small = 1.E-10_DP
    INTEGER, PARAMETER :: nintt = 3
    REAL(DP), PARAMETER :: left(3) = [ -1.5_DP, 2.E-10_DP, 1._DP ]
    REAL(DP), PARAMETER :: right(3) = [ 2.5_DP, 3.E-10_DP, 1.E8_DP ]
    !
    !* FIRST EXECUTABLE STATEMENT  DEVCHK
    machep = D1MACH(4)
    eps1 = four*machep
    eps2 = ten*machep
    !
    Fail = .FALSE.
    !
    IF( Kprint>=2 ) WRITE (Lout,99001)
    99001 FORMAT (//10X,'DEVCHK RESULTS'/10X,'--------------')
    !
    !  CYCLE OVER INTERVALS.
    !
    DO iintt = 1, nintt
      x1 = left(iintt)
      x2 = right(iintt)
      !
      fact = MAX(SQRT(x2-x1),one)
      tol1 = eps1*fact
      tol2 = eps2*fact
      !
      !  COMPUTE AND PRINT ENDPOINT VALUES.
      !
      CALL DFDTRU(x1,f1,d1)
      CALL DFDTRU(x2,f2,d2)
      !
      IF( Kprint>=3 ) THEN
        IF( iintt==1 ) WRITE (Lout,99002)
        !
        ! FORMATS.
        !
        99002 FORMAT (/10X,'DCHFDV ACCURACY TEST')
        WRITE (Lout,'(/)')
        WRITE (Lout,99017) 'X1', x1, 'X2', x2
        WRITE (Lout,99017) 'F1', f1, 'F2', f2
        WRITE (Lout,99017) 'D1', d1, 'D2', d2
      END IF
      !
      IF( Kprint>=2 ) WRITE (Lout,99003) x1, x2
      99003 FORMAT (/10X,'INTERVAL = (',1P,D12.5,',',D12.5,' ):')
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
      CALL DCHFDV(x1,x2,f1,f2,d1,d2,2,Xev,checkf,checkd,next,ierr)
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
      IF( Kprint>=3 ) THEN
        WRITE (Lout,99004) next, aef, aef2, aed, aed2
        99004 FORMAT (/' ERRORS AT ENDPOINTS:',40X,'(NEXT =',2I3,')'//1P,4X,'F1:',&
          D13.5,4X,'F2:',D13.5,4X,'D1:',D13.5,4X,'D2:',D13.5)
        WRITE (Lout,99005) ref, ref2, red, red2
        99005 FORMAT (1P,4(7X,D13.5))
      END IF
      !
      IF( failoc .AND. (Kprint>=2) ) WRITE (Lout,99006)
      99006 FORMAT (/' ***** DCHFDV FAILED TO REPRODUCE ENDPOINT VALUES.')
      !
      !  DCHFEV SHOULD AGREE EXACTLY WITH DCHFDV.
      !                     -------
      !     --------------------------------------------------------------
      CALL DCHFEV(x1,x2,f1,f2,d1,d2,2,Xev,check,next,ierr)
      !     --------------------------------------------------------------
      failoc = (check(1)/=checkf(1)) .OR. (check(2)/=checkf(2))
      Fail = Fail .OR. failoc
      !
      IF( failoc .AND. (Kprint>=2) ) WRITE (Lout,99007)
      99007 FORMAT (/' ***** DCHFEV DOES NOT AGREE WITH DCHFDV AT ENDPOINTS.')
      !
      !  EVALUATE AT NPTS 'UNIFORMLY RANDOM' POINTS IN (X1,X2).
      !     THIS VERSION EXTENDS EVALUATION DOMAIN BY ADDING 4 SUBINTERVALS
      !     TO LEFT AND 6 TO RIGHT OF [X1,X2].
      !
      dx = (x2-x1)/(Npts-10)
      DO i = 1, Npts
        Xev(i) = (x1+(i-5)*dx) + dx*RAND(0._SP)
      END DO
      !     --------------------------------------------------------
      CALL DCHFDV(x1,x2,f1,f2,d1,d2,Npts,Xev,Fev,Dev,next,ierr)
      !     --------------------------------------------------------
      IF( ierr/=0 ) THEN
        failoc = .TRUE.
        IF( Kprint>=2 ) WRITE (Lout,99008) ierr
        99008 FORMAT (/' ***** ERROR ***** DCHFDV RETURNED IERR =',I5)
      ELSE
        !
        !     CUMULATE LARGEST AND SMALLEST ERRORS FOR SUMMARY.
        !
        DO i = 1, Npts
          CALL DFDTRU(Xev(i),ftrue,dtrue)
          aef = Fev(i) - ftrue
          ref = RERR(aef,ftrue,floorf)
          aed = Dev(i) - dtrue
          red = RERR(aed,dtrue,floord)
          !
          IF( i==1 ) THEN
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
            IF( aef<aefmin ) THEN
              aefmin = aef
              xafmin = Xev(i)
            ELSEIF( aef>aefmax ) THEN
              aefmax = aef
              xafmax = Xev(i)
            END IF
            IF( aed<aedmin ) THEN
              aedmin = aed
              xadmin = Xev(i)
            ELSEIF( aed>aedmax ) THEN
              aedmax = aed
              xadmax = Xev(i)
            END IF
            IF( ref<refmin ) THEN
              refmin = ref
              xrfmin = Xev(i)
            ELSEIF( ref>refmax ) THEN
              refmax = ref
              xrfmax = Xev(i)
            END IF
            IF( red<redmin ) THEN
              redmin = red
              xrdmin = Xev(i)
            ELSEIF( red>redmax ) THEN
              redmax = red
              xrdmax = Xev(i)
            END IF
          END IF
        END DO
        !
        fermax = MAX(ABS(refmax),ABS(refmin))
        dermax = MAX(ABS(redmax),ABS(redmin))
        !
        failnx = (next(1)+next(2))/=10
        failoc = failnx .OR. (MAX(fermax,dermax)>tol2)
      END IF
      Fail = Fail .OR. failoc
      !
      !  PRINT SUMMARY.
      !
      IF( Kprint>=3 ) THEN
        WRITE (Lout,99009) Npts - 10, next
        99009 FORMAT (/' ERRORS AT ',I5,' INTERIOR POINTS + 10 OUTSIDE:',15X,&
          '(NEXT =',2I3,')'//30X,'FUNCTION',17X,'DERIVATIVE'/15X,&
          2(11X,'ABS',9X,'REL'))
        !
        WRITE (Lout,99018) 'MIN', aefmin, refmin, aedmin, redmin
        WRITE (Lout,99019) xafmin, xrfmin, xadmin, xrdmin
        WRITE (Lout,99018) 'MAX', aefmax, refmax, aedmax, redmax
        WRITE (Lout,99019) xafmax, xrfmax, xadmax, xrdmax
      END IF
      !
      IF( Kprint>=2 ) THEN
        IF( failoc ) THEN
          IF( fermax>tol2 ) WRITE (Lout,99020) 'F', fermax, tol2
          IF( dermax>tol2 ) WRITE (Lout,99020) 'D', dermax, tol2
          IF( failnx ) WRITE (Lout,99010) next
          99010 FORMAT (/' ***** REPORTED NEXT =',2I5,'   RATHER THAN    4    6')
        ELSE
          WRITE (Lout,99011)
          99011 FORMAT (/' DCHFDV RESULTS OK.')
        END IF
      END IF
      !
      !  CHECK THAT DCHFEV AGREES WITH DCHFDV.
      !
      !     -----------------------------------------------------------------
      CALL DCHFEV(x1,x2,f1,f2,d1,d2,Npts,Xev,Fev2,next2,ierr)
      !     -----------------------------------------------------------------
      IF( ierr/=0 ) THEN
        failoc = .TRUE.
        IF( Kprint>=2 ) WRITE (Lout,99012) ierr
        99012 FORMAT (/' ***** ERROR ***** DCHFEV RETURNED IERR =',I5)
      ELSE
        aefmax = ABS(Fev2(1)-Fev(1))
        xafmax = Xev(1)
        DO i = 2, Npts
          aef = ABS(Fev2(i)-Fev(i))
          IF( aef>aefmax ) THEN
            aefmax = aef
            xafmax = Xev(i)
          END IF
        END DO
        failnx = (next2(1)/=next(1)) .OR. (next2(2)/=next(2))
        failoc = failnx .OR. (aefmax/=zero)
        IF( Kprint>=2 ) THEN
          IF( failoc ) THEN
            WRITE (Lout,99013)
            99013 FORMAT (/' ***** DCHFEV DID NOT AGREE WITH DCHFDV:')
            IF( aefmax/=zero ) WRITE (Lout,99014) aefmax, xafmax
            99014 FORMAT (7X,'MAXIMUM DIFFERENCE ',1P,D12.5,'; OCCURRED AT X =',&
              D12.5)
            IF( failnx ) WRITE (Lout,99015) next2, next
            99015 FORMAT (7X,'REPORTED NEXT =',2I3,'   RATHER THAN ',2I3)
          ELSE
            WRITE (Lout,99016)
            99016 FORMAT (/' DCHFEV AGREES WITH DCHFDV.')
          END IF
        END IF
      END IF
      !
      Fail = Fail .OR. failoc
      !
      !  GO BACK FOR ANOTHER INTERVAL.
      !
    END DO
    !
    RETURN
    99017 FORMAT (10X,A2,' =',1P,D18.10,5X,A2,' =',D18.10)
    99018 FORMAT (/5X,A3,'IMUM ERROR:  ',1P,2D12.4,2X,2D12.4)
    99019 FORMAT (5X,'LOCATED AT X =  ',1P,2D12.4,2X,2D12.4)
    99020 FORMAT (/' ***** MAXIMUM RELATIVE ERROR IN ',A1,' =',1P,D12.5,','/17X,&
      'EXCEEDS TOLERANCE =',D12.5)
    !------------- LAST LINE OF DEVCHK FOLLOWS -----------------------------
  CONTAINS

    !
    !  DEFINE RELATIVE ERROR WITH FLOOR.
    !
    REAL(DP) FUNCTION RERR(err,value,floorr)
      REAL(DP), INTENT(IN) :: err, value, floorr
      RERR = err/MAX(ABS(value),floorr)
    END FUNCTION RERR
  END SUBROUTINE DEVCHK
  !** DEVERK
  SUBROUTINE DEVERK(Lout,Kprint,Fail)
    !> Test error returns from DPCHIP evaluators for DPCHQ1.
    !***
    ! **Library:**   SLATEC (PCHIP)
    !***
    ! **Type:**      DOUBLE PRECISION (EVERCK-S, DEVERK-D)
    !***
    ! **Keywords:**  PCHIP EVALUATOR QUICK CHECK
    !***
    ! **Author:**  Fritsch, F. N., (LLNL)
    !***
    ! **Description:**
    !
    ! --------- CODE TO TEST ERROR RETURNS FROM DPCHIP EVALUATORS. ---------
    !
    !
    !     FORTRAN LIBRARY ROUTINES USED:  (WRITE).
    !     SLATEC LIBRARY ROUTINES USED:  DCHFDV, DCHFEV, DPCHFD, DPCHFE,
    !                                    XERDMP, XGETF, XSETF.
    !     OTHER ROUTINES USED:  COMP.
    !
    !***
    ! **Routines called:**  COMP, DCHFDV, DCHFEV, DPCHFD, DPCHFE, XERDMP,
    !                    XGETF, XSETF

    !* REVISION HISTORY  (YYMMDD)
    !   820601  DATE WRITTEN
    !   820715  CONVERTED TO QUICK CHECK FOR SLATEC LIBRARY.
    !   890207  ADDED CALLS TO ERROR HANDLER.
    !   890316  Added call to XERDMP if KPRINT>2 (FNF).
    !   890706  Cosmetic changes to prologue.  (WRB)
    !   890911  Removed unnecessary intrinsics.  (WRB)
    !   891009  Removed unreferenced statement label.  (WRB)
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    !   900309  Added COMP to list of routines called.  (FNF)
    !   900315  Revised prologue and improved some output formats.  (FNF)
    !   900316  Deleted INCFD tests because some compilers object to them,
    !           and made additional minor cosmetic changes.  (FNF)
    !   900322  Made miscellaneous cosmetic changes.  (FNF)
    !   910708  Minor modifications in use of KPRINT.  (WRB)
    !   930504  Removed parens from constants in WRITE statements.  (FNF)
    USE slatec, ONLY : DCHFDV, DCHFEV, DPCHFD, DPCHFE, XERSVE, control_xer
    !
    !  Declare arguments.
    !
    INTEGER :: Lout, Kprint
    LOGICAL :: Fail
    !
    !  DECLARATIONS.
    !
    INTEGER :: i, ierr, kontrl, nerr, next(2)
    REAL(DP) :: d(10), dum(1), f(10), temp, x(10)
    LOGICAL :: skip
    !  INITIALIZE.
    INTEGER, PARAMETER :: N = 10
    !* FIRST EXECUTABLE STATEMENT  DEVERK
    nerr = 0
    !
    kontrl = control_xer
    IF( Kprint<=2 ) THEN
      control_xer = 0
    ELSE
      control_xer = 1
    END IF
    !
    IF( Kprint>=3 ) WRITE (Lout,99001)
    !
    ! FORMATS.
    !
    99001 FORMAT ('1'//10X,'TEST ERROR RETURNS')
    IF( Kprint>=2 ) WRITE (Lout,99002)
    99002 FORMAT (//10X,'DEVERK RESULTS'/10X,'--------------')
    !
    !  FIRST, TEST DCHFEV AND DCHFDV.
    !
    IF( Kprint>=3 ) WRITE (Lout,99005) -1
    CALL DCHFEV(0._DP,1._DP,3._DP,7._DP,3._DP,6._DP,0,dum,dum,next,ierr)
    IF( ierr /= -1 ) nerr = nerr + 1
    !
    IF( Kprint>=3 ) WRITE (Lout,99005) -2
    CALL DCHFEV(1._DP,1._DP,3._DP,7._DP,3._DP,6._DP,1,dum,dum,next,ierr)
    IF( ierr /= -2 ) nerr = nerr + 1
    !
    IF( Kprint>=3 ) WRITE (Lout,99005) -1
    CALL DCHFDV(0._DP,1._DP,3._DP,7._DP,3._DP,6._DP,0,dum,dum,dum,next,ierr)
    IF( ierr /= -1 ) nerr = nerr + 1
    !
    IF( Kprint>=3 ) WRITE (Lout,99005) -2
    CALL DCHFDV(1._DP,1._DP,3._DP,7._DP,3._DP,6._DP,1,dum,dum,dum,next,ierr)
    IF( ierr /= -2 ) nerr = nerr + 1
    !
    !  SET UP PCH DEFINITION.
    !
    DO i = 1, N
      x(i) = i
      f(i) = i + 2
      d(i) = 1._DP
    END DO
    !
    !  SWAP POINTS 4 AND 7, SO X-ARRAY IS OUT OF ORDER.
    !
    temp = x(4)
    x(4) = x(7)
    x(7) = temp
    !
    !  NOW, TEST DPCHFE AND DPCHFD.
    !
    IF( Kprint>=3 ) WRITE (Lout,99005) -1
    skip = .FALSE.
    CALL DPCHFE(1,x,f,d,1,skip,0,dum,dum,ierr)
    IF( ierr /= -1 ) nerr = nerr + 1
    !
    IF( Kprint>=3 ) WRITE (Lout,99005) -3
    skip = .FALSE.
    CALL DPCHFE(N,x,f,d,1,skip,0,dum,dum,ierr)
    IF( ierr /= -3 ) nerr = nerr + 1
    !
    IF( Kprint>=3 ) WRITE (Lout,99005) -4
    skip = .TRUE.
    CALL DPCHFE(N,x,f,d,1,skip,0,dum,dum,ierr)
    IF( ierr /= -4 ) nerr = nerr + 1
    !
    IF( Kprint>=3 ) WRITE (Lout,99005) -1
    skip = .FALSE.
    CALL DPCHFD(1,x,f,d,1,skip,0,dum,dum,dum,ierr)
    IF( ierr /= -1 ) nerr = nerr + 1
    !
    IF( Kprint>=3 ) WRITE (Lout,99005) -3
    skip = .FALSE.
    CALL DPCHFD(N,x,f,d,1,skip,0,dum,dum,dum,ierr)
    IF( ierr /= -3 ) nerr = nerr + 1
    !
    IF( Kprint>=3 ) WRITE (Lout,99005) -4
    skip = .TRUE.
    CALL DPCHFD(N,x,f,d,1,skip,0,dum,dum,dum,ierr)
    IF( ierr /= -4 ) nerr = nerr + 1
    !
    !  SUMMARIZE RESULTS.
    !
    IF( Kprint>2 ) CALL XERSVE(' ',' ',0,0,0,i)
    IF( nerr==0 ) THEN
      Fail = .FALSE.
      IF( Kprint>=2 ) WRITE (Lout,99003)
      99003 FORMAT (/' ALL ERROR RETURNS OK.')
    ELSE
      Fail = .TRUE.
      IF( Kprint>=2 ) WRITE (Lout,99004) nerr
      99004 FORMAT (//' ***** TROUBLE IN DEVERK *****'//5X,I5,&
        ' TESTS FAILED TO GIVE EXPECTED RESULTS.')
    END IF
    !
    !  TERMINATE.
    !
    control_xer = kontrl
    RETURN
    99005 FORMAT (/' THIS CALL SHOULD RETURN IERR =',I3)
    !------------- LAST LINE OF DEVERK FOLLOWS -----------------------------
  END SUBROUTINE DEVERK
  !** DEVPCK
  SUBROUTINE DEVPCK(Lout,Kprint,X,Y,F,Fx,Fy,Xe,Ye,Fe,De,Fe2,Fail)
    !> Test usage of increment argument in DPCHFD and DPCHFE for
    !            DPCHQ1.
    !***
    ! **Library:**   SLATEC (PCHIP)
    !***
    ! **Type:**      DOUBLE PRECISION (EVPCCK-S, DEVPCK-D)
    !***
    ! **Keywords:**  PCHIP EVALUATOR QUICK CHECK
    !***
    ! **Author:**  Fritsch, F. N., (LLNL)
    !***
    ! **Description:**
    !
    ! ---- CODE TO TEST USAGE OF INCREMENT ARGUMENT IN DPCHFD AND DPCHFE ---
    !
    !     EVALUATES A BICUBIC FUNCTION AND ITS FIRST PARTIAL DERIVATIVES
    !     ON A 4X6 MESH CONTAINED IN A 10X10 ARRAY.
    !
    !     INTERPOLATION OF THESE DATA ALONG MESH LINES IN EITHER DIMENSION
    !     SHOULD AGREE WITH CORRECT FUNCTION WITHIN ROUNDOFF ERROR.
    !
    !     ARRAYS ARE ARGUMENTS ONLY TO ALLOW SHARING STORAGE WITH OTHER
    !     TEST ROUTINES.
    !
    !     NOTE:  RUN WITH KPRINT=4 FOR FULL GORY DETAILS (10 PAGES WORTH).
    !
    !
    !     FORTRAN INTRINSICS USED:  ABS.
    !     FORTRAN LIBRARY ROUTINES USED:  (WRITE).
    !     SLATEC LIBRARY ROUTINES USED:  DPCHFD, DPCHFE, D1MACH.
    !
    !***
    ! **Routines called:**  D1MACH, DPCHFD, DPCHFE

    !* REVISION HISTORY  (YYMMDD)
    !   820601  DATE WRITTEN
    !   820714  CONVERTED TO QUICK CHECK FOR SLATEC LIBRARY.
    !   820715  1. CORRECTED SOME FORMATS.
    !           2. ADDED CALL TO D1MACH TO SET MACHEP.
    !   890406  1. Modified to make sure final elements of X and XE
    !             agree, to avoid possible failure due to roundoff
    !             error.
    !           2. Added printout of TOL in case of failure.
    !           3. Removed unnecessary IMPLICIT declaration.
    !           4. Corrected a few S.P. constants to D.P.
    !           5. Minor cosmetic changes.
    !   890706  Cosmetic changes to prologue.  (WRB)
    !   890911  Removed unnecessary intrinsics.  (WRB)
    !   891004  Cosmetic changes to prologue.  (WRB)
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    !   900315  Revised prologue and improved some output formats.  (FNF)
    !   900316  Additional minor cosmetic changes.  (FNF)
    !   900321  Made miscellaneous cosmetic changes.  (FNF)
    !   901130  Made many changes to output:  (FNF)
    !           1. Reduced amount of output for KPRINT=3.  (Now need to
    !              use KPRINT=4 for full output.)
    !           2. Added 1P's to formats and revised some to reduce maximum
    !              line length.
    !   910708  Minor modifications in use of KPRINT.  (WRB)
    !   930317  Improved output formats.  (FNF)
    USE slatec, ONLY : D1MACH, DPCHFD, DPCHFE
    !
    !  Declare arguments.
    !
    INTEGER :: Lout, Kprint
    LOGICAL :: Fail
    REAL(DP) :: X(10), Y(10), F(10,10), Fx(10,10), Fy(10,10), Xe(51), &
      Ye(51), Fe(51), De(51), Fe2(51)
    !
    !  DECLARATIONS.
    !
    INTEGER :: i, ier2, ierr, inc, j, k, nerr
    LOGICAL :: faild, faile, failoc, skip
    REAL(DP) :: dermax, derr, dtrue, dx, fdiff, fdifmx, fermax, &
      ferr, ftrue, machep, tol, pdermx, pdifmx, pfermx
    !
    INTEGER, PARAMETER :: nmax = 10, nx = 4, ny = 6
    INTEGER, PARAMETER :: ne = 51
    REAL(DP), PARAMETER :: zero = 0._DP
    !
    !  INITIALIZE.
    !
    !* FIRST EXECUTABLE STATEMENT  DEVPCK
    machep = D1MACH(4)
    !       Following tolerance is looser than S.P. version to avoid
    !       spurious failures on some systems.
    tol = 25._DP*machep
    !
    Fail = .FALSE.
    !
    !  SET UP 4-BY-6 MESH IN A 10-BY-10 ARRAY:
    !     X =  0.25(0.25)1.   ;
    !     Y = -0.75(0.5 )1.75 .
    !
    DO i = 1, nx - 1
      X(i) = 0.25_DP*i
    END DO
    X(nx) = 1._DP
    DO j = 1, ny
      Y(j) = 0.5_DP*j - 1.25_DP
      DO i = 1, nx
        F(i,j) = FCN(X(i),Y(j))
        Fx(i,j) = DFDX(X(i),Y(j))
        Fy(i,j) = DFDY(X(i),Y(j))
      END DO
    END DO
    !
    !  SET UP EVALUATION POINTS:
    !     XE =  0.(0.02)1. ;
    !     YE = -2.(0.08)2. .
    !
    dx = 1._DP/(ne-1)
    DO k = 1, ne - 1
      Xe(k) = dx*(k-1)
      Ye(k) = 4._DP*Xe(k) - 2._DP
    END DO
    Xe(ne) = 1._DP
    Ye(ne) = 2._DP
    !
    IF( Kprint>=3 ) WRITE (Lout,99001)
    !
    ! FORMATS.
    !
    99001 FORMAT ('1'//10X,'TEST DPCHFE AND DPCHFD')
    IF( Kprint>=2 ) WRITE (Lout,99002)
    99002 FORMAT (//10X,'DEVPCK RESULTS'/10X,'--------------')
    !
    !  EVALUATE ON HORIZONTAL MESH LINES (Y FIXED, X RUNNING) ..............
    !
    nerr = 0
    inc = 1
    skip = .FALSE.
    DO j = 1, ny
      !        --------------------------------------------------------------
      CALL DPCHFD(nx,X,F(1,j),Fx(1,j),inc,skip,ne,Xe,Fe,De,ierr)
      !        --------------------------------------------------------------
      IF( Kprint>=3 ) WRITE (Lout,99003) inc, 'J', j, 'Y', Y(j), ierr
      IF( ierr<0 ) THEN
        !
        failoc = .TRUE.
        IF( Kprint>=2 ) WRITE (Lout,99011) ierr
      ELSE
        IF( Kprint>3 ) WRITE (Lout,99004) 'X'
        !
        !        DPCHFE SHOULD AGREE EXACTLY WITH DPCHFD.
        !
        !        -----------------------------------------------------------
        CALL DPCHFE(nx,X,F(1,j),Fx(1,j),inc,skip,ne,Xe,Fe2,ier2)
        !        -----------------------------------------------------------
        !
        DO k = 1, ne
          ftrue = FCN(Xe(k),Y(j))
          ferr = Fe(k) - ftrue
          dtrue = DFDX(Xe(k),Y(j))
          derr = De(k) - dtrue
          IF( Kprint>3 ) WRITE (Lout,99005) Xe(k), ftrue, Fe(k), ferr, &
            dtrue, De(k), derr
          IF( k==1 ) THEN
            !              INITIALIZE.
            fermax = ABS(ferr)
            pfermx = Xe(1)
            dermax = ABS(derr)
            pdermx = Xe(1)
            fdifmx = ABS(Fe2(1)-Fe(1))
            pdifmx = Xe(1)
          ELSE
            !              SELECT.
            ferr = ABS(ferr)
            IF( ferr>fermax ) THEN
              fermax = ferr
              pfermx = Xe(k)
            END IF
            derr = ABS(derr)
            IF( derr>dermax ) THEN
              dermax = derr
              pdermx = Xe(k)
            END IF
            fdiff = ABS(Fe2(k)-Fe(k))
            IF( fdiff>fdifmx ) THEN
              fdifmx = fdiff
              pdifmx = Xe(k)
            END IF
          END IF
        END DO
        !
        faild = (fermax>tol) .OR. (dermax>tol)
        faile = fdifmx/=zero
        failoc = faild .OR. faile .OR. (ierr/=13) .OR. (ier2/=ierr)
        !
        IF( failoc .AND. (Kprint>=2) ) WRITE (Lout,99006) 'J', j, 'Y', Y(j)
        !
        IF( (Kprint>=3) .OR. (faild .AND. (Kprint==2)) ) WRITE (Lout,99007)&
          fermax, pfermx, dermax, pdermx
        IF( faild .AND. (Kprint>=2) ) WRITE (Lout,99010) tol
        !
        IF( (Kprint>=3) .OR. (faile .AND. (Kprint==2)) ) WRITE (Lout,99008)&
          fdifmx, pdifmx
        !
        IF( (ierr/=13) .AND. (Kprint>=2) ) WRITE (Lout,99009) 'D', ierr, 13
        !
        IF( (ier2/=ierr) .AND. (Kprint>=2) ) WRITE (Lout,99009) 'E', ier2, ierr
      END IF
      !
      IF( failoc ) nerr = nerr + 1
      Fail = Fail .OR. failoc
    END DO
    !
    IF( Kprint>=2 ) THEN
      IF( nerr>0 ) THEN
        WRITE (Lout,99012) nerr, 'J'
      ELSE
        WRITE (Lout,99013) 'J'
      END IF
    END IF
    !
    !  EVALUATE ON VERTICAL MESH LINES (X FIXED, Y RUNNING) ................
    !
    nerr = 0
    inc = nmax
    skip = .FALSE.
    DO i = 1, nx
      !        --------------------------------------------------------------
      CALL DPCHFD(ny,Y,F(i,1),Fy(i,1),inc,skip,ne,Ye,Fe,De,ierr)
      !        --------------------------------------------------------------
      IF( Kprint>=3 ) WRITE (Lout,99003) inc, 'I', i, 'X', X(i), ierr
      IF( ierr<0 ) THEN
        !
        failoc = .TRUE.
        IF( Kprint>=2 ) WRITE (Lout,99011) ierr
      ELSE
        IF( Kprint>3 ) WRITE (Lout,99004) 'Y'
        !
        !        DPCHFE SHOULD AGREE EXACTLY WITH DPCHFD.
        !
        !        -----------------------------------------------------------
        CALL DPCHFE(ny,Y,F(i,1),Fy(i,1),inc,skip,ne,Ye,Fe2,ier2)
        !        -----------------------------------------------------------
        !
        DO k = 1, ne
          ftrue = FCN(X(i),Ye(k))
          ferr = Fe(k) - ftrue
          dtrue = DFDY(X(i),Ye(k))
          derr = De(k) - dtrue
          IF( Kprint>3 ) WRITE (Lout,99005) Ye(k), ftrue, Fe(k), ferr, &
            dtrue, De(k), derr
          IF( k==1 ) THEN
            !              INITIALIZE.
            fermax = ABS(ferr)
            pfermx = Ye(1)
            dermax = ABS(derr)
            pdermx = Ye(1)
            fdifmx = ABS(Fe2(1)-Fe(1))
            pdifmx = Ye(1)
          ELSE
            !              SELECT.
            ferr = ABS(ferr)
            IF( ferr>fermax ) THEN
              fermax = ferr
              pfermx = Ye(k)
            END IF
            derr = ABS(derr)
            IF( derr>dermax ) THEN
              dermax = derr
              pdermx = Ye(k)
            END IF
            fdiff = ABS(Fe2(k)-Fe(k))
            IF( fdiff>fdifmx ) THEN
              fdifmx = fdiff
              pdifmx = Ye(k)
            END IF
          END IF
        END DO
        !
        faild = (fermax>tol) .OR. (dermax>tol)
        faile = fdifmx/=zero
        failoc = faild .OR. faile .OR. (ierr/=20) .OR. (ier2/=ierr)
        !
        IF( failoc .AND. (Kprint>=2) ) WRITE (Lout,99006) 'I', i, 'X', X(i)
        !
        IF( (Kprint>=3) .OR. (faild .AND. (Kprint==2)) ) WRITE (Lout,99007)&
          fermax, pfermx, dermax, pdermx
        IF( faild .AND. (Kprint>=2) ) WRITE (Lout,99010) tol
        !
        IF( (Kprint>=3) .OR. (faile .AND. (Kprint==2)) ) WRITE (Lout,99008)&
          fdifmx, pdifmx
        !
        IF( (ierr/=20) .AND. (Kprint>=2) ) WRITE (Lout,99009) 'D', ierr, 20
        !
        IF( (ier2/=ierr) .AND. (Kprint>=2) ) WRITE (Lout,99009) 'E', ier2, ierr
      END IF
      !
      IF( failoc ) nerr = nerr + 1
      Fail = Fail .OR. failoc
    END DO
    !
    IF( Kprint>=2 ) THEN
      IF( nerr>0 ) THEN
        WRITE (Lout,99012) nerr, 'I'
      ELSE
        WRITE (Lout,99013) 'I'
      END IF
    END IF
    !
    !  TERMINATE.
    !
    RETURN
    99003 FORMAT (//20X,'DPCHFD INCREMENT TEST -- INCFD = ',I2/15X,'ON ',A1,&
      '-LINE ',I2,',  ',A1,' =',F8.4,'  --  IERR =',I3)
    99004 FORMAT (/3X,A1,'E',10X,'F',8X,'FE',9X,'DIFF',13X,'D',8X,'DE',9X,'DIFF')
    99005 FORMAT (F7.2,2(2X,2F10.5,1P,D15.5,0P))
    99006 FORMAT (/' ***** DPCHFD AND/OR DPCHFE FAILED ON ',A1,'-LINE ',I1,',  ',A1,&
      ' =',F8.4)
    99007 FORMAT (/19X,'  MAXIMUM ERROR IN FUNCTION =',1P,1P,D13.5,0P,' (AT',F6.2,&
      '),'/33X,'IN DERIVATIVE =',1P,D13.5,0P,' (AT',F6.2,').')
    99008 FORMAT ('  MAXIMUM DIFFERENCE BETWEEN DPCHFE AND DPCHFD =',1P,D13.5,0P,&
      ' (AT',F6.2,').')
    99009 FORMAT (/'  DPCHF',A1,' RETURNED IERR = ',I2,' INSTEAD OF ',I2)
    99010 FORMAT ('  *** BOTH SHOULD BE <= TOL =',1P,D12.5,' ***')
    99011 FORMAT (//' ***** ERROR ***** DPCHFD RETURNED IERR =',I5//)
    99012 FORMAT (//' ***** ERROR ***** DPCHFD AND/OR DPCHFE FAILED ON',I2,1X,A1,&
      '-LINES.'//)
    99013 FORMAT (/' DPCHFD AND DPCHFE OK ON ',A1,'-LINES.')
    !------------- LAST LINE OF DEVPCK FOLLOWS -----------------------------
  CONTAINS

    !
    !  DEFINE TEST FUNCTION AND DERIVATIVES.
    !
    REAL(DP) FUNCTION FCN(ax,ay)
      REAL(DP) ax, ay
      FCN = ax*(ay*ay)*(ax*ax+1._SP)
    END FUNCTION FCN
    REAL(DP) FUNCTION DFDX(ax,ay)
      REAL(DP) ax, ay
      DFDX = (ay*ay)*(3._SP*ax*ax+1._SP)
    END FUNCTION DFDX
    REAL(DP) FUNCTION DFDY(ax,ay)
      REAL(DP) ax, ay
      DFDY = 2._SP*ax*ay*(ax*ax+1._SP)
    END FUNCTION DFDY
  END SUBROUTINE DEVPCK
  !** DPCHQ1
  SUBROUTINE DPCHQ1(Lun,Kprint,Ipass)
    !> Test the PCHIP evaluators DCHFDV, DCHFEV, DPCHFD, DPCHFE.
    !***
    ! **Library:**   SLATEC (PCHIP)
    !***
    ! **Type:**      DOUBLE PRECISION (PCHQK1-S, DPCHQ1-D)
    !***
    ! **Keywords:**  PCHIP EVALUATOR QUICK CHECK
    !***
    ! **Author:**  Fritsch, F. N., (LLNL)
    !***
    ! **Description:**
    !
    !             DPCHIP QUICK CHECK NUMBER 1
    !
    !     TESTS THE EVALUATORS:  DCHFDV, DCHFEV, DPCHFD, DPCHFE.
    !- Usage:
    !
    !        INTEGER  LUN, KPRINT, IPASS
    !
    !        CALL DPCHQ1 (LUN, KPRINT, IPASS)
    !
    !- Arguments:
    !
    !     LUN   :IN  is the unit number to which output is to be written.
    !
    !     KPRINT:IN  controls the amount of output, as specified in the
    !                SLATEC Guidelines.
    !
    !     IPASS:OUT  will contain a pass/fail flag.  IPASS=1 is good.
    !                IPASS=0 indicates one or more tests failed.
    !
    !- Description:
    !
    !   This routine carries out three tests of the PCH evaluators:
    !     DEVCHK tests the single-cubic evaluators.
    !     DEVPCK tests the full PCH evaluators.
    !     DEVERK exercises the error returns in all evaluators.
    !
    !***
    ! **Routines called:**  DEVCHK, DEVERK, DEVPCK

    !* REVISION HISTORY  (YYMMDD)
    !   820601  DATE WRITTEN
    !   890306  Changed IPASS to the more accurate name IFAIL.  (FNF)
    !   890307  Removed conditional on call to DEVERK.
    !   890706  Cosmetic changes to prologue.  (WRB)
    !   891004  Correction in prologue.  (WRB)
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    !   900309  Added DEVERK to list of routines called.  (FNF)
    !   900314  Improved some output formats.
    !   900315  Revised prologue and improved some output formats.  (FNF)
    !   900316  Additional minor cosmetic changes.  (FNF)
    !   900321  Removed IFAIL from call sequence for SLATEC standards and
    !           made miscellaneous cosmetic changes.  (FNF)
    !   930317  Improved output formats.  (FNF)

    !
    !  Declare arguments.
    !
    INTEGER :: Lun, Kprint, Ipass
    !
    !  DECLARE LOCAL VARIABLES.
    !
    INTEGER :: i1, i2, i3, i4, i5, i6, i7, i8, i9, ifail, npts
    REAL(DP) :: work(4000)
    LOGICAL :: fail
    !
    !* FIRST EXECUTABLE STATEMENT  DPCHQ1
    IF( Kprint>=2 ) WRITE (Lun,99001) Kprint
    !
    ! FORMATS.
    !
    99001 FORMAT ('1'/' ------------ DPCHIP QUICK CHECK OUTPUT',' ------------'//&
      20X,'( KPRINT =',I2,' )')
    !
    !  TEST DCHFDV AND DCHFEV.
    !
    ifail = 0
    npts = 1000
    i1 = 1 + npts
    i2 = i1 + npts
    i3 = i2 + npts
    CALL DEVCHK(Lun,Kprint,npts,work(1),work(i1),work(i2),work(i3),fail)
    IF( fail ) ifail = ifail + 1
    !
    !  TEST DPCHFD AND DPCHFE.
    !
    i1 = 1 + 10
    i2 = i1 + 10
    i3 = i2 + 100
    i4 = i3 + 100
    i5 = i4 + 100
    i6 = i5 + 51
    i7 = i6 + 51
    i8 = i7 + 51
    i9 = i8 + 51
    CALL DEVPCK(Lun,Kprint,work(1),work(i1),work(i2),work(i3),work(i4),&
      work(i5),work(i6),work(i7),work(i8),work(i9),fail)
    IF( fail ) ifail = ifail + 2
    !
    !  TEST ERROR RETURNS.
    !
    CALL DEVERK(Lun,Kprint,fail)
    IF( fail ) ifail = ifail + 4
    !
    !  PRINT SUMMARY AND TERMINATE.
    !     At this point, IFAIL has the following value:
    !        IFAIL = 0  IF ALL TESTS PASSED.
    !        IFAIL BETWEEN 1 AND 7 IS THE SUM OF:
    !           IFAIL=1  IF SINGLE CUBIC  TEST FAILED. (SEE DEVCHK OUTPUT.)
    !           IFAIL=2  IF DPCHFD/DPCHFE TEST FAILED. (SEE DEVPCK OUTPUT.)
    !           IFAIL=4  IF ERROR RETURN  TEST FAILED. (SEE DEVERK OUTPUT.)
    !
    IF( (Kprint>=2) .AND. (ifail/=0) ) WRITE (Lun,99002) ifail
    99002 FORMAT (/' *** TROUBLE ***',I5,' EVALUATION TESTS FAILED.')
    !
    IF( ifail==0 ) THEN
      Ipass = 1
      IF( Kprint>=2 ) WRITE (Lun,99003)
      99003 FORMAT (/' ------------ DPCHIP PASSED  ALL EVALUATION TESTS',&
        ' ------------')
    ELSE
      Ipass = 0
      IF( Kprint>=1 ) WRITE (Lun,99004)
      99004 FORMAT (/' ************ DPCHIP FAILED SOME EVALUATION TESTS',&
        ' ************')
    END IF
    !
    RETURN
    !------------- LAST LINE OF DPCHQ1 FOLLOWS -----------------------------
  END SUBROUTINE DPCHQ1
  !** DPCHQ2
  SUBROUTINE DPCHQ2(Lun,Kprint,Ipass)
    !> Test the PCHIP integrators DPCHIA and DPCHID.
    !***
    ! **Library:**   SLATEC (PCHIP)
    !***
    ! **Type:**      DOUBLE PRECISION (PCHQK2-S, DPCHQ2-D)
    !***
    ! **Keywords:**  PCHIP INTEGRATOR QUICK CHECK
    !***
    ! **Author:**  Fritsch, F. N., (LLNL)
    !***
    ! **Description:**
    !
    !             DPCHIP QUICK CHECK NUMBER 2
    !
    !     TESTS THE INTEGRATORS:  DPCHIA, DPCHID.
    !- Usage:
    !
    !        INTEGER  LUN, KPRINT, IPASS
    !
    !        CALL DPCHQ2 (LUN, KPRINT, IPASS)
    !
    !- Arguments:
    !
    !     LUN   :IN  is the unit number to which output is to be written.
    !
    !     KPRINT:IN  controls the amount of output, as specified in the
    !                SLATEC Guidelines.
    !
    !     IPASS:OUT  will contain a pass/fail flag.  IPASS=1 is good.
    !                IPASS=0 indicates one or more tests failed.
    !
    !- Description:
    !
    !   This routine constructs data from a cubic, integrates it with DPCHIA
    !   and compares the results with the correct answer.
    !   Since DPCHIA calls DPCHID, this tests both integrators.
    !
    !***
    ! **Routines called:**  D1MACH, DPCHIA

    !* REVISION HISTORY  (YYMMDD)
    !   820601  DATE WRITTEN
    !   890306  Changed IPASS to the more accurate name IFAIL.  (FNF)
    !   890316  1. Removed IMPLICIT statement.                  (FNF)
    !           2. Eliminated unnecessary variable N1.          (FNF)
    !           3. Miscellaneous cosmetic changes.              (FNF)
    !   891004  Cosmetic changes to prologue.  (WRB)
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    !   900314  Improved some output formats.  (FNF)
    !   900315  Revised prologue and improved some output formats.  (FNF)
    !   900316  Additional minor cosmetic changes.  (FNF)
    !   900321  Removed IFAIL from call sequence for SLATEC standards and
    !           made miscellaneous cosmetic changes.  (FNF)
    !   900323  Corrected list of routines called.  (FNF)
    !   901130  Added 1P's to formats; changed to allow KPRINT>3.  (FNF)
    !   910708  Minor modifications in use of KPRINT.  (WRB)
    !   930317  Improved output formats.  (FNF)
    USE slatec, ONLY : D1MACH, DPCHIA
    !
    !  Declare arguments.
    !
    INTEGER :: Lun, Kprint, Ipass
    !
    !  DECLARE VARIABLES.
    !
    INTEGER :: i, ierr, ifail
    REAL(DP) :: calc, d(7), errmax, error, f(7), machep, tol, true
    LOGICAL :: fail, skip
    !
    !  INITIALIZE.
    !
    REAL(DP), PARAMETER :: thrqtr = 0.75_DP, one = 1._DP, two = 2._DP, three = 3._DP
    INTEGER, PARAMETER :: n = 7
    REAL(DP), PARAMETER :: x(7) = [ -4._DP, -2._DP, -0.9_DP, 0._DP, 0.9_DP, 2._DP, 4._DP ]
    INTEGER, PARAMETER :: npairs = 17
    REAL(DP), PARAMETER :: a(17) = [ -3._DP, 3._DP, -0.5_DP, -0.5_DP, -0.5_DP, -4._DP, &
      -4._DP, 3._DP, -5._DP, -5._DP, -6._DP, 6._DP, -1.5_DP, -1.5_DP, -3._DP, 3._DP, 0.5_DP ]
    REAL(DP), PARAMETER :: b(17) =  [ 3._DP, -3._DP, 1._DP, 2._DP, 5._DP, -0.5_DP, &
      4._DP, 5._DP, -3._DP, 5._DP, -5._DP, 5._DP, -0.5_DP, -1._DP, -2.5_DP, 3.5_DP, 0.5_DP ]
    INTEGER, PARAMETER :: ierexp(17) = [ 0, 0, 0, 0, 2, 0, 0, 2, 1, 3, 3, 3, 0, &
      0, 0, 0, 0 ]
    !
    !  SET PASS/FAIL TOLERANCE.
    !
    !* FIRST EXECUTABLE STATEMENT  DPCHQ2
    machep = D1MACH(4)
    tol = 100._DP*machep
    !
    !  SET UP PCH FUNCTION DEFINITION.
    !
    DO i = 1, n
      f(i) = FCN(x(i))
      d(i) = DERIV(x(i))
    END DO
    !
    IF( Kprint>=3 ) WRITE (Lun,99001)
    !
    ! FORMATS.
    !
    99001 FORMAT ('1'//10X,'TEST DPCHIP INTEGRATORS')
    IF( Kprint>=2 ) WRITE (Lun,99002)
    99002 FORMAT (//10X,'DPCHQ2 RESULTS'/10X,'--------------')
    IF( Kprint>=3 ) WRITE (Lun,99003) (x(i),f(i),d(i),i=1,n)
    99003 FORMAT (//5X,'DATA:'//11X,'X',9X,'F',9X,'D'/(5X,3F10.3))
    !
    !  LOOP OVER (A,B)-PAIRS.
    !
    IF( Kprint>=3 ) WRITE (Lun,99004)
    99004 FORMAT (//5X,'TEST RESULTS:    A     B    ERR     TRUE',16X,'CALC',&
      15X,'ERROR')
    !
    ifail = 0
    !
    skip = .FALSE.
    DO i = 1, npairs
      !               ---------------------------------------------
      calc = DPCHIA(n,x,f,d,1,skip,a(i),b(i),ierr)
      !               ---------------------------------------------
      IF( ierr>=0 ) THEN
        fail = ierr/=ierexp(i)
        true = ANTDER(b(i)) - ANTDER(a(i))
        error = calc - true
        IF( Kprint>=3 ) THEN
          IF( fail ) THEN
            WRITE (Lun,99005) a(i), b(i), ierr, true, calc, error, ierexp(i)
            99005 FORMAT (2F6.1,I5,1P,2D20.10,D15.5,'  (',I1,') *****')
          ELSE
            WRITE (Lun,99010) a(i), b(i), ierr, true, calc, error
          END IF
        END IF
        !
        error = ABS(error)/MAX(one,ABS(true))
        IF( fail .OR. (error>tol) ) ifail = ifail + 1
        IF( i==1 ) THEN
          errmax = error
        ELSE
          errmax = MAX(errmax,error)
        END IF
      ELSE
        IF( Kprint>=3 ) WRITE (Lun,99010) a(i), b(i), ierr
        ifail = ifail + 1
      END IF
    END DO
    !
    !  PRINT SUMMARY.
    !
    IF( Kprint>=2 ) THEN
      WRITE (Lun,99006) errmax, tol
      99006 FORMAT (/'  MAXIMUM RELATIVE ERROR IS:',1P,D15.5,',   TOLERANCE:',1P,&
        D15.5)
      IF( ifail/=0 ) WRITE (Lun,99007) ifail
      99007 FORMAT (/' *** TROUBLE ***',I5,' INTEGRATION TESTS FAILED.')
    END IF
    !
    !  TERMINATE.
    !
    IF( ifail==0 ) THEN
      Ipass = 1
      IF( Kprint>=2 ) WRITE (Lun,99008)
      99008 FORMAT (/' ------------ DPCHIP PASSED  ALL INTEGRATION TESTS',&
        ' ------------')
    ELSE
      Ipass = 0
      IF( Kprint>=1 ) WRITE (Lun,99009)
      99009 FORMAT (/' ************ DPCHIP FAILED SOME INTEGRATION TESTS',&
        ' ************')
    END IF
    !
    RETURN
    99010 FORMAT (2F6.1,I5,1P,2D20.10,D15.5)
    !------------- LAST LINE OF DPCHQ2 FOLLOWS -----------------------------
  CONTAINS

    !
    !  DEFINE TEST FUNCTIONS.
    !
    REAL(DP) FUNCTION FCN(ax)
      REAL(DP), INTENT(IN) :: ax
      FCN = three*ax*ax*(ax-two)
    END FUNCTION FCN
    REAL(DP) FUNCTION DERIV(ax)
      REAL(DP), INTENT(IN) :: ax
      DERIV = three*ax*(two*(ax-two)+ax)
    END FUNCTION DERIV
    REAL(DP) FUNCTION ANTDER(ax)
      REAL(DP), INTENT(IN) :: ax
      ANTDER = ax**3*(thrqtr*ax-two)
    END FUNCTION ANTDER
  END SUBROUTINE DPCHQ2
  !** DPCHQ3
  SUBROUTINE DPCHQ3(Lun,Kprint,Ipass)
    !> Test the PCHIP interpolators DPCHIC, DPCHIM, DPCHSP.
    !***
    ! **Library:**   SLATEC (PCHIP)
    !***
    ! **Type:**      DOUBLE PRECISION (PCHQK3-S, DPCHQ3-D)
    !***
    ! **Keywords:**  PCHIP INTERPOLATOR QUICK CHECK
    !***
    ! **Author:**  Fritsch, F. N., (LLNL)
    !***
    ! **Description:**
    !
    !             DPCHIP QUICK CHECK NUMBER 3
    !
    !     TESTS THE INTERPOLATORS:  DPCHIC, DPCHIM, DPCHSP.
    !- Usage:
    !
    !        INTEGER  LUN, KPRINT, IPASS
    !
    !        CALL DPCHQ3 (LUN, KPRINT, IPASS)
    !
    !- Arguments:
    !
    !     LUN   :IN  is the unit number to which output is to be written.
    !
    !     KPRINT:IN  controls the amount of output, as specified in the
    !                SLATEC Guidelines.
    !
    !     IPASS:OUT  will contain a pass/fail flag.  IPASS=1 is good.
    !                IPASS=0 indicates one or more tests failed.
    !
    !- Description:
    !
    !   This routine interpolates a constructed data set with all three
    !   DPCHIP interpolators and compares the results with those obtained
    !   on a Cray X/MP.  Two different values of the DPCHIC parameter SWITCH
    !   are used.
    !
    !- Remarks:
    !     1. The Cray results are given only to nine significant figures,
    !        so don't expect them to match to more.
    !     2. The results will depend to some extent on the accuracy of
    !        the EXP function.
    !
    !***
    ! **Routines called:**  COMP, D1MACH, DPCHIC, DPCHIM, DPCHSP

    !* REVISION HISTORY  (YYMMDD)
    !   900309  DATE WRITTEN
    !   900314  Converted to a subroutine and added a SLATEC 4.0 prologue.
    !   900315  Revised prologue and improved some output formats.  (FNF)
    !   900316  Made TOLD machine-dependent and added extra output when
    !           KPRINT=3.  (FNF)
    !   900320  Added E0's to DATA statement for X to reduce single/double
    !           differences, and other minor cosmetic changes.
    !   900320  Converted to double precision.
    !   900321  Removed IFAIL from call sequence for SLATEC standards and
    !           made miscellaneous cosmetic changes.  (FNF)
    !   900322  Minor changes to reduce single/double differences.  (FNF)
    !   900530  Tolerance (TOLD) and argument to DPCHIC changed.  (WRB)
    !   900802  Modified TOLD formula and constants in DPCHIC calls to
    !           correct DPCHQ3 failures.  (FNF)
    !   901130  Several significant changes:  (FNF)
    !           1. Changed comparison between DPCHIM and DPCHIC to only
    !              require agreement to machine precision.
    !           2. Revised to print more output when KPRINT=3.
    !           3. Added 1P's to formats.
    !   910708  Minor modifications in use of KPRINT.  (WRB)
    !   930317  Improved output formats.  (FNF)
    USE slatec, ONLY : D1MACH, DPCHIC, DPCHIM, DPCHSP
    !
    !*Internal Notes:
    !
    !     TOLD is used to compare with stored Cray results.  Its value
    !          should be consistent with significance of stored values.
    !     TOLZ is used for cases in which exact equality is expected.
    !     TOL  is used for cases in which agreement to machine precision
    !          is expected.
    !**End
    !
    !  Declare arguments.
    !
    INTEGER :: Lun, Kprint, Ipass
    !
    !  Declare variables.
    !
    INTEGER :: i, ierr, ifail, nbad, nbadz
    INTEGER, PARAMETER :: N = 9, NWK = 2*N
    REAL(DP) :: d(N), dc(N), err, f(N), tol, told, tolz, vc(2), wk(NWK)
    REAL(DP), PARAMETER :: ZERO = 0._DP, MONE = -1._DP
    CHARACTER(6) :: result
    !
    !  Initialize.
    !
    !       Data.
    INTEGER, PARAMETER :: ic(2) = [ 0, 0 ]
    REAL(DP), PARAMETER :: x(N) = [ -2.2_DP, -1.2_DP, -1._DP, -0.5_DP, -0.01_DP, 0.5_DP, &
      1._DP, 2._DP, 2.2_DP ]
    !
    !       Results generated on Cray X/MP (9 sign. figs.)
    REAL(DP), PARAMETER :: dm(N) = [ 0._DP, 3.80027352E-01_DP, 7.17253009E-01_DP, &
      5.82014161E-01_DP, 0._DP, -5.68208031E-01_DP, -5.13501618E-01_DP, -7.77910977E-02_DP, &
      -2.45611117E-03_DP ]
    REAL(DP), PARAMETER :: dc5 = 1.76950158E-02_DP, dc6 = -5.69579814E-01_DP
    REAL(DP), PARAMETER :: ds(N) = [ -5.16830792E-02_DP, 5.71455855E-01_DP, &
      7.40530225E-01_DP, 7.63864934E-01_DP, 1.92614386E-02_DP, -7.65324380E-01_DP, &
      -7.28209035E-01_DP, -7.98445427E-02_DP, -2.85983446E-02_DP ]
    !
    !* FIRST EXECUTABLE STATEMENT  DPCHQ3
    ifail = 0
    !
    !        Set tolerances.
    tol = 10*D1MACH(4)
    told = SQRT(D1MACH(4))
    tolz = ZERO
    !
    IF( Kprint>=3 ) WRITE (Lun,99001)
    !
    ! FORMATS.
    !
    99001 FORMAT ('1'//10X,'TEST DPCHIP INTERPOLATORS')
    IF( Kprint>=2 ) WRITE (Lun,99002)
    99002 FORMAT (//10X,'DPCHQ3 RESULTS'/10X,'--------------')
    !
    !  Set up data.
    !
    DO i = 1, N
      f(i) = EXP(-x(i)**2)
    END DO
    !
    IF( Kprint>=3 ) THEN
      WRITE (Lun,99003)
      99003 FORMAT (//5X,'DATA:'/39X,'---------- EXPECTED D-VALUES ----------'/12X,&
        'X',9X,'F',18X,'DM',13X,'DC',13X,'DS')
      DO i = 1, 4
        WRITE (Lun,99009) x(i), f(i), dm(i), ds(i)
      END DO
      WRITE (Lun,99010) x(5), f(5), dm(5), dc5, ds(5)
      WRITE (Lun,99010) x(6), f(6), dm(6), dc6, ds(6)
      DO i = 7, N
        WRITE (Lun,99009) x(i), f(i), dm(i), ds(i)
      END DO
    END IF
    !
    !  Test DPCHIM.
    !
    IF( Kprint>=3 ) WRITE (Lun,99011) 'IM'
    !     --------------------------------
    CALL DPCHIM(N,x,f,d,1,ierr)
    !     --------------------------------
    !        Expect IERR=1 (one monotonicity switch).
    IF( Kprint>=3 ) WRITE (Lun,99012) 1
    IF( ierr /= 1 ) THEN
      ifail = ifail + 1
    ELSE
      IF( Kprint>=3 ) WRITE (Lun,99013)
      nbad = 0
      nbadz = 0
      DO i = 1, N
        result = '  OK'
        !             D-values should agree with stored values.
        !               (Zero values should agree exactly.)
        IF( dm(i)==ZERO ) THEN
          err = ABS(d(i))
          IF( err>tolz ) THEN
            nbadz = nbadz + 1
            result = '**BADZ'
          END IF
        ELSE
          err = ABS((d(i)-dm(i))/dm(i))
          IF( err>told ) THEN
            nbad = nbad + 1
            result = '**BAD'
          END IF
        END IF
        IF( Kprint>=3 ) WRITE (Lun,99014) i, x(i), d(i), err, result
      END DO
      IF( (nbadz/=0) .OR. (nbad/=0) ) THEN
        ifail = ifail + 1
        IF( (nbadz/=0) .AND. (Kprint>=2) ) WRITE (Lun,99004) nbad
        99004 FORMAT (/'    **',I5,' DPCHIM RESULTS FAILED TO BE EXACTLY ZERO.')
        IF( (nbad/=0) .AND. (Kprint>=2) ) WRITE (Lun,99015) nbad, 'IM', told
      ELSE
        IF( Kprint>=2 ) WRITE (Lun,99016) 'IM'
      END IF
    END IF
    !
    !  Test DPCHIC -- options set to reproduce DPCHIM.
    !
    IF( Kprint>=3 ) WRITE (Lun,99011) 'IC'
    !     --------------------------------------------------------
    CALL DPCHIC(ic,vc,ZERO,N,x,f,dc,1,wk,NWK,ierr)
    !     --------------------------------------------------------
    !        Expect IERR=0 .
    IF( Kprint>=3 ) WRITE (Lun,99012) 0
    IF( ierr /= 0 ) THEN
      ifail = ifail + 1
    ELSE
      IF( Kprint>=3 ) WRITE (Lun,99013)
      nbad = 0
      DO i = 1, N
        result = '  OK'
        !           D-values should agree exactly with those computed by DPCHIM.
        !            (To be generous, will only test to machine precision.)
        err = ABS(d(i)-dc(i))
        IF( err>tol ) THEN
          nbad = nbad + 1
          result = '**BAD'
        END IF
        IF( Kprint>=3 ) WRITE (Lun,99014) i, x(i), dc(i), err, result
      END DO
      IF( nbad/=0 ) THEN
        ifail = ifail + 1
        IF( Kprint>=2 ) WRITE (Lun,99015) nbad, 'IC', tol
      ELSE
        IF( Kprint>=2 ) WRITE (Lun,99016) 'IC'
      END IF
    END IF
    !
    !  Test DPCHIC -- default nonzero switch derivatives.
    !
    IF( Kprint>=3 ) WRITE (Lun,99011) 'IC'
    !     -------------------------------------------------------
    CALL DPCHIC(ic,vc,MONE,N,x,f,d,1,wk,NWK,ierr)
    !     -------------------------------------------------------
    !        Expect IERR=0 .
    IF( Kprint>=3 ) WRITE (Lun,99012) 0
    IF( ierr /= 0 ) THEN
      ifail = ifail + 1
    ELSE
      IF( Kprint>=3 ) WRITE (Lun,99013)
      nbad = 0
      nbadz = 0
      DO i = 1, N
        result = '  OK'
        !            D-values should agree exactly with those computed in
        !            previous call, except at points 5 and 6.
        IF( (i<5) .OR. (i>6) ) THEN
          err = ABS(d(i)-dc(i))
          IF( err>tolz ) THEN
            nbadz = nbadz + 1
            result = '**BADA'
          END IF
        ELSE
          IF( i==5 ) THEN
            err = ABS((d(i)-dc5)/dc5)
          ELSE
            err = ABS((d(i)-dc6)/dc6)
          END IF
          IF( err>told ) THEN
            nbad = nbad + 1
            result = '**BAD'
          END IF
        END IF
        IF( Kprint>=3 ) WRITE (Lun,99014) i, x(i), d(i), err, result
      END DO
      IF( (nbadz/=0) .OR. (nbad/=0) ) THEN
        ifail = ifail + 1
        IF( (nbadz/=0) .AND. (Kprint>=2) ) WRITE (Lun,99005) nbad
        99005 FORMAT (/'    **',I5,' DPCHIC RESULTS FAILED TO AGREE WITH',&
          ' PREVIOUS CALL.')
        IF( (nbad/=0) .AND. (Kprint>=2) ) WRITE (Lun,99015) nbad, 'IC', told
      ELSE
        IF( Kprint>=2 ) WRITE (Lun,99016) 'IC'
      END IF
    END IF
    !
    !  Test DPCHSP.
    !
    IF( Kprint>=3 ) WRITE (Lun,99011) 'SP'
    !     -------------------------------------------------
    CALL DPCHSP(ic,vc,N,x,f,d,1,wk,NWK,ierr)
    !     -------------------------------------------------
    !        Expect IERR=0 .
    IF( Kprint>=3 ) WRITE (Lun,99012) 0
    IF( ierr /= 0 ) THEN
      ifail = ifail + 1
    ELSE
      IF( Kprint>=3 ) WRITE (Lun,99013)
      nbad = 0
      DO i = 1, N
        result = '  OK'
        !             D-values should agree with stored values.
        err = ABS((d(i)-ds(i))/ds(i))
        IF( err>told ) THEN
          nbad = nbad + 1
          result = '**BAD'
        END IF
        IF( Kprint>=3 ) WRITE (Lun,99014) i, x(i), d(i), err, result
      END DO
      IF( nbad/=0 ) THEN
        ifail = ifail + 1
        IF( Kprint>=2 ) WRITE (Lun,99015) nbad, 'SP', told
      ELSE
        IF( Kprint>=2 ) WRITE (Lun,99016) 'SP'
      END IF
    END IF
    !
    !  PRINT SUMMARY AND TERMINATE.
    !
    IF( (Kprint>=2) .AND. (ifail/=0) ) WRITE (Lun,99006) ifail
    99006 FORMAT (/' *** TROUBLE ***',I5,' INTERPOLATION TESTS FAILED.')
    !
    IF( ifail==0 ) THEN
      Ipass = 1
      IF( Kprint>=2 ) WRITE (Lun,99007)
      99007 FORMAT (/' ------------ DPCHIP PASSED  ALL INTERPOLATION TESTS',&
        ' ------------')
    ELSE
      Ipass = 0
      IF( Kprint>=1 ) WRITE (Lun,99008)
      99008 FORMAT (/' ************ DPCHIP FAILED SOME INTERPOLATION TESTS',&
        ' ************')
    END IF
    !
    RETURN
    99009 FORMAT (5X,F10.2,1P,D15.5,4X,D15.5,15X,D15.5)
    99010 FORMAT (5X,F10.2,1P,D15.5,4X,3D15.5)
    99011 FORMAT (/5X,'DPCH',A2,' TEST:')
    99012 FORMAT (15X,'EXPECT  IERR =',I5)
    99013 FORMAT (/9X,'I',7X,'X',9X,'D',13X,'ERR')
    99014 FORMAT (5X,I5,F10.2,1P,2D15.5,2X,A)
    99015 FORMAT (/'    **',I5,' DPCH',A2,' RESULTS FAILED TOLERANCE TEST.',&
      '  TOL =',1P,D10.3)
    99016 FORMAT (/5X,'  ALL DPCH',A2,' RESULTS OK.')
    !------------- LAST LINE OF DPCHQ3 FOLLOWS -----------------------------
  END SUBROUTINE DPCHQ3
  !** DPCHQ4
  SUBROUTINE DPCHQ4(Lun,Kprint,Ipass)
    !> Test the PCHIP monotonicity checker DPCHCM.
    !***
    ! **Library:**   SLATEC (PCHIP)
    !***
    ! **Type:**      DOUBLE PRECISION (PCHQK4-S, DPCHQ4-D)
    !***
    ! **Keywords:**  PCHIP MONOTONICITY CHECKER QUICK CHECK
    !***
    ! **Author:**  Fritsch, F. N., (LLNL)
    !***
    ! **Description:**
    !
    !             DPCHIP QUICK CHECK NUMBER 4
    !
    !     TESTS THE MONOTONICITY CHECKER:  DPCHCM.
    !- Usage:
    !
    !        INTEGER  LUN, KPRINT, IPASS
    !
    !        CALL DPCHQ4 (LUN, KPRINT, IPASS)
    !
    !- Arguments:
    !
    !     LUN   :IN  is the unit number to which output is to be written.
    !
    !     KPRINT:IN  controls the amount of output, as specified in the
    !                SLATEC Guidelines.
    !
    !     IPASS:OUT  will contain a pass/fail flag.  IPASS=1 is good.
    !                IPASS=0 indicates one or more tests failed.
    !
    !- Description:
    !
    !   This routine tests a constructed data set with three different
    !   INCFD settings and compares with the expected results.  It then
    !   runs a special test to check for bug in overall monotonicity found
    !   in DPCHMC.  Finally, it reverses the data and repeats all tests.
    !
    !***
    ! **Routines called:**  DPCHCM

    !* REVISION HISTORY  (YYMMDD)
    !   890208  DATE WRITTEN
    !   890306  Changed LOUT to LUN and added it to call list.  (FNF)
    !   890316  Removed DATA statements to suit new quick check standards.
    !   890410  Changed PCHMC to PCHCM.
    !   890410  Added a SLATEC 4.0 format prologue.
    !   900314  Changed name from PCHQK3 to PCHQK4 and improved some output
    !           formats.
    !   900315  Revised prologue and improved some output formats.  (FNF)
    !   900320  Converted to double precision.
    !   900321  Removed IFAIL from call sequence for SLATEC standards and
    !           made miscellaneous cosmetic changes.  (FNF)
    !   900322  Added declarations so all variables are declared.  (FNF)
    !   910708  Minor modifications in use of KPRINT.  (WRB)
    !   930317  Improved output formats.  (FNF)
    USE slatec, ONLY : DPCHCM
    !
    !*Internal Notes:
    !
    !     Data set-up is done via assignment statements to avoid modifying
    !     DATA-loaded arrays, as required by the 1989 SLATEC Guidelines.
    !     Run with KPRINT=3 to display the data.
    !**End
    !
    !  Declare arguments.
    !
    INTEGER :: Lun, Kprint, Ipass
    !
    !  DECLARE VARIABLES.
    !
    INTEGER, PARAMETER :: MAXN = 16, MAXN2 = 8, MAXN3 = 6, NB = 7
    INTEGER :: i, ierr, ifail, incfd, ismon(MAXN), k, n
    REAL(DP) :: d(MAXN), db(NB), f(MAXN), fb(NB), x(MAXN)
    LOGICAL :: skip
    !
    !  DEFINE EXPECTED RESULTS.
    !
    INTEGER :: ismex1(MAXN) = [ 1, 1, -1, 1, 1, -1, 1, 1, -1, 1, 1, &
      -1, 1, 1, -1, 2 ]
    INTEGER :: ismex2(MAXN2) = [ 1, 2, 2, 1, 2, 2, 1, 2 ]
    INTEGER :: ismex3(MAXN3) = [ 1, 1, 1, 1, 1, 1 ]
    INTEGER :: ismexb(NB) = [ 1, 3, 1, -1, -3, -1, 2 ]
    !
    !  DEFINE TEST DATA.
    !
    INTEGER, PARAMETER :: ns(3) = [ 16, 8, 6 ]
    !
    !* FIRST EXECUTABLE STATEMENT  DPCHQ4
    IF( Kprint>=3 ) WRITE (Lun,99001)
    !
    ! FORMATS.
    !
    99001 FORMAT ('1'//10X,'TEST DPCHIP MONOTONICITY CHECKER')
    IF( Kprint>=2 ) WRITE (Lun,99002)
    99002 FORMAT (//10X,'DPCHQ4 RESULTS'/10X,'--------------')
    !
    !       Define X, F, D.
    DO i = 1, MAXN
      x(i) = i
      d(i) = 0._DP
    END DO
    DO i = 2, MAXN, 3
      d(i) = 2._DP
    END DO
    DO i = 1, 3
      f(i) = x(i)
      f(i+3) = f(i) + 1._DP
      f(i+6) = f(i+3) + 1._DP
      f(i+9) = f(i+6) + 1._DP
      f(i+12) = f(i+9) + 1._DP
    END DO
    f(16) = 6._DP
    !       Define FB, DB.
    fb(1) = 0._DP
    fb(2) = 2._DP
    fb(3) = 3._DP
    fb(4) = 5._DP
    db(1) = 1._DP
    db(2) = 3._DP
    db(3) = 3._DP
    db(4) = 0._DP
    DO i = 1, 3
      fb(NB-i+1) = fb(i)
      db(NB-i+1) = -db(i)
    END DO
    !
    !  INITIALIZE.
    !
    ifail = 0
    !
    IF( Kprint>=3 ) THEN
      WRITE (Lun,99003)
      99003 FORMAT (//5X,'DATA:'//9X,'I',4X,'X',5X,'F',5X,'D',5X,'FB',4X,'DB')
      DO i = 1, NB
        WRITE (Lun,99010) i, x(i), f(i), d(i), fb(i), db(i)
      END DO
      DO i = NB + 1, MAXN
        WRITE (Lun,99010) i, x(i), f(i), d(i)
      END DO
    END IF
    !
    !  TRANSFER POINT FOR SECOND SET OF TESTS.
    !
    !
    !  Loop over a series of values of INCFD.
    !
    100 CONTINUE
    DO incfd = 1, 3
      n = ns(incfd)
      skip = .FALSE.
      !        -------------------------------------------------
      CALL DPCHCM(n,x,f,d,incfd,skip,ismon,ierr)
      !        -------------------------------------------------
      IF( Kprint>=3 ) WRITE (Lun,99004) incfd, ierr, (ismon(i),i=1,n)
      99004 FORMAT (/4X,'INCFD =',I2,':  IERR =',I3/15X,'ISMON =',16I3)
      IF( ierr/=0 ) THEN
        ifail = ifail + 1
        IF( Kprint>=3 ) WRITE (Lun,99011)
      ELSE
        DO i = 1, n
          IF( incfd==1 ) THEN
            IF( ismon(i)/=ismex1(i) ) THEN
              ifail = ifail + 1
              IF( Kprint>=3 ) WRITE (Lun,99012) (ismex1(k),k=1,n)
              EXIT
            END IF
          ELSEIF( incfd==2 ) THEN
            IF( ismon(i)/=ismex2(i) ) THEN
              ifail = ifail + 1
              IF( Kprint>=3 ) WRITE (Lun,99012) (ismex2(k),k=1,n)
              EXIT
            END IF
          ELSEIF( ismon(i)/=ismex3(i) ) THEN
            ifail = ifail + 1
            IF( Kprint>=3 ) WRITE (Lun,99012) (ismex3(k),k=1,n)
            EXIT
          END IF
        END DO
      END IF
    END DO
    !
    !  Test for -1,3,1 bug.
    !
    skip = .FALSE.
    !     ------------------------------------------------
    CALL DPCHCM(NB,x,fb,db,1,skip,ismon,ierr)
    !     ------------------------------------------------
    IF( Kprint>=3 ) WRITE (Lun,99005) ierr, (ismon(i),i=1,NB)
    99005 FORMAT (/4X,' Bug test:  IERR =',I3/15X,'ISMON =',7I3)
    IF( ierr/=0 ) THEN
      ifail = ifail + 1
      IF( Kprint>=3 ) WRITE (Lun,99011)
    ELSE
      DO i = 1, NB
        IF( ismon(i)/=ismexb(i) ) THEN
          ifail = ifail + 1
          IF( Kprint>=3 ) WRITE (Lun,99012) (ismexb(k),k=1,NB)
          EXIT
        END IF
      END DO
    END IF
    !
    IF( f(1)<0._SP ) THEN
      !
      !  PRINT SUMMARY AND TERMINATE.
      !
      IF( (Kprint>=2) .AND. (ifail/=0) ) WRITE (Lun,99006) ifail
      99006 FORMAT (/' *** TROUBLE ***',I5,' MONOTONICITY TESTS FAILED.')
      !
      IF( ifail==0 ) THEN
        Ipass = 1
        IF( Kprint>=2 ) WRITE (Lun,99007)
        99007 FORMAT (/' ------------ DPCHIP PASSED  ALL MONOTONICITY TESTS',&
          ' ------------')
      ELSE
        Ipass = 0
        IF( Kprint>=1 ) WRITE (Lun,99008)
        99008 FORMAT (/' ************ DPCHIP FAILED SOME MONOTONICITY TESTS',&
          ' ************')
      END IF
      !
      RETURN
    ELSE
      !
      !  Change sign and do again.
      !
      IF( Kprint>=3 ) WRITE (Lun,99009)
      99009 FORMAT (/4X,'Changing sign of data.....')
      DO i = 1, MAXN
        f(i) = -f(i)
        d(i) = -d(i)
        IF( ismex1(i)/=2 ) ismex1(i) = -ismex1(i)
      END DO
      DO i = 1, MAXN2
        IF( ismex2(i)/=2 ) ismex2(i) = -ismex2(i)
      END DO
      DO i = 1, MAXN3
        IF( ismex3(i)/=2 ) ismex3(i) = -ismex3(i)
      END DO
      DO i = 1, NB
        fb(i) = -fb(i)
        db(i) = -db(i)
        IF( ismexb(i)/=2 ) ismexb(i) = -ismexb(i)
      END DO
      GOTO 100
    END IF
    99010 FORMAT (5X,I5,5F6.1)
    99011 FORMAT (' *** Failed -- bad IERR value.')
    99012 FORMAT (' *** Failed -- expect:',16I3)
    !------------- LAST LINE OF DPCHQ4 FOLLOWS -----------------------------
  END SUBROUTINE DPCHQ4
  !** DPCHQ5
  SUBROUTINE DPCHQ5(Lun,Kprint,Ipass)
    !> Test the PCH to B-spline conversion routine DPCHBS.
    !***
    ! **Library:**   SLATEC (PCHIP)
    !***
    ! **Type:**      DOUBLE PRECISION (PCHQK5-S, DPCHQ5-D)
    !***
    ! **Keywords:**  PCHIP CONVERSION ROUTINE QUICK CHECK
    !***
    ! **Author:**  Fritsch, F. N., (LLNL)
    !***
    ! **Description:**
    !
    !             DPCHIP QUICK CHECK NUMBER 5
    !
    !     TESTS THE CONVERSION ROUTINE:  DPCHBS.
    !- Usage:
    !
    !        INTEGER  LUN, KPRINT, IPASS
    !
    !        CALL DPCHQ5 (LUN, KPRINT, IPASS)
    !
    !- Arguments:
    !
    !     LUN   :IN  is the unit number to which output is to be written.
    !
    !     KPRINT:IN  controls the amount of output, as specified in the
    !                SLATEC Guidelines.
    !
    !     IPASS:OUT  will contain a pass/fail flag.  IPASS=1 is good.
    !                IPASS=0 indicates one or more tests failed.
    !
    !- Description:
    !
    !   This routine tests a constructed data set with four different
    !   KNOTYP settings.  It computes the function and derivatives of the
    !   resulting B-representation via DBVALU and compares with PCH data.
    !
    !- Caution:
    !   This routine assumes DBVALU has already been successfully tested.
    !
    !***
    ! **Routines called:**  DBVALU, DPCHBS, D1MACH

    !* REVISION HISTORY  (YYMMDD)
    !   900411  DATE WRITTEN
    !   900412  Corrected minor errors in initial implementation.
    !   900430  Produced double precision version.
    !   900501  Corrected declarations.
    !   930317  Improved output formats.  (FNF)
    USE slatec, ONLY : DBVALU, DPCHBS, D1MACH
    !
    !*Internal Notes:
    !  TOL  is the tolerance to use for quantities that should only
    !       theoretically be equal.
    !  TOLZ is the tolerance to use for quantities that should be exactly
    !       equal.
    !
    !**End
    !
    !  Declare arguments.
    !
    INTEGER :: Lun, Kprint, Ipass
    !
    !  Declare variables.
    !
    INTEGER :: i, ierr, ifail, inbv, j, knotyp, k, ndim, nknots
    INTEGER, PARAMETER :: N = 9
    REAL(DP) :: bcoef(2*N), dcalc, derr, dermax, fcalc, ferr, fermax, t(2*N+4), &
      terr, termax, tol, tolz, tsave(2*N+4), work(16*N)
    REAL(DP), PARAMETER :: ZERO = 0._DP
    LOGICAL :: fail
    !
    !  Define test data.
    !
    REAL(DP), PARAMETER :: x(N) = [ -2.2_DP, -1.2_DP, -1._DP, -0.5_DP, -0.01_DP, 0.5_DP, &
      1._DP, 2._DP, 2.2_DP ]
    REAL(DP), PARAMETER :: f(N) = [ 0.0079_DP, 0.2369_DP, 0.3679_DP, 0.7788_DP, 0.9999_DP, &
      0.7788_DP, 0.3679_DP, 0.1083_DP, 0.0079_DP ]
    REAL(DP), PARAMETER :: d(N) = [ 0.0000_DP, 0.3800_DP, 0.7173_DP, 0.5820_DP, 0.0177_DP, &
      -0.5696_DP, -0.5135_DP, -0.0778_DP, -0.0025_DP ]
    !
    !  Initialize.
    !
    !* FIRST EXECUTABLE STATEMENT  DPCHQ5
    ifail = 0
    tol = 100*D1MACH(4)
    tolz = ZERO
    !
    IF( Kprint>=3 ) WRITE (Lun,99001)
    !
    ! FORMATS.
    !
    99001 FORMAT ('1'//10X,'TEST PCH TO B-SPLINE CONVERTER')
    IF( Kprint>=2 ) WRITE (Lun,99002)
    99002 FORMAT (//10X,'DPCHQ5 RESULTS'/10X,'--------------')
    !
    !  Loop over a series of values of KNOTYP.
    !
    IF( Kprint>=3 ) WRITE (Lun,99003)
    99003 FORMAT (/4X,'(Results should be the same for all KNOTYP values.)')
    DO knotyp = 2, -1, -1
      !        ------------
      CALL DPCHBS(N,x,f,d,1,knotyp,nknots,t,bcoef,ndim,k,ierr)
      !        ------------
      IF( Kprint>=3 ) WRITE (Lun,99004) knotyp, nknots, ndim, k, ierr
      99004 FORMAT (/4X,'KNOTYP =',I2,':  NKNOTS =',I3,',  NDIM =',I3,',  K =',I2,&
        ',  IERR =',I3)
      IF( ierr/=0 ) THEN
        ifail = ifail + 1
        IF( Kprint>=3 ) WRITE (Lun,99005)
        99005 FORMAT (' *** Failed -- bad IERR value.')
      ELSE
        !             Compare evaluated results with inputs to DPCHBS.
        inbv = 1
        fermax = ZERO
        dermax = ZERO
        IF( Kprint>=3 ) THEN
          WRITE (Lun,99006)
          99006 FORMAT (/15X,'X',9X,'KNOTS',10X,'F',7X,'FERR',8X,'D',7X,'DERR')
          WRITE (Lun,99013) t(1), t(2)
          j = 1
        END IF
        DO i = 1, N
          fcalc = DBVALU(t,bcoef,ndim,k,0,x(i),inbv,work)
          ferr = f(i) - fcalc
          fermax = MAX(fermax,RELERR(ferr,f(i)))
          dcalc = DBVALU(t,bcoef,ndim,k,1,x(i),inbv,work)
          derr = d(i) - dcalc
          dermax = MAX(dermax,RELERR(derr,d(i)))
          IF( Kprint>=3 ) THEN
            j = j + 2
            WRITE (Lun,99007) x(i), t(j), t(j+1), f(i), ferr, d(i), derr
            99007 FORMAT (10X,3F8.2,F10.4,1P,D10.2,0P,F10.4,1P,D10.2)
          END IF
        END DO
        IF( Kprint>=3 ) THEN
          j = j + 2
          WRITE (Lun,99013) t(j), t(j+1)
        END IF
        fail = (fermax>tol) .OR. (dermax>tol)
        IF( fail ) ifail = ifail + 1
        IF( (Kprint>=3) .OR. (Kprint>=2) .AND. fail ) WRITE (Lun,99008) fermax, &
          dermax, tol
        99008 FORMAT (/5X,'Maximum relative errors:'/15X,'F-error =',1P,D13.5,5X,&
          'D-error =',D13.5/5X,'Both should be less than  TOL =',D13.5)
      END IF
      !
      !          Special check for KNOTYP=-1.
      IF( knotyp==0 ) THEN
        !             Save knot vector for next test.
        DO i = 1, nknots
          tsave(i) = t(i)
        END DO
      ELSEIF( knotyp==-1 ) THEN
        !             Check that knot vector is unchanged.
        termax = ZERO
        DO i = 1, nknots
          terr = ABS(t(i)-tsave(i))
          termax = MAX(termax,terr)
        END DO
        IF( termax>tolz ) THEN
          ifail = ifail + 1
          IF( Kprint>=2 ) WRITE (Lun,99009) termax, tolz
          99009 FORMAT (/' *** T-ARRAY MAXIMUM CHANGE =',1P,D13.5,&
            ';  SHOULD NOT EXCEED TOLZ =',D13.5)
        END IF
      END IF
    END DO
    !
    !  PRINT SUMMARY AND TERMINATE.
    !
    IF( (Kprint>=2) .AND. (ifail/=0) ) WRITE (Lun,99010) ifail
    99010 FORMAT (/' *** TROUBLE ***',I5,' CONVERSION TESTS FAILED.')
    !
    IF( ifail==0 ) THEN
      Ipass = 1
      IF( Kprint>=2 ) WRITE (Lun,99011)
      99011 FORMAT (/' ------------ DPCHIP PASSED  ALL CONVERSION TESTS',&
        ' ------------')
    ELSE
      Ipass = 0
      IF( Kprint>=1 ) WRITE (Lun,99012)
      99012 FORMAT (/' ************ DPCHIP FAILED SOME CONVERSION TESTS',&
        ' ************')
    END IF
    !
    RETURN
    99013 FORMAT (18X,2F8.2)
    !------------- LAST LINE OF DPCHQ5 FOLLOWS -----------------------------
  CONTAINS
    !
    !  Define relative error function.
    !
    REAL(DP) FUNCTION RELERR(err,ans)
      REAL(DP), INTENT(IN) :: ans, err
      RELERR = ABS(err)/MAX(1.E-5_DP,ABS(ans))
    END FUNCTION RELERR
  END SUBROUTINE DPCHQ5
END MODULE TEST33_MOD
!** TEST33
PROGRAM TEST33
  USE TEST33_MOD, ONLY : DPCHQ1, DPCHQ2, DPCHQ3, DPCHQ4, DPCHQ5
  USE slatec, ONLY : I1MACH, control_xer, max_xer
  USE common_mod, ONLY : GET_ARGUMENT
  IMPLICIT NONE
  !> Driver for testing SLATEC subprograms
  !***
  ! **Library:**   SLATEC
  !***
  ! **Category:**  E1A
  !***
  ! **Type:**      DOUBLE PRECISION (TEST32-S, TEST33-D)
  !***
  ! **Keywords:**  QUICK CHECK DRIVER
  !***
  ! **Author:**  SLATEC Common Mathematical Library Committee
  !***
  ! **Description:**
  !
  !- Usage:
  !     One input data record is required
  !         READ (LIN, '(I1)') KPRINT
  !
  !- Arguments:
  !     KPRINT = 0  Quick checks - No printing.
  !                 Driver       - Short pass or fail message printed.
  !              1  Quick checks - No message printed for passed tests,
  !                                short message printed for failed tests.
  !                 Driver       - Short pass or fail message printed.
  !              2  Quick checks - Print short message for passed tests,
  !                                fuller information for failed tests.
  !                 Driver       - Pass or fail message printed.
  !              3  Quick checks - Print complete quick check results.
  !                 Driver       - Pass or fail message printed.
  !
  !- Description:
  !     Driver for testing SLATEC subprograms
  !        DPCHIP
  !
  !***
  ! **References:**  Kirby W. Fong, Thomas H. Jefferson, Tokihiko Suyehiro
  !                 and Lee Walton, Guide to the SLATEC Common Mathema-
  !                 tical Library, April 10, 1990.
  !***
  ! **Routines called:**  DPCHQ1, DPCHQ2, DPCHQ3, DPCHQ4, DPCHQ5, I1MACH,
  !                    XERMAX, XSETF

  !* REVISION HISTORY  (YYMMDD)
  !   890618  DATE WRITTEN
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900319  Corrected category record.  (FNF)
  !   900320  Added new quick checks DPCHQ3, DPCHQ4.  (FNF)
  !   900321  Moved IPASS to call sequences for SLATEC standards.  (FNF)
  !   900322  Corrected list of routines called.  (FNF)
  !   900524  Cosmetic changes to code.  (WRB)
  !   930318  Added new quick check DPCHQ5.  (WRB,FNF)
  INTEGER :: ipass, kprint, lin, lun, nfail
  !* FIRST EXECUTABLE STATEMENT  TEST33
  lun = I1MACH(2)
  lin = I1MACH(1)
  nfail = 0
  !
  !     Read KPRINT parameter
  !
  CALL GET_ARGUMENT(kprint)
  max_xer = 1000
  IF( kprint<=1 ) THEN
    control_xer = 0
  ELSE
    control_xer = 1
  END IF
  !
  !     Test DPCHIP evaluators
  !
  CALL DPCHQ1(lun,kprint,ipass)
  IF( ipass==0 ) nfail = nfail + 1
  !
  !     Test DPCHIP integrators
  !
  CALL DPCHQ2(lun,kprint,ipass)
  IF( ipass==0 ) nfail = nfail + 1
  !
  !     Test DPCHIP interpolators
  !
  CALL DPCHQ3(lun,kprint,ipass)
  IF( ipass==0 ) nfail = nfail + 1
  !
  !     Test DPCHIP monotonicity checker
  !
  CALL DPCHQ4(lun,kprint,ipass)
  IF( ipass==0 ) nfail = nfail + 1
  !
  !     Test PCH to B-spline conversion.
  !
  CALL DPCHQ5(lun,kprint,ipass)
  IF( ipass==0 ) nfail = nfail + 1
  !
  !     Write PASS or FAIL message
  !
  IF( nfail==0 ) THEN
    WRITE (lun,99001)
    99001 FORMAT (/' --------------TEST33 PASSED ALL TESTS----------------')
  ELSE
    WRITE (lun,99002) nfail
    99002 FORMAT (/' ************* WARNING -- ',I5,&
      ' TEST(S) FAILED IN PROGRAM TEST33 *************')
  END IF
  STOP
END PROGRAM TEST33
