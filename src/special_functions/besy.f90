!** BESY
SUBROUTINE BESY(X,Fnu,N,Y)
  IMPLICIT NONE
  !>
  !***
  !  Implement forward recursion on the three term recursion
  !            relation for a sequence of non-negative order Bessel
  !            functions Y/SUB(FNU+I-1)/(X), I=1,...,N for real, positive
  !            X and non-negative orders FNU.
  !***
  ! **Library:**   SLATEC
  !***
  ! **Category:**  C10A3
  !***
  ! **Type:**      SINGLE PRECISION (BESY-S, DBESY-D)
  !***
  ! **Keywords:**  SPECIAL FUNCTIONS, Y BESSEL FUNCTION
  !***
  ! **Author:**  Amos, D. E., (SNLA)
  !***
  ! **Description:**
  !
  !     Abstract
  !         BESY implements forward recursion on the three term
  !         recursion relation for a sequence of non-negative order Bessel
  !         functions Y/sub(FNU+I-1)/(X), I=1,N for real X .GT. 0.0E0 and
  !         non-negative orders FNU.  If FNU .LT. NULIM, orders FNU and
  !         FNU+1 are obtained from BESYNU which computes by a power
  !         series for X .LE. 2, the K Bessel function of an imaginary
  !         argument for 2 .LT. X .LE. 20 and the asymptotic expansion for
  !         X .GT. 20.
  !
  !         If FNU .GE. NULIM, the uniform asymptotic expansion is coded
  !         in ASYJY for orders FNU and FNU+1 to start the recursion.
  !         NULIM is 70 or 100 depending on whether N=1 or N .GE. 2.  An
  !         overflow test is made on the leading term of the asymptotic
  !         expansion before any extensive computation is done.
  !
  !     Description of Arguments
  !
  !         Input
  !           X      - X .GT. 0.0E0
  !           FNU    - order of the initial Y function, FNU .GE. 0.0E0
  !           N      - number of members in the sequence, N .GE. 1
  !
  !         Output
  !           Y      - a vector whose first N components contain values
  !                    for the sequence Y(I)=Y/sub(FNU+I-1)/(X), I=1,N.
  !
  !     Error Conditions
  !         Improper input arguments - a fatal error
  !         Overflow - a fatal error
  !
  !***
  ! **References:**  F. W. J. Olver, Tables of Bessel Functions of Moderate
  !                 or Large Orders, NPL Mathematical Tables 6, Her
  !                 Majesty's Stationery Office, London, 1962.
  !               N. M. Temme, On the numerical evaluation of the modified
  !                 Bessel function of the third kind, Journal of
  !                 Computational Physics 19, (1975), pp. 324-337.
  !               N. M. Temme, On the numerical evaluation of the ordinary
  !                 Bessel function of the second kind, Journal of
  !                 Computational Physics 21, (1976), pp. 343-350.
  !***
  ! **Routines called:**  ASYJY, BESY0, BESY1, BESYNU, I1MACH, R1MACH,
  !                    XERMSG, YAIRY

  !* REVISION HISTORY  (YYMMDD)
  !   800501  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900326  Removed duplicate information from DESCRIPTION section.
  !           (WRB)
  !   920501  Reformatted the REFERENCES section.  (WRB)

  !
  INTEGER i, iflw, j, N, nb, nd, nn, nud
  REAL azn, cn, dnu, elim, flgjy, fn, Fnu, ran, s, s1, s2, tm, &
    trx, w(2), wk(7), w2n, X, xlim, xxn, Y(*)
  INTEGER, EXTERNAL :: I1MACH
  REAL, EXTERNAL :: BESY0, BESY1, R1MACH
  EXTERNAL :: YAIRY
  INTEGER, PARAMETER :: nulim(2) = [ 70, 100 ]
  !* FIRST EXECUTABLE STATEMENT  BESY
  nn = -I1MACH(12)
  elim = 2.303E0*(nn*R1MACH(5)-3.0E0)
  xlim = R1MACH(1)*1.0E+3
  IF ( Fnu<0.0E0 ) THEN
    !
    !
    !
    CALL XERMSG('SLATEC','BESY','ORDER, FNU, LESS THAN ZERO',2,1)
    RETURN
  ELSEIF ( X<=0.0E0 ) THEN
    CALL XERMSG('SLATEC','BESY','X LESS THAN OR EQUAL TO ZERO',2,1)
    RETURN
  ELSEIF ( X<xlim ) THEN
    CALL XERMSG('SLATEC','BESY',&
      'OVERFLOW, FNU OR N TOO LARGE OR X TOO SMALL',6,1)
    RETURN
  ELSEIF ( N<1 ) THEN
    CALL XERMSG('SLATEC','BESY','N LESS THAN ONE',2,1)
    RETURN
  ELSE
    !
    !     ND IS A DUMMY VARIABLE FOR N
    !
    nd = N
    nud = INT(Fnu)
    dnu = Fnu - nud
    nn = MIN(2,nd)
    fn = Fnu + N - 1
    IF ( fn<2.0E0 ) THEN
      !
      !     OVERFLOW TEST
      IF ( fn<=1.0E0 ) GOTO 200
      IF ( -fn*(LOG(X)-0.693E0)<=elim ) GOTO 200
      CALL XERMSG('SLATEC','BESY',&
        'OVERFLOW, FNU OR N TOO LARGE OR X TOO SMALL',6,1)
      RETURN
    ELSE
      !
      !     OVERFLOW TEST  (LEADING EXPONENTIAL OF ASYMPTOTIC EXPANSION)
      !     FOR THE LAST ORDER, FNU+N-1.GE.NULIM
      !
      xxn = X/fn
      w2n = 1.0E0 - xxn*xxn
      IF ( w2n>0.0E0 ) THEN
        ran = SQRT(w2n)
        azn = LOG((1.0E0+ran)/xxn) - ran
        cn = fn*azn
        IF ( cn>elim ) THEN
          CALL XERMSG('SLATEC','BESY',&
            'OVERFLOW, FNU OR N TOO LARGE OR X TOO SMALL',6,1)
          RETURN
        END IF
      END IF
      IF ( nud<nulim(nn) ) THEN
        !
        IF ( dnu/=0.0E0 ) THEN
          nb = 2
          IF ( nud==0.AND.nd==1 ) nb = 1
          CALL BESYNU(X,dnu,nb,w)
          s1 = w(1)
          IF ( nb==1 ) GOTO 20
          s2 = w(2)
        ELSE
          s1 = BESY0(X)
          IF ( nud==0.AND.nd==1 ) GOTO 20
          s2 = BESY1(X)
        END IF
        trx = 2.0E0/X
        tm = (dnu+dnu+2.0E0)/X
        !     FORWARD RECUR FROM DNU TO FNU+1 TO GET Y(1) AND Y(2)
        IF ( nd==1 ) nud = nud - 1
        IF ( nud>0 ) THEN
          DO i = 1, nud
            s = s2
            s2 = tm*s2 - s1
            s1 = s
            tm = tm + trx
          END DO
          IF ( nd==1 ) s1 = s2
        ELSEIF ( nd<=1 ) THEN
          s1 = s2
        END IF
      ELSE
        !
        !     ASYMPTOTIC EXPANSION FOR ORDERS FNU AND FNU+1.GE.NULIM
        !
        flgjy = -1.0E0
        CALL ASYJY(YAIRY,X,Fnu,flgjy,nn,Y,wk,iflw)
        IF ( iflw/=0 ) THEN
          CALL XERMSG('SLATEC','BESY',&
            'OVERFLOW, FNU OR N TOO LARGE OR X TOO SMALL',6,1)
          RETURN
        ELSE
          IF ( nn==1 ) RETURN
          trx = 2.0E0/X
          tm = (Fnu+Fnu+2.0E0)/X
          GOTO 100
        END IF
      END IF
      20  Y(1) = s1
      IF ( nd==1 ) RETURN
      Y(2) = s2
    END IF
  END IF
  100 CONTINUE
  IF ( nd==2 ) RETURN
  !     FORWARD RECUR FROM FNU+2 TO FNU+N-1
  DO i = 3, nd
    Y(i) = tm*Y(i-1) - Y(i-2)
    tm = tm + trx
  END DO
  RETURN
  200 CONTINUE
  IF ( dnu==0.0E0 ) THEN
    j = nud
    IF ( j/=1 ) THEN
      j = j + 1
      Y(j) = BESY0(X)
      IF ( nd==1 ) RETURN
      j = j + 1
    END IF
    Y(j) = BESY1(X)
    IF ( nd==1 ) RETURN
    trx = 2.0E0/X
    tm = trx
    GOTO 100
  ELSE
    CALL BESYNU(X,Fnu,nd,Y)
    RETURN
  END IF
  RETURN
END SUBROUTINE BESY
