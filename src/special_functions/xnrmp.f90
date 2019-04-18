!** XNRMP
SUBROUTINE XNRMP(Nu,Mu1,Mu2,Sarg,Mode,Spn,Ipn,Isig,Ierror)
  !>
  !***
  !  Compute normalized Legendre polynomials.
  !***
  ! **Library:**   SLATEC
  !***
  ! **Category:**  C3A2, C9
  !***
  ! **Type:**      SINGLE PRECISION (XNRMP-S, DXNRMP-D)
  !***
  ! **Keywords:**  LEGENDRE FUNCTIONS
  !***
  ! **Author:**  Lozier, Daniel W., (National Bureau of Standards)
  !           Smith, John M., (NBS and George Mason University)
  !***
  ! **Description:**
  !
  !        SUBROUTINE TO CALCULATE NORMALIZED LEGENDRE POLYNOMIALS
  !        (DXNRMP is double-precision version)
  !        XNRMP calculates normalized Legendre polynomials of varying
  !        order and fixed argument and degree. The order MU and degree
  !        NU are non-negative integers and the argument is real. Because
  !        the algorithm requires the use of numbers outside the normal
  !        machine range, this subroutine employs a special arithmetic
  !        called extended-range arithmetic. See J.M. Smith, F.W.J. Olver,
  !        and D.W. Lozier, Extended-Range Arithmetic and Normalized
  !        Legendre Polynomials, ACM Transactions on Mathematical Soft-
  !        ware, 93-105, March 1981, for a complete description of the
  !        algorithm and special arithmetic. Also see program comments
  !        in XSET.
  !
  !        The normalized Legendre polynomials are multiples of the
  !        associated Legendre polynomials of the first kind where the
  !        normalizing coefficients are chosen so as to make the integral
  !        from -1 to 1 of the square of each function equal to 1. See
  !        E. Jahnke, F. Emde and F. Losch, Tables of Higher Functions,
  !        McGraw-Hill, New York, 1960, p. 121.
  !
  !        The input values to XNRMP are NU, MU1, MU2, SARG, and MODE.
  !        These must satisfy
  !          1. NU .GE. 0 specifies the degree of the normalized Legendre
  !             polynomial that is wanted.
  !          2. MU1 .GE. 0 specifies the lowest-order normalized Legendre
  !             polynomial that is wanted.
  !          3. MU2 .GE. MU1 specifies the highest-order normalized Leg-
  !             endre polynomial that is wanted.
  !         4a. MODE = 1 and -1.0 .LE. SARG .LE. 1.0 specifies that
  !             Normalized Legendre(NU, MU, SARG) is wanted for MU = MU1,
  !             MU1 + 1, ..., MU2.
  !         4b. MODE = 2 and -3.14159... .LT. SARG .LT. 3.14159... spec-
  !             ifies that Normalized Legendre(NU, MU, COS(SARG)) is want-
  !             ed for MU = MU1, MU1 + 1, ..., MU2.
  !
  !        The output of XNRMP consists of the two vectors SPN and IPN
  !        and the error estimate ISIG. The computed values are stored as
  !        extended-range numbers such that
  !             (SPN(1),IPN(1))=NORMALIZED LEGENDRE(NU,MU1,X)
  !             (SPN(2),IPN(2))=NORMALIZED LEGENDRE(NU,MU1+1,X)
  !                .
  !                .
  !             (SPN(K),IPN(K))=NORMALIZED LEGENDRE(NU,MU2,X)
  !        where K = MU2 - MU1 + 1 and X = SARG or COS(SARG) according
  !        to whether MODE = 1 or 2. Finally, ISIG is an estimate of the
  !        number of decimal digits lost through rounding errors in the
  !        computation. For example if SARG is accurate to 12 significant
  !        decimals, then the computed function values are accurate to
  !        12 - ISIG significant decimals (except in neighborhoods of
  !        zeros).
  !
  !        The interpretation of (SPN(I),IPN(I)) is SPN(I)*(IR**IPN(I))
  !        where IR is the internal radix of the computer arithmetic. When
  !        IPN(I) = 0 the value of the normalized Legendre polynomial is
  !        contained entirely in SPN(I) and subsequent single-precision
  !        computations can be performed without further consideration of
  !        extended-range arithmetic. However, if IPN(I) .NE. 0 the corre-
  !        sponding value of the normalized Legendre polynomial cannot be
  !        represented in single-precision because of overflow or under-
  !        flow. THE USER MUST TEST IPN(I) IN HIS/HER PROGRAM. In the case
  !        that IPN(I) is nonzero, the user should try using double pre-
  !        cision if it has a wider exponent range. If double precision
  !        fails, the user could rewrite his/her program to use extended-
  !        range arithmetic.
  !
  !        The interpretation of (SPN(I),IPN(I)) can be changed to
  !        SPN(I)*(10**IPN(I)) by calling the extended-range subroutine
  !        XCON. This should be done before printing the computed values.
  !        As an example of usage, the Fortran coding
  !              J = K
  !              DO 20 I = 1, K
  !              CALL XCON(SPN(I), IPN(I),IERROR)
  !              IF (IERROR.NE.0) RETURN
  !              PRINT 10, SPN(I), IPN(I)
  !           10 FORMAT(1X, E30.8, I15)
  !              IF ((IPN(I) .EQ. 0) .OR. (J .LT. K)) GO TO 20
  !              J = I - 1
  !           20 CONTINUE
  !        will print all computed values and determine the largest J
  !        such that IPN(1) = IPN(2) = ... = IPN(J) = 0. Because of the
  !        change of representation caused by calling XCON, (SPN(I),
  !        IPN(I)) for I = J+1, J+2, ... cannot be used in subsequent
  !        extended-range computations.
  !
  !        IERROR is an error indicator. If no errors are detected,
  !        IERROR=0 when control returns to the calling routine. If
  !        an error is detected, IERROR is returned as nonzero. The
  !        calling routine must check the value of IERROR.
  !
  !        If IERROR=112 or 113, invalid input was provided to XNRMP.
  !        If IERROR=101,102,103, or 104, invalid input was provided
  !        to XSET.
  !        If IERROR=105 or 106, an internal consistency error occurred
  !        in XSET (probably due to a software malfunction in the
  !        library routine I1MACH).
  !        If IERROR=107, an overflow or underflow of an extended-range
  !        number was detected in XADJ.
  !        If IERROR=108, an overflow or underflow of an extended-range
  !        number was detected in XC210.
  !
  !***
  ! **See also:**  XSET
  !***
  ! **References:**  Smith, Olver and Lozier, Extended-Range Arithmetic and
  !                 Normalized Legendre Polynomials, ACM Trans on Math
  !                 Softw, v 7, n 1, March 1981, pp 93--105.
  !***
  ! **Routines called:**  XADD, XADJ, XERMSG, XRED, XSET

  !* REVISION HISTORY  (YYMMDD)
  !   820712  DATE WRITTEN
  !   881020  Revised to meet SLATEC CML recommendations.  (DWL and JMS)
  !   901019  Revisions to prologue.  (DWL and WRB)
  !   901106  Changed all specific intrinsics to generic.  (WRB)
  !           Corrected order of sections in prologue and added TYPE
  !           section.  (WRB)
  !           CALLs to XERROR changed to CALLs to XERMSG.  (WRB)
  !   920127  Revised PURPOSE section of prologue.  (DWL)
  USE service, ONLY : XERMSG
  INTEGER i, Ierror, ip, ip1, ip2, j, k, mu
  INTEGER Nu, Mu1, Mu2, Mode, Ipn(*), Isig
  REAL Sarg, Spn(*)
  REAL c1, c2, p, p1, p2, p3, s, sx, t, tx, x, rk
  ! CALL XSET TO INITIALIZE EXTENDED-RANGE ARITHMETIC (SEE XSET
  ! LISTING FOR DETAILS)
  !* FIRST EXECUTABLE STATEMENT  XNRMP
  Ierror = 0
  CALL XSET(0,0,0.0,0,Ierror)
  IF ( Ierror/=0 ) RETURN
  !
  !        TEST FOR PROPER INPUT VALUES.
  !
  IF ( Nu<0 ) GOTO 300
  IF ( Mu1<0 ) GOTO 300
  IF ( Mu1>Mu2 ) GOTO 300
  IF ( Nu==0 ) GOTO 200
  IF ( Mode<1.OR.Mode>2 ) GOTO 300
  IF ( Mode==2 ) THEN
    IF ( ABS(Sarg)>4.0*ATAN(1.0) ) GOTO 400
    IF ( Sarg==0.0 ) GOTO 200
    x = COS(Sarg)
    sx = ABS(SIN(Sarg))
    tx = x/sx
    Isig = INT( LOG10(2.0*Nu*(5.0+ABS(Sarg*tx))) )
  ELSE
    IF ( ABS(Sarg)>1.0 ) GOTO 400
    IF ( ABS(Sarg)==1.0 ) GOTO 200
    x = Sarg
    sx = SQRT((1.0+ABS(x))*((0.5-ABS(x))+0.5))
    tx = x/sx
    Isig = INT( LOG10(2.0*Nu*(5.0+tx**2)) )
  END IF
  !
  !        BEGIN CALCULATION
  !
  mu = Mu2
  i = Mu2 - Mu1 + 1
  !
  !        IF MU.GT.NU, NORMALIZED LEGENDRE(NU,MU,X)=0.
  !
  DO WHILE ( mu>Nu )
    Spn(i) = 0.0
    Ipn(i) = 0
    i = i - 1
    mu = mu - 1
    IF ( i<=0 ) THEN
      Isig = 0
      RETURN
    END IF
  END DO
  mu = Nu
  !
  !        P1 = 0. = NORMALIZED LEGENDRE(NU,NU+1,X)
  !
  p1 = 0.0
  ip1 = 0
  !
  !        CALCULATE P2 = NORMALIZED LEGENDRE(NU,NU,X)
  !
  p2 = 1.0
  ip2 = 0
  p3 = 0.5
  rk = 2.0
  DO j = 1, Nu
    p3 = ((rk+1.0)/rk)*p3
    p2 = p2*sx
    CALL XADJ(p2,ip2,Ierror)
    IF ( Ierror/=0 ) RETURN
    rk = rk + 2.0
  END DO
  p2 = p2*SQRT(p3)
  CALL XADJ(p2,ip2,Ierror)
  IF ( Ierror/=0 ) RETURN
  s = 2.0*tx
  t = 1.0/Nu
  IF ( Mu2>=Nu ) THEN
    Spn(i) = p2
    Ipn(i) = ip2
    i = i - 1
    IF ( i==0 ) GOTO 500
  END IF
  !
  !        RECURRENCE PROCESS
  !
  100  p = mu*t
  c1 = 1.0/SQRT((1.0-p+t)*(1.0+p))
  c2 = s*p*c1*p2
  c1 = -SQRT((1.0+p+t)*(1.0-p))*c1*p1
  CALL XADD(c2,ip2,c1,ip1,p,ip,Ierror)
  IF ( Ierror/=0 ) RETURN
  mu = mu - 1
  IF ( mu<=Mu2 ) THEN
    !
    !        STORE IN ARRAY SPN FOR RETURN TO CALLING ROUTINE.
    !
    Spn(i) = p
    Ipn(i) = ip
    i = i - 1
    IF ( i==0 ) GOTO 500
  END IF
  p1 = p2
  ip1 = ip2
  p2 = p
  ip2 = ip
  IF ( mu>Mu1 ) GOTO 100
  GOTO 500
  !
  !        SPECIAL CASE WHEN X=-1 OR +1, OR NU=0.
  !
  200  k = Mu2 - Mu1 + 1
  DO i = 1, k
    Spn(i) = 0.0
    Ipn(i) = 0
  END DO
  Isig = 0
  IF ( Mu1<=0 ) THEN
    Isig = 1
    Spn(1) = SQRT(Nu+0.5)
    Ipn(1) = 0
    IF ( MOD(Nu,2)/=0 ) THEN
      IF ( Mode/=1.OR.Sarg/=1.0 ) THEN
        IF ( Mode/=2 ) Spn(1) = -Spn(1)
      END IF
    END IF
  END IF
  RETURN
  !
  !          ERROR PRINTOUTS AND TERMINATION.
  !
  300  CALL XERMSG('SLATEC','XNRMP','NU, MU1, MU2 or MODE not valid',112,1)
  Ierror = 112
  RETURN
  400  CALL XERMSG('SLATEC','XNRMP','SARG out of range',113,1)
  Ierror = 113
  RETURN
  !
  !        RETURN TO CALLING PROGRAM
  !
  500  k = Mu2 - Mu1 + 1
  DO i = 1, k
    CALL XRED(Spn(i),Ipn(i),Ierror)
    IF ( Ierror/=0 ) RETURN
  END DO
  RETURN
END SUBROUTINE XNRMP
