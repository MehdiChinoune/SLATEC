!DECK DXNRMP
SUBROUTINE DXNRMP(Nu,Mu1,Mu2,Darg,Mode,Dpn,Ipn,Isig,Ierror)
  IMPLICIT NONE
  INTEGER i, Ierror, ip, ip1, ip2, j, k, mu
  !***BEGIN PROLOGUE  DXNRMP
  !***PURPOSE  Compute normalized Legendre polynomials.
  !***LIBRARY   SLATEC
  !***CATEGORY  C3A2, C9
  !***TYPE      DOUBLE PRECISION (XNRMP-S, DXNRMP-D)
  !***KEYWORDS  LEGENDRE FUNCTIONS
  !***AUTHOR  Lozier, Daniel W., (National Bureau of Standards)
  !           Smith, John M., (NBS and George Mason University)
  !***DESCRIPTION
  !
  !        SUBROUTINE TO CALCULATE NORMALIZED LEGENDRE POLYNOMIALS
  !        (XNRMP is single-precision version)
  !        DXNRMP calculates normalized Legendre polynomials of varying
  !        order and fixed argument and degree. The order MU and degree
  !        NU are non-negative integers and the argument is real. Because
  !        the algorithm requires the use of numbers outside the normal
  !        machine range, this subroutine employs a special arithmetic
  !        called extended-range arithmetic. See J.M. Smith, F.W.J. Olver,
  !        and D.W. Lozier, Extended-Range Arithmetic and Normalized
  !        Legendre Polynomials, ACM Transactions on Mathematical Soft-
  !        ware, 93-105, March 1981, for a complete description of the
  !        algorithm and special arithmetic. Also see program comments
  !        in DXSET.
  !
  !        The normalized Legendre polynomials are multiples of the
  !        associated Legendre polynomials of the first kind where the
  !        normalizing coefficients are chosen so as to make the integral
  !        from -1 to 1 of the square of each function equal to 1. See
  !        E. Jahnke, F. Emde and F. Losch, Tables of Higher Functions,
  !        McGraw-Hill, New York, 1960, p. 121.
  !
  !        The input values to DXNRMP are NU, MU1, MU2, DARG, and MODE.
  !        These must satisfy
  !          1. NU .GE. 0 specifies the degree of the normalized Legendre
  !             polynomial that is wanted.
  !          2. MU1 .GE. 0 specifies the lowest-order normalized Legendre
  !             polynomial that is wanted.
  !          3. MU2 .GE. MU1 specifies the highest-order normalized Leg-
  !             endre polynomial that is wanted.
  !         4a. MODE = 1 and -1.0D0 .LE. DARG .LE. 1.0D0 specifies that
  !             Normalized Legendre(NU, MU, DARG) is wanted for MU = MU1,
  !             MU1 + 1, ..., MU2.
  !         4b. MODE = 2 and -3.14159... .LT. DARG .LT. 3.14159... spec-
  !             ifies that Normalized Legendre(NU, MU, COS(DARG)) is
  !             wanted for MU = MU1, MU1 + 1, ..., MU2.
  !
  !        The output of DXNRMP consists of the two vectors DPN and IPN
  !        and the error estimate ISIG. The computed values are stored as
  !        extended-range numbers such that
  !             (DPN(1),IPN(1))=NORMALIZED LEGENDRE(NU,MU1,DX)
  !             (DPN(2),IPN(2))=NORMALIZED LEGENDRE(NU,MU1+1,DX)
  !                .
  !                .
  !             (DPN(K),IPN(K))=NORMALIZED LEGENDRE(NU,MU2,DX)
  !        where K = MU2 - MU1 + 1 and DX = DARG or COS(DARG) according
  !        to whether MODE = 1 or 2. Finally, ISIG is an estimate of the
  !        number of decimal digits lost through rounding errors in the
  !        computation. For example if DARG is accurate to 12 significant
  !        decimals, then the computed function values are accurate to
  !        12 - ISIG significant decimals (except in neighborhoods of
  !        zeros).
  !
  !        The interpretation of (DPN(I),IPN(I)) is DPN(I)*(IR**IPN(I))
  !        where IR is the internal radix of the computer arithmetic. When
  !        IPN(I) = 0 the value of the normalized Legendre polynomial is
  !        contained entirely in DPN(I) and subsequent double-precision
  !        computations can be performed without further consideration of
  !        extended-range arithmetic. However, if IPN(I) .NE. 0 the corre-
  !        sponding value of the normalized Legendre polynomial cannot be
  !        represented in double-precision because of overflow or under-
  !        flow. THE USER MUST TEST IPN(I) IN HIS/HER PROGRAM. In the case
  !        that IPN(I) is nonzero, the user could rewrite his/her program
  !        to use extended range arithmetic.
  !
  !
  !
  !        The interpretation of (DPN(I),IPN(I)) can be changed to
  !        DPN(I)*(10**IPN(I)) by calling the extended-range subroutine
  !        DXCON. This should be done before printing the computed values.
  !        As an example of usage, the Fortran coding
  !              J = K
  !              DO 20 I = 1, K
  !              CALL DXCON(DPN(I), IPN(I),IERROR)
  !              IF (IERROR.NE.0) RETURN
  !              PRINT 10, DPN(I), IPN(I)
  !           10 FORMAT(1X, D30.18, I15)
  !              IF ((IPN(I) .EQ. 0) .OR. (J .LT. K)) GO TO 20
  !              J = I - 1
  !           20 CONTINUE
  !        will print all computed values and determine the largest J
  !        such that IPN(1) = IPN(2) = ... = IPN(J) = 0. Because of the
  !        change of representation caused by calling DXCON, (DPN(I),
  !        IPN(I)) for I = J+1, J+2, ... cannot be used in subsequent
  !        extended-range computations.
  !
  !        IERROR is an error indicator. If no errors are detected,
  !        IERROR=0 when control returns to the calling routine. If
  !        an error is detected, IERROR is returned as nonzero. The
  !        calling routine must check the value of IERROR.
  !
  !        If IERROR=212 or 213, invalid input was provided to DXNRMP.
  !        If IERROR=201,202,203, or 204, invalid input was provided
  !        to DXSET.
  !        If IERROR=205 or 206, an internal consistency error occurred
  !        in DXSET (probably due to a software malfunction in the
  !        library routine I1MACH).
  !        If IERROR=207, an overflow or underflow of an extended-range
  !        number was detected in DXADJ.
  !        If IERROR=208, an overflow or underflow of an extended-range
  !        number was detected in DXC210.
  !
  !***SEE ALSO  DXSET
  !***REFERENCES  Smith, Olver and Lozier, Extended-Range Arithmetic and
  !                 Normalized Legendre Polynomials, ACM Trans on Math
  !                 Softw, v 7, n 1, March 1981, pp 93--105.
  !***ROUTINES CALLED  DXADD, DXADJ, DXRED, DXSET, XERMSG
  !***REVISION HISTORY  (YYMMDD)
  !   820712  DATE WRITTEN
  !   890126  Revised to meet SLATEC CML recommendations.  (DWL and JMS)
  !   901019  Revisions to prologue.  (DWL and WRB)
  !   901106  Changed all specific intrinsics to generic.  (WRB)
  !           Corrected order of sections in prologue and added TYPE
  !           section.  (WRB)
  !           CALLs to XERROR changed to CALLs to XERMSG.  (WRB)
  !   920127  Revised PURPOSE section of prologue.  (DWL)
  !***END PROLOGUE  DXNRMP
  INTEGER Nu, Mu1, Mu2, Mode, Ipn, Isig
  REAL(8) :: Darg, Dpn
  DIMENSION Dpn(*), Ipn(*)
  REAL(8) :: c1, c2, p, p1, p2, p3, s, sx, t, tx, x, dk
  ! CALL DXSET TO INITIALIZE EXTENDED-RANGE ARITHMETIC (SEE DXSET
  ! LISTING FOR DETAILS)
  !***FIRST EXECUTABLE STATEMENT  DXNRMP
  Ierror = 0
  CALL DXSET(0,0,0.0D0,0,Ierror)
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
    IF ( ABS(Darg)>4.0D0*ATAN(1.0D0) ) GOTO 400
    IF ( Darg==0.0D0 ) GOTO 200
    x = COS(Darg)
    sx = ABS(SIN(Darg))
    tx = x/sx
    Isig = INT( LOG10(2.0D0*Nu*(5.0D0+ABS(Darg*tx))) )
  ELSE
    IF ( ABS(Darg)>1.0D0 ) GOTO 400
    IF ( ABS(Darg)==1.0D0 ) GOTO 200
    x = Darg
    sx = SQRT((1.0D0+ABS(x))*((0.5D0-ABS(x))+0.5D0))
    tx = x/sx
    Isig = INT( LOG10(2.0D0*Nu*(5.0D0+tx**2)) )
  ENDIF
  !
  !        BEGIN CALCULATION
  !
  mu = Mu2
  i = Mu2 - Mu1 + 1
  !
  !        IF MU.GT.NU, NORMALIZED LEGENDRE(NU,MU,X)=0.
  !
  DO WHILE ( mu>Nu )
    Dpn(i) = 0.0D0
    Ipn(i) = 0
    i = i - 1
    mu = mu - 1
    IF ( i<=0 ) THEN
      Isig = 0
      RETURN
    ENDIF
  ENDDO
  mu = Nu
  !
  !        P1 = 0. = NORMALIZED LEGENDRE(NU,NU+1,X)
  !
  p1 = 0.0D0
  ip1 = 0
  !
  !        CALCULATE P2 = NORMALIZED LEGENDRE(NU,NU,X)
  !
  p2 = 1.0D0
  ip2 = 0
  p3 = 0.5D0
  dk = 2.0D0
  DO j = 1, Nu
    p3 = ((dk+1.0D0)/dk)*p3
    p2 = p2*sx
    CALL DXADJ(p2,ip2,Ierror)
    IF ( Ierror/=0 ) RETURN
    dk = dk + 2.0D0
  ENDDO
  p2 = p2*SQRT(p3)
  CALL DXADJ(p2,ip2,Ierror)
  IF ( Ierror/=0 ) RETURN
  s = 2.0D0*tx
  t = 1.0D0/Nu
  IF ( Mu2>=Nu ) THEN
    Dpn(i) = p2
    Ipn(i) = ip2
    i = i - 1
    IF ( i==0 ) GOTO 500
  ENDIF
  !
  !        RECURRENCE PROCESS
  !
  100  p = mu*t
  c1 = 1.0D0/SQRT((1.0D0-p+t)*(1.0D0+p))
  c2 = s*p*c1*p2
  c1 = -SQRT((1.0D0+p+t)*(1.0D0-p))*c1*p1
  CALL DXADD(c2,ip2,c1,ip1,p,ip,Ierror)
  IF ( Ierror/=0 ) RETURN
  mu = mu - 1
  IF ( mu<=Mu2 ) THEN
    !
    !        STORE IN ARRAY DPN FOR RETURN TO CALLING ROUTINE.
    !
    Dpn(i) = p
    Ipn(i) = ip
    i = i - 1
    IF ( i==0 ) GOTO 500
  ENDIF
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
    Dpn(i) = 0.0D0
    Ipn(i) = 0
  ENDDO
  Isig = 0
  IF ( Mu1<=0 ) THEN
    Isig = 1
    Dpn(1) = SQRT(Nu+0.5D0)
    Ipn(1) = 0
    IF ( MOD(Nu,2)/=0 ) THEN
      IF ( Mode/=1.OR.Darg/=1.0D0 ) THEN
        IF ( Mode/=2 ) Dpn(1) = -Dpn(1)
      ENDIF
    ENDIF
  ENDIF
  RETURN
  !
  !          ERROR PRINTOUTS AND TERMINATION.
  !
  300  CALL XERMSG('SLATEC','DXNRMP','NU, MU1, MU2 or MODE not valid',212,1)
  Ierror = 212
  RETURN
  400  CALL XERMSG('SLATEC','DXNRMP','DARG out of range',213,1)
  Ierror = 213
  RETURN
  !
  !        RETURN TO CALLING PROGRAM
  !
  500  k = Mu2 - Mu1 + 1
  DO i = 1, k
    CALL DXRED(Dpn(i),Ipn(i),Ierror)
    IF ( Ierror/=0 ) RETURN
  ENDDO
  RETURN
END SUBROUTINE DXNRMP
