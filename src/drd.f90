!DECK DRD
REAL(8) FUNCTION DRD(X,Y,Z,Ier)
  IMPLICIT NONE
  !***BEGIN PROLOGUE  DRD
  !***PURPOSE  Compute the incomplete or complete elliptic integral of
  !            the 2nd kind. For X and Y nonnegative, X+Y and Z positive,
  !            DRD(X,Y,Z) = Integral from zero to infinity of
  !                                -1/2     -1/2     -3/2
  !                      (3/2)(t+X)    (t+Y)    (t+Z)    dt.
  !            If X or Y is zero, the integral is complete.
  !***LIBRARY   SLATEC
  !***CATEGORY  C14
  !***TYPE      DOUBLE PRECISION (RD-S, DRD-D)
  !***KEYWORDS  COMPLETE ELLIPTIC INTEGRAL, DUPLICATION THEOREM,
  !             INCOMPLETE ELLIPTIC INTEGRAL, INTEGRAL OF THE SECOND KIND,
  !             TAYLOR SERIES
  !***AUTHOR  Carlson, B. C.
  !             Ames Laboratory-DOE
  !             Iowa State University
  !             Ames, IA  50011
  !           Notis, E. M.
  !             Ames Laboratory-DOE
  !             Iowa State University
  !             Ames, IA  50011
  !           Pexton, R. L.
  !             Lawrence Livermore National Laboratory
  !             Livermore, CA  94550
  !***DESCRIPTION
  !
  !   1.     DRD
  !          Evaluate an INCOMPLETE (or COMPLETE) ELLIPTIC INTEGRAL
  !          of the second kind
  !          Standard FORTRAN function routine
  !          Double precision version
  !          The routine calculates an approximation result to
  !          DRD(X,Y,Z) = Integral from zero to infinity of
  !                              -1/2     -1/2     -3/2
  !                    (3/2)(t+X)    (t+Y)    (t+Z)    dt,
  !          where X and Y are nonnegative, X + Y is positive, and Z is
  !          positive.  If X or Y is zero, the integral is COMPLETE.
  !          The duplication theorem is iterated until the variables are
  !          nearly equal, and the function is then expanded in Taylor
  !          series to fifth order.
  !
  !   2.     Calling Sequence
  !
  !          DRD( X, Y, Z, IER )
  !
  !          Parameters On Entry
  !          Values assigned by the calling routine
  !
  !          X      - Double precision, nonnegative variable
  !
  !          Y      - Double precision, nonnegative variable
  !
  !                   X + Y is positive
  !
  !          Z      - Double precision, positive variable
  !
  !
  !
  !          On Return    (values assigned by the DRD routine)
  !
  !          DRD     - Double precision approximation to the integral
  !
  !
  !          IER    - Integer
  !
  !                   IER = 0 Normal and reliable termination of the
  !                           routine. It is assumed that the requested
  !                           accuracy has been achieved.
  !
  !                   IER >  0 Abnormal termination of the routine
  !
  !
  !          X, Y, Z are unaltered.
  !
  !   3.    Error Messages
  !
  !         Value of IER assigned by the DRD routine
  !
  !                  Value assigned         Error message printed
  !                  IER = 1                MIN(X,Y) .LT. 0.0D0
  !                      = 2                MIN(X + Y, Z ) .LT. LOLIM
  !                      = 3                MAX(X,Y,Z) .GT. UPLIM
  !
  !
  !   4.     Control Parameters
  !
  !                  Values of LOLIM, UPLIM, and ERRTOL are set by the
  !                  routine.
  !
  !          LOLIM and UPLIM determine the valid range of X, Y, and Z
  !
  !          LOLIM  - Lower limit of valid arguments
  !
  !                    Not less  than 2 / (machine maximum) ** (2/3).
  !
  !          UPLIM  - Upper limit of valid arguments
  !
  !                 Not greater than (0.1D0 * ERRTOL / machine
  !                 minimum) ** (2/3), where ERRTOL is described below.
  !                 In the following table it is assumed that ERRTOL will
  !                 never be chosen smaller than 1.0D-5.
  !
  !
  !                    Acceptable values for:   LOLIM      UPLIM
  !                    IBM 360/370 SERIES   :   6.0D-51     1.0D+48
  !                    CDC 6000/7000 SERIES :   5.0D-215    2.0D+191
  !                    UNIVAC 1100 SERIES   :   1.0D-205    2.0D+201
  !                    CRAY                 :   3.0D-1644   1.69D+1640
  !                    VAX 11 SERIES        :   1.0D-25     4.5D+21
  !
  !
  !          ERRTOL determines the accuracy of the answer
  !
  !                 The value assigned by the routine will result
  !                 in solution precision within 1-2 decimals of
  !                 "machine precision".
  !
  !          ERRTOL    Relative error due to truncation is less than
  !                    3 * ERRTOL ** 6 / (1-ERRTOL) ** 3/2.
  !
  !
  !
  !        The accuracy of the computed approximation to the integral
  !        can be controlled by choosing the value of ERRTOL.
  !        Truncation of a Taylor series after terms of fifth order
  !        introduces an error less than the amount shown in the
  !        second column of the following table for each value of
  !        ERRTOL in the first column.  In addition to the truncation
  !        error there will be round-off error, but in practice the
  !        total error from both sources is usually less than the
  !        amount given in the table.
  !
  !
  !
  !
  !          Sample choices:  ERRTOL   Relative truncation
  !                                    error less than
  !                           1.0D-3    4.0D-18
  !                           3.0D-3    3.0D-15
  !                           1.0D-2    4.0D-12
  !                           3.0D-2    3.0D-9
  !                           1.0D-1    4.0D-6
  !
  !
  !                    Decreasing ERRTOL by a factor of 10 yields six more
  !                    decimal digits of accuracy at the expense of one or
  !                    two more iterations of the duplication theorem.
  !
  ! *Long Description:
  !
  !   DRD Special Comments
  !
  !
  !
  !          Check: DRD(X,Y,Z) + DRD(Y,Z,X) + DRD(Z,X,Y)
  !          = 3 / SQRT(X * Y * Z), where X, Y, and Z are positive.
  !
  !
  !          On Input:
  !
  !          X, Y, and Z are the variables in the integral DRD(X,Y,Z).
  !
  !
  !          On Output:
  !
  !
  !          X, Y, Z are unaltered.
  !
  !
  !
  !          ********************************************************
  !
  !          WARNING: Changes in the program may improve speed at the
  !                   expense of robustness.
  !
  !
  !
  !    -------------------------------------------------------------------
  !
  !
  !   Special double precision functions via DRD and DRF
  !
  !
  !                  Legendre form of ELLIPTIC INTEGRAL of 2nd kind
  !
  !                  -----------------------------------------
  !
  !
  !                                             2         2   2
  !                  E(PHI,K) = SIN(PHI) DRF(COS (PHI),1-K SIN (PHI),1) -
  !
  !                     2      3             2         2   2
  !                  -(K/3) SIN (PHI) DRD(COS (PHI),1-K SIN (PHI),1)
  !
  !
  !                                  2        2            2
  !                  E(K) = DRF(0,1-K ,1) - (K/3) DRD(0,1-K ,1)
  !
  !                         PI/2     2   2      1/2
  !                       = INT  (1-K SIN (PHI) )  D PHI
  !                          0
  !
  !                  Bulirsch form of ELLIPTIC INTEGRAL of 2nd kind
  !
  !                  -----------------------------------------
  !
  !                                               2 2    2
  !                  EL2(X,KC,A,B) = AX DRF(1,1+KC X ,1+X ) +
  !
  !                                              3          2 2    2
  !                                 +(1/3)(B-A) X DRD(1,1+KC X ,1+X )
  !
  !
  !
  !
  !                  Legendre form of alternative ELLIPTIC INTEGRAL
  !                  of 2nd kind
  !
  !                  -----------------------------------------
  !
  !
  !
  !                            Q     2       2   2  -1/2
  !                  D(Q,K) = INT SIN P  (1-K SIN P)     DP
  !                            0
  !
  !
  !
  !                                     3          2     2   2
  !                  D(Q,K) = (1/3) (SIN Q) DRD(COS Q,1-K SIN Q,1)
  !
  !
  !
  !
  !                  Lemniscate constant  B
  !
  !                  -----------------------------------------
  !
  !
  !
  !
  !                       1    2    4 -1/2
  !                  B = INT  S (1-S )    DS
  !                       0
  !
  !
  !                  B = (1/3) DRD (0,2,1)
  !
  !
  !                  Heuman's LAMBDA function
  !
  !                  -----------------------------------------
  !
  !
  !
  !                  (PI/2) LAMBDA0(A,B) =
  !
  !                                    2                2
  !                 = SIN(B) (DRF(0,COS (A),1)-(1/3) SIN (A) *
  !
  !                            2               2         2       2
  !                  *DRD(0,COS (A),1)) DRF(COS (B),1-COS (A) SIN (B),1)
  !
  !                            2       3             2
  !                  -(1/3) COS (A) SIN (B) DRF(0,COS (A),1) *
  !
  !                           2         2       2
  !                   *DRD(COS (B),1-COS (A) SIN (B),1)
  !
  !
  !
  !                  Jacobi ZETA function
  !
  !                  -----------------------------------------
  !
  !                             2                 2       2   2
  !                  Z(B,K) = (K/3) SIN(B) DRF(COS (B),1-K SIN (B),1)
  !
  !
  !                                       2             2
  !                             *DRD(0,1-K ,1)/DRF(0,1-K ,1)
  !
  !                               2       3           2       2   2
  !                            -(K /3) SIN (B) DRD(COS (B),1-K SIN (B),1)
  !
  !
  ! ---------------------------------------------------------------------
  !
  !***REFERENCES  B. C. Carlson and E. M. Notis, Algorithms for incomplete
  !                 elliptic integrals, ACM Transactions on Mathematical
  !                 Software 7, 3 (September 1981), pp. 398-403.
  !               B. C. Carlson, Computing elliptic integrals by
  !                 duplication, Numerische Mathematik 33, (1979),
  !                 pp. 1-16.
  !               B. C. Carlson, Elliptic integrals of the first kind,
  !                 SIAM Journal of Mathematical Analysis 8, (1977),
  !                 pp. 231-242.
  !***ROUTINES CALLED  D1MACH, XERMSG
  !***REVISION HISTORY  (YYMMDD)
  !   790801  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900326  Removed duplicate information from DESCRIPTION section.
  !           (WRB)
  !   900510  Modify calls to XERMSG to put in standard form.  (RWC)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  !***END PROLOGUE  DRD
  CHARACTER(16) :: xern3, xern4, xern5, xern6
  INTEGER Ier
  REAL(8) :: lolim, tuplim, uplim, epslon, errtol, D1MACH
  REAL(8) :: c1, c2, c3, c4, ea, eb, ec, ed, ef, lamda
  REAL(8) :: mu, power4, sigma, s1, s2, X, xn, xndev
  REAL(8) :: xnroot, Y, yn, yndev, ynroot, Z, zn, zndev ,&
    znroot
  LOGICAL first
  SAVE errtol, lolim, uplim, c1, c2, c3, c4, first
  DATA first/.TRUE./
  !
  !***FIRST EXECUTABLE STATEMENT  DRD
  IF ( first ) THEN
    errtol = (D1MACH(3)/3.0D0)**(1.0D0/6.0D0)
    lolim = 2.0D0/(D1MACH(2))**(2.0D0/3.0D0)
    tuplim = D1MACH(1)**(1.0E0/3.0E0)
    tuplim = (0.10D0*errtol)**(1.0E0/3.0E0)/tuplim
    uplim = tuplim**2.0D0
    !
    c1 = 3.0D0/14.0D0
    c2 = 1.0D0/6.0D0
    c3 = 9.0D0/22.0D0
    c4 = 3.0D0/26.0D0
  ENDIF
  first = .FALSE.
  !
  !         CALL ERROR HANDLER IF NECESSARY.
  !
  DRD = 0.0D0
  IF ( MIN(X,Y)<0.0D0 ) THEN
    Ier = 1
    WRITE (xern3,'(1PE15.6)') X
    WRITE (xern4,'(1PE15.6)') Y
    CALL XERMSG('SLATEC','DRD','MIN(X,Y).LT.0 WHERE X = '//xern3//&
      ' AND Y = '//xern4,1,1)
    RETURN
  ENDIF
  !
  IF ( MAX(X,Y,Z)>uplim ) THEN
    Ier = 3
    WRITE (xern3,'(1PE15.6)') X
    WRITE (xern4,'(1PE15.6)') Y
    WRITE (xern5,'(1PE15.6)') Z
    WRITE (xern6,'(1PE15.6)') uplim
    CALL XERMSG('SLATEC','DRD','MAX(X,Y,Z).GT.UPLIM WHERE X = '//xern3//&
      ' Y = '//xern4//' Z = '//xern5//' AND UPLIM = '//xern6,3,1)
    RETURN
  ENDIF
  !
  IF ( MIN(X+Y,Z)<lolim ) THEN
    Ier = 2
    WRITE (xern3,'(1PE15.6)') X
    WRITE (xern4,'(1PE15.6)') Y
    WRITE (xern5,'(1PE15.6)') Z
    WRITE (xern6,'(1PE15.6)') lolim
    CALL XERMSG('SLATEC','DRD','MIN(X+Y,Z).LT.LOLIM WHERE X = '//xern3//&
      ' Y = '//xern4//' Z = '//xern5//' AND LOLIM = '//xern6,2,1)
    RETURN
  ENDIF
  !
  Ier = 0
  xn = X
  yn = Y
  zn = Z
  sigma = 0.0D0
  power4 = 1.0D0
  DO
    !
    mu = (xn+yn+3.0D0*zn)*0.20D0
    xndev = (mu-xn)/mu
    yndev = (mu-yn)/mu
    zndev = (mu-zn)/mu
    epslon = MAX(ABS(xndev),ABS(yndev),ABS(zndev))
    IF ( epslon<errtol ) THEN
      !
      ea = xndev*yndev
      eb = zndev*zndev
      ec = ea - eb
      ed = ea - 6.0D0*eb
      ef = ed + ec + ec
      s1 = ed*(-c1+0.250D0*c3*ed-1.50D0*c4*zndev*ef)
      s2 = zndev*(c2*ef+zndev*(-c3*ec+zndev*c4*ea))
      DRD = 3.0D0*sigma + power4*(1.0D0+s1+s2)/(mu*SQRT(mu))
      EXIT
    ELSE
      xnroot = SQRT(xn)
      ynroot = SQRT(yn)
      znroot = SQRT(zn)
      lamda = xnroot*(ynroot+znroot) + ynroot*znroot
      sigma = sigma + power4/(znroot*(zn+lamda))
      power4 = power4*0.250D0
      xn = (xn+lamda)*0.250D0
      yn = (yn+lamda)*0.250D0
      zn = (zn+lamda)*0.250D0
    ENDIF
  ENDDO
  !
END FUNCTION DRD
