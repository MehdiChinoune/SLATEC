!** RD
REAL FUNCTION RD(X,Y,Z,Ier)
  !>
  !***
  !  Compute the incomplete or complete elliptic integral of the
  !            2nd kind.  For X and Y nonnegative, X+Y and Z positive,
  !             RD(X,Y,Z) = Integral from zero to infinity of
  !                                -1/2     -1/2     -3/2
  !                      (3/2)(t+X)    (t+Y)    (t+Z)    dt.
  !            If X or Y is zero, the integral is complete.
  !***
  ! **Library:**   SLATEC
  !***
  ! **Category:**  C14
  !***
  ! **Type:**      SINGLE PRECISION (RD-S, DRD-D)
  !***
  ! **Keywords:**  COMPLETE ELLIPTIC INTEGRAL, DUPLICATION THEOREM,
  !             INCOMPLETE ELLIPTIC INTEGRAL, INTEGRAL OF THE SECOND KIND,
  !             TAYLOR SERIES
  !***
  ! **Author:**  Carlson, B. C.
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
  !***
  ! **Description:**
  !
  !   1.     RD
  !          Evaluate an INCOMPLETE (or COMPLETE) ELLIPTIC INTEGRAL
  !          of the second kind
  !          Standard FORTRAN function routine
  !          Single precision version
  !          The routine calculates an approximation result to
  !          RD(X,Y,Z) = Integral from zero to infinity of
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
  !          RD( X, Y, Z, IER )
  !
  !          Parameters on Entry
  !          Values assigned by the calling routine
  !
  !          X      - Single precision, nonnegative variable
  !
  !          Y      - Single precision, nonnegative variable
  !
  !                   X + Y is positive
  !
  !          Z      - Real, positive variable
  !
  !
  !
  !          On Return     (values assigned by the RD routine)
  !
  !          RD     - Real approximation to the integral
  !
  !
  !          IER    - Integer
  !
  !                   IER = 0 Normal and reliable termination of the
  !                           routine.  It is assumed that the requested
  !                           accuracy has been achieved.
  !
  !                   IER >  0 Abnormal termination of the routine
  !
  !
  !          X, Y, Z are unaltered.
  !
  !   3.    Error Messages
  !
  !         Value of IER assigned by the RD routine
  !
  !                  Value Assigned         Error Message Printed
  !                  IER = 1                MIN(X,Y) .LT. 0.0E0
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
  !                    Not greater than (0.1E0 * ERRTOL / machine
  !                    minimum) ** (2/3), where ERRTOL is described below.
  !                    In the following table it is assumed that ERRTOL
  !                    will never be chosen smaller than 1.0E-5.
  !
  !
  !                    Acceptable Values For:   LOLIM      UPLIM
  !                    IBM 360/370 SERIES   :   6.0E-51     1.0E+48
  !                    CDC 6000/7000 SERIES :   5.0E-215    2.0E+191
  !                    UNIVAC 1100 SERIES   :   1.0E-25     2.0E+21
  !                    CRAY                 :   3.0E-1644   1.69E+1640
  !                    VAX 11 SERIES        :   1.0E-25     4.5E+21
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
  !              The accuracy of the computed approximation to the inte-
  !              gral can be controlled by choosing the value of ERRTOL.
  !              Truncation of a Taylor series after terms of fifth order
  !              introduces an error less than the amount shown in the
  !              second column of the following table for each value of
  !              ERRTOL in the first column.  In addition to the trunca-
  !              tion error there will be round-off error, but in prac-
  !              tice the total error from both sources is usually less
  !              than the amount given in the table.
  !
  !
  !
  !
  !          Sample Choices:  ERRTOL   Relative Truncation
  !                                    error less than
  !                           1.0E-3    4.0E-18
  !                           3.0E-3    3.0E-15
  !                           1.0E-2    4.0E-12
  !                           3.0E-2    3.0E-9
  !                           1.0E-1    4.0E-6
  !
  !
  !                    Decreasing ERRTOL by a factor of 10 yields six more
  !                    decimal digits of accuracy at the expense of one or
  !                    two more iterations of the duplication theorem.
  !
  !- Long Description:
  !
  !   RD Special Comments
  !
  !
  !
  !          Check: RD(X,Y,Z) + RD(Y,Z,X) + RD(Z,X,Y)
  !          = 3 /  SQRT(X * Y * Z), where X, Y, and Z are positive.
  !
  !
  !          On Input:
  !
  !          X, Y, and Z are the variables in the integral RD(X,Y,Z).
  !
  !
  !          On Output:
  !
  !
  !          X, Y, and Z are unaltered.
  !
  !
  !
  !          ********************************************************
  !
  !           WARNING: Changes in the program may improve speed at the
  !                    expense of robustness.
  !
  !
  !
  !    -------------------------------------------------------------------
  !
  !
  !   Special Functions via RD and RF
  !
  !
  !                  Legendre form of ELLIPTIC INTEGRAL of 2nd kind
  !                  ----------------------------------------------
  !
  !
  !                                            2         2   2
  !                  E(PHI,K) = SIN(PHI) RF(COS (PHI),1-K SIN (PHI),1) -
  !
  !                     2      3            2         2   2
  !                  -(K/3) SIN (PHI) RD(COS (PHI),1-K SIN (PHI),1)
  !
  !
  !                                 2        2           2
  !                  E(K) = RF(0,1-K ,1) - (K/3) RD(0,1-K ,1)
  !
  !
  !                         PI/2     2   2      1/2
  !                       = INT  (1-K SIN (PHI) )  D PHI
  !                          0
  !
  !
  !
  !                  Bulirsch form of ELLIPTIC INTEGRAL of 2nd kind
  !                  ----------------------------------------------
  !
  !                                              2 2    2
  !                  EL2(X,KC,A,B) = AX RF(1,1+KC X ,1+X ) +
  !
  !                                              3         2 2    2
  !                                 +(1/3)(B-A) X RD(1,1+KC X ,1+X )
  !
  !
  !
  !                  Legendre form of alternative ELLIPTIC INTEGRAL of 2nd
  !                  -----------------------------------------------------
  !                        kind
  !                        ----
  !
  !                            Q     2       2   2  -1/2
  !                  D(Q,K) = INT SIN P  (1-K SIN P)     DP
  !                            0
  !
  !
  !
  !                                   3          2     2   2
  !                  D(Q,K) =(1/3)(SIN Q)  RD(COS Q,1-K SIN Q,1)
  !
  !
  !
  !
  !
  !                  Lemniscate constant B
  !                  ---------------------
  !
  !
  !
  !                       1    2    4 -1/2
  !                  B = INT  S (1-S )    DS
  !                       0
  !
  !
  !                  B =(1/3)RD (0,2,1)
  !
  !
  !
  !
  !                  Heuman's LAMBDA function
  !                  ------------------------
  !
  !
  !
  !                  (PI/2) LAMBDA0(A,B) =
  !
  !                                       2                2
  !                     = SIN(B) (RF(0,COS (A),1)-(1/3) SIN (A) *
  !
  !                               2              2         2       2
  !                      *RD(0,COS (A),1)) RF(COS (B),1-COS (A) SIN (B),1)
  !
  !                               2       3            2
  !                     -(1/3) COS (A) SIN (B) RF(0,COS (A),1) *
  !
  !                             2         2       2
  !                      *RD(COS (B),1-COS (A) SIN (B),1)
  !
  !
  !
  !                  Jacobi ZETA function
  !                  --------------------
  !
  !
  !                             2                2       2   2
  !                  Z(B,K) = (K/3) SIN(B) RF(COS (B),1-K SIN (B),1)
  !
  !
  !                                      2            2
  !                             *RD(0,1-K ,1)/RF(0,1-K ,1)
  !
  !                               2       3          2       2   2
  !                            -(K /3) SIN (B) RD(COS (B),1-K SIN (B),1)
  !
  !
  !    -------------------------------------------------------------------
  !
  !***
  ! **References:**  B. C. Carlson and E. M. Notis, Algorithms for incomplete
  !                 elliptic integrals, ACM Transactions on Mathematical
  !                 Software 7, 3 (September 1981), pp. 398-403.
  !               B. C. Carlson, Computing elliptic integrals by
  !                 duplication, Numerische Mathematik 33, (1979),
  !                 pp. 1-16.
  !               B. C. Carlson, Elliptic integrals of the first kind,
  !                 SIAM Journal of Mathematical Analysis 8, (1977),
  !                 pp. 231-242.
  !***
  ! **Routines called:**  R1MACH, XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   790801  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900326  Removed duplicate information from DESCRIPTION section.
  !           (WRB)
  !   900510  Modify calls to XERMSG to put in standard form.  (RWC)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  USE service, ONLY : XERMSG, R1MACH
  CHARACTER(16) :: xern3, xern4, xern5, xern6
  INTEGER Ier
  REAL epslon, ea, eb, ec, ed, ef, lamda, mu, power4, sigma, s1, s2, X, xn, &
    xndev, xnroot, Y, yn, yndev, ynroot, Z, zn, zndev, znroot, tuplim
  REAL, SAVE :: errtol, lolim, uplim
  REAL, PARAMETER :: c1 = 3.0E0/14.0E0, c2 = 1.0E0/6.0E0, c3 = 9.0E0/22.0E0, &
    c4 = 3.0E0/26.0E0
  LOGICAL :: first = .TRUE.
  !
  !* FIRST EXECUTABLE STATEMENT  RD
  IF ( first ) THEN
    errtol = (R1MACH(3)/3.0E0)**(1.0E0/6.0E0)
    lolim = 2.0E0/(R1MACH(2))**(2.0E0/3.0E0)
    tuplim = R1MACH(1)**(1.0E0/3.0E0)
    tuplim = (0.10E0*errtol)**(1.0E0/3.0E0)/tuplim
    uplim = tuplim**2.0E0
    first = .FALSE.
  END IF
  !
  !         CALL ERROR HANDLER IF NECESSARY.
  !
  RD = 0.0E0
  IF ( MIN(X,Y)<0.0E0 ) THEN
    Ier = 1
    WRITE (xern3,'(1PE15.6)') X
    WRITE (xern4,'(1PE15.6)') Y
    CALL XERMSG('SLATEC','RD','MIN(X,Y).LT.0 WHERE X = '//xern3//&
      ' AND Y = '//xern4,1,1)
    RETURN
  END IF
  !
  IF ( MAX(X,Y,Z)>uplim ) THEN
    Ier = 3
    WRITE (xern3,'(1PE15.6)') X
    WRITE (xern4,'(1PE15.6)') Y
    WRITE (xern5,'(1PE15.6)') Z
    WRITE (xern6,'(1PE15.6)') uplim
    CALL XERMSG('SLATEC','RD','MAX(X,Y,Z).GT.UPLIM WHERE X = '//xern3//&
      ' Y = '//xern4//' Z = '//xern5//' AND UPLIM = '//xern6,3,1)
    RETURN
  END IF
  !
  IF ( MIN(X+Y,Z)<lolim ) THEN
    Ier = 2
    WRITE (xern3,'(1PE15.6)') X
    WRITE (xern4,'(1PE15.6)') Y
    WRITE (xern5,'(1PE15.6)') Z
    WRITE (xern6,'(1PE15.6)') lolim
    CALL XERMSG('SLATEC','RD','MIN(X+Y,Z).LT.LOLIM WHERE X = '//xern3//&
      ' Y = '//xern4//' Z = '//xern5//' AND LOLIM = '//xern6,2,1)
    RETURN
  END IF
  !
  Ier = 0
  xn = X
  yn = Y
  zn = Z
  sigma = 0.0E0
  power4 = 1.0E0
  DO
    !
    mu = (xn+yn+3.0E0*zn)*0.20E0
    xndev = (mu-xn)/mu
    yndev = (mu-yn)/mu
    zndev = (mu-zn)/mu
    epslon = MAX(ABS(xndev),ABS(yndev),ABS(zndev))
    IF ( epslon<errtol ) THEN
      !
      ea = xndev*yndev
      eb = zndev*zndev
      ec = ea - eb
      ed = ea - 6.0E0*eb
      ef = ed + ec + ec
      s1 = ed*(-c1+0.250E0*c3*ed-1.50E0*c4*zndev*ef)
      s2 = zndev*(c2*ef+zndev*(-c3*ec+zndev*c4*ea))
      RD = 3.0E0*sigma + power4*(1.0E0+s1+s2)/(mu*SQRT(mu))
      EXIT
    ELSE
      xnroot = SQRT(xn)
      ynroot = SQRT(yn)
      znroot = SQRT(zn)
      lamda = xnroot*(ynroot+znroot) + ynroot*znroot
      sigma = sigma + power4/(znroot*(zn+lamda))
      power4 = power4*0.250E0
      xn = (xn+lamda)*0.250E0
      yn = (yn+lamda)*0.250E0
      zn = (zn+lamda)*0.250E0
    END IF
  END DO
  !
END FUNCTION RD
