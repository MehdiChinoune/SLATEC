!** RJ
REAL FUNCTION RJ(X,Y,Z,P,Ier)
  IMPLICIT NONE
  !>
  !***
  !  Compute the incomplete or complete (X or Y or Z is zero)
  !            elliptic integral of the 3rd kind.  For X, Y, and Z non-
  !            negative, at most one of them zero, and P positive,
  !             RJ(X,Y,Z,P) = Integral from zero to infinity of
  !                                  -1/2     -1/2     -1/2     -1
  !                        (3/2)(t+X)    (t+Y)    (t+Z)    (t+P)  dt.
  !***
  ! **Library:**   SLATEC
  !***
  ! **Category:**  C14
  !***
  ! **Type:**      SINGLE PRECISION (RJ-S, DRJ-D)
  !***
  ! **Keywords:**  COMPLETE ELLIPTIC INTEGRAL, DUPLICATION THEOREM,
  !             INCOMPLETE ELLIPTIC INTEGRAL, INTEGRAL OF THE THIRD KIND,
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
  !   1.     RJ
  !          Standard FORTRAN function routine
  !          Single precision version
  !          The routine calculates an approximation result to
  !          RJ(X,Y,Z,P) = Integral from zero to infinity of
  !
  !                                -1/2     -1/2     -1/2     -1
  !                      (3/2)(t+X)    (t+Y)    (t+Z)    (t+P)  dt,
  !
  !          where X, Y, and Z are nonnegative, at most one of them is
  !          zero, and P is positive.  If X or Y or Z is zero, the
  !          integral is COMPLETE.  The duplication theorem is iterated
  !          until the variables are nearly equal, and the function is
  !          then expanded in Taylor series to fifth order.
  !
  !
  !   2.     Calling Sequence
  !          RJ( X, Y, Z, P, IER )
  !
  !          Parameters On Entry
  !          Values assigned by the calling routine
  !
  !          X      - Single precision, nonnegative variable
  !
  !          Y      - Single precision, nonnegative variable
  !
  !          Z      - Single precision, nonnegative variable
  !
  !          P      - Single precision, positive variable
  !
  !
  !          On  Return     (values assigned by the RJ routine)
  !
  !          RJ     - Single precision approximation to the integral
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
  !          X, Y, Z, P are unaltered.
  !
  !
  !   3.    Error Messages
  !
  !         Value of IER assigned by the RJ routine
  !
  !                  Value Assigned        Error Message Printed
  !                  IER = 1               MIN(X,Y,Z) .LT. 0.0E0
  !                      = 2               MIN(X+Y,X+Z,Y+Z,P) .LT. LOLIM
  !                      = 3               MAX(X,Y,Z,P) .GT. UPLIM
  !
  !
  !
  !   4.     Control Parameters
  !
  !                  Values of LOLIM, UPLIM, and ERRTOL are set by the
  !                  routine.
  !
  !
  !          LOLIM and UPLIM determine the valid range of X Y, Z, and P
  !
  !          LOLIM is not less than the cube root of the value
  !          of LOLIM used in the routine for RC.
  !
  !          UPLIM is not greater than 0.3 times the cube root of
  !          the value of UPLIM used in the routine for RC.
  !
  !
  !                     Acceptable Values For:   LOLIM      UPLIM
  !                     IBM 360/370 SERIES   :   2.0E-26     3.0E+24
  !                     CDC 6000/7000 SERIES :   5.0E-98     3.0E+106
  !                     UNIVAC 1100 SERIES   :   5.0E-13     6.0E+11
  !                     CRAY                 :   1.32E-822   1.4E+821
  !                     VAX 11 SERIES        :   2.5E-13     9.0E+11
  !
  !
  !
  !          ERRTOL determines the accuracy of the answer
  !
  !                 The value assigned by the routine will result
  !                 in solution precision within 1-2 decimals of
  !                 "machine precision".
  !
  !
  !
  !
  !          Relative error due to truncation of the series for RJ
  !          is less than 3 * ERRTOL ** 6 / (1 - ERRTOL) ** 3/2.
  !
  !
  !
  !              The accuracy of the computed approximation to the inte-
  !              gral can be controlled by choosing the value of ERRTOL.
  !              Truncation of a Taylor series after terms of fifth order
  !              Introduces an error less than the amount shown in the
  !              second column of the following table for each value of
  !              ERRTOL in the first column.  In addition to the trunca-
  !              tion error there will be round-off error, but in prac-
  !              tice the total error from both sources is usually less
  !              than the amount given in the table.
  !
  !
  !
  !          Sample choices:  ERRTOL   Relative Truncation
  !                                    error less than
  !                           1.0E-3    4.0E-18
  !                           3.0E-3    3.0E-15
  !                           1.0E-2    4.0E-12
  !                           3.0E-2    3.0E-9
  !                           1.0E-1    4.0E-6
  !
  !                    Decreasing ERRTOL by a factor of 10 yields six more
  !                    decimal digits of accuracy at the expense of one or
  !                    two more iterations of the duplication theorem.
  !
  !- Long Description:
  !
  !   RJ Special Comments
  !
  !
  !          Check by addition theorem: RJ(X,X+Z,X+W,X+P)
  !          + RJ(Y,Y+Z,Y+W,Y+P) + (A-B) * RJ(A,B,B,A) + 3 / SQRT(A)
  !          = RJ(0,Z,W,P), where X,Y,Z,W,P are positive and X * Y
  !          = Z * W,  A = P * P * (X+Y+Z+W),  B = P * (P+X) * (P+Y),
  !          and B - A = P * (P-Z) * (P-W).  The sum of the third and
  !          fourth terms on the left side is 3 * RC(A,B).
  !
  !
  !          On Input:
  !
  !          X, Y, Z, and P are the variables in the integral RJ(X,Y,Z,P).
  !
  !
  !          On Output:
  !
  !
  !          X, Y, Z, and P are unaltered.
  !
  !          ********************************************************
  !
  !          Warning: Changes in the program may improve speed at the
  !                   expense of robustness.
  !
  ! ------------------------------------------------------------
  !
  !
  !   Special Functions via RJ and RF
  !
  !
  !                  Legendre form of ELLIPTIC INTEGRAL of 3rd kind
  !                  ----------------------------------------------
  !
  !
  !                               PHI         2         -1
  !                  P(PHI,K,N) = INT (1+N SIN (THETA) )   *
  !                                0
  !
  !                                      2    2         -1/2
  !                                 *(1-K  SIN (THETA) )     D THETA
  !
  !
  !                                         2          2   2
  !                       = SIN (PHI) RF(COS (PHI), 1-K SIN (PHI),1)
  !
  !                                  3            2         2   2
  !                        -(N/3) SIN (PHI) RJ(COS (PHI),1-K SIN (PHI),
  !
  !                                 2
  !                        1,1+N SIN (PHI))
  !
  !
  !
  !                  Bulirsch form of ELLIPTIC INTEGRAL of 3rd kind
  !                  ----------------------------------------------
  !
  !
  !                                           2 2    2
  !                  EL3(X,KC,P) = X RF(1,1+KC X ,1+X ) +
  !
  !                                            3          2 2    2     2
  !                               +(1/3)(1-P) X  RJ(1,1+KC X ,1+X ,1+PX )
  !
  !
  !                                           2
  !                  CEL(KC,P,A,B) = A RF(0,KC ,1) +
  !
  !                                                     2
  !                                 +(1/3)(B-PA) RJ(0,KC ,1,P)
  !
  !
  !
  !
  !                  Heuman's LAMBDA function
  !                  ------------------------
  !
  !
  !                                 2                     2      2    1/2
  !                  L(A,B,P) = (COS(A)SIN(B)COS(B)/(1-COS (A)SIN (B))   )
  !
  !                                           2         2       2
  !                            *(SIN(P) RF(COS (P),1-SIN (A) SIN (P),1)
  !
  !                                 2       3            2       2
  !                            +(SIN (A) SIN (P)/(3(1-COS (A) SIN (B))))
  !
  !                                   2         2       2
  !                            *RJ(COS (P),1-SIN (A) SIN (P),1,1-
  !
  !                                2       2          2       2
  !                            -SIN (A) SIN (P)/(1-COS (A) SIN (B))))
  !
  !
  !
  !
  !                  (PI/2) LAMBDA0(A,B) =L(A,B,PI/2) =
  !
  !
  !                    2                         2       2    -1/2
  !               = COS (A)  SIN(B) COS(B) (1-COS (A) SIN (B))
  !
  !                           2                  2       2
  !                  *RF(0,COS (A),1) + (1/3) SIN (A) COS (A)
  !
  !                                       2       2    -3/2
  !                  *SIN(B) COS(B) (1-COS (A) SIN (B))
  !
  !                           2         2       2          2       2
  !                  *RJ(0,COS (A),1,COS (A) COS (B)/(1-COS (A) SIN (B)))
  !
  !
  !
  !                  Jacobi ZETA function
  !                  --------------------
  !
  !
  !                             2                     2   2    1/2
  !                  Z(B,K) = (K/3) SIN(B) COS(B) (1-K SIN (B))
  !
  !
  !                                      2      2   2                2
  !                             *RJ(0,1-K ,1,1-K SIN (B)) / RF (0,1-K ,1)
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
  ! **Routines called:**  R1MACH, RC, XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   790801  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   891009  Removed unreferenced statement labels.  (WRB)
  !   891009  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900326  Removed duplicate information from DESCRIPTION section.
  !           (WRB)
  !   900510  Changed calls to XERMSG to standard form, and some
  !           editorial changes.  (RWC)).
  !   920501  Reformatted the REFERENCES section.  (WRB)

  REAL R1MACH
  CHARACTER(16) :: xern3, xern4, xern5, xern6, xern7
  INTEGER Ier
  REAL alfa, beta, ea, eb, ec, e2, e3, epslon, lamda, mu, P, pn, pndev, power4, &
    RC, sigma, s1, s2, s3, X, xn, xndev, xnroot, Y, yn, yndev, ynroot, Z, zn, &
    zndev, znroot
  REAL, SAVE :: errtol, lolim, uplim
  REAL, PARAMETER :: c1 = 3.0E0/14.0E0, c2 = 1.0E0/3.0E0, c3 = 3.0E0/22.0E0, &
    c4 = 3.0E0/26.0E0
  LOGICAL :: first = .TRUE.
  !
  !* FIRST EXECUTABLE STATEMENT  RJ
  IF ( first ) THEN
    errtol = (R1MACH(3)/3.0E0)**(1.0E0/6.0E0)
    lolim = (5.0E0*R1MACH(1))**(1.0E0/3.0E0)
    uplim = 0.30E0*(R1MACH(2)/5.0E0)**(1.0E0/3.0E0)
    first = .FALSE.
  END IF
  !
  !         CALL ERROR HANDLER IF NECESSARY.
  !
  RJ = 0.0E0
  IF ( MIN(X,Y,Z)<0.0E0 ) THEN
    Ier = 1
    WRITE (xern3,'(1PE15.6)') X
    WRITE (xern4,'(1PE15.6)') Y
    WRITE (xern5,'(1PE15.6)') Z
    CALL XERMSG('SLATEC','RJ','MIN(X,Y,Z).LT.0 WHERE X = '//xern3//' Y = '//&
      xern4//' AND Z = '//xern5,1,1)
    RETURN
  END IF
  !
  IF ( MAX(X,Y,Z,P)>uplim ) THEN
    Ier = 3
    WRITE (xern3,'(1PE15.6)') X
    WRITE (xern4,'(1PE15.6)') Y
    WRITE (xern5,'(1PE15.6)') Z
    WRITE (xern6,'(1PE15.6)') P
    WRITE (xern7,'(1PE15.6)') uplim
    CALL XERMSG('SLATEC','RJ','MAX(X,Y,Z,P).GT.UPLIM WHERE X = '//xern3//&
      ' Y = '//xern4//' Z = '//xern5//' P = '//xern6//&
      ' AND UPLIM = '//xern7,3,1)
    RETURN
  END IF
  !
  IF ( MIN(X+Y,X+Z,Y+Z,P)<lolim ) THEN
    Ier = 2
    WRITE (xern3,'(1PE15.6)') X
    WRITE (xern4,'(1PE15.6)') Y
    WRITE (xern5,'(1PE15.6)') Z
    WRITE (xern6,'(1PE15.6)') P
    WRITE (xern7,'(1PE15.6)') lolim
    CALL XERMSG('SLATEC','RJ','MIN(X+Y,X+Z,Y+Z,P).LT.LOLIM WHERE X = '//&
      xern3//' Y = '//xern4//' Z = '//xern5//' P = '//xern6//' AND LOLIM = ',2,1)
    RETURN
  END IF
  !
  Ier = 0
  xn = X
  yn = Y
  zn = Z
  pn = P
  sigma = 0.0E0
  power4 = 1.0E0
  DO
    !
    mu = (xn+yn+zn+pn+pn)*0.20E0
    xndev = (mu-xn)/mu
    yndev = (mu-yn)/mu
    zndev = (mu-zn)/mu
    pndev = (mu-pn)/mu
    epslon = MAX(ABS(xndev),ABS(yndev),ABS(zndev),ABS(pndev))
    IF ( epslon<errtol ) THEN
      !
      ea = xndev*(yndev+zndev) + yndev*zndev
      eb = xndev*yndev*zndev
      ec = pndev*pndev
      e2 = ea - 3.0E0*ec
      e3 = eb + 2.0E0*pndev*(ea-ec)
      s1 = 1.0E0 + e2*(-c1+0.750E0*c3*e2-1.50E0*c4*e3)
      s2 = eb*(0.50E0*c2+pndev*(-c3-c3+pndev*c4))
      s3 = pndev*ea*(c2-pndev*c3) - c2*pndev*ec
      RJ = 3.0E0*sigma + power4*(s1+s2+s3)/(mu*SQRT(mu))
      EXIT
    ELSE
      xnroot = SQRT(xn)
      ynroot = SQRT(yn)
      znroot = SQRT(zn)
      lamda = xnroot*(ynroot+znroot) + ynroot*znroot
      alfa = pn*(xnroot+ynroot+znroot) + xnroot*ynroot*znroot
      alfa = alfa*alfa
      beta = pn*(pn+lamda)*(pn+lamda)
      sigma = sigma + power4*RC(alfa,beta,Ier)
      power4 = power4*0.250E0
      xn = (xn+lamda)*0.250E0
      yn = (yn+lamda)*0.250E0
      zn = (zn+lamda)*0.250E0
      pn = (pn+lamda)*0.250E0
    END IF
  END DO
END FUNCTION RJ
