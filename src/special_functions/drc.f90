!** DRC
REAL(DP) FUNCTION DRC(X,Y,Ier)
  !> Calculate a double precision approximation to
  !             DRC(X,Y) = Integral from zero to infinity of
  !                              -1/2     -1
  !                    (1/2)(t+X)    (t+Y)  dt,
  !            where X is nonnegative and Y is positive.
  !***
  ! **Library:**   SLATEC
  !***
  ! **Category:**  C14
  !***
  ! **Type:**      DOUBLE PRECISION (RC-S, DRC-D)
  !***
  ! **Keywords:**  DUPLICATION THEOREM, ELEMENTARY FUNCTIONS,
  !             ELLIPTIC INTEGRAL, TAYLOR SERIES
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
  !   1.     DRC
  !          Standard FORTRAN function routine
  !          Double precision version
  !          The routine calculates an approximation result to
  !          DRC(X,Y) = integral from zero to infinity of
  !
  !                              -1/2     -1
  !                    (1/2)(t+X)    (t+Y)  dt,
  !
  !          where X is nonnegative and Y is positive.  The duplication
  !          theorem is iterated until the variables are nearly equal,
  !          and the function is then expanded in Taylor series to fifth
  !          order.  Logarithmic, inverse circular, and inverse hyper-
  !          bolic functions can be expressed in terms of DRC.
  !
  !   2.     Calling Sequence
  !          DRC( X, Y, IER )
  !
  !          Parameters On Entry
  !          Values assigned by the calling routine
  !
  !          X      - Double precision, nonnegative variable
  !
  !          Y      - Double precision, positive variable
  !
  !
  !
  !          On Return  (values assigned by the DRC routine)
  !
  !          DRC    - Double precision approximation to the integral
  !
  !          IER    - Integer to indicate normal or abnormal termination.
  !
  !                     IER = 0 Normal and reliable termination of the
  !                             routine.  It is assumed that the requested
  !                             accuracy has been achieved.
  !
  !                     IER > 0 Abnormal termination of the routine
  !
  !          X and Y are unaltered.
  !
  !   3.    Error messages
  !
  !         Value of IER assigned by the DRC routine
  !
  !                  Value assigned         Error message printed
  !                  IER = 1                X<0.0D0 .OR. Y<=0.0D0
  !                      = 2                X+Y<LOLIM
  !                      = 3                MAX(X,Y) > UPLIM
  !
  !   4.     Control parameters
  !
  !                  Values of LOLIM, UPLIM, and ERRTOL are set by the
  !                  routine.
  !
  !          LOLIM and UPLIM determine the valid range of X and Y
  !
  !          LOLIM  - Lower limit of valid arguments
  !
  !                   Not less  than 5 * (machine minimum)  .
  !
  !          UPLIM  - Upper limit of valid arguments
  !
  !                   Not greater than (machine maximum) / 5 .
  !
  !
  !                     Acceptable values for:   LOLIM       UPLIM
  !                     IBM 360/370 SERIES   :   3.0D-78     1.0D+75
  !                     CDC 6000/7000 SERIES :   1.0D-292    1.0D+321
  !                     UNIVAC 1100 SERIES   :   1.0D-307    1.0D+307
  !                     CRAY                 :   2.3D-2466   1.0D+2465
  !                     VAX 11 SERIES        :   1.5D-38     3.0D+37
  !
  !          ERRTOL determines the accuracy of the answer
  !
  !                 The value assigned by the routine will result
  !                 in solution precision within 1-2 decimals of
  !                 "machine precision".
  !
  !
  !          ERRTOL  - relative error due to truncation is less than
  !                    16 * ERRTOL ** 6 / (1 - 2 * ERRTOL).
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
  !          Sample choices:  ERRTOL   Relative truncation
  !                                    error less than
  !                           1.0D-3    2.0D-17
  !                           3.0D-3    2.0D-14
  !                           1.0D-2    2.0D-11
  !                           3.0D-2    2.0D-8
  !                           1.0D-1    2.0D-5
  !
  !
  !                    Decreasing ERRTOL by a factor of 10 yields six more
  !                    decimal digits of accuracy at the expense of one or
  !                    two more iterations of the duplication theorem.
  !
  !- Long Description:
  !
  !   DRC special comments
  !
  !
  !
  !
  !                  Check: DRC(X,X+Z) + DRC(Y,Y+Z) = DRC(0,Z)
  !
  !                  where X, Y, and Z are positive and X * Y = Z * Z
  !
  !
  !          On Input:
  !
  !          X, and Y are the variables in the integral DRC(X,Y).
  !
  !          On Output:
  !
  !          X and Y are unaltered.
  !
  !
  !
  !                    DRC(0,1/4)=DRC(1/16,1/8)=PI=3.14159...
  !
  !                    DRC(9/4,2)=LN(2)
  !
  !
  !
  !          ********************************************************
  !
  !          WARNING: Changes in the program may improve speed at the
  !                   expense of robustness.
  !
  !
  !   --------------------------------------------------------------------
  !
  !   Special functions via DRC
  !
  !
  !
  !                  LN X                X > 0
  !
  !                                             2
  !                  LN(X) = (X-1) DRC(((1+X)/2) , X )
  !
  !
  !   --------------------------------------------------------------------
  !
  !                  ARCSIN X            -1 <= X <= 1
  !
  !                                       2
  !                  ARCSIN X = X DRC (1-X  ,1 )
  !
  !   --------------------------------------------------------------------
  !
  !                  ARCCOS X            0 <= X <= 1
  !
  !
  !                                     2       2
  !                  ARCCOS X = SQRT(1-X ) DRC(X  ,1 )
  !
  !   --------------------------------------------------------------------
  !
  !                  ARCTAN X            -INF < X < +INF
  !
  !                                        2
  !                  ARCTAN X = X DRC(1,1+X  )
  !
  !   --------------------------------------------------------------------
  !
  !                  ARCCOT X            0 <= X < INF
  !
  !                                  2   2
  !                  ARCCOT X = DRC(X  ,X +1 )
  !
  !   --------------------------------------------------------------------
  !
  !                  ARCSINH X           -INF < X < +INF
  !
  !                                       2
  !                  ARCSINH X = X DRC(1+X  ,1 )
  !
  !   --------------------------------------------------------------------
  !
  !                  ARCCOSH X           X >= 1
  !
  !                                    2         2
  !                  ARCCOSH X = SQRT(X -1) DRC(X  ,1 )
  !
  !   --------------------------------------------------------------------
  !
  !                  ARCTANH X           -1 < X < 1
  !
  !                                         2
  !                  ARCTANH X = X DRC(1,1-X  )
  !
  !   --------------------------------------------------------------------
  !
  !                  ARCCOTH X           X > 1
  !
  !                                   2   2
  !                  ARCCOTH X = DRC(X  ,X -1 )
  !
  !   --------------------------------------------------------------------
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
  ! **Routines called:**  D1MACH, XERMSG

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
  !           editorial changes.  (RWC))
  !   920501  Reformatted the REFERENCES section.  (WRB)
  USE service, ONLY : XERMSG, D1MACH
  INTEGER :: Ier
  REAL(DP) :: X, Y
  REAL(DP) :: mu, s, sn, xn, yn, lamda
  CHARACTER(16) :: xern3, xern4, xern5
  REAL(DP), PARAMETER :: errtol = (D1MACH(3)/16._DP)**(1._DP/6._DP), &
    lolim = 5._DP*D1MACH(1), uplim = D1MACH(2)/5._DP
  REAL(DP), PARAMETER :: c1 = 1._DP/7._DP, c2 = 9._DP/22._DP
  !* FIRST EXECUTABLE STATEMENT  DRC
  !         CALL ERROR HANDLER IF NECESSARY.
  !
  DRC = 0._DP
  IF( X<0._DP .OR. Y<=0._DP ) THEN
    Ier = 1
    WRITE (xern3,'(1PE15.6)') X
    WRITE (xern4,'(1PE15.6)') Y
    CALL XERMSG('DRC','X<0 .OR. Y<=0 WHERE X = '//xern3//&
      ' AND Y = '//xern4,1,1)
    RETURN
  END IF
  !
  IF( MAX(X,Y)>uplim ) THEN
    Ier = 3
    WRITE (xern3,'(1PE15.6)') X
    WRITE (xern4,'(1PE15.6)') Y
    WRITE (xern5,'(1PE15.6)') uplim
    CALL XERMSG('DRC','MAX(X,Y)>UPLIM WHERE X = '//xern3//&
      ' Y = '//xern4//' AND UPLIM = '//xern5,3,1)
    RETURN
  END IF
  !
  IF( X+Y<lolim ) THEN
    Ier = 2
    WRITE (xern3,'(1PE15.6)') X
    WRITE (xern4,'(1PE15.6)') Y
    WRITE (xern5,'(1PE15.6)') lolim
    CALL XERMSG('DRC','X+Y<LOLIM WHERE X = '//xern3//' Y = '//&
      xern4//' AND LOLIM = '//xern5,2,1)
    RETURN
  END IF
  !
  Ier = 0
  xn = X
  yn = Y
  DO
    !
    mu = (xn+yn+yn)/3._DP
    sn = (yn+mu)/mu - 2._DP
    IF( ABS(sn)<errtol ) THEN
      !
      s = sn*sn*(0.30_DP+sn*(c1+sn*(0.3750_DP+sn*c2)))
      DRC = (1._DP+s)/SQRT(mu)
      EXIT
    ELSE
      lamda = 2._DP*SQRT(xn)*SQRT(yn) + yn
      xn = (xn+lamda)*0.250_DP
      yn = (yn+lamda)*0.250_DP
    END IF
  END DO
END FUNCTION DRC
