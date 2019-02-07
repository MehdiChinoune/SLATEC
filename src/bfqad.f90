!*==BFQAD.f90  processed by SPAG 6.72Dc at 11:01 on  6 Feb 2019
!DECK BFQAD
SUBROUTINE BFQAD(F,T,Bcoef,N,K,Id,X1,X2,Tol,Quad,Ierr,Work)
  IMPLICIT NONE
  !*--BFQAD5
  !*** Start of declarations inserted by SPAG
  INTEGER inbv
  !*** End of declarations inserted by SPAG
  !***BEGIN PROLOGUE  BFQAD
  !***PURPOSE  Compute the integral of a product of a function and a
  !            derivative of a B-spline.
  !***LIBRARY   SLATEC
  !***CATEGORY  H2A2A1, E3, K6
  !***TYPE      SINGLE PRECISION (BFQAD-S, DBFQAD-D)
  !***KEYWORDS  INTEGRAL OF B-SPLINE, QUADRATURE
  !***AUTHOR  Amos, D. E., (SNLA)
  !***DESCRIPTION
  !
  !     Abstract
  !         BFQAD computes the integral on (X1,X2) of a product of a
  !         function F and the ID-th derivative of a K-th order B-spline,
  !         using the B-representation (T,BCOEF,N,K).  (X1,X2) must be
  !         a subinterval of T(K) .LE. X .le. T(N+1).  An integration
  !         routine BSGQ8 (a modification
  !         of GAUS8), integrates the product on sub-
  !         intervals of (X1,X2) formed by included (distinct) knots.
  !
  !     Description of Arguments
  !         Input
  !           F      - external function of one argument for the
  !                    integrand BF(X)=F(X)*BVALU(T,BCOEF,N,K,ID,X,INBV,
  !                    WORK)
  !           T      - knot array of length N+K
  !           BCOEF  - coefficient array of length N
  !           N      - length of coefficient array
  !           K      - order of B-spline, K .GE. 1
  !           ID     - order of the spline derivative, 0 .LE. ID .LE. K-1
  !                    ID=0 gives the spline function
  !           X1,X2  - end points of quadrature interval in
  !                    T(K) .LE. X .LE. T(N+1)
  !           TOL    - desired accuracy for the quadrature, suggest
  !                    10.*STOL .LT. TOL .LE. 0.1 where STOL is the single
  !                    precision unit roundoff for the machine = R1MACH(4)
  !
  !         Output
  !           QUAD   - integral of BF(X) on (X1,X2)
  !           IERR   - a status code
  !                    IERR=1  normal return
  !                         2  some quadrature on (X1,X2) does not meet
  !                            the requested tolerance.
  !           WORK   - work vector of length 3*K
  !
  !     Error Conditions
  !         X1 or X2 not in T(K) .LE. X .LE. T(N+1) is a fatal error.
  !         TOL not greater than the single precision unit roundoff or
  !         less than 0.1 is a fatal error.
  !         Some quadrature fails to meet the requested tolerance.
  !
  !***REFERENCES  D. E. Amos, Quadrature subroutines for splines and
  !                 B-splines, Report SAND79-1825, Sandia Laboratories,
  !                 December 1979.
  !***ROUTINES CALLED  BSGQ8, INTRV, R1MACH, XERMSG
  !***REVISION HISTORY  (YYMMDD)
  !   800901  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900326  Removed duplicate information from DESCRIPTION section.
  !           (WRB)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  !***END PROLOGUE  BFQAD
  !
  !
  INTEGER Id , Ierr , iflg , ilo , il1 , il2 , K , left , mflag , N , npk , &
    np1
  REAL a , aa , ans , b , bb , Bcoef , q , Quad , T , ta , tb , Tol , Work , &
    wtol , X1 , X2
  REAL R1MACH , F
  DIMENSION T(*) , Bcoef(*) , Work(*)
  EXTERNAL F
  !***FIRST EXECUTABLE STATEMENT  BFQAD
  Ierr = 1
  Quad = 0.0E0
  IF ( K<1 ) THEN
    CALL XERMSG('SLATEC','BFQAD','K DOES NOT SATISFY K.GE.1',2,1)
    RETURN
  ELSEIF ( N<K ) THEN
    CALL XERMSG('SLATEC','BFQAD','N DOES NOT SATISFY N.GE.K',2,1)
    RETURN
  ELSEIF ( Id<0.OR.Id>=K ) THEN
    CALL XERMSG('SLATEC','BFQAD','ID DOES NOT SATISFY 0 .LE. ID .LT. K',2,1)
    GOTO 99999
  ELSE
    wtol = R1MACH(4)
    IF ( Tol>=wtol.AND.Tol<=0.1E0 ) THEN
      aa = MIN(X1,X2)
      bb = MAX(X1,X2)
      IF ( aa>=T(K) ) THEN
        np1 = N + 1
        IF ( bb<=T(np1) ) THEN
          IF ( aa==bb ) RETURN
          npk = N + K
          !
          ilo = 1
          CALL INTRV(T,npk,aa,ilo,il1,mflag)
          CALL INTRV(T,npk,bb,ilo,il2,mflag)
          IF ( il2>=np1 ) il2 = N
          inbv = 1
          q = 0.0E0
          DO left = il1 , il2
            ta = T(left)
            tb = T(left+1)
            IF ( ta/=tb ) THEN
              a = MAX(aa,ta)
              b = MIN(bb,tb)
              CALL BSGQ8(F,T,Bcoef,N,K,Id,a,b,inbv,Tol,ans,iflg,Work)
              IF ( iflg>1 ) Ierr = 2
              q = q + ans
            ENDIF
          ENDDO
          IF ( X1>X2 ) q = -q
          Quad = q
          RETURN
        ENDIF
      ENDIF
      !
      !
      CALL XERMSG('SLATEC','BFQAD',&
        'X1 OR X2 OR BOTH DO NOT SATISFY T(K).LE.X.LE.T(N+1)',2,1)
      RETURN
    ENDIF
  ENDIF
  CALL XERMSG('SLATEC','BFQAD',&
    'TOL IS LESS THAN THE SINGLE PRECISION TOLERANCE OR '//&
    'GREATER THAN 0.1',2,1)
  RETURN
  99999 END SUBROUTINE BFQAD
