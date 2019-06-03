!** DPFQAD
SUBROUTINE DPFQAD(F,Ldc,C,Xi,Lxi,K,Id,X1,X2,Tol,Quad,Ierr)
  !>
  !  Compute the integral on (X1,X2) of a product of a
  !            function F and the ID-th derivative of a B-spline,
  !            (PP-representation).
  !***
  ! **Library:**   SLATEC
  !***
  ! **Category:**  H2A2A1, E3, K6
  !***
  ! **Type:**      DOUBLE PRECISION (PFQAD-S, DPFQAD-D)
  !***
  ! **Keywords:**  B-SPLINE, DATA FITTING, INTERPOLATION, QUADRATURE, SPLINES
  !***
  ! **Author:**  Amos, D. E., (SNLA)
  !***
  ! **Description:**
  !
  !     Abstract    **** a double precision routine ****
  !         DPFQAD computes the integral on (X1,X2) of a product of a
  !         function F and the ID-th derivative of a B-spline, using the
  !         PP-representation (C,XI,LXI,K).  (X1,X2) is normally a sub-
  !         interval of XI(1) .LE. X .LE. XI(LXI+1).  An integration
  !         routine, DPPGQ8 (a modification of GAUS8), integrates the
  !         product on subintervals of (X1,X2) formed by the included
  !         break points.  Integration outside of (XI(1),XI(LXI+1)) is
  !         permitted provided F is defined.
  !
  !         The maximum number of significant digits obtainable in
  !         DBSQAD is the smaller of 18 and the number of digits
  !         carried in double precision arithmetic.
  !
  !     Description of arguments
  !         Input      F,C,XI,X1,X2,TOL are double precision
  !           F      - external function of one argument for the
  !                    integrand PF(X)=F(X)*DPPVAL(LDC,C,XI,LXI,K,ID,X,
  !                    INPPV)
  !           LDC    - leading dimension of matrix C, LDC .GE. K
  !           C(I,J) - right Taylor derivatives at XI(J), I=1,K, J=1,LXI
  !           XI(*)  - break point array of length LXI+1
  !           LXI    - number of polynomial pieces
  !           K      - order of B-spline, K .GE. 1
  !           ID     - order of the spline derivative, 0 .LE. ID .LE. K-1
  !                    ID=0 gives the spline function
  !           X1,X2  - end points of quadrature interval, normally in
  !                    XI(1) .LE. X .LE. XI(LXI+1)
  !           TOL    - desired accuracy for the quadrature, suggest
  !                    10.*DTOL .LT. TOL .LE. 0.1 where DTOL is the
  !                    maximum of 1.0D-18 and double precision unit
  !                    roundoff for the machine = D1MACH(4)
  !
  !         Output     QUAD is double precision
  !           QUAD   - integral of PF(X) on (X1,X2)
  !           IERR   - a status code
  !                    IERR=1 normal return
  !                         2 some quadrature does not meet the
  !                           requested tolerance
  !
  !     Error Conditions
  !         Improper input is a fatal error.
  !         Some quadrature does not meet the requested tolerance.
  !
  !***
  ! **References:**  D. E. Amos, Quadrature subroutines for splines and
  !                 B-splines, Report SAND79-1825, Sandia Laboratories,
  !                 December 1979.
  !***
  ! **Routines called:**  D1MACH, DINTRV, DPPGQ8, XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   800901  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900326  Removed duplicate information from DESCRIPTION section.
  !           (WRB)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  USE service, ONLY : XERMSG, D1MACH
  USE interpolation, ONLY : DPPGQ8, DINTRV
  !
  INTERFACE
    REAL(DP) FUNCTION F(X)
      IMPORT DP
      REAL(DP), INTENT(IN) :: X
    END FUNCTION F
  END INTERFACE
  INTEGER :: Id, Ierr, K, Ldc, Lxi
  REAL(DP) :: Quad, Tol, X1, X2, C(Ldc,Lxi), Xi(Lxi+1)
  INTEGER :: iflg, ilo, il1, il2, inppv, left, mf1, mf2
  REAL(DP) :: a, aa, ans, b, bb, q, ta, tb, wtol
  !
  !* FIRST EXECUTABLE STATEMENT  DPFQAD
  Ierr = 1
  Quad = 0.0D0
  IF ( K<1 ) THEN
    CALL XERMSG('DPFQAD','K DOES NOT SATISFY K.GE.1',2,1)
    RETURN
  ELSEIF ( Ldc<K ) THEN
    CALL XERMSG('DPFQAD','LDC DOES NOT SATISFY LDC.GE.K',2,1)
    RETURN
  ELSEIF ( Id<0.OR.Id>=K ) THEN
    CALL XERMSG('DPFQAD','ID DOES NOT SATISFY 0.LE.ID.LT.K',2,1)
    RETURN
  ELSEIF ( Lxi<1 ) THEN
    CALL XERMSG('DPFQAD','LXI DOES NOT SATISFY LXI.GE.1',2,1)
    RETURN
  ELSE
    wtol = D1MACH(4)
    wtol = MAX(wtol,1.0D-18)
    IF ( Tol>=wtol.AND.Tol<=0.1D0 ) THEN
      aa = MIN(X1,X2)
      bb = MAX(X1,X2)
      IF ( aa==bb ) RETURN
      ilo = 1
      CALL DINTRV(Xi,Lxi,aa,ilo,il1,mf1)
      CALL DINTRV(Xi,Lxi,bb,ilo,il2,mf2)
      q = 0.0D0
      inppv = 1
      DO left = il1, il2
        ta = Xi(left)
        a = MAX(aa,ta)
        IF ( left==1 ) a = aa
        tb = bb
        IF ( left<Lxi ) tb = Xi(left+1)
        b = MIN(bb,tb)
        CALL DPPGQ8(F,Ldc,C,Xi,Lxi,K,Id,a,b,inppv,Tol,ans,iflg)
        IF ( iflg>1 ) Ierr = 2
        q = q + ans
      END DO
      IF ( X1>X2 ) q = -q
      Quad = q
      RETURN
    END IF
  END IF
  !
  CALL XERMSG('DPFQAD','TOL IS LESS DTOL OR GREATER THAN 0.1',2,1)
  RETURN
END SUBROUTINE DPFQAD
