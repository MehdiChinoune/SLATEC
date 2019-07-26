!** PFQAD
PURE SUBROUTINE PFQAD(F,Ldc,C,Xi,Lxi,K,Id,X1,X2,Tol,Quad,Ierr)
  !> Compute the integral on (X1,X2) of a product of a function F and the
  !  ID-th derivative of a B-spline, (PP-representation).
  !***
  ! **Library:**   SLATEC
  !***
  ! **Category:**  H2A2A1, E3, K6
  !***
  ! **Type:**      SINGLE PRECISION (PFQAD-S, DPFQAD-D)
  !***
  ! **Keywords:**  B-SPLINE, DATA FITTING, INTERPOLATION, QUADRATURE, SPLINES
  !***
  ! **Author:**  Amos, D. E., (SNLA)
  !***
  ! **Description:**
  !
  !     Abstract
  !         PFQAD computes the integral on (X1,X2) of a product of a
  !         function F and the ID-th derivative of a B-spline, using the
  !         PP-representation (C,XI,LXI,K). (X1,X2) is normally a sub-
  !         interval of XI(1) <= X <= XI(LXI+1).  An integration rou-
  !         tine, PPGQ8(a modification of GAUS8), integrates the product
  !         on sub-intervals of (X1,X2) formed by the included break
  !         points.  Integration outside of (XI(1),XI(LXI+1)) is permitted
  !         provided F is defined.
  !
  !     Description of Arguments
  !         Input
  !           F      - external function of one argument for the
  !                    integrand PF(X)=F(X)*PPVAL(LDC,C,XI,LXI,K,ID,X,
  !                    INPPV)
  !           LDC    - leading dimension of matrix C, LDC >= K
  !           C(I,J) - right Taylor derivatives at XI(J), I=1,K, J=1,LXI
  !           XI(*)  - break point array of length LXI+1
  !           LXI    - number of polynomial pieces
  !           K      - order of B-spline, K >= 1
  !           ID     - order of the spline derivative, 0 <= ID <= K-1
  !                    ID=0 gives the spline function
  !           X1,X2  - end points of quadrature interval, normally in
  !                    XI(1) <= X <= XI(LXI+1)
  !           TOL    - desired accuracy for the quadrature, suggest
  !                    10.*STOL < TOL <= 0.1 where STOL is the single
  !                    precision unit roundoff for the machine = eps_sp
  !
  !         Output
  !           QUAD   - integral of PF(X) on (X1,X2)
  !           IERR   - a status code
  !                    IERR=1 normal return
  !                         2 some quadrature does not meet the
  !                           requested tolerance
  !
  !     Error Conditions
  !         TOL not greater than the single precision unit roundoff or
  !         less than 0.1 is a fatal error.
  !         Some quadrature does not meet the requested tolerance.
  !
  !***
  ! **References:**  D. E. Amos, Quadrature subroutines for splines and
  !                 B-splines, Report SAND79-1825, Sandia Laboratories,
  !                 December 1979.
  !***
  ! **Routines called:**  INTRV, PPGQ8, R1MACH, XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   800901  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900326  Removed duplicate information from DESCRIPTIONsection.  (WRB)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  USE service, ONLY : eps_sp
  USE interpolation, ONLY : INTRV
  !
  INTERFACE
    PURE REAL(SP) FUNCTION F(X)
      IMPORT SP
      REAL(SP), INTENT(IN) :: X
    END FUNCTION F
  END INTERFACE
  INTEGER, INTENT(IN) :: Id, K, Ldc, Lxi
  INTEGER, INTENT(OUT) :: Ierr
  REAL(SP), INTENT(IN) :: C(Ldc,Lxi), Xi(Lxi+1), X1, X2
  REAL(SP), INTENT(INOUT) :: Tol
  REAL(SP), INTENT(OUT) :: Quad
  !
  INTEGER :: iflg, ilo, il1, il2, inppv, left, mf1, mf2
  REAL(SP) :: a, aa, ans, b, bb, q, ta, tb, wtol
  !
  !* FIRST EXECUTABLE STATEMENT  PFQAD
  Ierr = 1
  Quad = 0._SP
  IF( K<1 ) THEN
    ERROR STOP 'PFQAD : K DOES NOT SATISFY K>=1'
    RETURN
  ELSEIF( Ldc<K ) THEN
    ERROR STOP 'PFQAD : LDC DOES NOT SATISFY LDC>=K'
    RETURN
  ELSEIF( Id<0 .OR. Id>=K ) THEN
    ERROR STOP 'PFQAD : ID DOES NOT SATISFY 0<=ID<K'
    RETURN
  ELSEIF( Lxi<1 ) THEN
    ERROR STOP 'PFQAD : LXI DOES NOT SATISFY LXI>=1'
    RETURN
  ELSE
    wtol = eps_sp
    IF( Tol>=wtol .AND. Tol<=0.1_SP ) THEN
      aa = MIN(X1,X2)
      bb = MAX(X1,X2)
      IF( aa==bb ) RETURN
      ilo = 1
      CALL INTRV(Xi,Lxi,aa,ilo,il1,mf1)
      CALL INTRV(Xi,Lxi,bb,ilo,il2,mf2)
      q = 0._SP
      inppv = 1
      DO left = il1, il2
        ta = Xi(left)
        a = MAX(aa,ta)
        IF( left==1 ) a = aa
        tb = bb
        IF( left<Lxi ) tb = Xi(left+1)
        b = MIN(bb,tb)
        CALL PPGQ8(F,Ldc,C,Xi,Lxi,K,Id,a,b,Tol,ans,iflg)
        IF( iflg>1 ) Ierr = 2
        q = q + ans
      END DO
      IF( X1>X2 ) q = -q
      Quad = q
      RETURN
    END IF
  END IF
  !
  ERROR STOP 'PFQAD : TOL IS LESS THAN THE SINGLE PRECISION TOLERANCE OR GREATER THAN 0.1'
  !
  RETURN
END SUBROUTINE PFQAD