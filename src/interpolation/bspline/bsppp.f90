!** BSPPP
SUBROUTINE BSPPP(T,A,N,K,Ldc,C,Xi,Lxi,Work)
  !>
  !***
  !  Convert the B-representation of a B-spline to the piecewise
  !            polynomial (PP) form.
  !***
  ! **Library:**   SLATEC
  !***
  ! **Category:**  E3, K6
  !***
  ! **Type:**      SINGLE PRECISION (BSPPP-S, DBSPPP-D)
  !***
  ! **Keywords:**  B-SPLINE, PIECEWISE POLYNOMIAL
  !***
  ! **Author:**  Amos, D. E., (SNLA)
  !***
  ! **Description:**
  !
  !     Written by Carl de Boor and modified by D. E. Amos
  !
  !     Abstract
  !         BSPPP is the BSPLPP routine of the reference.
  !
  !         BSPPP converts the B-representation (T,A,N,K) to the
  !         piecewise polynomial (PP) form (C,XI,LXI,K) for use with
  !         PPVAL.  Here XI(*), the break point array of length LXI, is
  !         the knot array T(*) with multiplicities removed.  The columns
  !         of the matrix C(I,J) contain the right Taylor derivatives
  !         for the polynomial expansion about XI(J) for the intervals
  !         XI(J) .LE. X .LE. XI(J+1), I=1,K, J=1,LXI.  Function PPVAL
  !         makes this evaluation at a specified point X in
  !         XI(1) .LE. X .LE. XI(LXI(1) .LE. X .LE. XI+1)
  !
  !     Description of Arguments
  !         Input
  !          T       - knot vector of length N+K
  !          A       - B-spline coefficient vector of length N
  !          N       - number of B-spline coefficients
  !                    N = sum of knot multiplicities-K
  !          K       - order of the B-spline, K .GE. 1
  !          LDC     - leading dimension of C, LDC .GE. K
  !
  !         Output
  !          C       - matrix of dimension at least (K,LXI) containing
  !                    right derivatives at break points
  !          XI      - XI break point vector of length LXI+1
  !          LXI     - number of break points, LXI .LE. N-K+1
  !          WORK    - work vector of length K*(N+3)
  !
  !     Error Conditions
  !         Improper input is a fatal error
  !
  !***
  ! **References:**  Carl de Boor, Package for calculating with B-splines,
  !                 SIAM Journal on Numerical Analysis 14, 3 (June 1977),
  !                 pp. 441-472.
  !***
  ! **Routines called:**  BSPDR, BSPEV, XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   800901  DATE WRITTEN
  !   890831  Modified array declarations.  (WRB)
  !   890831  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900326  Removed duplicate information from DESCRIPTION section.
  !           (WRB)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  USE service, ONLY : XERMSG
  !
  INTEGER ileft, inev, K, Ldc, Lxi, N, nk
  REAL A(*), C(Ldc,*), T(*), Work(*), Xi(*)
  !     DIMENSION T(N+K),XI(LXI+1),C(LDC,*)
  !     HERE, * = THE FINAL VALUE OF THE OUTPUT PARAMETER LXI.
  !* FIRST EXECUTABLE STATEMENT  BSPPP
  IF ( K<1 ) THEN
    CALL XERMSG('SLATEC','BSPPP','K DOES NOT SATISFY K.GE.1',2,1)
    RETURN
  ELSEIF ( N<K ) THEN
    CALL XERMSG('SLATEC','BSPPP','N DOES NOT SATISFY N.GE.K',2,1)
    RETURN
  ELSEIF ( Ldc<K ) THEN
    CALL XERMSG('SLATEC','BSPPP','LDC DOES NOT SATISFY LDC.GE.K',2,1)
    RETURN
  END IF
  CALL BSPDR(T,A,N,K,K,Work)
  Lxi = 0
  Xi(1) = T(K)
  inev = 1
  nk = N*K + 1
  DO ileft = K, N
    IF ( T(ileft+1)/=T(ileft) ) THEN
      Lxi = Lxi + 1
      Xi(Lxi+1) = T(ileft+1)
      CALL BSPEV(T,Work(1),N,K,K,Xi(Lxi),inev,C(1,Lxi),Work(nk))
    END IF
  END DO
  RETURN
END SUBROUTINE BSPPP
