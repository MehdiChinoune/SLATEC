!** PPQAD
SUBROUTINE PPQAD(Ldc,C,Xi,Lxi,K,X1,X2,Pquad)
  !> Compute the integral on (X1,X2) of a K-th order B-spline
  !            using the piecewise polynomial (PP) representation.
  !***
  ! **Library:**   SLATEC
  !***
  ! **Category:**  H2A2A1, E3, K6
  !***
  ! **Type:**      SINGLE PRECISION (PPQAD-S, DPPQAD-D)
  !***
  ! **Keywords:**  B-SPLINE, DATA FITTING, INTERPOLATION, QUADRATURE, SPLINES
  !***
  ! **Author:**  Amos, D. E., (SNLA)
  !***
  ! **Description:**
  !
  !     Abstract
  !         PPQAD computes the integral on (X1,X2) of a K-th order
  !         B-spline using the piecewise polynomial representation
  !         (C,XI,LXI,K).  Here the Taylor expansion about the left
  !         end point XI(J) of the J-th interval is integrated and
  !         evaluated on subintervals of (X1,X2) which are formed by
  !         included break points.  Integration outside (XI(1),XI(LXI+1))
  !         is permitted.
  !
  !     Description of Arguments
  !         Input
  !           LDC    - leading dimension of matrix C, LDC >= K
  !           C(I,J) - right Taylor derivatives at XI(J), I=1,K, J=1,LXI
  !           XI(*)  - break point array of length LXI+1
  !           LXI    - number of polynomial pieces
  !           K      - order of B-spline, K >= 1
  !           X1,X2  - end points of quadrature interval, normally in
  !                    XI(1) <= X <= XI(LXI+1)
  !
  !         Output
  !           PQUAD  - integral of the PP representation over (X1,X2)
  !
  !     Error Conditions
  !         Improper input is a fatal error
  !
  !***
  ! **References:**  D. E. Amos, Quadrature subroutines for splines and
  !                 B-splines, Report SAND79-1825, Sandia Laboratories,
  !                 December 1979.
  !***
  ! **Routines called:**  INTRV, XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   800901  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  USE service, ONLY : XERMSG
  !
  INTEGER :: K, Ldc, Lxi
  REAL(SP) :: C(Ldc,Lxi), Pquad, Xi(Lxi+1), X1, X2
  INTEGER :: i, ii, il, ilo, il1, il2, im, left, mf1, mf2
  REAL(SP) :: a, aa, bb, dx, flk, q, s, ss(2), ta, tb, x
  !
  !* FIRST EXECUTABLE STATEMENT  PPQAD
  Pquad = 0._SP
  IF( K<1 ) THEN
    !
    !
    CALL XERMSG('PPQAD','K DOES NOT SATISFY K>=1',2,1)
    RETURN
  ELSEIF( Lxi<1 ) THEN
    CALL XERMSG('PPQAD','LXI DOES NOT SATISFY LXI>=1',2,1)
    RETURN
  ELSEIF( Ldc<K ) THEN
    CALL XERMSG('PPQAD','LDC DOES NOT SATISFY LDC>=K',2,1)
    RETURN
  END IF
  aa = MIN(X1,X2)
  bb = MAX(X1,X2)
  IF( aa==bb ) RETURN
  ilo = 1
  CALL INTRV(Xi,Lxi,aa,ilo,il1,mf1)
  CALL INTRV(Xi,Lxi,bb,ilo,il2,mf2)
  q = 0._SP
  DO left = il1, il2
    ta = Xi(left)
    a = MAX(aa,ta)
    IF( left==1 ) a = aa
    tb = bb
    IF( left<Lxi ) tb = Xi(left+1)
    x = MIN(bb,tb)
    DO ii = 1, 2
      ss(ii) = 0._SP
      dx = x - Xi(left)
      IF( dx/=0._SP ) THEN
        s = C(K,left)
        flk = K
        im = K - 1
        il = im
        DO i = 1, il
          s = s*dx/flk + C(im,left)
          im = im - 1
          flk = flk - 1._SP
        END DO
        ss(ii) = s*dx
      END IF
      x = a
    END DO
    q = q + (ss(1)-ss(2))
  END DO
  IF( X1>X2 ) q = -q
  Pquad = q
  RETURN
END SUBROUTINE PPQAD
