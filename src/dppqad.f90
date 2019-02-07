!*==DPPQAD.f90  processed by SPAG 6.72Dc at 11:01 on  6 Feb 2019
!DECK DPPQAD
SUBROUTINE DPPQAD(Ldc,C,Xi,Lxi,K,X1,X2,Pquad)
  IMPLICIT NONE
  !*--DPPQAD5
  !***BEGIN PROLOGUE  DPPQAD
  !***PURPOSE  Compute the integral on (X1,X2) of a K-th order B-spline
  !            using the piecewise polynomial (PP) representation.
  !***LIBRARY   SLATEC
  !***CATEGORY  H2A2A1, E3, K6
  !***TYPE      DOUBLE PRECISION (PPQAD-S, DPPQAD-D)
  !***KEYWORDS  B-SPLINE, DATA FITTING, INTERPOLATION, QUADRATURE, SPLINES
  !***AUTHOR  Amos, D. E., (SNLA)
  !***DESCRIPTION
  !
  !     Abstract    **** a double precision routine ****
  !         DPPQAD computes the integral on (X1,X2) of a K-th order
  !         B-spline using the piecewise polynomial representation
  !         (C,XI,LXI,K).  Here the Taylor expansion about the left
  !         end point XI(J) of the J-th interval is integrated and
  !         evaluated on subintervals of (X1,X2) which are formed by
  !         included break points.  Integration outside (XI(1),XI(LXI+1))
  !         is permitted.
  !
  !     Description of Arguments
  !         Input      C,XI,X1,X2 are double precision
  !           LDC    - leading dimension of matrix C, LDC .GE. K
  !           C(I,J) - right Taylor derivatives at XI(J), I=1,K , J=1,LXI
  !           XI(*)  - break point array of length LXI+1
  !           LXI    - number of polynomial pieces
  !           K      - order of B-spline, K .GE. 1
  !           X1,X2  - end points of quadrature interval, normally in
  !                    XI(1) .LE. X .LE. XI(LXI+1)
  !
  !         Output     PQUAD is double precision
  !           PQUAD  - integral of the PP representation over (X1,X2)
  !
  !     Error Conditions
  !         Improper input is a fatal error
  !
  !***REFERENCES  D. E. Amos, Quadrature subroutines for splines and
  !                 B-splines, Report SAND79-1825, Sandia Laboratories,
  !                 December 1979.
  !***ROUTINES CALLED  DINTRV, XERMSG
  !***REVISION HISTORY  (YYMMDD)
  !   800901  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890911  Removed unnecessary intrinsics.  (WRB)
  !   890911  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  !***END PROLOGUE  DPPQAD
  !
  INTEGER i , ii , il , ilo , il1 , il2 , im , K , Ldc , left , Lxi , mf1 , &
    mf2
  DOUBLE PRECISION a , aa , bb , C , dx , flk , Pquad , q , s , ss , ta , &
    tb , x , Xi , X1 , X2
  DIMENSION Xi(*) , C(Ldc,*) , ss(2)
  !
  !***FIRST EXECUTABLE STATEMENT  DPPQAD
  Pquad = 0.0D0
  IF ( K<1 ) THEN
    !
    !
    CALL XERMSG('SLATEC','DPPQAD','K DOES NOT SATISFY K.GE.1',2,1)
    RETURN
  ELSEIF ( Lxi<1 ) THEN
    CALL XERMSG('SLATEC','DPPQAD','LXI DOES NOT SATISFY LXI.GE.1',2,1)
    RETURN
  ELSEIF ( Ldc<K ) THEN
    CALL XERMSG('SLATEC','DPPQAD','LDC DOES NOT SATISFY LDC.GE.K',2,1)
    GOTO 99999
  ENDIF
  aa = MIN(X1,X2)
  bb = MAX(X1,X2)
  IF ( aa==bb ) RETURN
  ilo = 1
  CALL DINTRV(Xi,Lxi,aa,ilo,il1,mf1)
  CALL DINTRV(Xi,Lxi,bb,ilo,il2,mf2)
  q = 0.0D0
  DO left = il1 , il2
    ta = Xi(left)
    a = MAX(aa,ta)
    IF ( left==1 ) a = aa
    tb = bb
    IF ( left<Lxi ) tb = Xi(left+1)
    x = MIN(bb,tb)
    DO ii = 1 , 2
      ss(ii) = 0.0D0
      dx = x - Xi(left)
      IF ( dx/=0.0D0 ) THEN
        s = C(K,left)
        flk = K
        im = K - 1
        il = im
        DO i = 1 , il
          s = s*dx/flk + C(im,left)
          im = im - 1
          flk = flk - 1.0D0
        ENDDO
        ss(ii) = s*dx
      ENDIF
      x = a
    ENDDO
    q = q + (ss(1)-ss(2))
  ENDDO
  IF ( X1>X2 ) q = -q
  Pquad = q
  RETURN
  99999 END SUBROUTINE DPPQAD
