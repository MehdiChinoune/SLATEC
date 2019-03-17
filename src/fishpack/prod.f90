!DECK PROD
SUBROUTINE PROD(Nd,Bd,Nm1,Bm1,Nm2,Bm2,Na,Aa,X,Y,M,A,B,C,D,W,U)
  IMPLICIT NONE
  !***BEGIN PROLOGUE  PROD
  !***SUBSIDIARY
  !***PURPOSE  Subsidiary to BLKTRI
  !***LIBRARY   SLATEC
  !***TYPE      SINGLE PRECISION (PROD-S, PROC-C)
  !***AUTHOR  (UNKNOWN)
  !***DESCRIPTION
  !
  ! PROD applies a sequence of matrix operations to the vector X and
  ! stores the result in Y.
  !
  ! BD,BM1,BM2 are arrays containing roots of certain B polynomials.
  ! ND,NM1,NM2 are the lengths of the arrays BD,BM1,BM2 respectively.
  ! AA         Array containing scalar multipliers of the vector X.
  ! NA         is the length of the array AA.
  ! X,Y        The matrix operations are applied to X and the result is Y.
  ! A,B,C      are arrays which contain the tridiagonal matrix.
  ! M          is the order of the matrix.
  ! D,W,U      are working arrays.
  ! IS         determines whether or not a change in sign is made.
  !
  !***SEE ALSO  BLKTRI
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   801001  DATE WRITTEN
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900402  Added TYPE section.  (WRB)
  !***END PROLOGUE  PROD
  REAL A, Aa, B, Bd, Bm1, Bm2, C, D, den, rt, U, W, X, Y
  INTEGER ia, ibr, id, j, k, M, m1, m2, mm, Na, Nd, Nm1, Nm2
  DIMENSION A(*), B(*), C(*), X(*), Y(*), D(*), W(*), Bd(*), Bm1(*)&
    , Bm2(*), Aa(*), U(*)
  !***FIRST EXECUTABLE STATEMENT  PROD
  DO j = 1, M
    W(j) = X(j)
    Y(j) = W(j)
  ENDDO
  mm = M - 1
  id = Nd
  ibr = 0
  m1 = Nm1
  m2 = Nm2
  ia = Na
  100 CONTINUE
  DO
    IF ( ia>0 ) THEN
      rt = Aa(ia)
      IF ( Nd==0 ) rt = -rt
      ia = ia - 1
      !
      ! SCALAR MULTIPLICATION
      !
      DO j = 1, M
        Y(j) = rt*W(j)
      ENDDO
    ENDIF
    IF ( id<=0 ) RETURN
    rt = Bd(id)
    id = id - 1
    IF ( id==0 ) ibr = 1
    !
    ! BEGIN SOLUTION TO SYSTEM
    !
    D(M) = A(M)/(B(M)-rt)
    W(M) = Y(M)/(B(M)-rt)
    DO j = 2, mm
      k = M - j
      den = B(k+1) - rt - C(k+1)*D(k+2)
      D(k+1) = A(k+1)/den
      W(k+1) = (Y(k+1)-C(k+1)*W(k+2))/den
    ENDDO
    den = B(1) - rt - C(1)*D(2)
    W(1) = 1.
    IF ( den/=0 ) W(1) = (Y(1)-C(1)*W(2))/den
    DO j = 2, M
      W(j) = W(j) - D(j)*W(j-1)
    ENDDO
    IF ( Na<=0 ) THEN
      IF ( m1<=0 ) THEN
        IF ( m2>0 ) GOTO 400
        EXIT
      ELSE
        IF ( m2<=0 ) GOTO 300
        IF ( ABS(Bm1(m1))>ABS(Bm2(m2)) ) GOTO 300
        GOTO 400
      ENDIF
    ENDIF
  ENDDO
  200 CONTINUE
  DO j = 1, M
    Y(j) = W(j)
  ENDDO
  ibr = 1
  GOTO 100
  300 CONTINUE
  IF ( ibr<=0 ) THEN
    IF ( ABS(Bm1(m1)-Bd(id))<ABS(Bm1(m1)-rt) ) GOTO 200
  ENDIF
  rt = rt - Bm1(m1)
  m1 = m1 - 1
  GOTO 500
  400 CONTINUE
  IF ( ibr<=0 ) THEN
    IF ( ABS(Bm2(m2)-Bd(id))<ABS(Bm2(m2)-rt) ) GOTO 200
  ENDIF
  rt = rt - Bm2(m2)
  m2 = m2 - 1
  500 CONTINUE
  DO j = 1, M
    Y(j) = Y(j) + rt*W(j)
  ENDDO
  GOTO 100
  RETURN
END SUBROUTINE PROD
