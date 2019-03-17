!DECK PROCP
SUBROUTINE PROCP(Nd,Bd,Nm1,Bm1,Nm2,Bm2,Na,Aa,X,Y,M,A,B,C,D,U,W)
  IMPLICIT NONE
  !***BEGIN PROLOGUE  PROCP
  !***SUBSIDIARY
  !***PURPOSE  Subsidiary to CBLKTR
  !***LIBRARY   SLATEC
  !***TYPE      COMPLEX (PRODP-C, PROCP-C)
  !***AUTHOR  (UNKNOWN)
  !***DESCRIPTION
  !
  ! PROCP applies a sequence of matrix operations to the vector X and
  ! stores the result in Y (periodic boundary conditions).
  !
  ! BD,BM1,BM2 are arrays containing roots of certain B polynomials.
  ! ND,NM1,NM2 are the lengths of the arrays BD,BM1,BM2 respectively.
  ! AA         Array containing scalar multipliers of the vector X.
  ! NA         is the length of the array AA.
  ! X,Y        The matrix operations are applied to X and the result is Y.
  ! A,B,C      are arrays which contain the tridiagonal matrix.
  ! M          is the order of the matrix.
  ! D,U,W      are working arrays.
  ! IS         determines whether or not a change in sign is made.
  !
  !***SEE ALSO  CBLKTR
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   801001  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900402  Added TYPE section.  (WRB)
  !***END PROLOGUE  PROCP
  REAL Aa, Bd, Bm1, Bm2, rt
  INTEGER ia, ibr, id, j, k, M, m1, m2, mm, mm2, Na, Nd, Nm1, Nm2
  DIMENSION A(*), B(*), C(*), X(*), Y(*), D(*), U(*), Bd(*), Bm1(*)&
    , Bm2(*), Aa(*), W(*)
  COMPLEX X, Y, A, B, C, D, U, W, den, ym, v, bh, am
  !***FIRST EXECUTABLE STATEMENT  PROCP
  DO j = 1, M
    Y(j) = X(j)
    W(j) = Y(j)
  ENDDO
  mm = M - 1
  mm2 = M - 2
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
    bh = B(M) - rt
    ym = Y(M)
    den = B(1) - rt
    D(1) = C(1)/den
    U(1) = A(1)/den
    W(1) = Y(1)/den
    v = C(M)
    IF ( mm2>=2 ) THEN
      DO j = 2, mm2
        den = B(j) - rt - A(j)*D(j-1)
        D(j) = C(j)/den
        U(j) = -A(j)*U(j-1)/den
        W(j) = (Y(j)-A(j)*W(j-1))/den
        bh = bh - v*U(j-1)
        ym = ym - v*W(j-1)
        v = -v*D(j-1)
      ENDDO
    ENDIF
    den = B(M-1) - rt - A(M-1)*D(M-2)
    D(M-1) = (C(M-1)-A(M-1)*U(M-2))/den
    W(M-1) = (Y(M-1)-A(M-1)*W(M-2))/den
    am = A(M) - v*D(M-2)
    bh = bh - v*U(M-2)
    ym = ym - v*W(M-2)
    den = bh - am*D(M-1)
    IF ( ABS(den)/=0 ) THEN
      W(M) = (ym-am*W(M-1))/den
    ELSE
      W(M) = (1.,0.)
    ENDIF
    W(M-1) = W(M-1) - D(M-1)*W(M)
    DO j = 2, mm
      k = M - j
      W(k) = W(k) - D(k)*W(k+1) - U(k)*W(M)
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
END SUBROUTINE PROCP
