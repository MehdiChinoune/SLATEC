!DECK CPROC
SUBROUTINE CPROC(Nd,Bd,Nm1,Bm1,Nm2,Bm2,Na,Aa,X,Y,M,A,B,C,D,W,Yy)
  IMPLICIT NONE
  REAL Aa, Bm1, Bm2, rt, Yy
  INTEGER ia, id, iflg, j, k, M, m1, m2, mm, Na, Nd, Nm1, Nm2
  !***BEGIN PROLOGUE  CPROC
  !***SUBSIDIARY
  !***PURPOSE  Subsidiary to CBLKTR
  !***LIBRARY   SLATEC
  !***TYPE      COMPLEX (CPROD-S, CPROC-C)
  !***AUTHOR  (UNKNOWN)
  !***DESCRIPTION
  !
  ! PROC applies a sequence of matrix operations to the vector X and
  ! stores the result in Y.
  ! AA     Array containing scalar multipliers of the vector X.
  ! ND,NM1,NM2 are the lengths of the arrays BD,BM1,BM2 respectively.
  ! BD,BM1,BM2 are arrays containing roots of certain B polynomials.
  ! NA     is the length of the array AA.
  ! X,Y    The matrix operations are applied to X and the result is Y.
  ! A,B,C  are arrays which contain the tridiagonal matrix.
  ! M      is the order of the matrix.
  ! D,W    are work arrays.
  ! ISGN   determines whether or not a change in sign is made.
  !
  !***SEE ALSO  CBLKTR
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   801001  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900402  Added TYPE section.  (WRB)
  !***END PROLOGUE  CPROC
  !
  COMPLEX Y, D, W, Bd, crt, den, y1, y2, X, A, B, C
  DIMENSION A(*), B(*), C(*), X(*), Y(*), D(*), W(*), Bd(*), Bm1(*)&
    , Bm2(*), Aa(*), Yy(*)
  !***FIRST EXECUTABLE STATEMENT  CPROC
  DO j = 1, M
    Y(j) = X(j)
  ENDDO
  mm = M - 1
  id = Nd
  m1 = Nm1
  m2 = Nm2
  ia = Na
  100 CONTINUE
  IFlg = 0
  IF ( id>0 ) THEN
    crt = Bd(id)
    id = id - 1
    !
    ! BEGIN SOLUTION TO SYSTEM
    !
    D(M) = A(M)/(B(M)-crt)
    W(M) = Y(M)/(B(M)-crt)
    DO j = 2, mm
      k = M - j
      den = B(k+1) - crt - C(k+1)*D(k+2)
      D(k+1) = A(k+1)/den
      W(k+1) = (Y(k+1)-C(k+1)*W(k+2))/den
    ENDDO
    den = B(1) - crt - C(1)*D(2)
    IF ( ABS(den)/=0 ) THEN
      Y(1) = (Y(1)-C(1)*W(2))/den
    ELSE
      Y(1) = (1.,0.)
    ENDIF
    DO j = 2, M
      Y(j) = W(j) - D(j)*Y(j-1)
    ENDDO
  ENDIF
  IF ( m1<=0 ) THEN
    IF ( m2<=0 ) THEN
      IF ( ia>0 ) THEN
        rt = Aa(ia)
        ia = ia - 1
        iflg = 1
        !
        ! SCALAR MULTIPLICATION
        !
        DO j = 1, M
          Y(j) = rt*Y(j)
        ENDDO
      ENDIF
      IF ( iflg>0 ) GOTO 100
      RETURN
    ELSE
      rt = Bm2(m2)
      m2 = m2 - 1
    ENDIF
  ELSEIF ( m2<=0 ) THEN
    rt = Bm1(m1)
    m1 = m1 - 1
  ELSEIF ( ABS(Bm1(m1))<=ABS(Bm2(m2)) ) THEN
    rt = Bm2(m2)
    m2 = m2 - 1
  ELSE
    rt = Bm1(m1)
    m1 = m1 - 1
  ENDIF
  y1 = (B(1)-rt)*Y(1) + C(1)*Y(2)
  IF ( mm>=2 ) THEN
    !
    ! MATRIX MULTIPLICATION
    !
    DO j = 2, mm
      y2 = A(j)*Y(j-1) + (B(j)-rt)*Y(j) + C(j)*Y(j+1)
      Y(j-1) = y1
      y1 = y2
    ENDDO
  ENDIF
  Y(M) = A(M)*Y(M-1) + (B(M)-rt)*Y(M)
  Y(M-1) = y1
  iflg = 1
  GOTO 100
  RETURN
END SUBROUTINE CPROC
