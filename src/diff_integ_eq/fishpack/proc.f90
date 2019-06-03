!** PROC
SUBROUTINE PROC(Nd,Bd,Nm1,Bm1,Nm2,Bm2,Na,Aa,X,Y,M,A,B,C,D,W,U)
  !>
  !  Subsidiary to CBLKTR
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      COMPLEX (PROD-S, PROC-C)
  !***
  ! **Author:**  (UNKNOWN)
  !***
  ! **Description:**
  !
  ! PROC applies a sequence of matrix operations to the vector X and
  !  stores the result in Y.
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
  !***
  ! **See also:**  CBLKTR
  !***
  ! **Routines called:**  (NONE)

  !* REVISION HISTORY  (YYMMDD)
  !   801001  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900402  Added TYPE section.  (WRB)

  INTEGER :: M, Na, Nd, Nm1, Nm2
  REAL(SP) :: Aa(Na), Bd(Nd), Bm1(Nm1), Bm2(Nm2)
  COMPLEX(SP) :: X(M), Y(M), A(M), B(M), C(M), D(M), W(M), U(M)
  INTEGER :: ia, ibr, id, j, k, m1, m2, mm
  REAL(SP) :: rt
  COMPLEX(SP) :: den
  !* FIRST EXECUTABLE STATEMENT  PROC
  DO j = 1, M
    W(j) = X(j)
    Y(j) = W(j)
  END DO
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
      END DO
    END IF
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
    END DO
    den = B(1) - rt - C(1)*D(2)
    W(1) = (1.,0.)
    IF ( ABS(den)/=0 ) W(1) = (Y(1)-C(1)*W(2))/den
    DO j = 2, M
      W(j) = W(j) - D(j)*W(j-1)
    END DO
    IF ( Na<=0 ) THEN
      IF ( m1<=0 ) THEN
        IF ( m2>0 ) GOTO 400
        EXIT
      ELSE
        IF ( m2<=0 ) GOTO 300
        IF ( ABS(Bm1(m1))>ABS(Bm2(m2)) ) GOTO 300
        GOTO 400
      END IF
    END IF
  END DO
  200 CONTINUE
  DO j = 1, M
    Y(j) = W(j)
  END DO
  ibr = 1
  GOTO 100
  300 CONTINUE
  IF ( ibr<=0 ) THEN
    IF ( ABS(Bm1(m1)-Bd(id))<ABS(Bm1(m1)-rt) ) GOTO 200
  END IF
  rt = rt - Bm1(m1)
  m1 = m1 - 1
  GOTO 500
  400 CONTINUE
  IF ( ibr<=0 ) THEN
    IF ( ABS(Bm2(m2)-Bd(id))<ABS(Bm2(m2)-rt) ) GOTO 200
  END IF
  rt = rt - Bm2(m2)
  m2 = m2 - 1
  500 CONTINUE
  DO j = 1, M
    Y(j) = Y(j) + rt*W(j)
  END DO
  GOTO 100
  RETURN
END SUBROUTINE PROC
