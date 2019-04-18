!** CPROCP
SUBROUTINE CPROCP(Nd,Bd,Nm1,Bm1,Nm2,Bm2,Na,Aa,X,Y,M,A,B,C,D,U,Yy)
  !>
  !  Subsidiary to CBLKTR
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      COMPLEX (CPRODP-S, CPROCP-C)
  !***
  ! **Author:**  (UNKNOWN)
  !***
  ! **Description:**
  !
  ! CPROCP applies a sequence of matrix operations to the vector X and
  ! stores the result in Y.
  !
  ! BD,BM1,BM2  are arrays containing roots of certain B polynomials.
  ! ND,NM1,NM2  are the lengths of the arrays BD,BM1,BM2 respectively.
  ! AA          Array containing scalar multipliers of the vector X.
  ! NA          is the length of the array AA.
  ! X,Y        The matrix operations are applied to X and the result is Y.
  ! A,B,C       are arrays which contain the tridiagonal matrix.
  ! M           is the order of the matrix.
  ! D,U         are work arrays.
  ! ISGN        determines whether or not a change in sign is made.
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

  REAL Aa(*), Bm1(*), Bm2(*), rt, Yy(*)
  INTEGER ia, id, iflg, j, k, M, m1, m2, mm, mm2, Na, Nd, Nm1, Nm2
  COMPLEX Y(*), D(*), U(*), v, den, bh, ym, am, y1, y2, yh, Bd(*), crt, X(*), &
    A(*), B(*), C(*)
  !* FIRST EXECUTABLE STATEMENT  CPROCP
  DO j = 1, M
    Y(j) = X(j)
  END DO
  mm = M - 1
  mm2 = M - 2
  id = Nd
  m1 = Nm1
  m2 = Nm2
  ia = Na
  100 CONTINUE
  IFlg = 0
  IF ( id>0 ) THEN
    crt = Bd(id)
    id = id - 1
    iflg = 1
    !
    ! BEGIN SOLUTION TO SYSTEM
    !
    bh = B(M) - crt
    ym = Y(M)
    den = B(1) - crt
    D(1) = C(1)/den
    U(1) = A(1)/den
    Y(1) = Y(1)/den
    v = C(M)
    IF ( mm2>=2 ) THEN
      DO j = 2, mm2
        den = B(j) - crt - A(j)*D(j-1)
        D(j) = C(j)/den
        U(j) = -A(j)*U(j-1)/den
        Y(j) = (Y(j)-A(j)*Y(j-1))/den
        bh = bh - v*U(j-1)
        ym = ym - v*Y(j-1)
        v = -v*D(j-1)
      END DO
    END IF
    den = B(M-1) - crt - A(M-1)*D(M-2)
    D(M-1) = (C(M-1)-A(M-1)*U(M-2))/den
    Y(M-1) = (Y(M-1)-A(M-1)*Y(M-2))/den
    am = A(M) - v*D(M-2)
    bh = bh - v*U(M-2)
    ym = ym - v*Y(M-2)
    den = bh - am*D(M-1)
    IF ( ABS(den)/=0 ) THEN
      Y(M) = (ym-am*Y(M-1))/den
    ELSE
      Y(M) = (1.,0.)
    END IF
    Y(M-1) = Y(M-1) - D(M-1)*Y(M)
    DO j = 2, mm
      k = M - j
      Y(k) = Y(k) - D(k)*Y(k+1) - U(k)*Y(M)
    END DO
  END IF
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
        END DO
      END IF
      IF ( iflg>0 ) GOTO 100
      RETURN
    ELSE
      rt = Bm2(m2)
      m2 = m2 - 1
    END IF
  ELSEIF ( m2<=0 ) THEN
    rt = Bm1(m1)
    m1 = m1 - 1
  ELSEIF ( ABS(Bm1(m1))<=ABS(Bm2(m2)) ) THEN
    rt = Bm2(m2)
    m2 = m2 - 1
  ELSE
    rt = Bm1(m1)
    m1 = m1 - 1
  END IF
  !
  ! MATRIX MULTIPLICATION
  !
  yh = Y(1)
  y1 = (B(1)-rt)*Y(1) + C(1)*Y(2) + A(1)*Y(M)
  IF ( mm>=2 ) THEN
    DO j = 2, mm
      y2 = A(j)*Y(j-1) + (B(j)-rt)*Y(j) + C(j)*Y(j+1)
      Y(j-1) = y1
      y1 = y2
    END DO
  END IF
  Y(M) = A(M)*Y(M-1) + (B(M)-rt)*Y(M) + C(M)*yh
  Y(M-1) = y1
  iflg = 1
  GOTO 100
  RETURN
END SUBROUTINE CPROCP
