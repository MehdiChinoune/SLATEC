!** CPROD
SUBROUTINE CPROD(Nd,Bd,Nm1,Bm1,Nm2,Bm2,Na,Aa,X,Yy,M,A,B,C,D,W,Y)
  !>
  !  Subsidiary to BLKTRI
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      SINGLE PRECISION (CPROD-S, CPROC-C)
  !***
  ! **Author:**  (UNKNOWN)
  !***
  ! **Description:**
  !
  ! PROD applies a sequence of matrix operations to the vector X and
  ! stores the result in YY.   (COMPLEX case)
  ! AA         array containing scalar multipliers of the vector X.
  ! ND,NM1,NM2 are the lengths of the arrays BD,BM1,BM2 respectively.
  ! BD,BM1,BM2 are arrays containing roots of certain B polynomials.
  ! NA         is the length of the array AA.
  ! X,YY      The matrix operations are applied to X and the result is YY.
  ! A,B,C      are arrays which contain the tridiagonal matrix.
  ! M          is the order of the matrix.
  ! D,W,Y      are working arrays.
  ! ISGN       determines whether or not a change in sign is made.
  !
  !***
  ! **See also:**  BLKTRI
  !***
  ! **Routines called:**  (NONE)

  !* REVISION HISTORY  (YYMMDD)
  !   801001  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900402  Added TYPE section.  (WRB)

  INTEGER :: M, Na, Nd, Nm1, Nm2
  REAL(SP) :: A(M), Aa(Na), B(M), Bm1(Nm1), Bm2(Nm2), C(M), X(M), Yy(M)
  COMPLEX(SP) :: Y(M), D(M), W(M), Bd(Nd)
  INTEGER :: ia, id, iflg, j, k, m1, m2, mm
  REAL(SP) :: rt
  COMPLEX(SP) :: crt, den, y1, y2
  !* FIRST EXECUTABLE STATEMENT  CPROD
  DO j = 1, M
    Y(j) = CMPLX(X(j),0.)
  END DO
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
    END DO
    den = B(1) - crt - C(1)*D(2)
    IF ( ABS(den)/=0 ) THEN
      Y(1) = (Y(1)-C(1)*W(2))/den
    ELSE
      Y(1) = (1.,0.)
    END IF
    DO j = 2, M
      Y(j) = W(j) - D(j)*Y(j-1)
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
      DO j = 1, M
        Yy(j) = REAL(Y(j))
      END DO
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
  y1 = (B(1)-rt)*Y(1) + C(1)*Y(2)
  IF ( mm>=2 ) THEN
    !
    ! MATRIX MULTIPLICATION
    !
    DO j = 2, mm
      y2 = A(j)*Y(j-1) + (B(j)-rt)*Y(j) + C(j)*Y(j+1)
      Y(j-1) = y1
      y1 = y2
    END DO
  END IF
  Y(M) = A(M)*Y(M-1) + (B(M)-rt)*Y(M)
  Y(M-1) = y1
  iflg = 1
  GOTO 100
  RETURN
END SUBROUTINE CPROD
