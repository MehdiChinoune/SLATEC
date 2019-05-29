!** DLPDP
SUBROUTINE DLPDP(A,Mda,M,N1,N2,Prgopt,X,Wnorm,Mode,Ws,Is)
  !>
  !  Subsidiary to DLSEI
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      DOUBLE PRECISION (LPDP-S, DLPDP-D)
  !***
  ! **Author:**  Hanson, R. J., (SNLA)
  !           Haskell, K. H., (SNLA)
  !***
  ! **Description:**
  !
  !  **** Double Precision version of LPDP ****
  !     DIMENSION A(MDA,N+1),PRGOPT(*),X(N),WS((M+2)*(N+7)),IS(M+N+1),
  !     where N=N1+N2.  This is a slight overestimate for WS(*).
  !
  !     Determine an N1-vector W, and
  !               an N2-vector Z
  !     which minimizes the Euclidean length of W
  !     subject to G*W+H*Z .GE. Y.
  !     This is the least projected distance problem, LPDP.
  !     The matrices G and H are of respective
  !     dimensions M by N1 and M by N2.
  !
  !     Called by subprogram DLSI( ).
  !
  !     The matrix
  !                (G H Y)
  !
  !     occupies rows 1,...,M and cols 1,...,N1+N2+1 of A(*,*).
  !
  !     The solution (W) is returned in X(*).
  !                  (Z)
  !
  !     The value of MODE indicates the status of
  !     the computation after returning to the user.
  !
  !          MODE=1  The solution was successfully obtained.
  !
  !          MODE=2  The inequalities are inconsistent.
  !
  !***
  ! **See also:**  DLSEI
  !***
  ! **Routines called:**  DCOPY, DDOT, DNRM2, DSCAL, DWNNLS

  !* REVISION HISTORY  (YYMMDD)
  !   790701  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900328  Added TYPE section.  (WRB)
  !   910408  Updated the AUTHOR section.  (WRB)
  !
  INTEGER :: M, Mda, Mode, N1, N2
  INTEGER :: Is(M+N1+N2+1)
  REAL(8) :: Wnorm, A(Mda,N1+N2+1), Prgopt(:), Ws((M+2)*(N1+N2+7)), X(N1+N2)
  INTEGER :: i, iw, ix, j, l, modew, n, np1
  REAL(8) :: rnorm, sc, ynorm
  REAL(8), PARAMETER :: zero = 0.0D0, one = 1.0D0, fac = 0.1D0
  !* FIRST EXECUTABLE STATEMENT  DLPDP
  n = N1 + N2
  Mode = 1
  IF ( M>0 ) THEN
    !        BEGIN BLOCK PERMITTING ...EXITS TO 190
    np1 = n + 1
    !
    !           SCALE NONZERO ROWS OF INEQUALITY MATRIX TO HAVE LENGTH ONE.
    DO i = 1, M
      sc = NORM2(A(i,1:n))
      IF ( sc/=zero ) THEN
        sc = one/sc
        A(i,1:np1) = A(i,1:np1)*sc
      END IF
    END DO
    !
    !           SCALE RT.-SIDE VECTOR TO HAVE LENGTH ONE (OR ZERO).
    ynorm = NORM2(A(1:M,np1))
    IF ( ynorm/=zero ) THEN
      sc = one/ynorm
      A(1:M,np1) = A(1:M,np1)*sc
    END IF
    !
    !           SCALE COLS OF MATRIX H.
    j = N1 + 1
    DO WHILE ( j<=n )
      sc = NORM2(A(1:M,j))
      IF ( sc/=zero ) sc = one/sc
      A(1:M,j) = A(1:M,j)*sc
      X(j) = sc
      j = j + 1
    END DO
    IF ( N1>0 ) THEN
      !
      !              COPY TRANSPOSE OF (H G Y) TO WORK ARRAY WS(*).
      iw = 0
      DO i = 1, M
        !
        !                 MOVE COL OF TRANSPOSE OF H INTO WORK ARRAY.
        Ws(iw+1:iw+N2) = A(i,N1+1:N1+N2)
        iw = iw + N2
        !
        !                 MOVE COL OF TRANSPOSE OF G INTO WORK ARRAY.
        Ws(iw+1:iw+N1) = A(i,1:N1)
        iw = iw + N1
        !
        !                 MOVE COMPONENT OF VECTOR Y INTO WORK ARRAY.
        Ws(iw+1) = A(i,np1)
        iw = iw + 1
      END DO
      Ws(iw+1:iw+n) = zero
      iw = iw + n
      Ws(iw+1) = one
      iw = iw + 1
      !
      !              SOLVE EU=F SUBJECT TO (TRANSPOSE OF H)U=0, U.GE.0.  THE
      !              MATRIX E = TRANSPOSE OF (G Y), AND THE (N+1)-VECTOR
      !              F = TRANSPOSE OF (0,...,0,1).
      ix = iw + 1
      iw = iw + M
      !
      !              DO NOT CHECK LENGTHS OF WORK ARRAYS IN THIS USAGE OF
      !              DWNNLS( ).
      Is(1) = 0
      Is(2) = 0
      CALL DWNNLS(Ws,np1,N2,np1-N2,M,0,Prgopt,Ws(ix),rnorm,modew,Is,Ws(iw+1))
      !
      !              COMPUTE THE COMPONENTS OF THE SOLN DENOTED ABOVE BY W.
      sc = one - DOT_PRODUCT(A(1:M,np1),Ws(ix:ix+M-1))
      IF ( one+fac*ABS(sc)==one.OR.rnorm<=zero ) THEN
        Mode = 2
        !        .........EXIT
        RETURN
      ELSE
        sc = one/sc
        DO j = 1, N1
          X(j) = sc*DOT_PRODUCT(A(1:M,j),Ws(ix:ix+M-1))
        END DO
        !
        !                 COMPUTE THE VECTOR Q=Y-GW.  OVERWRITE Y WITH THIS
        !                 VECTOR.
        DO i = 1, M
          A(i,np1) = A(i,np1) - DOT_PRODUCT(A(i,1:N1),X(1:N1))
        END DO
      END IF
    END IF
    IF ( N2>0 ) THEN
      !
      !              COPY TRANSPOSE OF (H Q) TO WORK ARRAY WS(*).
      iw = 0
      DO i = 1, M
        Ws(iw+1:iw+N2) = A(i,N1+1:N1+N2)
        iw = iw + N2
        Ws(iw+1) = A(i,np1)
        iw = iw + 1
      END DO
      Ws(iw+1:iw+N2) = zero
      iw = iw + N2
      Ws(iw+1) = one
      iw = iw + 1
      ix = iw + 1
      iw = iw + M
      !
      !              SOLVE RV=S SUBJECT TO V.GE.0.  THE MATRIX R =(TRANSPOSE
      !              OF (H Q)), WHERE Q=Y-GW.  THE (N2+1)-VECTOR S =(TRANSPOSE
      !              OF (0,...,0,1)).
      !
      !              DO NOT CHECK LENGTHS OF WORK ARRAYS IN THIS USAGE OF
      !              DWNNLS( ).
      Is(1) = 0
      Is(2) = 0
      CALL DWNNLS(Ws,N2+1,0,N2+1,M,0,Prgopt,Ws(ix),rnorm,modew,Is,Ws(iw+1))
      !
      !              COMPUTE THE COMPONENTS OF THE SOLN DENOTED ABOVE BY Z.
      sc = one - DOT_PRODUCT(A(1:M,np1),Ws(ix:ix+M-1))
      IF ( one+fac*ABS(sc)==one.OR.rnorm<=zero ) THEN
        Mode = 2
        !        .........EXIT
        RETURN
      ELSE
        sc = one/sc
        DO j = 1, N2
          l = N1 + j
          X(l) = sc*DOT_PRODUCT(A(1:M,l),Ws(ix:ix+M-1))*X(l)
        END DO
      END IF
    END IF
    !
    !           ACCOUNT FOR SCALING OF RT.-SIDE VECTOR IN SOLUTION.
    X(1:n) = X(1:n)*ynorm
    Wnorm = NORM2(X(1:N1))
  ELSE
    IF ( n>0 ) THEN
      X(1:n) = zero
    END IF
    Wnorm = zero
  END IF
  RETURN
END SUBROUTINE DLPDP
