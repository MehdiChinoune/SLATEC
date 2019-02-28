!DECK DLPDP
SUBROUTINE DLPDP(A,Mda,M,N1,N2,Prgopt,X,Wnorm,Mode,Ws,Is)
  IMPLICIT NONE
  !***BEGIN PROLOGUE  DLPDP
  !***SUBSIDIARY
  !***PURPOSE  Subsidiary to DLSEI
  !***LIBRARY   SLATEC
  !***TYPE      DOUBLE PRECISION (LPDP-S, DLPDP-D)
  !***AUTHOR  Hanson, R. J., (SNLA)
  !           Haskell, K. H., (SNLA)
  !***DESCRIPTION
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
  !***SEE ALSO  DLSEI
  !***ROUTINES CALLED  DCOPY, DDOT, DNRM2, DSCAL, DWNNLS
  !***REVISION HISTORY  (YYMMDD)
  !   790701  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900328  Added TYPE section.  (WRB)
  !   910408  Updated the AUTHOR section.  (WRB)
  !***END PROLOGUE  DLPDP
  !
  INTEGER i, Is(*), iw, ix, j, l, M, Mda, Mode, modew, n, N1, &
    N2, np1
  REAL(8) :: A(Mda,*), DDOT, DNRM2, fac, one, Prgopt(*), rnorm, &
    sc, Wnorm, Ws(*), X(*), ynorm, zero
  SAVE zero, one, fac
  DATA zero, one/0.0D0, 1.0D0/, fac/0.1D0/
  !***FIRST EXECUTABLE STATEMENT  DLPDP
  n = N1 + N2
  Mode = 1
  IF ( M>0 ) THEN
    !        BEGIN BLOCK PERMITTING ...EXITS TO 190
    np1 = n + 1
    !
    !           SCALE NONZERO ROWS OF INEQUALITY MATRIX TO HAVE LENGTH ONE.
    DO i = 1, M
      sc = DNRM2(n,A(i,1),Mda)
      IF ( sc/=zero ) THEN
        sc = one/sc
        CALL DSCAL(np1,sc,A(i,1),Mda)
      ENDIF
    ENDDO
    !
    !           SCALE RT.-SIDE VECTOR TO HAVE LENGTH ONE (OR ZERO).
    ynorm = DNRM2(M,A(1,np1),1)
    IF ( ynorm/=zero ) THEN
      sc = one/ynorm
      CALL DSCAL(M,sc,A(1,np1),1)
    ENDIF
    !
    !           SCALE COLS OF MATRIX H.
    j = N1 + 1
    DO WHILE ( j<=n )
      sc = DNRM2(M,A(1,j),1)
      IF ( sc/=zero ) sc = one/sc
      CALL DSCAL(M,sc,A(1,j),1)
      X(j) = sc
      j = j + 1
    ENDDO
    IF ( N1>0 ) THEN
      !
      !              COPY TRANSPOSE OF (H G Y) TO WORK ARRAY WS(*).
      iw = 0
      DO i = 1, M
        !
        !                 MOVE COL OF TRANSPOSE OF H INTO WORK ARRAY.
        CALL DCOPY(N2,A(i,N1+1),Mda,Ws(iw+1),1)
        iw = iw + N2
        !
        !                 MOVE COL OF TRANSPOSE OF G INTO WORK ARRAY.
        CALL DCOPY(N1,A(i,1),Mda,Ws(iw+1),1)
        iw = iw + N1
        !
        !                 MOVE COMPONENT OF VECTOR Y INTO WORK ARRAY.
        Ws(iw+1) = A(i,np1)
        iw = iw + 1
      ENDDO
      Ws(iw+1) = zero
      CALL DCOPY(n,Ws(iw+1),0,Ws(iw+1),1)
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
      CALL DWNNLS(Ws,np1,N2,np1-N2,M,0,Prgopt,Ws(ix),rnorm,modew,Is,Ws(iw+1)&
        )
      !
      !              COMPUTE THE COMPONENTS OF THE SOLN DENOTED ABOVE BY W.
      sc = one - DDOT(M,A(1,np1),1,Ws(ix),1)
      IF ( one+fac*ABS(sc)==one.OR.rnorm<=zero ) THEN
        Mode = 2
        !        .........EXIT
        RETURN
      ELSE
        sc = one/sc
        DO j = 1, N1
          X(j) = sc*DDOT(M,A(1,j),1,Ws(ix),1)
        ENDDO
        !
        !                 COMPUTE THE VECTOR Q=Y-GW.  OVERWRITE Y WITH THIS
        !                 VECTOR.
        DO i = 1, M
          A(i,np1) = A(i,np1) - DDOT(N1,A(i,1),Mda,X,1)
        ENDDO
      ENDIF
    ENDIF
    IF ( N2>0 ) THEN
      !
      !              COPY TRANSPOSE OF (H Q) TO WORK ARRAY WS(*).
      iw = 0
      DO i = 1, M
        CALL DCOPY(N2,A(i,N1+1),Mda,Ws(iw+1),1)
        iw = iw + N2
        Ws(iw+1) = A(i,np1)
        iw = iw + 1
      ENDDO
      Ws(iw+1) = zero
      CALL DCOPY(N2,Ws(iw+1),0,Ws(iw+1),1)
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
      sc = one - DDOT(M,A(1,np1),1,Ws(ix),1)
      IF ( one+fac*ABS(sc)==one.OR.rnorm<=zero ) THEN
        Mode = 2
        !        .........EXIT
        RETURN
      ELSE
        sc = one/sc
        DO j = 1, N2
          l = N1 + j
          X(l) = sc*DDOT(M,A(1,l),1,Ws(ix),1)*X(l)
        ENDDO
      ENDIF
    ENDIF
    !
    !           ACCOUNT FOR SCALING OF RT.-SIDE VECTOR IN SOLUTION.
    CALL DSCAL(n,ynorm,X,1)
    Wnorm = DNRM2(N1,X,1)
  ELSE
    IF ( n>0 ) THEN
      X(1) = zero
      CALL DCOPY(n,X,0,X,1)
    ENDIF
    Wnorm = zero
  ENDIF
  RETURN
END SUBROUTINE DLPDP
