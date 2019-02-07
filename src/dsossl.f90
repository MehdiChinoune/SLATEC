!*==DSOSSL.f90  processed by SPAG 6.72Dc at 11:01 on  6 Feb 2019
!DECK DSOSSL
SUBROUTINE DSOSSL(K,N,L,X,C,B,M)
  IMPLICIT NONE
  !*--DSOSSL5
  !***BEGIN PROLOGUE  DSOSSL
  !***SUBSIDIARY
  !***PURPOSE  Subsidiary to DSOS
  !***LIBRARY   SLATEC
  !***TYPE      DOUBLE PRECISION (SOSSOL-S, DSOSSL-D)
  !***AUTHOR  (UNKNOWN)
  !***DESCRIPTION
  !
  !     DSOSSL solves an upper triangular type of linear system by back
  !     substitution.
  !
  !     The matrix C is upper trapezoidal and stored as a linear array by
  !     rows. The equations have been normalized so that the diagonal
  !     entries of C are understood to be unity. The off diagonal entries
  !     and the elements of the constant right hand side vector B have
  !     already been stored as the negatives of the corresponding equation
  !     values.
  !     With each call to DSOSSL a (K-1) by (K-1) triangular system is
  !     resolved. For L greater than K, column L of C is included in the
  !     right hand side vector.
  !
  !***SEE ALSO  DSOS
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   801001  DATE WRITTEN
  !   890831  Modified array declarations.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900328  Added TYPE section.  (WRB)
  !***END PROLOGUE  DSOSSL
  !
  !
  INTEGER j , jkm , K , kj , km , km1 , kmm1 , kn , L , lk , M , N , np1
  REAL(8) :: B(*) , C(*) , X(*) , xmax
  !
  !***FIRST EXECUTABLE STATEMENT  DSOSSL
  np1 = N + 1
  km1 = K - 1
  lk = km1
  IF ( L==K ) lk = K
  kn = M
  !
  !
  DO kj = 1 , km1
    kmm1 = K - kj
    km = kmm1 + 1
    xmax = 0.0D0
    kn = kn - np1 + kmm1
    IF ( km<=lk ) THEN
      jkm = kn
      !
      DO j = km , lk
        jkm = jkm + 1
        xmax = xmax + C(jkm)*X(j)
      ENDDO
    ENDIF
    !
    IF ( L>K ) THEN
      jkm = kn + L - kmm1
      xmax = xmax + C(jkm)*X(L)
    ENDIF
    X(kmm1) = xmax + B(kmm1)
  ENDDO
  !
END SUBROUTINE DSOSSL
