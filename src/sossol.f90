!*==SOSSOL.f90  processed by SPAG 6.72Dc at 11:02 on  6 Feb 2019
!DECK SOSSOL
SUBROUTINE SOSSOL(K,N,L,X,C,B,M)
  IMPLICIT NONE
  !*--SOSSOL5
  !*** Start of declarations inserted by SPAG
  REAL B , C , X , xmax
  INTEGER j , jkm , K , kj , km , km1 , kmm1 , kn , L , lk , M , N , np1
  !*** End of declarations inserted by SPAG
  !***BEGIN PROLOGUE  SOSSOL
  !***SUBSIDIARY
  !***PURPOSE  Subsidiary to SOS
  !***LIBRARY   SLATEC
  !***TYPE      SINGLE PRECISION (SOSSOL-S, DSOSSL-D)
  !***AUTHOR  (UNKNOWN)
  !***DESCRIPTION
  !
  !     SOSSOL solves an upper triangular type of linear system by back
  !     substitution.
  !
  !     The matrix C is upper trapezoidal and stored as a linear array by
  !     rows. The equations have been normalized so that the diagonal
  !     entries of C are understood to be unity. The off diagonal entries
  !     and the elements of the constant right hand side vector B have
  !     already been stored as the negatives of the corresponding equation
  !     values.
  !     with each call to SOSSOL a (K-1) by (K-1) triangular system is
  !     resolved. For L greater than K, column L of C is included in the
  !     right hand side vector.
  !
  !***SEE ALSO  SOS
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   801001  DATE WRITTEN
  !   890831  Modified array declarations.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900328  Added TYPE section.  (WRB)
  !***END PROLOGUE  SOSSOL
  !
  !
  DIMENSION X(*) , C(*) , B(*)
  !
  !***FIRST EXECUTABLE STATEMENT  SOSSOL
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
    xmax = 0.
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
END SUBROUTINE SOSSOL
