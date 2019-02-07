!*==LDERES.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK LDERES
LOGICAL FUNCTION LDERES(Type,Uplo,M,N,Aa,As,Lda)
  IMPLICIT NONE
  !*--LDERES5
  !***BEGIN PROLOGUE  LDERES
  !***SUBSIDIARY
  !***PURPOSE  Test if selected elements in two arrays are equal.
  !***LIBRARY   SLATEC (BLAS)
  !***AUTHOR  Du Croz, J. J., (NAG)
  !           Hanson, R. J., (SNLA)
  !***DESCRIPTION
  !
  !  Tests if selected elements in two arrays are equal.
  !
  !  TYPE is 'GE', 'SY' or 'SP'.
  !
  !  Auxiliary routine for test program for Level 2 Blas.
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   870810  DATE WRITTEN
  !   910620  Modified to meet SLATEC code and prologue standards.  (BKS)
  !***END PROLOGUE  LDERES
  !     .. Scalar Arguments ..
  INTEGER Lda , M , N
  CHARACTER :: Uplo
  CHARACTER(2) :: Type
  !     .. Array Arguments ..
  DOUBLE PRECISION Aa(Lda,*) , As(Lda,*)
  !     .. Local Scalars ..
  INTEGER i , ibeg , iend , j
  LOGICAL upper
  !***FIRST EXECUTABLE STATEMENT  LDERES
  upper = Uplo=='U'
  IF ( Type=='GE' ) THEN
    DO j = 1 , N
      DO i = M + 1 , Lda
        IF ( Aa(i,j)/=As(i,j) ) GOTO 100
      ENDDO
    ENDDO
  ELSEIF ( Type=='SY' ) THEN
    DO j = 1 , N
      IF ( upper ) THEN
        ibeg = 1
        iend = j
      ELSE
        ibeg = j
        iend = N
      ENDIF
      DO i = 1 , ibeg - 1
        IF ( Aa(i,j)/=As(i,j) ) GOTO 100
      ENDDO
      DO i = iend + 1 , Lda
        IF ( Aa(i,j)/=As(i,j) ) GOTO 100
      ENDDO
    ENDDO
  ENDIF
  !
  LDERES = .TRUE.
  GOTO 99999
  100  LDERES = .FALSE.
  !
  !     End of LDERES.
  !
  99999 END FUNCTION LDERES
