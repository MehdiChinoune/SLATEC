!*==LCERES.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK LCERES
LOGICAL FUNCTION LCERES(Type,Uplo,M,N,Aa,As,Lda)
  IMPLICIT NONE
  !*--LCERES5
  !***BEGIN PROLOGUE  LCERES
  !***SUBSIDIARY
  !***PURPOSE  Test if selected elements in two arrays are equal.
  !***LIBRARY   SLATEC (BLAS)
  !***AUTHOR  Du Croz, J. J., (NAG)
  !           Hanson, R. J., (SNLA)
  !***DESCRIPTION
  !
  !  Tests if selected elements in two arrays are equal.
  !
  !  TYPE is 'GE', 'HE' or 'HP'.
  !
  !  Auxiliary routine for test program for Level 2 Blas.
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   870810  DATE WRITTEN
  !   910620  Modified to meet SLATEC code and prologue standards.  (BKS)
  !***END PROLOGUE  LCERES
  !     .. Scalar Arguments ..
  INTEGER Lda, M, N
  CHARACTER :: Uplo
  CHARACTER(2) :: Type
  !     .. Array Arguments ..
  COMPLEX Aa(Lda,*), As(Lda,*)
  !     .. Local Scalars ..
  INTEGER i, ibeg, iend, j
  LOGICAL upper
  !***FIRST EXECUTABLE STATEMENT  LCERES
  upper = Uplo=='U'
  IF ( Type=='GE' ) THEN
    DO j = 1, N
      DO i = M + 1, Lda
        IF ( Aa(i,j)/=As(i,j) ) GOTO 100
      ENDDO
    ENDDO
  ELSEIF ( Type=='HE' ) THEN
    DO j = 1, N
      IF ( upper ) THEN
        ibeg = 1
        iend = j
      ELSE
        ibeg = j
        iend = N
      ENDIF
      DO i = 1, ibeg - 1
        IF ( Aa(i,j)/=As(i,j) ) GOTO 100
      ENDDO
      DO i = iend + 1, Lda
        IF ( Aa(i,j)/=As(i,j) ) GOTO 100
      ENDDO
    ENDDO
  ENDIF
  !
  LCERES = .TRUE.
  GOTO 99999
  100  LCERES = .FALSE.
  !
  !     End of LCERES.
  !
  99999 CONTINUE
  END FUNCTION LCERES
