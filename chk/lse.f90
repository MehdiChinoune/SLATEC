!DECK LSE
LOGICAL FUNCTION LSE(Ri,Rj,Lr)
  IMPLICIT NONE
  !***BEGIN PROLOGUE  LSE
  !***SUBSIDIARY
  !***PURPOSE  Test if two arrays are identical.
  !***LIBRARY   SLATEC (BLAS)
  !***AUTHOR  Du Croz, J. J., (NAG)
  !           Hanson, R. J., (SNLA)
  !***DESCRIPTION
  !
  !  Tests if two arrays are identical.
  !
  !  Auxiliary routine for test program for Level 2 Blas.
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   870810  DATE WRITTEN
  !   910620  Modified to meet SLATEC code and prologue standards.  (BKS)
  !***END PROLOGUE  LSE
  !     .. Scalar Arguments ..
  INTEGER Lr
  !     .. Array Arguments ..
  REAL Ri(*), Rj(*)
  !     .. Local Scalars ..
  INTEGER i
  !***FIRST EXECUTABLE STATEMENT  LSE
  LSE = .TRUE.
  DO i = 1, Lr
    IF ( Ri(i)/=Rj(i) ) THEN
      LSE = .FALSE.
      EXIT
    ENDIF
  ENDDO
  !
  !     End of LSE.
  !
END FUNCTION LSE
