!DECK LCE
LOGICAL FUNCTION LCE(Ri,Rj,Lr)
  IMPLICIT NONE
  !***BEGIN PROLOGUE  LCE
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
  !***END PROLOGUE  LCE
  !     .. Scalar Arguments ..
  INTEGER Lr
  !     .. Array Arguments ..
  COMPLEX Ri(*), Rj(*)
  !     .. Local Scalars ..
  INTEGER i
  !***FIRST EXECUTABLE STATEMENT  LCE
  LCE = .TRUE.
  DO i = 1, Lr
    IF ( Ri(i)/=Rj(i) ) THEN
      LCE = .FALSE.
      EXIT
    ENDIF
  ENDDO
  !
  !     End of LCE.
  !
END FUNCTION LCE
