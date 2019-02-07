!*==LDE.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK LDE
LOGICAL FUNCTION LDE(Ri,Rj,Lr)
  IMPLICIT NONE
  !*--LDE5
  !***BEGIN PROLOGUE  LDE
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
  !***END PROLOGUE  LDE
  !     .. Scalar Arguments ..
  INTEGER Lr
  !     .. Array Arguments ..
  REAL(8) :: Ri(*), Rj(*)
  !     .. Local Scalars ..
  INTEGER i
  !***FIRST EXECUTABLE STATEMENT  LDE
  LDE = .TRUE.
  DO i = 1, Lr
    IF ( Ri(i)/=Rj(i) ) THEN
      LDE = .FALSE.
      EXIT
    ENDIF
  ENDDO
  !
  !     End of LDE.
  !
END FUNCTION LDE
