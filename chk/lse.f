*DECK LSE
      LOGICAL FUNCTION LSE (RI, RJ, LR)
C***BEGIN PROLOGUE  LSE
C***SUBSIDIARY
C***PURPOSE  Test if two arrays are identical.
C***LIBRARY   SLATEC (BLAS)
C***AUTHOR  Du Croz, J. J., (NAG)
C           Hanson, R. J., (SNLA)
C***DESCRIPTION
C
C  Tests if two arrays are identical.
C
C  Auxiliary routine for test program for Level 2 Blas.
C***REFERENCES  (NONE)
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   870810  DATE WRITTEN
C   910620  Modified to meet SLATEC code and prologue standards.  (BKS)
C***END PROLOGUE  LSE
C     .. Scalar Arguments ..
      INTEGER            LR
C     .. Array Arguments ..
      REAL               RI( * ), RJ( * )
C     .. Local Scalars ..
      INTEGER            I
C***FIRST EXECUTABLE STATEMENT  LSE
      LSE = .TRUE.
      DO 10 I = 1, LR
         IF( RI( I ).NE.RJ( I ) ) THEN
           LSE = .FALSE.
           GO TO 30
         ENDIF
   10 CONTINUE
   30 RETURN
C
C     End of LSE.
C
      END
