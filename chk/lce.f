*DECK LCE
      LOGICAL FUNCTION LCE (RI, RJ, LR)
C***BEGIN PROLOGUE  LCE
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
C***END PROLOGUE  LCE
C     .. Scalar Arguments ..
      INTEGER            LR
C     .. Array Arguments ..
      COMPLEX            RI( * ), RJ( * )
C     .. Local Scalars ..
      INTEGER            I
C***FIRST EXECUTABLE STATEMENT  LCE
      LCE = .TRUE.
      DO 10 I = 1, LR
         IF( RI( I ).NE.RJ( I ) ) THEN
           LCE = .FALSE.
           GO TO 30
         ENDIF
   10 CONTINUE
   30 RETURN
C
C     End of LCE.
C
      END
