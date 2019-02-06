*DECK SBEG
      REAL FUNCTION SBEG (RESET)
C***BEGIN PROLOGUE  SBEG
C***SUBSIDIARY
C***PURPOSE  Generate random numbers.
C***LIBRARY   SLATEC (BLAS)
C***AUTHOR  Du Croz, J. (NAG)
C           Hanson, R. J. (SNLA)
C***DESCRIPTION
C
C  Generates random numbers uniformly distributed between -0.5 and 0.5.
C
C  Auxiliary routine for test program for Level 2 Blas.
C***REFERENCES  (NONE)
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   870810  DATE WRITTEN
C   910619  Modified to meet SLATEC code and prologue standards.  (BKS)
C***END PROLOGUE  SBEG
C     .. Scalar Arguments ..
      LOGICAL            RESET
C     .. Local Scalars ..
      INTEGER            I, IC, MI
C     .. Save statement ..
      SAVE               I, IC, MI
C     .. Intrinsic Functions ..
      INTRINSIC          REAL
C ***FIRST EXECUTABLE STATEMENT  SBEG
      IF( RESET )THEN
C        Initialize local variables.
         MI = 891
         I = 7
         IC = 0
         RESET = .FALSE.
      END IF
C
C     The sequence of values of I is bounded between 1 and 999.
C     If initial I = 1,2,3,6,7 or 9, the period will be 50.
C     If initial I = 4 or 8, the period will be 25.
C     If initial I = 5, the period will be 10.
C     IC is used to break up the period by skipping 1 value of I in 6.
C
      IC = IC + 1
   10 I = I*MI
      I = I - 1000*( I/1000 )
      IF( IC.GE.5 )THEN
         IC = 0
         GO TO 10
      END IF
      SBEG = REAL( I - 500 )/1001.0
      RETURN
C
C     End of SBEG.
C
      END
