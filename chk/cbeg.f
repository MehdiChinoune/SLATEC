*DECK CBEG
      COMPLEX FUNCTION CBEG (RESET)
C***BEGIN PROLOGUE  CBEG
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
C***END PROLOGUE  CBEG
C     .. Scalar Arguments ..
      LOGICAL            RESET
C     .. Local Scalars ..
      INTEGER            I, IC, J, MI, MJ
C     .. Save statement ..
      SAVE               I, IC, J, MI, MJ
C     .. Intrinsic Functions ..
      INTRINSIC          CMPLX
C***FIRST EXECUTABLE STATEMENT  CBEG
      IF( RESET )THEN
C        Initialize local variables.
         MI = 891
         MJ = 457
         I = 7
         J = 7
         IC = 0
         RESET = .FALSE.
      END IF
C
C     The sequence of values of I or J is bounded between 1 and 999.
C     If initial I or J = 1,2,3,6,7 or 9, the period will be 50.
C     If initial I or J = 4 or 8, the period will be 25.
C     If initial I or J = 5, the period will be 10.
C     IC is used to break up the period by skipping 1 value of I or J
C     in 6.
C
      IC = IC + 1
   10 I = I*MI
      J = J*MJ
      I = I - 1000*( I/1000 )
      J = J - 1000*( J/1000 )
      IF( IC.GE.5 )THEN
         IC = 0
         GO TO 10
      END IF
      CBEG = CMPLX( ( I - 500 )/1001.0, ( J - 500 )/1001.0 )
      RETURN
C
C     End of CBEG.
C
      END
