*DECK T4
      REAL FUNCTION T4 (X)
C***BEGIN PROLOGUE  T4
C***PURPOSE  Subsidiary to
C***LIBRARY   SLATEC
C***AUTHOR  (UNKNOWN)
C***ROUTINES CALLED  F4S
C***REVISION HISTORY  (YYMMDD)
C   ??????  DATE WRITTEN
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C***END PROLOGUE  T4
      REAL A,B,F4S,X,X1,Y
C***FIRST EXECUTABLE STATEMENT  T4
      A = 0.0E+00
      B = 0.1E+01
      X1 = X+0.1E+01
      Y = (B-A)/X1+A
      T4 = (B-A)*F4S(Y)/X1/X1
      RETURN
      END
