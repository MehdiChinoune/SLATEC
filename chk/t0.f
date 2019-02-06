*DECK T0
      REAL FUNCTION T0 (X)
C***BEGIN PROLOGUE  T0
C***PURPOSE  Subsidiary to
C***LIBRARY   SLATEC
C***AUTHOR  (UNKNOWN)
C***ROUTINES CALLED  F0S
C***REVISION HISTORY  (YYMMDD)
C   ??????  DATE WRITTEN
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C***END PROLOGUE  T0
      REAL A,B,F0S,X,X1,Y
C***FIRST EXECUTABLE STATEMENT  T0
      A = 0.0E+00
      B = 0.1E+01
      X1 = X+0.1E+01
      Y = (B-A)/X1+A
      T0 = (B-A)*F0S(Y)/X1/X1
      RETURN
      END
