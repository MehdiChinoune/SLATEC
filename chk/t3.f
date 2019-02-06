*DECK T3
      REAL FUNCTION T3 (X)
C***BEGIN PROLOGUE  T3
C***PURPOSE  Subsidiary to
C***LIBRARY   SLATEC
C***AUTHOR  (UNKNOWN)
C***ROUTINES CALLED  F3S
C***REVISION HISTORY  (YYMMDD)
C   ??????  DATE WRITTEN
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C***END PROLOGUE  T3
      REAL A,B,F3S,X,X1,Y
C***FIRST EXECUTABLE STATEMENT  T3
      A = 0.0E+00
      B = 0.5E+01
      X1 = X+0.1E+01
      Y = (B-A)/X1+A
      T3 = (B-A)*F3S(Y)/X1/X1
      RETURN
      END
