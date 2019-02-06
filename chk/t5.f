*DECK T5
      REAL FUNCTION T5 (X)
C***BEGIN PROLOGUE  T5
C***PURPOSE  Subsidiary to
C***LIBRARY   SLATEC
C***AUTHOR  (UNKNOWN)
C***ROUTINES CALLED  F5S
C***REVISION HISTORY  (YYMMDD)
C   ??????  DATE WRITTEN
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C***END PROLOGUE  T5
      REAL A,B,F5S,X,X1,Y
C***FIRST EXECUTABLE STATEMENT  T5
      A = 0.0E+00
      B = 0.1E+01
      X1 = X+0.1E+01
      Y = (B-A)/X1+A
      T5 = (B-A)*F5S(Y)/X1/X1
      RETURN
      END
