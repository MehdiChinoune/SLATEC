*DECK T2
      REAL FUNCTION T2 (X)
C***BEGIN PROLOGUE  T2
C***PURPOSE  Subsidiary to
C***LIBRARY   SLATEC
C***AUTHOR  (UNKNOWN)
C***ROUTINES CALLED  F2S
C***REVISION HISTORY  (YYMMDD)
C   ??????  DATE WRITTEN
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C***END PROLOGUE  T2
      REAL A,B,F2S,X,X1,Y
C***FIRST EXECUTABLE STATEMENT  T2
      A = 0.1E+00
      B = 0.1E+01
      X1 = X+0.1E+01
      Y = (B-A)/X1+A
      T2 = (B-A)*F2S(Y)/X1/X1
      RETURN
      END
