*DECK T1
      REAL FUNCTION T1 (X)
C***BEGIN PROLOGUE  T1
C***PURPOSE  Subsidiary to
C***LIBRARY   SLATEC
C***AUTHOR  (UNKNOWN)
C***ROUTINES CALLED  F1S
C***REVISION HISTORY  (YYMMDD)
C   ??????  DATE WRITTEN
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C***END PROLOGUE  T1
      REAL A,B,F1S,X,X1,Y
C***FIRST EXECUTABLE STATEMENT  T1
      A = 0.0E+00
      B = 0.1E+01
      X1 = X+0.1E+01
      Y = (B-A)/X1+A
      T1 = (B-A)*F1S(Y)/X1/X1
      RETURN
      END
