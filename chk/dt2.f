*DECK DT2
      DOUBLE PRECISION FUNCTION DT2 (X)
C***BEGIN PROLOGUE  DT2
C***PURPOSE  Subsidiary to
C***LIBRARY   SLATEC
C***AUTHOR  (UNKNOWN)
C***ROUTINES CALLED  DF2S
C***REVISION HISTORY  (YYMMDD)
C   ??????  DATE WRITTEN
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C***END PROLOGUE  DT2
      DOUBLE PRECISION A,B,DF2S,X,X1,Y
C***FIRST EXECUTABLE STATEMENT  DT2
      A = 0.1D+00
      B = 0.1D+01
      X1 = X+0.1D+01
      Y = (B-A)/X1+A
      DT2 = (B-A)*DF2S(Y)/X1/X1
      RETURN
      END
