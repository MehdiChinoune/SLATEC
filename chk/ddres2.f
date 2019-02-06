*DECK DDRES2
      SUBROUTINE DDRES2 (T, Y, YPRIME, DELTA, IRES, RPAR, IPAR)
C***BEGIN PROLOGUE  DDRES2
C***SUBSIDIARY
C***PURPOSE  Second residual evaluator for DDASQC.
C***LIBRARY   SLATEC (DASSL)
C***TYPE      DOUBLE PRECISION (SDRES2-S, DDRES2-D)
C***AUTHOR  PETZOLD, LINDA R., (LLNL)
C***SEE ALSO  DDASQC
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   891013  DATE WRITTEN
C   901001  Converted prologue to 4.0 format and made all argument
C           declarations explicit.  (FNF)
C   901030  Made all local declarations explicit.  (FNF)
C***END PROLOGUE  DDRES2
      INTEGER  IRES, IPAR(*)
      DOUBLE PRECISION  T, Y(*), YPRIME(*), DELTA(*), RPAR(*)
      INTEGER  I, J, K, NG
      DOUBLE PRECISION  ALPH1, ALPH2, D
      DATA ALPH1/1.0D0/, ALPH2/1.0D0/, NG/5/
C***FIRST EXECUTABLE STATEMENT  DDRES2
      DO 10 J = 1,NG
      DO 10 I = 1,NG
        K = I + (J - 1)*NG
        D = -2.0D0*Y(K)
        IF (I .NE. 1) D = D + Y(K-1)*ALPH1
        IF (J .NE. 1) D = D + Y(K-NG)*ALPH2
 10     DELTA(K) = D - YPRIME(K)
      RETURN
      END
