*DECK SDRES2
      SUBROUTINE SDRES2 (T, Y, YPRIME, DELTA, IRES, RPAR, IPAR)
C***BEGIN PROLOGUE  SDRES2
C***SUBSIDIARY
C***PURPOSE  Second residual evaluator for SDASQC.
C***LIBRARY   SLATEC (DASSL)
C***TYPE      SINGLE PRECISION (SDRES2-S, DDRES2-D)
C***AUTHOR  PETZOLD, LINDA R., (LLNL)
C***SEE ALSO  SDASQC
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   891013  DATE WRITTEN
C   901001  Converted prologue to 4.0 format and made all argument
C           declarations explicit.  (FNF)
C   901030  Made all local declarations explicit.  (FNF)
C***END PROLOGUE  SDRES2
      INTEGER  IRES, IPAR(*)
      REAL  T, Y(*), YPRIME(*), DELTA(*), RPAR(*)
      INTEGER  I, J, K, NG
      REAL  ALPH1, ALPH2, D
      DATA ALPH1/1.0E0/, ALPH2/1.0E0/, NG/5/
C***FIRST EXECUTABLE STATEMENT  SDRES2
      DO 10 J = 1,NG
      DO 10 I = 1,NG
        K = I + (J - 1)*NG
        D = -2.0E0*Y(K)
        IF (I .NE. 1) D = D + Y(K-1)*ALPH1
        IF (J .NE. 1) D = D + Y(K-NG)*ALPH2
 10     DELTA(K) = D - YPRIME(K)
      RETURN
      END
