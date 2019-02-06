*DECK DDJAC2
      SUBROUTINE DDJAC2 (T, Y, YPRIME, PD, CJ, RPAR, IPAR)
C***BEGIN PROLOGUE  DDJAC2
C***SUBSIDIARY
C***PURPOSE  Second Jacobian evaluator for DDASQC.
C***LIBRARY   SLATEC (DASSL)
C***TYPE      DOUBLE PRECISION (SDJAC2-S, DDJAC2-D)
C***AUTHOR  PETZOLD, LINDA R., (LLNL)
C***SEE ALSO  DDASQC
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   891013  DATE WRITTEN
C   901001  Converted prologue to 4.0 format and made all argument
C           declarations explicit.  (FNF)
C   901001  Eliminated 7-character variable names MBANDPn by explicitly
C           including MBAND+n in expressions.  (FNF)
C   901030  Made all local declarations explicit.  (FNF)
C***END PROLOGUE  DDJAC2
      INTEGER  IPAR(*)
      DOUBLE PRECISION  T, Y(*), YPRIME(*), PD(11,25), CJ, RPAR(*)
      INTEGER  J, MBAND, ML, MU, NEQ, NG
      DOUBLE PRECISION  ALPH1, ALPH2
      DATA ALPH1/1.0D0/, ALPH2/1.0D0/, NG/5/
      DATA ML/5/, MU/0/, NEQ/25/
C***FIRST EXECUTABLE STATEMENT  DDJAC2
      MBAND = ML + MU + 1
      DO 10 J = 1,NEQ
        PD(MBAND,J) = -2.0D0 - CJ
        PD(MBAND+1,J) = ALPH1
        PD(MBAND+2,J) = 0.0D0
        PD(MBAND+3,J) = 0.0D0
        PD(MBAND+4,J) = 0.0D0
 10     PD(MBAND+5,J) = ALPH2
      DO 20 J = 1,NEQ,NG
 20     PD(MBAND+1,J) = 0.0D0
      RETURN
      END
