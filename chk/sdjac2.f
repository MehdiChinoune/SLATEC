*DECK SDJAC2
      SUBROUTINE SDJAC2 (T, Y, YPRIME, PD, CJ, RPAR, IPAR)
C***BEGIN PROLOGUE  SDJAC2
C***SUBSIDIARY
C***PURPOSE  Second Jacobian evaluator for SDASQC.
C***LIBRARY   SLATEC (DASSL)
C***TYPE      SINGLE PRECISION (SDJAC2-S, DDJAC2-D)
C***AUTHOR  PETZOLD, LINDA R., (LLNL)
C***SEE ALSO  SDASQC
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   891013  DATE WRITTEN
C   901001  Converted prologue to 4.0 format and made all argument
C           declarations explicit.  (FNF)
C   901001  Eliminated 7-character variable names MBANDPn by explicitly
C           including MBAND+n in expressions.  (FNF)
C   901030  Made all local declarations explicit.  (FNF)
C***END PROLOGUE  SDJAC2
      INTEGER  IPAR(*)
      REAL  T, Y(*), YPRIME(*), PD(11,25), CJ, RPAR(*)
      INTEGER  J, MBAND, ML, MU, NEQ, NG
      REAL  ALPH1, ALPH2
      DATA ALPH1/1.0E0/, ALPH2/1.0E0/, NG/5/
      DATA ML/5/, MU/0/, NEQ/25/
C***FIRST EXECUTABLE STATEMENT  SDJAC2
      MBAND = ML + MU + 1
      DO 10 J = 1,NEQ
        PD(MBAND,J) = -2.0E0 - CJ
        PD(MBAND+1,J) = ALPH1
        PD(MBAND+2,J) = 0.0E0
        PD(MBAND+3,J) = 0.0E0
        PD(MBAND+4,J) = 0.0E0
 10     PD(MBAND+5,J) = ALPH2
      DO 20 J = 1,NEQ,NG
 20     PD(MBAND+1,J) = 0.0E0
      RETURN
      END
