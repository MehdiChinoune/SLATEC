*DECK EDIT2
      SUBROUTINE EDIT2 (Y, T, ERM)
C***BEGIN PROLOGUE  EDIT2
C***SUBSIDIARY
C***PURPOSE  Subsidiary to SDASQC.
C***LIBRARY   SLATEC (DASSL)
C***TYPE      SINGLE PRECISION (EDIT2-S, DEDIT2-D)
C***AUTHOR  PETZOLD, LINDA R., (LLNL)
C***SEE ALSO  SDASQC
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   891013  DATE WRITTEN
C   901001  Converted prologue to 4.0 format and made all argument
C           declarations explicit.  (FNF)
C   901009  Changed AMAX1 to MAX.  (FNF)
C   901030  Removed FLOAT's; made all local declarations explicit. (FNF)
C***END PROLOGUE  EDIT2
      REAL  Y(*), T, ERM
      INTEGER  I, J, K, NG
      REAL  ALPH1, ALPH2, A1, A2, ER, EX, YT
      DATA ALPH1/1.0E0/, ALPH2/1.0E0/, NG/5/
C***FIRST EXECUTABLE STATEMENT  EDIT2
      ERM = 0.0E0
      IF (T .EQ. 0.0E0) RETURN
      EX = 0.0E0
      IF (T .LE. 30.0E0) EX = EXP(-2.0E0*T)
      A2 = 1.0E0
      DO 60 J = 1,NG
        A1 = 1.0E0
        DO 50 I = 1,NG
          K = I + (J - 1)*NG
          YT = T**(I+J-2)*EX*A1*A2
          ER = ABS(Y(K)-YT)
          ERM = MAX(ERM,ER)
          A1 = A1*ALPH1/I
 50       CONTINUE
        A2 = A2*ALPH2/J
 60     CONTINUE
      RETURN
      END
