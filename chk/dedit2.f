*DECK DEDIT2
      SUBROUTINE DEDIT2 (Y, T, ERM)
C***BEGIN PROLOGUE  DEDIT2
C***SUBSIDIARY
C***PURPOSE  Subsidiary to DDASQC.
C***LIBRARY   SLATEC (DASSL)
C***TYPE      DOUBLE PRECISION (EDIT2-S, DEDIT2-D)
C***AUTHOR  PETZOLD, LINDA R., (LLNL)
C***SEE ALSO  DDASQC
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   891013  DATE WRITTEN
C   901001  Converted prologue to 4.0 format and made all argument
C           declarations explicit.  (FNF)
C   901009  Changed AMAX1 to MAX.  (FNF)
C   901030  Removed FLOAT's; made all local declarations explicit. (FNF)
C***END PROLOGUE  DEDIT2
      DOUBLE PRECISION  Y(*), T, ERM
      INTEGER  I, J, K, NG
      DOUBLE PRECISION  ALPH1, ALPH2, A1, A2, ER, EX, YT
      DATA ALPH1/1.0D0/, ALPH2/1.0D0/, NG/5/
C***FIRST EXECUTABLE STATEMENT  DEDIT2
      ERM = 0.0D0
      IF (T .EQ. 0.0D0) RETURN
      EX = 0.0D0
      IF (T .LE. 30.0D0) EX = EXP(-2.0D0*T)
      A2 = 1.0D0
      DO 60 J = 1,NG
        A1 = 1.0D0
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
