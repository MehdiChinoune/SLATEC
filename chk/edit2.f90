!DECK EDIT2
SUBROUTINE EDIT2(Y,T,Erm)
  IMPLICIT NONE
  !***BEGIN PROLOGUE  EDIT2
  !***SUBSIDIARY
  !***PURPOSE  Subsidiary to SDASQC.
  !***LIBRARY   SLATEC (DASSL)
  !***TYPE      SINGLE PRECISION (EDIT2-S, DEDIT2-D)
  !***AUTHOR  PETZOLD, LINDA R., (LLNL)
  !***SEE ALSO  SDASQC
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   891013  DATE WRITTEN
  !   901001  Converted prologue to 4.0 format and made all argument
  !           declarations explicit.  (FNF)
  !   901009  Changed AMAX1 to MAX.  (FNF)
  !   901030  Removed FLOAT's; made all local declarations explicit. (FNF)
  !***END PROLOGUE  EDIT2
  REAL Y(*), T, Erm
  INTEGER i, j, k, ng
  REAL alph1, alph2, a1, a2, er, ex, yt
  DATA alph1/1.0E0/, alph2/1.0E0/, ng/5/
  !***FIRST EXECUTABLE STATEMENT  EDIT2
  Erm = 0.0E0
  IF ( T==0.0E0 ) RETURN
  ex = 0.0E0
  IF ( T<=30.0E0 ) ex = EXP(-2.0E0*T)
  a2 = 1.0E0
  DO j = 1, ng
    a1 = 1.0E0
    DO i = 1, ng
      k = i + (j-1)*ng
      yt = T**(i+j-2)*ex*a1*a2
      er = ABS(Y(k)-yt)
      Erm = MAX(Erm,er)
      a1 = a1*alph1/i
    ENDDO
    a2 = a2*alph2/j
  ENDDO
END SUBROUTINE EDIT2
