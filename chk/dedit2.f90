!*==DEDIT2.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK DEDIT2
SUBROUTINE DEDIT2(Y,T,Erm)
  IMPLICIT NONE
  !*--DEDIT25
  !***BEGIN PROLOGUE  DEDIT2
  !***SUBSIDIARY
  !***PURPOSE  Subsidiary to DDASQC.
  !***LIBRARY   SLATEC (DASSL)
  !***TYPE      DOUBLE PRECISION (EDIT2-S, DEDIT2-D)
  !***AUTHOR  PETZOLD, LINDA R., (LLNL)
  !***SEE ALSO  DDASQC
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   891013  DATE WRITTEN
  !   901001  Converted prologue to 4.0 format and made all argument
  !           declarations explicit.  (FNF)
  !   901009  Changed AMAX1 to MAX.  (FNF)
  !   901030  Removed FLOAT's; made all local declarations explicit. (FNF)
  !***END PROLOGUE  DEDIT2
  DOUBLE PRECISION Y(*) , T , Erm
  INTEGER i , j , k , ng
  DOUBLE PRECISION alph1 , alph2 , a1 , a2 , er , ex , yt
  DATA alph1/1.0D0/ , alph2/1.0D0/ , ng/5/
  !***FIRST EXECUTABLE STATEMENT  DEDIT2
  Erm = 0.0D0
  IF ( T==0.0D0 ) RETURN
  ex = 0.0D0
  IF ( T<=30.0D0 ) ex = EXP(-2.0D0*T)
  a2 = 1.0D0
  DO j = 1 , ng
    a1 = 1.0D0
    DO i = 1 , ng
      k = i + (j-1)*ng
      yt = T**(i+j-2)*ex*a1*a2
      er = ABS(Y(k)-yt)
      Erm = MAX(Erm,er)
      a1 = a1*alph1/i
    ENDDO
    a2 = a2*alph2/j
  ENDDO
END SUBROUTINE DEDIT2
