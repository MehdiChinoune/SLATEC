!DECK BESKES
SUBROUTINE BESKES(Xnu,X,Nin,Bke)
  IMPLICIT NONE
  REAL alnbig, Bke, bknu1, direct, R1MACH, v, vend, vincr, X, Xnu
  INTEGER i, iswtch, n, Nin
  !***BEGIN PROLOGUE  BESKES
  !***PURPOSE  Compute a sequence of exponentially scaled modified Bessel
  !            functions of the third kind of fractional order.
  !***LIBRARY   SLATEC (FNLIB)
  !***CATEGORY  C10B3
  !***TYPE      SINGLE PRECISION (BESKES-S, DBSKES-D)
  !***KEYWORDS  EXPONENTIALLY SCALED, FNLIB, FRACTIONAL ORDER,
  !             MODIFIED BESSEL FUNCTION, SEQUENCE OF BESSEL FUNCTIONS,
  !             SPECIAL FUNCTIONS, THIRD KIND
  !***AUTHOR  Fullerton, W., (LANL)
  !***DESCRIPTION
  !
  ! BESKES computes a sequence of exponentially scaled
  ! (i.e., multipled by EXP(X)) modified Bessel
  ! functions of the third kind of order XNU + I at X, where X .GT. 0,
  ! XNU lies in (-1,1), and I = 0, 1, ..., NIN - 1, if NIN is positive
  ! and I = 0, -1, ..., NIN + 1, if NIN is negative.  On return, the
  ! vector BKE(.) contains the results at X for order starting at XNU.
  !
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  R1MACH, R9KNUS, XERMSG
  !***REVISION HISTORY  (YYMMDD)
  !   770601  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890911  Removed unnecessary intrinsics.  (WRB)
  !   890911  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900326  Removed duplicate information from DESCRIPTION section.
  !           (WRB)
  !***END PROLOGUE  BESKES
  DIMENSION Bke(*)
  SAVE alnbig
  DATA alnbig/0./
  !***FIRST EXECUTABLE STATEMENT  BESKES
  IF ( alnbig==0. ) alnbig = LOG(R1MACH(2))
  !
  v = ABS(Xnu)
  n = ABS(Nin)
  !
  IF ( v>=1. ) CALL XERMSG('SLATEC','BESKES','ABS(XNU) MUST BE LT 1',2,2)
  IF ( X<=0. ) CALL XERMSG('SLATEC','BESKES','X IS LE 0',3,2)
  IF ( n==0 ) CALL XERMSG('SLATEC','BESKES',&
    'N THE NUMBER IN THE SEQUENCE IS 0',4,2)
  !
  CALL R9KNUS(v,X,Bke(1),bknu1,iswtch)
  IF ( n==1 ) RETURN
  !
  vincr = SIGN(1.0,REAL(Nin))
  direct = vincr
  IF ( Xnu/=0. ) direct = vincr*SIGN(1.0,Xnu)
  IF ( iswtch==1.AND.direct>0. ) CALL XERMSG('SLATEC','BESKES',&
    'X SO SMALL BESSEL K-SUB-XNU+1 OVERFLOWS',5,2)
  Bke(2) = bknu1
  !
  IF ( direct<0. ) CALL R9KNUS(ABS(Xnu+vincr),X,Bke(2),bknu1,iswtch)
  IF ( n==2 ) RETURN
  !
  vend = ABS(Xnu+Nin) - 1.0
  IF ( (vend-0.5)*LOG(vend)+0.27-vend*(LOG(X)-.694)>alnbig )&
    CALL XERMSG('SLATEC','BESKES',&
    'X SO SMALL OR ABS(NU) SO BIG THAT BESSEL K-SUB-NU OVERFLOWS',5,2)
  !
  v = Xnu
  DO i = 3, n
    v = v + vincr
    Bke(i) = 2.0*v*Bke(i-1)/X + Bke(i-2)
  ENDDO
  !
END SUBROUTINE BESKES
