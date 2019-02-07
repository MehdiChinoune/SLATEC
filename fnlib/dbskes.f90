!*==DBSKES.f90  processed by SPAG 6.72Dc at 10:56 on  6 Feb 2019
!DECK DBSKES
SUBROUTINE DBSKES(Xnu,X,Nin,Bke)
  IMPLICIT NONE
  !*--DBSKES5
  !*** Start of declarations inserted by SPAG
  INTEGER i , iswtch , n , Nin
  !*** End of declarations inserted by SPAG
  !***BEGIN PROLOGUE  DBSKES
  !***PURPOSE  Compute a sequence of exponentially scaled modified Bessel
  !            functions of the third kind of fractional order.
  !***LIBRARY   SLATEC (FNLIB)
  !***CATEGORY  C10B3
  !***TYPE      DOUBLE PRECISION (BESKES-S, DBSKES-D)
  !***KEYWORDS  EXPONENTIALLY SCALED, FNLIB, FRACTIONAL ORDER,
  !             MODIFIED BESSEL FUNCTION, SEQUENCE OF BESSEL FUNCTIONS,
  !             SPECIAL FUNCTIONS, THIRD KIND
  !***AUTHOR  Fullerton, W., (LANL)
  !***DESCRIPTION
  !
  ! DBSKES(XNU,X,NIN,BKE) computes a double precision sequence
  ! of exponentially scaled modified Bessel functions
  ! of the third kind of order XNU + I at X, where X .GT. 0,
  ! XNU lies in (-1,1), and I = 0, 1, ... , NIN - 1, if NIN is positive
  ! and I = 0, -1, ... , NIN + 1, if NIN is negative.  On return, the
  ! vector BKE(.) contains the results at X for order starting at XNU.
  ! XNU, X, and BKE are double precision.  NIN is integer.
  !
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  D1MACH, D9KNUS, XERMSG
  !***REVISION HISTORY  (YYMMDD)
  !   770601  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890911  Removed unnecessary intrinsics.  (WRB)
  !   890911  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !***END PROLOGUE  DBSKES
  DOUBLE PRECISION Xnu , X , Bke(*) , bknu1 , v , vincr , vend , alnbig , &
    D1MACH , direct
  SAVE alnbig
  DATA alnbig/0.D0/
  !***FIRST EXECUTABLE STATEMENT  DBSKES
  IF ( alnbig==0.D0 ) alnbig = LOG(D1MACH(2))
  !
  v = ABS(Xnu)
  n = ABS(Nin)
  !
  IF ( v>=1.D0 ) CALL XERMSG('SLATEC','DBSKES','ABS(XNU) MUST BE LT 1',2,2)
  IF ( X<=0.D0 ) CALL XERMSG('SLATEC','DBSKES','X IS LE 0',3,2)
  IF ( n==0 ) CALL XERMSG('SLATEC','DBSKES',&
    'N THE NUMBER IN THE SEQUENCE IS 0',4,2)
  !
  CALL D9KNUS(v,X,Bke(1),bknu1,iswtch)
  IF ( n==1 ) RETURN
  !
  vincr = SIGN(1.0,REAL(Nin))
  direct = vincr
  IF ( Xnu/=0.D0 ) direct = vincr*SIGN(1.D0,Xnu)
  IF ( iswtch==1.AND.direct>0. ) CALL XERMSG('SLATEC','DBSKES',&
    'X SO SMALL BESSEL K-SUB-XNU+1 OVERFLOWS',5,2)
  Bke(2) = bknu1
  !
  IF ( direct<0. ) CALL D9KNUS(ABS(Xnu+vincr),X,Bke(2),bknu1,iswtch)
  IF ( n==2 ) RETURN
  !
  vend = ABS(Xnu+Nin) - 1.0D0
  IF ( (vend-.5D0)*LOG(vend)+0.27D0-vend*(LOG(X)-.694D0)>alnbig )&
    CALL XERMSG('SLATEC','DBSKES',&
    'X SO SMALL OR ABS(NU) SO BIG THAT BESSEL K-SUB-NU '//&
    'OVERFLOWS',5,2)
  !
  v = Xnu
  DO i = 3 , n
    v = v + vincr
    Bke(i) = 2.0D0*v*Bke(i-1)/X + Bke(i-2)
  ENDDO
  !
END SUBROUTINE DBSKES
