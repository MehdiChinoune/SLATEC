!** DBSKES
PURE SUBROUTINE DBSKES(Xnu,X,Nin,Bke)
  !> Compute a sequence of exponentially scaled modified Bessel
  !  functions of the third kind of fractional order.
  !***
  ! **Library:**   SLATEC (FNLIB)
  !***
  ! **Category:**  C10B3
  !***
  ! **Type:**      DOUBLE PRECISION (BESKES-S, DBSKES-D)
  !***
  ! **Keywords:**  EXPONENTIALLY SCALED, FNLIB, FRACTIONAL ORDER,
  !             MODIFIED BESSEL FUNCTION, SEQUENCE OF BESSEL FUNCTIONS,
  !             SPECIAL FUNCTIONS, THIRD KIND
  !***
  ! **Author:**  Fullerton, W., (LANL)
  !***
  ! **Description:**
  !
  ! DBSKES(XNU,X,NIN,BKE) computes a double precision sequence
  ! of exponentially scaled modified Bessel functions
  ! of the third kind of order XNU + I at X, where X > 0,
  ! XNU lies in (-1,1), and I = 0, 1, ..., NIN - 1, if NIN is positive
  ! and I = 0, -1, ..., NIN + 1, if NIN is negative.  On return, the
  ! vector BKE(.) contains the results at X for order starting at XNU.
  ! XNU, X, and BKE are double precision.  NIN is integer.
  !
  !***
  ! **References:**  (NONE)
  !***
  ! **Routines called:**  D1MACH, D9KNUS, XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   770601  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890911  Removed unnecessary intrinsics.  (WRB)
  !   890911  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  USE service, ONLY : D1MACH
  INTEGER, INTENT(IN) :: Nin
  REAL(DP), INTENT(IN) :: Xnu, X
  REAL(DP), INTENT(OUT) :: Bke(Nin)
  INTEGER :: i, iswtch, n
  REAL(DP) :: bknu1, v, vincr, vend, direct
  REAL(DP), PARAMETER :: alnbig = LOG(D1MACH(2))
  !* FIRST EXECUTABLE STATEMENT  DBSKES
  !
  v = ABS(Xnu)
  n = ABS(Nin)
  !
  IF( v>=1._DP ) THEN
    ERROR STOP 'DBSKES : ABS(XNU) MUST BE LT 1'
  ELSEIF( X<=0._DP ) THEN
    ERROR STOP 'DBSKES : X IS LE 0'
  ELSEIF( n==0 ) THEN
    ERROR STOP 'DBSKES : N THE NUMBER IN THE SEQUENCE IS 0'
  END IF
  !
  CALL D9KNUS(v,X,Bke(1),bknu1,iswtch)
  IF( n==1 ) RETURN
  !
  vincr = SIGN(1._SP,REAL(Nin,SP))
  direct = vincr
  IF( Xnu/=0._DP ) direct = vincr*SIGN(1._DP,Xnu)
  IF( iswtch==1 .AND. direct>0. ) THEN
    ERROR STOP 'DBSKES : X SO SMALL BESSEL K-SUB-XNU+1 OVERFLOWS'
  END IF
  Bke(2) = bknu1
  !
  IF( direct<0. ) CALL D9KNUS(ABS(Xnu+vincr),X,Bke(2),bknu1,iswtch)
  IF( n==2 ) RETURN
  !
  vend = ABS(Xnu+Nin) - 1._DP
  IF( (vend-.5_DP)*LOG(vend)+0.27_DP-vend*(LOG(X)-.694_DP)>alnbig ) THEN
    ERROR STOP 'DBSKES : X SO SMALL OR ABS(NU) SO BIG THAT BESSEL K-SUB-NU OVERFLOWS'
  END IF
  !
  v = Xnu
  DO i = 3, n
    v = v + vincr
    Bke(i) = 2._DP*v*Bke(i-1)/X + Bke(i-2)
  END DO
  !
END SUBROUTINE DBSKES