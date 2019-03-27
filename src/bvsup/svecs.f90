!** SVECS
SUBROUTINE SVECS(Ncomp,Lnfc,Yhp,Work,Iwork,Inhomo,Iflag)
  IMPLICIT NONE
  !>
  !***
  !  Subsidiary to BVSUP
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      SINGLE PRECISION (SVECS-S, DVECS-D)
  !***
  ! **Author:**  Watts, H. A., (SNLA)
  !***
  ! **Description:**
  !
  !  This subroutine is used for the special structure of complex valued
  !  problems. MGSBV is called upon to obtain LNFC vectors from an
  !  original set of 2*LNFC independent vectors so that the resulting
  !  LNFC vectors together with their imaginary product or mate vectors
  !  form an independent set.
  !
  !***
  ! **See also:**  BVSUP
  !***
  ! **Routines called:**  MGSBV
  !***
  ! COMMON BLOCKS    ML18JR

  !* REVISION HISTORY  (YYMMDD)
  !   750601  DATE WRITTEN
  !   890921  Realigned order of variables in certain COMMON blocks.
  !           (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900328  Added TYPE section.  (WRB)
  !   910722  Updated AUTHOR section.  (ALS)
  
  INTEGER ICOco, idp, Iflag, INDpvt, Inhomo, INTeg, Iwork(*), k, kp, &
    Lnfc, LNFcc, MXNon, Ncomp, NDIsk, NEQ, NEQivp, NIC, niv, NOPg, NPS
  REAL AE, dum, RE, TOL, Work(*), Yhp(Ncomp,*)
  INTEGER NTApe, NTP, NUMort, NXPts
  COMMON /ML18JR/ AE, RE, TOL, NXPts, NIC, NOPg, MXNon, NDIsk, &
    NTApe, NEQ, INDpvt, INTeg, NPS, NTP, NEQivp, NUMort, LNFcc, ICOco
  !* FIRST EXECUTABLE STATEMENT  SVECS
  IF ( Lnfc/=1 ) THEN
    niv = Lnfc
    Lnfc = 2*Lnfc
    LNFcc = 2*LNFcc
    kp = Lnfc + 2 + LNFcc
    idp = INDpvt
    INDpvt = 0
    CALL MGSBV(Ncomp,Lnfc,Yhp,Ncomp,niv,Iflag,Work(1),Work(kp),Iwork(1),&
      Inhomo,Yhp(1,Lnfc+1),Work(Lnfc+2),dum)
    Lnfc = Lnfc/2
    LNFcc = LNFcc/2
    INDpvt = idp
    IF ( Iflag/=0.OR.niv/=Lnfc ) THEN
      Iflag = 99
      RETURN
    ENDIF
  ENDIF
  DO k = 1, Ncomp
    Yhp(k,Lnfc+1) = Yhp(k,LNFcc+1)
  ENDDO
  Iflag = 1
END SUBROUTINE SVECS
