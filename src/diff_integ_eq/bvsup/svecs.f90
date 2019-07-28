!** SVECS
SUBROUTINE SVECS(Ncomp,Lnfc,Yhp,Work,Iwork,Inhomo,Iflag)
  !> Subsidiary to BVSUP
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
  !   890921  Realigned order of variables in certain COMMON blocks.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900328  Added TYPE section.  (WRB)
  !   910722  Updated AUTHOR section.  (ALS)
  USE ML, ONLY: indpvt_com, nfcc_com
  !
  INTEGER, INTENT(IN) :: Inhomo, Ncomp
  INTEGER, INTENT(INOUT) :: Iflag, Lnfc, Iwork(*)
  REAL(SP), INTENT(INOUT) :: Work(*), Yhp(:,:)
  !
  INTEGER :: idp, k, kp, niv
  REAL(SP) :: dum
  !* FIRST EXECUTABLE STATEMENT  SVECS
  IF( Lnfc/=1 ) THEN
    niv = Lnfc
    Lnfc = 2*Lnfc
    nfcc_com = 2*nfcc_com
    kp = Lnfc + 2 + nfcc_com
    idp = indpvt_com
    indpvt_com = 0
    CALL MGSBV(Ncomp,Lnfc,Yhp,Ncomp,niv,Iflag,Work(1),Work(kp),Iwork(1),&
      Inhomo,Yhp(:,Lnfc+1),Work(Lnfc+2),dum)
    Lnfc = Lnfc/2
    nfcc_com = nfcc_com/2
    indpvt_com = idp
    IF( Iflag/=0 .OR. niv/=Lnfc ) THEN
      Iflag = 99
      RETURN
    END IF
  END IF
  DO k = 1, Ncomp
    Yhp(k,Lnfc+1) = Yhp(k,nfcc_com+1)
  END DO
  Iflag = 1
  !
END SUBROUTINE SVECS