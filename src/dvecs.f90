!*==DVECS.f90  processed by SPAG 6.72Dc at 11:01 on  6 Feb 2019
!DECK DVECS
      SUBROUTINE DVECS(Ncomp,Lnfc,Yhp,Work,Iwork,Inhomo,Iflag)
      IMPLICIT NONE
!*--DVECS5
!***BEGIN PROLOGUE  DVECS
!***SUBSIDIARY
!***PURPOSE  Subsidiary to DBVSUP
!***LIBRARY   SLATEC
!***TYPE      DOUBLE PRECISION (SVECS-S, DVECS-D)
!***AUTHOR  Watts, H. A., (SNLA)
!***DESCRIPTION
!
!  This subroutine is used for the special structure of COMPLEX*16
!  valued problems. DMGSBV is called upon to obtain LNFC vectors from an
!  original set of 2*LNFC independent vectors so that the resulting
!  LNFC vectors together with their imaginary product or mate vectors
!  form an independent set.
!
!***SEE ALSO  DBVSUP
!***ROUTINES CALLED  DMGSBV
!***COMMON BLOCKS    DML18J
!***REVISION HISTORY  (YYMMDD)
!   750601  DATE WRITTEN
!   890831  Modified array declarations.  (WRB)
!   890921  Realigned order of variables in certain COMMON blocks.
!           (WRB)
!   891009  Removed unreferenced statement label.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900328  Added TYPE section.  (WRB)
!   910722  Updated AUTHOR section.  (ALS)
!***END PROLOGUE  DVECS
!
      INTEGER ICOco , idp , Iflag , INDpvt , Inhomo , INTeg , Iwork(*) , k , 
     &        kp , Lnfc , LNFcc , MXNon , Ncomp , NDIsk , NEQ , NEQivp , NIC , 
     &        niv , NOPg , NPS , NTApe , NTP , NUMort , NXPts
      DOUBLE PRECISION AE , dum , RE , TOL , Work(*) , Yhp(Ncomp,*)
      COMMON /DML18J/ AE , RE , TOL , NXPts , NIC , NOPg , MXNon , NDIsk , 
     &                NTApe , NEQ , INDpvt , INTeg , NPS , NTP , NEQivp , 
     &                NUMort , LNFcc , ICOco
!***FIRST EXECUTABLE STATEMENT  DVECS
      IF ( Lnfc/=1 ) THEN
        niv = Lnfc
        Lnfc = 2*Lnfc
        LNFcc = 2*LNFcc
        kp = Lnfc + 2 + LNFcc
        idp = INDpvt
        INDpvt = 0
        CALL DMGSBV(Ncomp,Lnfc,Yhp,Ncomp,niv,Iflag,Work(1),Work(kp),Iwork(1),
     &              Inhomo,Yhp(1,Lnfc+1),Work(Lnfc+2),dum)
        Lnfc = Lnfc/2
        LNFcc = LNFcc/2
        INDpvt = idp
        IF ( Iflag/=0.OR.niv/=Lnfc ) THEN
          Iflag = 99
        ELSE
          DO k = 1 , Ncomp
            Yhp(k,Lnfc+1) = Yhp(k,LNFcc+1)
          ENDDO
          Iflag = 1
        ENDIF
      ELSE
        DO k = 1 , Ncomp
          Yhp(k,Lnfc+1) = Yhp(k,LNFcc+1)
        ENDDO
        Iflag = 1
      ENDIF
      END SUBROUTINE DVECS
