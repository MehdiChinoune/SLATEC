!*==REORT.f90  processed by SPAG 6.72Dc at 11:01 on  6 Feb 2019
!DECK REORT
SUBROUTINE REORT(Ncomp,Y,Yp,Yhp,Niv,W,S,P,Ip,Stowa,Iflag)
  IMPLICIT NONE
  !*--REORT5
  !*** Start of declarations inserted by SPAG
  REAL AE , C , dnd , dndt , dx , P , PWCnd , PX , RE , S , SDOT , srp , &
    Stowa , TND , TOL , vnorm , W , wcnd , X , XBEg
  REAL XENd , XOP , XOT , XSAv , Y , Yhp , Yp , ypnm
  INTEGER ICOco , Iflag , IGOfx , ijk , INDpvt , INFo , INHomo , INTeg , &
    Ip , ISTkop , IVP , j , k , kk , KNSwot , KOP , l , LOTjp , &
    mflag , MNSwot
  INTEGER MXNon , Ncomp , NCOmpd , NDIsk , NEQ , NEQivp , NFC , NFCc , &
    nfcp , NIC , Niv , NOPg , NPS , NSWot , NTApe , NTP , NUMort , &
    NXPts
  !*** End of declarations inserted by SPAG
  !***BEGIN PROLOGUE  REORT
  !***SUBSIDIARY
  !***PURPOSE  Subsidiary to BVSUP
  !***LIBRARY   SLATEC
  !***TYPE      SINGLE PRECISION (REORT-S, DREORT-D)
  !***AUTHOR  Watts, H. A., (SNLA)
  !***DESCRIPTION
  !
  ! **********************************************************************
  !   INPUT
  ! *********
  !     Y, YP and YHP = homogeneous solution matrix and particular
  !                     solution vector to be orthonormalized.
  !     IFLAG = 1 --  store YHP into Y and YP, test for
  !                   reorthonormalization, orthonormalize if needed,
  !                   save restart data.
  !             2 --  store YHP into Y and YP, reorthonormalization,
  !                   no restarts.
  !                   (preset orthonormalization mode)
  !             3 --  store YHP into Y and YP, reorthonormalization
  !                   (when INHOMO=3 and X=XEND).
  ! **********************************************************************
  !   OUTPUT
  ! *********
  !     Y, YP = orthonormalized solutions.
  !     NIV = number of independent vectors returned from DMGSBV.
  !     IFLAG = 0 --  reorthonormalization was performed.
  !            10 --  solution process must be restarted at the last
  !                   orthonormalization point.
  !            30 --  solutions are linearly dependent, problem must
  !                   be restarted from the beginning.
  !     W, P, IP = orthonormalization information.
  ! **********************************************************************
  !
  !***SEE ALSO  BVSUP
  !***ROUTINES CALLED  MGSBV, SDOT, STOR1, STWAY
  !***COMMON BLOCKS    ML15TO, ML18JR, ML8SZ
  !***REVISION HISTORY  (YYMMDD)
  !   750601  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890831  Modified array declarations.  (WRB)
  !   890921  Realigned order of variables in certain COMMON blocks.
  !           (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900328  Added TYPE section.  (WRB)
  !   910722  Updated AUTHOR section.  (ALS)
  !***END PROLOGUE  REORT
  !
  DIMENSION Y(Ncomp,*) , Yp(*) , W(*) , S(*) , P(*) , Ip(*) , Stowa(*) , &
    Yhp(Ncomp,*)
  !
  ! **********************************************************************
  !
  COMMON /ML8SZ / C , XSAv , IGOfx , INHomo , IVP , NCOmpd , NFC
  COMMON /ML15TO/ PX , PWCnd , TND , X , XBEg , XENd , XOT , XOP , INFo(15)&
    , ISTkop , KNSwot , KOP , LOTjp , MNSwot , NSWot
  COMMON /ML18JR/ AE , RE , TOL , NXPts , NIC , NOPg , MXNon , NDIsk , &
    NTApe , NEQ , INDpvt , INTeg , NPS , NTP , NEQivp , &
    NUMort , NFCc , ICOco
  !
  ! **********************************************************************
  !***FIRST EXECUTABLE STATEMENT  REORT
  nfcp = NFC + 1
  !
  !     CHECK TO SEE IF ORTHONORMALIZATION TEST IS TO BE PERFORMED
  !
  IF ( Iflag==1 ) THEN
    KNSwot = KNSwot + 1
    IF ( KNSwot<NSWot ) THEN
      IF ( (XENd-X)*(X-XOT)<0. ) RETURN
    ENDIF
  ENDIF
  CALL STOR1(Y,Yhp,Yp,Yhp(1,nfcp),1,0,0)
  !
  !     ****************************************
  !
  !     ORTHOGONALIZE THE HOMOGENEOUS SOLUTIONS Y
  !     AND PARTICULAR SOLUTION YP.
  !
  Niv = NFC
  CALL MGSBV(Ncomp,NFC,Y,Ncomp,Niv,mflag,S,P,Ip,INHomo,Yp,W,wcnd)
  !
  !     ****************************************
  !
  !  CHECK FOR LINEAR DEPENDENCE OF THE SOLUTIONS.
  !
  IF ( mflag==0 ) THEN
    !
    !     ****************************************
    !
    IF ( Iflag==1 ) THEN
      !
      !     TEST FOR ORTHONORMALIZATION
      !
      IF ( wcnd>=50.*TOL ) THEN
        DO ijk = 1 , nfcp
          IF ( S(ijk)>1.0E+20 ) GOTO 50
        ENDDO
        !
        !     USE LINEAR EXTRAPOLATION ON LOGARITHMIC VALUES OF THE NORM
        !     DECREMENTS TO DETERMINE NEXT ORTHONORMALIZATION CHECKPOINT.
        !     OTHER CONTROLS ON THE NUMBER OF STEPS TO THE NEXT CHECKPOINT
        !     ARE ADDED FOR SAFETY PURPOSES.
        !
        NSWot = KNSwot
        KNSwot = 0
        LOTjp = 0
        wcnd = LOG10(wcnd)
        IF ( wcnd>TND+3. ) NSWot = 2*NSWot
        IF ( wcnd>=PWCnd ) THEN
          XOT = XENd
        ELSE
          dx = X - PX
          dnd = PWCnd - wcnd
          IF ( dnd>=4 ) NSWot = NSWot/2
          dndt = wcnd - TND
          IF ( ABS(dx*dndt)>dnd*ABS(XENd-X) ) THEN
            XOT = XENd
          ELSE
            XOT = X + dx*dndt/dnd
          ENDIF
        ENDIF
        NSWot = MIN(MNSwot,NSWot)
        PWCnd = wcnd
        PX = X
        RETURN
      ENDIF
    ENDIF
    !
    !     ****************************************
    !
    !     ORTHONORMALIZATION NECESSARY SO WE NORMALIZE THE HOMOGENEOUS
    !     SOLUTION VECTORS AND CHANGE W ACCORDINGLY.
    !
    50     NSWot = 1
    KNSwot = 0
    LOTjp = 1
    kk = 1
    l = 1
    DO k = 1 , NFCc
      srp = SQRT(P(kk))
      IF ( INHomo==1 ) W(k) = srp*W(k)
      vnorm = 1./srp
      P(kk) = vnorm
      kk = kk + NFCc + 1 - k
      IF ( NFC/=NFCc ) THEN
        IF ( l/=k/2 ) CYCLE
      ENDIF
      DO j = 1 , Ncomp
        Y(j,l) = Y(j,l)*vnorm
      ENDDO
      l = l + 1
    ENDDO
    !
    IF ( INHomo==1.AND.NPS/=1 ) THEN
      !
      !     NORMALIZE THE PARTICULAR SOLUTION
      !
      ypnm = SDOT(Ncomp,Yp,1,Yp,1)
      IF ( ypnm==0.0 ) ypnm = 1.0
      ypnm = SQRT(ypnm)
      S(nfcp) = ypnm
      DO j = 1 , Ncomp
        Yp(j) = Yp(j)/ypnm
      ENDDO
      DO j = 1 , NFCc
        W(j) = C*W(j)
      ENDDO
    ENDIF
    !
    IF ( Iflag==1 ) CALL STWAY(Y,Yp,Yhp,0,Stowa)
    Iflag = 0
  ELSE
    IF ( Iflag/=2 ) THEN
      IF ( NSWot>1.OR.LOTjp==0 ) THEN
        !
        !     RETRIEVE DATA FOR A RESTART AT LAST ORTHONORMALIZATION POINT
        !
        CALL STWAY(Y,Yp,Yhp,1,Stowa)
        LOTjp = 1
        NSWot = 1
        KNSwot = 0
        MNSwot = MNSwot/2
        TND = TND + 1.
        Iflag = 10
        RETURN
      ENDIF
    ENDIF
    Iflag = 30
    RETURN
  ENDIF
END SUBROUTINE REORT
