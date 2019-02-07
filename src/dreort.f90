!*==DREORT.f90  processed by SPAG 6.72Dc at 11:01 on  6 Feb 2019
!DECK DREORT
SUBROUTINE DREORT(Ncomp,Y,Yp,Yhp,Niv,W,S,P,Ip,Stowa,Iflag)
  IMPLICIT NONE
  !*--DREORT5
  !***BEGIN PROLOGUE  DREORT
  !***SUBSIDIARY
  !***PURPOSE  Subsidiary to DBVSUP
  !***LIBRARY   SLATEC
  !***TYPE      DOUBLE PRECISION (REORT-S, DREORT-D)
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
  !***SEE ALSO  DBVSUP
  !***ROUTINES CALLED  DDOT, DMGSBV, DSTOR1, DSTWAY
  !***COMMON BLOCKS    DML15T, DML18J, DML8SZ
  !***REVISION HISTORY  (YYMMDD)
  !   750601  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890831  Modified array declarations.  (WRB)
  !   890921  Realigned order of variables in certain COMMON blocks.
  !           (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900328  Added TYPE section.  (WRB)
  !   910722  Updated AUTHOR section.  (ALS)
  !***END PROLOGUE  DREORT
  !
  REAL(8) :: DDOT
  INTEGER ICOco , Iflag , IGOfx , ijk , INDpvt , INFo , INHomo , INTeg , &
    Ip(*) , ISTkop , IVP , j , k , kk , KNSwot , KOP , l , LOTjp , &
    mflag , MNSwot , MXNon , Ncomp , NCOmpd , NDIsk , NEQ , NEQivp , &
    NFC , NFCc , nfcp , NIC , Niv , NOPg , NPS , NSWot , NTApe , NTP , &
    NUMort , NXPts
  REAL(8) :: AE , C , dnd , dndt , dx , P(*) , PWCnd , PX , RE , &
    S(*) , srp , Stowa(*) , TND , TOL , vnorm , W(*) , wcnd , &
    X , XBEg , XENd , XOP , XOT , XSAv , Y(Ncomp,*) , &
    Yhp(Ncomp,*) , Yp(*) , ypnm
  !
  !     ******************************************************************
  !
  COMMON /DML8SZ/ C , XSAv , IGOfx , INHomo , IVP , NCOmpd , NFC
  COMMON /DML15T/ PX , PWCnd , TND , X , XBEg , XENd , XOT , XOP , INFo(15)&
    , ISTkop , KNSwot , KOP , LOTjp , MNSwot , NSWot
  COMMON /DML18J/ AE , RE , TOL , NXPts , NIC , NOPg , MXNon , NDIsk , &
    NTApe , NEQ , INDpvt , INTeg , NPS , NTP , NEQivp , &
    NUMort , NFCc , ICOco
  !
  ! **********************************************************************
  !     BEGIN BLOCK PERMITTING ...EXITS TO 210
  !        BEGIN BLOCK PERMITTING ...EXITS TO 10
  !***FIRST EXECUTABLE STATEMENT  DREORT
  nfcp = NFC + 1
  !
  !           CHECK TO SEE IF ORTHONORMALIZATION TEST IS TO BE PERFORMED
  !
  !        ...EXIT
  IF ( Iflag==1 ) THEN
    KNSwot = KNSwot + 1
    !        ...EXIT
    IF ( KNSwot<NSWot ) THEN
      !     ......EXIT
      IF ( (XENd-X)*(X-XOT)<0.0D0 ) GOTO 99999
    ENDIF
  ENDIF
  CALL DSTOR1(Y,Yhp,Yp,Yhp(1,nfcp),1,0,0)
  !
  !        ***************************************************************
  !
  !        ORTHOGONALIZE THE HOMOGENEOUS SOLUTIONS Y
  !        AND PARTICULAR SOLUTION YP.
  !
  Niv = NFC
  CALL DMGSBV(Ncomp,NFC,Y,Ncomp,Niv,mflag,S,P,Ip,INHomo,Yp,W,wcnd)
  !
  !           ************************************************************
  !
  !        CHECK FOR LINEAR DEPENDENCE OF THE SOLUTIONS.
  !
  IF ( mflag==0 ) THEN
    !           BEGIN BLOCK PERMITTING ...EXITS TO 190
    !              BEGIN BLOCK PERMITTING ...EXITS TO 110
    !
    !                 ******************************************************
    !
    !              ...EXIT
    IF ( Iflag==1 ) THEN
      !
      !                 TEST FOR ORTHONORMALIZATION
      !
      !              ...EXIT
      IF ( wcnd>=50.0D0*TOL ) THEN
        DO ijk = 1 , nfcp
          !              ......EXIT
          IF ( S(ijk)>1.0D20 ) GOTO 50
        ENDDO
        !
        !                 USE LINEAR EXTRAPOLATION ON LOGARITHMIC VALUES OF THE
        !                 NORM DECREMENTS TO DETERMINE NEXT ORTHONORMALIZATION
        !                 CHECKPOINT.  OTHER CONTROLS ON THE NUMBER OF STEPS TO
        !                 THE NEXT CHECKPOINT ARE ADDED FOR SAFETY PURPOSES.
        !
        NSWot = KNSwot
        KNSwot = 0
        LOTjp = 0
        wcnd = LOG10(wcnd)
        IF ( wcnd>TND+3.0D0 ) NSWot = 2*NSWot
        IF ( wcnd<PWCnd ) THEN
          dx = X - PX
          dnd = PWCnd - wcnd
          IF ( dnd>=4 ) NSWot = NSWot/2
          dndt = wcnd - TND
          IF ( ABS(dx*dndt)<=dnd*ABS(XENd-X) ) THEN
            XOT = X + dx*dndt/dnd
            NSWot = MIN(MNSwot,NSWot)
            PWCnd = wcnd
            !           ......EXIT
            PX = X
          ELSE
            XOT = XENd
            NSWot = MIN(MNSwot,NSWot)
            PWCnd = wcnd
            PX = X
          ENDIF
        ELSE
          XOT = XENd
          NSWot = MIN(MNSwot,NSWot)
          PWCnd = wcnd
          PX = X
        ENDIF
        GOTO 99999
      ENDIF
    ENDIF
    !
    !              *********************************************************
    !
    !              ORTHONORMALIZATION NECESSARY SO WE NORMALIZE THE
    !              HOMOGENEOUS SOLUTION VECTORS AND CHANGE W ACCORDINGLY.
    !
    50     NSWot = 1
    KNSwot = 0
    LOTjp = 1
    kk = 1
    l = 1
    DO k = 1 , NFCc
      !                 BEGIN BLOCK PERMITTING ...EXITS TO 140
      srp = SQRT(P(kk))
      IF ( INHomo==1 ) W(k) = srp*W(k)
      vnorm = 1.0D0/srp
      P(kk) = vnorm
      kk = kk + NFCc + 1 - k
      IF ( NFC/=NFCc ) THEN
        !                 ......EXIT
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
      !                 NORMALIZE THE PARTICULAR SOLUTION
      !
      ypnm = DDOT(Ncomp,Yp,1,Yp,1)
      IF ( ypnm==0.0D0 ) ypnm = 1.0D0
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
    IF ( Iflag==1 ) CALL DSTWAY(Y,Yp,Yhp,0,Stowa)
    Iflag = 0
    !           BEGIN BLOCK PERMITTING ...EXITS TO 40
  ELSEIF ( Iflag==2 ) THEN
    Iflag = 30
  ELSEIF ( NSWot<=1.AND.LOTjp/=0 ) THEN
    Iflag = 30
  ELSE
    !
    !                    RETRIEVE DATA FOR A RESTART AT LAST
    !                    ORTHONORMALIZATION POINT
    !
    CALL DSTWAY(Y,Yp,Yhp,1,Stowa)
    LOTjp = 1
    NSWot = 1
    KNSwot = 0
    MNSwot = MNSwot/2
    TND = TND + 1.0D0
    !           .........EXIT
    Iflag = 10
  ENDIF
  99999 END SUBROUTINE DREORT
