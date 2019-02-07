!*==DSTWAY.f90  processed by SPAG 6.72Dc at 11:01 on  6 Feb 2019
!DECK DSTWAY
SUBROUTINE DSTWAY(U,V,Yhp,Inout,Stowa)
  IMPLICIT NONE
  !*--DSTWAY5
  !***BEGIN PROLOGUE  DSTWAY
  !***SUBSIDIARY
  !***PURPOSE  Subsidiary to DBVSUP
  !***LIBRARY   SLATEC
  !***TYPE      DOUBLE PRECISION (STWAY-S, DSTWAY-D)
  !***AUTHOR  Watts, H. A., (SNLA)
  !***DESCRIPTION
  !
  !  This subroutine stores (recalls) integration data in the event
  !  that a restart is needed (the homogeneous solution vectors become
  !  too dependent to continue).
  !
  !***SEE ALSO  DBVSUP
  !***ROUTINES CALLED  DSTOR1
  !***COMMON BLOCKS    DML15T, DML18J, DML8SZ
  !***REVISION HISTORY  (YYMMDD)
  !   750601  DATE WRITTEN
  !   890831  Modified array declarations.  (WRB)
  !   890921  Realigned order of variables in certain COMMON blocks.
  !           (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900328  Added TYPE section.  (WRB)
  !   910722  Updated AUTHOR section.  (ALS)
  !***END PROLOGUE  DSTWAY
  !
  INTEGER ICOco , IGOfx , INDpvt , INFo , INHomo , Inout , INTeg , ISTkop , &
    IVP , j , k , KNSwot , ko , KOP , ks , ksj , LOTjp , MNSwot , &
    MXNon , NCOmp , NDIsk , NEQ , NEQivp , NFC , NFCc , NIC , NOPg , &
    NPS , NSWot , NTApe , NTP , NUMort , NXPts
  REAL(8) :: AE , C , PWCnd , PX , RE , Stowa(*) , TND , TOL , U(*) , &
    V(*) , X , XBEg , XENd , XOP , XOT , XSAv , Yhp(*)
  !
  COMMON /DML8SZ/ C , XSAv , IGOfx , INHomo , IVP , NCOmp , NFC
  COMMON /DML15T/ PX , PWCnd , TND , X , XBEg , XENd , XOT , XOP , INFo(15)&
    , ISTkop , KNSwot , KOP , LOTjp , MNSwot , NSWot
  COMMON /DML18J/ AE , RE , TOL , NXPts , NIC , NOPg , MXNon , NDIsk , &
    NTApe , NEQ , INDpvt , INTeg , NPS , NTP , NEQivp , &
    NUMort , NFCc , ICOco
  !
  !***FIRST EXECUTABLE STATEMENT  DSTWAY
  IF ( Inout==1 ) THEN
    !
    !        RECALL FROM STOWA ARRAY AND ISTKOP
    !
    ks = NFC*NCOmp
    CALL DSTOR1(Yhp,Stowa,Yhp(ks+1),Stowa(ks+1),1,0,0)
    ks = ks + NCOmp
    IF ( NEQivp>=1 ) THEN
      DO j = 1 , NEQivp
        ksj = ks + j
        Yhp(ksj) = Stowa(ksj)
      ENDDO
    ENDIF
    ks = ks + NEQivp
    X = Stowa(ks+1)
    INFo(1) = 0
    ko = KOP - ISTkop
    KOP = ISTkop
    IF ( NDIsk/=0.AND.ko/=0 ) THEN
      DO k = 1 , ko
        BACKSPACE NTApe
      ENDDO
    ENDIF
  ELSE
    !
    !        SAVE IN STOWA ARRAY AND ISTKOP
    !
    ks = NFC*NCOmp
    CALL DSTOR1(Stowa,U,Stowa(ks+1),V,1,0,0)
    ks = ks + NCOmp
    IF ( NEQivp>=1 ) THEN
      DO j = 1 , NEQivp
        ksj = ks + j
        Stowa(ksj) = Yhp(ksj)
      ENDDO
    ENDIF
    ks = ks + NEQivp
    Stowa(ks+1) = X
    ISTkop = KOP
    IF ( XOP==X ) ISTkop = KOP + 1
  ENDIF
END SUBROUTINE DSTWAY
