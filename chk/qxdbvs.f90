!*==QXDBVS.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK QXDBVS
SUBROUTINE QXDBVS(Lun,Kprint,Ipass)
  IMPLICIT NONE
  !*--QXDBVS5
  !*** Start of declarations inserted by SPAG
  INTEGER i , iflag , igofx , Ipass , ipss , j , kont , kount , Kprint , l ,&
    Lun , ncomp , ndiw , ndw , neqivp , nfc , nic , nrowa , nrowb ,&
    nrowy
  INTEGER numort , nxpts
  !*** End of declarations inserted by SPAG
  !***BEGIN PROLOGUE  QXDBVS
  !***PURPOSE  Quick check for DBVSUP.
  !***LIBRARY   SLATEC
  !***TYPE      DOUBLE PRECISION (QXBVSP-S, QXDBVS-D)
  !***AUTHOR  (UNKNOWN)
  !***ROUTINES CALLED  DBVSUP, PASS
  !***COMMON BLOCKS    DSAVEX
  !***REVISION HISTORY  (YYMMDD)
  !   ??????  DATE WRITTEN
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   901014  Made editorial changes and added correct result to
  !           output. (RWC)
  !   910708  Minor modifications in use of KPRINT.  (WRB)
  !***END PROLOGUE  QXDBVS
  INTEGER itmp(9) , iwork(100)
  DOUBLE PRECISION work(1000) , ae , re , XSAve , sve , TERm , tol
  DOUBLE PRECISION y(4,15) , xpts(15) , a(2,4) , alpha(2) , b(2,4) , beta(2)&
    , yans(2,15) , reler , abser
  CHARACTER(4) :: msg
  COMMON /DSAVEX/ XSAve , TERm
  DATA yans(1,1) , yans(2,1) , yans(1,2) , yans(2,2) , yans(1,3) , yans(2,3)&
    , yans(1,4) , yans(2,4) , yans(1,5) , yans(2,5) , yans(1,6) ,&
    yans(2,6) , yans(1,7) , yans(2,7) , yans(1,8) , yans(2,8) , yans(1,9)&
    , yans(2,9) , yans(1,10) , yans(2,10) , yans(1,11) , yans(2,11) ,&
    yans(1,12) , yans(2,12) , yans(1,13) , yans(2,13) , yans(1,14) ,&
    yans(2,14) , yans(1,15) , yans(2,15)/5.000000000D+00 ,&
    -6.888880126D-01 , 8.609248635D+00 , -1.083092311D+00 ,&
    1.674923836D+01 , -2.072210073D+00 , 3.351098494D+01 ,&
    -4.479263780D+00 , 6.601103894D+01 , -8.909222513D+00 ,&
    8.579580988D+01 , -1.098742758D+01 , 1.106536877D+02 ,&
    -1.402469444D+01 , 1.421228220D+02 , -1.742236546D+01 ,&
    1.803383474D+02 , -2.086465851D+01 , 2.017054332D+02 ,&
    -1.990879843D+01 , 2.051622475D+02 , -1.324886978D+01 ,&
    2.059197452D+02 , 1.051529813D+01 , 1.972191446D+02 ,&
    9.320592785D+01 , 1.556894846D+02 , 3.801682434D+02 ,&
    1.818989404D-12 , 1.379853993D+03/
  DATA xpts(1) , xpts(2) , xpts(3) , xpts(4) , xpts(5) , xpts(6) , xpts(7) ,&
    xpts(8) , xpts(9) , xpts(10) , xpts(11) , xpts(12) , xpts(13) ,&
    xpts(14) , xpts(15)/60.0D+00 , 55.0D+00 , 50.0D+00 , 45.0D+00 ,&
    40.0D+00 , 38.0D+00 , 36.0D+00 , 34.0D+00 , 32.0D+00 , 31.0D+00 ,&
    30.8D+00 , 30.6D+00 , 30.4D+00 , 30.2D+00 , 30.0D+00/
  !***FIRST EXECUTABLE STATEMENT  QXDBVS
  IF ( Kprint>=2 ) THEN
    WRITE (Lun,99001)
    !
    99001   FORMAT ('1')
    WRITE (Lun,99002)
    99002   FORMAT (/' DBVSUP QUICK CHECK')
  ENDIF
  !
  !-----INITIALIZE VARIABLES FOR TEST PROBLEM.
  !
  DO i = 1 , 9
    itmp(i) = 0
  ENDDO
  !
  tol = 1.0D-03
  XSAve = 0.0D+00
  nrowy = 4
  ncomp = 2
  nxpts = 15
  a(1,1) = 1.0D+00
  a(1,2) = 0.0D+00
  nrowa = 2
  alpha(1) = 5.0D+00
  nic = 1
  b(1,1) = 1.0D+00
  b(1,2) = 0.0D+00
  nrowb = 2
  beta(1) = 0.0D+00
  nfc = 1
  igofx = 1
  re = 1.0D-05
  ae = 1.0D-05
  ndw = 1000
  ndiw = 100
  neqivp = 0
  Ipass = 1
  !
  DO i = 1 , 15
    iwork(i) = 0
  ENDDO
  !
  CALL DBVSUP(y,nrowy,ncomp,xpts,nxpts,a,nrowa,alpha,nic,b,nrowb,beta,nfc,&
    igofx,re,ae,iflag,work,ndw,iwork,ndiw,neqivp)
  !
  !-----IF IFLAG = 0, WE HAVE A SUCCESSFUL SOLUTION; OTHERWISE, SKIP
  !     THE ARGUMENT CHECKING AND GO TO THE END.
  !
  IF ( iflag/=0 ) THEN
    Ipass = 0
    IF ( Kprint>1 ) WRITE (Lun,99003) iflag
    99003   FORMAT (10X,'IFLAG =',I2)
    GOTO 300
  ENDIF
  !
  !-----CHECK THE ACCURACY OF THE SOLUTION.
  !
  numort = iwork(1)
  DO j = 1 , nxpts
    DO l = 1 , 2
      abser = ABS(yans(l,j)-y(l,j))
      reler = abser/ABS(yans(l,j))
      IF ( reler>tol.AND.abser>tol ) Ipass = 0
    ENDDO
  ENDDO
  !
  !-----CHECK FOR SUPPRESSION OF PRINTING.
  !
  IF ( Kprint==0.OR.(Kprint==1.AND.Ipass==1) ) GOTO 400
  !
  IF ( Kprint/=1.OR.Ipass/=0 ) THEN
    IF ( Kprint>=3.OR.Ipass==0 ) THEN
      WRITE (Lun,99004)
      99004     FORMAT (/' ACCURACY TEST')
      WRITE (Lun,99005) numort
      99005     FORMAT (/' NUMBER OF ORTHONORMALIZATIONS =',I3)
      WRITE (Lun,99006) (work(j),j=1,numort)
      99006     FORMAT (/' ORTHONORMALIZATION POINTS ARE'/(1X,4F10.2))
      WRITE (Lun,99007)
      99007     FORMAT (//20X,'CALCULATION',30X,'TRUE SOLUTION'/2X,'X',14X,'Y',17X,&
        'Y-PRIME',15X,'Y',17X,'Y-PRIME'/)
      DO j = 1 , nxpts
        msg = 'PASS'
        abser = ABS(yans(1,j)-y(1,j))
        reler = abser/ABS(yans(1,j))
        IF ( reler>tol.AND.abser>tol ) msg = 'FAIL'
        abser = ABS(yans(2,j)-y(2,j))
        reler = abser/ABS(yans(2,j))
        IF ( reler>tol.AND.abser>tol ) msg = 'FAIL'
        WRITE (Lun,99008) xpts(j) , y(1,j) , y(2,j) , yans(1,j) , yans(2,j)&
          , msg
        99008       FORMAT (F5.1,4E20.7,5X,A)
      ENDDO
    ENDIF
  ENDIF
  !
  !-----SEND MESSAGE INDICATING PASSAGE OR FAILURE OF TESTS.
  !
  CALL PASS(Lun,1,Ipass)
  !
  !-----ERROR MESSAGE TESTS.
  !
  IF ( Kprint==1 ) GOTO 400
  kont = 1
  WRITE (Lun,99009)
  99009 FORMAT (/' (7) TESTS OF IFLAG VALUES')
  !
  !-----NROWY LESS THAN NCOMP
  !
  kount = 1
  nrowy = 1
  100  DO
  DO i = 1 , 15
    iwork(i) = 0
  ENDDO
  CALL DBVSUP(y,nrowy,ncomp,xpts,nxpts,a,nrowa,alpha,nic,b,nrowb,beta,nfc,&
    igofx,re,ae,iflag,work,ndw,iwork,ndiw,neqivp)
  SELECT CASE (kount)
    CASE (2)
      !
      WRITE (Lun,99013) iflag
      IF ( iflag==-2 ) itmp(kont) = 1
      kont = kont + 1
      !
      !-----RE OR AE NEGATIVE
      !
      kount = 3
      igofx = 1
      re = -1.0D+00
      ae = -2.0D+00
    CASE (3)
      !
      WRITE (Lun,99013) iflag
      IF ( iflag==-2 ) itmp(kont) = 1
      kont = kont + 1
      !
      !-----NROWA LESS THAN NIC
      !
      kount = 4
      re = 1.0D-05
      ae = 1.0D-05
      nrowa = 0
      EXIT
    CASE (4)
      EXIT
    CASE (5)
      GOTO 200
    CASE (6)
      !
      WRITE (Lun,99010) iflag
      99010     FORMAT (/' IFLAG SHOULD BE -1, IFLAG =',I3)
      IF ( iflag==-1 ) itmp(kont) = 1
      kont = kont + 1
      !-----INCORRECT ORDERING OF XPTS
      kount = 7
      ndiw = 100
      sve = xpts(1)
      xpts(1) = xpts(4)
      xpts(4) = sve
    CASE (7)
      !
      WRITE (Lun,99013) iflag
      IF ( iflag==-2 ) itmp(kont) = 1
      GOTO 300
    CASE DEFAULT
      !
      WRITE (Lun,99013) iflag
      IF ( iflag==-2 ) itmp(kont) = 1
      kont = kont + 1
      !
      !-----IGOFX NOT EQUAL TO 0 OR 1
      !
      kount = 2
      nrowy = 2
      igofx = 3
  END SELECT
ENDDO
!
WRITE (Lun,99013) iflag
IF ( iflag==-2 ) itmp(kont) = 1
kont = kont + 1
!-----NROWB LESS THAN NFC
kount = 5
nrowa = 2
nrowb = 0
!
200  WRITE (Lun,99013) iflag
IF ( iflag==-2 ) itmp(kont) = 1
kont = kont + 1
!-----STORAGE ALLOCATION IS INSUFFICIENT
kount = 6
nrowb = 2
ndiw = 17
GOTO 100
!
!-----SEE IF IFLAG TESTS PASSED
!
300  ipss = 1
DO i = 1 , kont
  ipss = ipss*itmp(i)
ENDDO
!
CALL PASS(Lun,2,ipss)
!
!-----SEE IF ALL TESTS PASSED.
!
Ipass = Ipass*ipss
!
400  IF ( Ipass==1.AND.Kprint>1 ) WRITE (Lun,99011)
99011 FORMAT (/' ***************DBVSUP PASSED ALL TESTS***************')
IF ( Ipass==0.AND.Kprint/=0 ) WRITE (Lun,99012)
99012 FORMAT (/' ***************DBVSUP FAILED SOME TESTS**************')
RETURN
99013 FORMAT (/' IFLAG SHOULD BE -2, IFLAG =',I3)
END SUBROUTINE QXDBVS
