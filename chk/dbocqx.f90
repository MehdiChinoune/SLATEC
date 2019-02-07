!*==DBOCQX.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK DBOCQX
SUBROUTINE DBOCQX(Lun,Kprint,Ipass)
  !***BEGIN PROLOGUE  DBOCQX
  !***PURPOSE  Quick check for DBOCLS.
  !***LIBRARY   SLATEC
  !***TYPE      DOUBLE PRECISION (SBOCQX-S, DBOCQX-D)
  !***KEYWORDS  QUICK CHECK
  !***AUTHOR  (UNKNOWN)
  !***DESCRIPTION
  !
  !     MINIMAL TEST DRIVER FOR DBOCLS, BOUNDED CONSTRAINED LEAST
  !     SQUARES SOLVER.  DELIVERS THE VALUE IPASS=1 IF 8 TESTS WERE
  !     PASSED.  DELIVER THE VALUE IPASS=0 IF ANY ONE OF THEM FAILED.
  !
  !     RUN FOUR BOUNDED LEAST SQUARES PROBLEMS THAT COME FROM THE
  !     DIPLOME WORK OF P. ZIMMERMANN.
  !
  !***ROUTINES CALLED  D1MACH, DBOCLS, DBOLS, DCOPY, DNRM2
  !***REVISION HISTORY  (YYMMDD)
  !   850310  DATE WRITTEN
  !   890618  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   901010  Added PASS/FAIL message.  (RWC)
  !***END PROLOGUE  DBOCQX
  IMPLICIT NONE
  !*--DBOCQX27
  !*** Start of declarations inserted by SPAG
  DOUBLE PRECISION D1MACH , DNRM2 , rnorm , rnormc , sr
  INTEGER i , ib , Ipass , irhs , itest , j , Kprint , Lun , mcon , mdw ,&
    mode , mpass , mrows , ncols
  !*** End of declarations inserted by SPAG
  DOUBLE PRECISION d(6,5) , w(11,11) , bl(5,2) , bu(5,2) , x(30) , rw(55) ,&
    xtrue(9)
  DOUBLE PRECISION c(5,5)
  DOUBLE PRECISION bl1(10) , bu1(10)
  INTEGER ind(10) , iw(20) , iopt(40)
  DOUBLE PRECISION rhs(6,2)
  CHARACTER(4) :: msg
  !
  DATA ((c(i,j),i=1,5),j=1,5)/1.D0 , 10.D0 , 4.D0 , 8.D0 , 1.D0 , 1.D0 ,&
    10.D0 , 2.D0 , -1.D0 , 1.D0 , 1.D0 , -3.D0 , -3.D0 , 2.D0 , 1.D0 ,&
    1.D0 , 5.D0 , 5.D0 , 5.D0 , 1.D0 , 1.D0 , 4.D0 , -1.D0 , -3.D0 ,&
    1.D0/
  DATA ((d(i,j),i=1,6),j=1,5)/ - 74.D0 , 14.D0 , 66.D0 , -12.D0 , 3.D0 ,&
    4.D0 , 80.D0 , -69.D0 , -72.D0 , 66.D0 , 8.D0 , -12.D0 , 18.D0 ,&
    21.D0 , -5.D0 , -30.D0 , -7.D0 , 4.D0 , -11.D0 , 28.D0 , 7.D0 ,&
    -23.D0 , -4.D0 , 4.D0 , -4.D0 , 0.D0 , 1.D0 , 3.D0 , 1.D0 , 0.D0/
  DATA ((bl(i,j),i=1,5),j=1,2)/1.D0 , 0.D0 , -1.D0 , 1.D0 , -4.D0 , -1.D0 ,&
    0.D0 , -3.D0 , 1.D0 , -6.D0/
  DATA ((bu(i,j),i=1,5),j=1,2)/3.D0 , 2.D0 , 1.D0 , 3.D0 , -2.D0 , 3.D0 ,&
    4.D0 , 1.D0 , 5.D0 , -2.D0/
  DATA ((rhs(i,j),i=1,6),j=1,2)/51.D0 , -61.D0 , -56.D0 , 69.D0 , 10.D0 ,&
    -12.D0 , -5.D0 , -9.D0 , 708.D0 , 4165.D0 , -13266.D0 , 8409.D0/
  DATA (xtrue(j),j=1,9)/1.D0 , 2.D0 , -1.D0 , 3.D0 , -4.D0 , 1.D0 , 32.D0 ,&
    30.D0 , 31.D0/
  !***FIRST EXECUTABLE STATEMENT  DBOCQX
  mdw = 11
  mrows = 6
  ncols = 5
  mcon = 4
  iopt(1) = 99
  Ipass = 1
  itest = 0
  !
  IF ( Kprint>=2 ) WRITE (Lun,99001)
  99001 FORMAT (' TEST   IB IRHS             SR')
  !
  DO ib = 1 , 2
    DO irhs = 1 , 2
      !
      !           TRANSFER DATA TO WORKING ARRAY W(*,*).
      !
      DO j = 1 , ncols
        CALL DCOPY(mrows,d(1,j),1,w(1,j),1)
      ENDDO
      !
      CALL DCOPY(mrows,rhs(1,irhs),1,w(1,ncols+1),1)
      !
      !             SET BOUND INDICATOR FLAGS.
      !
      DO j = 1 , ncols
        ind(j) = 3
      ENDDO
      !
      CALL DBOLS(w,mdw,mrows,ncols,bl(1,ib),bu(1,ib),ind,iopt,x,rnorm,mode,&
        rw,iw)
      DO j = 1 , ncols
        x(j) = x(j) - xtrue(j)
      ENDDO
      !
      sr = DNRM2(ncols,x,1)
      mpass = 1
      IF ( sr>10.D2*SQRT(D1MACH(4)) ) mpass = 0
      Ipass = Ipass*mpass
      IF ( Kprint>=2 ) THEN
        msg = 'PASS'
        IF ( mpass==0 ) msg = 'FAIL'
        itest = itest + 1
        WRITE (Lun,99003) itest , ib , irhs , sr , msg
      ENDIF
    ENDDO
  ENDDO
  !
  !     RUN STOER'S PROBLEM FROM 1971 SIAM J. N. ANAL. PAPER.
  !
  DO ib = 1 , 2
    DO irhs = 1 , 2
      CALL DCOPY(11*10,0.D0,0,w,1)
      CALL DCOPY(ncols,bl(1,ib),1,bl1,1)
      CALL DCOPY(ncols,bu(1,ib),1,bu1,1)
      ind(ncols+1) = 2
      ind(ncols+2) = 1
      ind(ncols+3) = 2
      ind(ncols+4) = 3
      bu1(ncols+1) = 5.
      bl1(ncols+2) = 20.
      bu1(ncols+3) = 30.
      bl1(ncols+4) = 11.
      bu1(ncols+4) = 40.
      DO j = 1 , ncols
        CALL DCOPY(mcon,c(1,j),1,w(1,j),1)
        CALL DCOPY(mrows,d(1,j),1,w(mcon+1,j),1)
      ENDDO
      !
      CALL DCOPY(mrows,rhs(1,irhs),1,w(mcon+1,ncols+1),1)
      !
      !           CHECK LENGTHS OF REQD. ARRAYS.
      !
      iopt(01) = 2
      iopt(02) = 11
      iopt(03) = 11
      iopt(04) = 10
      iopt(05) = 30
      iopt(06) = 55
      iopt(07) = 20
      iopt(08) = 40
      iopt(09) = 99
      CALL DBOCLS(w,mdw,mcon,mrows,ncols,bl1,bu1,ind,iopt,x,rnormc,rnorm,&
        mode,rw,iw)
      DO j = 1 , ncols + mcon
        x(j) = x(j) - xtrue(j)
      ENDDO
      !
      sr = DNRM2(ncols+mcon,x,1)
      mpass = 1
      IF ( sr>10.D2*SQRT(D1MACH(4)) ) mpass = 0
      Ipass = Ipass*mpass
      IF ( Kprint>=2 ) THEN
        msg = 'PASS'
        IF ( mpass==0 ) msg = 'FAIL'
        itest = itest + 1
        WRITE (Lun,99003) itest , ib , irhs , sr , msg
      ENDIF
    ENDDO
  ENDDO
  !
  !     HERE THE VALUE OF IPASS=1 SAYS THAT DBOCLS HAS PASSED ITS TESTS.
  !          THE VALUE OF IPASS=0 SAYS THAT DBOCLS HAS NOT PASSED.
  !
  IF ( Kprint>=3 ) WRITE (Lun,&
    '('' IPASS VALUE. (A 1 IS GOOD, 0 IS BAD.)'',I4)')&
    Ipass
  IF ( Kprint>=2.AND.Ipass==0 ) WRITE (Lun,99002)
  !
  99002 FORMAT (' ERROR IN DBOCLS OR DBOLS')
  RETURN
  99003 FORMAT (3I5,1P,E20.6,' TEST ',A,'ED.')
END SUBROUTINE DBOCQX
