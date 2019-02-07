!*==CBLAT2.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK CBLAT2
SUBROUTINE CBLAT2(Nout,Kprint,Ipass)
  IMPLICIT NONE
  !*--CBLAT25
  !***BEGIN PROLOGUE  CBLAT2
  !***PURPOSE  Driver for testing Level 2 BLAS complex subroutines.
  !***LIBRARY   SLATEC (BLAS)
  !***CATEGORY  A4
  !***TYPE      COMPLEX (SBLAT2-S, DBLAT2-D, CBLAT2-C)
  !***KEYWORDS  BLAS, QUICK CHECK DRIVER
  !***AUTHOR  Du Croz, J. J., (NAG)
  !           Hanson, R. J., (SNLA)
  !***DESCRIPTION
  !
  !  Test program for the COMPLEX              Level 2 Blas.
  !
  !***REFERENCES  Dongarra, J. J., Du Croz, J. J., Hammarling, S. and
  !                 Hanson, R. J.  An  extended  set of Fortran Basic
  !                 Linear Algebra Subprograms. ACM TOMS, Vol. 14, No. 1,
  !                 pp. 1-17, March 1988.
  !***ROUTINES CALLED  CCHK12, CCHK22, CCHK32, CCHK42, CCHK52, CCHK62,
  !                    CCHKE2, CMVCH, LCE, R1MACH, XERCLR
  !***REVISION HISTORY  (YYMMDD)
  !   870810  DATE WRITTEN
  !   910619  Modified to meet SLATEC code and prologue standards.  (BKS)
  !   930315  Removed unused variables.  (WRB)
  !   930618  Code modified to improve PASS/FAIL reporting.  (BKS, WRB)
  !***END PROLOGUE  CBLAT2
  !     .. Parameters ..
  INTEGER NSUBS
  PARAMETER (NSUBS=17)
  COMPLEX ZERO , ONE
  PARAMETER (ZERO=(0.0,0.0),ONE=(1.0,0.0))
  REAL RZERO
  PARAMETER (RZERO=0.0)
  INTEGER NMAX , INCMAX
  PARAMETER (NMAX=65,INCMAX=2)
  !     .. Scalar Arguments ..
  INTEGER Ipass , Kprint
  !     .. Local Scalars ..
  REAL eps , err , thresh
  INTEGER i , isnum , j , n , NALF , NBET , NIDIM , NINC , NKB , Nout
  PARAMETER (NIDIM=6,NKB=4,NINC=4,NALF=3,NBET=3)
  LOGICAL same , tsterr , ftl , ftl1 , ftl2
  CHARACTER :: trans
  !     .. Local Arrays ..
  COMPLEX a(NMAX,NMAX) , aa(NMAX*NMAX) , alf(NALF) , as(NMAX*NMAX) ,&
    bet(NBET) , x(NMAX) , xs(NMAX*INCMAX) , xx(NMAX*INCMAX) , y(NMAX)&
    , ys(NMAX*INCMAX) , yt(NMAX) , yy(NMAX*INCMAX) , z(2*NMAX)
  REAL g(NMAX)
  INTEGER idim(NIDIM) , inc(NINC) , kb(NKB)
  LOGICAL ltest(NSUBS)
  CHARACTER(6) :: snames(NSUBS)
  !     .. External Functions ..
  REAL R1MACH
  LOGICAL LCE
  EXTERNAL LCE , R1MACH
  !     .. External Subroutines ..
  EXTERNAL CCHK12 , CCHK22 , CCHK32 , CCHK42 , CCHK52 , CCHK62 , CCHKE2 ,&
    CMVCH
  !     .. Intrinsic Functions ..
  INTRINSIC ABS , MAX , MIN
  !     .. Data statements ..
  DATA snames/'CGEMV ' , 'CGBMV ' , 'CHEMV ' , 'CHBMV ' , 'CHPMV ' ,&
    'CTRMV ' , 'CTBMV ' , 'CTPMV ' , 'CTRSV ' , 'CTBSV ' , 'CTPSV ' ,&
    'CGERC ' , 'CGERU ' , 'CHER  ' , 'CHPR  ' , 'CHER2 ' , 'CHPR2 '/
  DATA idim/0 , 1 , 2 , 3 , 5 , 9/
  DATA kb/0 , 1 , 2 , 4/
  DATA inc/1 , 2 , -1 , -2/
  DATA alf/(0.0,0.0) , (1.0,0.0) , (0.7,-0.9)/
  DATA bet/(0.0,0.0) , (1.0,0.0) , (1.3,-1.1)/
  !***FIRST EXECUTABLE STATEMENT  CBLAT2
  !
  !     Set the flag that indicates whether error exits are to be tested.
  tsterr = .TRUE.
  !     Set the threshold value of the test ratio
  thresh = 16.0
  !
  !     Set IPASS = 1 assuming all tests will pass.
  !
  Ipass = 1
  !
  !     Report values of parameters.
  !
  IF ( Kprint>=3 ) THEN
    WRITE (Nout,FMT=99002)
    WRITE (Nout,FMT=99003) (idim(i),i=1,NIDIM)
    WRITE (Nout,FMT=99004) (kb(i),i=1,NKB)
    WRITE (Nout,FMT=99005) (inc(i),i=1,NINC)
    WRITE (Nout,FMT=99006) (alf(i),i=1,NALF)
    WRITE (Nout,FMT=99007) (bet(i),i=1,NBET)
    IF ( .NOT.tsterr ) WRITE (Nout,FMT=99010)
    WRITE (Nout,FMT=99001) thresh
  ENDIF
  !
  !     Set names of subroutines and flags which indicate
  !     whether they are to be tested.
  !
  DO i = 1 , NSUBS
    ltest(i) = .TRUE.
  ENDDO
  !
  !     Set EPS (the machine precision).
  !
  eps = R1MACH(4)
  !
  !     Check the reliability of CMVCH using exact data.
  !
  n = MIN(32,NMAX)
  DO j = 1 , n
    DO i = 1 , n
      a(i,j) = MAX(i-j+1,0)
    ENDDO
    x(j) = j
    y(j) = ZERO
  ENDDO
  DO j = 1 , n
    yy(j) = j*((j+1)*j)/2 - ((j+1)*j*(j-1))/3
  ENDDO
  !     YY holds the exact result. On exit from CMVCH YT holds
  !     the result computed by CMVCH.
  trans = 'N'
  ftl = .FALSE.
  CALL CMVCH(trans,n,n,ONE,a,NMAX,x,1,ZERO,y,1,yt,g,yy,eps,err,ftl,Nout,&
    .TRUE.,Kprint)
  same = LCE(yy,yt,n)
  IF ( .NOT.same.OR.err/=RZERO ) THEN
    Ipass = 0
    IF ( Kprint>=2 ) WRITE (Nout,FMT=99008) trans , same , err
  ENDIF
  trans = 'T'
  ftl = .FALSE.
  CALL CMVCH(trans,n,n,ONE,a,NMAX,x,-1,ZERO,y,-1,yt,g,yy,eps,err,ftl,Nout,&
    .TRUE.,Kprint)
  same = LCE(yy,yt,n)
  IF ( .NOT.same.OR.err/=RZERO ) THEN
    Ipass = 0
    IF ( Kprint>=2 ) WRITE (Nout,FMT=99008) trans , same , err
  ENDIF
  !
  !     Test each subroutine in turn.
  !
  DO isnum = 1 , NSUBS
    IF ( .NOT.ltest(isnum) ) THEN
      !           Subprogram is not to be tested.
      WRITE (Nout,FMT=99009) snames(isnum)
    ELSE
      !           Test error exits.
      ftl1 = .FALSE.
      IF ( tsterr ) CALL CCHKE2(isnum,snames(isnum),Nout,Kprint,ftl1)
      !           Test computations.
      ftl2 = .FALSE.
      CALL XERCLR
      SELECT CASE (isnum)
        CASE (3,4,5)
          !           Test CHEMV, 03, CHBMV, 04, and CHPMV, 05.
          CALL CCHK22(snames(isnum),eps,thresh,Nout,Kprint,ftl2,NIDIM,idim,&
            NKB,kb,NALF,alf,NBET,bet,NINC,inc,NMAX,INCMAX,a,aa,as,x,&
            xx,xs,y,yy,ys,yt,g)
        CASE (6,7,8,9,10,11)
          !           Test CTRMV, 06, CTBMV, 07, CTPMV, 08,
          !           CTRSV, 09, CTBSV, 10, and CTPSV, 11.
          CALL CCHK32(snames(isnum),eps,thresh,Nout,Kprint,ftl2,NIDIM,idim,&
            NKB,kb,NINC,inc,NMAX,INCMAX,a,aa,as,y,yy,ys,yt,g,z)
        CASE (12,13)
          !           Test CGERC, 12, CGERU, 13.
          CALL CCHK42(snames(isnum),eps,thresh,Nout,Kprint,ftl2,NIDIM,idim,&
            NALF,alf,NINC,inc,NMAX,INCMAX,a,aa,as,x,xx,xs,y,yy,ys,&
            yt,g,z)
        CASE (14,15)
          !           Test CHER, 14, and CHPR, 15.
          CALL CCHK52(snames(isnum),eps,thresh,Nout,Kprint,ftl2,NIDIM,idim,&
            NALF,alf,NINC,inc,NMAX,INCMAX,a,aa,as,x,xx,xs,y,yy,ys,&
            yt,g,z)
        CASE (16,17)
          !           Test CHER2, 16, and CHPR2, 17.
          CALL CCHK62(snames(isnum),eps,thresh,Nout,Kprint,ftl2,NIDIM,idim,&
            NALF,alf,NINC,inc,NMAX,INCMAX,a,aa,as,x,xx,xs,y,yy,ys,&
            yt,g,z)
        CASE DEFAULT
          !           Test CGEMV, 01, and CGBMV, 02.
          CALL CCHK12(snames(isnum),eps,thresh,Nout,Kprint,ftl2,NIDIM,idim,&
            NKB,kb,NALF,alf,NBET,bet,NINC,inc,NMAX,INCMAX,a,aa,as,x,&
            xx,xs,y,yy,ys,yt,g)
      END SELECT
      !
      IF ( ftl1.OR.ftl2 ) Ipass = 0
    ENDIF
  ENDDO
  RETURN
  !
  99001 FORMAT (' ROUTINES PASS COMPUTATIONAL TESTS IF TEST RATIO IS LES',&
    'S THAN',F8.2)
  99002 FORMAT (' TESTS OF THE COMPLEX          LEVEL 2 BLAS',//' THE F',&
    'OLLOWING PARAMETER VALUES WILL BE USED:')
  99003 FORMAT ('   FOR N              ',9I6)
  99004 FORMAT ('   FOR K              ',7I6)
  99005 FORMAT ('   FOR INCX AND INCY  ',7I6)
  99006 FORMAT ('   FOR ALPHA          ',7('(',F4.1,',',F4.1,')  ',:))
  99007 FORMAT ('   FOR BETA           ',7('(',F4.1,',',F4.1,')  ',:))
  99008 FORMAT (' ERROR IN CMVCH -  IN-LINE DOT PRODUCTS ARE BEING EVALU',&
    'ATED WRONGLY.',/' CMVCH WAS CALLED WITH TRANS = ',A1,&
    ' AND RETURNED SAME = ',L1,' AND ERR = ',F12.3,'. ',&
    /'THIS MAY BE DUE TO FAULTS IN THE ARITHMETIC OR THE COMPILER.')
  99009 FORMAT (1X,A6,' WAS NOT TESTED')
  99010 FORMAT (' ERROR-EXITS WILL NOT BE TESTED')
  !
  !     End of CBLAT2.
  !
END SUBROUTINE CBLAT2
