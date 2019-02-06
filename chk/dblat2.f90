!*==DBLAT2.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK DBLAT2
      SUBROUTINE DBLAT2(Nout,Kprint,Ipass)
      IMPLICIT NONE
!*--DBLAT25
!***BEGIN PROLOGUE  DBLAT2
!***PURPOSE  Driver for testing Level 2 BLAS double precision
!            subroutines.
!***LIBRARY   SLATEC (BLAS)
!***CATEGORY  A3B
!***TYPE      DOUBLE PRECISION (SBLAT2-S, DBLAT2-D, CBLAT2-C)
!***KEYWORDS  BLAS, QUICK CHECK DRIVER
!***AUTHOR  Du Croz, J. J., (NAG)
!           Hanson, R. J., (SNLA)
!***DESCRIPTION
!
!  Test program for the DOUBLE           Level 2 Blas.
!
!***REFERENCES  Dongarra, J. J., Du Croz, J. J., Hammarling, S. and
!                 Hanson, R. J.  An  extended  set of Fortran Basic
!                 Linear Algebra Subprograms. ACM TOMS, Vol. 14, No. 1,
!                 pp. 1-17, March 1988.
!***ROUTINES CALLED  DCHK12, DCHK22, DCHK32, DCHK42, DCHK52, DCHK62,
!                    DCHKE2, DMVCH, LDE, R1MACH, XERCLR
!***REVISION HISTORY  (YYMMDD)
!   870810  DATE WRITTEN
!   910619  Modified to meet SLATEC code and prologue standards. (BKS)
!   930315  Removed unused variables.  (WRB)
!   930618  Code modified to improve PASS/FAIL reporting.  (BKS, WRB)
!***END PROLOGUE  DBLAT2
!     .. Parameters ..
      INTEGER NSUBS
      PARAMETER (NSUBS=16)
      DOUBLE PRECISION ZERO , ONE
      PARAMETER (ZERO=0.0D0,ONE=1.0D0)
      INTEGER NMAX , INCMAX
      PARAMETER (NMAX=65,INCMAX=2)
!     .. Scalar Arguments ..
      INTEGER Ipass , Kprint
!     .. Local Scalars ..
      DOUBLE PRECISION eps , err , thresh
      INTEGER i , isnum , j , n , NALF , NBET , NIDIM , NINC , NKB , Nout
      PARAMETER (NIDIM=6,NKB=4,NINC=4,NALF=3,NBET=3)
      LOGICAL same , tsterr , ftl , ftl1 , ftl2
      CHARACTER :: trans
!     .. Local Arrays ..
      DOUBLE PRECISION a(NMAX,NMAX) , aa(NMAX*NMAX) , alf(NALF) , as(NMAX*NMAX)
     &                 , bet(NBET) , g(NMAX) , x(NMAX) , xs(NMAX*INCMAX) ,
     &                 xx(NMAX*INCMAX) , y(NMAX) , ys(NMAX*INCMAX) , yt(NMAX) ,
     &                 yy(NMAX*INCMAX) , z(2*NMAX)
      INTEGER idim(NIDIM) , inc(NINC) , kb(NKB)
      LOGICAL ltest(NSUBS)
      CHARACTER(6) :: snames(NSUBS)
!     .. External Functions ..
      REAL R1MACH
      LOGICAL LDE
      EXTERNAL LDE , R1MACH
!     .. External Subroutines ..
      EXTERNAL DCHK12 , DCHK22 , DCHK32 , DCHK42 , DCHK52 , DCHK62 , DCHKE2 ,
     &         DMVCH
!     .. Intrinsic Functions ..
      INTRINSIC ABS , MAX , MIN
!     .. Data statements ..
      DATA snames/'DGEMV ' , 'DGBMV ' , 'DSYMV ' , 'DSBMV ' , 'DSPMV ' ,
     &     'DTRMV ' , 'DTBMV ' , 'DTPMV ' , 'DTRSV ' , 'DTBSV ' , 'DTPSV ' ,
     &     'DGER  ' , 'DSYR  ' , 'DSPR  ' , 'DSYR2 ' , 'DSPR2 '/
      DATA idim/0 , 1 , 2 , 3 , 5 , 9/
      DATA kb/0 , 1 , 2 , 4/
      DATA inc/1 , 2 , -1 , -2/
      DATA alf/0.0 , 1.0 , 0.7/
      DATA bet/0.0 , 1.0 , 0.9/
!***FIRST EXECUTABLE STATEMENT  DBLAT2
!     Set the flag that indicates whether error exits are to be tested.
      tsterr = .TRUE.
!     Set the threshold value of the test ratio
      thresh = 16.0
!
!     Set IPASS to 1 assuming it will pass.
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
!     Check the reliability of DMVCH using exact data.
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
!     YY holds the exact result. On exit from DMVCH YT holds
!     the result computed by DMVCH.
      trans = 'N'
      ftl = .FALSE.
      CALL DMVCH(trans,n,n,ONE,a,NMAX,x,1,ZERO,y,1,yt,g,yy,eps,err,ftl,Nout,
     &           .TRUE.,Kprint)
      same = LDE(yy,yt,n)
      IF ( .NOT.same.OR.err/=ZERO ) THEN
        Ipass = 0
        IF ( Kprint>=2 ) WRITE (Nout,FMT=99008) trans , same , err
      ENDIF
      trans = 'T'
      ftl = .FALSE.
      CALL DMVCH(trans,n,n,ONE,a,NMAX,x,-1,ZERO,y,-1,yt,g,yy,eps,err,ftl,Nout,
     &           .TRUE.,Kprint)
      same = LDE(yy,yt,n)
      IF ( .NOT.same.OR.err/=ZERO ) THEN
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
          IF ( tsterr ) CALL DCHKE2(isnum,snames(isnum),Nout,Kprint,ftl1)
!           Test computations.
          CALL XERCLR
          ftl2 = .FALSE.
          SELECT CASE (isnum)
          CASE (3,4,5)
!           Test DSYMV, 03, DSBMV, 04, and DSPMV, 05.
            CALL DCHK22(snames(isnum),eps,thresh,Nout,Kprint,ftl2,NIDIM,idim,
     &                  NKB,kb,NALF,alf,NBET,bet,NINC,inc,NMAX,INCMAX,a,aa,as,x,
     &                  xx,xs,y,yy,ys,yt,g)
          CASE (6,7,8,9,10,11)
!           Test DTRMV, 06, DTBMV, 07, DTPMV, 08,
!           DTRSV, 09, DTBSV, 10, and DTPSV, 11.
            CALL DCHK32(snames(isnum),eps,thresh,Nout,Kprint,ftl2,NIDIM,idim,
     &                  NKB,kb,NINC,inc,NMAX,INCMAX,a,aa,as,y,yy,ys,yt,g,z)
          CASE (12)
!           Test DGER, 12.
            CALL DCHK42(snames(isnum),eps,thresh,Nout,Kprint,ftl2,NIDIM,idim,
     &                  NALF,alf,NINC,inc,NMAX,INCMAX,a,aa,as,x,xx,xs,y,yy,ys,
     &                  yt,g,z)
          CASE (13,14)
!           Test DSYR, 13, and DSPR, 14.
            CALL DCHK52(snames(isnum),eps,thresh,Nout,Kprint,ftl2,NIDIM,idim,
     &                  NALF,alf,NINC,inc,NMAX,INCMAX,a,aa,as,x,xx,xs,y,yy,ys,
     &                  yt,g,z)
          CASE (15,16)
!           Test DSYR2, 15, and DSPR2, 16.
            CALL DCHK62(snames(isnum),eps,thresh,Nout,Kprint,ftl2,NIDIM,idim,
     &                  NALF,alf,NINC,inc,NMAX,INCMAX,a,aa,as,x,xx,xs,y,yy,ys,
     &                  yt,g,z)
          CASE DEFAULT
!           Test DGEMV, 01, and DGBMV, 02.
            CALL DCHK12(snames(isnum),eps,thresh,Nout,Kprint,ftl2,NIDIM,idim,
     &                  NKB,kb,NALF,alf,NBET,bet,NINC,inc,NMAX,INCMAX,a,aa,as,x,
     &                  xx,xs,y,yy,ys,yt,g)
          END SELECT
!
          IF ( ftl1.OR.ftl2 ) Ipass = 0
        ENDIF
      ENDDO
      RETURN
!
99001 FORMAT (' ROUTINES PASS COMPUTATIONAL TESTS IF TEST RATIO IS LES',
     &        'S THAN',F8.2)
99002 FORMAT (' TESTS OF THE DOUBLE PRECISION LEVEL 2 BLAS',//' THE F',
     &        'OLLOWING PARAMETER VALUES WILL BE USED:')
99003 FORMAT ('   FOR N              ',9I6)
99004 FORMAT ('   FOR K              ',7I6)
99005 FORMAT ('   FOR INCX AND INCY  ',7I6)
99006 FORMAT ('   FOR ALPHA          ',7F6.1)
99007 FORMAT ('   FOR BETA           ',7F6.1)
99008 FORMAT (' ERROR IN DMVCH -  IN-LINE DOT PRODUCTS ARE BEING EVALU',
     &        'ATED WRONGLY.',/' DMVCH WAS CALLED WITH TRANS = ',A1,
     &        ' AND RETURNED SAME = ',L1,' AND ERR = ',F12.3,'.',
     &        /' THIS MAY BE DUE TO FAULTS IN THE ARITHMETIC OR THE',
     &        ' COMPILER.')
99009 FORMAT (1X,A6,' WAS NOT TESTED')
99010 FORMAT (' ERROR-EXITS WILL NOT BE TESTED')
!
!     End of DBLAT2.
!
      END SUBROUTINE DBLAT2
