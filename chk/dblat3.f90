!*==DBLAT3.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK DBLAT3
      SUBROUTINE DBLAT3(Nout,Kprint,Ipass)
      IMPLICIT NONE
!*--DBLAT35
!***BEGIN PROLOGUE  DBLAT3
!***PURPOSE  Driver for testing Level 3 BLAS double precision
!            subroutines.
!***LIBRARY   SLATEC (BLAS)
!***CATEGORY  A3B
!***TYPE      DOUBLE PRECISION (SBLAT3-S, DBLAT3-D, CBLAT3-C)
!***KEYWORDS  BLAS, QUICK CHECK DRIVER
!***AUTHOR  Dongarra, J. J., (ANL)
!           Duff, I., (AERE)
!           Du Croz, J., (NAG)
!           Hammarling, S., (NAG)
!***DESCRIPTION
!
!  Test program for the DOUBLE           Level 3 Blas.
!
!***REFERENCES  Dongarra, J., Du Croz, J., Duff, I., and Hammarling, S.
!                 A set of level 3 basic linear algebra subprograms.
!                 ACM TOMS, Vol. 16, No. 1, pp. 1-17, March 1990.
!***ROUTINES CALLED  DCHK13, DCHK23, DCHK33, DCHK43, DCHK53, DCHKE3,
!                    DMMCH, LDE, R1MACH, XERCLR
!***REVISION HISTORY  (YYMMDD)
!   890208  DATE WRITTEN
!   910619  Modified to meet SLATEC code and prologue standards. (BKS)
!   930315  Removed unused variables.  (WRB)
!   930618  Code modified to improve PASS/FAIL reporting.  (BKS, WRB)
!***END PROLOGUE  DBLAT3
!     .. Parameters ..
      INTEGER NSUBS
      PARAMETER (NSUBS=6)
      DOUBLE PRECISION ZERO , ONE
      PARAMETER (ZERO=0.0D0,ONE=1.0D0)
      INTEGER NMAX
      PARAMETER (NMAX=65)
!     .. Scalar Arguments ..
      INTEGER Ipass , Kprint
!     .. Local Scalars ..
      DOUBLE PRECISION eps , err , thresh
      INTEGER i , isnum , j , n , NALF , NBET , NIDIM , Nout
      PARAMETER (NIDIM=6,NALF=3,NBET=3)
      LOGICAL same , tsterr , ftl , ftl1 , ftl2
      CHARACTER :: transa , transb
!     .. Local Arrays ..
      DOUBLE PRECISION ab(NMAX,2*NMAX) , aa(NMAX*NMAX) , alf(NALF) ,
     &                 as(NMAX*NMAX) , bet(NBET) , g(NMAX) , bb(NMAX*NMAX) ,
     &                 bs(NMAX*NMAX) , c(NMAX,NMAX) , cc(NMAX*NMAX) ,
     &                 cs(NMAX*NMAX) , ct(NMAX) , w(2*NMAX)
      INTEGER idim(NIDIM)
      LOGICAL ltest(NSUBS)
      CHARACTER(6) :: snames(NSUBS)
!     .. External Functions ..
      REAL R1MACH
      LOGICAL LDE
      EXTERNAL LDE , R1MACH
!     .. External Subroutines ..
      EXTERNAL DCHK13 , DCHK23 , DCHK33 , DCHK43 , DCHK53 , DCHKE3 , DMMCH
!     .. Intrinsic Functions ..
      INTRINSIC ABS , MAX , MIN
!     .. Data statements ..
      DATA snames/'DGEMM ' , 'DSYMM ' , 'DTRMM ' , 'DTRSM ' , 'DSYRK ' ,
     &     'DSYR2K'/
      DATA idim/0 , 1 , 2 , 3 , 5 , 9/
      DATA alf/0.0 , 1.0 , 0.7/
      DATA bet/0.0 , 1.0 , 1.3/
!***FIRST EXECUTABLE STATEMENT  DBLAT3
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
        WRITE (Nout,FMT=99004) (alf(i),i=1,NALF)
        WRITE (Nout,FMT=99005) (bet(i),i=1,NBET)
        IF ( .NOT.tsterr ) WRITE (Nout,FMT=99008)
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
!     Check the reliability of DMMCH using exact data.
!
      n = MIN(32,NMAX)
      DO j = 1 , n
        DO i = 1 , n
          ab(i,j) = MAX(i-j+1,0)
        ENDDO
        ab(j,NMAX+1) = j
        ab(1,NMAX+j) = j
        c(j,1) = ZERO
      ENDDO
      DO j = 1 , n
        cc(j) = j*((j+1)*j)/2 - ((j+1)*j*(j-1))/3
      ENDDO
!     CC holds the exact result. On exit from DMMCH CT holds
!     the result computed by DMMCH.
      transa = 'N'
      transb = 'N'
      ftl = .FALSE.
      CALL DMMCH(transa,transb,n,1,n,ONE,ab,NMAX,ab(1,NMAX+1),NMAX,ZERO,c,NMAX,
     &           ct,g,cc,NMAX,eps,err,ftl,Nout,.TRUE.,Kprint)
      same = LDE(cc,ct,n)
      IF ( .NOT.same.OR.err/=ZERO ) THEN
        Ipass = 0
        IF ( Kprint>=2 ) WRITE (Nout,FMT=99006) transa , transb , same , err
      ENDIF
      transb = 'T'
      ftl = .FALSE.
      CALL DMMCH(transa,transb,n,1,n,ONE,ab,NMAX,ab(1,NMAX+1),NMAX,ZERO,c,NMAX,
     &           ct,g,cc,NMAX,eps,err,ftl,Nout,.TRUE.,Kprint)
      same = LDE(cc,ct,n)
      IF ( .NOT.same.OR.err/=ZERO ) THEN
        Ipass = 0
        IF ( Kprint>=2 ) WRITE (Nout,FMT=99006) transa , transb , same , err
      ENDIF
      DO j = 1 , n
        ab(j,NMAX+1) = n - j + 1
        ab(1,NMAX+j) = n - j + 1
      ENDDO
      DO j = 1 , n
        cc(n-j+1) = j*((j+1)*j)/2 - ((j+1)*j*(j-1))/3
      ENDDO
      transa = 'T'
      transb = 'N'
      ftl = .FALSE.
      CALL DMMCH(transa,transb,n,1,n,ONE,ab,NMAX,ab(1,NMAX+1),NMAX,ZERO,c,NMAX,
     &           ct,g,cc,NMAX,eps,err,ftl,Nout,.TRUE.,Kprint)
      same = LDE(cc,ct,n)
      IF ( .NOT.same.OR.err/=ZERO ) THEN
        Ipass = 0
        IF ( Kprint>=2 ) WRITE (Nout,FMT=99006) transa , transb , same , err
      ENDIF
      transb = 'T'
      ftl = .FALSE.
      CALL DMMCH(transa,transb,n,1,n,ONE,ab,NMAX,ab(1,NMAX+1),NMAX,ZERO,c,NMAX,
     &           ct,g,cc,NMAX,eps,err,ftl,Nout,.TRUE.,Kprint)
      same = LDE(cc,ct,n)
      IF ( .NOT.same.OR.err/=ZERO ) THEN
        Ipass = 0
        IF ( Kprint>=2 ) WRITE (Nout,FMT=99006) transa , transb , same , err
      ENDIF
!
!     Test each subroutine in turn.
!
      DO isnum = 1 , NSUBS
        IF ( .NOT.ltest(isnum) ) THEN
!           Subprogram is not to be tested.
          WRITE (Nout,FMT=99007) snames(isnum)
        ELSE
!           Test error exits.
          ftl1 = .FALSE.
          IF ( tsterr ) CALL DCHKE3(isnum,snames(isnum),Nout,Kprint,ftl1)
!           Test computations.
          ftl2 = .FALSE.
          CALL XERCLR
          SELECT CASE (isnum)
          CASE (2)
!           Test DSYMM, 02.
            CALL DCHK23(snames(isnum),eps,thresh,Nout,Kprint,ftl2,NIDIM,idim,
     &                  NALF,alf,NBET,bet,NMAX,ab,aa,as,ab(1,NMAX+1),bb,bs,c,cc,
     &                  cs,ct,g)
          CASE (3,4)
!           Test DTRMM, 03, DTRSM, 04.
            CALL DCHK33(snames(isnum),eps,thresh,Nout,Kprint,ftl2,NIDIM,idim,
     &                  NALF,alf,NMAX,ab,aa,as,ab(1,NMAX+1),bb,bs,ct,g,c)
          CASE (5)
!           Test DSYRK, 05.
            CALL DCHK43(snames(isnum),eps,thresh,Nout,Kprint,ftl2,NIDIM,idim,
     &                  NALF,alf,NBET,bet,NMAX,ab,aa,as,ab(1,NMAX+1),bb,bs,c,cc,
     &                  cs,ct,g)
          CASE (6)
!           Test DSYR2K, 06.
            CALL DCHK53(snames(isnum),eps,thresh,Nout,Kprint,ftl2,NIDIM,idim,
     &                  NALF,alf,NBET,bet,NMAX,ab,aa,as,bb,bs,c,cc,cs,ct,g,w)
          CASE DEFAULT
!           Test DGEMM, 01.
            CALL DCHK13(snames(isnum),eps,thresh,Nout,Kprint,ftl2,NIDIM,idim,
     &                  NALF,alf,NBET,bet,NMAX,ab,aa,as,ab(1,NMAX+1),bb,bs,c,cc,
     &                  cs,ct,g)
          END SELECT
!
          IF ( ftl1.OR.ftl2 ) Ipass = 0
        ENDIF
      ENDDO
      RETURN
!
99001 FORMAT (' ROUTINES PASS COMPUTATIONAL TESTS IF TEST RATIO IS LES',
     &        'S THAN',F8.2)
99002 FORMAT (' TESTS OF THE DOUBLE PRECISION LEVEL 3 BLAS',//' THE F',
     &        'OLLOWING PARAMETER VALUES WILL BE USED:')
99003 FORMAT ('   FOR N              ',9I6)
99004 FORMAT ('   FOR ALPHA          ',7F6.1)
99005 FORMAT ('   FOR BETA           ',7F6.1)
99006 FORMAT (' ERROR IN DMMCH -  IN-LINE DOT PRODUCTS ARE BEING EVALU',
     &        'ATED WRONGLY.',/' DMMCH WAS CALLED WITH TRANSA = ',A1,
     &        ' AND TRANSB = ',A1,/' AND RETURNED SAME = ',L1,' AND ','ERR = ',
     &        F12.3,'.',/' THIS MAY BE DUE TO FAULTS IN THE ',
     &        'ARITHMETIC OR THE COMPILER.')
99007 FORMAT (1X,A6,' WAS NOT TESTED')
99008 FORMAT (' ERROR-EXITS WILL NOT BE TESTED')
!
!     End of DBLAT3.
!
      END SUBROUTINE DBLAT3
