!*==SBLAT3.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK SBLAT3
      SUBROUTINE SBLAT3(Nout,Kprint,Ipass)
      IMPLICIT NONE
!*--SBLAT35
!***BEGIN PROLOGUE  SBLAT3
!***PURPOSE  Driver for testing Level 3 BLAS single precision
!            subroutines.
!***LIBRARY   SLATEC (BLAS)
!***CATEGORY  A3A
!***TYPE      SINGLE PRECISION (SBLAT3-S, DBLAT3-D, CBLAT3-C)
!***KEYWORDS  BLAS, QUICK CHECK DRIVER
!***AUTHOR  Dongarra, J. J., (ANL)
!           Duff, I., (AERE)
!           Du Croz, J., (NAG)
!           Hammarling, S., (NAG)
!***DESCRIPTION
!
!  Test program for the REAL             Level 3 Blas.
!
!***REFERENCES  Dongarra, J., Du Croz, J., Duff, I., and Hammarling, S.
!                 A set of level 3 basic linear algebra subprograms.
!                 ACM TOMS, Vol. 16, No. 1, pp. 1-17, March 1990.
!***ROUTINES CALLED  LSE, R1MACH, SCHK13, SCHK23, SCHK33, SCHK43,
!                    SCHK53, SCHKE3, SMMCH, XERCLR
!***REVISION HISTORY  (YYMMDD)
!   890208  DATE WRITTEN
!   910619  Modified to meet SLATEC code and prologue standards.  (BKS)
!   930315  Removed unused variables.  (WRB)
!   930618  Code modified to improve PASS/FAIL reporting.  (BKS, WRB)
!   930701  Call to SCHKE5 changed to call to SCHKE3.  (BKS)
!***END PROLOGUE  SBLAT3
!     .. Parameters ..
      INTEGER NSUBS
      PARAMETER (NSUBS=6)
      REAL ZERO , ONE
      PARAMETER (ZERO=0.0,ONE=1.0)
      INTEGER NMAX , INCMAX
      PARAMETER (NMAX=65,INCMAX=2)
!     .. Scalar Arguments ..
      INTEGER Ipass , Kprint
!     .. Local Scalars ..
      REAL eps , err , thresh
      INTEGER i , isnum , j , n , NALF , NBET , NIDIM , Nout
      PARAMETER (NIDIM=6,NALF=3,NBET=3)
      LOGICAL same , tsterr , ftl , ftl1 , ftl2
      CHARACTER :: transa , transb
!     .. Local Arrays ..
      REAL ab(NMAX,2*NMAX) , aa(NMAX*NMAX) , alf(NALF) , as(NMAX*NMAX) ,
     &     bet(NBET) , g(NMAX) , bb(NMAX*NMAX) , bs(NMAX*NMAX) , c(NMAX,NMAX) ,
     &     cc(NMAX*NMAX) , cs(NMAX*NMAX) , ct(NMAX) , w(2*NMAX)
      INTEGER idim(NIDIM)
      LOGICAL ltest(NSUBS)
      CHARACTER(6) :: snames(NSUBS)
!     .. External Functions ..
      REAL R1MACH
      LOGICAL LSE
      EXTERNAL LSE , R1MACH
!     .. External Subroutines ..
      EXTERNAL SCHK13 , SCHK23 , SCHK33 , SCHK43 , SCHK53 , SCHKE3 , SMMCH
!     .. Intrinsic Functions ..
      INTRINSIC MAX , MIN
!     .. Data statements ..
      DATA snames/'SGEMM ' , 'SSYMM ' , 'STRMM ' , 'STRSM ' , 'SSYRK ' ,
     &     'SSYR2K'/
      DATA idim/0 , 1 , 2 , 3 , 5 , 9/
      DATA alf/0.0 , 1.0 , 0.7/
      DATA bet/0.0 , 1.0 , 1.3/
!***FIRST EXECUTABLE STATEMENT  SBLAT3
!
!     Set the flag that indicates whether error exits are to be tested.
!
      tsterr = .TRUE.
!
!     Set the threshold value of the test ratio
!
      thresh = 16.0
!
!     Initialize IPASS to 1 assuming everything will pass.
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
!     Check the reliability of SMMCH using exact data.
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
!     CC holds the exact result. On exit from SMMCH CT holds
!     the result computed by SMMCH.
      transa = 'N'
      transb = 'N'
      ftl = .FALSE.
      CALL SMMCH(transa,transb,n,1,n,ONE,ab,NMAX,ab(1,NMAX+1),NMAX,ZERO,c,NMAX,
     &           ct,g,cc,NMAX,eps,err,ftl,Nout,.TRUE.,Kprint)
      same = LSE(cc,ct,n)
      IF ( .NOT.same.OR.err/=ZERO ) THEN
        Ipass = 0
        IF ( Kprint>=2 ) WRITE (Nout,FMT=99006) transa , transb , same , err
      ENDIF
      transb = 'T'
      ftl = .FALSE.
      CALL SMMCH(transa,transb,n,1,n,ONE,ab,NMAX,ab(1,NMAX+1),NMAX,ZERO,c,NMAX,
     &           ct,g,cc,NMAX,eps,err,ftl,Nout,.TRUE.,Kprint)
      same = LSE(cc,ct,n)
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
      CALL SMMCH(transa,transb,n,1,n,ONE,ab,NMAX,ab(1,NMAX+1),NMAX,ZERO,c,NMAX,
     &           ct,g,cc,NMAX,eps,err,ftl,Nout,.TRUE.,Kprint)
      same = LSE(cc,ct,n)
      IF ( .NOT.same.OR.err/=ZERO ) THEN
        Ipass = 0
        IF ( Kprint>=2 ) WRITE (Nout,FMT=99006) transa , transb , same , err
      ENDIF
      transb = 'T'
      ftl = .FALSE.
      CALL SMMCH(transa,transb,n,1,n,ONE,ab,NMAX,ab(1,NMAX+1),NMAX,ZERO,c,NMAX,
     &           ct,g,cc,NMAX,eps,err,ftl,Nout,.TRUE.,Kprint)
      same = LSE(cc,ct,n)
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
          IF ( tsterr ) CALL SCHKE3(isnum,snames(isnum),Nout,Kprint,ftl1)
!           Test computations.
          ftl2 = .FALSE.
          CALL XERCLR
          SELECT CASE (isnum)
          CASE (2)
!           Test SSYMM, 02.
            CALL SCHK23(snames(isnum),eps,thresh,Nout,Kprint,ftl2,NIDIM,idim,
     &                  NALF,alf,NBET,bet,NMAX,ab,aa,as,ab(1,NMAX+1),bb,bs,c,cc,
     &                  cs,ct,g)
          CASE (3,4)
!           Test STRMM, 03, STRSM, 04.
            CALL SCHK33(snames(isnum),eps,thresh,Nout,Kprint,ftl2,NIDIM,idim,
     &                  NALF,alf,NMAX,ab,aa,as,ab(1,NMAX+1),bb,bs,ct,g,c)
          CASE (5)
!           Test SSYRK, 05.
            CALL SCHK43(snames(isnum),eps,thresh,Nout,Kprint,ftl2,NIDIM,idim,
     &                  NALF,alf,NBET,bet,NMAX,ab,aa,as,ab(1,NMAX+1),bb,bs,c,cc,
     &                  cs,ct,g)
          CASE (6)
!           Test SSYR2K, 06.
            CALL SCHK53(snames(isnum),eps,thresh,Nout,Kprint,ftl2,NIDIM,idim,
     &                  NALF,alf,NBET,bet,NMAX,ab,aa,as,bb,bs,c,cc,cs,ct,g,w)
          CASE DEFAULT
!           Test SGEMM, 01.
            CALL SCHK13(snames(isnum),eps,thresh,Nout,Kprint,ftl2,NIDIM,idim,
     &                  NALF,alf,NBET,bet,NMAX,ab,aa,as,ab(1,NMAX+1),bb,bs,c,cc,
     &                  cs,ct,g)
          END SELECT
          IF ( ftl1.OR.ftl2 ) Ipass = 0
        ENDIF
      ENDDO
      RETURN
!
99001 FORMAT (' ROUTINES PASS COMPUTATIONAL TESTS IF TEST RATIO IS LES',
     &        'S THAN',F8.2)
99002 FORMAT (' TESTS OF THE REAL             LEVEL 3 BLAS',//' THE F',
     &        'OLLOWING PARAMETER VALUES WILL BE USED:')
99003 FORMAT ('   FOR N              ',9I6)
99004 FORMAT ('   FOR ALPHA          ',7F6.1)
99005 FORMAT ('   FOR BETA           ',7F6.1)
99006 FORMAT (' ERROR IN SMMCH -  IN-LINE DOT PRODUCTS ARE BEING EVALU',
     &        'ATED WRONGLY.',/' SMMCH WAS CALLED WITH TRANSA = ',A1,
     &        ' AND TRANSB = ',A1,' AND RETURNED SAME = ',L1,' AND ERR = ',
     &        F12.3,'.',/' THIS MAY BE DUE TO FAULTS IN THE ARITHMETIC OR THE',
     &        ' COMPILER.')
99007 FORMAT (1X,A6,' WAS NOT TESTED')
99008 FORMAT (' ERROR-EXITS WILL NOT BE TESTED')
!
!     End of SBLAT3.
!
      END SUBROUTINE SBLAT3
