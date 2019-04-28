MODULE TEST19_MOD
  IMPLICIT NONE

CONTAINS
  !** DBEG
  REAL(8) FUNCTION DBEG(Reset)
    !>
    !  Generate random numbers.
    !***
    ! **Library:**   SLATEC (BLAS)
    !***
    ! **Author:**  Du Croz, J. (NAG)
    !           Hanson, R. J. (SNLA)
    !***
    ! **Description:**
    !
    !  Generates random numbers uniformly distributed between -0.5 and 0.5.
    !
    !  Auxiliary routine for test program for Level 2 Blas.
    !***
    ! **References:**  (NONE)
    !***
    ! **Routines called:**  (NONE)

    !* REVISION HISTORY  (YYMMDD)
    !   870810  DATE WRITTEN
    !   910619  Modified to meet SLATEC code and prologue standards.  (BKS)

    !     .. Scalar Arguments ..
    LOGICAL Reset
    !     .. Local Scalars ..
    INTEGER, SAVE :: i, ic
    INTEGER, PARAMETER :: mi = 891
    !     .. Intrinsic Functions ..
    INTRINSIC REAL
    !* FIRST EXECUTABLE STATEMENT  DBEG
    IF ( Reset ) THEN
      !        Initialize local variables.
      i = 7
      ic = 0
      Reset = .FALSE.
    END IF
    !
    !     The sequence of values of I is bounded between 1 and 999.
    !     If initial I = 1,2,3,6,7 or 9, the period will be 50.
    !     If initial I = 4 or 8, the period will be 25.
    !     If initial I = 5, the period will be 10.
    !     IC is used to break up the period by skipping 1 value of I in 6.
    !
    ic = ic + 1
    DO
      i = i*mi
      i = i - 1000*(i/1000)
      IF ( ic>=5 ) THEN
        ic = 0
        CYCLE
      END IF
      DBEG = REAL(i-500, 8)/1001.0D0
      EXIT
    END DO
    !
    !     End of DBEG.
    !
  END FUNCTION DBEG
  !** DBLAT2
  SUBROUTINE DBLAT2(Nout,Kprint,Ipass)
    !>
    !  Driver for testing Level 2 BLAS double precision
    !            subroutines.
    !***
    ! **Library:**   SLATEC (BLAS)
    !***
    ! **Category:**  A3B
    !***
    ! **Type:**      DOUBLE PRECISION (SBLAT2-S, DBLAT2-D, CBLAT2-C)
    !***
    ! **Keywords:**  BLAS, QUICK CHECK DRIVER
    !***
    ! **Author:**  Du Croz, J. J., (NAG)
    !           Hanson, R. J., (SNLA)
    !***
    ! **Description:**
    !
    !  Test program for the DOUBLE           Level 2 Blas.
    !
    !***
    ! **References:**  Dongarra, J. J., Du Croz, J. J., Hammarling, S. and
    !                 Hanson, R. J.  An  extended  set of Fortran Basic
    !                 Linear Algebra Subprograms. ACM TOMS, Vol. 14, No. 1,
    !                 pp. 1-17, March 1988.
    !***
    ! **Routines called:**  DCHK12, DCHK22, DCHK32, DCHK42, DCHK52, DCHK62,
    !                    DCHKE2, DMVCH, LDE, R1MACH, XERCLR

    !* REVISION HISTORY  (YYMMDD)
    !   870810  DATE WRITTEN
    !   910619  Modified to meet SLATEC code and prologue standards. (BKS)
    !   930315  Removed unused variables.  (WRB)
    !   930618  Code modified to improve PASS/FAIL reporting.  (BKS, WRB)
    USE slatec, ONLY : R1MACH, XERCLR
    !     .. Parameters ..
    INTEGER, PARAMETER :: NSUBS = 16
    REAL(8), PARAMETER :: ZERO = 0.0D0, ONE = 1.0D0
    INTEGER, PARAMETER :: NMAX = 65, INCMAX = 2
    !     .. Scalar Arguments ..
    INTEGER Ipass, Kprint
    !     .. Local Scalars ..
    REAL(8) :: eps, err, thresh
    INTEGER i, isnum, j, n, Nout
    INTEGER, PARAMETER :: NIDIM = 6, NKB = 4, NINC = 4, NALF = 3, NBET = 3
    LOGICAL same, tsterr, ftl, ftl1, ftl2
    CHARACTER :: trans
    !     .. Local Arrays ..
    REAL(8) :: a(NMAX,NMAX), aa(NMAX*NMAX), as(NMAX*NMAX), &
      g(NMAX), x(NMAX), xs(NMAX*INCMAX), xx(NMAX*INCMAX), y(NMAX), ys(NMAX*INCMAX), &
      yt(NMAX), yy(NMAX*INCMAX), z(2*NMAX)
    LOGICAL ltest(NSUBS)
    !     .. Intrinsic Functions ..
    INTRINSIC ABS, MAX, MIN
    !     .. Data statements ..
    CHARACTER(6), PARAMETER :: snames(NSUBS) = [ 'DGEMV ', 'DGBMV ', 'DSYMV ', &
      'DSBMV ', 'DSPMV ', 'DTRMV ', 'DTBMV ', 'DTPMV ', 'DTRSV ', 'DTBSV ', &
      'DTPSV ', 'DGER  ', 'DSYR  ', 'DSPR  ', 'DSYR2 ', 'DSPR2 ' ]
    INTEGER, PARAMETER :: idimm(NIDIM) = [ 0, 1, 2, 3, 5, 9 ]
    INTEGER, PARAMETER :: kb(NKB) = [ 0, 1, 2, 4 ]
    INTEGER, PARAMETER :: inc(NINC) = [ 1, 2, -1, -2 ]
    REAL(8), PARAMETER :: alf(NALF) = [ 0.0, 1.0, 0.7 ]
    REAL(8), PARAMETER :: bet(NBET) = [ 0.0, 1.0, 0.9 ]
    !* FIRST EXECUTABLE STATEMENT  DBLAT2
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
      WRITE (Nout,FMT=99003) (idimm(i),i=1,NIDIM)
      WRITE (Nout,FMT=99004) (kb(i),i=1,NKB)
      WRITE (Nout,FMT=99005) (inc(i),i=1,NINC)
      WRITE (Nout,FMT=99006) (alf(i),i=1,NALF)
      WRITE (Nout,FMT=99007) (bet(i),i=1,NBET)
      IF ( .NOT.tsterr ) WRITE (Nout,FMT=99010)
      WRITE (Nout,FMT=99001) thresh
    END IF
    !
    !     Set names of subroutines and flags which indicate
    !     whether they are to be tested.
    !
    DO i = 1, NSUBS
      ltest(i) = .TRUE.
    END DO
    !
    !     Set EPS (the machine precision).
    !
    eps = R1MACH(4)
    !
    !     Check the reliability of DMVCH using exact data.
    !
    n = MIN(32,NMAX)
    DO j = 1, n
      DO i = 1, n
        a(i,j) = MAX(i-j+1,0)
      END DO
      x(j) = j
      y(j) = ZERO
    END DO
    DO j = 1, n
      yy(j) = j*((j+1)*j)/2 - ((j+1)*j*(j-1))/3
    END DO
    !     YY holds the exact result. On exit from DMVCH YT holds
    !     the result computed by DMVCH.
    trans = 'N'
    ftl = .FALSE.
    CALL DMVCH(trans,n,n,ONE,a,NMAX,x,1,ZERO,y,1,yt,g,yy,eps,err,ftl,Nout,&
      .TRUE.,Kprint)
    same = LDE(yy,yt,n)
    IF ( .NOT.same.OR.err/=ZERO ) THEN
      Ipass = 0
      IF ( Kprint>=2 ) WRITE (Nout,FMT=99008) trans, same, err
    END IF
    trans = 'T'
    ftl = .FALSE.
    CALL DMVCH(trans,n,n,ONE,a,NMAX,x,-1,ZERO,y,-1,yt,g,yy,eps,err,ftl,Nout,&
      .TRUE.,Kprint)
    same = LDE(yy,yt,n)
    IF ( .NOT.same.OR.err/=ZERO ) THEN
      Ipass = 0
      IF ( Kprint>=2 ) WRITE (Nout,FMT=99008) trans, same, err
    END IF
    !
    !     Test each subroutine in turn.
    !
    DO isnum = 1, NSUBS
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
            CALL DCHK22(snames(isnum),eps,thresh,Nout,Kprint,ftl2,NIDIM,idimm,&
              NKB,kb,NALF,alf,NBET,bet,NINC,inc,NMAX,INCMAX,a,aa,as,x,&
              xx,xs,y,yy,ys,yt,g)
          CASE (6,7,8,9,10,11)
            !           Test DTRMV, 06, DTBMV, 07, DTPMV, 08,
            !           DTRSV, 09, DTBSV, 10, and DTPSV, 11.
            CALL DCHK32(snames(isnum),eps,thresh,Nout,Kprint,ftl2,NIDIM,idimm,&
              NKB,kb,NINC,inc,NMAX,INCMAX,a,aa,as,y,yy,ys,yt,g,z)
          CASE (12)
            !           Test DGER, 12.
            CALL DCHK42(snames(isnum),eps,thresh,Nout,Kprint,ftl2,NIDIM,idimm,&
              NALF,alf,NINC,inc,NMAX,INCMAX,a,aa,as,x,xx,xs,y,yy,ys,yt,g,z)
          CASE (13,14)
            !           Test DSYR, 13, and DSPR, 14.
            CALL DCHK52(snames(isnum),eps,thresh,Nout,Kprint,ftl2,NIDIM,idimm,&
              NALF,alf,NINC,inc,NMAX,INCMAX,a,aa,as,x,xx,xs,yt,g,z)
          CASE (15,16)
            !           Test DSYR2, 15, and DSPR2, 16.
            CALL DCHK62(snames(isnum),eps,thresh,Nout,Kprint,ftl2,NIDIM,idimm,&
              NALF,alf,NINC,inc,NMAX,INCMAX,a,aa,as,x,xx,xs,y,yy,ys,yt,g,z)
          CASE DEFAULT
            !           Test DGEMV, 01, and DGBMV, 02.
            CALL DCHK12(snames(isnum),eps,thresh,Nout,Kprint,ftl2,NIDIM,idimm,&
              NKB,kb,NALF,alf,NBET,bet,NINC,inc,NMAX,INCMAX,a,aa,as,x,&
              xx,xs,y,yy,ys,yt,g)
        END SELECT
        !
        IF ( ftl1.OR.ftl2 ) Ipass = 0
      END IF
    END DO
    RETURN
    !
    99001 FORMAT (' ROUTINES PASS COMPUTATIONAL TESTS IF TEST RATIO IS LES',&
      'S THAN',F8.2)
    99002 FORMAT (' TESTS OF THE DOUBLE PRECISION LEVEL 2 BLAS',//' THE F',&
      'OLLOWING PARAMETER VALUES WILL BE USED:')
    99003 FORMAT ('   FOR N              ',9I6)
    99004 FORMAT ('   FOR K              ',7I6)
    99005 FORMAT ('   FOR INCX AND INCY  ',7I6)
    99006 FORMAT ('   FOR ALPHA          ',7F6.1)
    99007 FORMAT ('   FOR BETA           ',7F6.1)
    99008 FORMAT (' ERROR IN DMVCH -  IN-LINE DOT PRODUCTS ARE BEING EVALU',&
      'ATED WRONGLY.',/' DMVCH WAS CALLED WITH TRANS = ',A1,&
      ' AND RETURNED SAME = ',L1,' AND ERR = ',F12.3,'.',&
      /' THIS MAY BE DUE TO FAULTS IN THE ARITHMETIC OR THE COMPILER.')
    99009 FORMAT (1X,A6,' WAS NOT TESTED')
    99010 FORMAT (' ERROR-EXITS WILL NOT BE TESTED')
    !
    !     End of DBLAT2.
    !
  END SUBROUTINE DBLAT2
  !** DBLAT3
  SUBROUTINE DBLAT3(Nout,Kprint,Ipass)
    !>
    !  Driver for testing Level 3 BLAS double precision
    !            subroutines.
    !***
    ! **Library:**   SLATEC (BLAS)
    !***
    ! **Category:**  A3B
    !***
    ! **Type:**      DOUBLE PRECISION (SBLAT3-S, DBLAT3-D, CBLAT3-C)
    !***
    ! **Keywords:**  BLAS, QUICK CHECK DRIVER
    !***
    ! **Author:**  Dongarra, J. J., (ANL)
    !           Duff, I., (AERE)
    !           Du Croz, J., (NAG)
    !           Hammarling, S., (NAG)
    !***
    ! **Description:**
    !
    !  Test program for the DOUBLE           Level 3 Blas.
    !
    !***
    ! **References:**  Dongarra, J., Du Croz, J., Duff, I., and Hammarling, S.
    !                 A set of level 3 basic linear algebra subprograms.
    !                 ACM TOMS, Vol. 16, No. 1, pp. 1-17, March 1990.
    !***
    ! **Routines called:**  DCHK13, DCHK23, DCHK33, DCHK43, DCHK53, DCHKE3,
    !                    DMMCH, LDE, R1MACH, XERCLR

    !* REVISION HISTORY  (YYMMDD)
    !   890208  DATE WRITTEN
    !   910619  Modified to meet SLATEC code and prologue standards. (BKS)
    !   930315  Removed unused variables.  (WRB)
    !   930618  Code modified to improve PASS/FAIL reporting.  (BKS, WRB)
    USE slatec, ONLY : R1MACH, XERCLR
    !     .. Parameters ..
    INTEGER, PARAMETER :: NSUBS = 6
    REAL(8), PARAMETER :: ZERO = 0.0D0, ONE = 1.0D0
    INTEGER, PARAMETER :: NMAX = 65
    !     .. Scalar Arguments ..
    INTEGER Ipass, Kprint
    !     .. Local Scalars ..
    REAL(8) :: eps, err, thresh
    INTEGER i, isnum, j, n, Nout
    INTEGER, PARAMETER :: NIDIM = 6, NALF = 3, NBET = 3
    LOGICAL same, tsterr, ftl, ftl1, ftl2
    CHARACTER :: transa, transb
    !     .. Local Arrays ..
    REAL(8) :: ab(NMAX,2*NMAX), aa(NMAX*NMAX), as(NMAX*NMAX), g(NMAX), bb(NMAX*NMAX), &
      bs(NMAX*NMAX), c(NMAX,NMAX), cc(NMAX*NMAX), cs(NMAX*NMAX), ct(NMAX), w(2*NMAX)
    LOGICAL ltest(NSUBS)
    !     .. Intrinsic Functions ..
    INTRINSIC ABS, MAX, MIN
    !     .. Data statements ..
    CHARACTER(6), PARAMETER :: snames(NSUBS)  = [ 'DGEMM ', 'DSYMM ', 'DTRMM ', &
      'DTRSM ', 'DSYRK ', 'DSYR2K' ]
    INTEGER, PARAMETER :: idimm(NIDIM) = [ 0, 1, 2, 3, 5, 9 ]
    REAL(8), PARAMETER :: alf(NALF) = [ 0.0, 1.0, 0.7 ]
    REAL(8), PARAMETER :: bet(NBET) = [ 0.0, 1.0, 1.3 ]
    !* FIRST EXECUTABLE STATEMENT  DBLAT3
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
      WRITE (Nout,FMT=99003) (idimm(i),i=1,NIDIM)
      WRITE (Nout,FMT=99004) (alf(i),i=1,NALF)
      WRITE (Nout,FMT=99005) (bet(i),i=1,NBET)
      IF ( .NOT.tsterr ) WRITE (Nout,FMT=99008)
      WRITE (Nout,FMT=99001) thresh
    END IF
    !
    !     Set names of subroutines and flags which indicate
    !     whether they are to be tested.
    !
    DO i = 1, NSUBS
      ltest(i) = .TRUE.
    END DO
    !
    !     Set EPS (the machine precision).
    !
    eps = R1MACH(4)
    !
    !     Check the reliability of DMMCH using exact data.
    !
    n = MIN(32,NMAX)
    DO j = 1, n
      DO i = 1, n
        ab(i,j) = MAX(i-j+1,0)
      END DO
      ab(j,NMAX+1) = j
      ab(1,NMAX+j) = j
      c(j,1) = ZERO
    END DO
    DO j = 1, n
      cc(j) = j*((j+1)*j)/2 - ((j+1)*j*(j-1))/3
    END DO
    !     CC holds the exact result. On exit from DMMCH CT holds
    !     the result computed by DMMCH.
    transa = 'N'
    transb = 'N'
    ftl = .FALSE.
    CALL DMMCH(transa,transb,n,1,n,ONE,ab,NMAX,ab(1,NMAX+1),NMAX,ZERO,c,NMAX,&
      ct,g,cc,NMAX,eps,err,ftl,Nout,.TRUE.,Kprint)
    same = LDE(cc,ct,n)
    IF ( .NOT.same.OR.err/=ZERO ) THEN
      Ipass = 0
      IF ( Kprint>=2 ) WRITE (Nout,FMT=99006) transa, transb, same, err
    END IF
    transb = 'T'
    ftl = .FALSE.
    CALL DMMCH(transa,transb,n,1,n,ONE,ab,NMAX,ab(1,NMAX+1),NMAX,ZERO,c,NMAX,&
      ct,g,cc,NMAX,eps,err,ftl,Nout,.TRUE.,Kprint)
    same = LDE(cc,ct,n)
    IF ( .NOT.same.OR.err/=ZERO ) THEN
      Ipass = 0
      IF ( Kprint>=2 ) WRITE (Nout,FMT=99006) transa, transb, same, err
    END IF
    DO j = 1, n
      ab(j,NMAX+1) = n - j + 1
      ab(1,NMAX+j) = n - j + 1
    END DO
    DO j = 1, n
      cc(n-j+1) = j*((j+1)*j)/2 - ((j+1)*j*(j-1))/3
    END DO
    transa = 'T'
    transb = 'N'
    ftl = .FALSE.
    CALL DMMCH(transa,transb,n,1,n,ONE,ab,NMAX,ab(1,NMAX+1),NMAX,ZERO,c,NMAX,&
      ct,g,cc,NMAX,eps,err,ftl,Nout,.TRUE.,Kprint)
    same = LDE(cc,ct,n)
    IF ( .NOT.same.OR.err/=ZERO ) THEN
      Ipass = 0
      IF ( Kprint>=2 ) WRITE (Nout,FMT=99006) transa, transb, same, err
    END IF
    transb = 'T'
    ftl = .FALSE.
    CALL DMMCH(transa,transb,n,1,n,ONE,ab,NMAX,ab(1,NMAX+1),NMAX,ZERO,c,NMAX,&
      ct,g,cc,NMAX,eps,err,ftl,Nout,.TRUE.,Kprint)
    same = LDE(cc,ct,n)
    IF ( .NOT.same.OR.err/=ZERO ) THEN
      Ipass = 0
      IF ( Kprint>=2 ) WRITE (Nout,FMT=99006) transa, transb, same, err
    END IF
    !
    !     Test each subroutine in turn.
    !
    DO isnum = 1, NSUBS
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
            CALL DCHK23(snames(isnum),eps,thresh,Nout,Kprint,ftl2,NIDIM,idimm,&
              NALF,alf,NBET,bet,NMAX,ab,aa,as,ab(1,NMAX+1),bb,bs,c,cc,cs,ct,g)
          CASE (3,4)
            !           Test DTRMM, 03, DTRSM, 04.
            CALL DCHK33(snames(isnum),eps,thresh,Nout,Kprint,ftl2,NIDIM,idimm,&
              NALF,alf,NMAX,ab,aa,as,ab(1,NMAX+1),bb,bs,ct,g,c)
          CASE (5)
            !           Test DSYRK, 05.
            CALL DCHK43(snames(isnum),eps,thresh,Nout,Kprint,ftl2,NIDIM,idimm,&
              NALF,alf,NBET,bet,NMAX,ab,aa,as,c,cc,cs,ct,g)
          CASE (6)
            !           Test DSYR2K, 06.
            CALL DCHK53(snames(isnum),eps,thresh,Nout,Kprint,ftl2,NIDIM,idimm,&
              NALF,alf,NBET,bet,NMAX,ab,aa,as,bb,bs,c,cc,cs,ct,g,w)
          CASE DEFAULT
            !           Test DGEMM, 01.
            CALL DCHK13(snames(isnum),eps,thresh,Nout,Kprint,ftl2,NIDIM,idimm,&
              NALF,alf,NBET,bet,NMAX,ab,aa,as,ab(1,NMAX+1),bb,bs,c,cc,cs,ct,g)
        END SELECT
        !
        IF ( ftl1.OR.ftl2 ) Ipass = 0
      END IF
    END DO
    RETURN
    !
    99001 FORMAT (' ROUTINES PASS COMPUTATIONAL TESTS IF TEST RATIO IS LES',&
      'S THAN',F8.2)
    99002 FORMAT (' TESTS OF THE DOUBLE PRECISION LEVEL 3 BLAS',//' THE F',&
      'OLLOWING PARAMETER VALUES WILL BE USED:')
    99003 FORMAT ('   FOR N              ',9I6)
    99004 FORMAT ('   FOR ALPHA          ',7F6.1)
    99005 FORMAT ('   FOR BETA           ',7F6.1)
    99006 FORMAT (' ERROR IN DMMCH -  IN-LINE DOT PRODUCTS ARE BEING EVALU',&
      'ATED WRONGLY.',/' DMMCH WAS CALLED WITH TRANSA = ',A1,&
      ' AND TRANSB = ',A1,/' AND RETURNED SAME = ',L1,' AND ','ERR = ',&
      F12.3,'.',/' THIS MAY BE DUE TO FAULTS IN THE ',&
      'ARITHMETIC OR THE COMPILER.')
    99007 FORMAT (1X,A6,' WAS NOT TESTED')
    99008 FORMAT (' ERROR-EXITS WILL NOT BE TESTED')
    !
    !     End of DBLAT3.
    !
  END SUBROUTINE DBLAT3
  !** DMAKE2
  SUBROUTINE DMAKE2(Type,Uplo,Diag,M,N,A,Nmax,Aa,Lda,Kl,Ku,Reset,Transl)
    !>
    !  Generate values for an M by N matrix A.
    !***
    ! **Library:**   SLATEC (BLAS)
    !***
    ! **Author:**  Du Croz, J. J., (NAG)
    !           Hanson, R. J., (SNLA)
    !***
    ! **Description:**
    !
    !  Generates values for an M by N matrix A within the bandwidth
    !  defined by KL and KU.
    !  Stores the values in the array AA in the data structure required
    !  by the routine, with unwanted elements set to rogue value.
    !
    !  TYPE is 'GE', 'GB', 'SY', 'SB', 'SP', 'TR', 'TB' OR 'TP'.
    !
    !  Auxiliary routine for test program for Level 2 Blas.
    !***
    ! **References:**  (NONE)
    !***
    ! **Routines called:**  DBEG

    !* REVISION HISTORY  (YYMMDD)
    !   870810  DATE WRITTEN
    !   910620  Modified to meet SLATEC code and prologue standards.  (BKS)

    !     .. Parameters ..
    REAL(8), PARAMETER :: ZERO = 0.0D0, ONE = 1.0D0
    REAL(8), PARAMETER :: ROGUE = -1.0D10
    !     .. Scalar Arguments ..
    REAL(8) :: Transl
    INTEGER Kl, Ku, Lda, M, N, Nmax
    LOGICAL Reset
    CHARACTER :: Diag, Uplo
    CHARACTER(2) :: Type
    !     .. Array Arguments ..
    REAL(8) :: A(Nmax,*), Aa(*)
    !     .. Local Scalars ..
    INTEGER i, i1, i2, i3, ibeg, iend, ioff, j, kk
    LOGICAL gen, lower, sym, tri, unit, upper
    !     .. Intrinsic Functions ..
    INTRINSIC MAX, MIN
    !* FIRST EXECUTABLE STATEMENT  DMAKE2
    gen = Type(1:1)=='G'
    sym = Type(1:1)=='S'
    tri = Type(1:1)=='T'
    upper = (sym.OR.tri) .AND. Uplo=='U'
    lower = (sym.OR.tri) .AND. Uplo=='L'
    unit = tri .AND. Diag=='U'
    !
    !     Generate data in array A.
    !
    DO j = 1, N
      DO i = 1, M
        IF ( gen.OR.(upper.AND.i<=j).OR.(lower.AND.i>=j) ) THEN
          IF ( (i<=j.AND.j-i<=Ku).OR.(i>=j.AND.i-j<=Kl) ) THEN
            A(i,j) = DBEG(Reset) + Transl
          ELSE
            A(i,j) = ZERO
          END IF
          IF ( i/=j ) THEN
            IF ( sym ) THEN
              A(j,i) = A(i,j)
            ELSEIF ( tri ) THEN
              A(j,i) = ZERO
            END IF
          END IF
        END IF
      END DO
      IF ( tri ) A(j,j) = A(j,j) + ONE
      IF ( unit ) A(j,j) = ONE
    END DO
    !
    !     Store elements in array AS in data structure required by routine.
    !
    IF ( Type=='GE' ) THEN
      DO j = 1, N
        DO i = 1, M
          Aa(i+(j-1)*Lda) = A(i,j)
        END DO
        DO i = M + 1, Lda
          Aa(i+(j-1)*Lda) = ROGUE
        END DO
      END DO
    ELSEIF ( Type=='GB' ) THEN
      DO j = 1, N
        DO i1 = 1, Ku + 1 - j
          Aa(i1+(j-1)*Lda) = ROGUE
        END DO
        DO i2 = i1, MIN(Kl+Ku+1,Ku+1+M-j)
          Aa(i2+(j-1)*Lda) = A(i2+j-Ku-1,j)
        END DO
        DO i3 = i2, Lda
          Aa(i3+(j-1)*Lda) = ROGUE
        END DO
      END DO
    ELSEIF ( Type=='SY'.OR.Type=='TR' ) THEN
      DO j = 1, N
        IF ( upper ) THEN
          ibeg = 1
          IF ( unit ) THEN
            iend = j - 1
          ELSE
            iend = j
          END IF
        ELSE
          IF ( unit ) THEN
            ibeg = j + 1
          ELSE
            ibeg = j
          END IF
          iend = N
        END IF
        DO i = 1, ibeg - 1
          Aa(i+(j-1)*Lda) = ROGUE
        END DO
        DO i = ibeg, iend
          Aa(i+(j-1)*Lda) = A(i,j)
        END DO
        DO i = iend + 1, Lda
          Aa(i+(j-1)*Lda) = ROGUE
        END DO
      END DO
    ELSEIF ( Type=='SB'.OR.Type=='TB' ) THEN
      DO j = 1, N
        IF ( upper ) THEN
          kk = Kl + 1
          ibeg = MAX(1,Kl+2-j)
          IF ( unit ) THEN
            iend = Kl
          ELSE
            iend = Kl + 1
          END IF
        ELSE
          kk = 1
          IF ( unit ) THEN
            ibeg = 2
          ELSE
            ibeg = 1
          END IF
          iend = MIN(Kl+1,1+M-j)
        END IF
        DO i = 1, ibeg - 1
          Aa(i+(j-1)*Lda) = ROGUE
        END DO
        DO i = ibeg, iend
          Aa(i+(j-1)*Lda) = A(i+j-kk,j)
        END DO
        DO i = iend + 1, Lda
          Aa(i+(j-1)*Lda) = ROGUE
        END DO
      END DO
    ELSEIF ( Type=='SP'.OR.Type=='TP' ) THEN
      ioff = 0
      DO j = 1, N
        IF ( upper ) THEN
          ibeg = 1
          iend = j
        ELSE
          ibeg = j
          iend = N
        END IF
        DO i = ibeg, iend
          ioff = ioff + 1
          Aa(ioff) = A(i,j)
          IF ( i==j ) THEN
            IF ( unit ) Aa(ioff) = ROGUE
          END IF
        END DO
      END DO
    END IF
    !
    !     End of DMAKE2.
    !
  END SUBROUTINE DMAKE2
  !** DMAKE3
  SUBROUTINE DMAKE3(Type,Uplo,Diag,M,N,A,Nmax,Aa,Lda,Reset,Transl)
    !>
    !  Generate values for an M by N matrix A.
    !***
    ! **Library:**   SLATEC (BLAS)
    !***
    ! **Author:**  Dongarra, J. J., (ANL)
    !           Duff, I., (AERE)
    !           Du Croz, J., (NAG)
    !           Hammarling, S., (NAG)
    !***
    ! **Description:**
    !
    !  Generates values for an M by N matrix A within the bandwidth
    !  defined by KL and KU.
    !  Stores the values in the array AA in the data structure required
    !  by the routine, with unwanted elements set to rogue value.
    !
    !  TYPE is 'GE', 'SY' or 'TR'.
    !
    !  Auxiliary routine for test program for Level 3 Blas.
    !***
    ! **References:**  (NONE)
    !***
    ! **Routines called:**  DBEG

    !* REVISION HISTORY  (YYMMDD)
    !   890208  DATE WRITTEN
    !   910620  Modified to meet SLATEC code and prologue standards.  (BKS)

    !     .. Parameters ..
    REAL(8), PARAMETER :: ZERO = 0.0D0, ONE = 1.0D0
    REAL(8), PARAMETER :: ROGUE = -1.0D10
    !     .. Scalar Arguments ..
    REAL(8) :: Transl
    INTEGER Lda, M, N, Nmax
    LOGICAL Reset
    CHARACTER :: Diag, Uplo
    CHARACTER(2) :: Type
    !     .. Array Arguments ..
    REAL(8) :: A(Nmax,*), Aa(*)
    !     .. Local Scalars ..
    INTEGER i, ibeg, iend, j
    LOGICAL gen, lower, sym, tri, unit, upper
    !     .. Intrinsic Functions ..
    INTRINSIC MAX, MIN
    !* FIRST EXECUTABLE STATEMENT  DMAKE3
    gen = Type=='GE'
    sym = Type=='SY'
    tri = Type=='TR'
    upper = (sym.OR.tri) .AND. Uplo=='U'
    lower = (sym.OR.tri) .AND. Uplo=='L'
    unit = tri .AND. Diag=='U'
    !
    !     Generate data in array A.
    !
    DO j = 1, N
      DO i = 1, M
        IF ( gen.OR.(upper.AND.i<=j).OR.(lower.AND.i>=j) ) THEN
          A(i,j) = DBEG(Reset) + Transl
          IF ( i/=j ) THEN
            !                 Set some elements to zero
            IF ( N>3.AND.j==N/2 ) A(i,j) = ZERO
            IF ( sym ) THEN
              A(j,i) = A(i,j)
            ELSEIF ( tri ) THEN
              A(j,i) = ZERO
            END IF
          END IF
        END IF
      END DO
      IF ( tri ) A(j,j) = A(j,j) + ONE
      IF ( unit ) A(j,j) = ONE
    END DO
    !
    !     Store elements in array AS in data structure required by routine.
    !
    IF ( Type=='GE' ) THEN
      DO j = 1, N
        DO i = 1, M
          Aa(i+(j-1)*Lda) = A(i,j)
        END DO
        DO i = M + 1, Lda
          Aa(i+(j-1)*Lda) = ROGUE
        END DO
      END DO
    ELSEIF ( Type=='SY'.OR.Type=='TR' ) THEN
      DO j = 1, N
        IF ( upper ) THEN
          ibeg = 1
          IF ( unit ) THEN
            iend = j - 1
          ELSE
            iend = j
          END IF
        ELSE
          IF ( unit ) THEN
            ibeg = j + 1
          ELSE
            ibeg = j
          END IF
          iend = N
        END IF
        DO i = 1, ibeg - 1
          Aa(i+(j-1)*Lda) = ROGUE
        END DO
        DO i = ibeg, iend
          Aa(i+(j-1)*Lda) = A(i,j)
        END DO
        DO i = iend + 1, Lda
          Aa(i+(j-1)*Lda) = ROGUE
        END DO
      END DO
    END IF
    !
    !     End of DMAKE3.
    !
  END SUBROUTINE DMAKE3
  !** DMMCH
  SUBROUTINE DMMCH(Transa,Transb,M,N,Kk,Alpha,A,Lda,B,Ldb,Beta,C,Ldc,Ct,G,&
      Cc,Ldcc,Eps,Err,Ftl,Nout,Mv,Kprint)
    !>
    !  Check the results of the computational tests.
    !***
    ! **Library:**   SLATEC (BLAS)
    !***
    ! **Author:**  Dongarra, J. J., (ANL)
    !           Duff, I., (AERE)
    !           Du Croz, J., (NAG)
    !           Hammarling, S., (NAG)
    !***
    ! **Description:**
    !
    !  Checks the results of the computational tests.
    !
    !  Auxiliary routine for test program for Level 3 Blas.
    !***
    ! **References:**  (NONE)
    !***
    ! **Routines called:**  (NONE)

    !* REVISION HISTORY  (YYMMDD)
    !   890208  DATE WRITTEN
    !   910620  Modified to meet SLATEC code and prologue standards.  (BKS)

    !     .. Parameters ..
    REAL(8), PARAMETER :: ZERO = 0.0D0, ONE = 1.0D0
    !     .. Scalar Arguments ..
    LOGICAL Ftl
    REAL(8) :: Alpha, Beta, Eps, Err
    INTEGER Kk, Kprint, Lda, Ldb, Ldc, Ldcc, M, N, Nout
    LOGICAL Mv
    CHARACTER :: Transa, Transb
    !     .. Array Arguments ..
    REAL(8) :: A(Lda,*), B(Ldb,*), C(Ldc,*), Cc(Ldcc,*), Ct(*), G(*)
    !     .. Local Scalars ..
    REAL(8) :: erri
    INTEGER i, j, k
    LOGICAL trana, tranb
    !     .. Intrinsic Functions ..
    INTRINSIC ABS, MAX, SQRT
    !* FIRST EXECUTABLE STATEMENT  DMMCH
    trana = Transa=='T' .OR. Transa=='C'
    tranb = Transb=='T' .OR. Transb=='C'
    !
    !     Compute expected result, one column at a time, in CT using data
    !     in A, B and C.
    !     Compute gauges in G.
    !
    DO j = 1, N
      !
      DO i = 1, M
        Ct(i) = ZERO
        G(i) = ZERO
      END DO
      IF ( .NOT.trana.AND..NOT.tranb ) THEN
        DO k = 1, Kk
          DO i = 1, M
            Ct(i) = Ct(i) + A(i,k)*B(k,j)
            G(i) = G(i) + ABS(A(i,k))*ABS(B(k,j))
          END DO
        END DO
      ELSEIF ( trana.AND..NOT.tranb ) THEN
        DO k = 1, Kk
          DO i = 1, M
            Ct(i) = Ct(i) + A(k,i)*B(k,j)
            G(i) = G(i) + ABS(A(k,i))*ABS(B(k,j))
          END DO
        END DO
      ELSEIF ( .NOT.trana.AND.tranb ) THEN
        DO k = 1, Kk
          DO i = 1, M
            Ct(i) = Ct(i) + A(i,k)*B(j,k)
            G(i) = G(i) + ABS(A(i,k))*ABS(B(j,k))
          END DO
        END DO
      ELSEIF ( trana.AND.tranb ) THEN
        DO k = 1, Kk
          DO i = 1, M
            Ct(i) = Ct(i) + A(k,i)*B(j,k)
            G(i) = G(i) + ABS(A(k,i))*ABS(B(j,k))
          END DO
        END DO
      END IF
      DO i = 1, M
        Ct(i) = Alpha*Ct(i) + Beta*C(i,j)
        G(i) = ABS(Alpha)*G(i) + ABS(Beta)*ABS(C(i,j))
      END DO
      !
      !        Compute the error ratio for this result.
      !
      Err = ZERO
      DO i = 1, M
        erri = ABS(Ct(i)-Cc(i,j))/Eps
        IF ( G(i)/=ZERO ) erri = erri/G(i)
        Err = MAX(Err,erri)
        IF ( Err*SQRT(Eps)>=ONE ) THEN
          Ftl = .TRUE.
          IF ( Kprint>=2 ) THEN
            WRITE (Nout,FMT=99001)
            DO k = 1, M
              IF ( Mv ) THEN
                WRITE (Nout,FMT=99002) k, Ct(k), Cc(k,j)
              ELSE
                WRITE (Nout,FMT=99002) k, Cc(k,j), Ct(k)
              END IF
            END DO
          END IF
        END IF
      END DO
    END DO
    RETURN
    !
    99001 FORMAT (' ******* FATAL ERROR - COMPUTED RESULT IS LESS THAN HAL',&
      'F ACCURATE *******',/'           EXPECTED RESULT   COMPUTED RESULT')
    99002 FORMAT (1X,I7,2G18.6)
    !
    !     End of DMMCH.
    !
  END SUBROUTINE DMMCH
  !** DMVCH
  SUBROUTINE DMVCH(Trans,M,N,Alpha,A,Nmax,X,Incx,Beta,Y,Incy,Yt,G,Yy,Eps,&
      Err,Ftl,Nout,Mv,Kprint)
    !>
    !  Check the results of the computational tests.
    !***
    ! **Library:**   SLATEC (BLAS)
    !***
    ! **Author:**  Du Croz, J. J., (NAG)
    !           Hanson, R. J., (SNLA)
    !***
    ! **Description:**
    !
    !  Checks the results of the computational tests.
    !
    !  Auxiliary routine for test program for Level 2 Blas.
    !***
    ! **References:**  (NONE)
    !***
    ! **Routines called:**  (NONE)

    !* REVISION HISTORY  (YYMMDD)
    !   870810  DATE WRITTEN
    !   910620  Modified to meet SLATEC code and prologue standards.  (BKS)

    !     .. Parameters ..
    REAL(8), PARAMETER :: ZERO = 0.0D0, ONE = 1.0D0
    !     .. Scalar Arguments ..
    REAL(8) :: Alpha, Beta, Eps, Err
    INTEGER Incx, Incy, Kprint, M, N, Nmax, Nout
    LOGICAL Mv, Ftl
    CHARACTER :: Trans
    !     .. Array Arguments ..
    REAL(8) :: A(Nmax,*), G(*), X(*), Y(*), Yt(*), Yy(*)
    !     .. Local Scalars ..
    REAL(8) :: erri
    INTEGER i, incxl, incyl, iy, j, jx, k, kx, ky, ml, nl
    LOGICAL tran
    !     .. Intrinsic Functions ..
    INTRINSIC ABS, MAX, SQRT
    !* FIRST EXECUTABLE STATEMENT  DMVCH
    tran = Trans=='T' .OR. Trans=='C'
    IF ( tran ) THEN
      ml = N
      nl = M
    ELSE
      ml = M
      nl = N
    END IF
    IF ( Incx<0 ) THEN
      kx = nl
      incxl = -1
    ELSE
      kx = 1
      incxl = 1
    END IF
    IF ( Incy<0 ) THEN
      ky = ml
      incyl = -1
    ELSE
      ky = 1
      incyl = 1
    END IF
    !
    !     Compute expected result in YT using data in A, X and Y.
    !     Compute gauges in G.
    !
    iy = ky
    DO i = 1, ml
      Yt(iy) = ZERO
      G(iy) = ZERO
      jx = kx
      IF ( tran ) THEN
        DO j = 1, nl
          Yt(iy) = Yt(iy) + A(j,i)*X(jx)
          G(iy) = G(iy) + ABS(A(j,i)*X(jx))
          jx = jx + incxl
        END DO
      ELSE
        DO j = 1, nl
          Yt(iy) = Yt(iy) + A(i,j)*X(jx)
          G(iy) = G(iy) + ABS(A(i,j)*X(jx))
          jx = jx + incxl
        END DO
      END IF
      Yt(iy) = Alpha*Yt(iy) + Beta*Y(iy)
      G(iy) = ABS(Alpha)*G(iy) + ABS(Beta*Y(iy))
      iy = iy + incyl
    END DO
    !
    !     Compute the error ratio for this result.
    !
    Err = ZERO
    DO i = 1, ml
      erri = ABS(Yt(i)-Yy(1+(i-1)*ABS(Incy)))/Eps
      IF ( G(i)/=ZERO ) erri = erri/G(i)
      Err = MAX(Err,erri)
      IF ( Err*SQRT(Eps)>=ONE ) THEN
        Ftl = .TRUE.
        IF ( Kprint>=2 ) THEN
          WRITE (Nout,FMT=99001)
          DO k = 1, ml
            IF ( Mv ) THEN
              WRITE (Nout,FMT=99002) k, Yt(k), Yy(1+(k-1)*ABS(Incy))
            ELSE
              WRITE (Nout,FMT=99002) k, Yy(1+(k-1)*ABS(Incy)), Yt(k)
            END IF
          END DO
        END IF
      END IF
    END DO
    RETURN
    !
    99001 FORMAT (' ******* FATAL ERROR - COMPUTED RESULT IS LESS THAN HAL',&
      'F ACCURATE *******',/'           EXPECTED RESULT   COMPUTED RESULT')
    99002 FORMAT (1X,I7,2G18.6)
    !
    !     End of DMVCH.
    !
  END SUBROUTINE DMVCH
  !** LDE
  LOGICAL FUNCTION LDE(Ri,Rj,Lr)
    !>
    !  Test if two arrays are identical.
    !***
    ! **Library:**   SLATEC (BLAS)
    !***
    ! **Author:**  Du Croz, J. J., (NAG)
    !           Hanson, R. J., (SNLA)
    !***
    ! **Description:**
    !
    !  Tests if two arrays are identical.
    !
    !  Auxiliary routine for test program for Level 2 Blas.
    !***
    ! **References:**  (NONE)
    !***
    ! **Routines called:**  (NONE)

    !* REVISION HISTORY  (YYMMDD)
    !   870810  DATE WRITTEN
    !   910620  Modified to meet SLATEC code and prologue standards.  (BKS)

    !     .. Scalar Arguments ..
    INTEGER Lr
    !     .. Array Arguments ..
    REAL(8) :: Ri(*), Rj(*)
    !     .. Local Scalars ..
    INTEGER i
    !* FIRST EXECUTABLE STATEMENT  LDE
    LDE = .TRUE.
    DO i = 1, Lr
      IF ( Ri(i)/=Rj(i) ) THEN
        LDE = .FALSE.
        EXIT
      END IF
    END DO
    !
    !     End of LDE.
    !
  END FUNCTION LDE
  !** LDERES
  LOGICAL FUNCTION LDERES(Type,Uplo,M,N,Aa,As,Lda)
    !>
    !  Test if selected elements in two arrays are equal.
    !***
    ! **Library:**   SLATEC (BLAS)
    !***
    ! **Author:**  Du Croz, J. J., (NAG)
    !           Hanson, R. J., (SNLA)
    !***
    ! **Description:**
    !
    !  Tests if selected elements in two arrays are equal.
    !
    !  TYPE is 'GE', 'SY' or 'SP'.
    !
    !  Auxiliary routine for test program for Level 2 Blas.
    !***
    ! **References:**  (NONE)
    !***
    ! **Routines called:**  (NONE)

    !* REVISION HISTORY  (YYMMDD)
    !   870810  DATE WRITTEN
    !   910620  Modified to meet SLATEC code and prologue standards.  (BKS)

    !     .. Scalar Arguments ..
    INTEGER Lda, M, N
    CHARACTER :: Uplo
    CHARACTER(2) :: Type
    !     .. Array Arguments ..
    REAL(8) :: Aa(Lda,*), As(Lda,*)
    !     .. Local Scalars ..
    INTEGER i, ibeg, iend, j
    LOGICAL upper
    !* FIRST EXECUTABLE STATEMENT  LDERES
    upper = Uplo=='U'
    IF ( Type=='GE' ) THEN
      DO j = 1, N
        DO i = M + 1, Lda
          IF ( Aa(i,j)/=As(i,j) ) GOTO 100
        END DO
      END DO
    ELSEIF ( Type=='SY' ) THEN
      DO j = 1, N
        IF ( upper ) THEN
          ibeg = 1
          iend = j
        ELSE
          ibeg = j
          iend = N
        END IF
        DO i = 1, ibeg - 1
          IF ( Aa(i,j)/=As(i,j) ) GOTO 100
        END DO
        DO i = iend + 1, Lda
          IF ( Aa(i,j)/=As(i,j) ) GOTO 100
        END DO
      END DO
    END IF
    !
    LDERES = .TRUE.
    RETURN
    100  LDERES = .FALSE.
    !
    !     End of LDERES.
    !
    RETURN
  END FUNCTION LDERES
  !** DCHK12
  SUBROUTINE DCHK12(Sname,Eps,Thresh,Nout,Kprint,Fatal,Nidim,Idimm,Nkb,Kb,&
      Nalf,Alf,Nbet,Bet,Ninc,Inc,Nmax,Incmax,A,Aa,As,X,Xx,Xs,Y,Yy,Ys,Yt,G)
    !>
    !  Test DGEMV and DGBMV.
    !***
    ! **Library:**   SLATEC (BLAS)
    !***
    ! **Keywords:**  BLAS, QUICK CHECK SERVICE ROUTINE
    !***
    ! **Author:**  Du Croz, J. (NAG)
    !           Hanson, R. J. (SNLA)
    !***
    ! **Description:**
    !
    !  Quick check for DGEMV and DGBMV.
    !
    !  Auxiliary routine for test program for Level 2 Blas.
    !***
    ! **References:**  (NONE)
    !***
    ! **Routines called:**  DGBMV, DGEMV, DMAKE2, DMVCH, LDE, LDERES, NUMXER

    !* REVISION HISTORY  (YYMMDD)
    !   870810  DATE WRITTEN
    !   910619  Modified to meet SLATEC code and prologue standards. (BKS)
    USE slatec, ONLY : DGBMV, DGEMV, NUMXER
    !     .. Parameters ..
    REAL(8), PARAMETER :: ZERO = 0.0D0, HALF = 0.5D0
    !     .. Scalar Arguments ..
    LOGICAL Fatal
    REAL(8) :: Eps, Thresh
    INTEGER Incmax, Kprint, Nalf, Nbet, Nidim, Ninc, Nkb, Nmax, Nout
    CHARACTER(6) :: Sname
    !     .. Array Arguments ..
    REAL(8) :: A(Nmax,Nmax), Aa(Nmax*Nmax), Alf(Nalf), As(Nmax*Nmax), &
      Bet(Nbet), G(Nmax), X(Nmax), Xs(Nmax*Incmax), &
      Xx(Nmax*Incmax), Y(Nmax), Ys(Nmax*Incmax), Yt(Nmax), Yy(Nmax*Incmax)
    INTEGER Idimm(Nidim), Inc(Ninc), Kb(Nkb)
    !     .. Local Scalars ..
    REAL(8) :: alpha, als, beta, bls, err, errmax, transl
    INTEGER i, ia, ib, ic, iku, im, in, incx, incxs, incy, incys, &
      ix, iy, kl, kls, ku, kus, laa, lda, ldas, lx, ly, m, &
      ml, ms, n, nargs, nc, nd, nk, nl, ns, nerr
    LOGICAL banded, ftl, full, nul, reset, tran
    CHARACTER :: trans, transs
    !     .. Local Arrays ..
    LOGICAL isame(13)
    !     .. Intrinsic Functions ..
    INTRINSIC ABS, MAX, MIN
    !     .. Data statements ..
    CHARACTER(3), PARAMETER :: ich = 'NTC'
    !* FIRST EXECUTABLE STATEMENT  DCHK12
    full = Sname(3:3)=='E'
    banded = Sname(3:3)=='B'
    !     Define the number of arguments.
    IF ( full ) THEN
      nargs = 11
    ELSEIF ( banded ) THEN
      nargs = 13
    END IF
    !
    nc = 0
    reset = .TRUE.
    errmax = ZERO
    !
    DO in = 1, Nidim
      n = Idimm(in)
      nd = n/2 + 1
      !
      DO im = 1, 2
        IF ( im==1 ) m = MAX(n-nd,0)
        IF ( im==2 ) m = MIN(n+nd,Nmax)
        !
        IF ( banded ) THEN
          nk = Nkb
        ELSE
          nk = 1
        END IF
        DO iku = 1, nk
          IF ( banded ) THEN
            ku = Kb(iku)
            kl = MAX(ku-1,0)
          ELSE
            ku = n - 1
            kl = m - 1
          END IF
          !              Set LDA to 1 more than minimum value if room.
          IF ( banded ) THEN
            lda = kl + ku + 1
          ELSE
            lda = m
          END IF
          IF ( lda<Nmax ) lda = lda + 1
          !              Skip tests if not enough room.
          IF ( lda<=Nmax ) THEN
            laa = lda*n
            nul = n<=0 .OR. m<=0
            !
            !              Generate the matrix A.
            !
            transl = ZERO
            CALL DMAKE2(Sname(2:3),' ',' ',m,n,A,Nmax,Aa,lda,kl,ku,reset,transl)
            !
            DO ic = 1, 3
              trans = ich(ic:ic)
              tran = trans=='T' .OR. trans=='C'
              !
              IF ( tran ) THEN
                ml = n
                nl = m
              ELSE
                ml = m
                nl = n
              END IF
              !
              DO ix = 1, Ninc
                incx = Inc(ix)
                lx = ABS(incx)*nl
                !
                !                    Generate the vector X.
                !
                transl = HALF
                CALL DMAKE2('GE',' ',' ',1,nl,X,1,Xx,ABS(incx),0,nl-1,reset,&
                  transl)
                IF ( nl>1 ) THEN
                  X(nl/2) = ZERO
                  Xx(1+ABS(incx)*(nl/2-1)) = ZERO
                END IF
                !
                DO iy = 1, Ninc
                  incy = Inc(iy)
                  ly = ABS(incy)*ml
                  !
                  DO ia = 1, Nalf
                    alpha = Alf(ia)
                    !
                    DO ib = 1, Nbet
                      beta = Bet(ib)
                      !
                      !                             Generate the vector Y.
                      !
                      transl = ZERO
                      CALL DMAKE2('GE',' ',' ',1,ml,Y,1,Yy,ABS(incy),0,ml-1,&
                        reset,transl)
                      !
                      nc = nc + 1
                      !
                      !                             Save every datum before calling the
                      !                             subroutine.
                      !
                      transs = trans
                      ms = m
                      ns = n
                      kls = kl
                      kus = ku
                      als = alpha
                      DO i = 1, laa
                        As(i) = Aa(i)
                      END DO
                      ldas = lda
                      DO i = 1, lx
                        Xs(i) = Xx(i)
                      END DO
                      incxs = incx
                      bls = beta
                      DO i = 1, ly
                        Ys(i) = Yy(i)
                      END DO
                      incys = incy
                      !
                      !                             Call the subroutine.
                      !
                      IF ( full ) THEN
                        CALL DGEMV(trans,m,n,alpha,Aa,lda,Xx,incx,beta,Yy,incy)
                      ELSEIF ( banded ) THEN
                        CALL DGBMV(trans,m,n,kl,ku,alpha,Aa,lda,Xx,incx,beta,&
                          Yy,incy)
                      END IF
                      !
                      !                             Check if error-exit was taken incorrectly.
                      !
                      IF ( NUMXER(nerr)/=0 ) THEN
                        IF ( Kprint>=2 ) WRITE (Nout,FMT=99007)
                        Fatal = .TRUE.
                      END IF
                      !
                      !                             See what data changed inside subroutines.
                      !
                      isame(1) = trans==transs
                      isame(2) = ms==m
                      isame(3) = ns==n
                      IF ( full ) THEN
                        isame(4) = als==alpha
                        isame(5) = LDE(As,Aa,laa)
                        isame(6) = ldas==lda
                        isame(7) = LDE(Xs,Xx,lx)
                        isame(8) = incxs==incx
                        isame(9) = bls==beta
                        IF ( nul ) THEN
                          isame(10) = LDE(Ys,Yy,ly)
                        ELSE
                          isame(10) = LDERES('GE',' ',1,ml,Ys,Yy,ABS(incy))
                        END IF
                        isame(11) = incys==incy
                      ELSEIF ( banded ) THEN
                        isame(4) = kls==kl
                        isame(5) = kus==ku
                        isame(6) = als==alpha
                        isame(7) = LDE(As,Aa,laa)
                        isame(8) = ldas==lda
                        isame(9) = LDE(Xs,Xx,lx)
                        isame(10) = incxs==incx
                        isame(11) = bls==beta
                        IF ( nul ) THEN
                          isame(12) = LDE(Ys,Yy,ly)
                        ELSE
                          isame(12) = LDERES('GE',' ',1,ml,Ys,Yy,ABS(incy))
                        END IF
                        isame(13) = incys==incy
                      END IF
                      !
                      !                             If data was incorrectly changed, report
                      !                             and return.
                      !
                      DO i = 1, nargs
                        IF ( .NOT.isame(i) ) THEN
                          Fatal = .TRUE.
                          IF ( Kprint>=2 ) WRITE (Nout,FMT=99002) i
                        END IF
                      END DO
                      !
                      ftl = .FALSE.
                      IF ( .NOT.nul ) THEN
                        !
                        !                                Check the result.
                        !
                        CALL DMVCH(trans,m,n,alpha,A,Nmax,X,incx,beta,Y,incy,&
                          Yt,G,Yy,Eps,err,ftl,Nout,.TRUE.,Kprint)
                        errmax = MAX(errmax,err)
                      END IF
                      IF ( ftl ) THEN
                        Fatal = .TRUE.
                        IF ( Kprint>=3 ) THEN
                          WRITE (Nout,FMT=99004) Sname
                          IF ( full ) THEN
                            WRITE (Nout,FMT=99006) nc, Sname, trans, m, &
                              n, alpha, lda, incx, beta, incy
                          ELSEIF ( banded ) THEN
                            WRITE (Nout,FMT=99005) nc, Sname, trans, m, &
                              n, kl, ku, alpha, lda, incx, beta, incy
                          END IF
                        END IF
                      END IF
                      !
                    END DO
                    !
                  END DO
                  !
                END DO
                !
              END DO
              !
            END DO
          END IF
          !
        END DO
        !
      END DO
      !
    END DO
    !
    !     Report result.
    !
    IF ( .NOT.Fatal ) THEN
      IF ( Kprint>=3 ) THEN
        IF ( errmax<Thresh ) THEN
          WRITE (Nout,FMT=99001) Sname, nc
        ELSE
          WRITE (Nout,FMT=99003) Sname, nc, errmax
        END IF
      END IF
    END IF
    RETURN
    !
    99001 FORMAT (' ',A6,' PASSED THE COMPUTATIONAL TESTS (',I6,' CALL','S)')
    99002 FORMAT (' ******* FATAL ERROR - PARAMETER NUMBER ',I2,' WAS CH',&
      'ANGED INCORRECTLY *******')
    99003 FORMAT (' ',A6,' COMPLETED THE COMPUTATIONAL TESTS (',I6,' C','ALLS)',&
      /' ******* BUT WITH MAXIMUM TEST RATIO',F8.2,' - SUSPECT *******')
    99004 FORMAT (' ******* ',A6,' FAILED ON CALL NUMBER:')
    99005 FORMAT (1X,I6,': ',A6,'(''',A1,''',',4(I3,','),F4.1,', A,',I3,', X,',I2,&
      ',',F4.1,', Y,',I2,') .')
    99006 FORMAT (1X,I6,': ',A6,'(''',A1,''',',2(I3,','),F4.1,', A,',I3,', X,',I2,&
      ',',F4.1,', Y,',I2,')         .')
    99007 FORMAT (' ******* FATAL ERROR - ERROR-EXIT TAKEN ON VALID CALL *',&
      '******')
    !
    !     End of DCHK12.
    !
  END SUBROUTINE DCHK12
  !** DCHK13
  SUBROUTINE DCHK13(Sname,Eps,Thresh,Nout,Kprint,Fatal,Nidim,Idimm,Nalf,Alf,&
      Nbet,Bet,Nmax,A,Aa,As,B,Bb,Bs,C,Cc,Cs,Ct,G)
    !>
    !  Test DGEMM.
    !***
    ! **Library:**   SLATEC (BLAS)
    !***
    ! **Keywords:**  BLAS, QUICK CHECK SERVICE ROUTINE
    !***
    ! **Author:**  Dongarra, J. J., (ANL)
    !           Duff, I., (AERE)
    !           Du Croz, J., (NAG)
    !           Hammarling, S., (NAG)
    !***
    ! **Description:**
    !
    !  Quick check for DGEMM.
    !
    !  Auxiliary routine for test program for Level 3 Blas.
    !***
    ! **References:**  (NONE)
    !***
    ! **Routines called:**  DGEMM, DMAKE3, DMMCH, LDE, LDERES, NUMXER

    !* REVISION HISTORY  (YYMMDD)
    !   890208  DATE WRITTEN
    !   910619  Modified to meet SLATEC code and prologue standards. (BKS)
    USE slatec, ONLY : DGEMM, NUMXER
    !     .. Parameters ..
    REAL(8), PARAMETER :: ZERO = 0.0D0
    !     .. Scalar Arguments ..
    LOGICAL Fatal
    REAL(8) :: Eps, Thresh
    INTEGER Kprint, Nalf, Nbet, Nidim, Nmax, Nout
    CHARACTER(6) :: Sname
    !     .. Array Arguments ..
    REAL(8) :: A(Nmax,Nmax), Aa(Nmax*Nmax), Alf(Nalf), As(Nmax*Nmax), &
      G(Nmax), Bb(Nmax*Nmax), Bet(Nbet), Bs(Nmax*Nmax), &
      C(Nmax,Nmax), Cc(Nmax*Nmax), Cs(Nmax*Nmax), Ct(Nmax), B(Nmax,Nmax)
    INTEGER Idimm(Nidim)
    !     .. Local Scalars ..
    REAL(8) :: alpha, als, beta, bls, err, errmax
    INTEGER i, ia, ib, ica, icb, ik, im, in, k, ks, laa, lbb, &
      lcc, lda, ldas, ldb, ldbs, ldc, ldcs, m, ma, mb, ms, &
      n, na, nargs, nb, nc, nerr, ns
    LOGICAL ftl, nul, reset, trana, tranb
    CHARACTER :: tranas, tranbs, transa, transb
    !     .. Local Arrays ..
    LOGICAL isame(13)
    !     .. Intrinsic Functions ..
    INTRINSIC ABS, MAX, MIN
    !     .. Data statements ..
    CHARACTER(3), PARAMETER :: ich = 'NTC'
    !* FIRST EXECUTABLE STATEMENT  DCHK13
    nargs = 13
    nc = 0
    reset = .TRUE.
    errmax = ZERO
    !
    DO im = 1, Nidim
      m = Idimm(im)
      !
      DO in = 1, Nidim
        n = Idimm(in)
        !           Set LDC to 1 more than minimum value if room.
        ldc = m
        IF ( ldc<Nmax ) ldc = ldc + 1
        !           Skip tests if not enough room.
        IF ( ldc<=Nmax ) THEN
          lcc = ldc*n
          nul = n<=0 .OR. m<=0
          !
          DO ik = 1, Nidim
            k = Idimm(ik)
            !
            DO ica = 1, 3
              transa = ich(ica:ica)
              trana = transa=='T' .OR. transa=='C'
              !
              IF ( trana ) THEN
                ma = k
                na = m
              ELSE
                ma = m
                na = k
              END IF
              !                 Set LDA to 1 more than minimum value if room.
              lda = ma
              IF ( lda<Nmax ) lda = lda + 1
              !                 Skip tests if not enough room.
              IF ( lda<=Nmax ) THEN
                laa = lda*na
                !
                !                 Generate the matrix A.
                !
                CALL DMAKE3('GE',' ',' ',ma,na,A,Nmax,Aa,lda,reset,ZERO)
                !
                DO icb = 1, 3
                  transb = ich(icb:icb)
                  tranb = transb=='T' .OR. transb=='C'
                  !
                  IF ( tranb ) THEN
                    mb = n
                    nb = k
                  ELSE
                    mb = k
                    nb = n
                  END IF
                  !                    Set LDB to 1 more than minimum value if room.
                  ldb = mb
                  IF ( ldb<Nmax ) ldb = ldb + 1
                  !                    Skip tests if not enough room.
                  IF ( ldb<=Nmax ) THEN
                    lbb = ldb*nb
                    !
                    !                    Generate the matrix B.
                    !
                    CALL DMAKE3('GE',' ',' ',mb,nb,B,Nmax,Bb,ldb,reset,ZERO)
                    !
                    DO ia = 1, Nalf
                      alpha = Alf(ia)
                      !
                      DO ib = 1, Nbet
                        beta = Bet(ib)
                        !
                        !                          Generate the matrix C.
                        !
                        CALL DMAKE3('GE',' ',' ',m,n,C,Nmax,Cc,ldc,reset,ZERO)
                        !
                        nc = nc + 1
                        !
                        !                          Save every datum before calling the
                        !                          subroutine.
                        !
                        tranas = transa
                        tranbs = transb
                        ms = m
                        ns = n
                        ks = k
                        als = alpha
                        DO i = 1, laa
                          As(i) = Aa(i)
                        END DO
                        ldas = lda
                        DO i = 1, lbb
                          Bs(i) = Bb(i)
                        END DO
                        ldbs = ldb
                        bls = beta
                        DO i = 1, lcc
                          Cs(i) = Cc(i)
                        END DO
                        ldcs = ldc
                        !
                        !                          Call the subroutine.
                        !
                        CALL DGEMM(transa,transb,m,n,k,alpha,Aa,lda,Bb,ldb,&
                          beta,Cc,ldc)
                        !
                        !                          Check if error-exit was taken incorrectly.
                        !
                        IF ( NUMXER(nerr)/=0 ) THEN
                          IF ( Kprint>=2 ) WRITE (Nout,FMT=99006)
                          Fatal = .TRUE.
                        END IF
                        !
                        !                          See what data changed inside subroutines.
                        !
                        isame(1) = transa==tranas
                        isame(2) = transb==tranbs
                        isame(3) = ms==m
                        isame(4) = ns==n
                        isame(5) = ks==k
                        isame(6) = als==alpha
                        isame(7) = LDE(As,Aa,laa)
                        isame(8) = ldas==lda
                        isame(9) = LDE(Bs,Bb,lbb)
                        isame(10) = ldbs==ldb
                        isame(11) = bls==beta
                        IF ( nul ) THEN
                          isame(12) = LDE(Cs,Cc,lcc)
                        ELSE
                          isame(12) = LDERES('GE',' ',m,n,Cs,Cc,ldc)
                        END IF
                        isame(13) = ldcs==ldc
                        !
                        !                          If data was incorrectly changed, report
                        !                          and return.
                        !
                        DO i = 1, nargs
                          IF ( .NOT.isame(i) ) THEN
                            Fatal = .TRUE.
                            IF ( Kprint>=2 ) WRITE (Nout,FMT=99002) i
                          END IF
                        END DO
                        !
                        ftl = .FALSE.
                        IF ( .NOT.nul ) THEN
                          !
                          !                             Check the result.
                          !
                          CALL DMMCH(transa,transb,m,n,k,alpha,A,Nmax,B,Nmax,&
                            beta,C,Nmax,Ct,G,Cc,ldc,Eps,err,ftl,Nout,.TRUE.,Kprint)
                          errmax = MAX(errmax,err)
                        END IF
                        IF ( ftl ) THEN
                          Fatal = .TRUE.
                          IF ( Kprint>=3 ) THEN
                            WRITE (Nout,FMT=99004) Sname
                            WRITE (Nout,FMT=99005) nc, Sname, transa, &
                              transb, m, n, k, alpha, lda, ldb, beta, ldc
                          END IF
                        END IF
                        !
                      END DO
                      !
                    END DO
                  END IF
                  !
                END DO
              END IF
              !
            END DO
            !
          END DO
        END IF
        !
      END DO
      !
    END DO
    !
    !     Report result.
    !
    IF ( .NOT.Fatal ) THEN
      IF ( Kprint>=3 ) THEN
        IF ( errmax<Thresh ) THEN
          WRITE (Nout,FMT=99001) Sname, nc
        ELSE
          WRITE (Nout,FMT=99003) Sname, nc, errmax
        END IF
      END IF
    END IF
    RETURN
    !
    99001 FORMAT (' ',A6,' PASSED THE COMPUTATIONAL TESTS (',I6,' CALL','S)')
    99002 FORMAT (' ******* FATAL ERROR - PARAMETER NUMBER ',I2,' WAS CH',&
      'ANGED INCORRECTLY *******')
    99003 FORMAT (' ',A6,' COMPLETED THE COMPUTATIONAL TESTS (',I6,' C','ALLS)',&
      /' ******* BUT WITH MAXIMUM TEST RATIO',F8.2,' - SUSPECT *******')
    99004 FORMAT (' ******* ',A6,' FAILED ON CALL NUMBER:')
    99005 FORMAT (1X,I6,': ',A6,'(''',A1,''',''',A1,''',',3(I3,','),F4.1,', A,',I3,&
      ', B,',I3,',',F4.1,', ','C,',I3,').')
    99006 FORMAT (' ******* FATAL ERROR - ERROR-EXIT TAKEN ON VALID CALL *',&
      '******')
    !
    !     End of DCHK13.
    !
  END SUBROUTINE DCHK13
  !** DCHK22
  SUBROUTINE DCHK22(Sname,Eps,Thresh,Nout,Kprint,Fatal,Nidim,Idimm,Nkb,Kb,&
      Nalf,Alf,Nbet,Bet,Ninc,Inc,Nmax,Incmax,A,Aa,As,X,Xx,Xs,Y,Yy,Ys,Yt,G)
    !>
    !  Test DSYMV, DSBMV and DSPMV.
    !***
    ! **Library:**   SLATEC (BLAS)
    !***
    ! **Keywords:**  BLAS, QUICK CHECK SERVICE ROUTINE
    !***
    ! **Author:**  Du Croz, J. (NAG)
    !           Hanson, R. J. (SNLA)
    !***
    ! **Description:**
    !
    !  Quick check for DSYMV, DSBMV and DSPMV.
    !
    !  Auxiliary routine for test program for Level 2 Blas.
    !***
    ! **References:**  (NONE)
    !***
    ! **Routines called:**  DMAKE2, DMVCH, DSBMV, DSPMV, DSYMV, LDE, LDERES,
    !                    NUMXER

    !* REVISION HISTORY  (YYMMDD)
    !   870810  DATE WRITTEN
    !   910619  Modified to meet SLATEC code and prologue standards. (BKS)
    USE slatec, ONLY : DSBMV, DSPMV, DSYMV, NUMXER
    !     .. Parameters ..
    REAL(8), PARAMETER :: ZERO = 0.0D0, HALF = 0.5D0
    !     .. Scalar Arguments ..
    LOGICAL Fatal
    REAL(8) :: Eps, Thresh
    INTEGER Incmax, Kprint, Nalf, Nbet, Nidim, Ninc, Nkb, Nmax, Nout
    CHARACTER(6) :: Sname
    !     .. Array Arguments ..
    REAL(8) :: A(Nmax,Nmax), Aa(Nmax*Nmax), Alf(Nalf), As(Nmax*Nmax), &
      Bet(Nbet), G(Nmax), X(Nmax), Xs(Nmax*Incmax), &
      Xx(Nmax*Incmax), Y(Nmax), Ys(Nmax*Incmax), Yt(Nmax), Yy(Nmax*Incmax)
    INTEGER Idimm(Nidim), Inc(Ninc), Kb(Nkb)
    !     .. Local Scalars ..
    REAL(8) :: alpha, als, beta, bls, err, errmax, transl
    INTEGER i, ia, ib, ic, ik, in, incx, incxs, incy, incys, ix, &
      iy, k, ks, laa, lda, ldas, lx, ly, n, nargs, nc, nk, ns, nerr
    LOGICAL banded, ftl, full, nul, packed, reset
    CHARACTER :: uplo, uplos
    !     .. Local Arrays ..
    LOGICAL isame(13)
    !     .. Intrinsic Functions ..
    INTRINSIC ABS, MAX
    !     .. Data statements ..
    CHARACTER(2), PARAMETER :: ich = 'UL'
    !* FIRST EXECUTABLE STATEMENT  DCHK22
    full = Sname(3:3)=='Y'
    banded = Sname(3:3)=='B'
    packed = Sname(3:3)=='P'
    !     Define the number of arguments.
    IF ( full ) THEN
      nargs = 10
    ELSEIF ( banded ) THEN
      nargs = 11
    ELSEIF ( packed ) THEN
      nargs = 9
    END IF
    !
    nc = 0
    reset = .TRUE.
    errmax = ZERO
    !
    DO in = 1, Nidim
      n = Idimm(in)
      !
      IF ( banded ) THEN
        nk = Nkb
      ELSE
        nk = 1
      END IF
      DO ik = 1, nk
        IF ( banded ) THEN
          k = Kb(ik)
        ELSE
          k = n - 1
        END IF
        !           Set LDA to 1 more than minimum value if room.
        IF ( banded ) THEN
          lda = k + 1
        ELSE
          lda = n
        END IF
        IF ( lda<Nmax ) lda = lda + 1
        !           Skip tests if not enough room.
        IF ( lda<=Nmax ) THEN
          IF ( packed ) THEN
            laa = (n*(n+1))/2
          ELSE
            laa = lda*n
          END IF
          nul = n<=0
          !
          DO ic = 1, 2
            uplo = ich(ic:ic)
            !
            !              Generate the matrix A.
            !
            transl = ZERO
            CALL DMAKE2(Sname(2:3),uplo,' ',n,n,A,Nmax,Aa,lda,k,k,reset,transl)
            !
            DO ix = 1, Ninc
              incx = Inc(ix)
              lx = ABS(incx)*n
              !
              !                 Generate the vector X.
              !
              transl = HALF
              CALL DMAKE2('GE',' ',' ',1,n,X,1,Xx,ABS(incx),0,n-1,reset,transl)
              IF ( n>1 ) THEN
                X(n/2) = ZERO
                Xx(1+ABS(incx)*(n/2-1)) = ZERO
              END IF
              !
              DO iy = 1, Ninc
                incy = Inc(iy)
                ly = ABS(incy)*n
                !
                DO ia = 1, Nalf
                  alpha = Alf(ia)
                  !
                  DO ib = 1, Nbet
                    beta = Bet(ib)
                    !
                    !                          Generate the vector Y.
                    !
                    transl = ZERO
                    CALL DMAKE2('GE',' ',' ',1,n,Y,1,Yy,ABS(incy),0,n-1,reset,&
                      transl)
                    !
                    nc = nc + 1
                    !
                    !                          Save every datum before calling the
                    !                          subroutine.
                    !
                    uplos = uplo
                    ns = n
                    ks = k
                    als = alpha
                    DO i = 1, laa
                      As(i) = Aa(i)
                    END DO
                    ldas = lda
                    DO i = 1, lx
                      Xs(i) = Xx(i)
                    END DO
                    incxs = incx
                    bls = beta
                    DO i = 1, ly
                      Ys(i) = Yy(i)
                    END DO
                    incys = incy
                    !
                    !                          Call the subroutine.
                    !
                    IF ( full ) THEN
                      CALL DSYMV(uplo,n,alpha,Aa,lda,Xx,incx,beta,Yy,incy)
                    ELSEIF ( banded ) THEN
                      CALL DSBMV(uplo,n,k,alpha,Aa,lda,Xx,incx,beta,Yy,incy)
                    ELSEIF ( packed ) THEN
                      CALL DSPMV(uplo,n,alpha,Aa,Xx,incx,beta,Yy,incy)
                    END IF
                    !
                    !                          Check if error-exit was taken incorrectly.
                    !
                    IF ( NUMXER(nerr)/=0 ) THEN
                      IF ( Kprint>=2 ) WRITE (Nout,FMT=99008)
                      Fatal = .TRUE.
                    END IF
                    !
                    !                          See what data changed inside subroutines.
                    !
                    isame(1) = uplo==uplos
                    isame(2) = ns==n
                    IF ( full ) THEN
                      isame(3) = als==alpha
                      isame(4) = LDE(As,Aa,laa)
                      isame(5) = ldas==lda
                      isame(6) = LDE(Xs,Xx,lx)
                      isame(7) = incxs==incx
                      isame(8) = bls==beta
                      IF ( nul ) THEN
                        isame(9) = LDE(Ys,Yy,ly)
                      ELSE
                        isame(9) = LDERES('GE',' ',1,n,Ys,Yy,ABS(incy))
                      END IF
                      isame(10) = incys==incy
                    ELSEIF ( banded ) THEN
                      isame(3) = ks==k
                      isame(4) = als==alpha
                      isame(5) = LDE(As,Aa,laa)
                      isame(6) = ldas==lda
                      isame(7) = LDE(Xs,Xx,lx)
                      isame(8) = incxs==incx
                      isame(9) = bls==beta
                      IF ( nul ) THEN
                        isame(10) = LDE(Ys,Yy,ly)
                      ELSE
                        isame(10) = LDERES('GE',' ',1,n,Ys,Yy,ABS(incy))
                      END IF
                      isame(11) = incys==incy
                    ELSEIF ( packed ) THEN
                      isame(3) = als==alpha
                      isame(4) = LDE(As,Aa,laa)
                      isame(5) = LDE(Xs,Xx,lx)
                      isame(6) = incxs==incx
                      isame(7) = bls==beta
                      IF ( nul ) THEN
                        isame(8) = LDE(Ys,Yy,ly)
                      ELSE
                        isame(8) = LDERES('GE',' ',1,n,Ys,Yy,ABS(incy))
                      END IF
                      isame(9) = incys==incy
                    END IF
                    !
                    !                          If data was incorrectly changed, report and
                    !                          return.
                    !
                    DO i = 1, nargs
                      IF ( .NOT.isame(i) ) THEN
                        Fatal = .TRUE.
                        IF ( Kprint>=2 ) WRITE (Nout,FMT=99002) i
                      END IF
                    END DO
                    !
                    ftl = .FALSE.
                    IF ( .NOT.nul ) THEN
                      !
                      !                             Check the result.
                      !
                      CALL DMVCH('N',n,n,alpha,A,Nmax,X,incx,beta,Y,incy,Yt,G,&
                        Yy,Eps,err,ftl,Nout,.TRUE.,Kprint)
                      errmax = MAX(errmax,err)
                    END IF
                    IF ( ftl ) THEN
                      Fatal = .TRUE.
                      IF ( Kprint>=3 ) THEN
                        WRITE (Nout,FMT=99004) Sname
                        IF ( full ) THEN
                          WRITE (Nout,FMT=99007) nc, Sname, uplo, n, &
                            alpha, lda, incx, beta, incy
                        ELSEIF ( banded ) THEN
                          WRITE (Nout,FMT=99006) nc, Sname, uplo, n, &
                            alpha, incx, beta, incy
                        ELSEIF ( packed ) THEN
                          WRITE (Nout,FMT=99005) nc, Sname, uplo, n, &
                            alpha, incx, beta, incy
                        END IF
                      END IF
                    END IF
                    !
                  END DO
                  !
                END DO
                !
              END DO
              !
            END DO
            !
          END DO
        END IF
        !
      END DO
      !
    END DO
    !
    !     Report result.
    !
    IF ( .NOT.Fatal ) THEN
      IF ( Kprint>=3 ) THEN
        IF ( errmax<Thresh ) THEN
          WRITE (Nout,FMT=99001) Sname, nc
        ELSE
          WRITE (Nout,FMT=99003) Sname, nc, errmax
        END IF
      END IF
    END IF
    RETURN
    !
    99001 FORMAT (' ',A6,' PASSED THE COMPUTATIONAL TESTS (',I6,' CALL','S)')
    99002 FORMAT (' ******* FATAL ERROR - PARAMETER NUMBER ',I2,' WAS CH',&
      'ANGED INCORRECTLY *******')
    99003 FORMAT (' ',A6,' COMPLETED THE COMPUTATIONAL TESTS (',I6,' C','ALLS)',&
      /' ******* BUT WITH MAXIMUM TEST RATIO',F8.2,' - SUSPECT *******')
    99004 FORMAT (' ******* ',A6,' FAILED ON CALL NUMBER:')
    99005 FORMAT (1X,I6,': ',A6,'(''',A1,''',',I3,',',F4.1,', AP',', X,',I2,',',&
      F4.1,', Y,',I2,')                .')
    99006 FORMAT (1X,I6,': ',A6,'(''',A1,''',',2(I3,','),F4.1,', A,',I3,', X,',I2,&
      ',',F4.1,', Y,',I2,')         .')
    99007 FORMAT (1X,I6,': ',A6,'(''',A1,''',',I3,',',F4.1,', A,',I3,', X,',I2,',',&
      F4.1,', Y,',I2,')             .')
    99008 FORMAT (' ******* FATAL ERROR - ERROR-EXIT TAKEN ON VALID CALL *',&
      '******')
    !
    !     End of DCHK22.
    !
  END SUBROUTINE DCHK22
  !** DCHK23
  SUBROUTINE DCHK23(Sname,Eps,Thresh,Nout,Kprint,Fatal,Nidim,Idimm,Nalf,Alf,&
      Nbet,Bet,Nmax,A,Aa,As,B,Bb,Bs,C,Cc,Cs,Ct,G)
    !>
    !  Test DSYMM.
    !***
    ! **Library:**   SLATEC (BLAS)
    !***
    ! **Keywords:**  BLAS, QUICK CHECK SERVICE ROUTINE
    !***
    ! **Author:**  Dongarra, J. J., (ANL)
    !           Duff, I., (AERE)
    !           Du Croz, J., (NAG)
    !           Hammarling, S., (NAG)
    !***
    ! **Description:**
    !
    !  Quick check for DSYMM.
    !
    !  Auxiliary routine for test program for Level 3 Blas.
    !***
    ! **References:**  (NONE)
    !***
    ! **Routines called:**  DMAKE3, DMMCH, DSYMM, LDE, LDERES, NUMXER

    !* REVISION HISTORY  (YYMMDD)
    !   890208  DATE WRITTEN
    !   910619  Modified to meet SLATEC code and prologue standards. (BKS)
    USE slatec, ONLY : DSYMM, NUMXER
    !     .. Parameters ..
    REAL(8), PARAMETER :: ZERO = 0.0D0
    !     .. Scalar Arguments ..
    LOGICAL Fatal
    REAL(8) :: Eps, Thresh
    INTEGER Kprint, Nalf, Nbet, Nidim, Nmax, Nout
    CHARACTER(6) :: Sname
    !     .. Array Arguments ..
    REAL(8) :: A(Nmax,Nmax), Aa(Nmax*Nmax), Alf(Nalf), As(Nmax*Nmax), &
      G(Nmax), Bb(Nmax*Nmax), Bet(Nbet), Bs(Nmax*Nmax), &
      C(Nmax,Nmax), Cc(Nmax*Nmax), Cs(Nmax*Nmax), Ct(Nmax), B(Nmax,Nmax)
    INTEGER Idimm(Nidim)
    !     .. Local Scalars ..
    REAL(8) :: alpha, als, beta, bls, err, errmax
    INTEGER i, ia, ib, ics, icu, im, in, laa, lbb, lcc, lda, ldas, &
      ldb, ldbs, ldc, ldcs, m, ms, n, na, nargs, nc, nerr, ns
    LOGICAL ftl, left, nul, reset
    CHARACTER :: side, sides, uplo, uplos
    !     .. Local Arrays ..
    LOGICAL isame(13)
    !     .. Intrinsic Functions ..
    INTRINSIC ABS, MAX, MIN
    !     .. Data statements ..
    CHARACTER(2), PARAMETER :: ichs = 'LR', ichu = 'UL'
    !* FIRST EXECUTABLE STATEMENT  DCHK23
    nargs = 12
    nc = 0
    reset = .TRUE.
    errmax = ZERO
    !
    DO im = 1, Nidim
      m = Idimm(im)
      !
      DO in = 1, Nidim
        n = Idimm(in)
        !           Set LDC to 1 more than minimum value if room.
        ldc = m
        IF ( ldc<Nmax ) ldc = ldc + 1
        !           Skip tests if not enough room.
        IF ( ldc<=Nmax ) THEN
          lcc = ldc*n
          nul = n<=0 .OR. m<=0
          !
          !           Set LDB to 1 more than minimum value if room.
          ldb = m
          IF ( ldb<Nmax ) ldb = ldb + 1
          !           Skip tests if not enough room.
          IF ( ldb<=Nmax ) THEN
            lbb = ldb*n
            !
            !           Generate the matrix B.
            !
            CALL DMAKE3('GE',' ',' ',m,n,B,Nmax,Bb,ldb,reset,ZERO)
            !
            DO ics = 1, 2
              side = ichs(ics:ics)
              left = side=='L'
              !
              IF ( left ) THEN
                na = m
              ELSE
                na = n
              END IF
              !              Set LDA to 1 more than minimum value if room.
              lda = na
              IF ( lda<Nmax ) lda = lda + 1
              !              Skip tests if not enough room.
              IF ( lda<=Nmax ) THEN
                laa = lda*na
                !
                DO icu = 1, 2
                  uplo = ichu(icu:icu)
                  !
                  !                 Generate the symmetric matrix A.
                  !
                  CALL DMAKE3('SY',uplo,' ',na,na,A,Nmax,Aa,lda,reset,ZERO)
                  !
                  DO ia = 1, Nalf
                    alpha = Alf(ia)
                    !
                    DO ib = 1, Nbet
                      beta = Bet(ib)
                      !
                      !                       Generate the matrix C.
                      !
                      CALL DMAKE3('GE',' ',' ',m,n,C,Nmax,Cc,ldc,reset,ZERO)
                      !
                      nc = nc + 1
                      !
                      !                       Save every datum before calling the
                      !                       subroutine.
                      !
                      sides = side
                      uplos = uplo
                      ms = m
                      ns = n
                      als = alpha
                      DO i = 1, laa
                        As(i) = Aa(i)
                      END DO
                      ldas = lda
                      DO i = 1, lbb
                        Bs(i) = Bb(i)
                      END DO
                      ldbs = ldb
                      bls = beta
                      DO i = 1, lcc
                        Cs(i) = Cc(i)
                      END DO
                      ldcs = ldc
                      !
                      !                       Call the subroutine.
                      !
                      CALL DSYMM(side,uplo,m,n,alpha,Aa,lda,Bb,ldb,beta,Cc,ldc)
                      !
                      !                       Check if error-exit was taken incorrectly.
                      !
                      IF ( NUMXER(nerr)/=0 ) THEN
                        IF ( Kprint>=2 ) WRITE (Nout,FMT=99006)
                        Fatal = .TRUE.
                      END IF
                      !
                      !                       See what data changed inside subroutines.
                      !
                      isame(1) = sides==side
                      isame(2) = uplos==uplo
                      isame(3) = ms==m
                      isame(4) = ns==n
                      isame(5) = als==alpha
                      isame(6) = LDE(As,Aa,laa)
                      isame(7) = ldas==lda
                      isame(8) = LDE(Bs,Bb,lbb)
                      isame(9) = ldbs==ldb
                      isame(10) = bls==beta
                      IF ( nul ) THEN
                        isame(11) = LDE(Cs,Cc,lcc)
                      ELSE
                        isame(11) = LDERES('GE',' ',m,n,Cs,Cc,ldc)
                      END IF
                      isame(12) = ldcs==ldc
                      !
                      !                       If data was incorrectly changed, report and
                      !                       return.
                      !
                      DO i = 1, nargs
                        IF ( .NOT.isame(i) ) THEN
                          Fatal = .TRUE.
                          IF ( Kprint>=2 ) WRITE (Nout,FMT=99002) i
                        END IF
                      END DO
                      !
                      ftl = .FALSE.
                      IF ( .NOT.nul ) THEN
                        !
                        !                          Check the result.
                        !
                        IF ( left ) THEN
                          CALL DMMCH('N','N',m,n,m,alpha,A,Nmax,B,Nmax,beta,C,&
                            Nmax,Ct,G,Cc,ldc,Eps,err,ftl,Nout,.TRUE.,Kprint)
                        ELSE
                          CALL DMMCH('N','N',m,n,n,alpha,B,Nmax,A,Nmax,beta,C,&
                            Nmax,Ct,G,Cc,ldc,Eps,err,ftl,Nout,.TRUE.,Kprint)
                        END IF
                        errmax = MAX(errmax,err)
                      END IF
                      IF ( ftl ) THEN
                        Fatal = .TRUE.
                        IF ( Kprint>=3 ) THEN
                          WRITE (Nout,FMT=99004) Sname
                          WRITE (Nout,FMT=99005) nc, Sname, side, uplo, &
                            m, n, alpha, lda, ldb, beta, ldc
                        END IF
                      END IF
                      !
                    END DO
                    !
                  END DO
                  !
                END DO
              END IF
              !
            END DO
          END IF
        END IF
        !
      END DO
      !
    END DO
    !
    !     Report result.
    !
    IF ( .NOT.Fatal ) THEN
      IF ( Kprint>=3 ) THEN
        IF ( errmax<Thresh ) THEN
          WRITE (Nout,FMT=99001) Sname, nc
        ELSE
          WRITE (Nout,FMT=99003) Sname, nc, errmax
        END IF
      END IF
    END IF
    RETURN
    !
    99001 FORMAT (' ',A6,' PASSED THE COMPUTATIONAL TESTS (',I6,' CALL','S)')
    99002 FORMAT (' ******* FATAL ERROR - PARAMETER NUMBER ',I2,' WAS CH',&
      'ANGED INCORRECTLY *******')
    99003 FORMAT (' ',A6,' COMPLETED THE COMPUTATIONAL TESTS (',I6,' C','ALLS)',&
      /' ******* BUT WITH MAXIMUM TEST RATIO',F8.2,' - SUSPECT *******')
    99004 FORMAT (' ******* ',A6,' FAILED ON CALL NUMBER:')
    99005 FORMAT (1X,I6,': ',A6,'(',2('''',A1,''','),2(I3,','),F4.1,', A,',I3,&
      ', B,',I3,',',F4.1,', C,',I3,')   ',' .')
    99006 FORMAT (' ******* FATAL ERROR - ERROR-EXIT TAKEN ON VALID CALL *',&
      '******')
    !
    !     End of DCHK23.
    !
  END SUBROUTINE DCHK23
  !** DCHK32
  SUBROUTINE DCHK32(Sname,Eps,Thresh,Nout,Kprint,Fatal,Nidim,Idimm,Nkb,Kb,&
      Ninc,Inc,Nmax,Incmax,A,Aa,As,X,Xx,Xs,Xt,G,Z)
    !>
    !  Test DTRMV, DTBMV, DTPMV, DTRSV, DTBSV and DTPSV.
    !***
    ! **Library:**   SLATEC (BLAS)
    !***
    ! **Keywords:**  BLAS, QUICK CHECK SERVICE ROUTINE
    !***
    ! **Author:**  Du Croz, J. (NAG)
    !           Hanson, R. J. (SNLA)
    !***
    ! **Description:**
    !
    !  Quick check for DTRMV, DTBMV, DTPMV, DTRSV, DTBSV and DTPSV.
    !
    !  Auxiliary routine for test program for Level 2 Blas.
    !***
    ! **References:**  (NONE)
    !***
    ! **Routines called:**  DMAKE2, DMVCH, DTBMV, DTBSV, DTPMV, DTPSV, DTRMV,
    !                    DTRSV, LDE, LDERES, NUMXER

    !* REVISION HISTORY  (YYMMDD)
    !   870810  DATE WRITTEN
    !   910619  Modified to meet SLATEC code and prologue standards. (BKS)
    USE slatec, ONLY : DTBMV, DTBSV, DTPMV, DTPSV, DTRMV, DTRSV, NUMXER
    !     .. Parameters ..
    REAL(8), PARAMETER :: ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0
    !     .. Scalar Arguments ..
    LOGICAL Fatal
    REAL(8) :: Eps, Thresh
    INTEGER Incmax, Kprint, Nidim, Ninc, Nkb, Nmax, Nout
    CHARACTER(6) :: Sname
    !     .. Array Arguments ..
    REAL(8) :: A(Nmax,Nmax), Aa(Nmax*Nmax), As(Nmax*Nmax), G(Nmax), &
      X(Nmax), Xs(Nmax*Incmax), Xt(Nmax), Xx(Nmax*Incmax), Z(Nmax)
    INTEGER Idimm(Nidim), Inc(Ninc), Kb(Nkb)
    !     .. Local Scalars ..
    REAL(8) :: err, errmax, transl
    INTEGER i, icd, ict, icu, ik, in, incx, incxs, ix, k, ks, laa, &
      lda, ldas, lx, n, nargs, nc, nk, ns, nerr
    LOGICAL banded, ftl, full, nul, packed, reset
    CHARACTER :: diag, diags, trans, transs, uplo, uplos
    !     .. Local Arrays ..
    LOGICAL isame(13)
    !     .. Intrinsic Functions ..
    INTRINSIC ABS, MAX
    !     .. Data statements ..
    CHARACTER(2), PARAMETER :: ichu = 'UL', ichd = 'UN'
    CHARACTER(3), PARAMETER :: icht = 'NTC'
    !* FIRST EXECUTABLE STATEMENT  DCHK32
    full = Sname(3:3)=='R'
    banded = Sname(3:3)=='B'
    packed = Sname(3:3)=='P'
    !     Define the number of arguments.
    IF ( full ) THEN
      nargs = 8
    ELSEIF ( banded ) THEN
      nargs = 9
    ELSEIF ( packed ) THEN
      nargs = 7
    END IF
    !
    nc = 0
    reset = .TRUE.
    errmax = ZERO
    !     Set up zero vector for DMVCH.
    DO i = 1, Nmax
      Z(i) = ZERO
    END DO
    !
    DO in = 1, Nidim
      n = Idimm(in)
      !
      IF ( banded ) THEN
        nk = Nkb
      ELSE
        nk = 1
      END IF
      DO ik = 1, nk
        IF ( banded ) THEN
          k = Kb(ik)
        ELSE
          k = n - 1
        END IF
        !           Set LDA to 1 more than minimum value if room.
        IF ( banded ) THEN
          lda = k + 1
        ELSE
          lda = n
        END IF
        IF ( lda<Nmax ) lda = lda + 1
        !           Skip tests if not enough room.
        IF ( lda<=Nmax ) THEN
          IF ( packed ) THEN
            laa = (n*(n+1))/2
          ELSE
            laa = lda*n
          END IF
          nul = n<=0
          !
          DO icu = 1, 2
            uplo = ichu(icu:icu)
            !
            DO ict = 1, 3
              trans = icht(ict:ict)
              !
              DO icd = 1, 2
                diag = ichd(icd:icd)
                !
                !                    Generate the matrix A.
                !
                transl = ZERO
                CALL DMAKE2(Sname(2:3),uplo,diag,n,n,A,Nmax,Aa,lda,k,k,reset,&
                  transl)
                !
                DO ix = 1, Ninc
                  incx = Inc(ix)
                  lx = ABS(incx)*n
                  !
                  !                       Generate the vector X.
                  !
                  transl = HALF
                  CALL DMAKE2('GE',' ',' ',1,n,X,1,Xx,ABS(incx),0,n-1,reset,&
                    transl)
                  IF ( n>1 ) THEN
                    X(n/2) = ZERO
                    Xx(1+ABS(incx)*(n/2-1)) = ZERO
                  END IF
                  !
                  nc = nc + 1
                  !
                  !                       Save every datum before calling the subroutine.
                  !
                  uplos = uplo
                  transs = trans
                  diags = diag
                  ns = n
                  ks = k
                  DO i = 1, laa
                    As(i) = Aa(i)
                  END DO
                  ldas = lda
                  DO i = 1, lx
                    Xs(i) = Xx(i)
                  END DO
                  incxs = incx
                  !
                  !                       Call the subroutine.
                  !
                  IF ( Sname(4:5)=='MV' ) THEN
                    IF ( full ) THEN
                      CALL DTRMV(uplo,trans,diag,n,Aa,lda,Xx,incx)
                    ELSEIF ( banded ) THEN
                      CALL DTBMV(uplo,trans,diag,n,k,Aa,lda,Xx,incx)
                    ELSEIF ( packed ) THEN
                      CALL DTPMV(uplo,trans,diag,n,Aa,Xx,incx)
                    END IF
                  ELSEIF ( Sname(4:5)=='SV' ) THEN
                    IF ( full ) THEN
                      CALL DTRSV(uplo,trans,diag,n,Aa,lda,Xx,incx)
                    ELSEIF ( banded ) THEN
                      CALL DTBSV(uplo,trans,diag,n,k,Aa,lda,Xx,incx)
                    ELSEIF ( packed ) THEN
                      CALL DTPSV(uplo,trans,diag,n,Aa,Xx,incx)
                    END IF
                  END IF
                  !
                  !                       Check if error-exit was taken incorrectly.
                  !
                  IF ( NUMXER(nerr)/=0 ) THEN
                    IF ( Kprint>=2 ) WRITE (Nout,FMT=99008)
                    Fatal = .TRUE.
                  END IF
                  !
                  !                       See what data changed inside subroutines.
                  !
                  isame(1) = uplo==uplos
                  isame(2) = trans==transs
                  isame(3) = diag==diags
                  isame(4) = ns==n
                  IF ( full ) THEN
                    isame(5) = LDE(As,Aa,laa)
                    isame(6) = ldas==lda
                    IF ( nul ) THEN
                      isame(7) = LDE(Xs,Xx,lx)
                    ELSE
                      isame(7) = LDERES('GE',' ',1,n,Xs,Xx,ABS(incx))
                    END IF
                    isame(8) = incxs==incx
                  ELSEIF ( banded ) THEN
                    isame(5) = ks==k
                    isame(6) = LDE(As,Aa,laa)
                    isame(7) = ldas==lda
                    IF ( nul ) THEN
                      isame(8) = LDE(Xs,Xx,lx)
                    ELSE
                      isame(8) = LDERES('GE',' ',1,n,Xs,Xx,ABS(incx))
                    END IF
                    isame(9) = incxs==incx
                  ELSEIF ( packed ) THEN
                    isame(5) = LDE(As,Aa,laa)
                    IF ( nul ) THEN
                      isame(6) = LDE(Xs,Xx,lx)
                    ELSE
                      isame(6) = LDERES('GE',' ',1,n,Xs,Xx,ABS(incx))
                    END IF
                    isame(7) = incxs==incx
                  END IF
                  !
                  !                       If data was incorrectly changed, report and
                  !                       return.
                  !
                  DO i = 1, nargs
                    IF ( .NOT.isame(i) ) THEN
                      Fatal = .TRUE.
                      IF ( Kprint>=2 ) WRITE (Nout,FMT=99002) i
                    END IF
                  END DO
                  !
                  ftl = .FALSE.
                  IF ( .NOT.nul ) THEN
                    IF ( Sname(4:5)=='MV' ) THEN
                      !
                      !                             Check the result.
                      !
                      CALL DMVCH(trans,n,n,ONE,A,Nmax,X,incx,ZERO,Z,incx,Xt,G,&
                        Xx,Eps,err,ftl,Nout,.TRUE.,Kprint)
                    ELSEIF ( Sname(4:5)=='SV' ) THEN
                      !
                      !                             Compute approximation to original vector.
                      !
                      DO i = 1, n
                        Z(i) = Xx(1+(i-1)*ABS(incx))
                        Xx(1+(i-1)*ABS(incx)) = X(i)
                      END DO
                      CALL DMVCH(trans,n,n,ONE,A,Nmax,Z,incx,ZERO,X,incx,Xt,G,&
                        Xx,Eps,err,ftl,Nout,.FALSE.,Kprint)
                    END IF
                    errmax = MAX(errmax,err)
                  END IF
                  IF ( ftl ) THEN
                    Fatal = .TRUE.
                    IF ( Kprint>=3 ) THEN
                      WRITE (Nout,FMT=99004) Sname
                      IF ( full ) THEN
                        WRITE (Nout,FMT=99007) nc, Sname, uplo, trans, &
                          diag, n, lda, incx
                      ELSEIF ( banded ) THEN
                        WRITE (Nout,FMT=99006) nc, Sname, uplo, trans, &
                          diag, n, k, lda, incx
                      ELSEIF ( packed ) THEN
                        WRITE (Nout,FMT=99005) nc, Sname, uplo, trans, &
                          diag, n, incx
                      END IF
                    END IF
                  END IF
                  !
                END DO
                !
              END DO
              !
            END DO
            !
          END DO
        END IF
        !
      END DO
      !
    END DO
    !
    !     Report result.
    !
    IF ( .NOT.Fatal ) THEN
      IF ( Kprint>=3 ) THEN
        IF ( errmax<Thresh ) THEN
          WRITE (Nout,FMT=99001) Sname, nc
        ELSE
          WRITE (Nout,FMT=99003) Sname, nc, errmax
        END IF
      END IF
    END IF
    RETURN
    !
    99001 FORMAT (' ',A6,' PASSED THE COMPUTATIONAL TESTS (',I6,' CALL','S)')
    99002 FORMAT (' ******* FATAL ERROR - PARAMETER NUMBER ',I2,' WAS CH',&
      'ANGED INCORRECTLY *******')
    99003 FORMAT (' ',A6,' COMPLETED THE COMPUTATIONAL TESTS (',I6,' C','ALLS)',&
      /' ******* BUT WITH MAXIMUM TEST RATIO',F8.2,' - SUSPECT *******')
    99004 FORMAT (' ******* ',A6,' FAILED ON CALL NUMBER:')
    99005 FORMAT (1X,I6,': ',A6,'(',3('''',A1,''','),I3,', AP, ','X,',I2,&
      ')                        .')
    99006 FORMAT (1X,I6,': ',A6,'(',3('''',A1,''','),2(I3,','),' A,',I3,', X,',I2,&
      ')                 .')
    99007 FORMAT (1X,I6,': ',A6,'(',3('''',A1,''','),I3,', A,',I3,', X,',I2,&
      ')                     .')
    99008 FORMAT (' ******* FATAL ERROR - ERROR-EXIT TAKEN ON VALID CALL *',&
      '******')
    !
    !     End of DCHK32.
    !
  END SUBROUTINE DCHK32
  !** DCHK33
  SUBROUTINE DCHK33(Sname,Eps,Thresh,Nout,Kprint,Fatal,Nidim,Idimm,Nalf,Alf,&
      Nmax,A,Aa,As,B,Bb,Bs,Ct,G,C)
    !>
    !  Test DTRMM and DTRSM.
    !***
    ! **Library:**   SLATEC (BLAS)
    !***
    ! **Keywords:**  BLAS, QUICK CHECK SERVICE ROUTINE
    !***
    ! **Author:**  Dongarra, J. J., (ANL)
    !           Duff, I., (AERE)
    !           Du Croz, J., (NAG)
    !           Hammarling, S., (NAG)
    !***
    ! **Description:**
    !
    !  Quick check for DTRMM and DTRSM.
    !
    !  Auxiliary routine for test program for Level 3 Blas.
    !***
    ! **References:**  (NONE)
    !***
    ! **Routines called:**  DMAKE3, DMMCH, DTRMM, DTRSM, LDE, LDERES, NUMXER

    !* REVISION HISTORY  (YYMMDD)
    !   890208  DATE WRITTEN
    !   910619  Modified to meet SLATEC code and prologue standards. (BKS)
    USE slatec, ONLY : DTRMM, DTRSM, NUMXER
    !     .. Parameters ..
    REAL(8), PARAMETER :: ZERO = 0.0D0, ONE = 1.0D0
    !     .. Scalar Arguments ..
    LOGICAL Fatal
    REAL(8) :: Eps, Thresh
    INTEGER Kprint, Nalf, Nidim, Nmax, Nout
    CHARACTER(6) :: Sname
    !     .. Array Arguments ..
    REAL(8) :: A(Nmax,Nmax), Aa(Nmax*Nmax), Alf(Nalf), As(Nmax*Nmax), &
      G(Nmax), Bb(Nmax*Nmax), Bs(Nmax*Nmax), C(Nmax,Nmax), &
      Ct(Nmax), B(Nmax,Nmax)
    INTEGER Idimm(Nidim)
    !     .. Local Scalars ..
    REAL(8) :: alpha, als, err, errmax
    INTEGER i, ia, icd, ics, ict, icu, im, in, j, laa, lbb, lda, &
      ldas, ldb, ldbs, m, ms, n, na, nargs, nc, nerr, ns
    LOGICAL ftl, left, nul, reset
    CHARACTER :: diag, diags, side, sides, tranas, transa, uplo, uplos
    !     .. Local Arrays ..
    LOGICAL isame(13)
    !     .. Intrinsic Functions ..
    INTRINSIC ABS, MAX, MIN
    !     .. Data statements ..
    CHARACTER(2), PARAMETER :: ichu = 'UL', ichd = 'UN', ichs = 'LR'
    CHARACTER(3), PARAMETER :: icht = 'NTC'
    !* FIRST EXECUTABLE STATEMENT  DCHK33
    nargs = 11
    nc = 0
    reset = .TRUE.
    errmax = ZERO
    !     Set up zero matrix for DMMCH.
    DO j = 1, Nmax
      DO i = 1, Nmax
        C(i,j) = ZERO
      END DO
    END DO
    !
    DO im = 1, Nidim
      m = Idimm(im)
      !
      DO in = 1, Nidim
        n = Idimm(in)
        !           Set LDB to 1 more than minimum value if room.
        ldb = m
        IF ( ldb<Nmax ) ldb = ldb + 1
        !           Skip tests if not enough room.
        IF ( ldb<=Nmax ) THEN
          lbb = ldb*n
          nul = m<=0 .OR. n<=0
          !
          DO ics = 1, 2
            side = ichs(ics:ics)
            left = side=='L'
            IF ( left ) THEN
              na = m
            ELSE
              na = n
            END IF
            !              Set LDA to 1 more than minimum value if room.
            lda = na
            IF ( lda<Nmax ) lda = lda + 1
            !              Skip tests if not enough room.
            IF ( lda>Nmax ) EXIT
            laa = lda*na
            !
            DO icu = 1, 2
              uplo = ichu(icu:icu)
              !
              DO ict = 1, 3
                transa = icht(ict:ict)
                !
                DO icd = 1, 2
                  diag = ichd(icd:icd)
                  !
                  DO ia = 1, Nalf
                    alpha = Alf(ia)
                    !
                    !                          Generate the matrix A.
                    !
                    CALL DMAKE3('TR',uplo,diag,na,na,A,Nmax,Aa,lda,reset,ZERO)
                    !
                    !                          Generate the matrix B.
                    !
                    CALL DMAKE3('GE',' ',' ',m,n,B,Nmax,Bb,ldb,reset,ZERO)
                    !
                    nc = nc + 1
                    !
                    !                          Save every datum before calling the
                    !                          subroutine.
                    !
                    sides = side
                    uplos = uplo
                    tranas = transa
                    diags = diag
                    ms = m
                    ns = n
                    als = alpha
                    DO i = 1, laa
                      As(i) = Aa(i)
                    END DO
                    ldas = lda
                    DO i = 1, lbb
                      Bs(i) = Bb(i)
                    END DO
                    ldbs = ldb
                    !
                    !                          Call the subroutine.
                    !
                    IF ( Sname(4:5)=='MM' ) THEN
                      CALL DTRMM(side,uplo,transa,diag,m,n,alpha,Aa,lda,Bb,ldb)
                    ELSEIF ( Sname(4:5)=='SM' ) THEN
                      CALL DTRSM(side,uplo,transa,diag,m,n,alpha,Aa,lda,Bb,ldb)
                    END IF
                    !
                    !                          Check if error-exit was taken incorrectly.
                    !
                    IF ( NUMXER(nerr)/=0 ) THEN
                      IF ( Kprint>=2 ) WRITE (Nout,FMT=99006)
                      Fatal = .TRUE.
                    END IF
                    !
                    !                          See what data changed inside subroutines.
                    !
                    isame(1) = sides==side
                    isame(2) = uplos==uplo
                    isame(3) = tranas==transa
                    isame(4) = diags==diag
                    isame(5) = ms==m
                    isame(6) = ns==n
                    isame(7) = als==alpha
                    isame(8) = LDE(As,Aa,laa)
                    isame(9) = ldas==lda
                    IF ( nul ) THEN
                      isame(10) = LDE(Bs,Bb,lbb)
                    ELSE
                      isame(10) = LDERES('GE',' ',m,n,Bs,Bb,ldb)
                    END IF
                    isame(11) = ldbs==ldb
                    !
                    !                          If data was incorrectly changed, report and
                    !                          return.
                    !
                    DO i = 1, nargs
                      IF ( .NOT.isame(i) ) THEN
                        Fatal = .TRUE.
                        IF ( Kprint>=2 ) WRITE (Nout,FMT=99002) i
                      END IF
                    END DO
                    !
                    ftl = .FALSE.
                    IF ( .NOT.nul ) THEN
                      IF ( Sname(4:5)=='MM' ) THEN
                        !
                        !                                Check the result.
                        !
                        IF ( left ) THEN
                          CALL DMMCH(transa,'N',m,n,m,alpha,A,Nmax,B,Nmax,&
                            ZERO,C,Nmax,Ct,G,Bb,ldb,Eps,err,ftl,Nout,.TRUE.,Kprint)
                        ELSE
                          CALL DMMCH('N',transa,m,n,n,alpha,B,Nmax,A,Nmax,&
                            ZERO,C,Nmax,Ct,G,Bb,ldb,Eps,err,ftl,Nout,.TRUE.,Kprint)
                        END IF
                      ELSEIF ( Sname(4:5)=='SM' ) THEN
                        !
                        !                                Compute approximation to original
                        !                                matrix.
                        !
                        DO j = 1, n
                          DO i = 1, m
                            C(i,j) = Bb(i+(j-1)*ldb)
                            Bb(i+(j-1)*ldb) = alpha*B(i,j)
                          END DO
                        END DO
                        !
                        IF ( left ) THEN
                          CALL DMMCH(transa,'N',m,n,m,ONE,A,Nmax,C,Nmax,ZERO,&
                            B,Nmax,Ct,G,Bb,ldb,Eps,err,ftl,Nout,.FALSE.,Kprint)
                        ELSE
                          CALL DMMCH('N',transa,m,n,n,ONE,C,Nmax,A,Nmax,ZERO,&
                            B,Nmax,Ct,G,Bb,ldb,Eps,err,ftl,Nout,.FALSE.,Kprint)
                        END IF
                      END IF
                      errmax = MAX(errmax,err)
                    END IF
                    IF ( ftl ) THEN
                      Fatal = .TRUE.
                      IF ( Kprint>=3 ) THEN
                        WRITE (Nout,FMT=99004) Sname
                        WRITE (Nout,FMT=99005) nc, Sname, side, uplo, &
                          transa, diag, m, n, alpha, lda, ldb
                      END IF
                    END IF
                    !
                  END DO
                  !
                END DO
                !
              END DO
              !
            END DO
            !
          END DO
        END IF
        !
      END DO
      !
    END DO
    !
    !     Report result.
    !
    IF ( .NOT.Fatal ) THEN
      IF ( Kprint>=3 ) THEN
        IF ( errmax<Thresh ) THEN
          WRITE (Nout,FMT=99001) Sname, nc
        ELSE
          WRITE (Nout,FMT=99003) Sname, nc, errmax
        END IF
      END IF
    END IF
    RETURN
    !
    99001 FORMAT (' ',A6,' PASSED THE COMPUTATIONAL TESTS (',I6,' CALL','S)')
    99002 FORMAT (' ******* FATAL ERROR - PARAMETER NUMBER ',I2,' WAS CH',&
      'ANGED INCORRECTLY *******')
    99003 FORMAT (' ',A6,' COMPLETED THE COMPUTATIONAL TESTS (',I6,' C','ALLS)',&
      /' ******* BUT WITH MAXIMUM TEST RATIO',F8.2,' - SUSPECT *******')
    99004 FORMAT (' ******* ',A6,' FAILED ON CALL NUMBER:')
    99005 FORMAT (1X,I6,': ',A6,'(',4('''',A1,''','),2(I3,','),F4.1,', A,',I3,&
      ', B,',I3,')        .')
    99006 FORMAT (' ******* FATAL ERROR - ERROR-EXIT TAKEN ON VALID CALL *',&
      '******')
    !
    !     End of DCHK33.
    !
  END SUBROUTINE DCHK33
  !** DCHK42
  SUBROUTINE DCHK42(Sname,Eps,Thresh,Nout,Kprint,Fatal,Nidim,Idimm,Nalf,Alf,&
      Ninc,Inc,Nmax,Incmax,A,Aa,As,X,Xx,Xs,Y,Yy,Ys,Yt,G,Z)
    !>
    !  Test DGER.
    !***
    ! **Library:**   SLATEC (BLAS)
    !***
    ! **Keywords:**  BLAS, QUICK CHECK SERVICE ROUTINE
    !***
    ! **Author:**  Du Croz, J. (NAG)
    !           Hanson, R. J. (SNLA)
    !***
    ! **Description:**
    !
    !  Quick check for DGER.
    !
    !  Auxiliary routine for test program for Level 2 Blas.
    !***
    ! **References:**  (NONE)
    !***
    ! **Routines called:**  DGER, DMAKE2, DMVCH, LDE, LDERES, NUMXER

    !* REVISION HISTORY  (YYMMDD)
    !   870810  DATE WRITTEN
    !   910619  Modified to meet SLATEC code and prologue standards. (BKS)
    USE slatec, ONLY : DGER, NUMXER
    !     .. Parameters ..
    REAL(8), PARAMETER :: ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0
    !     .. Scalar Arguments ..
    LOGICAL Fatal
    REAL(8) :: Eps, Thresh
    INTEGER Incmax, Kprint, Nalf, Nidim, Ninc, Nmax, Nout
    CHARACTER(6) :: Sname
    !     .. Array Arguments ..
    REAL(8) :: A(Nmax,Nmax), Aa(Nmax*Nmax), Alf(Nalf), As(Nmax*Nmax), &
      G(Nmax), X(Nmax), Xs(Nmax*Incmax), Xx(Nmax*Incmax), &
      Y(Nmax), Ys(Nmax*Incmax), Yt(Nmax), Yy(Nmax*Incmax), Z(Nmax)
    INTEGER Idimm(Nidim), Inc(Ninc)
    !     .. Local Scalars ..
    REAL(8) :: alpha, als, err, errmax, transl
    INTEGER i, ia, im, in, incx, incxs, incy, incys, ix, iy, j, &
      laa, lda, ldas, lx, ly, m, ms, n, nargs, nc, nd, ns, nerr
    LOGICAL ftl, nul, reset
    !     .. Local Arrays ..
    REAL(8) :: w(1)
    LOGICAL isame(13)
    !     .. Intrinsic Functions ..
    INTRINSIC ABS, MAX, MIN
    !* FIRST EXECUTABLE STATEMENT  DCHK42
    !     Define the number of arguments.
    nargs = 9
    !
    nc = 0
    reset = .TRUE.
    errmax = ZERO
    !
    DO in = 1, Nidim
      n = Idimm(in)
      nd = n/2 + 1
      !
      DO im = 1, 2
        IF ( im==1 ) m = MAX(n-nd,0)
        IF ( im==2 ) m = MIN(n+nd,Nmax)
        !
        !           Set LDA to 1 more than minimum value if room.
        lda = m
        IF ( lda<Nmax ) lda = lda + 1
        !           Skip tests if not enough room.
        IF ( lda<=Nmax ) THEN
          laa = lda*n
          nul = n<=0 .OR. m<=0
          !
          DO ix = 1, Ninc
            incx = Inc(ix)
            lx = ABS(incx)*m
            !
            !              Generate the vector X.
            !
            transl = HALF
            CALL DMAKE2('GE',' ',' ',1,m,X,1,Xx,ABS(incx),0,m-1,reset,transl)
            IF ( m>1 ) THEN
              X(m/2) = ZERO
              Xx(1+ABS(incx)*(m/2-1)) = ZERO
            END IF
            !
            DO iy = 1, Ninc
              incy = Inc(iy)
              ly = ABS(incy)*n
              !
              !                 Generate the vector Y.
              !
              transl = ZERO
              CALL DMAKE2('GE',' ',' ',1,n,Y,1,Yy,ABS(incy),0,n-1,reset,transl)
              IF ( n>1 ) THEN
                Y(n/2) = ZERO
                Yy(1+ABS(incy)*(n/2-1)) = ZERO
              END IF
              !
              DO ia = 1, Nalf
                alpha = Alf(ia)
                !
                !                    Generate the matrix A.
                !
                transl = ZERO
                CALL DMAKE2(Sname(2:3),' ',' ',m,n,A,Nmax,Aa,lda,m-1,n-1,&
                  reset,transl)
                !
                nc = nc + 1
                !
                !                    Save every datum before calling the subroutine.
                !
                ms = m
                ns = n
                als = alpha
                DO i = 1, laa
                  As(i) = Aa(i)
                END DO
                ldas = lda
                DO i = 1, lx
                  Xs(i) = Xx(i)
                END DO
                incxs = incx
                DO i = 1, ly
                  Ys(i) = Yy(i)
                END DO
                incys = incy
                !
                !                    Call the subroutine.
                !
                CALL DGER(m,n,alpha,Xx,incx,Yy,incy,Aa,lda)
                !
                !                    Check if error-exit was taken incorrectly.
                !
                IF ( NUMXER(nerr)/=0 ) THEN
                  IF ( Kprint>=2 ) WRITE (Nout,FMT=99007)
                  Fatal = .TRUE.
                END IF
                !
                !                    See what data changed inside subroutine.
                !
                isame(1) = ms==m
                isame(2) = ns==n
                isame(3) = als==alpha
                isame(4) = LDE(Xs,Xx,lx)
                isame(5) = incxs==incx
                isame(6) = LDE(Ys,Yy,ly)
                isame(7) = incys==incy
                IF ( nul ) THEN
                  isame(8) = LDE(As,Aa,laa)
                ELSE
                  isame(8) = LDERES('GE',' ',m,n,As,Aa,lda)
                END IF
                isame(9) = ldas==lda
                !
                !                    If data was incorrectly changed, report and return.
                !
                DO i = 1, nargs
                  IF ( .NOT.isame(i) ) THEN
                    Fatal = .TRUE.
                    IF ( Kprint>=2 ) WRITE (Nout,FMT=99002) i
                  END IF
                END DO
                !
                ftl = .FALSE.
                IF ( .NOT.nul ) THEN
                  !
                  !                       Check the result column by column.
                  !
                  IF ( incx>0 ) THEN
                    DO i = 1, m
                      Z(i) = X(i)
                    END DO
                  ELSE
                    DO i = 1, m
                      Z(i) = X(m-i+1)
                    END DO
                  END IF
                  DO j = 1, n
                    IF ( incy>0 ) THEN
                      w(1) = Y(j)
                    ELSE
                      w(1) = Y(n-j+1)
                    END IF
                    CALL DMVCH('N',m,1,alpha,Z,Nmax,w,1,ONE,A(1,j),1,Yt,G,&
                      Aa(1+(j-1)*lda),Eps,err,ftl,Nout,.TRUE.,Kprint)
                    errmax = MAX(errmax,err)
                  END DO
                END IF
                IF ( ftl ) THEN
                  Fatal = .TRUE.
                  IF ( Kprint>=3 ) THEN
                    WRITE (Nout,FMT=99005) j
                    WRITE (Nout,FMT=99004) Sname
                    WRITE (Nout,FMT=99006) nc, Sname, m, n, alpha, incx, incy, lda
                  END IF
                END IF
                !
              END DO
              !
            END DO
            !
          END DO
        END IF
        !
      END DO
      !
    END DO
    !
    !     Report result.
    !
    IF ( .NOT.Fatal ) THEN
      IF ( Kprint>=3 ) THEN
        IF ( errmax<Thresh ) THEN
          WRITE (Nout,FMT=99001) Sname, nc
        ELSE
          WRITE (Nout,FMT=99003) Sname, nc, errmax
        END IF
      END IF
    END IF
    RETURN
    !
    99001 FORMAT (' ',A6,' PASSED THE COMPUTATIONAL TESTS (',I6,' CALL','S)')
    99002 FORMAT (' ******* FATAL ERROR - PARAMETER NUMBER ',I2,' WAS CH',&
      'ANGED INCORRECTLY *******')
    99003 FORMAT (' ',A6,' COMPLETED THE COMPUTATIONAL TESTS (',I6,' C','ALLS)',&
      /' ******* BUT WITH MAXIMUM TEST RATIO',F8.2,' - SUSPECT *******')
    99004 FORMAT (' ******* ',A6,' FAILED ON CALL NUMBER:')
    99005 FORMAT ('      THESE ARE THE RESULTS FOR COLUMN ',I3)
    99006 FORMAT (1X,I6,': ',A6,'(',2(I3,','),F4.1,', X,',I2,', Y,',I2,', A,',I3,&
      ')                  .')
    99007 FORMAT (' ******* FATAL ERROR - ERROR-EXIT TAKEN ON VALID CALL *',&
      '******')
    !
    !     End of DCHK42.
    !
  END SUBROUTINE DCHK42
  !** DCHK43
  SUBROUTINE DCHK43(Sname,Eps,Thresh,Nout,Kprint,Fatal,Nidim,Idimm,Nalf,Alf,&
      Nbet,Bet,Nmax,A,Aa,As,C,Cc,Cs,Ct,G)
    !>
    !  Test DSYRK.
    !***
    ! **Library:**   SLATEC (BLAS)
    !***
    ! **Keywords:**  BLAS, QUICK CHECK SERVICE ROUTINE
    !***
    ! **Author:**  Dongarra, J. J., (ANL)
    !           Duff, I., (AERE)
    !           Du Croz, J., (NAG)
    !           Hammarling, S., (NAG)
    !***
    ! **Description:**
    !
    !  Quick check for DSYRK.
    !
    !  Auxiliary routine for test program for Level 3 Blas.
    !***
    ! **References:**  (NONE)
    !***
    ! **Routines called:**  DMAKE3, DMMCH, DSYRK, LDE, LDERES, NUMXER

    !* REVISION HISTORY  (YYMMDD)
    !   890208  DATE WRITTEN
    !   910619  Modified to meet SLATEC code and prologue standards. (BKS)
    USE slatec, ONLY : DSYRK, NUMXER
    !     .. Parameters ..
    REAL(8), PARAMETER :: ZERO = 0.0D0
    !     .. Scalar Arguments ..
    LOGICAL Fatal
    REAL(8) :: Eps, Thresh
    INTEGER Kprint, Nalf, Nbet, Nidim, Nmax, Nout
    CHARACTER(6) :: Sname
    !     .. Array Arguments ..
    REAL(8) :: A(Nmax,Nmax), Aa(Nmax*Nmax), Alf(Nalf), As(Nmax*Nmax), &
      G(Nmax), Bet(Nbet), C(Nmax,Nmax), Cc(Nmax*Nmax), Ct(Nmax), Cs(Nmax*Nmax)
    INTEGER Idimm(Nidim)
    !     .. Local Scalars ..
    REAL(8) :: alpha, als, beta, bets, err, errmax
    INTEGER i, ia, ib, ict, icu, ik, in, j, jc, jj, k, ks, laa, &
      lcc, lda, ldas, ldc, ldcs, lj, ma, n, na, nerr, ns, nargs, nc
    LOGICAL ftl, nul, reset, tran, upper
    CHARACTER :: trans, transs, uplo, uplos
    !     .. Local Arrays ..
    LOGICAL isame(13)
    !     .. Intrinsic Functions ..
    INTRINSIC ABS, MAX, MIN
    !     .. Data statements ..
    CHARACTER(2), PARAMETER :: ichu = 'UL'
    CHARACTER(3), PARAMETER :: icht = 'NTC'
    !* FIRST EXECUTABLE STATEMENT  DCHK43
    nargs = 10
    nc = 0
    reset = .TRUE.
    errmax = ZERO
    !
    DO in = 1, Nidim
      n = Idimm(in)
      !        Set LDC to 1 more than minimum value if room.
      ldc = n
      IF ( ldc<Nmax ) ldc = ldc + 1
      !        Skip tests if not enough room.
      IF ( ldc<=Nmax ) THEN
        lcc = ldc*n
        nul = n<=0
        !
        DO ik = 1, Nidim
          k = Idimm(ik)
          !
          DO ict = 1, 3
            trans = icht(ict:ict)
            tran = trans=='T' .OR. trans=='C'
            IF ( tran ) THEN
              ma = k
              na = n
            ELSE
              ma = n
              na = k
            END IF
            !              Set LDA to 1 more than minimum value if room.
            lda = ma
            IF ( lda<Nmax ) lda = lda + 1
            !              Skip tests if not enough room.
            IF ( lda<=Nmax ) THEN
              laa = lda*na
              !
              !              Generate the matrix A.
              !
              CALL DMAKE3('GE',' ',' ',ma,na,A,Nmax,Aa,lda,reset,ZERO)
              !
              DO icu = 1, 2
                uplo = ichu(icu:icu)
                upper = uplo=='U'
                !
                DO ia = 1, Nalf
                  alpha = Alf(ia)
                  !
                  DO ib = 1, Nbet
                    beta = Bet(ib)
                    !
                    !                       Generate the matrix C.
                    !
                    CALL DMAKE3('SY',uplo,' ',n,n,C,Nmax,Cc,ldc,reset,ZERO)
                    !
                    nc = nc + 1
                    !
                    !                       Save every datum before calling the subroutine.
                    !
                    uplos = uplo
                    transs = trans
                    ns = n
                    ks = k
                    als = alpha
                    DO i = 1, laa
                      As(i) = Aa(i)
                    END DO
                    ldas = lda
                    bets = beta
                    DO i = 1, lcc
                      Cs(i) = Cc(i)
                    END DO
                    ldcs = ldc
                    !
                    !                       Call the subroutine.
                    !
                    CALL DSYRK(uplo,trans,n,k,alpha,Aa,lda,beta,Cc,ldc)
                    !
                    !                       Check if error-exit was taken incorrectly.
                    !
                    IF ( NUMXER(nerr)/=0 ) THEN
                      IF ( Kprint>=2 ) WRITE (Nout,FMT=99006)
                      Fatal = .TRUE.
                    END IF
                    !
                    !                       See what data changed inside subroutines.
                    !
                    isame(1) = uplos==uplo
                    isame(2) = transs==trans
                    isame(3) = ns==n
                    isame(4) = ks==k
                    isame(5) = als==alpha
                    isame(6) = LDE(As,Aa,laa)
                    isame(7) = ldas==lda
                    isame(8) = bets==beta
                    IF ( nul ) THEN
                      isame(9) = LDE(Cs,Cc,lcc)
                    ELSE
                      isame(9) = LDERES('SY',uplo,n,n,Cs,Cc,ldc)
                    END IF
                    isame(10) = ldcs==ldc
                    !
                    !                       If data was incorrectly changed, report and
                    !                       return.
                    !
                    DO i = 1, nargs
                      IF ( .NOT.isame(i) ) THEN
                        Fatal = .TRUE.
                        IF ( Kprint>=2 ) WRITE (Nout,FMT=99002) i
                      END IF
                    END DO
                    !
                    ftl = .FALSE.
                    IF ( .NOT.nul ) THEN
                      !
                      !                          Check the result column by column.
                      !
                      jc = 1
                      DO j = 1, n
                        IF ( upper ) THEN
                          jj = 1
                          lj = j
                        ELSE
                          jj = j
                          lj = n - j + 1
                        END IF
                        IF ( tran ) THEN
                          CALL DMMCH('T','N',lj,1,k,alpha,A(1,jj),Nmax,A(1,j),&
                            Nmax,beta,C(jj,j),Nmax,Ct,G,Cc(jc),ldc,&
                            Eps,err,ftl,Nout,.TRUE.,Kprint)
                        ELSE
                          CALL DMMCH('N','T',lj,1,k,alpha,A(jj,1),Nmax,A(j,1),&
                            Nmax,beta,C(jj,j),Nmax,Ct,G,Cc(jc),ldc,&
                            Eps,err,ftl,Nout,.TRUE.,Kprint)
                        END IF
                        IF ( upper ) THEN
                          jc = jc + ldc
                        ELSE
                          jc = jc + ldc + 1
                        END IF
                        errmax = MAX(errmax,err)
                      END DO
                    END IF
                    IF ( ftl ) THEN
                      Fatal = .TRUE.
                      IF ( Kprint>=3 ) THEN
                        WRITE (Nout,FMT=99004) Sname
                        WRITE (Nout,FMT=99005) nc, Sname, uplo, trans, n, &
                          k, alpha, lda, beta, ldc
                      END IF
                    END IF
                    !
                  END DO
                  !
                END DO
                !
              END DO
            END IF
            !
          END DO
          !
        END DO
      END IF
      !
    END DO
    !
    !     Report result.
    !
    IF ( .NOT.Fatal ) THEN
      IF ( Kprint>=3 ) THEN
        IF ( errmax<Thresh ) THEN
          WRITE (Nout,FMT=99001) Sname, nc
        ELSE
          WRITE (Nout,FMT=99003) Sname, nc, errmax
        END IF
      END IF
    END IF
    RETURN
    !
    99001 FORMAT (' ',A6,' PASSED THE COMPUTATIONAL TESTS (',I6,' CALL','S)')
    99002 FORMAT (' ******* FATAL ERROR - PARAMETER NUMBER ',I2,' WAS CH',&
      'ANGED INCORRECTLY *******')
    99003 FORMAT (' ',A6,' COMPLETED THE COMPUTATIONAL TESTS (',I6,' C','ALLS)',&
      /' ******* BUT WITH MAXIMUM TEST RATIO',F8.2,' - SUSPECT *******')
    99004 FORMAT (' ******* ',A6,' FAILED ON CALL NUMBER:')
    99005 FORMAT (1X,I6,': ',A6,'(',2('''',A1,''','),2(I3,','),F4.1,', A,',I3,',',&
      F4.1,', C,',I3,')           .')
    99006 FORMAT (' ******* FATAL ERROR - ERROR-EXIT TAKEN ON VALID CALL *',&
      '******')
    !
    !     End of DCHK43.
    !
  END SUBROUTINE DCHK43
  !** DCHK52
  SUBROUTINE DCHK52(Sname,Eps,Thresh,Nout,Kprint,Fatal,Nidim,Idimm,Nalf,Alf,&
      Ninc,Inc,Nmax,Incmax,A,Aa,As,X,Xx,Xs,Yt,G,Z)
    !>
    !  Quick check for DSYR and DSPR.
    !***
    ! **Library:**   SLATEC (BLAS)
    !***
    ! **Keywords:**  BLAS, QUICK CHECK SERVICE ROUTINE
    !***
    ! **Author:**  Du Croz, J. (NAG)
    !           Hanson, R. J. (SNLA)
    !***
    ! **Description:**
    !
    !  Quick check for DSYR and DSPR.
    !
    !  Auxiliary routine for test program for Level 2 Blas.
    !***
    ! **References:**  (NONE)
    !***
    ! **Routines called:**  DMAKE2, DMVCH, DSPR, DSYR, LDE, LDERES, NUMXER

    !* REVISION HISTORY  (YYMMDD)
    !   870810  DATE WRITTEN
    !   910619  Modified to meet SLATEC code and prologue standards. (BKS)
    USE slatec, ONLY : DSPR, DSYR, NUMXER
    !     .. Parameters ..
    REAL(8), PARAMETER :: ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0
    !     .. Scalar Arguments ..
    LOGICAL Fatal
    REAL(8) :: Eps, Thresh
    INTEGER Incmax, Kprint, Nalf, Nidim, Ninc, Nmax, Nout
    CHARACTER(6) :: Sname
    !     .. Array Arguments ..
    REAL(8) :: A(Nmax,Nmax), Aa(Nmax*Nmax), Alf(Nalf), As(Nmax*Nmax), &
      G(Nmax), X(Nmax), Xs(Nmax*Incmax), Xx(Nmax*Incmax), Yt(Nmax), Z(Nmax)
    INTEGER Idimm(Nidim), Inc(Ninc)
    !     .. Local Scalars ..
    REAL(8) :: alpha, als, err, errmax, transl
    INTEGER i, ia, ic, in, incx, incxs, ix, j, ja, jj, laa, lda, &
      ldas, lj, lx, n, nargs, nc, ns, nerr
    LOGICAL ftl, full, nul, packed, reset, upper
    CHARACTER :: uplo, uplos
    !     .. Local Arrays ..
    REAL(8) :: w(1)
    LOGICAL isame(13)
    !     .. Intrinsic Functions ..
    INTRINSIC ABS, MAX
    !     .. Data statements ..
    CHARACTER(2), PARAMETER :: ich = 'UL'
    !* FIRST EXECUTABLE STATEMENT  DCHK52
    full = Sname(3:3)=='Y'
    packed = Sname(3:3)=='P'
    !     Define the number of arguments.
    IF ( full ) THEN
      nargs = 7
    ELSEIF ( packed ) THEN
      nargs = 6
    END IF
    !
    nc = 0
    reset = .TRUE.
    errmax = ZERO
    !
    DO in = 1, Nidim
      n = Idimm(in)
      !        Set LDA to 1 more than minimum value if room.
      lda = n
      IF ( lda<Nmax ) lda = lda + 1
      !        Skip tests if not enough room.
      IF ( lda<=Nmax ) THEN
        IF ( packed ) THEN
          laa = (n*(n+1))/2
        ELSE
          laa = lda*n
        END IF
        !
        DO ic = 1, 2
          uplo = ich(ic:ic)
          upper = uplo=='U'
          !
          DO ix = 1, Ninc
            incx = Inc(ix)
            lx = ABS(incx)*n
            !
            !              Generate the vector X.
            !
            transl = HALF
            CALL DMAKE2('GE',' ',' ',1,n,X,1,Xx,ABS(incx),0,n-1,reset,transl)
            IF ( n>1 ) THEN
              X(n/2) = ZERO
              Xx(1+ABS(incx)*(n/2-1)) = ZERO
            END IF
            !
            DO ia = 1, Nalf
              alpha = Alf(ia)
              nul = n<=0 .OR. alpha==ZERO
              !
              !                 Generate the matrix A.
              !
              transl = ZERO
              CALL DMAKE2(Sname(2:3),uplo,' ',n,n,A,Nmax,Aa,lda,n-1,n-1,reset,&
                transl)
              !
              nc = nc + 1
              !
              !                 Save every datum before calling the subroutine.
              !
              uplos = uplo
              ns = n
              als = alpha
              DO i = 1, laa
                As(i) = Aa(i)
              END DO
              ldas = lda
              DO i = 1, lx
                Xs(i) = Xx(i)
              END DO
              incxs = incx
              !
              !                 Call the subroutine.
              !
              IF ( full ) THEN
                CALL DSYR(uplo,n,alpha,Xx,incx,Aa,lda)
              ELSEIF ( packed ) THEN
                CALL DSPR(uplo,n,alpha,Xx,incx,Aa)
              END IF
              !
              !                 Check if error-exit was taken incorrectly.
              !
              IF ( NUMXER(nerr)/=0 ) THEN
                IF ( Kprint>=2 ) WRITE (Nout,FMT=99008)
                Fatal = .TRUE.
              END IF
              !
              !                 See what data changed inside subroutines.
              !
              isame(1) = uplo==uplos
              isame(2) = ns==n
              isame(3) = als==alpha
              isame(4) = LDE(Xs,Xx,lx)
              isame(5) = incxs==incx
              IF ( nul ) THEN
                isame(6) = LDE(As,Aa,laa)
              ELSE
                isame(6) = LDERES(Sname(2:3),uplo,n,n,As,Aa,lda)
              END IF
              IF ( .NOT.packed ) isame(7) = ldas==lda
              !
              !                 If data was incorrectly changed, report and return.
              !
              DO i = 1, nargs
                IF ( .NOT.isame(i) ) THEN
                  Fatal = .TRUE.
                  IF ( Kprint>=2 ) WRITE (Nout,FMT=99002) i
                END IF
              END DO
              !
              ftl = .FALSE.
              IF ( .NOT.nul ) THEN
                !
                !                    Check the result column by column.
                !
                IF ( incx>0 ) THEN
                  DO i = 1, n
                    Z(i) = X(i)
                  END DO
                ELSE
                  DO i = 1, n
                    Z(i) = X(n-i+1)
                  END DO
                END IF
                ja = 1
                DO j = 1, n
                  w(1) = Z(j)
                  IF ( upper ) THEN
                    jj = 1
                    lj = j
                  ELSE
                    jj = j
                    lj = n - j + 1
                  END IF
                  ftl = .FALSE.
                  CALL DMVCH('N',lj,1,alpha,Z(jj),lj,w,1,ONE,A(jj,j),1,Yt,G,&
                    Aa(ja),Eps,err,ftl,Nout,.TRUE.,Kprint)
                  IF ( .NOT.(full) ) THEN
                    ja = ja + lj
                  ELSEIF ( upper ) THEN
                    ja = ja + lda
                  ELSE
                    ja = ja + lda + 1
                  END IF
                  errmax = MAX(errmax,err)
                END DO
              END IF
              IF ( ftl ) THEN
                Fatal = .TRUE.
                IF ( Kprint>=3 ) THEN
                  WRITE (Nout,FMT=99005) j
                  WRITE (Nout,FMT=99004) Sname
                  IF ( full ) THEN
                    WRITE (Nout,FMT=99007) nc, Sname, uplo, n, alpha, incx, lda
                  ELSEIF ( packed ) THEN
                    WRITE (Nout,FMT=99006) nc, Sname, uplo, n, alpha, incx
                  END IF
                END IF
              END IF
              !
            END DO
            !
          END DO
          !
        END DO
      END IF
      !
    END DO
    !
    !     Report result.
    !
    IF ( .NOT.Fatal ) THEN
      IF ( Kprint>=3 ) THEN
        IF ( errmax<Thresh ) THEN
          WRITE (Nout,FMT=99001) Sname, nc
        ELSE
          WRITE (Nout,FMT=99003) Sname, nc, errmax
        END IF
      END IF
    END IF
    RETURN
    !
    99001 FORMAT (' ',A6,' PASSED THE COMPUTATIONAL TESTS (',I6,' CALL','S)')
    99002 FORMAT (' ******* FATAL ERROR - PARAMETER NUMBER ',I2,' WAS CH',&
      'ANGED INCORRECTLY *******')
    99003 FORMAT (' ',A6,' COMPLETED THE COMPUTATIONAL TESTS (',I6,' C','ALLS)',&
      /' ******* BUT WITH MAXIMUM TEST RATIO',F8.2,' - SUSPECT *******')
    99004 FORMAT (' ******* ',A6,' FAILED ON CALL NUMBER:')
    99005 FORMAT ('      THESE ARE THE RESULTS FOR COLUMN ',I3)
    99006 FORMAT (1X,I6,': ',A6,'(''',A1,''',',I3,',',F4.1,', X,',I2,&
      ', AP)                           .')
    99007 FORMAT (1X,I6,': ',A6,'(''',A1,''',',I3,',',F4.1,', X,',I2,', A,',I3,&
      ')                        .')
    99008 FORMAT (' ******* FATAL ERROR - ERROR-EXIT TAKEN ON VALID CALL *',&
      '******')
    !
    !     End of DCHK52.
    !
  END SUBROUTINE DCHK52
  !** DCHK53
  SUBROUTINE DCHK53(Sname,Eps,Thresh,Nout,Kprint,Fatal,Nidim,Idimm,Nalf,Alf,&
      Nbet,Bet,Nmax,Ab,Aa,As,Bb,Bs,C,Cc,Cs,Ct,G,W)
    !>
    !  Test DSYR2K.
    !***
    ! **Library:**   SLATEC (BLAS)
    !***
    ! **Keywords:**  BLAS, QUICK CHECK SERVICE ROUTINE
    !***
    ! **Author:**  Dongarra, J. J., (ANL)
    !           Duff, I., (AERE)
    !           Du Croz, J., (NAG)
    !           Hammarling, S., (NAG)
    !***
    ! **Description:**
    !
    !  Quick check for DSYR2K.
    !
    !  Auxiliary routine for test program for Level 3 Blas.
    !***
    ! **References:**  (NONE)
    !***
    ! **Routines called:**  DMAKE3, DMMCH, DSYR2K, LDE, LDERES, NUMXER

    !* REVISION HISTORY  (YYMMDD)
    !   890208  DATE WRITTEN
    !   910619  Modified to meet SLATEC code and prologue standards. (BKS)
    USE slatec, ONLY : DSYR2K, NUMXER
    !     .. Parameters ..
    REAL(8), PARAMETER :: ZERO = 0.0D0
    !     .. Scalar Arguments ..
    LOGICAL Fatal
    REAL(8) :: Eps, Thresh
    INTEGER Kprint, Nalf, Nbet, Nidim, Nmax, Nout
    CHARACTER(6) :: Sname
    !     .. Array Arguments ..
    REAL(8) :: Aa(Nmax*Nmax), Ab(2*Nmax*Nmax), Alf(Nalf), &
      As(Nmax*Nmax), Bb(Nmax*Nmax), Bet(Nbet), Bs(Nmax*Nmax), &
      C(Nmax,Nmax), Cc(Nmax*Nmax), Cs(Nmax*Nmax), Ct(Nmax), &
      G(Nmax), W(2*Nmax)
    INTEGER Idimm(Nidim)
    !     .. Local Scalars ..
    REAL(8) :: alpha, als, beta, bets, err, errmax
    INTEGER i, ia, ib, ict, icu, ik, in, j, jc, jj, jjab, k, ks, laa, lbb, lcc, &
      lda, ldas, ldb, ldbs, ldc, ldcs, lj, ma, n, na, nargs, nc, nerr, ns
    LOGICAL ftl, nul, reset, tran, upper
    CHARACTER :: trans, transs, uplo, uplos
    !     .. Local Arrays ..
    LOGICAL isame(13)
    !     .. Intrinsic Functions ..
    INTRINSIC ABS, MAX, MIN
    !     .. Data statements ..
    CHARACTER(2), PARAMETER :: ichu = 'UL'
    CHARACTER(3), PARAMETER :: icht = 'NTC'
    !* FIRST EXECUTABLE STATEMENT  DCHK53
    nargs = 12
    nc = 0
    reset = .TRUE.
    errmax = ZERO
    !
    DO in = 1, Nidim
      n = Idimm(in)
      !        Set LDC to 1 more than minimum value if room.
      ldc = n
      IF ( ldc<Nmax ) ldc = ldc + 1
      !        Skip tests if not enough room.
      IF ( ldc<=Nmax ) THEN
        lcc = ldc*n
        nul = n<=0
        !
        DO ik = 1, Nidim
          k = Idimm(ik)
          !
          DO ict = 1, 3
            trans = icht(ict:ict)
            tran = trans=='T' .OR. trans=='C'
            IF ( tran ) THEN
              ma = k
              na = n
            ELSE
              ma = n
              na = k
            END IF
            !              Set LDA to 1 more than minimum value if room.
            lda = ma
            IF ( lda<Nmax ) lda = lda + 1
            !              Skip tests if not enough room.
            IF ( lda<=Nmax ) THEN
              laa = lda*na
              !
              !              Generate the matrix A.
              !
              IF ( tran ) THEN
                CALL DMAKE3('GE',' ',' ',ma,na,Ab,2*Nmax,Aa,lda,reset,ZERO)
              ELSE
                CALL DMAKE3('GE',' ',' ',ma,na,Ab,Nmax,Aa,lda,reset,ZERO)
              END IF
              !
              !              Generate the matrix B.
              !
              ldb = lda
              lbb = laa
              IF ( tran ) THEN
                CALL DMAKE3('GE',' ',' ',ma,na,Ab(k+1),2*Nmax,Bb,ldb,reset,ZERO)
              ELSE
                CALL DMAKE3('GE',' ',' ',ma,na,Ab(k*Nmax+1),Nmax,Bb,ldb,reset,&
                  ZERO)
              END IF
              !
              DO icu = 1, 2
                uplo = ichu(icu:icu)
                upper = uplo=='U'
                !
                DO ia = 1, Nalf
                  alpha = Alf(ia)
                  !
                  DO ib = 1, Nbet
                    beta = Bet(ib)
                    !
                    !                       Generate the matrix C.
                    !
                    CALL DMAKE3('SY',uplo,' ',n,n,C,Nmax,Cc,ldc,reset,ZERO)
                    !
                    nc = nc + 1
                    !
                    !                       Save every datum before calling the subroutine.
                    !
                    uplos = uplo
                    transs = trans
                    ns = n
                    ks = k
                    als = alpha
                    DO i = 1, laa
                      As(i) = Aa(i)
                    END DO
                    ldas = lda
                    DO i = 1, lbb
                      Bs(i) = Bb(i)
                    END DO
                    ldbs = ldb
                    bets = beta
                    DO i = 1, lcc
                      Cs(i) = Cc(i)
                    END DO
                    ldcs = ldc
                    !
                    !                       Call the subroutine.
                    !
                    CALL DSYR2K(uplo,trans,n,k,alpha,Aa,lda,Bb,ldb,beta,Cc,ldc)
                    !
                    !                       Check if error-exit was taken incorrectly.
                    !
                    IF ( NUMXER(nerr)/=0 ) THEN
                      IF ( Kprint>=2 ) WRITE (Nout,FMT=99006)
                      Fatal = .TRUE.
                    END IF
                    !
                    !                       See what data changed inside subroutines.
                    !
                    isame(1) = uplos==uplo
                    isame(2) = transs==trans
                    isame(3) = ns==n
                    isame(4) = ks==k
                    isame(5) = als==alpha
                    isame(6) = LDE(As,Aa,laa)
                    isame(7) = ldas==lda
                    isame(8) = LDE(Bs,Bb,lbb)
                    isame(9) = ldbs==ldb
                    isame(10) = bets==beta
                    IF ( nul ) THEN
                      isame(11) = LDE(Cs,Cc,lcc)
                    ELSE
                      isame(11) = LDERES('SY',uplo,n,n,Cs,Cc,ldc)
                    END IF
                    isame(12) = ldcs==ldc
                    !
                    !                       If data was incorrectly changed, report and
                    !                       return.
                    !
                    DO i = 1, nargs
                      IF ( .NOT.isame(i) ) THEN
                        Fatal = .TRUE.
                        IF ( Kprint>=2 ) WRITE (Nout,FMT=99002) i
                      END IF
                    END DO
                    !
                    ftl = .FALSE.
                    IF ( .NOT.nul ) THEN
                      !
                      !                          Check the result column by column.
                      !
                      jjab = 1
                      jc = 1
                      DO j = 1, n
                        IF ( upper ) THEN
                          jj = 1
                          lj = j
                        ELSE
                          jj = j
                          lj = n - j + 1
                        END IF
                        IF ( tran ) THEN
                          DO i = 1, k
                            W(i) = Ab((j-1)*2*Nmax+k+i)
                            W(k+i) = Ab((j-1)*2*Nmax+i)
                          END DO
                          CALL DMMCH('T','N',lj,1,2*k,alpha,Ab(jjab),2*Nmax,W,&
                            2*Nmax,beta,C(jj,j),Nmax,Ct,G,Cc(jc),ldc,&
                            Eps,err,ftl,Nout,.TRUE.,Kprint)
                        ELSE
                          DO i = 1, k
                            W(i) = Ab((k+i-1)*Nmax+j)
                            W(k+i) = Ab((i-1)*Nmax+j)
                          END DO
                          CALL DMMCH('N','N',lj,1,2*k,alpha,Ab(jj),Nmax,W,&
                            2*Nmax,beta,C(jj,j),Nmax,Ct,G,Cc(jc),ldc,&
                            Eps,err,ftl,Nout,.TRUE.,Kprint)
                        END IF
                        IF ( upper ) THEN
                          jc = jc + ldc
                        ELSE
                          jc = jc + ldc + 1
                          IF ( tran ) jjab = jjab + 2*Nmax
                        END IF
                        errmax = MAX(errmax,err)
                      END DO
                    END IF
                    IF ( ftl ) THEN
                      Fatal = .TRUE.
                      IF ( Kprint>=3 ) THEN
                        WRITE (Nout,FMT=99004) Sname
                        WRITE (Nout,FMT=99005) nc, Sname, uplo, trans, n, &
                          k, alpha, lda, ldb, beta, ldc
                      END IF
                    END IF
                    !
                  END DO
                  !
                END DO
                !
              END DO
            END IF
            !
          END DO
          !
        END DO
      END IF
      !
    END DO
    !
    !     Report result.
    !
    IF ( .NOT.Fatal ) THEN
      IF ( Kprint>=3 ) THEN
        IF ( errmax<Thresh ) THEN
          WRITE (Nout,FMT=99001) Sname, nc
        ELSE
          WRITE (Nout,FMT=99003) Sname, nc, errmax
        END IF
      END IF
    END IF
    RETURN
    !
    99001 FORMAT (' ',A6,' PASSED THE COMPUTATIONAL TESTS (',I6,' CALL','S)')
    99002 FORMAT (' ******* FATAL ERROR - PARAMETER NUMBER ',I2,' WAS CH',&
      'ANGED INCORRECTLY *******')
    99003 FORMAT (' ',A6,' COMPLETED THE COMPUTATIONAL TESTS (',I6,' C','ALLS)',&
      /' ******* BUT WITH MAXIMUM TEST RATIO',F8.2,' - SUSPECT *******')
    99004 FORMAT (' ******* ',A6,' FAILED ON CALL NUMBER:')
    99005 FORMAT (1X,I6,': ',A6,'(',2('''',A1,''','),2(I3,','),F4.1,', A,',I3,&
      ', B,',I3,',',F4.1,', C,',I3,')   ',' .')
    99006 FORMAT (' ******* FATAL ERROR - ERROR-EXIT TAKEN ON VALID CALL *',&
      '******')
    !
    !     End of DCHK53.
    !
  END SUBROUTINE DCHK53
  !** DCHK62
  SUBROUTINE DCHK62(Sname,Eps,Thresh,Nout,Kprint,Fatal,Nidim,Idimm,Nalf,Alf,&
      Ninc,Inc,Nmax,Incmax,A,Aa,As,X,Xx,Xs,Y,Yy,Ys,Yt,G,Z)
    !>
    !  Test DSYR2 and DSPR2.
    !***
    ! **Library:**   SLATEC (BLAS)
    !***
    ! **Keywords:**  BLAS, QUICK CHECK SERVICE ROUTINE
    !***
    ! **Author:**  Du Croz, J. (NAG)
    !           Hanson, R. J. (SNLA)
    !***
    ! **Description:**
    !
    !  Quick check for DSYR2 and DSPR2.
    !
    !  Auxiliary routine for test program for Level 2 Blas.
    !***
    ! **References:**  (NONE)
    !***
    ! **Routines called:**  DMAKE2, DMVCH, DSPR2, DSYR2, LDE, LDERES, NUMXER

    !* REVISION HISTORY  (YYMMDD)
    !   870810  DATE WRITTEN
    !   910619  Modified to meet SLATEC code and prologue standards. (BKS)
    USE slatec, ONLY : DSPR2, DSYR2, NUMXER
    !     .. Parameters ..
    REAL(8), PARAMETER :: ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0
    !     .. Scalar Arguments ..
    LOGICAL Fatal
    REAL(8) :: Eps, Thresh
    INTEGER Incmax, Kprint, Nalf, Nidim, Ninc, Nmax, Nout
    CHARACTER(6) :: Sname
    !     .. Array Arguments ..
    REAL(8) :: A(Nmax,Nmax), Aa(Nmax*Nmax), Alf(Nalf), As(Nmax*Nmax), &
      G(Nmax), X(Nmax), Xs(Nmax*Incmax), Xx(Nmax*Incmax), &
      Y(Nmax), Ys(Nmax*Incmax), Yt(Nmax), Yy(Nmax*Incmax), Z(Nmax,2)
    INTEGER Idimm(Nidim), Inc(Ninc)
    !     .. Local Scalars ..
    REAL(8) :: alpha, als, err, errmax, transl
    INTEGER i, ia, ic, in, incx, incxs, incy, incys, ix, iy, j, &
      ja, jj, laa, lda, ldas, lj, lx, ly, n, nargs, nc, ns, nerr
    LOGICAL ftl, full, nul, packed, reset, upper
    CHARACTER :: uplo, uplos
    !     .. Local Arrays ..
    REAL(8) :: w(2)
    LOGICAL isame(13)
    !     .. Intrinsic Functions ..
    INTRINSIC ABS, MAX
    !     .. Data statements ..
    CHARACTER(2), PARAMETER :: ich = 'UL'
    !* FIRST EXECUTABLE STATEMENT  DCHK62
    full = Sname(3:3)=='Y'
    packed = Sname(3:3)=='P'
    !     Define the number of arguments.
    IF ( full ) THEN
      nargs = 9
    ELSEIF ( packed ) THEN
      nargs = 8
    END IF
    !
    nc = 0
    reset = .TRUE.
    errmax = ZERO
    !
    DO in = 1, Nidim
      n = Idimm(in)
      !        Set LDA to 1 more than minimum value if room.
      lda = n
      IF ( lda<Nmax ) lda = lda + 1
      !        Skip tests if not enough room.
      IF ( lda<=Nmax ) THEN
        IF ( packed ) THEN
          laa = (n*(n+1))/2
        ELSE
          laa = lda*n
        END IF
        !
        DO ic = 1, 2
          uplo = ich(ic:ic)
          upper = uplo=='U'
          !
          DO ix = 1, Ninc
            incx = Inc(ix)
            lx = ABS(incx)*n
            !
            !              Generate the vector X.
            !
            transl = HALF
            CALL DMAKE2('GE',' ',' ',1,n,X,1,Xx,ABS(incx),0,n-1,reset,transl)
            IF ( n>1 ) THEN
              X(n/2) = ZERO
              Xx(1+ABS(incx)*(n/2-1)) = ZERO
            END IF
            !
            DO iy = 1, Ninc
              incy = Inc(iy)
              ly = ABS(incy)*n
              !
              !                 Generate the vector Y.
              !
              transl = ZERO
              CALL DMAKE2('GE',' ',' ',1,n,Y,1,Yy,ABS(incy),0,n-1,reset,transl)
              IF ( n>1 ) THEN
                Y(n/2) = ZERO
                Yy(1+ABS(incy)*(n/2-1)) = ZERO
              END IF
              !
              DO ia = 1, Nalf
                alpha = Alf(ia)
                nul = n<=0 .OR. alpha==ZERO
                !
                !                    Generate the matrix A.
                !
                transl = ZERO
                CALL DMAKE2(Sname(2:3),uplo,' ',n,n,A,Nmax,Aa,lda,n-1,n-1,&
                  reset,transl)
                !
                nc = nc + 1
                !
                !                    Save every datum before calling the subroutine.
                !
                uplos = uplo
                ns = n
                als = alpha
                DO i = 1, laa
                  As(i) = Aa(i)
                END DO
                ldas = lda
                DO i = 1, lx
                  Xs(i) = Xx(i)
                END DO
                incxs = incx
                DO i = 1, ly
                  Ys(i) = Yy(i)
                END DO
                incys = incy
                !
                !                    Call the subroutine.
                !
                IF ( full ) THEN
                  CALL DSYR2(uplo,n,alpha,Xx,incx,Yy,incy,Aa,lda)
                ELSEIF ( packed ) THEN
                  CALL DSPR2(uplo,n,alpha,Xx,incx,Yy,incy,Aa)
                END IF
                !
                !                    Check if error-exit was taken incorrectly.
                !
                IF ( NUMXER(nerr)/=0 ) THEN
                  IF ( Kprint>=2 ) WRITE (Nout,FMT=99008)
                  Fatal = .TRUE.
                END IF
                !
                !                    See what data changed inside subroutines.
                !
                isame(1) = uplo==uplos
                isame(2) = ns==n
                isame(3) = als==alpha
                isame(4) = LDE(Xs,Xx,lx)
                isame(5) = incxs==incx
                isame(6) = LDE(Ys,Yy,ly)
                isame(7) = incys==incy
                IF ( nul ) THEN
                  isame(8) = LDE(As,Aa,laa)
                ELSE
                  isame(8) = LDERES(Sname(2:3),uplo,n,n,As,Aa,lda)
                END IF
                IF ( .NOT.packed ) isame(9) = ldas==lda
                !
                !                    If data was incorrectly changed, report and return.
                !
                DO i = 1, nargs
                  IF ( .NOT.isame(i) ) THEN
                    Fatal = .TRUE.
                    IF ( Kprint>=2 ) WRITE (Nout,FMT=99002) i
                  END IF
                END DO
                !
                IF ( .NOT.nul ) THEN
                  !
                  !                       Check the result column by column.
                  !
                  IF ( incx>0 ) THEN
                    DO i = 1, n
                      Z(i,1) = X(i)
                    END DO
                  ELSE
                    DO i = 1, n
                      Z(i,1) = X(n-i+1)
                    END DO
                  END IF
                  IF ( incy>0 ) THEN
                    DO i = 1, n
                      Z(i,2) = Y(i)
                    END DO
                  ELSE
                    DO i = 1, n
                      Z(i,2) = Y(n-i+1)
                    END DO
                  END IF
                  ja = 1
                  DO j = 1, n
                    w(1) = Z(j,2)
                    w(2) = Z(j,1)
                    IF ( upper ) THEN
                      jj = 1
                      lj = j
                    ELSE
                      jj = j
                      lj = n - j + 1
                    END IF
                    ftl = .FALSE.
                    CALL DMVCH('N',lj,2,alpha,Z(jj,1),Nmax,w,1,ONE,A(jj,j),1,&
                      Yt,G,Aa(ja),Eps,err,ftl,Nout,.TRUE.,Kprint)
                    IF ( .NOT.(full) ) THEN
                      ja = ja + lj
                    ELSEIF ( upper ) THEN
                      ja = ja + lda
                    ELSE
                      ja = ja + lda + 1
                    END IF
                    errmax = MAX(errmax,err)
                    IF ( ftl ) THEN
                      Fatal = .TRUE.
                      IF ( Kprint>=3 ) THEN
                        WRITE (Nout,FMT=99005) j
                        WRITE (Nout,FMT=99004) Sname
                        IF ( full ) THEN
                          WRITE (Nout,FMT=99007) nc, Sname, uplo, n, &
                            alpha, incx, incy, lda
                        ELSEIF ( packed ) THEN
                          WRITE (Nout,FMT=99006) nc, Sname, uplo, n, &
                            alpha, incx, incy
                        END IF
                      END IF
                    END IF
                  END DO
                END IF
                !
              END DO
              !
            END DO
            !
          END DO
          !
        END DO
      END IF
      !
    END DO
    !
    !     Report result.
    !
    IF ( .NOT.Fatal ) THEN
      IF ( Kprint>=3 ) THEN
        IF ( errmax<Thresh ) THEN
          WRITE (Nout,FMT=99001) Sname, nc
        ELSE
          WRITE (Nout,FMT=99003) Sname, nc, errmax
        END IF
      END IF
    END IF
    RETURN
    !
    99001 FORMAT (' ',A6,' PASSED THE COMPUTATIONAL TESTS (',I6,' CALL','S)')
    99002 FORMAT (' ******* FATAL ERROR - PARAMETER NUMBER ',I2,' WAS CH',&
      'ANGED INCORRECTLY *******')
    99003 FORMAT (' ',A6,' COMPLETED THE COMPUTATIONAL TESTS (',I6,' C','ALLS)',&
      /' ******* BUT WITH MAXIMUM TEST RATIO',F8.2,' - SUSPECT *******')
    99004 FORMAT (' ******* ',A6,' FAILED ON CALL NUMBER:')
    99005 FORMAT ('      THESE ARE THE RESULTS FOR COLUMN ',I3)
    99006 FORMAT (1X,I6,': ',A6,'(''',A1,''',',I3,',',F4.1,', X,',I2,', Y,',I2,&
      ', AP)                     .')
    99007 FORMAT (1X,I6,': ',A6,'(''',A1,''',',I3,',',F4.1,', X,',I2,', Y,',I2,&
      ', A,',I3,')                  .')
    99008 FORMAT (' ******* FATAL ERROR - ERROR-EXIT TAKEN ON VALID CALL *',&
      '******')
    !
    !     End of DCHK62.
    !
  END SUBROUTINE DCHK62
  !** DCHKE2
  SUBROUTINE DCHKE2(Isnum,Srnamt,Nout,Kprint,Fatal)
    !>
    !  Test the error exits from the Level 2 Blas.
    !***
    ! **Library:**   SLATEC (BLAS)
    !***
    ! **Author:**  Du Croz, J. J., (NAG)
    !           Hanson, R. J., (SNLA)
    !***
    ! **Description:**
    !
    !  Tests the error exits from the Level 2 Blas.
    !  ALPHA, BETA, A, X and Y should not need to be defined.
    !
    !  Auxiliary routine for test program for Level 2 Blas.
    !***
    ! **References:**  (NONE)
    !***
    ! **Routines called:**  CHKXER, DGBMV, DGEMV, DGER, DSBMV, DSPMV, DSPR,
    !                    DSPR2, DSYMV, DSYR, DSYR2, DTBMV, DTBSV, DTPMV,
    !                    DTPSV, DTRMV, DTRSV, XERCLR, XERDMP, XGETF, XSETF

    !* REVISION HISTORY  (YYMMDD)
    !   870810  DATE WRITTEN
    !   910620  Modified to meet SLATEC code and prologue standards.  (BKS)
    USE slatec, ONLY : DGBMV, DGEMV, DGER, DSBMV, DSPMV, DSPR, DSPR2, DSYMV, DSYR, &
      DSYR2, DTBMV, DTBSV, DTPMV, DTPSV, DTRMV, DTRSV, XERCLR, XERDMP, XGETF, XSETF
    USE common_mod, ONLY : CHKXER
    !     .. Scalar Arguments ..
    INTEGER Isnum, Nout
    LOGICAL Fatal
    CHARACTER(6) :: Srnamt
    INTEGER infot, Kprint
    !     .. Local Scalars ..
    REAL(8) :: alpha, beta
    INTEGER kontrl
    !     .. Local Arrays ..
    REAL(8) :: a(1,1), x(1), y(1)
    !* FIRST EXECUTABLE STATEMENT  DCHKE2
    CALL XGETF(kontrl)
    IF ( Kprint<=2 ) THEN
      CALL XSETF(0)
    ELSE
      CALL XSETF(1)
    END IF
    SELECT CASE (Isnum)
      CASE (2)
        infot = 1
        CALL XERCLR
        CALL DGBMV('/',0,0,0,0,alpha,a,1,x,1,beta,y,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 2
        CALL XERCLR
        CALL DGBMV('N',-1,0,0,0,alpha,a,1,x,1,beta,y,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 3
        CALL XERCLR
        CALL DGBMV('N',0,-1,0,0,alpha,a,1,x,1,beta,y,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 4
        CALL XERCLR
        CALL DGBMV('N',0,0,-1,0,alpha,a,1,x,1,beta,y,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 5
        CALL XERCLR
        CALL DGBMV('N',2,0,0,-1,alpha,a,1,x,1,beta,y,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 8
        CALL XERCLR
        CALL DGBMV('N',0,0,1,0,alpha,a,1,x,1,beta,y,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 10
        CALL XERCLR
        CALL DGBMV('N',0,0,0,0,alpha,a,1,x,0,beta,y,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 13
        CALL XERCLR
        CALL DGBMV('N',0,0,0,0,alpha,a,1,x,1,beta,y,0)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      CASE (3)
        infot = 1
        CALL XERCLR
        CALL DSYMV('/',0,alpha,a,1,x,1,beta,y,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 2
        CALL XERCLR
        CALL DSYMV('U',-1,alpha,a,1,x,1,beta,y,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 5
        CALL XERCLR
        CALL DSYMV('U',2,alpha,a,1,x,1,beta,y,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 7
        CALL XERCLR
        CALL DSYMV('U',0,alpha,a,1,x,0,beta,y,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 10
        CALL XERCLR
        CALL DSYMV('U',0,alpha,a,1,x,1,beta,y,0)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      CASE (4)
        infot = 1
        CALL XERCLR
        CALL DSBMV('/',0,0,alpha,a,1,x,1,beta,y,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 2
        CALL XERCLR
        CALL DSBMV('U',-1,0,alpha,a,1,x,1,beta,y,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 3
        CALL XERCLR
        CALL DSBMV('U',0,-1,alpha,a,1,x,1,beta,y,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 6
        CALL XERCLR
        CALL DSBMV('U',0,1,alpha,a,1,x,1,beta,y,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 8
        CALL XERCLR
        CALL DSBMV('U',0,0,alpha,a,1,x,0,beta,y,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 11
        CALL XERCLR
        CALL DSBMV('U',0,0,alpha,a,1,x,1,beta,y,0)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      CASE (5)
        infot = 1
        CALL XERCLR
        CALL DSPMV('/',0,alpha,a,x,1,beta,y,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 2
        CALL XERCLR
        CALL DSPMV('U',-1,alpha,a,x,1,beta,y,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 6
        CALL XERCLR
        CALL DSPMV('U',0,alpha,a,x,0,beta,y,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 9
        CALL XERCLR
        CALL DSPMV('U',0,alpha,a,x,1,beta,y,0)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      CASE (6)
        infot = 1
        CALL XERCLR
        CALL DTRMV('/','N','N',0,a,1,x,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 2
        CALL XERCLR
        CALL DTRMV('U','/','N',0,a,1,x,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 3
        CALL XERCLR
        CALL DTRMV('U','N','/',0,a,1,x,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 4
        CALL XERCLR
        CALL DTRMV('U','N','N',-1,a,1,x,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 6
        CALL XERCLR
        CALL DTRMV('U','N','N',2,a,1,x,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 8
        CALL XERCLR
        CALL DTRMV('U','N','N',0,a,1,x,0)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      CASE (7)
        infot = 1
        CALL XERCLR
        CALL DTBMV('/','N','N',0,0,a,1,x,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 2
        CALL XERCLR
        CALL DTBMV('U','/','N',0,0,a,1,x,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 3
        CALL XERCLR
        CALL DTBMV('U','N','/',0,0,a,1,x,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 4
        CALL XERCLR
        CALL DTBMV('U','N','N',-1,0,a,1,x,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 5
        CALL XERCLR
        CALL DTBMV('U','N','N',0,-1,a,1,x,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 7
        CALL XERCLR
        CALL DTBMV('U','N','N',0,1,a,1,x,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 9
        CALL XERCLR
        CALL DTBMV('U','N','N',0,0,a,1,x,0)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      CASE (8)
        infot = 1
        CALL XERCLR
        CALL DTPMV('/','N','N',0,a,x,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 2
        CALL XERCLR
        CALL DTPMV('U','/','N',0,a,x,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 3
        CALL XERCLR
        CALL DTPMV('U','N','/',0,a,x,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 4
        CALL XERCLR
        CALL DTPMV('U','N','N',-1,a,x,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 7
        CALL XERCLR
        CALL DTPMV('U','N','N',0,a,x,0)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      CASE (9)
        infot = 1
        CALL XERCLR
        CALL DTRSV('/','N','N',0,a,1,x,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 2
        CALL XERCLR
        CALL DTRSV('U','/','N',0,a,1,x,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 3
        CALL XERCLR
        CALL DTRSV('U','N','/',0,a,1,x,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 4
        CALL XERCLR
        CALL DTRSV('U','N','N',-1,a,1,x,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 6
        CALL XERCLR
        CALL DTRSV('U','N','N',2,a,1,x,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 8
        CALL XERCLR
        CALL DTRSV('U','N','N',0,a,1,x,0)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      CASE (10)
        infot = 1
        CALL XERCLR
        CALL DTBSV('/','N','N',0,0,a,1,x,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 2
        CALL XERCLR
        CALL DTBSV('U','/','N',0,0,a,1,x,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 3
        CALL XERCLR
        CALL DTBSV('U','N','/',0,0,a,1,x,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 4
        CALL XERCLR
        CALL DTBSV('U','N','N',-1,0,a,1,x,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 5
        CALL XERCLR
        CALL DTBSV('U','N','N',0,-1,a,1,x,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 7
        CALL XERCLR
        CALL DTBSV('U','N','N',0,1,a,1,x,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 9
        CALL XERCLR
        CALL DTBSV('U','N','N',0,0,a,1,x,0)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      CASE (11)
        infot = 1
        CALL XERCLR
        CALL DTPSV('/','N','N',0,a,x,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 2
        CALL XERCLR
        CALL DTPSV('U','/','N',0,a,x,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 3
        CALL XERCLR
        CALL DTPSV('U','N','/',0,a,x,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 4
        CALL XERCLR
        CALL DTPSV('U','N','N',-1,a,x,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 7
        CALL XERCLR
        CALL DTPSV('U','N','N',0,a,x,0)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      CASE (12)
        infot = 1
        CALL XERCLR
        CALL DGER(-1,0,alpha,x,1,y,1,a,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 2
        CALL XERCLR
        CALL DGER(0,-1,alpha,x,1,y,1,a,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 5
        CALL XERCLR
        CALL DGER(0,0,alpha,x,0,y,1,a,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 7
        CALL XERCLR
        CALL DGER(0,0,alpha,x,1,y,0,a,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 9
        CALL XERCLR
        CALL DGER(2,0,alpha,x,1,y,1,a,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      CASE (13)
        infot = 1
        CALL XERCLR
        CALL DSYR('/',0,alpha,x,1,a,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 2
        CALL XERCLR
        CALL DSYR('U',-1,alpha,x,1,a,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 5
        CALL XERCLR
        CALL DSYR('U',0,alpha,x,0,a,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 7
        CALL XERCLR
        CALL DSYR('U',2,alpha,x,1,a,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      CASE (14)
        infot = 1
        CALL XERCLR
        CALL DSPR('/',0,alpha,x,1,a)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 2
        CALL XERCLR
        CALL DSPR('U',-1,alpha,x,1,a)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 5
        CALL XERCLR
        CALL DSPR('U',0,alpha,x,0,a)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      CASE (15)
        infot = 1
        CALL XERCLR
        CALL DSYR2('/',0,alpha,x,1,y,1,a,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 2
        CALL XERCLR
        CALL DSYR2('U',-1,alpha,x,1,y,1,a,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 5
        CALL XERCLR
        CALL DSYR2('U',0,alpha,x,0,y,1,a,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 7
        CALL XERCLR
        CALL DSYR2('U',0,alpha,x,1,y,0,a,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 9
        CALL XERCLR
        CALL DSYR2('U',2,alpha,x,1,y,1,a,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      CASE (16)
        infot = 1
        CALL XERCLR
        CALL DSPR2('/',0,alpha,x,1,y,1,a)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 2
        CALL XERCLR
        CALL DSPR2('U',-1,alpha,x,1,y,1,a)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 5
        CALL XERCLR
        CALL DSPR2('U',0,alpha,x,0,y,1,a)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 7
        CALL XERCLR
        CALL DSPR2('U',0,alpha,x,1,y,0,a)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      CASE DEFAULT
        infot = 1
        CALL XERCLR
        CALL DGEMV('/',0,0,alpha,a,1,x,1,beta,y,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 2
        CALL XERCLR
        CALL DGEMV('N',-1,0,alpha,a,1,x,1,beta,y,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 3
        CALL XERCLR
        CALL DGEMV('N',0,-1,alpha,a,1,x,1,beta,y,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 6
        CALL XERCLR
        CALL DGEMV('N',2,0,alpha,a,1,x,1,beta,y,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 8
        CALL XERCLR
        CALL DGEMV('N',0,0,alpha,a,1,x,0,beta,y,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 11
        CALL XERCLR
        CALL DGEMV('N',0,0,alpha,a,1,x,1,beta,y,0)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
    END SELECT
    !
    IF ( Kprint>=2 ) THEN
      CALL XERDMP
      IF ( .NOT.Fatal ) THEN
        WRITE (Nout,FMT=99001) Srnamt
      ELSE
        WRITE (Nout,FMT=99002) Srnamt
      END IF
    END IF
    CALL XSETF(kontrl)
    RETURN
    !
    99001 FORMAT (' ',A6,' PASSED THE TESTS OF ERROR-EXITS')
    99002 FORMAT (' ******* ',A6,' FAILED THE TESTS OF ERROR-EXITS *****','**')
    !
    !     End of DCHKE2.
    !
  END SUBROUTINE DCHKE2
  !** DCHKE3
  SUBROUTINE DCHKE3(Isnum,Srnamt,Nout,Kprint,Fatal)
    !>
    !  Test the error exits from the Level 3 Blas.
    !***
    ! **Library:**   SLATEC (BLAS)
    !***
    ! **Author:**  Dongarra, J. J., (ANL)
    !           Duff, I., (AERE)
    !           Du Croz, J., (NAG)
    !           Hammarling, S., (NAG)
    !***
    ! **Description:**
    !
    !  Tests the error exits from the Level 3 Blas.
    !  ALPHA, BETA, A, X and Y should not need to be defined.
    !
    !  Auxiliary routine for test program for Level 3 Blas.
    !***
    ! **References:**  (NONE)
    !***
    ! **Routines called:**  CHKXER, DGEMM, DSYMM, DSYR2K, DSYRK, DTRMM, DTRSM,
    !                    XERCLR, XERDMP, XGETF, XSETF

    !* REVISION HISTORY  (YYMMDD)
    !   890208  DATE WRITTEN
    !   910620  Modified to meet SLATEC code and prologue standards.  (BKS)
    USE slatec, ONLY : DGEMM, DSYMM, DSYR2K, DSYRK, DTRMM, DTRSM, XERCLR, XERDMP, &
      XGETF, XSETF
    USE common_mod, ONLY : CHKXER
    !     .. Scalar Arguments ..
    LOGICAL Fatal
    INTEGER Isnum, Nout
    CHARACTER(6) :: Srnamt
    INTEGER infot, Kprint
    !     .. Local Scalars ..
    REAL(8) :: alpha, beta
    INTEGER kontrl
    !     .. Local Arrays ..
    REAL(8) :: a(1,1), b(1,1), c(1,1)
    !* FIRST EXECUTABLE STATEMENT  DCHKE3
    CALL XGETF(kontrl)
    IF ( Kprint<=2 ) THEN
      CALL XSETF(0)
    ELSE
      CALL XSETF(1)
    END IF
    SELECT CASE (Isnum)
      CASE (2)
        infot = 1
        CALL XERCLR
        CALL DSYMM('/','U',0,0,alpha,a,1,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 2
        CALL XERCLR
        CALL DSYMM('L','/',0,0,alpha,a,1,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 3
        CALL XERCLR
        CALL DSYMM('L','U',-1,0,alpha,a,1,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 3
        CALL XERCLR
        CALL DSYMM('R','U',-1,0,alpha,a,1,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 3
        CALL XERCLR
        CALL DSYMM('L','L',-1,0,alpha,a,1,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 3
        CALL XERCLR
        CALL DSYMM('R','L',-1,0,alpha,a,1,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 4
        CALL XERCLR
        CALL DSYMM('L','U',0,-1,alpha,a,1,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 4
        CALL XERCLR
        CALL DSYMM('R','U',0,-1,alpha,a,1,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 4
        CALL XERCLR
        CALL DSYMM('L','L',0,-1,alpha,a,1,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 4
        CALL XERCLR
        CALL DSYMM('R','L',0,-1,alpha,a,1,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 7
        CALL XERCLR
        CALL DSYMM('L','U',2,0,alpha,a,1,b,2,beta,c,2)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 7
        CALL XERCLR
        CALL DSYMM('R','U',0,2,alpha,a,1,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 7
        CALL XERCLR
        CALL DSYMM('L','L',2,0,alpha,a,1,b,2,beta,c,2)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 7
        CALL XERCLR
        CALL DSYMM('R','L',0,2,alpha,a,1,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 9
        CALL XERCLR
        CALL DSYMM('L','U',2,0,alpha,a,2,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 9
        CALL XERCLR
        CALL DSYMM('R','U',2,0,alpha,a,1,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 9
        CALL XERCLR
        CALL DSYMM('L','L',2,0,alpha,a,2,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 9
        CALL XERCLR
        CALL DSYMM('R','L',2,0,alpha,a,1,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 12
        CALL XERCLR
        CALL DSYMM('L','U',2,0,alpha,a,2,b,2,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 12
        CALL XERCLR
        CALL DSYMM('R','U',2,0,alpha,a,1,b,2,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 12
        CALL XERCLR
        CALL DSYMM('L','L',2,0,alpha,a,2,b,2,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 12
        CALL XERCLR
        CALL DSYMM('R','L',2,0,alpha,a,1,b,2,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      CASE (3)
        infot = 1
        CALL XERCLR
        CALL DTRMM('/','U','N','N',0,0,alpha,a,1,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 2
        CALL XERCLR
        CALL DTRMM('L','/','N','N',0,0,alpha,a,1,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 3
        CALL XERCLR
        CALL DTRMM('L','U','/','N',0,0,alpha,a,1,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 4
        CALL XERCLR
        CALL DTRMM('L','U','N','/',0,0,alpha,a,1,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 5
        CALL XERCLR
        CALL DTRMM('L','U','N','N',-1,0,alpha,a,1,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 5
        CALL XERCLR
        CALL DTRMM('L','U','T','N',-1,0,alpha,a,1,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 5
        CALL XERCLR
        CALL DTRMM('R','U','N','N',-1,0,alpha,a,1,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 5
        CALL XERCLR
        CALL DTRMM('R','U','T','N',-1,0,alpha,a,1,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 5
        CALL XERCLR
        CALL DTRMM('L','L','N','N',-1,0,alpha,a,1,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 5
        CALL XERCLR
        CALL DTRMM('L','L','T','N',-1,0,alpha,a,1,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 5
        CALL XERCLR
        CALL DTRMM('R','L','N','N',-1,0,alpha,a,1,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 5
        CALL XERCLR
        CALL DTRMM('R','L','T','N',-1,0,alpha,a,1,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 6
        CALL XERCLR
        CALL DTRMM('L','U','N','N',0,-1,alpha,a,1,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 6
        CALL XERCLR
        CALL DTRMM('L','U','T','N',0,-1,alpha,a,1,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 6
        CALL XERCLR
        CALL DTRMM('R','U','N','N',0,-1,alpha,a,1,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 6
        CALL XERCLR
        CALL DTRMM('R','U','T','N',0,-1,alpha,a,1,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 6
        CALL XERCLR
        CALL DTRMM('L','L','N','N',0,-1,alpha,a,1,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 6
        CALL XERCLR
        CALL DTRMM('L','L','T','N',0,-1,alpha,a,1,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 6
        CALL XERCLR
        CALL DTRMM('R','L','N','N',0,-1,alpha,a,1,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 6
        CALL XERCLR
        CALL DTRMM('R','L','T','N',0,-1,alpha,a,1,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 9
        CALL XERCLR
        CALL DTRMM('L','U','N','N',2,0,alpha,a,1,b,2)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 9
        CALL XERCLR
        CALL DTRMM('L','U','T','N',2,0,alpha,a,1,b,2)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 9
        CALL XERCLR
        CALL DTRMM('R','U','N','N',0,2,alpha,a,1,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 9
        CALL XERCLR
        CALL DTRMM('R','U','T','N',0,2,alpha,a,1,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 9
        CALL XERCLR
        CALL DTRMM('L','L','N','N',2,0,alpha,a,1,b,2)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 9
        CALL XERCLR
        CALL DTRMM('L','L','T','N',2,0,alpha,a,1,b,2)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 9
        CALL XERCLR
        CALL DTRMM('R','L','N','N',0,2,alpha,a,1,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 9
        CALL XERCLR
        CALL DTRMM('R','L','T','N',0,2,alpha,a,1,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 11
        CALL XERCLR
        CALL DTRMM('L','U','N','N',2,0,alpha,a,2,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 11
        CALL XERCLR
        CALL DTRMM('L','U','T','N',2,0,alpha,a,2,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 11
        CALL XERCLR
        CALL DTRMM('R','U','N','N',2,0,alpha,a,1,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 11
        CALL XERCLR
        CALL DTRMM('R','U','T','N',2,0,alpha,a,1,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 11
        CALL XERCLR
        CALL DTRMM('L','L','N','N',2,0,alpha,a,2,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 11
        CALL XERCLR
        CALL DTRMM('L','L','T','N',2,0,alpha,a,2,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 11
        CALL XERCLR
        CALL DTRMM('R','L','N','N',2,0,alpha,a,1,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 11
        CALL XERCLR
        CALL DTRMM('R','L','T','N',2,0,alpha,a,1,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      CASE (4)
        infot = 1
        CALL XERCLR
        CALL DTRSM('/','U','N','N',0,0,alpha,a,1,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 2
        CALL XERCLR
        CALL DTRSM('L','/','N','N',0,0,alpha,a,1,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 3
        CALL XERCLR
        CALL DTRSM('L','U','/','N',0,0,alpha,a,1,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 4
        CALL XERCLR
        CALL DTRSM('L','U','N','/',0,0,alpha,a,1,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 5
        CALL XERCLR
        CALL DTRSM('L','U','N','N',-1,0,alpha,a,1,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 5
        CALL XERCLR
        CALL DTRSM('L','U','T','N',-1,0,alpha,a,1,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 5
        CALL XERCLR
        CALL DTRSM('R','U','N','N',-1,0,alpha,a,1,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 5
        CALL XERCLR
        CALL DTRSM('R','U','T','N',-1,0,alpha,a,1,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 5
        CALL XERCLR
        CALL DTRSM('L','L','N','N',-1,0,alpha,a,1,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 5
        CALL XERCLR
        CALL DTRSM('L','L','T','N',-1,0,alpha,a,1,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 5
        CALL XERCLR
        CALL DTRSM('R','L','N','N',-1,0,alpha,a,1,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 5
        CALL XERCLR
        CALL DTRSM('R','L','T','N',-1,0,alpha,a,1,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 6
        CALL XERCLR
        CALL DTRSM('L','U','N','N',0,-1,alpha,a,1,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 6
        CALL XERCLR
        CALL DTRSM('L','U','T','N',0,-1,alpha,a,1,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 6
        CALL XERCLR
        CALL DTRSM('R','U','N','N',0,-1,alpha,a,1,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 6
        CALL XERCLR
        CALL DTRSM('R','U','T','N',0,-1,alpha,a,1,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 6
        CALL XERCLR
        CALL DTRSM('L','L','N','N',0,-1,alpha,a,1,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 6
        CALL XERCLR
        CALL DTRSM('L','L','T','N',0,-1,alpha,a,1,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 6
        CALL XERCLR
        CALL DTRSM('R','L','N','N',0,-1,alpha,a,1,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 6
        CALL XERCLR
        CALL DTRSM('R','L','T','N',0,-1,alpha,a,1,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 9
        CALL XERCLR
        CALL DTRSM('L','U','N','N',2,0,alpha,a,1,b,2)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 9
        CALL XERCLR
        CALL DTRSM('L','U','T','N',2,0,alpha,a,1,b,2)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 9
        CALL XERCLR
        CALL DTRSM('R','U','N','N',0,2,alpha,a,1,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 9
        CALL XERCLR
        CALL DTRSM('R','U','T','N',0,2,alpha,a,1,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 9
        CALL XERCLR
        CALL DTRSM('L','L','N','N',2,0,alpha,a,1,b,2)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 9
        CALL XERCLR
        CALL DTRSM('L','L','T','N',2,0,alpha,a,1,b,2)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 9
        CALL XERCLR
        CALL DTRSM('R','L','N','N',0,2,alpha,a,1,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 9
        CALL XERCLR
        CALL DTRSM('R','L','T','N',0,2,alpha,a,1,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 11
        CALL XERCLR
        CALL DTRSM('L','U','N','N',2,0,alpha,a,2,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 11
        CALL XERCLR
        CALL DTRSM('L','U','T','N',2,0,alpha,a,2,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 11
        CALL XERCLR
        CALL DTRSM('R','U','N','N',2,0,alpha,a,1,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 11
        CALL XERCLR
        CALL DTRSM('R','U','T','N',2,0,alpha,a,1,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 11
        CALL XERCLR
        CALL DTRSM('L','L','N','N',2,0,alpha,a,2,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 11
        CALL XERCLR
        CALL DTRSM('L','L','T','N',2,0,alpha,a,2,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 11
        CALL XERCLR
        CALL DTRSM('R','L','N','N',2,0,alpha,a,1,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 11
        CALL XERCLR
        CALL DTRSM('R','L','T','N',2,0,alpha,a,1,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      CASE (5)
        infot = 1
        CALL XERCLR
        CALL DSYRK('/','N',0,0,alpha,a,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 2
        CALL XERCLR
        CALL DSYRK('U','/',0,0,alpha,a,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 3
        CALL XERCLR
        CALL DSYRK('U','N',-1,0,alpha,a,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 3
        CALL XERCLR
        CALL DSYRK('U','T',-1,0,alpha,a,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 3
        CALL XERCLR
        CALL DSYRK('L','N',-1,0,alpha,a,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 3
        CALL XERCLR
        CALL DSYRK('L','T',-1,0,alpha,a,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 4
        CALL XERCLR
        CALL DSYRK('U','N',0,-1,alpha,a,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 4
        CALL XERCLR
        CALL DSYRK('U','T',0,-1,alpha,a,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 4
        CALL XERCLR
        CALL DSYRK('L','N',0,-1,alpha,a,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 4
        CALL XERCLR
        CALL DSYRK('L','T',0,-1,alpha,a,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 7
        CALL XERCLR
        CALL DSYRK('U','N',2,0,alpha,a,1,beta,c,2)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 7
        CALL XERCLR
        CALL DSYRK('U','T',0,2,alpha,a,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 7
        CALL XERCLR
        CALL DSYRK('L','N',2,0,alpha,a,1,beta,c,2)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 7
        CALL XERCLR
        CALL DSYRK('L','T',0,2,alpha,a,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 10
        CALL XERCLR
        CALL DSYRK('U','N',2,0,alpha,a,2,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 10
        CALL XERCLR
        CALL DSYRK('U','T',2,0,alpha,a,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 10
        CALL XERCLR
        CALL DSYRK('L','N',2,0,alpha,a,2,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 10
        CALL XERCLR
        CALL DSYRK('L','T',2,0,alpha,a,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      CASE (6)
        infot = 1
        CALL XERCLR
        CALL DSYR2K('/','N',0,0,alpha,a,1,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 2
        CALL XERCLR
        CALL DSYR2K('U','/',0,0,alpha,a,1,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 3
        CALL XERCLR
        CALL DSYR2K('U','N',-1,0,alpha,a,1,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 3
        CALL XERCLR
        CALL DSYR2K('U','T',-1,0,alpha,a,1,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 3
        CALL XERCLR
        CALL DSYR2K('L','N',-1,0,alpha,a,1,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 3
        CALL XERCLR
        CALL DSYR2K('L','T',-1,0,alpha,a,1,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 4
        CALL XERCLR
        CALL DSYR2K('U','N',0,-1,alpha,a,1,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 4
        CALL XERCLR
        CALL DSYR2K('U','T',0,-1,alpha,a,1,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 4
        CALL XERCLR
        CALL DSYR2K('L','N',0,-1,alpha,a,1,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 4
        CALL XERCLR
        CALL DSYR2K('L','T',0,-1,alpha,a,1,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 7
        CALL XERCLR
        CALL DSYR2K('U','N',2,0,alpha,a,1,b,1,beta,c,2)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 7
        CALL XERCLR
        CALL DSYR2K('U','T',0,2,alpha,a,1,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 7
        CALL XERCLR
        CALL DSYR2K('L','N',2,0,alpha,a,1,b,1,beta,c,2)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 7
        CALL XERCLR
        CALL DSYR2K('L','T',0,2,alpha,a,1,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 9
        CALL XERCLR
        CALL DSYR2K('U','N',2,0,alpha,a,2,b,1,beta,c,2)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 9
        CALL XERCLR
        CALL DSYR2K('U','T',0,2,alpha,a,2,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 9
        CALL XERCLR
        CALL DSYR2K('L','N',2,0,alpha,a,2,b,1,beta,c,2)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 9
        CALL XERCLR
        CALL DSYR2K('L','T',0,2,alpha,a,2,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 12
        CALL XERCLR
        CALL DSYR2K('U','N',2,0,alpha,a,2,b,2,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 12
        CALL XERCLR
        CALL DSYR2K('U','T',2,0,alpha,a,1,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 12
        CALL XERCLR
        CALL DSYR2K('L','N',2,0,alpha,a,2,b,2,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 12
        CALL XERCLR
        CALL DSYR2K('L','T',2,0,alpha,a,1,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      CASE DEFAULT
        infot = 1
        CALL XERCLR
        CALL DGEMM('/','N',0,0,0,alpha,a,1,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 1
        CALL XERCLR
        CALL DGEMM('/','T',0,0,0,alpha,a,1,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 2
        CALL XERCLR
        CALL DGEMM('N','/',0,0,0,alpha,a,1,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 2
        CALL XERCLR
        CALL DGEMM('T','/',0,0,0,alpha,a,1,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 3
        CALL XERCLR
        CALL DGEMM('N','N',-1,0,0,alpha,a,1,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 3
        CALL XERCLR
        CALL DGEMM('N','T',-1,0,0,alpha,a,1,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 3
        CALL XERCLR
        CALL DGEMM('T','N',-1,0,0,alpha,a,1,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 3
        CALL XERCLR
        CALL DGEMM('T','T',-1,0,0,alpha,a,1,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 4
        CALL XERCLR
        CALL DGEMM('N','N',0,-1,0,alpha,a,1,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 4
        CALL XERCLR
        CALL DGEMM('N','T',0,-1,0,alpha,a,1,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 4
        CALL XERCLR
        CALL DGEMM('T','N',0,-1,0,alpha,a,1,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 4
        CALL XERCLR
        CALL DGEMM('T','T',0,-1,0,alpha,a,1,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 5
        CALL XERCLR
        CALL DGEMM('N','N',0,0,-1,alpha,a,1,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 5
        CALL XERCLR
        CALL DGEMM('N','T',0,0,-1,alpha,a,1,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 5
        CALL XERCLR
        CALL DGEMM('T','N',0,0,-1,alpha,a,1,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 5
        CALL XERCLR
        CALL DGEMM('T','T',0,0,-1,alpha,a,1,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 8
        CALL XERCLR
        CALL DGEMM('N','N',2,0,0,alpha,a,1,b,1,beta,c,2)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 8
        CALL XERCLR
        CALL DGEMM('N','T',2,0,0,alpha,a,1,b,1,beta,c,2)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 8
        CALL XERCLR
        CALL DGEMM('T','N',0,0,2,alpha,a,1,b,2,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 8
        CALL XERCLR
        CALL DGEMM('T','T',0,0,2,alpha,a,1,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 10
        CALL XERCLR
        CALL DGEMM('N','N',0,0,2,alpha,a,1,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 10
        CALL XERCLR
        CALL DGEMM('T','N',0,0,2,alpha,a,2,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 10
        CALL XERCLR
        CALL DGEMM('N','T',0,2,0,alpha,a,1,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 10
        CALL XERCLR
        CALL DGEMM('T','T',0,2,0,alpha,a,1,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 13
        CALL XERCLR
        CALL DGEMM('N','N',2,0,0,alpha,a,2,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 13
        CALL XERCLR
        CALL DGEMM('N','T',2,0,0,alpha,a,2,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 13
        CALL XERCLR
        CALL DGEMM('T','N',2,0,0,alpha,a,1,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 13
        CALL XERCLR
        CALL DGEMM('T','T',2,0,0,alpha,a,1,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
    END SELECT
    !
    IF ( Kprint>=2 ) THEN
      CALL XERDMP
      IF ( .NOT.Fatal ) THEN
        WRITE (Nout,FMT=99001) Srnamt
      ELSE
        WRITE (Nout,FMT=99002) Srnamt
      END IF
    END IF
    CALL XSETF(kontrl)
    RETURN
    !
    99001 FORMAT (' ',A6,' PASSED THE TESTS OF ERROR-EXITS')
    99002 FORMAT (' ******* ',A6,' FAILED THE TESTS OF ERROR-EXITS *****','**')
    !
    !     End of DCHKE3.
    !
  END SUBROUTINE DCHKE3
END MODULE TEST19_MOD
!** TEST19
PROGRAM TEST19
  USE TEST19_MOD, ONLY : DBLAT2, DBLAT3
  USE slatec, ONLY : I1MACH, XSETF, XSETUN, XERMAX
  USE common_mod, ONLY : GET_ARGUMENT
  IMPLICIT NONE
  !>
  !  Driver for testing SLATEC subprograms
  !***
  ! **Library:**   SLATEC
  !***
  ! **Category:**  D1B
  !***
  ! **Keywords:**  QUICK CHECK DRIVER
  !***
  ! **Type:**      DOUBLE PRECISION (TEST18-S, TEST19-D, TEST20-C)
  !***
  ! **Author:**  SLATEC Common Mathematical Library Committee
  !***
  ! **Description:**
  !
  !- Usage:
  !     One input data record is required
  !         READ (LIN, '(I1)') KPRINT
  !
  !- Arguments:
  !     KPRINT = 0  Quick checks - No printing.
  !                 Driver       - Short pass or fail message printed.
  !              1  Quick checks - No message printed for passed tests,
  !                                short message printed for failed tests.
  !                 Driver       - Short pass or fail message printed.
  !              2  Quick checks - Print short message for passed tests,
  !                                fuller information for failed tests.
  !                 Driver       - Pass or fail message printed.
  !              3  Quick checks - Print complete quick check results.
  !                 Driver       - Pass or fail message printed.
  !
  !- Description:
  !     Driver for testing SLATEC subprograms
  !        double precision Levels 2 and 3 BLAS routines
  !
  !***
  ! **References:**  Kirby W. Fong,  Thomas H. Jefferson, Tokihiko Suyehiro
  !                 and Lee Walton, Guide to the SLATEC Common Mathema-
  !                 tical Library, April 10, 1990.
  !***
  ! **Routines called:**  I1MACH, DBLAT2, DBLAT3, XERMAX, XSETF, XSETUN

  !* REVISION HISTORY  (YYMMDD)
  !   920601  DATE WRITTEN

  INTEGER ipass, kprint, lin, lun, nfail
  !* FIRST EXECUTABLE STATEMENT  TEST19
  lun = I1MACH(2)
  lin = I1MACH(1)
  nfail = 0
  !
  !     Read KPRINT parameter
  !
  CALL GET_ARGUMENT(kprint)
  CALL XERMAX(1000)
  CALL XSETUN(lun)
  IF ( kprint<=1 ) THEN
    CALL XSETF(0)
  ELSE
    CALL XSETF(1)
  END IF
  !
  !     Test double precision Level 2 BLAS routines
  !
  CALL DBLAT2(lun,kprint,ipass)
  IF ( ipass==0 ) nfail = nfail + 1
  !
  !     Test double precision Level 3 BLAS routines
  !
  CALL DBLAT3(lun,kprint,ipass)
  IF ( ipass==0 ) nfail = nfail + 1
  !
  !     Write PASS or FAIL message
  !
  IF ( nfail==0 ) THEN
    WRITE (lun,99001)
    99001 FORMAT (/' --------------TEST19 PASSED ALL TESTS----------------')
  ELSE
    WRITE (lun,99002) nfail
    99002 FORMAT (/' ************* WARNING -- ',I5,&
      ' TEST(S) FAILED IN PROGRAM TEST19 *************')
  END IF
  STOP
END PROGRAM TEST19
