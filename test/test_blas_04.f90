MODULE TEST20_MOD
  IMPLICIT NONE

CONTAINS
  !** CBEG
  COMPLEX FUNCTION CBEG(Reset)
    IMPLICIT NONE
    !>
    !***
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
    INTEGER i, ic, j, mi, mj
    !     .. Save statement ..
    SAVE i, ic, j, mi, mj
    !     .. Intrinsic Functions ..
    INTRINSIC CMPLX
    !* FIRST EXECUTABLE STATEMENT  CBEG
    IF ( Reset ) THEN
      !        Initialize local variables.
      mi = 891
      mj = 457
      i = 7
      j = 7
      ic = 0
      Reset = .FALSE.
    ENDIF
    !
    !     The sequence of values of I or J is bounded between 1 and 999.
    !     If initial I or J = 1,2,3,6,7 or 9, the period will be 50.
    !     If initial I or J = 4 or 8, the period will be 25.
    !     If initial I or J = 5, the period will be 10.
    !     IC is used to break up the period by skipping 1 value of I or J
    !     in 6.
    !
    ic = ic + 1
    DO
      i = i*mi
      j = j*mj
      i = i - 1000*(i/1000)
      j = j - 1000*(j/1000)
      IF ( ic>=5 ) THEN
        ic = 0
        CYCLE
      ENDIF
      CBEG = CMPLX((i-500)/1001.0,(j-500)/1001.0)
      EXIT
    ENDDO
    !
    !     End of CBEG.
    !
  END FUNCTION CBEG
  !** CBLAT2
  SUBROUTINE CBLAT2(Nout,Kprint,Ipass)
    IMPLICIT NONE
    !>
    !***
    !  Driver for testing Level 2 BLAS complex subroutines.
    !***
    ! **Library:**   SLATEC (BLAS)
    !***
    ! **Category:**  A4
    !***
    ! **Type:**      COMPLEX (SBLAT2-S, DBLAT2-D, CBLAT2-C)
    !***
    ! **Keywords:**  BLAS, QUICK CHECK DRIVER
    !***
    ! **Author:**  Du Croz, J. J., (NAG)
    !           Hanson, R. J., (SNLA)
    !***
    ! **Description:**
    !
    !  Test program for the COMPLEX              Level 2 Blas.
    !
    !***
    ! **References:**  Dongarra, J. J., Du Croz, J. J., Hammarling, S. and
    !                 Hanson, R. J.  An  extended  set of Fortran Basic
    !                 Linear Algebra Subprograms. ACM TOMS, Vol. 14, No. 1,
    !                 pp. 1-17, March 1988.
    !***
    ! **Routines called:**  CCHK12, CCHK22, CCHK32, CCHK42, CCHK52, CCHK62,
    !                    CCHKE2, CMVCH, LCE, R1MACH, XERCLR

    !* REVISION HISTORY  (YYMMDD)
    !   870810  DATE WRITTEN
    !   910619  Modified to meet SLATEC code and prologue standards.  (BKS)
    !   930315  Removed unused variables.  (WRB)
    !   930618  Code modified to improve PASS/FAIL reporting.  (BKS, WRB)

    !     .. Parameters ..
    INTEGER, PARAMETER :: NSUBS = 17
    COMPLEX, PARAMETER :: ZERO = (0.0,0.0), ONE = (1.0,0.0)
    REAL, PARAMETER :: RZERO = 0.0
    INTEGER, PARAMETER :: NMAX = 65, INCMAX = 2
    !     .. Scalar Arguments ..
    INTEGER Ipass, Kprint
    !     .. Local Scalars ..
    REAL eps, err, thresh
    INTEGER i, isnum, j, n, Nout
    INTEGER, PARAMETER :: NIDIM = 6, NKB = 4, NINC = 4, NALF = 3, NBET = 3
    LOGICAL same, tsterr, ftl, ftl1, ftl2
    CHARACTER :: trans
    !     .. Local Arrays ..
    COMPLEX a(NMAX,NMAX), aa(NMAX*NMAX), as(NMAX*NMAX), x(NMAX), xs(NMAX*INCMAX), &
      xx(NMAX*INCMAX), y(NMAX), ys(NMAX*INCMAX), yt(NMAX), yy(NMAX*INCMAX), z(2*NMAX)
    REAL g(NMAX)
    LOGICAL ltest(NSUBS)
    !     .. External Functions ..
    REAL, EXTERNAL :: R1MACH
    !     .. Intrinsic Functions ..
    INTRINSIC ABS, MAX, MIN
    !     .. Data statements ..
    CHARACTER(6), PARAMETER :: snames(NSUBS) = [ 'CGEMV ', 'CGBMV ', 'CHEMV ', &
      'CHBMV ', 'CHPMV ', 'CTRMV ', 'CTBMV ', 'CTPMV ', 'CTRSV ', 'CTBSV ', &
      'CTPSV ', 'CGERC ', 'CGERU ', 'CHER  ', 'CHPR  ', 'CHER2 ', 'CHPR2 ' ]
    INTEGER, PARAMETER :: idim(NIDIM) = [ 0, 1, 2, 3, 5, 9 ]
    INTEGER, PARAMETER :: kb(NKB) = [ 0, 1, 2, 4 ]
    INTEGER, PARAMETER :: inc(NINC) = [ 1, 2, -1, -2 ]
    COMPLEX, PARAMETER :: alf(NALF) = [ (0.0,0.0), (1.0,0.0), (0.7,-0.9) ]
    COMPLEX, PARAMETER :: bet(NBET) = [ (0.0,0.0), (1.0,0.0), (1.3,-1.1) ]
    !* FIRST EXECUTABLE STATEMENT  CBLAT2
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
    DO i = 1, NSUBS
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
    DO j = 1, n
      DO i = 1, n
        a(i,j) = MAX(i-j+1,0)
      ENDDO
      x(j) = j
      y(j) = ZERO
    ENDDO
    DO j = 1, n
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
      IF ( Kprint>=2 ) WRITE (Nout,FMT=99008) trans, same, err
    ENDIF
    trans = 'T'
    ftl = .FALSE.
    CALL CMVCH(trans,n,n,ONE,a,NMAX,x,-1,ZERO,y,-1,yt,g,yy,eps,err,ftl,Nout,&
      .TRUE.,Kprint)
    same = LCE(yy,yt,n)
    IF ( .NOT.same.OR.err/=RZERO ) THEN
      Ipass = 0
      IF ( Kprint>=2 ) WRITE (Nout,FMT=99008) trans, same, err
    ENDIF
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
              NALF,alf,NINC,inc,NMAX,INCMAX,a,aa,as,x,xx,xs,y,yy,ys,yt,g,z)
          CASE (14,15)
            !           Test CHER, 14, and CHPR, 15.
            CALL CCHK52(snames(isnum),eps,thresh,Nout,Kprint,ftl2,NIDIM,idim,&
              NALF,alf,NINC,inc,NMAX,INCMAX,a,aa,as,x,xx,xs,y,yy,ys,yt,g,z)
          CASE (16,17)
            !           Test CHER2, 16, and CHPR2, 17.
            CALL CCHK62(snames(isnum),eps,thresh,Nout,Kprint,ftl2,NIDIM,idim,&
              NALF,alf,NINC,inc,NMAX,INCMAX,a,aa,as,x,xx,xs,y,yy,ys,yt,g,z)
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
  !** CBLAT3
  SUBROUTINE CBLAT3(Nout,Kprint,Ipass)
    IMPLICIT NONE
    !>
    !***
    !  Driver for testing Level 3 BLAS complex subroutines.
    !***
    ! **Library:**   SLATEC (BLAS)
    !***
    ! **Category:**  A4
    !***
    ! **Type:**      COMPLEX (SBLAT3-S, DBLAT3-D, CBLAT3-C)
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
    !  Test program for the COMPLEX              Level 3 Blas.
    !
    !***
    ! **References:**  Dongarra, J., Du Croz, J., Duff, I., and Hammarling, S.
    !                 A set of level 3 basic linear algebra subprograms.
    !                 ACM TOMS, Vol. 16, No. 1, pp. 1-17, March 1990.
    !***
    ! **Routines called:**  CCHK13, CCHK23, CCHK33, CCHK43, CCHK53, CCHKE3,
    !                    CMMCH, LCE, R1MACH, XERCLR

    !* REVISION HISTORY  (YYMMDD)
    !   890208  DATE WRITTEN
    !   910619  Modified to meet SLATEC code and prologue standards.  (BKS)
    !   930315  Removed unused variables.  (WRB)
    !   930618  Code modified to improve PASS/FAIL reporting.  (BKS, WRB)

    !     .. Parameters ..
    INTEGER, PARAMETER :: NSUBS = 9
    COMPLEX, PARAMETER :: ZERO = (0.0,0.0), ONE = (1.0,0.0)
    REAL, PARAMETER :: RZERO = 0.0
    INTEGER, PARAMETER :: NMAX = 65
    !     .. Scalar Arguments ..
    INTEGER Ipass, Kprint
    !     .. Local Scalars ..
    REAL eps, err, thresh
    INTEGER i, isnum, j, n, Nout
    INTEGER, PARAMETER :: NIDIM = 6, NALF = 3, NBET = 3
    LOGICAL same, tsterr, ftl, ftl1, ftl2
    CHARACTER :: transa, transb
    !     .. Local Arrays ..
    COMPLEX aa(NMAX*NMAX), ab(NMAX,2*NMAX), as(NMAX*NMAX), bb(NMAX*NMAX), &
      bs(NMAX*NMAX), c(NMAX,NMAX), cc(NMAX*NMAX), cs(NMAX*NMAX), ct(NMAX), w(2*NMAX)
    REAL g(NMAX)
    LOGICAL ltest(NSUBS)
    !     .. External Functions ..
    REAL, EXTERNAL :: R1MACH
    !     .. Intrinsic Functions ..
    INTRINSIC ABS, MAX, MIN
    !     .. Data statements ..
    CHARACTER(6), PARAMETER :: snames(NSUBS) = [ 'CGEMM ', 'CHEMM ', 'CSYMM ', &
      'CTRMM ', 'CTRSM ', 'CHERK ', 'CSYRK ', 'CHER2K', 'CSYR2K' ]
    INTEGER, PARAMETER :: idim(NIDIM) = [ 0, 1, 2, 3, 5, 9 ]
    COMPLEX, PARAMETER :: alf(NALF) = [ (0.0,0.0), (1.0,0.0), (0.7,-0.9) ]
    COMPLEX, PARAMETER :: bet(NBET) = [ (0.0,0.0), (1.0,0.0), (1.3,-1.1) ]
    !* FIRST EXECUTABLE STATEMENT  CBLAT3
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
      WRITE (Nout,FMT=99004) (alf(i),i=1,NALF)
      WRITE (Nout,FMT=99005) (bet(i),i=1,NBET)
      IF ( .NOT.tsterr ) WRITE (Nout,FMT=99008)
      WRITE (Nout,FMT=99001) thresh
    ENDIF
    !
    !     Set names of subroutines and flags which indicate
    !     whether they are to be tested.
    !
    DO i = 1, NSUBS
      ltest(i) = .TRUE.
    ENDDO
    !
    !     Set EPS (the machine precision).
    !
    eps = R1MACH(4)
    !
    !     Check the reliability of CMMCH using exact data.
    !
    n = MIN(32,NMAX)
    DO j = 1, n
      DO i = 1, n
        ab(i,j) = MAX(i-j+1,0)
      ENDDO
      ab(j,NMAX+1) = j
      ab(1,NMAX+j) = j
      c(j,1) = ZERO
    ENDDO
    DO j = 1, n
      cc(j) = j*((j+1)*j)/2 - ((j+1)*j*(j-1))/3
    ENDDO
    !     CC holds the exact result. On exit from CMMCH CT holds
    !     the result computed by CMMCH.
    transa = 'N'
    transb = 'N'
    ftl = .FALSE.
    CALL CMMCH(transa,transb,n,1,n,ONE,ab,NMAX,ab(1,NMAX+1),NMAX,ZERO,c,NMAX,&
      ct,g,cc,NMAX,eps,err,ftl,Nout,.TRUE.,Kprint)
    same = LCE(cc,ct,n)
    IF ( .NOT.same.OR.err/=RZERO ) THEN
      Ipass = 0
      IF ( Kprint>=2 ) WRITE (Nout,FMT=99006) transa, transb, same, err
    ENDIF
    transb = 'C'
    ftl = .FALSE.
    CALL CMMCH(transa,transb,n,1,n,ONE,ab,NMAX,ab(1,NMAX+1),NMAX,ZERO,c,NMAX,&
      ct,g,cc,NMAX,eps,err,ftl,Nout,.TRUE.,Kprint)
    same = LCE(cc,ct,n)
    IF ( .NOT.same.OR.err/=RZERO ) THEN
      Ipass = 0
      IF ( Kprint>=2 ) WRITE (Nout,FMT=99006) transa, transb, same, err
    ENDIF
    DO j = 1, n
      ab(j,NMAX+1) = n - j + 1
      ab(1,NMAX+j) = n - j + 1
    ENDDO
    DO j = 1, n
      cc(n-j+1) = j*((j+1)*j)/2 - ((j+1)*j*(j-1))/3
    ENDDO
    transa = 'C'
    transb = 'N'
    ftl = .FALSE.
    CALL CMMCH(transa,transb,n,1,n,ONE,ab,NMAX,ab(1,NMAX+1),NMAX,ZERO,c,NMAX,&
      ct,g,cc,NMAX,eps,err,ftl,Nout,.TRUE.,Kprint)
    same = LCE(cc,ct,n)
    IF ( .NOT.same.OR.err/=RZERO ) THEN
      Ipass = 0
      IF ( Kprint>=2 ) WRITE (Nout,FMT=99006) transa, transb, same, err
    ENDIF
    transb = 'C'
    ftl = .FALSE.
    CALL CMMCH(transa,transb,n,1,n,ONE,ab,NMAX,ab(1,NMAX+1),NMAX,ZERO,c,NMAX,&
      ct,g,cc,NMAX,eps,err,ftl,Nout,.TRUE.,Kprint)
    same = LCE(cc,ct,n)
    IF ( .NOT.same.OR.err/=RZERO ) THEN
      Ipass = 0
      IF ( Kprint>=2 ) WRITE (Nout,FMT=99006) transa, transb, same, err
    ENDIF
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
        IF ( tsterr ) CALL CCHKE3(isnum,snames(isnum),Nout,Kprint,ftl1)
        !           Test computations.
        ftl2 = .FALSE.
        CALL XERCLR
        SELECT CASE (isnum)
          CASE (2,3)
            !           Test CHEMM, 02, CSYMM, 03.
            CALL CCHK23(snames(isnum),eps,thresh,Nout,Kprint,ftl2,NIDIM,idim,&
              NALF,alf,NBET,bet,NMAX,ab,aa,as,ab(1,NMAX+1),bb,bs,c,cc,cs,ct,g)
          CASE (4,5)
            !           Test CTRMM, 04, CTRSM, 05.
            CALL CCHK33(snames(isnum),eps,thresh,Nout,Kprint,ftl2,NIDIM,idim,&
              NALF,alf,NMAX,ab,aa,as,ab(1,NMAX+1),bb,bs,ct,g,c)
          CASE (6,7)
            !           Test CHERK, 06, CSYRK, 07.
            CALL CCHK43(snames(isnum),eps,thresh,Nout,Kprint,ftl2,NIDIM,idim,&
              NALF,alf,NBET,bet,NMAX,ab,aa,as,ab(1,NMAX+1),bb,bs,c,cc,cs,ct,g)
          CASE (8,9)
            !           Test CHER2K, 08, CSYR2K, 09.
            CALL CCHK53(snames(isnum),eps,thresh,Nout,Kprint,ftl2,NIDIM,idim,&
              NALF,alf,NBET,bet,NMAX,ab,aa,as,bb,bs,c,cc,cs,ct,g,w)
          CASE DEFAULT
            !           Test CGEMM, 01.
            CALL CCHK13(snames(isnum),eps,thresh,Nout,Kprint,ftl2,NIDIM,idim,&
              NALF,alf,NBET,bet,NMAX,ab,aa,as,ab(1,NMAX+1),bb,bs,c,cc,cs,ct,g)
        END SELECT
        !
        IF ( ftl1.OR.ftl2 ) Ipass = 0
      ENDIF
    ENDDO
    RETURN
    !
    99001 FORMAT (' ROUTINES PASS COMPUTATIONAL TESTS IF TEST RATIO IS LES',&
      'S THAN',F8.2)
    99002 FORMAT (' TESTS OF THE COMPLEX          LEVEL 3 BLAS',//' THE F',&
      'OLLOWING PARAMETER VALUES WILL BE USED:')
    99003 FORMAT ('   FOR N              ',9I6)
    99004 FORMAT ('   FOR ALPHA          ',7('(',F4.1,',',F4.1,')  ',:))
    99005 FORMAT ('   FOR BETA           ',7('(',F4.1,',',F4.1,')  ',:))
    99006 FORMAT (' ERROR IN CMMCH -  IN-LINE DOT PRODUCTS ARE BEING EVALU',&
      'ATED WRONGLY.',/' CMMCH WAS CALLED WITH TRANSA = ',A1,&
      ' AND TRANSB = ',A1,/' AND RETURNED SAME = ',L1,' AND ','ERR = ',&
      F12.3,'.',/' THIS MAY BE DUE TO FAULTS IN THE ',&
      'ARITHMETIC OR THE COMPILER.')
    99007 FORMAT (1X,A6,' WAS NOT TESTED')
    99008 FORMAT (' ERROR-EXITS WILL NOT BE TESTED')
    !
    !     End of CBLAT3.
    !
  END SUBROUTINE CBLAT3
  !** CMAKE2
  SUBROUTINE CMAKE2(Type,Uplo,Diag,M,N,A,Nmax,Aa,Lda,Kl,Ku,Reset,Transl)
    IMPLICIT NONE
    !>
    !***
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
    ! **Routines called:**  CBEG

    !* REVISION HISTORY  (YYMMDD)
    !   870810  DATE WRITTEN
    !   910620  Modified to meet SLATEC code and prologue standards.  (BKS)

    !     .. Parameters ..
    COMPLEX, PARAMETER :: ZERO = (0.0,0.0), ONE = (1.0,0.0)
    COMPLEX, PARAMETER :: ROGUE=(-1.0E10,1.0E10)
    REAL, PARAMETER :: RZERO = 0.0
    REAL, PARAMETER :: RROGUE = -1.0E10
    !     .. Scalar Arguments ..
    COMPLEX Transl
    INTEGER Kl, Ku, Lda, M, N, Nmax
    LOGICAL Reset
    CHARACTER :: Diag, Uplo
    CHARACTER(2) :: Type
    !     .. Array Arguments ..
    COMPLEX A(Nmax,*), Aa(*)
    !     .. Local Scalars ..
    INTEGER i, i1, i2, i3, ibeg, iend, ioff, j, jj, kk
    LOGICAL gen, lower, sym, tri, unit, upper
    !     .. Intrinsic Functions ..
    INTRINSIC CMPLX, CONJG, MAX, MIN, REAL
    !* FIRST EXECUTABLE STATEMENT  CMAKE2
    gen = Type(1:1)=='G'
    sym = Type(1:1)=='H'
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
            A(i,j) = CBEG(Reset) + Transl
          ELSE
            A(i,j) = ZERO
          ENDIF
          IF ( i/=j ) THEN
            IF ( sym ) THEN
              A(j,i) = CONJG(A(i,j))
            ELSEIF ( tri ) THEN
              A(j,i) = ZERO
            ENDIF
          ENDIF
        ENDIF
      ENDDO
      IF ( sym ) A(j,j) = CMPLX(REAL(A(j,j)),RZERO)
      IF ( tri ) A(j,j) = A(j,j) + ONE
      IF ( unit ) A(j,j) = ONE
    ENDDO
    !
    !     Store elements in array AS in data structure required by routine.
    !
    IF ( Type=='GE' ) THEN
      DO j = 1, N
        DO i = 1, M
          Aa(i+(j-1)*Lda) = A(i,j)
        ENDDO
        DO i = M + 1, Lda
          Aa(i+(j-1)*Lda) = ROGUE
        ENDDO
      ENDDO
    ELSEIF ( Type=='GB' ) THEN
      DO j = 1, N
        DO i1 = 1, Ku + 1 - j
          Aa(i1+(j-1)*Lda) = ROGUE
        ENDDO
        DO i2 = i1, MIN(Kl+Ku+1,Ku+1+M-j)
          Aa(i2+(j-1)*Lda) = A(i2+j-Ku-1,j)
        ENDDO
        DO i3 = i2, Lda
          Aa(i3+(j-1)*Lda) = ROGUE
        ENDDO
      ENDDO
    ELSEIF ( Type=='HE'.OR.Type=='TR' ) THEN
      DO j = 1, N
        IF ( upper ) THEN
          ibeg = 1
          IF ( unit ) THEN
            iend = j - 1
          ELSE
            iend = j
          ENDIF
        ELSE
          IF ( unit ) THEN
            ibeg = j + 1
          ELSE
            ibeg = j
          ENDIF
          iend = N
        ENDIF
        DO i = 1, ibeg - 1
          Aa(i+(j-1)*Lda) = ROGUE
        ENDDO
        DO i = ibeg, iend
          Aa(i+(j-1)*Lda) = A(i,j)
        ENDDO
        DO i = iend + 1, Lda
          Aa(i+(j-1)*Lda) = ROGUE
        ENDDO
        IF ( sym ) THEN
          jj = j + (j-1)*Lda
          Aa(jj) = CMPLX(REAL(Aa(jj)),RROGUE)
        ENDIF
      ENDDO
    ELSEIF ( Type=='HB'.OR.Type=='TB' ) THEN
      DO j = 1, N
        IF ( upper ) THEN
          kk = Kl + 1
          ibeg = MAX(1,Kl+2-j)
          IF ( unit ) THEN
            iend = Kl
          ELSE
            iend = Kl + 1
          ENDIF
        ELSE
          kk = 1
          IF ( unit ) THEN
            ibeg = 2
          ELSE
            ibeg = 1
          ENDIF
          iend = MIN(Kl+1,1+M-j)
        ENDIF
        DO i = 1, ibeg - 1
          Aa(i+(j-1)*Lda) = ROGUE
        ENDDO
        DO i = ibeg, iend
          Aa(i+(j-1)*Lda) = A(i+j-kk,j)
        ENDDO
        DO i = iend + 1, Lda
          Aa(i+(j-1)*Lda) = ROGUE
        ENDDO
        IF ( sym ) THEN
          jj = kk + (j-1)*Lda
          Aa(jj) = CMPLX(REAL(Aa(jj)),RROGUE)
        ENDIF
      ENDDO
    ELSEIF ( Type=='HP'.OR.Type=='TP' ) THEN
      ioff = 0
      DO j = 1, N
        IF ( upper ) THEN
          ibeg = 1
          iend = j
        ELSE
          ibeg = j
          iend = N
        ENDIF
        DO i = ibeg, iend
          ioff = ioff + 1
          Aa(ioff) = A(i,j)
          IF ( i==j ) THEN
            IF ( unit ) Aa(ioff) = ROGUE
            IF ( sym ) Aa(ioff) = CMPLX(REAL(Aa(ioff)),RROGUE)
          ENDIF
        ENDDO
      ENDDO
    ENDIF
    !
    !     End of CMAKE2.
    !
  END SUBROUTINE CMAKE2
  !** CMAKE3
  SUBROUTINE CMAKE3(Type,Uplo,Diag,M,N,A,Nmax,Aa,Lda,Reset,Transl)
    IMPLICIT NONE
    !>
    !***
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
    !  TYPE is 'GE', 'HE', 'SY', OR 'TR'.
    !
    !  Auxiliary routine for test program for Level 3 Blas.
    !***
    ! **References:**  (NONE)
    !***
    ! **Routines called:**  CBEG

    !* REVISION HISTORY  (YYMMDD)
    !   890208  DATE WRITTEN
    !   910620  Modified to meet SLATEC code and prologue standards.  (BKS)

    !     .. Parameters ..
    COMPLEX, PARAMETER :: ZERO = (0.0,0.0), ONE = (1.0,0.0)
    COMPLEX, PARAMETER :: ROGUE = (-1.0E10,1.0E10)
    REAL, PARAMETER :: RZERO = 0.0
    REAL, PARAMETER :: RROGUE = -1.0E10
    !     .. Scalar Arguments ..
    COMPLEX Transl
    INTEGER Lda, M, N, Nmax
    LOGICAL Reset
    CHARACTER :: Diag, Uplo
    CHARACTER(2) :: Type
    !     .. Array Arguments ..
    COMPLEX A(Nmax,*), Aa(*)
    !     .. Local Scalars ..
    INTEGER i, ibeg, iend, j, jj
    LOGICAL gen, lower, sym, tri, unit, upper, her
    !     .. Intrinsic Functions ..
    INTRINSIC CMPLX, CONJG, REAL
    !* FIRST EXECUTABLE STATEMENT  CMAKE3
    gen = Type=='GE'
    her = Type=='HE'
    sym = Type=='SY'
    tri = Type=='TR'
    upper = (her.OR.sym.OR.tri) .AND. Uplo=='U'
    lower = (her.OR.sym.OR.tri) .AND. Uplo=='L'
    unit = tri .AND. Diag=='U'
    !
    !     Generate data in array A.
    !
    DO j = 1, N
      DO i = 1, M
        IF ( gen.OR.(upper.AND.i<=j).OR.(lower.AND.i>=j) ) THEN
          A(i,j) = CBEG(Reset) + Transl
          IF ( i/=j ) THEN
            !                 Set some elements to zero
            IF ( N>3.AND.j==N/2 ) A(i,j) = ZERO
            IF ( her ) THEN
              A(j,i) = CONJG(A(i,j))
            ELSEIF ( sym ) THEN
              A(j,i) = A(i,j)
            ELSEIF ( tri ) THEN
              A(j,i) = ZERO
            ENDIF
          ENDIF
        ENDIF
      ENDDO
      IF ( her ) A(j,j) = CMPLX(REAL(A(j,j)),RZERO)
      IF ( tri ) A(j,j) = A(j,j) + ONE
      IF ( unit ) A(j,j) = ONE
    ENDDO
    !
    !     Store elements in array AS in data structure required by routine.
    !
    IF ( Type=='GE' ) THEN
      DO j = 1, N
        DO i = 1, M
          Aa(i+(j-1)*Lda) = A(i,j)
        ENDDO
        DO i = M + 1, Lda
          Aa(i+(j-1)*Lda) = ROGUE
        ENDDO
      ENDDO
    ELSEIF ( Type=='HE'.OR.Type=='SY'.OR.Type=='TR' ) THEN
      DO j = 1, N
        IF ( upper ) THEN
          ibeg = 1
          IF ( unit ) THEN
            iend = j - 1
          ELSE
            iend = j
          ENDIF
        ELSE
          IF ( unit ) THEN
            ibeg = j + 1
          ELSE
            ibeg = j
          ENDIF
          iend = N
        ENDIF
        DO i = 1, ibeg - 1
          Aa(i+(j-1)*Lda) = ROGUE
        ENDDO
        DO i = ibeg, iend
          Aa(i+(j-1)*Lda) = A(i,j)
        ENDDO
        DO i = iend + 1, Lda
          Aa(i+(j-1)*Lda) = ROGUE
        ENDDO
        IF ( her ) THEN
          jj = j + (j-1)*Lda
          Aa(jj) = CMPLX(REAL(Aa(jj)),RROGUE)
        ENDIF
      ENDDO
    ENDIF
    !
    !     End of CMAKE3.
    !
  END SUBROUTINE CMAKE3
  !** CMMCH
  SUBROUTINE CMMCH(Transa,Transb,M,N,Kk,Alpha,A,Lda,B,Ldb,Beta,C,Ldc,Ct,G,&
      Cc,Ldcc,Eps,Err,Ftl,Nout,Mv,Kprint)
    IMPLICIT NONE
    !>
    !***
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
    COMPLEX, PARAMETER :: ZERO = (0.0,0.0)
    REAL, PARAMETER :: RZERO = 0.0, RONE = 1.0
    !     .. Scalar Arguments ..
    LOGICAL Ftl
    COMPLEX Alpha, Beta
    REAL Eps, Err
    INTEGER Kk, Kprint, Lda, Ldb, Ldc, Ldcc, M, N, Nout
    LOGICAL Mv
    CHARACTER :: Transa, Transb
    !     .. Array Arguments ..
    COMPLEX A(Lda,*), B(Ldb,*), C(Ldc,*), Cc(Ldcc,*), Ct(*)
    REAL G(*)
    !     .. Local Scalars ..
    REAL erri
    INTEGER i, j, k
    LOGICAL ctrana, ctranb, trana, tranb
    !     .. Intrinsic Functions ..
    INTRINSIC ABS, AIMAG, CONJG, MAX, REAL, SQRT
    REAL CABS1
    !* FIRST EXECUTABLE STATEMENT  CMMCH
    trana = Transa=='T' .OR. Transa=='C'
    tranb = Transb=='T' .OR. Transb=='C'
    ctrana = Transa=='C'
    ctranb = Transb=='C'
    !
    !     Compute expected result, one column at a time, in CT using data
    !     in A, B and C.
    !     Compute gauges in G.
    !
    DO j = 1, N
      !
      DO i = 1, M
        Ct(i) = ZERO
        G(i) = RZERO
      ENDDO
      IF ( .NOT.trana.AND..NOT.tranb ) THEN
        DO k = 1, Kk
          DO i = 1, M
            Ct(i) = Ct(i) + A(i,k)*B(k,j)
            G(i) = G(i) + CABS1(A(i,k))*CABS1(B(k,j))
          ENDDO
        ENDDO
      ELSEIF ( trana.AND..NOT.tranb ) THEN
        IF ( ctrana ) THEN
          DO k = 1, Kk
            DO i = 1, M
              Ct(i) = Ct(i) + CONJG(A(k,i))*B(k,j)
              G(i) = G(i) + CABS1(A(k,i))*CABS1(B(k,j))
            ENDDO
          ENDDO
        ELSE
          DO k = 1, Kk
            DO i = 1, M
              Ct(i) = Ct(i) + A(k,i)*B(k,j)
              G(i) = G(i) + CABS1(A(k,i))*CABS1(B(k,j))
            ENDDO
          ENDDO
        ENDIF
      ELSEIF ( .NOT.trana.AND.tranb ) THEN
        IF ( ctranb ) THEN
          DO k = 1, Kk
            DO i = 1, M
              Ct(i) = Ct(i) + A(i,k)*CONJG(B(j,k))
              G(i) = G(i) + CABS1(A(i,k))*CABS1(B(j,k))
            ENDDO
          ENDDO
        ELSE
          DO k = 1, Kk
            DO i = 1, M
              Ct(i) = Ct(i) + A(i,k)*B(j,k)
              G(i) = G(i) + CABS1(A(i,k))*CABS1(B(j,k))
            ENDDO
          ENDDO
        ENDIF
      ELSEIF ( trana.AND.tranb ) THEN
        IF ( ctrana ) THEN
          IF ( ctranb ) THEN
            DO k = 1, Kk
              DO i = 1, M
                Ct(i) = Ct(i) + CONJG(A(k,i))*CONJG(B(j,k))
                G(i) = G(i) + CABS1(A(k,i))*CABS1(B(j,k))
              ENDDO
            ENDDO
          ELSE
            DO k = 1, Kk
              DO i = 1, M
                Ct(i) = Ct(i) + CONJG(A(k,i))*B(j,k)
                G(i) = G(i) + CABS1(A(k,i))*CABS1(B(j,k))
              ENDDO
            ENDDO
          ENDIF
        ELSEIF ( ctranb ) THEN
          DO k = 1, Kk
            DO i = 1, M
              Ct(i) = Ct(i) + A(k,i)*CONJG(B(j,k))
              G(i) = G(i) + CABS1(A(k,i))*CABS1(B(j,k))
            ENDDO
          ENDDO
        ELSE
          DO k = 1, Kk
            DO i = 1, M
              Ct(i) = Ct(i) + A(k,i)*B(j,k)
              G(i) = G(i) + CABS1(A(k,i))*CABS1(B(j,k))
            ENDDO
          ENDDO
        ENDIF
      ENDIF
      DO i = 1, M
        Ct(i) = Alpha*Ct(i) + Beta*C(i,j)
        G(i) = CABS1(Alpha)*G(i) + CABS1(Beta)*CABS1(C(i,j))
      ENDDO
      !
      !        Compute the error ratio for this result.
      !
      Err = ZERO
      DO i = 1, M
        erri = CABS1(Ct(i)-Cc(i,j))/Eps
        IF ( G(i)/=RZERO ) erri = erri/G(i)
        Err = MAX(Err,erri)
        IF ( Err*SQRT(Eps)>=RONE ) THEN
          Ftl = .TRUE.
          IF ( Kprint>=2 ) THEN
            WRITE (Nout,FMT=99001)
            DO k = 1, M
              IF ( Mv ) THEN
                WRITE (Nout,FMT=99002) k, Ct(k), Cc(k,j)
              ELSE
                WRITE (Nout,FMT=99002) k, Cc(k,j), Ct(k)
              ENDIF
            ENDDO
          ENDIF
        ENDIF
      ENDDO
    ENDDO
    RETURN
    !
    99001 FORMAT (' ******* FATAL ERROR - COMPUTED RESULT IS LESS THAN HAL',&
      'F ACCURATE *******',/'                       EXPECTED RE',&
      'SULT                    COMPUTED RESULT')
    99002 FORMAT (1X,I7,2('  (',G15.6,',',G15.6,')'))
    !
    !     End of CMMCH.
    !
  END SUBROUTINE CMMCH
  !** CMVCH
  SUBROUTINE CMVCH(Trans,M,N,Alpha,A,Nmax,X,Incx,Beta,Y,Incy,Yt,G,Yy,Eps,&
      Err,Ftl,Nout,Mv,Kprint)
    IMPLICIT NONE
    !>
    !***
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
    COMPLEX, PARAMETER :: ZERO = (0.0,0.0)
    REAL, PARAMETER :: RZERO = 0.0, RONE = 1.0
    !     .. Scalar Arguments ..
    COMPLEX Alpha, Beta
    REAL Eps, Err
    INTEGER Incx, Incy, Kprint, M, N, Nmax, Nout
    LOGICAL Mv, Ftl
    CHARACTER :: Trans
    !     .. Array Arguments ..
    COMPLEX A(Nmax,*), X(*), Y(*), Yt(*), Yy(*)
    REAL G(*)
    !     .. Local Scalars ..
    REAL erri
    INTEGER i, incxl, incyl, iy, j, jx, k, kx, ky, ml, nl
    LOGICAL ctran, tran
    !     .. Intrinsic Functions ..
    INTRINSIC ABS, AIMAG, CONJG, MAX, REAL, SQRT
    REAL CABS1
    !* FIRST EXECUTABLE STATEMENT  CMVCH
    tran = Trans=='T'
    ctran = Trans=='C'
    IF ( tran.OR.ctran ) THEN
      ml = N
      nl = M
    ELSE
      ml = M
      nl = N
    ENDIF
    IF ( Incx<0 ) THEN
      kx = nl
      incxl = -1
    ELSE
      kx = 1
      incxl = 1
    ENDIF
    IF ( Incy<0 ) THEN
      ky = ml
      incyl = -1
    ELSE
      ky = 1
      incyl = 1
    ENDIF
    !
    !     Compute expected result in YT using data in A, X and Y.
    !     Compute gauges in G.
    !
    iy = ky
    DO i = 1, ml
      Yt(iy) = ZERO
      G(iy) = RZERO
      jx = kx
      IF ( tran ) THEN
        DO j = 1, nl
          Yt(iy) = Yt(iy) + A(j,i)*X(jx)
          G(iy) = G(iy) + CABS1(A(j,i))*CABS1(X(jx))
          jx = jx + incxl
        ENDDO
      ELSEIF ( ctran ) THEN
        DO j = 1, nl
          Yt(iy) = Yt(iy) + CONJG(A(j,i))*X(jx)
          G(iy) = G(iy) + CABS1(A(j,i))*CABS1(X(jx))
          jx = jx + incxl
        ENDDO
      ELSE
        DO j = 1, nl
          Yt(iy) = Yt(iy) + A(i,j)*X(jx)
          G(iy) = G(iy) + CABS1(A(i,j))*CABS1(X(jx))
          jx = jx + incxl
        ENDDO
      ENDIF
      Yt(iy) = Alpha*Yt(iy) + Beta*Y(iy)
      G(iy) = CABS1(Alpha)*G(iy) + CABS1(Beta)*CABS1(Y(iy))
      iy = iy + incyl
    ENDDO
    !
    !     Compute the error ratio for this result.
    !
    Err = ZERO
    DO i = 1, ml
      erri = ABS(Yt(i)-Yy(1+(i-1)*ABS(Incy)))/Eps
      IF ( G(i)/=RZERO ) erri = erri/G(i)
      Err = MAX(Err,erri)
      IF ( Err*SQRT(Eps)>=RONE ) THEN
        Ftl = .TRUE.
        IF ( Kprint>=2 ) THEN
          WRITE (Nout,FMT=99001)
          DO k = 1, ml
            IF ( Mv ) THEN
              WRITE (Nout,FMT=99002) k, Yt(k), Yy(1+(k-1)*ABS(Incy))
            ELSE
              WRITE (Nout,FMT=99002) i, Yy(1+(k-1)*ABS(Incy)), Yt(k)
            ENDIF
          ENDDO
        ENDIF
      ENDIF
      !
    ENDDO
    RETURN
    !
    99001 FORMAT (' ******* FATAL ERROR - COMPUTED RESULT IS LESS THAN HAL',&
      'F ACCURATE *******',/'                       EXPECTED RE',&
      'SULT                    COMPUTED RESULT')
    99002 FORMAT (1X,I7,2('  (',G15.6,',',G15.6,')'))
    !
    !     End of CMVCH.
    !
  END SUBROUTINE CMVCH
  !** LCE
  LOGICAL FUNCTION LCE(Ri,Rj,Lr)
    IMPLICIT NONE
    !>
    !***
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
    COMPLEX Ri(*), Rj(*)
    !     .. Local Scalars ..
    INTEGER i
    !* FIRST EXECUTABLE STATEMENT  LCE
    LCE = .TRUE.
    DO i = 1, Lr
      IF ( Ri(i)/=Rj(i) ) THEN
        LCE = .FALSE.
        EXIT
      ENDIF
    ENDDO
    !
    !     End of LCE.
    !
  END FUNCTION LCE
  !** LCERES
  LOGICAL FUNCTION LCERES(Type,Uplo,M,N,Aa,As,Lda)
    IMPLICIT NONE
    !>
    !***
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
    !  TYPE is 'GE', 'HE' or 'HP'.
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
    COMPLEX Aa(Lda,*), As(Lda,*)
    !     .. Local Scalars ..
    INTEGER i, ibeg, iend, j
    LOGICAL upper
    !* FIRST EXECUTABLE STATEMENT  LCERES
    upper = Uplo=='U'
    IF ( Type=='GE' ) THEN
      DO j = 1, N
        DO i = M + 1, Lda
          IF ( Aa(i,j)/=As(i,j) ) GOTO 100
        ENDDO
      ENDDO
    ELSEIF ( Type=='HE' ) THEN
      DO j = 1, N
        IF ( upper ) THEN
          ibeg = 1
          iend = j
        ELSE
          ibeg = j
          iend = N
        ENDIF
        DO i = 1, ibeg - 1
          IF ( Aa(i,j)/=As(i,j) ) GOTO 100
        ENDDO
        DO i = iend + 1, Lda
          IF ( Aa(i,j)/=As(i,j) ) GOTO 100
        ENDDO
      ENDDO
    ENDIF
    !
    LCERES = .TRUE.
    RETURN
    100  LCERES = .FALSE.
    !
    !     End of LCERES.
    !
    RETURN
  END FUNCTION LCERES
  !** CCHK12
  SUBROUTINE CCHK12(Sname,Eps,Thresh,Nout,Kprint,Fatal,Nidim,Idim,Nkb,Kb,&
      Nalf,Alf,Nbet,Bet,Ninc,Inc,Nmax,Incmax,A,Aa,As,X,Xx,Xs,Y,Yy,Ys,Yt,G)
    IMPLICIT NONE
    !>
    !***
    !  Quick check for CGEMV and CGBMV.
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
    !  Quick check for CGEMV and CGBMV.
    !
    !  Auxiliary routine for test program for Level 2 Blas.
    !***
    ! **References:**  (NONE)
    !***
    ! **Routines called:**  CGBMV, CGEMV, CMAKE2, CMVCH, LCE, LCERES, NUMXER

    !* REVISION HISTORY  (YYMMDD)
    !   870810  DATE WRITTEN
    !   910619  Modified to meet SLATEC code and prologue standards.  (BKS)

    !     .. Parameters ..
    COMPLEX, PARAMETER :: ZERO = (0.0,0.0), HALF = (0.5,0.0)
    REAL, PARAMETER :: RZERO = 0.0
    !     .. Scalar Arguments ..
    LOGICAL Fatal
    REAL Eps, Thresh
    INTEGER Incmax, Kprint, Nalf, Nbet, Nidim, Ninc, Nkb, Nmax, Nout
    CHARACTER(6) :: Sname
    !     .. Array Arguments ..
    COMPLEX A(Nmax,Nmax), Aa(Nmax*Nmax), Alf(Nalf), As(Nmax*Nmax), &
      Bet(Nbet), X(Nmax), Xs(Nmax*Incmax), Xx(Nmax*Incmax), Y(Nmax), &
      Ys(Nmax*Incmax), Yt(Nmax), Yy(Nmax*Incmax)
    REAL G(Nmax)
    INTEGER Idim(Nidim), Inc(Ninc), Kb(Nkb)
    !     .. Local Scalars ..
    COMPLEX alpha, als, beta, bls, transl
    REAL err, errmax
    INTEGER i, ia, ib, ic, iku, im, in, incx, incxs, incy, incys, &
      ix, iy, kl, kls, ku, kus, laa, lda, ldas, lx, ly, m, &
      ml, ms, n, nargs, nc, nd, nk, nerr, nl, ns
    LOGICAL banded, ftl, full, null, reset, tran
    CHARACTER :: trans, transs
    !     .. Local Arrays ..
    LOGICAL isame(13)
    !     .. External Functions ..
    INTEGER, EXTERNAL :: NUMXER
    !     .. External Subroutines ..
    EXTERNAL :: CGBMV, CGEMV
    !     .. Intrinsic Functions ..
    INTRINSIC ABS, MAX, MIN
    !     .. Data statements ..
    CHARACTER(3), PARAMETER :: ich = 'NTC'
    !* FIRST EXECUTABLE STATEMENT  CCHK12
    full = Sname(3:3)=='E'
    banded = Sname(3:3)=='B'
    !     Define the number of arguments.
    IF ( full ) THEN
      nargs = 11
    ELSEIF ( banded ) THEN
      nargs = 13
    ENDIF
    !
    nc = 0
    reset = .TRUE.
    errmax = RZERO
    !
    DO in = 1, Nidim
      n = Idim(in)
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
        ENDIF
        DO iku = 1, nk
          IF ( banded ) THEN
            ku = Kb(iku)
            kl = MAX(ku-1,0)
          ELSE
            ku = n - 1
            kl = m - 1
          ENDIF
          !              Set LDA to 1 more than minimum value if room.
          IF ( banded ) THEN
            lda = kl + ku + 1
          ELSE
            lda = m
          ENDIF
          IF ( lda<Nmax ) lda = lda + 1
          !              Skip tests if not enough room.
          IF ( lda<=Nmax ) THEN
            laa = lda*n
            null = n<=0 .OR. m<=0
            !
            !              Generate the matrix A.
            !
            transl = ZERO
            CALL CMAKE2(Sname(2:3),' ',' ',m,n,A,Nmax,Aa,lda,kl,ku,reset,transl)
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
              ENDIF
              !
              DO ix = 1, Ninc
                incx = Inc(ix)
                lx = ABS(incx)*nl
                !
                !                    Generate the vector X.
                !
                transl = HALF
                CALL CMAKE2('GE',' ',' ',1,nl,X,1,Xx,ABS(incx),0,nl-1,reset,&
                  transl)
                IF ( nl>1 ) THEN
                  X(nl/2) = ZERO
                  Xx(1+ABS(incx)*(nl/2-1)) = ZERO
                ENDIF
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
                      CALL CMAKE2('GE',' ',' ',1,ml,Y,1,Yy,ABS(incy),0,ml-1,&
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
                      ENDDO
                      ldas = lda
                      DO i = 1, lx
                        Xs(i) = Xx(i)
                      ENDDO
                      incxs = incx
                      bls = beta
                      DO i = 1, ly
                        Ys(i) = Yy(i)
                      ENDDO
                      incys = incy
                      !
                      !                             Call the subroutine.
                      !
                      IF ( full ) THEN
                        CALL CGEMV(trans,m,n,alpha,Aa,lda,Xx,incx,beta,Yy,incy)
                      ELSEIF ( banded ) THEN
                        CALL CGBMV(trans,m,n,kl,ku,alpha,Aa,lda,Xx,incx,beta,&
                          Yy,incy)
                      ENDIF
                      !
                      !                             Check if error-exit was taken incorrectly.
                      !
                      IF ( NUMXER(nerr)/=0 ) THEN
                        IF ( Kprint>=2 ) WRITE (Nout,FMT=99007)
                        Fatal = .TRUE.
                      ENDIF
                      !
                      !                             See what data changed inside subroutines.
                      !
                      isame(1) = trans==transs
                      isame(2) = ms==m
                      isame(3) = ns==n
                      IF ( full ) THEN
                        isame(4) = als==alpha
                        isame(5) = LCE(As,Aa,laa)
                        isame(6) = ldas==lda
                        isame(7) = LCE(Xs,Xx,lx)
                        isame(8) = incxs==incx
                        isame(9) = bls==beta
                        IF ( null ) THEN
                          isame(10) = LCE(Ys,Yy,ly)
                        ELSE
                          isame(10) = LCERES('GE',' ',1,ml,Ys,Yy,ABS(incy))
                        ENDIF
                        isame(11) = incys==incy
                      ELSEIF ( banded ) THEN
                        isame(4) = kls==kl
                        isame(5) = kus==ku
                        isame(6) = als==alpha
                        isame(7) = LCE(As,Aa,laa)
                        isame(8) = ldas==lda
                        isame(9) = LCE(Xs,Xx,lx)
                        isame(10) = incxs==incx
                        isame(11) = bls==beta
                        IF ( null ) THEN
                          isame(12) = LCE(Ys,Yy,ly)
                        ELSE
                          isame(12) = LCERES('GE',' ',1,ml,Ys,Yy,ABS(incy))
                        ENDIF
                        isame(13) = incys==incy
                      ENDIF
                      !
                      !                             If data was incorrectly changed, report
                      !                             and return.
                      !
                      DO i = 1, nargs
                        IF ( .NOT.isame(i) ) THEN
                          Fatal = .TRUE.
                          IF ( Kprint>=2 ) WRITE (Nout,FMT=99002) i
                        ENDIF
                      ENDDO
                      !
                      ftl = .FALSE.
                      IF ( .NOT.null ) THEN
                        !
                        !                                Check the result.
                        !
                        CALL CMVCH(trans,m,n,alpha,A,Nmax,X,incx,beta,Y,incy,&
                          Yt,G,Yy,Eps,err,ftl,Nout,.TRUE.,Kprint)
                        errmax = MAX(errmax,err)
                      ENDIF
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
                          ENDIF
                        ENDIF
                      ENDIF
                      !
                    ENDDO
                    !
                  ENDDO
                  !
                ENDDO
                !
              ENDDO
              !
            ENDDO
          ENDIF
          !
        ENDDO
        !
      ENDDO
      !
    ENDDO
    !
    !     Report result.
    !
    IF ( .NOT.Fatal ) THEN
      IF ( Kprint>=3 ) THEN
        IF ( errmax<Thresh ) THEN
          WRITE (Nout,FMT=99001) Sname, nc
        ELSE
          WRITE (Nout,FMT=99003) Sname, nc, errmax
        ENDIF
      ENDIF
    ENDIF
    RETURN
    !
    99001 FORMAT (' ',A6,' PASSED THE COMPUTATIONAL TESTS (',I6,' CALL','S)')
    99002 FORMAT (' ******* FATAL ERROR - PARAMETER NUMBER ',I2,' WAS CH',&
      'ANGED INCORRECTLY *******')
    99003 FORMAT (' ',A6,' COMPLETED THE COMPUTATIONAL TESTS (',I6,' C','ALLS)',&
      /' ******* BUT WITH MAXIMUM TEST RATIO',F8.2,' - SUSPECT *******')
    99004 FORMAT (' ******* ',A6,' FAILED ON CALL NUMBER:')
    99005 FORMAT (1X,I6,': ',A6,'(''',A1,''',',4(I3,','),'(',F4.1,',',F4.1,'), A,',&
      I3,', X,',I2,',(',F4.1,',',F4.1,'), Y,',I2,') .')
    99006 FORMAT (1X,I6,': ',A6,'(''',A1,''',',2(I3,','),'(',F4.1,',',F4.1,'), A,',&
      I3,', X,',I2,',(',F4.1,',',F4.1,'), Y,',I2,')         .')
    99007 FORMAT (' ******* FATAL ERROR - ERROR-EXIT TAKEN ON VALID CALL *',&
      '******')
    !
    !     End of CCHK12.
    !
  END SUBROUTINE CCHK12
  !** CCHK13
  SUBROUTINE CCHK13(Sname,Eps,Thresh,Nout,Kprint,Fatal,Nidim,Idim,Nalf,Alf,&
      Nbet,Bet,Nmax,A,Aa,As,B,Bb,Bs,C,Cc,Cs,Ct,G)
    IMPLICIT NONE
    !>
    !***
    !  Quick check for CGEMM.
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
    !  Quick check for CGEMM.
    !
    !  Auxiliary routine for test program for Level 3 Blas.
    !***
    ! **References:**  (NONE)
    !***
    ! **Routines called:**  CGEMM, CMAKE3, CMMCH, LCE, LCERES, NUMXER

    !* REVISION HISTORY  (YYMMDD)
    !   890208  DATE WRITTEN
    !   910619  Modified to meet SLATEC code and prologue standards.  (BKS)

    !     .. Parameters ..
    COMPLEX, PARAMETER :: ZERO = (0.0,0.0)
    REAL, PARAMETER :: RZERO = 0.0
    !     .. Scalar Arguments ..
    LOGICAL Fatal
    REAL Eps, Thresh
    INTEGER Kprint, Nalf, Nbet, Nidim, Nmax, Nout
    CHARACTER(6) :: Sname
    !     .. Array Arguments ..
    COMPLEX A(Nmax,Nmax), Aa(Nmax*Nmax), Alf(Nalf), As(Nmax*Nmax), &
      B(Nmax,Nmax), Bb(Nmax*Nmax), Bet(Nbet), Bs(Nmax*Nmax), &
      C(Nmax,Nmax), Cc(Nmax*Nmax), Cs(Nmax*Nmax), Ct(Nmax)
    REAL G(Nmax)
    INTEGER Idim(Nidim)
    !     .. Local Scalars ..
    COMPLEX alpha, als, beta, bls
    REAL err, errmax
    INTEGER i, ia, ib, ica, icb, ik, im, in, k, ks, laa, lbb, &
      lcc, lda, ldas, ldb, ldbs, ldc, ldcs, m, ma, mb, ms, &
      n, na, nargs, nb, nc, nerr, ns
    LOGICAL ftl, null, reset, trana, tranb
    CHARACTER :: tranas, tranbs, transa, transb
    !     .. Local Arrays ..
    LOGICAL isame(13)
    !     .. External Functions ..
    INTEGER, EXTERNAL :: NUMXER
    !     .. External Subroutines ..
    EXTERNAL :: CGEMM
    !     .. Intrinsic Functions ..
    INTRINSIC ABS, MAX, MIN
    !     .. Data statements ..
    CHARACTER(3), PARAMETER :: ich = 'NTC'
    !* FIRST EXECUTABLE STATEMENT  CCHK13
    nargs = 13
    nc = 0
    reset = .TRUE.
    errmax = RZERO
    !
    DO im = 1, Nidim
      m = Idim(im)
      !
      DO in = 1, Nidim
        n = Idim(in)
        !           Set LDC to 1 more than minimum value if room.
        ldc = m
        IF ( ldc<Nmax ) ldc = ldc + 1
        !           Skip tests if not enough room.
        IF ( ldc<=Nmax ) THEN
          lcc = ldc*n
          null = n<=0 .OR. m<=0
          !
          DO ik = 1, Nidim
            k = Idim(ik)
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
              ENDIF
              !                 Set LDA to 1 more than minimum value if room.
              lda = ma
              IF ( lda<Nmax ) lda = lda + 1
              !                 Skip tests if not enough room.
              IF ( lda<=Nmax ) THEN
                laa = lda*na
                !
                !                 Generate the matrix A.
                !
                CALL CMAKE3('GE',' ',' ',ma,na,A,Nmax,Aa,lda,reset,ZERO)
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
                  ENDIF
                  !                    Set LDB to 1 more than minimum value if room.
                  ldb = mb
                  IF ( ldb<Nmax ) ldb = ldb + 1
                  !                    Skip tests if not enough room.
                  IF ( ldb<=Nmax ) THEN
                    lbb = ldb*nb
                    !
                    !                    Generate the matrix B.
                    !
                    CALL CMAKE3('GE',' ',' ',mb,nb,B,Nmax,Bb,ldb,reset,ZERO)
                    !
                    DO ia = 1, Nalf
                      alpha = Alf(ia)
                      !
                      DO ib = 1, Nbet
                        beta = Bet(ib)
                        !
                        !                          Generate the matrix C.
                        !
                        CALL CMAKE3('GE',' ',' ',m,n,C,Nmax,Cc,ldc,reset,ZERO)
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
                        ENDDO
                        ldas = lda
                        DO i = 1, lbb
                          Bs(i) = Bb(i)
                        ENDDO
                        ldbs = ldb
                        bls = beta
                        DO i = 1, lcc
                          Cs(i) = Cc(i)
                        ENDDO
                        ldcs = ldc
                        !
                        !                          Call the subroutine.
                        !
                        CALL CGEMM(transa,transb,m,n,k,alpha,Aa,lda,Bb,ldb,&
                          beta,Cc,ldc)
                        !
                        !                          Check if error-exit was taken incorrectly.
                        !
                        IF ( NUMXER(nerr)/=0 ) THEN
                          IF ( Kprint>=2 ) WRITE (Nout,FMT=99006)
                          Fatal = .TRUE.
                        ENDIF
                        !
                        !                          See what data changed inside subroutines.
                        !
                        isame(1) = transa==tranas
                        isame(2) = transb==tranbs
                        isame(3) = ms==m
                        isame(4) = ns==n
                        isame(5) = ks==k
                        isame(6) = als==alpha
                        isame(7) = LCE(As,Aa,laa)
                        isame(8) = ldas==lda
                        isame(9) = LCE(Bs,Bb,lbb)
                        isame(10) = ldbs==ldb
                        isame(11) = bls==beta
                        IF ( null ) THEN
                          isame(12) = LCE(Cs,Cc,lcc)
                        ELSE
                          isame(12) = LCERES('GE',' ',m,n,Cs,Cc,ldc)
                        ENDIF
                        isame(13) = ldcs==ldc
                        !
                        !                          If data was incorrectly changed, report
                        !
                        DO i = 1, nargs
                          IF ( .NOT.isame(i) ) THEN
                            Fatal = .TRUE.
                            IF ( Kprint>=2 ) WRITE (Nout,FMT=99002) i
                          ENDIF
                        ENDDO
                        !
                        ftl = .FALSE.
                        IF ( .NOT.null ) THEN
                          !
                          !                             Check the result.
                          !
                          CALL CMMCH(transa,transb,m,n,k,alpha,A,Nmax,B,Nmax,&
                            beta,C,Nmax,Ct,G,Cc,ldc,Eps,err,ftl,Nout,.TRUE.,Kprint)
                          errmax = MAX(errmax,err)
                        ENDIF
                        IF ( ftl ) THEN
                          Fatal = .TRUE.
                          IF ( Kprint>=3 ) THEN
                            WRITE (Nout,FMT=99004) Sname
                            WRITE (Nout,FMT=99005) nc, Sname, transa, &
                              transb, m, n, k, alpha, lda, ldb, beta, ldc
                          ENDIF
                        ENDIF
                      ENDDO
                      !
                    ENDDO
                  ENDIF
                  !
                ENDDO
              ENDIF
              !
            ENDDO
            !
          ENDDO
        ENDIF
        !
      ENDDO
      !
    ENDDO
    !
    !     Report result.
    !
    IF ( .NOT.Fatal ) THEN
      IF ( Kprint>=3 ) THEN
        IF ( errmax<Thresh ) THEN
          WRITE (Nout,FMT=99001) Sname, nc
        ELSE
          WRITE (Nout,FMT=99003) Sname, nc, errmax
        ENDIF
      ENDIF
    ENDIF
    RETURN
    !
    99001 FORMAT (' ',A6,' PASSED THE COMPUTATIONAL TESTS (',I6,' CALL','S)')
    99002 FORMAT (' ******* FATAL ERROR - PARAMETER NUMBER ',I2,' WAS CH',&
      'ANGED INCORRECTLY *******')
    99003 FORMAT (' ',A6,' COMPLETED THE COMPUTATIONAL TESTS (',I6,' C','ALLS)',&
      /' ******* BUT WITH MAXIMUM TEST RATIO',F8.2,' - SUSPECT *******')
    99004 FORMAT (' ******* ',A6,' FAILED ON CALL NUMBER:')
    99005 FORMAT (1X,I6,': ',A6,'(''',A1,''',''',A1,''',',3(I3,','),'(',F4.1,',',&
      F4.1,'), A,',I3,', B,',I3,',(',F4.1,',',F4.1,'), C,',I3,').')
    99006 FORMAT (' ******* FATAL ERROR - ERROR-EXIT TAKEN ON VALID CALL *',&
      '******')
    !
    !     End of CCHK13.
    !
  END SUBROUTINE CCHK13
  !** CCHK22
  SUBROUTINE CCHK22(Sname,Eps,Thresh,Nout,Kprint,Fatal,Nidim,Idim,Nkb,Kb,&
      Nalf,Alf,Nbet,Bet,Ninc,Inc,Nmax,Incmax,A,Aa,As,X,Xx,Xs,Y,Yy,Ys,Yt,G)
    IMPLICIT NONE
    !>
    !***
    !  Quick check for CHEMV, CHBMV, CHPMV.
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
    !  Quick check for CHEMV, CHBMV and CHPMV.
    !
    !  Auxiliary routine for test program for Level 2 Blas.
    !***
    ! **References:**  (NONE)
    !***
    ! **Routines called:**  CHBMV, CHEMV, CHPMV, CMAKE2, CMVCH, LCE, LCERES,
    !                    NUMXER

    !* REVISION HISTORY  (YYMMDD)
    !   870810  DATE WRITTEN
    !   910619  Modified to meet SLATEC code and prologue standards.  (BKS)

    !     .. Parameters ..
    COMPLEX, PARAMETER :: ZERO = (0.0,0.0), HALF = (0.5,0.0)
    REAL, PARAMETER :: RZERO = 0.0
    !     .. Scalar Arguments ..
    LOGICAL Fatal
    REAL Eps, Thresh
    INTEGER Incmax, Kprint, Nalf, Nbet, Nidim, Ninc, Nkb, Nmax, Nout
    CHARACTER(6) :: Sname
    !     .. Array Arguments ..
    COMPLEX A(Nmax,Nmax), Aa(Nmax*Nmax), Alf(Nalf), As(Nmax*Nmax), &
      Bet(Nbet), X(Nmax), Xs(Nmax*Incmax), Xx(Nmax*Incmax), Y(Nmax), &
      Ys(Nmax*Incmax), Yt(Nmax), Yy(Nmax*Incmax)
    REAL G(Nmax)
    INTEGER Idim(Nidim), Inc(Ninc), Kb(Nkb)
    !     .. Local Scalars ..
    COMPLEX alpha, als, beta, bls, transl
    REAL err, errmax
    INTEGER i, ia, ib, ic, ik, in, incx, incxs, incy, incys, ix, &
      iy, k, ks, laa, lda, ldas, lx, ly, n, nargs, nc, nerr, nk, ns
    LOGICAL banded, ftl, full, null, packed, reset
    CHARACTER :: uplo, uplos
    !     .. Local Arrays ..
    LOGICAL isame(13)
    !     .. External Functions ..
    INTEGER, EXTERNAL :: NUMXER
    !     .. External Subroutines ..
    EXTERNAL :: CHBMV, CHEMV, CHPMV
    !     .. Intrinsic Functions ..
    INTRINSIC ABS, MAX
    !     .. Data statements ..
    CHARACTER(2), PARAMETER :: ich = 'UL'
    !* FIRST EXECUTABLE STATEMENT  CCHK22
    full = Sname(3:3)=='E'
    banded = Sname(3:3)=='B'
    packed = Sname(3:3)=='P'
    !     Define the number of arguments.
    IF ( full ) THEN
      nargs = 10
    ELSEIF ( banded ) THEN
      nargs = 11
    ELSEIF ( packed ) THEN
      nargs = 9
    ENDIF
    !
    nc = 0
    reset = .TRUE.
    errmax = RZERO
    !
    DO in = 1, Nidim
      n = Idim(in)
      !
      IF ( banded ) THEN
        nk = Nkb
      ELSE
        nk = 1
      ENDIF
      DO ik = 1, nk
        IF ( banded ) THEN
          k = Kb(ik)
        ELSE
          k = n - 1
        ENDIF
        !           Set LDA to 1 more than minimum value if room.
        IF ( banded ) THEN
          lda = k + 1
        ELSE
          lda = n
        ENDIF
        IF ( lda<Nmax ) lda = lda + 1
        !           Skip tests if not enough room.
        IF ( lda<=Nmax ) THEN
          IF ( packed ) THEN
            laa = (n*(n+1))/2
          ELSE
            laa = lda*n
          ENDIF
          null = n<=0
          !
          DO ic = 1, 2
            uplo = ich(ic:ic)
            !
            !              Generate the matrix A.
            !
            transl = ZERO
            CALL CMAKE2(Sname(2:3),uplo,' ',n,n,A,Nmax,Aa,lda,k,k,reset,transl)
            !
            DO ix = 1, Ninc
              incx = Inc(ix)
              lx = ABS(incx)*n
              !
              !                 Generate the vector X.
              !
              transl = HALF
              CALL CMAKE2('GE',' ',' ',1,n,X,1,Xx,ABS(incx),0,n-1,reset,transl)
              IF ( n>1 ) THEN
                X(n/2) = ZERO
                Xx(1+ABS(incx)*(n/2-1)) = ZERO
              ENDIF
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
                    CALL CMAKE2('GE',' ',' ',1,n,Y,1,Yy,ABS(incy),0,n-1,reset,&
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
                    ENDDO
                    ldas = lda
                    DO i = 1, lx
                      Xs(i) = Xx(i)
                    ENDDO
                    incxs = incx
                    bls = beta
                    DO i = 1, ly
                      Ys(i) = Yy(i)
                    ENDDO
                    incys = incy
                    !
                    !                          Call the subroutine.
                    !
                    IF ( full ) THEN
                      CALL CHEMV(uplo,n,alpha,Aa,lda,Xx,incx,beta,Yy,incy)
                    ELSEIF ( banded ) THEN
                      CALL CHBMV(uplo,n,k,alpha,Aa,lda,Xx,incx,beta,Yy,incy)
                    ELSEIF ( packed ) THEN
                      CALL CHPMV(uplo,n,alpha,Aa,Xx,incx,beta,Yy,incy)
                    ENDIF
                    !
                    !                          Check if error-exit was taken incorrectly.
                    !
                    IF ( NUMXER(nerr)/=0 ) THEN
                      IF ( Kprint>=2 ) WRITE (Nout,FMT=99007)
                      Fatal = .TRUE.
                    ENDIF
                    !
                    !                          See what data changed inside subroutines.
                    !
                    isame(1) = uplo==uplos
                    isame(2) = ns==n
                    IF ( full ) THEN
                      isame(3) = als==alpha
                      isame(4) = LCE(As,Aa,laa)
                      isame(5) = ldas==lda
                      isame(6) = LCE(Xs,Xx,lx)
                      isame(7) = incxs==incx
                      isame(8) = bls==beta
                      IF ( null ) THEN
                        isame(9) = LCE(Ys,Yy,ly)
                      ELSE
                        isame(9) = LCERES('GE',' ',1,n,Ys,Yy,ABS(incy))
                      ENDIF
                      isame(10) = incys==incy
                    ELSEIF ( banded ) THEN
                      isame(3) = ks==k
                      isame(4) = als==alpha
                      isame(5) = LCE(As,Aa,laa)
                      isame(6) = ldas==lda
                      isame(7) = LCE(Xs,Xx,lx)
                      isame(8) = incxs==incx
                      isame(9) = bls==beta
                      IF ( null ) THEN
                        isame(10) = LCE(Ys,Yy,ly)
                      ELSE
                        isame(10) = LCERES('GE',' ',1,n,Ys,Yy,ABS(incy))
                      ENDIF
                      isame(11) = incys==incy
                    ELSEIF ( packed ) THEN
                      isame(3) = als==alpha
                      isame(4) = LCE(As,Aa,laa)
                      isame(5) = LCE(Xs,Xx,lx)
                      isame(6) = incxs==incx
                      isame(7) = bls==beta
                      IF ( null ) THEN
                        isame(8) = LCE(Ys,Yy,ly)
                      ELSE
                        isame(8) = LCERES('GE',' ',1,n,Ys,Yy,ABS(incy))
                      ENDIF
                      isame(9) = incys==incy
                    ENDIF
                    !
                    !                          If data was incorrectly changed, report and
                    !                          return.
                    !
                    DO i = 1, nargs
                      IF ( .NOT.isame(i) ) THEN
                        Fatal = .TRUE.
                        IF ( Kprint>=2 ) WRITE (Nout,FMT=99002) i
                      ENDIF
                    ENDDO
                    !
                    ftl = .FALSE.
                    IF ( .NOT.null ) THEN
                      !
                      !                             Check the result.
                      !
                      CALL CMVCH('N',n,n,alpha,A,Nmax,X,incx,beta,Y,incy,Yt,G,&
                        Yy,Eps,err,ftl,Nout,.TRUE.,Kprint)
                      errmax = MAX(errmax,err)
                    ENDIF
                    IF ( ftl ) THEN
                      Fatal = .TRUE.
                      IF ( Kprint>=3 ) THEN
                        WRITE (Nout,FMT=99004) Sname
                        IF ( full ) THEN
                          WRITE (Nout,FMT=99006) nc, Sname, uplo, n, &
                            alpha, lda, incx, beta, incy
                        ELSEIF ( banded ) THEN
                          WRITE (Nout,FMT=99005) nc, Sname, uplo, n, &
                            alpha, lda, incx, beta, incy
                        ELSEIF ( packed ) THEN
                          WRITE (Nout,FMT=99005) nc, Sname, uplo, n, &
                            alpha, incx, beta, incy
                        ENDIF
                      ENDIF
                    ENDIF
                    !
                  ENDDO
                  !
                ENDDO
                !
              ENDDO
              !
            ENDDO
            !
          ENDDO
        ENDIF
        !
      ENDDO
      !
    ENDDO
    !
    !     Report result.
    !
    IF ( .NOT.Fatal ) THEN
      IF ( Kprint>=3 ) THEN
        IF ( errmax<Thresh ) THEN
          WRITE (Nout,FMT=99001) Sname, nc
        ELSE
          WRITE (Nout,FMT=99003) Sname, nc, errmax
        ENDIF
      ENDIF
    ENDIF
    RETURN
    !
    99001 FORMAT (' ',A6,' PASSED THE COMPUTATIONAL TESTS (',I6,' CALL','S)')
    99002 FORMAT (' ******* FATAL ERROR - PARAMETER NUMBER ',I2,' WAS CH',&
      'ANGED INCORRECTLY *******')
    99003 FORMAT (' ',A6,' COMPLETED THE COMPUTATIONAL TESTS (',I6,' C','ALLS)',&
      /' ******* BUT WITH MAXIMUM TEST RATIO',F8.2,' - SUSPECT *******')
    99004 FORMAT (' ******* ',A6,' FAILED ON CALL NUMBER:')
    99005 FORMAT (1X,I6,': ',A6,'(''',A1,''',',I3,',(',F4.1,',',F4.1,'), AP, X,',I2,&
      ',(',F4.1,',',F4.1,'), Y,',I2,')                .')
    99006 FORMAT (1X,I6,': ',A6,'(''',A1,''',',2(I3,','),'(',F4.1,',',F4.1,'), A,',&
      I3,', X,',I2,',(',F4.1,',',F4.1,'), Y,',I2,')         .')
    99007 FORMAT (' ******* FATAL ERROR - ERROR-EXIT TAKEN ON VALID CALL *',&
      '******')
    !
    !     End of CCHK22.
    !
  END SUBROUTINE CCHK22
  !** CCHK23
  SUBROUTINE CCHK23(Sname,Eps,Thresh,Nout,Kprint,Fatal,Nidim,Idim,Nalf,Alf,&
      Nbet,Bet,Nmax,A,Aa,As,B,Bb,Bs,C,Cc,Cs,Ct,G)
    IMPLICIT NONE
    !>
    !***
    !  Quick check for CHEMM and CSYMM.
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
    !  Quick check for CHEMM and CSYMM.
    !
    !  Auxiliary routine for test program for Level 3 Blas.
    !***
    ! **References:**  (NONE)
    !***
    ! **Routines called:**  CHEMM, CMAKE3, CMMCH, CSYMM, LCE, LCERES, NUMXER

    !* REVISION HISTORY  (YYMMDD)
    !   890208  DATE WRITTEN
    !   910619  Modified to meet SLATEC code and prologue standards.  (BKS)

    !     .. Parameters ..
    COMPLEX, PARAMETER :: ZERO = (0.0,0.0)
    REAL, PARAMETER :: RZERO = 0.0
    !     .. Scalar Arguments ..
    LOGICAL Fatal
    REAL Eps, Thresh
    INTEGER Kprint, Nalf, Nbet, Nidim, Nmax, Nout
    CHARACTER(6) :: Sname
    !     .. Array Arguments ..
    COMPLEX A(Nmax,Nmax), Aa(Nmax*Nmax), Alf(Nalf), As(Nmax*Nmax), &
      B(Nmax,Nmax), Bb(Nmax*Nmax), Bet(Nbet), Bs(Nmax*Nmax), &
      C(Nmax,Nmax), Cc(Nmax*Nmax), Cs(Nmax*Nmax), Ct(Nmax)
    REAL G(Nmax)
    INTEGER Idim(Nidim)
    !     .. Local Scalars ..
    COMPLEX alpha, als, beta, bls
    REAL err, errmax
    INTEGER i, ia, ib, ics, icu, im, in, laa, lbb, lcc, lda, ldas, &
      ldb, ldbs, ldc, ldcs, m, ms, n, na, nargs, nc, nerr, ns
    LOGICAL conj, ftl, left, null, reset
    CHARACTER :: side, sides, uplo, uplos
    !     .. Local Arrays ..
    LOGICAL isame(13)
    !     .. External Functions ..
    INTEGER, EXTERNAL :: NUMXER
    !     .. External Subroutines ..
    EXTERNAL :: CHEMM, CSYMM
    !     .. Intrinsic Functions ..
    INTRINSIC ABS, MAX, MIN
    !     .. Data statements ..
    CHARACTER(2), PARAMETER :: ichs = 'LR', ichu = 'UL'
    !* FIRST EXECUTABLE STATEMENT  CCHK23
    conj = Sname(2:3)=='HE'
    !
    nargs = 12
    nc = 0
    reset = .TRUE.
    errmax = RZERO
    !
    DO im = 1, Nidim
      m = Idim(im)
      !
      DO in = 1, Nidim
        n = Idim(in)
        !           Set LDC to 1 more than minimum value if room.
        ldc = m
        IF ( ldc<Nmax ) ldc = ldc + 1
        !           Skip tests if not enough room.
        IF ( ldc<=Nmax ) THEN
          lcc = ldc*n
          null = n<=0 .OR. m<=0
          !           Set LDB to 1 more than minimum value if room.
          ldb = m
          IF ( ldb<Nmax ) ldb = ldb + 1
          !           Skip tests if not enough room.
          IF ( ldb<=Nmax ) THEN
            lbb = ldb*n
            !
            !           Generate the matrix B.
            !
            CALL CMAKE3('GE',' ',' ',m,n,B,Nmax,Bb,ldb,reset,ZERO)
            !
            DO ics = 1, 2
              side = ichs(ics:ics)
              left = side=='L'
              !
              IF ( left ) THEN
                na = m
              ELSE
                na = n
              ENDIF
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
                  !                 Generate the hermitian or symmetric matrix A.
                  !
                  CALL CMAKE3(Sname(2:3),uplo,' ',na,na,A,Nmax,Aa,lda,reset,ZERO)
                  !
                  DO ia = 1, Nalf
                    alpha = Alf(ia)
                    !
                    DO ib = 1, Nbet
                      beta = Bet(ib)
                      !
                      !                       Generate the matrix C.
                      !
                      CALL CMAKE3('GE',' ',' ',m,n,C,Nmax,Cc,ldc,reset,ZERO)
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
                      ENDDO
                      ldas = lda
                      DO i = 1, lbb
                        Bs(i) = Bb(i)
                      ENDDO
                      ldbs = ldb
                      bls = beta
                      DO i = 1, lcc
                        Cs(i) = Cc(i)
                      ENDDO
                      ldcs = ldc
                      !
                      !                       Call the subroutine.
                      !
                      IF ( conj ) THEN
                        CALL CHEMM(side,uplo,m,n,alpha,Aa,lda,Bb,ldb,beta,Cc,ldc)
                      ELSE
                        CALL CSYMM(side,uplo,m,n,alpha,Aa,lda,Bb,ldb,beta,Cc,ldc)
                      ENDIF
                      !
                      !                       Check if error-exit was taken incorrectly.
                      !
                      IF ( NUMXER(nerr)/=0 ) THEN
                        IF ( Kprint>=2 ) WRITE (Nout,FMT=99006)
                        Fatal = .TRUE.
                      ENDIF
                      !
                      !                       See what data changed inside subroutines.
                      !
                      isame(1) = sides==side
                      isame(2) = uplos==uplo
                      isame(3) = ms==m
                      isame(4) = ns==n
                      isame(5) = als==alpha
                      isame(6) = LCE(As,Aa,laa)
                      isame(7) = ldas==lda
                      isame(8) = LCE(Bs,Bb,lbb)
                      isame(9) = ldbs==ldb
                      isame(10) = bls==beta
                      IF ( null ) THEN
                        isame(11) = LCE(Cs,Cc,lcc)
                      ELSE
                        isame(11) = LCERES('GE',' ',m,n,Cs,Cc,ldc)
                      ENDIF
                      isame(12) = ldcs==ldc
                      !
                      !                       If data was incorrectly changed, report and
                      !                       return.
                      !
                      DO i = 1, nargs
                        IF ( .NOT.isame(i) ) THEN
                          Fatal = .TRUE.
                          IF ( Kprint>=2 ) WRITE (Nout,FMT=99002) i
                        ENDIF
                      ENDDO
                      !
                      ftl = .FALSE.
                      IF ( .NOT.null ) THEN
                        !
                        !                          Check the result.
                        !
                        IF ( left ) THEN
                          CALL CMMCH('N','N',m,n,m,alpha,A,Nmax,B,Nmax,beta,C,&
                            Nmax,Ct,G,Cc,ldc,Eps,err,ftl,Nout,.TRUE.,Kprint)
                        ELSE
                          CALL CMMCH('N','N',m,n,n,alpha,B,Nmax,A,Nmax,beta,C,&
                            Nmax,Ct,G,Cc,ldc,Eps,err,ftl,Nout,.TRUE.,Kprint)
                        ENDIF
                        errmax = MAX(errmax,err)
                      ENDIF
                      !
                      IF ( ftl ) THEN
                        Fatal = .TRUE.
                        IF ( Kprint>=3 ) THEN
                          WRITE (Nout,FMT=99004) Sname
                          WRITE (Nout,FMT=99005) nc, Sname, side, uplo, &
                            m, n, alpha, lda, ldb, beta, ldc
                        ENDIF
                      ENDIF
                    ENDDO
                    !
                  ENDDO
                  !
                ENDDO
              ENDIF
              !
            ENDDO
          ENDIF
        ENDIF
        !
      ENDDO
      !
    ENDDO
    !
    !     Report result.
    !
    IF ( .NOT.Fatal ) THEN
      IF ( Kprint>=3 ) THEN
        IF ( errmax<Thresh ) THEN
          WRITE (Nout,FMT=99001) Sname, nc
        ELSE
          WRITE (Nout,FMT=99003) Sname, nc, errmax
        ENDIF
      ENDIF
    ENDIF
    RETURN
    !
    99001 FORMAT (' ',A6,' PASSED THE COMPUTATIONAL TESTS (',I6,' CALL','S)')
    99002 FORMAT (' ******* FATAL ERROR - PARAMETER NUMBER ',I2,' WAS CH',&
      'ANGED INCORRECTLY *******')
    99003 FORMAT (' ',A6,' COMPLETED THE COMPUTATIONAL TESTS (',I6,' C','ALLS)',&
      /' ******* BUT WITH MAXIMUM TEST RATIO',F8.2,' - SUSPECT *******')
    99004 FORMAT (' ******* ',A6,' FAILED ON CALL NUMBER:')
    99005 FORMAT (1X,I6,': ',A6,'(',2('''',A1,''','),2(I3,','),'(',F4.1,',',F4.1,&
      '), A,',I3,', B,',I3,',(',F4.1,',',F4.1,'), C,',I3,')    .')
    99006 FORMAT (' ******* FATAL ERROR - ERROR-EXIT TAKEN ON VALID CALL *',&
      '******')
    !
    !     End of CCHK23.
    !
  END SUBROUTINE CCHK23
  !** CCHK32
  SUBROUTINE CCHK32(Sname,Eps,Thresh,Nout,Kprint,Fatal,Nidim,Idim,Nkb,Kb,&
      Ninc,Inc,Nmax,Incmax,A,Aa,As,X,Xx,Xs,Xt,G,Z)
    IMPLICIT NONE
    !>
    !***
    !  Quick check for CTRMV, CTBMV, CTPMV, CTRSV, CTBSV and
    !            CTPSV.
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
    !  Quick check for CTRMV, CTBMV, CTPMV, CTRSV, CTBSV and CTPSV.
    !
    !  Auxiliary routine for test program for Level 2 Blas.
    !***
    ! **References:**  (NONE)
    !***
    ! **Routines called:**  CMAKE2, CMVCH, CTBMV, CTBSV, CTPMV, CTPSV, CTRMV,
    !                    CTRSV, LCE, LCERES, NUMXER

    !* REVISION HISTORY  (YYMMDD)
    !   870810  DATE WRITTEN
    !   910619  Modified to meet SLATEC code and prologue standards.  (BKS)

    !     .. Parameters ..
    COMPLEX, PARAMETER :: ZERO = (0.0,0.0), HALF = (0.5,0.0), ONE = (1.0,0.0)
    REAL, PARAMETER :: RZERO = 0.0
    !     .. Scalar Arguments ..
    LOGICAL Fatal
    REAL Eps, Thresh
    INTEGER Incmax, Kprint, Nidim, Ninc, Nkb, Nmax, Nout
    CHARACTER(6) :: Sname
    !     .. Array Arguments ..
    COMPLEX A(Nmax,Nmax), Aa(Nmax*Nmax), As(Nmax*Nmax), X(Nmax), &
      Xs(Nmax*Incmax), Xt(Nmax), Xx(Nmax*Incmax), Z(Nmax)
    REAL G(Nmax)
    INTEGER Idim(Nidim), Inc(Ninc), Kb(Nkb)
    !     .. Local Scalars ..
    COMPLEX transl
    REAL err, errmax
    INTEGER i, icd, ict, icu, ik, in, incx, incxs, ix, k, ks, laa, &
      lda, ldas, lx, n, nargs, nc, nerr, nk, ns
    LOGICAL banded, ftl, full, null, packed, reset
    CHARACTER :: diag, diags, trans, transs, uplo, uplos
    !     .. Local Arrays ..
    LOGICAL isame(13)
    !     .. External Functions ..
    INTEGER, EXTERNAL :: NUMXER
    !     .. External Subroutines ..
    EXTERNAL :: CTBMV, CTBSV, CTPMV, CTPSV, CTRMV, CTRSV
    !     .. Intrinsic Functions ..
    INTRINSIC ABS, MAX
    !     .. Data statements ..
    CHARACTER(2), PARAMETER :: ichu = 'UL', ichd = 'UN'
    CHARACTER(3), PARAMETER :: icht = 'NTC'
    !* FIRST EXECUTABLE STATEMENT  CCHK32
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
    ENDIF
    !
    nc = 0
    reset = .TRUE.
    errmax = RZERO
    !     Set up zero vector for CMVCH.
    DO i = 1, Nmax
      Z(i) = ZERO
    ENDDO
    !
    DO in = 1, Nidim
      n = Idim(in)
      !
      IF ( banded ) THEN
        nk = Nkb
      ELSE
        nk = 1
      ENDIF
      DO ik = 1, nk
        IF ( banded ) THEN
          k = Kb(ik)
        ELSE
          k = n - 1
        ENDIF
        !           Set LDA to 1 more than minimum value if room.
        IF ( banded ) THEN
          lda = k + 1
        ELSE
          lda = n
        ENDIF
        IF ( lda<Nmax ) lda = lda + 1
        !           Skip tests if not enough room.
        IF ( lda<=Nmax ) THEN
          IF ( packed ) THEN
            laa = (n*(n+1))/2
          ELSE
            laa = lda*n
          ENDIF
          null = n<=0
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
                CALL CMAKE2(Sname(2:3),uplo,diag,n,n,A,Nmax,Aa,lda,k,k,reset,&
                  transl)
                !
                DO ix = 1, Ninc
                  incx = Inc(ix)
                  lx = ABS(incx)*n
                  !
                  !                       Generate the vector X.
                  !
                  transl = HALF
                  CALL CMAKE2('GE',' ',' ',1,n,X,1,Xx,ABS(incx),0,n-1,reset,&
                    transl)
                  IF ( n>1 ) THEN
                    X(n/2) = ZERO
                    Xx(1+ABS(incx)*(n/2-1)) = ZERO
                  ENDIF
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
                  ENDDO
                  ldas = lda
                  DO i = 1, lx
                    Xs(i) = Xx(i)
                  ENDDO
                  incxs = incx
                  !
                  !                       Call the subroutine.
                  !
                  IF ( Sname(4:5)=='MV' ) THEN
                    IF ( full ) THEN
                      CALL CTRMV(uplo,trans,diag,n,Aa,lda,Xx,incx)
                    ELSEIF ( banded ) THEN
                      CALL CTBMV(uplo,trans,diag,n,k,Aa,lda,Xx,incx)
                    ELSEIF ( packed ) THEN
                      CALL CTPMV(uplo,trans,diag,n,Aa,Xx,incx)
                    ENDIF
                  ELSEIF ( Sname(4:5)=='SV' ) THEN
                    IF ( full ) THEN
                      CALL CTRSV(uplo,trans,diag,n,Aa,lda,Xx,incx)
                    ELSEIF ( banded ) THEN
                      CALL CTBSV(uplo,trans,diag,n,k,Aa,lda,Xx,incx)
                    ELSEIF ( packed ) THEN
                      CALL CTPSV(uplo,trans,diag,n,Aa,Xx,incx)
                    ENDIF
                  ENDIF
                  !
                  !                       Check if error-exit was taken incorrectly.
                  !
                  IF ( NUMXER(nerr)/=0 ) THEN
                    IF ( Kprint>=2 ) WRITE (Nout,FMT=99007)
                    Fatal = .TRUE.
                  ENDIF
                  !
                  !                       See what data changed inside subroutines.
                  !
                  isame(1) = uplo==uplos
                  isame(2) = trans==transs
                  isame(3) = diag==diags
                  isame(4) = ns==n
                  IF ( full ) THEN
                    isame(5) = LCE(As,Aa,laa)
                    isame(6) = ldas==lda
                    IF ( null ) THEN
                      isame(7) = LCE(Xs,Xx,lx)
                    ELSE
                      isame(7) = LCERES('GE',' ',1,n,Xs,Xx,ABS(incx))
                    ENDIF
                    isame(8) = incxs==incx
                  ELSEIF ( banded ) THEN
                    isame(5) = ks==k
                    isame(6) = LCE(As,Aa,laa)
                    isame(7) = ldas==lda
                    IF ( null ) THEN
                      isame(8) = LCE(Xs,Xx,lx)
                    ELSE
                      isame(8) = LCERES('GE',' ',1,n,Xs,Xx,ABS(incx))
                    ENDIF
                    isame(9) = incxs==incx
                  ELSEIF ( packed ) THEN
                    isame(5) = LCE(As,Aa,laa)
                    IF ( null ) THEN
                      isame(6) = LCE(Xs,Xx,lx)
                    ELSE
                      isame(6) = LCERES('GE',' ',1,n,Xs,Xx,ABS(incx))
                    ENDIF
                    isame(7) = incxs==incx
                  ENDIF
                  !
                  !                       If data was incorrectly changed, report and
                  !                       return.
                  !
                  DO i = 1, nargs
                    IF ( .NOT.isame(i) ) THEN
                      Fatal = .TRUE.
                      IF ( Kprint>=2 ) WRITE (Nout,FMT=99002) i
                    ENDIF
                  ENDDO
                  !
                  ftl = .FALSE.
                  IF ( .NOT.null ) THEN
                    IF ( Sname(4:5)=='MV' ) THEN
                      !
                      !                             Check the result.
                      !
                      CALL CMVCH(trans,n,n,ONE,A,Nmax,X,incx,ZERO,Z,incx,Xt,G,&
                        Xx,Eps,err,ftl,Nout,.TRUE.,Kprint)
                    ELSEIF ( Sname(4:5)=='SV' ) THEN
                      !
                      !                             Compute approximation to original vector.
                      !
                      DO i = 1, n
                        Z(i) = Xx(1+(i-1)*ABS(incx))
                        Xx(1+(i-1)*ABS(incx)) = X(i)
                      ENDDO
                      CALL CMVCH(trans,n,n,ONE,A,Nmax,Z,incx,ZERO,X,incx,Xt,G,&
                        Xx,Eps,err,ftl,Nout,.FALSE.,Kprint)
                    ENDIF
                    errmax = MAX(errmax,err)
                  ENDIF
                  IF ( ftl ) THEN
                    Fatal = .TRUE.
                    IF ( Kprint>=3 ) THEN
                      WRITE (Nout,FMT=99004) Sname
                      IF ( full ) THEN
                        WRITE (Nout,FMT=99006) nc, Sname, uplo, trans, &
                          diag, n, lda, incx
                      ELSEIF ( banded ) THEN
                        WRITE (Nout,FMT=99005) nc, Sname, uplo, trans, &
                          diag, n, k, lda, incx
                      ELSEIF ( packed ) THEN
                        WRITE (Nout,FMT=99005) nc, Sname, uplo, trans, diag, n, incx
                      ENDIF
                    ENDIF
                  ENDIF
                  !
                ENDDO
                !
              ENDDO
              !
            ENDDO
            !
          ENDDO
        ENDIF
        !
      ENDDO
      !
    ENDDO
    !
    !     Report result.
    !
    IF ( .NOT.Fatal ) THEN
      IF ( Kprint>=3 ) THEN
        IF ( errmax<Thresh ) THEN
          WRITE (Nout,FMT=99001) Sname, nc
        ELSE
          WRITE (Nout,FMT=99003) Sname, nc, errmax
        ENDIF
      ENDIF
    ENDIF
    RETURN
    !
    99001 FORMAT (' ',A6,' PASSED THE COMPUTATIONAL TESTS (',I6,' CALL','S)')
    99002 FORMAT (' ******* FATAL ERROR - PARAMETER NUMBER ',I2,' WAS CH',&
      'ANGED INCORRECTLY *******')
    99003 FORMAT (' ',A6,' COMPLETED THE COMPUTATIONAL TESTS (',I6,' C','ALLS)',&
      /' ******* BUT WITH MAXIMUM TEST RATIO',F8.2,' - SUSPECT *******')
    99004 FORMAT (' ******* ',A6,' FAILED ON CALL NUMBER:')
    99005 FORMAT (1X,I6,': ',A6,'(',3('''',A1,''','),I3,', AP, ','X,',I2,&
      ')                                      .')
    99006 FORMAT (1X,I6,': ',A6,'(',3('''',A1,''','),2(I3,','),' A,',I3,', X,',I2,&
      ')                               .')
    99007 FORMAT (' ******* FATAL ERROR - ERROR-EXIT TAKEN ON VALID CALL *',&
      '******')
    !
    !     End of CCHK32.
    !
  END SUBROUTINE CCHK32
  !** CCHK33
  SUBROUTINE CCHK33(Sname,Eps,Thresh,Nout,Kprint,Fatal,Nidim,Idim,Nalf,Alf,&
      Nmax,A,Aa,As,B,Bb,Bs,Ct,G,C)
    IMPLICIT NONE
    !>
    !***
    !  Quick check for CTRMM and CTRSM.
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
    !  Quick check for CTRMM and CTRSM.
    !
    !  Auxiliary routine for test program for Level 3 Blas.
    !***
    ! **References:**  (NONE)
    !***
    ! **Routines called:**  CMAKE3, CMMCH, CTRMM, CTRSM, LCE, LCERES, NUMXER

    !* REVISION HISTORY  (YYMMDD)
    !   890208  DATE WRITTEN
    !   910619  Modified to meet SLATEC code and prologue standards.  (BKS)

    !     .. Parameters ..
    COMPLEX, PARAMETER :: ZERO = (0.0,0.0), ONE = (1.0,0.0)
    REAL, PARAMETER :: RZERO = 0.0
    !     .. Scalar Arguments ..
    LOGICAL Fatal
    REAL Eps, Thresh
    INTEGER Kprint, Nalf, Nidim, Nmax, Nout
    CHARACTER(6) :: Sname
    !     .. Array Arguments ..
    COMPLEX A(Nmax,Nmax), Aa(Nmax*Nmax), Alf(Nalf), As(Nmax*Nmax), &
      B(Nmax,Nmax), Bb(Nmax*Nmax), Bs(Nmax*Nmax), C(Nmax,Nmax), Ct(Nmax)
    REAL G(Nmax)
    INTEGER Idim(Nidim)
    !     .. Local Scalars ..
    COMPLEX alpha, als
    REAL err, errmax
    INTEGER i, ia, icd, ics, ict, icu, im, in, j, laa, lbb, lda, &
      ldas, ldb, ldbs, m, ms, n, na, nargs, nc, nerr, ns
    LOGICAL ftl, left, null, reset
    CHARACTER :: diag, diags, side, sides, tranas, transa, uplo, uplos
    !     .. Local Arrays ..
    LOGICAL isame(13)
    !     .. External Functions ..
    INTEGER, EXTERNAL :: NUMXER
    !     .. External Subroutines ..
    EXTERNAL :: CTRMM, CTRSM
    !     .. Intrinsic Functions ..
    INTRINSIC ABS, MAX, MIN
    !     .. Data statements ..
    CHARACTER(2), PARAMETER :: ichs = 'LR', ichu = 'UL', ichd = 'UN'
    CHARACTER(3), PARAMETER :: icht = 'NTC'
    !* FIRST EXECUTABLE STATEMENT  CCHK33
    nargs = 11
    nc = 0
    reset = .TRUE.
    errmax = RZERO
    !     Set up zero matrix for CMMCH.
    DO j = 1, Nmax
      DO i = 1, Nmax
        C(i,j) = ZERO
      ENDDO
    ENDDO
    !
    DO im = 1, Nidim
      m = Idim(im)
      !
      DO in = 1, Nidim
        n = Idim(in)
        !           Set LDB to 1 more than minimum value if room.
        ldb = m
        IF ( ldb<Nmax ) ldb = ldb + 1
        !           Skip tests if not enough room.
        IF ( ldb<=Nmax ) THEN
          lbb = ldb*n
          null = m<=0 .OR. n<=0
          !
          DO ics = 1, 2
            side = ichs(ics:ics)
            left = side=='L'
            IF ( left ) THEN
              na = m
            ELSE
              na = n
            ENDIF
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
                    CALL CMAKE3('TR',uplo,diag,na,na,A,Nmax,Aa,lda,reset,ZERO)
                    !
                    !                          Generate the matrix B.
                    !
                    CALL CMAKE3('GE',' ',' ',m,n,B,Nmax,Bb,ldb,reset,ZERO)
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
                    ENDDO
                    ldas = lda
                    DO i = 1, lbb
                      Bs(i) = Bb(i)
                    ENDDO
                    ldbs = ldb
                    !
                    !                          Call the subroutine.
                    !
                    IF ( Sname(4:5)=='MM' ) THEN
                      CALL CTRMM(side,uplo,transa,diag,m,n,alpha,Aa,lda,Bb,ldb)
                    ELSEIF ( Sname(4:5)=='SM' ) THEN
                      CALL CTRSM(side,uplo,transa,diag,m,n,alpha,Aa,lda,Bb,ldb)
                    ENDIF
                    !
                    !                          Check if error-exit was taken incorrectly.
                    !
                    IF ( NUMXER(nerr)/=0 ) THEN
                      IF ( Kprint>=2 ) WRITE (Nout,FMT=99006)
                      Fatal = .TRUE.
                    ENDIF
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
                    isame(8) = LCE(As,Aa,laa)
                    isame(9) = ldas==lda
                    IF ( null ) THEN
                      isame(10) = LCE(Bs,Bb,lbb)
                    ELSE
                      isame(10) = LCERES('GE',' ',m,n,Bs,Bb,ldb)
                    ENDIF
                    isame(11) = ldbs==ldb
                    !
                    !                          If data was incorrectly changed, report and
                    !                          return.
                    !
                    DO i = 1, nargs
                      IF ( .NOT.isame(i) ) THEN
                        Fatal = .TRUE.
                        IF ( Kprint>=2 ) WRITE (Nout,FMT=99002) i
                      ENDIF
                    ENDDO
                    !
                    ftl = .FALSE.
                    IF ( .NOT.null ) THEN
                      IF ( Sname(4:5)=='MM' ) THEN
                        !
                        !                                Check the result.
                        !
                        IF ( left ) THEN
                          CALL CMMCH(transa,'N',m,n,m,alpha,A,Nmax,B,Nmax,&
                            ZERO,C,Nmax,Ct,G,Bb,ldb,Eps,err,ftl,Nout,.TRUE.,Kprint)
                        ELSE
                          CALL CMMCH('N',transa,m,n,n,alpha,B,Nmax,A,Nmax,&
                            ZERO,C,Nmax,Ct,G,Bb,ldb,Eps,err,ftl,Nout,.TRUE.,Kprint)
                        ENDIF
                      ELSEIF ( Sname(4:5)=='SM' ) THEN
                        !
                        !                                Compute approximation to original
                        !                                matrix.
                        !
                        DO j = 1, n
                          DO i = 1, m
                            C(i,j) = Bb(i+(j-1)*ldb)
                            Bb(i+(j-1)*ldb) = alpha*B(i,j)
                          ENDDO
                        ENDDO
                        !
                        IF ( left ) THEN
                          CALL CMMCH(transa,'N',m,n,m,ONE,A,Nmax,C,Nmax,ZERO,&
                            B,Nmax,Ct,G,Bb,ldb,Eps,err,ftl,Nout,.FALSE.,Kprint)
                        ELSE
                          CALL CMMCH('N',transa,m,n,n,ONE,C,Nmax,A,Nmax,ZERO,&
                            B,Nmax,Ct,G,Bb,ldb,Eps,err,ftl,Nout,.FALSE.,Kprint)
                        ENDIF
                      ENDIF
                      errmax = MAX(errmax,err)
                    ENDIF
                    IF ( ftl ) THEN
                      Fatal = .TRUE.
                      IF ( Kprint>=3 ) THEN
                        WRITE (Nout,FMT=99004) Sname
                        WRITE (Nout,FMT=99005) nc, Sname, side, uplo, &
                          transa, diag, m, n, alpha, lda, ldb
                      ENDIF
                    ENDIF
                    !
                  ENDDO
                  !
                ENDDO
                !
              ENDDO
              !
            ENDDO
            !
          ENDDO
        ENDIF
        !
      ENDDO
      !
    ENDDO
    !
    !     Report result.
    !
    IF ( .NOT.Fatal ) THEN
      IF ( Kprint>=3 ) THEN
        IF ( errmax<Thresh ) THEN
          WRITE (Nout,FMT=99001) Sname, nc
        ELSE
          WRITE (Nout,FMT=99003) Sname, nc, errmax
        ENDIF
      ENDIF
    ENDIF
    RETURN
    !
    99001 FORMAT (' ',A6,' PASSED THE COMPUTATIONAL TESTS (',I6,' CALL','S)')
    99002 FORMAT (' ******* FATAL ERROR - PARAMETER NUMBER ',I2,' WAS CH',&
      'ANGED INCORRECTLY *******')
    99003 FORMAT (' ',A6,' COMPLETED THE COMPUTATIONAL TESTS (',I6,' C','ALLS)',&
      /' ******* BUT WITH MAXIMUM TEST RATIO',F8.2,' - SUSPECT *******')
    99004 FORMAT (' ******* ',A6,' FAILED ON CALL NUMBER:')
    99005 FORMAT (1X,I6,': ',A6,'(',4('''',A1,''','),2(I3,','),'(',F4.1,',',F4.1,&
      '), A,',I3,', B,',I3,')         ','      .')
    99006 FORMAT (' ******* FATAL ERROR - ERROR-EXIT TAKEN ON VALID CALL *',&
      '******')
    !
    !     End of CCHK33.
    !
  END SUBROUTINE CCHK33
  !** CCHK42
  SUBROUTINE CCHK42(Sname,Eps,Thresh,Nout,Kprint,Fatal,Nidim,Idim,Nalf,Alf,&
      Ninc,Inc,Nmax,Incmax,A,Aa,As,X,Xx,Xs,Y,Yy,Ys,Yt,G,Z)
    IMPLICIT NONE
    !>
    !***
    !  Quick check for CGERC and CGERU.
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
    !  Quick check for CGERC and CGERU.
    !
    !  Auxiliary routine for test program for Level 2 Blas.
    !***
    ! **References:**  (NONE)
    !***
    ! **Routines called:**  CGERC, CGERU, CMAKE2, CMVCH, LCE, LCERES, NUMXER

    !* REVISION HISTORY  (YYMMDD)
    !   870810  DATE WRITTEN
    !   910619  Modified to meet SLATEC code and prologue standards.  (BKS)

    !     .. Parameters ..
    COMPLEX, PARAMETER :: ZERO = (0.0,0.0), HALF = (0.5,0.0), ONE = (1.0,0.0)
    REAL, PARAMETER :: RZERO = 0.0
    !     .. Scalar Arguments ..
    LOGICAL Fatal
    REAL Eps, Thresh
    INTEGER Incmax, Kprint, Nalf, Nidim, Ninc, Nmax, Nout
    CHARACTER(6) :: Sname
    !     .. Array Arguments ..
    COMPLEX A(Nmax,Nmax), Aa(Nmax*Nmax), Alf(Nalf), As(Nmax*Nmax), X(Nmax), &
      Xs(Nmax*Incmax), Xx(Nmax*Incmax), Y(Nmax), Ys(Nmax*Incmax), &
      Yt(Nmax), Yy(Nmax*Incmax), Z(Nmax)
    REAL G(Nmax)
    INTEGER Idim(Nidim), Inc(Ninc)
    !     .. Local Scalars ..
    COMPLEX alpha, als, transl
    REAL err, errmax
    INTEGER i, ia, im, in, incx, incxs, incy, incys, ix, iy, j, &
      laa, lda, ldas, lx, ly, m, ms, n, nargs, nc, nd, nerr, ns
    LOGICAL conj, ftl, null, reset
    !     .. Local Arrays ..
    COMPLEX w(1)
    LOGICAL isame(13)
    !     .. External Functions ..
    INTEGER, EXTERNAL :: NUMXER
    !     .. External Subroutines ..
    EXTERNAL :: CGERC, CGERU
    !     .. Intrinsic Functions ..
    INTRINSIC ABS, CONJG, MAX, MIN
    !* FIRST EXECUTABLE STATEMENT  CCHK42
    conj = Sname(5:5)=='C'
    !     Define the number of arguments.
    nargs = 9
    !
    nc = 0
    reset = .TRUE.
    errmax = RZERO
    !
    DO in = 1, Nidim
      n = Idim(in)
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
          null = n<=0 .OR. m<=0
          !
          DO ix = 1, Ninc
            incx = Inc(ix)
            lx = ABS(incx)*m
            !
            !              Generate the vector X.
            !
            transl = HALF
            CALL CMAKE2('GE',' ',' ',1,m,X,1,Xx,ABS(incx),0,m-1,reset,transl)
            IF ( m>1 ) THEN
              X(m/2) = ZERO
              Xx(1+ABS(incx)*(m/2-1)) = ZERO
            ENDIF
            !
            DO iy = 1, Ninc
              incy = Inc(iy)
              ly = ABS(incy)*n
              !
              !                 Generate the vector Y.
              !
              transl = ZERO
              CALL CMAKE2('GE',' ',' ',1,n,Y,1,Yy,ABS(incy),0,n-1,reset,transl)
              IF ( n>1 ) THEN
                Y(n/2) = ZERO
                Yy(1+ABS(incy)*(n/2-1)) = ZERO
              ENDIF
              !
              DO ia = 1, Nalf
                alpha = Alf(ia)
                !
                !                    Generate the matrix A.
                !
                transl = ZERO
                CALL CMAKE2(Sname(2:3),' ',' ',m,n,A,Nmax,Aa,lda,m-1,n-1,&
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
                ENDDO
                ldas = lda
                DO i = 1, lx
                  Xs(i) = Xx(i)
                ENDDO
                incxs = incx
                DO i = 1, ly
                  Ys(i) = Yy(i)
                ENDDO
                incys = incy
                !
                !                    Call the subroutine.
                !
                IF ( conj ) THEN
                  CALL CGERC(m,n,alpha,Xx,incx,Yy,incy,Aa,lda)
                ELSE
                  CALL CGERU(m,n,alpha,Xx,incx,Yy,incy,Aa,lda)
                ENDIF
                !
                !                    Check if error-exit was taken incorrectly.
                !
                IF ( NUMXER(nerr)/=0 ) THEN
                  IF ( Kprint>=2 ) WRITE (Nout,FMT=99007)
                  Fatal = .TRUE.
                ENDIF
                !
                !                    See what data changed inside subroutine.
                !
                isame(1) = ms==m
                isame(2) = ns==n
                isame(3) = als==alpha
                isame(4) = LCE(Xs,Xx,lx)
                isame(5) = incxs==incx
                isame(6) = LCE(Ys,Yy,ly)
                isame(7) = incys==incy
                IF ( null ) THEN
                  isame(8) = LCE(As,Aa,laa)
                ELSE
                  isame(8) = LCERES('GE',' ',m,n,As,Aa,lda)
                ENDIF
                isame(9) = ldas==lda
                !
                !                    If data was incorrectly changed, report and return.
                !
                DO i = 1, nargs
                  IF ( .NOT.isame(i) ) THEN
                    Fatal = .TRUE.
                    IF ( Kprint>=2 ) WRITE (Nout,FMT=99002) i
                  ENDIF
                ENDDO
                !
                ftl = .FALSE.
                IF ( .NOT.null ) THEN
                  !
                  !                       Check the result column by column.
                  !
                  IF ( incx>0 ) THEN
                    DO i = 1, m
                      Z(i) = X(i)
                    ENDDO
                  ELSE
                    DO i = 1, m
                      Z(i) = X(m-i+1)
                    ENDDO
                  ENDIF
                  DO j = 1, n
                    IF ( incy>0 ) THEN
                      w(1) = Y(j)
                    ELSE
                      w(1) = Y(n-j+1)
                    ENDIF
                    IF ( conj ) w(1) = CONJG(w(1))
                    CALL CMVCH('N',m,1,alpha,Z,Nmax,w,1,ONE,A(1,j),1,Yt,G,&
                      Aa(1+(j-1)*lda),Eps,err,ftl,Nout,.TRUE.,Kprint)
                    errmax = MAX(errmax,err)
                  ENDDO
                ENDIF
                IF ( ftl ) THEN
                  Fatal = .TRUE.
                  IF ( Kprint>=3 ) THEN
                    WRITE (Nout,FMT=99005) j
                    WRITE (Nout,FMT=99004) Sname
                    WRITE (Nout,FMT=99006) nc, Sname, m, n, alpha, incx, incy, lda
                  ENDIF
                ENDIF
                !
              ENDDO
              !
            ENDDO
            !
          ENDDO
        ENDIF
        !
      ENDDO
      !
    ENDDO
    !
    !     Report result.
    !
    IF ( .NOT.Fatal ) THEN
      IF ( Kprint>=3 ) THEN
        IF ( errmax<Thresh ) THEN
          WRITE (Nout,FMT=99001) Sname, nc
        ELSE
          WRITE (Nout,FMT=99003) Sname, nc, errmax
        ENDIF
      ENDIF
    ENDIF
    RETURN
    !
    99001 FORMAT (' ',A6,' PASSED THE COMPUTATIONAL TESTS (',I6,' CALL','S)')
    99002 FORMAT (' ******* FATAL ERROR - PARAMETER NUMBER ',I2,' WAS CH',&
      'ANGED INCORRECTLY *******')
    99003 FORMAT (' ',A6,' COMPLETED THE COMPUTATIONAL TESTS (',I6,' C','ALLS)',&
      /' ******* BUT WITH MAXIMUM TEST RATIO',F8.2,' - SUSPECT *******')
    99004 FORMAT (' ******* ',A6,' FAILED ON CALL NUMBER:')
    99005 FORMAT ('      THESE ARE THE RESULTS FOR COLUMN ',I3)
    99006 FORMAT (1X,I6,': ',A6,'(',2(I3,','),'(',F4.1,',',F4.1,'), X,',I2,', Y,',&
      I2,', A,',I3,')                   ','      .')
    99007 FORMAT (' ******* FATAL ERROR - ERROR-EXIT TAKEN ON VALID CALL *',&
      '******')
    !
    !     End of CCHK42.
    !
  END SUBROUTINE CCHK42
  !** CCHK43
  SUBROUTINE CCHK43(Sname,Eps,Thresh,Nout,Kprint,Fatal,Nidim,Idim,Nalf,Alf,&
      Nbet,Bet,Nmax,A,Aa,As,B,Bb,Bs,C,Cc,Cs,Ct,G)
    IMPLICIT NONE
    !>
    !***
    !  Quick check for CHERK and CSYRK.
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
    !  Quick check for CHERK and CSYRK.
    !
    !  Auxiliary routine for test program for Level 3 Blas.
    !***
    ! **References:**  (NONE)
    !***
    ! **Routines called:**  CHERK, CMAKE3, CMMCH, CSYRK, LCE, LCERES, NUMXER

    !* REVISION HISTORY  (YYMMDD)
    !   890208  DATE WRITTEN
    !   910619  Modified to meet SLATEC code and prologue standards.  (BKS)

    !     .. Parameters ..
    COMPLEX, PARAMETER :: ZERO = (0.0,0.0)
    REAL, PARAMETER :: RZERO = 0.0, RONE = 1.0
    !     .. Scalar Arguments ..
    LOGICAL Fatal
    REAL Eps, Thresh
    INTEGER Kprint, Nalf, Nbet, Nidim, Nmax, Nout
    CHARACTER(6) :: Sname
    !     .. Array Arguments ..
    COMPLEX A(Nmax,Nmax), Aa(Nmax*Nmax), Alf(Nalf), As(Nmax*Nmax), &
      B(Nmax,Nmax), Bb(Nmax*Nmax), Bs(Nmax*Nmax), Bet(Nbet), &
      C(Nmax,Nmax), Cc(Nmax*Nmax), Ct(Nmax), Cs(Nmax*Nmax)
    REAL G(Nmax)
    INTEGER Idim(Nidim)
    !     .. Local Scalars ..
    COMPLEX alpha, als, beta, bets
    REAL err, errmax, ralpha, rals, rbeta, rbets
    INTEGER i, ia, ib, ict, icu, ik, in, j, jc, jj, k, ks, laa, &
      lcc, lda, ldas, ldc, ldcs, lj, ma, n, na, nargs, nc, nerr, ns
    LOGICAL conj, ftl, null, reset, tran, upper
    CHARACTER :: trans, transs, uplo, transt, uplos
    !     .. Local Arrays ..
    LOGICAL isame(13)
    !     .. External Functions ..
    INTEGER, EXTERNAL :: NUMXER
    !     .. External Subroutines ..
    EXTERNAL :: CHERK, CSYRK
    !     .. Intrinsic Functions ..
    INTRINSIC ABS, MAX, MIN
    !     .. Data statements ..
    CHARACTER(2), PARAMETER :: ichu = 'UL', icht = 'NC'
    !* FIRST EXECUTABLE STATEMENT  CCHK43
    conj = Sname(2:3)=='HE'
    !
    nargs = 10
    nc = 0
    reset = .TRUE.
    errmax = RZERO
    !
    DO in = 1, Nidim
      n = Idim(in)
      !        Set LDC to 1 more than minimum value if room.
      ldc = n
      IF ( ldc<Nmax ) ldc = ldc + 1
      !        Skip tests if not enough room.
      IF ( ldc<=Nmax ) THEN
        lcc = ldc*n
        !
        DO ik = 1, Nidim
          k = Idim(ik)
          !
          DO ict = 1, 2
            trans = icht(ict:ict)
            tran = trans=='C'
            IF ( tran.AND..NOT.conj ) trans = 'T'
            IF ( tran ) THEN
              ma = k
              na = n
            ELSE
              ma = n
              na = k
            ENDIF
            !              Set LDA to 1 more than minimum value if room.
            lda = ma
            IF ( lda<Nmax ) lda = lda + 1
            !              Skip tests if not enough room.
            IF ( lda<=Nmax ) THEN
              laa = lda*na
              !
              !              Generate the matrix A.
              !
              CALL CMAKE3('GE',' ',' ',ma,na,A,Nmax,Aa,lda,reset,ZERO)
              !
              DO icu = 1, 2
                uplo = ichu(icu:icu)
                upper = uplo=='U'
                !
                DO ia = 1, Nalf
                  alpha = Alf(ia)
                  IF ( conj ) THEN
                    ralpha = REAL(alpha)
                    alpha = CMPLX(ralpha,RZERO)
                  ENDIF
                  !
                  DO ib = 1, Nbet
                    beta = Bet(ib)
                    IF ( conj ) THEN
                      rbeta = REAL(beta)
                      beta = CMPLX(rbeta,RZERO)
                    ENDIF
                    null = n<=0
                    IF ( conj ) null = null .OR.&
                      ((k<=0.OR.ralpha==RZERO).AND.rbeta==RONE)
                    !
                    !                       Generate the matrix C.
                    !
                    CALL CMAKE3(Sname(2:3),uplo,' ',n,n,C,Nmax,Cc,ldc,reset,ZERO)
                    !
                    nc = nc + 1
                    !
                    !                       Save every datum before calling the subroutine.
                    !
                    uplos = uplo
                    transs = trans
                    ns = n
                    ks = k
                    IF ( conj ) THEN
                      rals = ralpha
                    ELSE
                      als = alpha
                    ENDIF
                    DO i = 1, laa
                      As(i) = Aa(i)
                    ENDDO
                    ldas = lda
                    IF ( conj ) THEN
                      rbets = rbeta
                    ELSE
                      bets = beta
                    ENDIF
                    DO i = 1, lcc
                      Cs(i) = Cc(i)
                    ENDDO
                    ldcs = ldc
                    !
                    !                       Call the subroutine.
                    !
                    IF ( conj ) THEN
                      CALL CHERK(uplo,trans,n,k,ralpha,Aa,lda,rbeta,Cc,ldc)
                    ELSE
                      CALL CSYRK(uplo,trans,n,k,alpha,Aa,lda,beta,Cc,ldc)
                    ENDIF
                    !
                    !                       Check if error-exit was taken incorrectly.
                    !
                    IF ( NUMXER(nerr)/=0 ) THEN
                      IF ( Kprint>=2 ) WRITE (Nout,FMT=99007)
                      Fatal = .TRUE.
                    ENDIF
                    !
                    !                       See what data changed inside subroutines.
                    !
                    isame(1) = uplos==uplo
                    isame(2) = transs==trans
                    isame(3) = ns==n
                    isame(4) = ks==k
                    IF ( conj ) THEN
                      isame(5) = rals==ralpha
                    ELSE
                      isame(5) = als==alpha
                    ENDIF
                    isame(6) = LCE(As,Aa,laa)
                    isame(7) = ldas==lda
                    IF ( conj ) THEN
                      isame(8) = rbets==rbeta
                    ELSE
                      isame(8) = bets==beta
                    ENDIF
                    IF ( null ) THEN
                      isame(9) = LCE(Cs,Cc,lcc)
                    ELSE
                      isame(9) = LCERES(Sname(2:3),uplo,n,n,Cs,Cc,ldc)
                    ENDIF
                    isame(10) = ldcs==ldc
                    !
                    !                       If data was incorrectly changed, report and
                    !                       return.
                    !
                    DO i = 1, nargs
                      IF ( .NOT.isame(i) ) THEN
                        Fatal = .TRUE.
                        IF ( Kprint>=2 ) WRITE (Nout,FMT=99002) i
                      ENDIF
                    ENDDO
                    !
                    IF ( .NOT.null ) THEN
                      !
                      !                          Check the result column by column.
                      !
                      IF ( conj ) THEN
                        transt = 'C'
                      ELSE
                        transt = 'T'
                      ENDIF
                      jc = 1
                      DO j = 1, n
                        IF ( upper ) THEN
                          jj = 1
                          lj = j
                        ELSE
                          jj = j
                          lj = n - j + 1
                        ENDIF
                        IF ( tran ) THEN
                          ftl = .FALSE.
                          CALL CMMCH(transt,'N',lj,1,k,alpha,A(1,jj),Nmax,&
                            A(1,j),Nmax,beta,C(jj,j),Nmax,Ct,G,Cc(jc)&
                            ,ldc,Eps,err,ftl,Nout,.TRUE.,Kprint)
                        ELSE
                          ftl = .FALSE.
                          CALL CMMCH('N',transt,lj,1,k,alpha,A(jj,1),Nmax,&
                            A(j,1),Nmax,beta,C(jj,j),Nmax,Ct,G,Cc(jc)&
                            ,ldc,Eps,err,ftl,Nout,.TRUE.,Kprint)
                        ENDIF
                        IF ( upper ) THEN
                          jc = jc + ldc
                        ELSE
                          jc = jc + ldc + 1
                        ENDIF
                        errmax = MAX(errmax,err)
                        IF ( ftl ) THEN
                          Fatal = .TRUE.
                          IF ( Kprint>=3 ) THEN
                            WRITE (Nout,FMT=99004) Sname
                            IF ( conj ) THEN
                              WRITE (Nout,FMT=99005) nc, Sname, uplo, &
                                trans, n, k, ralpha, lda, rbeta, ldc
                            ELSE
                              WRITE (Nout,FMT=99006) nc, Sname, uplo, &
                                trans, n, k, alpha, lda, beta, ldc
                            ENDIF
                          ENDIF
                        ENDIF
                      ENDDO
                    ENDIF
                    !
                  ENDDO
                  !
                ENDDO
                !
              ENDDO
            ENDIF
            !
          ENDDO
          !
        ENDDO
      ENDIF
      !
    ENDDO
    !
    !     Report result.
    !
    IF ( .NOT.Fatal ) THEN
      IF ( Kprint>=3 ) THEN
        IF ( errmax<Thresh ) THEN
          WRITE (Nout,FMT=99001) Sname, nc
        ELSE
          WRITE (Nout,FMT=99003) Sname, nc, errmax
        ENDIF
      ENDIF
    ENDIF
    RETURN
    !
    99001 FORMAT (' ',A6,' PASSED THE COMPUTATIONAL TESTS (',I6,' CALL','S)')
    99002 FORMAT (' ******* FATAL ERROR - PARAMETER NUMBER ',I2,' WAS CH',&
      'ANGED INCORRECTLY *******')
    99003 FORMAT (' ',A6,' COMPLETED THE COMPUTATIONAL TESTS (',I6,' C','ALLS)',&
      /' ******* BUT WITH MAXIMUM TEST RATIO',F8.2,' - SUSPECT *******')
    99004 FORMAT (' ******* ',A6,' FAILED ON CALL NUMBER:')
    99005 FORMAT (1X,I6,': ',A6,'(',2('''',A1,''','),2(I3,','),F4.1,', A,',I3,',',&
      F4.1,', C,',I3,')               ','          .')
    99006 FORMAT (1X,I6,': ',A6,'(',2('''',A1,''','),2(I3,','),'(',F4.1,',',F4.1,&
      '), A,',I3,',(',F4.1,',',F4.1,'), C,',I3,')          .')
    99007 FORMAT (' ******* FATAL ERROR - ERROR-EXIT TAKEN ON VALID CALL *',&
      '******')
    !
    !     End of CCHK43.
    !
  END SUBROUTINE CCHK43
  !** CCHK52
  SUBROUTINE CCHK52(Sname,Eps,Thresh,Nout,Kprint,Fatal,Nidim,Idim,Nalf,Alf,&
      Ninc,Inc,Nmax,Incmax,A,Aa,As,X,Xx,Xs,Y,Yy,Ys,Yt,G,Z)
    IMPLICIT NONE
    !>
    !***
    !  Quick check for CHER and CHPR.
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
    !  Quick check for CHER and CHPR.
    !
    !  Auxiliary routine for test program for Level 2 Blas.
    !***
    ! **References:**  (NONE)
    !***
    ! **Routines called:**  CHER, CHPR, CMAKE2, CMVCH, LCE, LCERES, NUMXER

    !* REVISION HISTORY  (YYMMDD)
    !   870810  DATE WRITTEN
    !   910619  Modified to meet SLATEC code and prologue standards.  (BKS)

    !     .. Parameters ..
    COMPLEX, PARAMETER :: ZERO = (0.0,0.0), HALF = (0.5,0.0), ONE = (1.0,0.0)
    REAL, PARAMETER :: RZERO = 0.0
    !     .. Scalar Arguments ..
    LOGICAL Fatal
    REAL Eps, Thresh
    INTEGER Incmax, Kprint, Nalf, Nidim, Ninc, Nmax, Nout
    CHARACTER(6) :: Sname
    !     .. Array Arguments ..
    COMPLEX A(Nmax,Nmax), Aa(Nmax*Nmax), Alf(Nalf), As(Nmax*Nmax), X(Nmax), &
      Xs(Nmax*Incmax), Xx(Nmax*Incmax), Y(Nmax), Ys(Nmax*Incmax), &
      Yt(Nmax), Yy(Nmax*Incmax), Z(Nmax)
    REAL G(Nmax)
    INTEGER Idim(Nidim), Inc(Ninc)
    !     .. Local Scalars ..
    COMPLEX alpha, transl
    REAL err, errmax, ralpha, rals
    INTEGER i, ia, ic, in, incx, incxs, ix, j, ja, jj, laa, lda, &
      ldas, lj, lx, n, nargs, nc, nerr, ns
    LOGICAL ftl, full, null, packed, reset, upper
    CHARACTER :: uplo, uplos
    !     .. Local Arrays ..
    COMPLEX w(1)
    LOGICAL isame(13)
    !     .. External Functions ..
    INTEGER, EXTERNAL :: NUMXER
    !     .. External Subroutines ..
    EXTERNAL :: CHER, CHPR
    !     .. Intrinsic Functions ..
    INTRINSIC ABS, CMPLX, CONJG, MAX, REAL
    !     .. Data statements ..
    CHARACTER(2), PARAMETER :: ich = 'UL'
    !* FIRST EXECUTABLE STATEMENT  CCHK52
    full = Sname(3:3)=='E'
    packed = Sname(3:3)=='P'
    !     Define the number of arguments.
    IF ( full ) THEN
      nargs = 7
    ELSEIF ( packed ) THEN
      nargs = 6
    ENDIF
    !
    nc = 0
    reset = .TRUE.
    errmax = RZERO
    !
    DO in = 1, Nidim
      n = Idim(in)
      !        Set LDA to 1 more than minimum value if room.
      lda = n
      IF ( lda<Nmax ) lda = lda + 1
      !        Skip tests if not enough room.
      IF ( lda<=Nmax ) THEN
        IF ( packed ) THEN
          laa = (n*(n+1))/2
        ELSE
          laa = lda*n
        ENDIF
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
            CALL CMAKE2('GE',' ',' ',1,n,X,1,Xx,ABS(incx),0,n-1,reset,transl)
            IF ( n>1 ) THEN
              X(n/2) = ZERO
              Xx(1+ABS(incx)*(n/2-1)) = ZERO
            ENDIF
            !
            DO ia = 1, Nalf
              ralpha = REAL(Alf(ia))
              alpha = CMPLX(ralpha,RZERO)
              null = n<=0 .OR. ralpha==RZERO
              !
              !                 Generate the matrix A.
              !
              transl = ZERO
              CALL CMAKE2(Sname(2:3),uplo,' ',n,n,A,Nmax,Aa,lda,n-1,n-1,reset,&
                transl)
              !
              nc = nc + 1
              !
              !                 Save every datum before calling the subroutine.
              !
              uplos = uplo
              ns = n
              rals = ralpha
              DO i = 1, laa
                As(i) = Aa(i)
              ENDDO
              ldas = lda
              DO i = 1, lx
                Xs(i) = Xx(i)
              ENDDO
              incxs = incx
              !
              !                 Call the subroutine.
              !
              IF ( full ) THEN
                CALL CHER(uplo,n,ralpha,Xx,incx,Aa,lda)
              ELSEIF ( packed ) THEN
                CALL CHPR(uplo,n,ralpha,Xx,incx,Aa)
              ENDIF
              !
              !                 Check if error-exit was taken incorrectly.
              !
              IF ( NUMXER(nerr)/=0 ) THEN
                IF ( Kprint>=2 ) WRITE (Nout,FMT=99007)
                Fatal = .TRUE.
              ENDIF
              !
              !                 See what data changed inside subroutines.
              !
              isame(1) = uplo==uplos
              isame(2) = ns==n
              isame(3) = rals==ralpha
              isame(4) = LCE(Xs,Xx,lx)
              isame(5) = incxs==incx
              IF ( null ) THEN
                isame(6) = LCE(As,Aa,laa)
              ELSE
                isame(6) = LCERES(Sname(2:3),uplo,n,n,As,Aa,lda)
              ENDIF
              IF ( .NOT.packed ) isame(7) = ldas==lda
              !
              !                 If data was incorrectly changed, report and return.
              !
              DO i = 1, nargs
                IF ( .NOT.isame(i) ) THEN
                  Fatal = .TRUE.
                  IF ( Kprint>=2 ) WRITE (Nout,FMT=99002) i
                ENDIF
              ENDDO
              !
              ftl = .FALSE.
              IF ( .NOT.null ) THEN
                !
                !                    Check the result column by column.
                !
                IF ( incx>0 ) THEN
                  DO i = 1, n
                    Z(i) = X(i)
                  ENDDO
                ELSE
                  DO i = 1, n
                    Z(i) = X(n-i+1)
                  ENDDO
                ENDIF
                ja = 1
                DO j = 1, n
                  w(1) = CONJG(Z(j))
                  IF ( upper ) THEN
                    jj = 1
                    lj = j
                  ELSE
                    jj = j
                    lj = n - j + 1
                  ENDIF
                  CALL CMVCH('N',lj,1,alpha,Z(jj),lj,w,1,ONE,A(jj,j),1,Yt,G,&
                    Aa(ja),Eps,err,ftl,Nout,.TRUE.,Kprint)
                  IF ( .NOT.(full) ) THEN
                    ja = ja + lj
                  ELSEIF ( upper ) THEN
                    ja = ja + lda
                  ELSE
                    ja = ja + lda + 1
                  ENDIF
                  errmax = MAX(errmax,err)
                ENDDO
              ENDIF
              IF ( ftl ) THEN
                Fatal = .TRUE.
                IF ( Kprint>=3 ) THEN
                  WRITE (Nout,FMT=99004) Sname
                  IF ( full ) THEN
                    WRITE (Nout,FMT=99006) nc, Sname, uplo, n, ralpha, incx, lda
                  ELSEIF ( packed ) THEN
                    WRITE (Nout,FMT=99005) nc, Sname, uplo, n, ralpha, incx
                  ENDIF
                ENDIF
              ENDIF
              !
            ENDDO
            !
          ENDDO
          !
        ENDDO
      ENDIF
      !
    ENDDO
    !
    !     Report result.
    !
    IF ( .NOT.Fatal ) THEN
      IF ( Kprint>=3 ) THEN
        IF ( errmax<Thresh ) THEN
          WRITE (Nout,FMT=99001) Sname, nc
        ELSE
          WRITE (Nout,FMT=99003) Sname, nc, errmax
        ENDIF
      ENDIF
    ENDIF
    RETURN
    !
    99001 FORMAT (' ',A6,' PASSED THE COMPUTATIONAL TESTS (',I6,' CALL','S)')
    99002 FORMAT (' ******* FATAL ERROR - PARAMETER NUMBER ',I2,' WAS CH',&
      'ANGED INCORRECTLY *******')
    99003 FORMAT (' ',A6,' COMPLETED THE COMPUTATIONAL TESTS (',I6,' C','ALLS)',&
      /' ******* BUT WITH MAXIMUM TEST RATIO',F8.2,' - SUSPECT *******')
    99004 FORMAT (' ******* ',A6,' FAILED ON CALL NUMBER:')
    99005 FORMAT (1X,I6,': ',A6,'(''',A1,''',',I3,',',F4.1,', X,',I2,&
      ', AP)                                         .')
    99006 FORMAT (1X,I6,': ',A6,'(''',A1,''',',I3,',',F4.1,', X,',I2,', A,',I3,&
      ')                                      .')
    99007 FORMAT (' ******* FATAL ERROR - ERROR-EXIT TAKEN ON VALID CALL *',&
      '******')
    !
    !     End of CCHK52.
    !
  END SUBROUTINE CCHK52
  !** CCHK53
  SUBROUTINE CCHK53(Sname,Eps,Thresh,Nout,Kprint,Fatal,Nidim,Idim,Nalf,Alf,&
      Nbet,Bet,Nmax,Ab,Aa,As,Bb,Bs,C,Cc,Cs,Ct,G,W)
    IMPLICIT NONE
    !>
    !***
    !  Quick check for CHER2K and CSYR2K.
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
    !  Quick check for CHER2K and CSYR2K.
    !
    !  Auxiliary routine for test program for Level 3 Blas.
    !***
    ! **References:**  (NONE)
    !***
    ! **Routines called:**  CHER2K, CMAKE3, CMMCH, CSYR2K, LCE, LCERES, NUMXER

    !* REVISION HISTORY  (YYMMDD)
    !   890208  DATE WRITTEN
    !   910619  Modified to meet SLATEC code and prologue standards.  (BKS)

    !     .. Parameters ..
    COMPLEX, PARAMETER :: ZERO = (0.0,0.0), ONE = (1.0,0.0)
    REAL, PARAMETER :: RZERO = 0.0, RONE = 1.0
    !     .. Scalar Arguments ..
    LOGICAL Fatal
    REAL Eps, Thresh
    INTEGER Kprint, Nalf, Nbet, Nidim, Nmax, Nout
    CHARACTER(6) :: Sname
    !     .. Array Arguments ..
    COMPLEX Aa(Nmax*Nmax), Ab(2*Nmax*Nmax), Alf(Nalf), As(Nmax*Nmax), &
      Bb(Nmax*Nmax), Bet(Nbet), Bs(Nmax*Nmax), C(Nmax,Nmax), &
      Cc(Nmax*Nmax), Cs(Nmax*Nmax), Ct(Nmax), W(2*Nmax)
    REAL G(Nmax)
    INTEGER Idim(Nidim)
    !     .. Local Scalars ..
    COMPLEX alpha, als, beta, bets
    REAL err, errmax, rbeta, rbets
    INTEGER i, ia, ib, ict, icu, ik, in, j, jc, jj, jjab, k, ks, &
      laa, lbb, lcc, lda, ldas, ldb, ldbs, ldc, ldcs, lj, ma, &
      n, na, nargs, nc, nerr, ns
    LOGICAL conj, ftl, null, reset, tran, upper
    CHARACTER :: trans, transs, uplo, transt, uplos
    !     .. Local Arrays ..
    LOGICAL isame(13)
    !     .. External Functions ..
    INTEGER, EXTERNAL :: NUMXER
    !     .. External Subroutines ..
    EXTERNAL :: CHER2K, CSYR2K
    !     .. Intrinsic Functions ..
    INTRINSIC ABS, MAX, MIN
    !     .. Data statements ..
    CHARACTER(2), PARAMETER :: ichu = 'UL', icht = 'NC'
    !* FIRST EXECUTABLE STATEMENT  CCHK53
    conj = Sname(2:3)=='HE'
    !
    nargs = 12
    nc = 0
    reset = .TRUE.
    errmax = RZERO
    !
    DO in = 1, Nidim
      n = Idim(in)
      !        Set LDC to 1 more than minimum value if room.
      ldc = n
      IF ( ldc<Nmax ) ldc = ldc + 1
      !        Skip tests if not enough room.
      IF ( ldc<=Nmax ) THEN
        lcc = ldc*n
        !
        DO ik = 1, Nidim
          k = Idim(ik)
          !
          DO ict = 1, 2
            trans = icht(ict:ict)
            tran = trans=='C'
            IF ( tran.AND..NOT.conj ) trans = 'T'
            IF ( tran ) THEN
              ma = k
              na = n
            ELSE
              ma = n
              na = k
            ENDIF
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
                CALL CMAKE3('GE',' ',' ',ma,na,Ab,2*Nmax,Aa,lda,reset,ZERO)
              ELSE
                CALL CMAKE3('GE',' ',' ',ma,na,Ab,Nmax,Aa,lda,reset,ZERO)
              ENDIF
              !
              !              Generate the matrix B.
              !
              ldb = lda
              lbb = laa
              IF ( tran ) THEN
                CALL CMAKE3('GE',' ',' ',ma,na,Ab(k+1),2*Nmax,Bb,ldb,reset,ZERO)
              ELSE
                CALL CMAKE3('GE',' ',' ',ma,na,Ab(k*Nmax+1),Nmax,Bb,ldb,reset,&
                  ZERO)
              ENDIF
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
                    IF ( conj ) THEN
                      rbeta = REAL(beta)
                      beta = CMPLX(rbeta,RZERO)
                    ENDIF
                    null = n<=0
                    IF ( conj ) null = null .OR.&
                      ((k<=0.OR.alpha==ZERO).AND.rbeta==RONE)
                    !
                    !                       Generate the matrix C.
                    !
                    CALL CMAKE3(Sname(2:3),uplo,' ',n,n,C,Nmax,Cc,ldc,reset,ZERO)
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
                    ENDDO
                    ldas = lda
                    DO i = 1, lbb
                      Bs(i) = Bb(i)
                    ENDDO
                    ldbs = ldb
                    IF ( conj ) THEN
                      rbets = rbeta
                    ELSE
                      bets = beta
                    ENDIF
                    DO i = 1, lcc
                      Cs(i) = Cc(i)
                    ENDDO
                    ldcs = ldc
                    !
                    !                       Call the subroutine.
                    !
                    IF ( conj ) THEN
                      CALL CHER2K(uplo,trans,n,k,alpha,Aa,lda,Bb,ldb,rbeta,Cc,ldc)
                    ELSE
                      CALL CSYR2K(uplo,trans,n,k,alpha,Aa,lda,Bb,ldb,beta,Cc,ldc)
                    ENDIF
                    !
                    !                       Check if error-exit was taken incorrectly.
                    !
                    IF ( NUMXER(nerr)/=0 ) THEN
                      IF ( Kprint>=2 ) WRITE (Nout,FMT=99007)
                      Fatal = .TRUE.
                    ENDIF
                    !
                    !                       See what data changed inside subroutines.
                    !
                    isame(1) = uplos==uplo
                    isame(2) = transs==trans
                    isame(3) = ns==n
                    isame(4) = ks==k
                    isame(5) = als==alpha
                    isame(6) = LCE(As,Aa,laa)
                    isame(7) = ldas==lda
                    isame(8) = LCE(Bs,Bb,lbb)
                    isame(9) = ldbs==ldb
                    IF ( conj ) THEN
                      isame(10) = rbets==rbeta
                    ELSE
                      isame(10) = bets==beta
                    ENDIF
                    IF ( null ) THEN
                      isame(11) = LCE(Cs,Cc,lcc)
                    ELSE
                      isame(11) = LCERES('HE',uplo,n,n,Cs,Cc,ldc)
                    ENDIF
                    isame(12) = ldcs==ldc
                    !
                    !                       If data was incorrectly changed, report and
                    !                       return.
                    !
                    DO i = 1, nargs
                      IF ( .NOT.isame(i) ) THEN
                        Fatal = .TRUE.
                        IF ( Kprint>=2 ) WRITE (Nout,FMT=99002) i
                      ENDIF
                    ENDDO
                    !
                    IF ( .NOT.null ) THEN
                      !
                      !                          Check the result column by column.
                      !
                      IF ( conj ) THEN
                        transt = 'C'
                      ELSE
                        transt = 'T'
                      ENDIF
                      jjab = 1
                      jc = 1
                      DO j = 1, n
                        IF ( upper ) THEN
                          jj = 1
                          lj = j
                        ELSE
                          jj = j
                          lj = n - j + 1
                        ENDIF
                        IF ( tran ) THEN
                          DO i = 1, k
                            W(i) = alpha*Ab((j-1)*2*Nmax+k+i)
                            IF ( conj ) THEN
                              W(k+i) = CONJG(alpha)*Ab((j-1)*2*Nmax+i)
                            ELSE
                              W(k+i) = alpha*Ab((j-1)*2*Nmax+i)
                            ENDIF
                          ENDDO
                          ftl = .FALSE.
                          CALL CMMCH(transt,'N',lj,1,2*k,ONE,Ab(jjab),2*Nmax,&
                            W,2*Nmax,beta,C(jj,j),Nmax,Ct,G,Cc(jc),&
                            ldc,Eps,err,ftl,Nout,.TRUE.,Kprint)
                        ELSE
                          DO i = 1, k
                            IF ( conj ) THEN
                              W(i) = alpha*CONJG(Ab((k+i-1)*Nmax+j))
                              W(k+i) = CONJG(alpha*Ab((i-1)*Nmax+j))
                            ELSE
                              W(i) = alpha*Ab((k+i-1)*Nmax+j)
                              W(k+i) = alpha*Ab((i-1)*Nmax+j)
                            ENDIF
                          ENDDO
                          ftl = .FALSE.
                          CALL CMMCH('N','N',lj,1,2*k,ONE,Ab(jj),Nmax,W,&
                            2*Nmax,beta,C(jj,j),Nmax,Ct,G,Cc(jc),ldc,&
                            Eps,err,ftl,Nout,.TRUE.,Kprint)
                        ENDIF
                        IF ( upper ) THEN
                          jc = jc + ldc
                        ELSE
                          jc = jc + ldc + 1
                          IF ( tran ) jjab = jjab + 2*Nmax
                        ENDIF
                        errmax = MAX(errmax,err)
                        IF ( ftl ) THEN
                          Fatal = .TRUE.
                          IF ( Kprint>=3 ) THEN
                            WRITE (Nout,FMT=99004) Sname
                            IF ( conj ) THEN
                              WRITE (Nout,FMT=99005) nc, Sname, uplo, &
                                trans, n, k, alpha, lda, ldb, rbeta, ldc
                            ELSE
                              WRITE (Nout,FMT=99006) nc, Sname, uplo, &
                                trans, n, k, alpha, lda, ldb, beta, ldc
                            ENDIF
                          ENDIF
                        ENDIF
                      ENDDO
                    ENDIF
                    !
                  ENDDO
                  !
                ENDDO
                !
              ENDDO
            ENDIF
            !
          ENDDO
          !
        ENDDO
      ENDIF
      !
    ENDDO
    !
    !     Report result.
    !
    IF ( .NOT.Fatal ) THEN
      IF ( Kprint>=3 ) THEN
        IF ( errmax<Thresh ) THEN
          WRITE (Nout,FMT=99001) Sname, nc
        ELSE
          WRITE (Nout,FMT=99003) Sname, nc, errmax
        ENDIF
      ENDIF
    ENDIF
    RETURN
    !
    99001 FORMAT (' ',A6,' PASSED THE COMPUTATIONAL TESTS (',I6,' CALL','S)')
    99002 FORMAT (' ******* FATAL ERROR - PARAMETER NUMBER ',I2,' WAS CH',&
      'ANGED INCORRECTLY *******')
    99003 FORMAT (' ',A6,' COMPLETED THE COMPUTATIONAL TESTS (',I6,' C','ALLS)',&
      /' ******* BUT WITH MAXIMUM TEST RATIO',F8.2,' - SUSPECT *******')
    99004 FORMAT (' ******* ',A6,' FAILED ON CALL NUMBER:')
    99005 FORMAT (1X,I6,': ',A6,'(',2('''',A1,''','),2(I3,','),'(',F4.1,',',F4.1,&
      '), A,',I3,', B,',I3,',',F4.1,', C,',I3,')           .')
    99006 FORMAT (1X,I6,': ',A6,'(',2('''',A1,''','),2(I3,','),'(',F4.1,',',F4.1,&
      '), A,',I3,', B,',I3,',(',F4.1,',',F4.1,'), C,',I3,')    .')
    99007 FORMAT (' ******* FATAL ERROR - ERROR-EXIT TAKEN ON VALID CALL *',&
      '******')
    !
    !     End of CCHK53.
    !
  END SUBROUTINE CCHK53
  !** CCHK62
  SUBROUTINE CCHK62(Sname,Eps,Thresh,Nout,Kprint,Fatal,Nidim,Idim,Nalf,Alf,&
      Ninc,Inc,Nmax,Incmax,A,Aa,As,X,Xx,Xs,Y,Yy,Ys,Yt,G,Z)
    IMPLICIT NONE
    !>
    !***
    !  Quick check for CHER2 and CHPR2.
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
    !  Quick check for CHER2 and CHPR2.
    !
    !  Auxiliary routine for test program for Level 2 Blas.
    !***
    ! **References:**  (NONE)
    !***
    ! **Routines called:**  CHER2, CHPR2, CMAKE2, CMVCH, LCE, LCERES, NUMXER

    !* REVISION HISTORY  (YYMMDD)
    !   870810  DATE WRITTEN
    !   910619  Modified to meet SLATEC code and prologue standards.  (BKS)

    !     .. Parameters ..
    COMPLEX, PARAMETER :: ZERO = (0.0,0.0), HALF = (0.5,0.0), ONE = (1.0,0.0)
    REAL, PARAMETER :: RZERO = 0.0
    !     .. Scalar Arguments ..
    LOGICAL Fatal
    REAL Eps, Thresh
    INTEGER Incmax, Kprint, Nalf, Nidim, Ninc, Nmax, Nout
    CHARACTER(6) :: Sname
    !     .. Array Arguments ..
    COMPLEX A(Nmax,Nmax), Aa(Nmax*Nmax), Alf(Nalf), As(Nmax*Nmax), X(Nmax), &
      Xs(Nmax*Incmax), Xx(Nmax*Incmax), Y(Nmax), Ys(Nmax*Incmax), &
      Yt(Nmax), Yy(Nmax*Incmax), Z(Nmax,2)
    REAL G(Nmax)
    INTEGER Idim(Nidim), Inc(Ninc)
    !     .. Local Scalars ..
    COMPLEX alpha, als, transl
    REAL err, errmax
    INTEGER i, ia, ic, in, incx, incxs, incy, incys, ix, iy, j, &
      ja, jj, laa, lda, ldas, lj, lx, ly, n, nargs, nc, nerr, ns
    LOGICAL ftl, full, null, packed, reset, upper
    CHARACTER :: uplo, uplos
    !     .. Local Arrays ..
    COMPLEX w(2)
    LOGICAL isame(13)
    !     .. External Functions ..
    INTEGER, EXTERNAL :: NUMXER
    !     .. External Subroutines ..
    EXTERNAL :: CHER2, CHPR2
    !     .. Intrinsic Functions ..
    INTRINSIC ABS, CONJG, MAX
    !     .. Data statements ..
    CHARACTER(2), PARAMETER :: ich = 'UL'
    !* FIRST EXECUTABLE STATEMENT  CCHK62
    full = Sname(3:3)=='E'
    packed = Sname(3:3)=='P'
    !     Define the number of arguments.
    IF ( full ) THEN
      nargs = 9
    ELSEIF ( packed ) THEN
      nargs = 8
    ENDIF
    !
    nc = 0
    reset = .TRUE.
    errmax = RZERO
    !
    DO in = 1, Nidim
      n = Idim(in)
      !        Set LDA to 1 more than minimum value if room.
      lda = n
      IF ( lda<Nmax ) lda = lda + 1
      !        Skip tests if not enough room.
      IF ( lda<=Nmax ) THEN
        IF ( packed ) THEN
          laa = (n*(n+1))/2
        ELSE
          laa = lda*n
        ENDIF
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
            CALL CMAKE2('GE',' ',' ',1,n,X,1,Xx,ABS(incx),0,n-1,reset,transl)
            IF ( n>1 ) THEN
              X(n/2) = ZERO
              Xx(1+ABS(incx)*(n/2-1)) = ZERO
            ENDIF
            !
            DO iy = 1, Ninc
              incy = Inc(iy)
              ly = ABS(incy)*n
              !
              !                 Generate the vector Y.
              !
              transl = ZERO
              CALL CMAKE2('GE',' ',' ',1,n,Y,1,Yy,ABS(incy),0,n-1,reset,transl)
              IF ( n>1 ) THEN
                Y(n/2) = ZERO
                Yy(1+ABS(incy)*(n/2-1)) = ZERO
              ENDIF
              !
              DO ia = 1, Nalf
                alpha = Alf(ia)
                null = n<=0 .OR. alpha==ZERO
                !
                !                    Generate the matrix A.
                !
                transl = ZERO
                CALL CMAKE2(Sname(2:3),uplo,' ',n,n,A,Nmax,Aa,lda,n-1,n-1,&
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
                ENDDO
                ldas = lda
                DO i = 1, lx
                  Xs(i) = Xx(i)
                ENDDO
                incxs = incx
                DO i = 1, ly
                  Ys(i) = Yy(i)
                ENDDO
                incys = incy
                !
                !                    Call the subroutine.
                !
                IF ( full ) THEN
                  CALL CHER2(uplo,n,alpha,Xx,incx,Yy,incy,Aa,lda)
                ELSEIF ( packed ) THEN
                  CALL CHPR2(uplo,n,alpha,Xx,incx,Yy,incy,Aa)
                ENDIF
                !
                !                    Check if error-exit was taken incorrectly.
                !
                IF ( NUMXER(nerr)/=0 ) THEN
                  IF ( Kprint>=2 ) WRITE (Nout,FMT=99008)
                  Fatal = .TRUE.
                ENDIF
                !
                !                    See what data changed inside subroutines.
                !
                isame(1) = uplo==uplos
                isame(2) = ns==n
                isame(3) = als==alpha
                isame(4) = LCE(Xs,Xx,lx)
                isame(5) = incxs==incx
                isame(6) = LCE(Ys,Yy,ly)
                isame(7) = incys==incy
                IF ( null ) THEN
                  isame(8) = LCE(As,Aa,laa)
                ELSE
                  isame(8) = LCERES(Sname(2:3),uplo,n,n,As,Aa,lda)
                ENDIF
                IF ( .NOT.packed ) isame(9) = ldas==lda
                !
                !                    If data was incorrectly changed, report and return.
                !
                DO i = 1, nargs
                  IF ( .NOT.isame(i) ) THEN
                    Fatal = .TRUE.
                    IF ( Kprint>=2 ) WRITE (Nout,FMT=99002) i
                  ENDIF
                ENDDO
                !
                ftl = .FALSE.
                IF ( .NOT.null ) THEN
                  !
                  !                       Check the result column by column.
                  !
                  IF ( incx>0 ) THEN
                    DO i = 1, n
                      Z(i,1) = X(i)
                    ENDDO
                  ELSE
                    DO i = 1, n
                      Z(i,1) = X(n-i+1)
                    ENDDO
                  ENDIF
                  IF ( incy>0 ) THEN
                    DO i = 1, n
                      Z(i,2) = Y(i)
                    ENDDO
                  ELSE
                    DO i = 1, n
                      Z(i,2) = Y(n-i+1)
                    ENDDO
                  ENDIF
                  ja = 1
                  DO j = 1, n
                    w(1) = alpha*CONJG(Z(j,2))
                    w(2) = CONJG(alpha)*CONJG(Z(j,1))
                    IF ( upper ) THEN
                      jj = 1
                      lj = j
                    ELSE
                      jj = j
                      lj = n - j + 1
                    ENDIF
                    CALL CMVCH('N',lj,2,ONE,Z(jj,1),Nmax,w,1,ONE,A(jj,j),1,Yt,&
                      G,Aa(ja),Eps,err,ftl,Nout,.TRUE.,Kprint)
                    IF ( .NOT.(full) ) THEN
                      ja = ja + lj
                    ELSEIF ( upper ) THEN
                      ja = ja + lda
                    ELSE
                      ja = ja + lda + 1
                    ENDIF
                    errmax = MAX(errmax,err)
                  ENDDO
                ENDIF
                IF ( ftl ) THEN
                  Fatal = .TRUE.
                  IF ( Kprint>=3 ) THEN
                    WRITE (Nout,FMT=99005) j
                    WRITE (Nout,FMT=99004) Sname
                    IF ( full ) THEN
                      WRITE (Nout,FMT=99007) nc, Sname, uplo, n, alpha, &
                        incx, incy, lda
                    ELSEIF ( packed ) THEN
                      WRITE (Nout,FMT=99006) nc, Sname, uplo, n, alpha, incx, incy
                    ENDIF
                  ENDIF
                ENDIF
                !
              ENDDO
              !
            ENDDO
            !
          ENDDO
          !
        ENDDO
      ENDIF
      !
    ENDDO
    !
    !     Report result.
    !
    IF ( .NOT.Fatal ) THEN
      IF ( Kprint>=3 ) THEN
        IF ( errmax<Thresh ) THEN
          WRITE (Nout,FMT=99001) Sname, nc
        ELSE
          WRITE (Nout,FMT=99003) Sname, nc, errmax
        ENDIF
      ENDIF
    ENDIF
    RETURN
    !
    99001 FORMAT (' ',A6,' PASSED THE COMPUTATIONAL TESTS (',I6,' CALL','S)')
    99002 FORMAT (' ******* FATAL ERROR - PARAMETER NUMBER ',I2,' WAS CH',&
      'ANGED INCORRECTLY *******')
    99003 FORMAT (' ',A6,' COMPLETED THE COMPUTATIONAL TESTS (',I6,' C','ALLS)',&
      /' ******* BUT WITH MAXIMUM TEST RATIO',F8.2,' - SUSPECT *******')
    99004 FORMAT (' ******* ',A6,' FAILED ON CALL NUMBER:')
    99005 FORMAT ('      THESE ARE THE RESULTS FOR COLUMN ',I3)
    99006 FORMAT (1X,I6,': ',A6,'(''',A1,''',',I3,',(',F4.1,',',F4.1,'), X,',I2,&
      ', Y,',I2,', AP)                     ','       .')
    99007 FORMAT (1X,I6,': ',A6,'(''',A1,''',',I3,',(',F4.1,',',F4.1,'), X,',I2,&
      ', Y,',I2,', A,',I3,')             ','            .')
    99008 FORMAT (' ******* FATAL ERROR - ERROR-EXIT TAKEN ON VALID CALL *',&
      '******')
    !
    !     End of CCHK62.
    !
  END SUBROUTINE CCHK62
  !** CCHKE2
  SUBROUTINE CCHKE2(Isnum,Srnamt,Nout,Kprint,Fatal)
    IMPLICIT NONE
    !>
    !***
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
    ! **Routines called:**  CGBMV, CGEMV, CGERC, CGERU, CHBMV, CHEMV, CHER,
    !                    CHER2, CHKXER, CHPMV, CHPR, CHPR2, CTBMV, CTBSV,
    !                    CTPMV, CTPSV, CTRMV, CTRSV, XERCLR, XERDMP, XGETF,
    !                    XSETF

    !* REVISION HISTORY  (YYMMDD)
    !   870810  DATE WRITTEN
    !   910620  Modified to meet SLATEC code and prologue standards.  (BKS)

    !     .. Scalar Arguments ..
    LOGICAL Fatal
    INTEGER Isnum, Kprint, Nout
    CHARACTER(6) :: Srnamt
    !     .. Scalars in Common ..
    INTEGER infot
    !     .. Local Scalars ..
    COMPLEX alpha, beta
    REAL ralpha
    INTEGER kontrl
    !     .. Local Arrays ..
    COMPLEX a(1,1), x(1), y(1)
    !     .. External Subroutines ..
    EXTERNAL :: CGBMV, CGEMV, CGERC, CGERU, CHBMV, CHEMV, CHER, CHER2, &
      CHKXER, CHPMV, CHPR, CHPR2, CTBMV, CTBSV, CTPMV, CTPSV, CTRMV, CTRSV
    !* FIRST EXECUTABLE STATEMENT  CCHKE2
    CALL XGETF(kontrl)
    IF ( Kprint<=2 ) THEN
      CALL XSETF(0)
    ELSE
      CALL XSETF(1)
    ENDIF
    SELECT CASE (Isnum)
      CASE (2)
        infot = 1
        CALL XERCLR
        CALL CGBMV('/',0,0,0,0,alpha,a,1,x,1,beta,y,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 2
        CALL XERCLR
        CALL CGBMV('N',-1,0,0,0,alpha,a,1,x,1,beta,y,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 3
        CALL XERCLR
        CALL CGBMV('N',0,-1,0,0,alpha,a,1,x,1,beta,y,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 4
        CALL XERCLR
        CALL CGBMV('N',0,0,-1,0,alpha,a,1,x,1,beta,y,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 5
        CALL XERCLR
        CALL CGBMV('N',2,0,0,-1,alpha,a,1,x,1,beta,y,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 8
        CALL XERCLR
        CALL CGBMV('N',0,0,1,0,alpha,a,1,x,1,beta,y,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 10
        CALL XERCLR
        CALL CGBMV('N',0,0,0,0,alpha,a,1,x,0,beta,y,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 13
        CALL XERCLR
        CALL CGBMV('N',0,0,0,0,alpha,a,1,x,1,beta,y,0)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      CASE (3)
        infot = 1
        CALL XERCLR
        CALL CHEMV('/',0,alpha,a,1,x,1,beta,y,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 2
        CALL XERCLR
        CALL CHEMV('U',-1,alpha,a,1,x,1,beta,y,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 5
        CALL XERCLR
        CALL CHEMV('U',2,alpha,a,1,x,1,beta,y,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 7
        CALL XERCLR
        CALL CHEMV('U',0,alpha,a,1,x,0,beta,y,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 10
        CALL XERCLR
        CALL CHEMV('U',0,alpha,a,1,x,1,beta,y,0)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      CASE (4)
        infot = 1
        CALL XERCLR
        CALL CHBMV('/',0,0,alpha,a,1,x,1,beta,y,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 2
        CALL XERCLR
        CALL CHBMV('U',-1,0,alpha,a,1,x,1,beta,y,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 3
        CALL XERCLR
        CALL CHBMV('U',0,-1,alpha,a,1,x,1,beta,y,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 6
        CALL XERCLR
        CALL CHBMV('U',0,1,alpha,a,1,x,1,beta,y,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 8
        CALL XERCLR
        CALL CHBMV('U',0,0,alpha,a,1,x,0,beta,y,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 11
        CALL XERCLR
        CALL CHBMV('U',0,0,alpha,a,1,x,1,beta,y,0)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      CASE (5)
        infot = 1
        CALL XERCLR
        CALL CHPMV('/',0,alpha,a,x,1,beta,y,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 2
        CALL XERCLR
        CALL CHPMV('U',-1,alpha,a,x,1,beta,y,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 6
        CALL XERCLR
        CALL CHPMV('U',0,alpha,a,x,0,beta,y,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 9
        CALL XERCLR
        CALL CHPMV('U',0,alpha,a,x,1,beta,y,0)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      CASE (6)
        infot = 1
        CALL XERCLR
        CALL CTRMV('/','N','N',0,a,1,x,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 2
        CALL XERCLR
        CALL CTRMV('U','/','N',0,a,1,x,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 3
        CALL XERCLR
        CALL CTRMV('U','N','/',0,a,1,x,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 4
        CALL XERCLR
        CALL CTRMV('U','N','N',-1,a,1,x,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 6
        CALL XERCLR
        CALL CTRMV('U','N','N',2,a,1,x,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 8
        CALL XERCLR
        CALL CTRMV('U','N','N',0,a,1,x,0)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      CASE (7)
        infot = 1
        CALL XERCLR
        CALL CTBMV('/','N','N',0,0,a,1,x,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 2
        CALL XERCLR
        CALL CTBMV('U','/','N',0,0,a,1,x,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 3
        CALL XERCLR
        CALL CTBMV('U','N','/',0,0,a,1,x,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 4
        CALL XERCLR
        CALL CTBMV('U','N','N',-1,0,a,1,x,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 5
        CALL XERCLR
        CALL CTBMV('U','N','N',0,-1,a,1,x,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 7
        CALL XERCLR
        CALL CTBMV('U','N','N',0,1,a,1,x,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 9
        CALL XERCLR
        CALL CTBMV('U','N','N',0,0,a,1,x,0)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      CASE (8)
        infot = 1
        CALL XERCLR
        CALL CTPMV('/','N','N',0,a,x,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 2
        CALL XERCLR
        CALL CTPMV('U','/','N',0,a,x,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 3
        CALL XERCLR
        CALL CTPMV('U','N','/',0,a,x,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 4
        CALL XERCLR
        CALL CTPMV('U','N','N',-1,a,x,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 7
        CALL XERCLR
        CALL CTPMV('U','N','N',0,a,x,0)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      CASE (9)
        infot = 1
        CALL XERCLR
        CALL CTRSV('/','N','N',0,a,1,x,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 2
        CALL XERCLR
        CALL CTRSV('U','/','N',0,a,1,x,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 3
        CALL XERCLR
        CALL CTRSV('U','N','/',0,a,1,x,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 4
        CALL XERCLR
        CALL CTRSV('U','N','N',-1,a,1,x,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 6
        CALL XERCLR
        CALL CTRSV('U','N','N',2,a,1,x,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 8
        CALL XERCLR
        CALL CTRSV('U','N','N',0,a,1,x,0)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      CASE (10)
        infot = 1
        CALL XERCLR
        CALL CTBSV('/','N','N',0,0,a,1,x,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 2
        CALL XERCLR
        CALL CTBSV('U','/','N',0,0,a,1,x,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 3
        CALL XERCLR
        CALL CTBSV('U','N','/',0,0,a,1,x,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 4
        CALL XERCLR
        CALL CTBSV('U','N','N',-1,0,a,1,x,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 5
        CALL XERCLR
        CALL CTBSV('U','N','N',0,-1,a,1,x,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 7
        CALL XERCLR
        CALL CTBSV('U','N','N',0,1,a,1,x,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 9
        CALL XERCLR
        CALL CTBSV('U','N','N',0,0,a,1,x,0)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      CASE (11)
        infot = 1
        CALL XERCLR
        CALL CTPSV('/','N','N',0,a,x,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 2
        CALL XERCLR
        CALL CTPSV('U','/','N',0,a,x,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 3
        CALL XERCLR
        CALL CTPSV('U','N','/',0,a,x,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 4
        CALL XERCLR
        CALL CTPSV('U','N','N',-1,a,x,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 7
        CALL XERCLR
        CALL CTPSV('U','N','N',0,a,x,0)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      CASE (12)
        infot = 1
        CALL XERCLR
        CALL CGERC(-1,0,alpha,x,1,y,1,a,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 2
        CALL XERCLR
        CALL CGERC(0,-1,alpha,x,1,y,1,a,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 5
        CALL XERCLR
        CALL CGERC(0,0,alpha,x,0,y,1,a,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 7
        CALL XERCLR
        CALL CGERC(0,0,alpha,x,1,y,0,a,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 9
        CALL XERCLR
        CALL CGERC(2,0,alpha,x,1,y,1,a,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      CASE (13)
        infot = 1
        CALL XERCLR
        CALL CGERU(-1,0,alpha,x,1,y,1,a,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 2
        CALL XERCLR
        CALL CGERU(0,-1,alpha,x,1,y,1,a,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 5
        CALL XERCLR
        CALL CGERU(0,0,alpha,x,0,y,1,a,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 7
        CALL XERCLR
        CALL CGERU(0,0,alpha,x,1,y,0,a,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 9
        CALL XERCLR
        CALL CGERU(2,0,alpha,x,1,y,1,a,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      CASE (14)
        infot = 1
        CALL XERCLR
        CALL CHER('/',0,ralpha,x,1,a,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 2
        CALL XERCLR
        CALL CHER('U',-1,ralpha,x,1,a,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 5
        CALL XERCLR
        CALL CHER('U',0,ralpha,x,0,a,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 7
        CALL XERCLR
        CALL CHER('U',2,ralpha,x,1,a,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      CASE (15)
        infot = 1
        CALL XERCLR
        CALL CHPR('/',0,ralpha,x,1,a)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 2
        CALL XERCLR
        CALL CHPR('U',-1,ralpha,x,1,a)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 5
        CALL XERCLR
        CALL CHPR('U',0,ralpha,x,0,a)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      CASE (16)
        infot = 1
        CALL XERCLR
        CALL CHER2('/',0,alpha,x,1,y,1,a,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 2
        CALL XERCLR
        CALL CHER2('U',-1,alpha,x,1,y,1,a,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 5
        CALL XERCLR
        CALL CHER2('U',0,alpha,x,0,y,1,a,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 7
        CALL XERCLR
        CALL CHER2('U',0,alpha,x,1,y,0,a,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 9
        CALL XERCLR
        CALL CHER2('U',2,alpha,x,1,y,1,a,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      CASE (17)
        infot = 1
        CALL XERCLR
        CALL CHPR2('/',0,alpha,x,1,y,1,a)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 2
        CALL XERCLR
        CALL CHPR2('U',-1,alpha,x,1,y,1,a)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 5
        CALL XERCLR
        CALL CHPR2('U',0,alpha,x,0,y,1,a)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 7
        CALL XERCLR
        CALL CHPR2('U',0,alpha,x,1,y,0,a)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      CASE DEFAULT
        infot = 1
        CALL XERCLR
        CALL CGEMV('/',0,0,alpha,a,1,x,1,beta,y,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 2
        CALL XERCLR
        CALL CGEMV('N',-1,0,alpha,a,1,x,1,beta,y,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 3
        CALL XERCLR
        CALL CGEMV('N',0,-1,alpha,a,1,x,1,beta,y,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 6
        CALL XERCLR
        CALL CGEMV('N',2,0,alpha,a,1,x,1,beta,y,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 8
        CALL XERCLR
        CALL CGEMV('N',0,0,alpha,a,1,x,0,beta,y,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 11
        CALL XERCLR
        CALL CGEMV('N',0,0,alpha,a,1,x,1,beta,y,0)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
    END SELECT
    !
    IF ( Kprint>=2 ) THEN
      CALL XERDMP
      IF ( .NOT.Fatal ) THEN
        WRITE (Nout,FMT=99001) Srnamt
      ELSE
        WRITE (Nout,FMT=99002) Srnamt
      ENDIF
    ENDIF
    CALL XSETF(kontrl)
    RETURN
    !
    99001 FORMAT (' ',A6,' PASSED THE TESTS OF ERROR-EXITS')
    99002 FORMAT (' ******* ',A6,' FAILED THE TESTS OF ERROR-EXITS *****','**')
    !
    !     End of CCHKE2.
    !
  END SUBROUTINE CCHKE2
  !** CCHKE3
  SUBROUTINE CCHKE3(Isnum,Srnamt,Nout,Kprint,Fatal)
    IMPLICIT NONE
    !>
    !***
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
    ! **Routines called:**  CGEMM, CHEMM, CHER2K, CHERK, CHKXER, CSYMM, CSYR2K,
    !                    CSYRK, CTRMM, CTRSM, XERCLR, XERDMP, XGETF, XSETF

    !* REVISION HISTORY  (YYMMDD)
    !   890208  DATE WRITTEN
    !   910620  Modified to meet SLATEC code and prologue standards.  (BKS)

    !     .. Scalar Arguments ..
    LOGICAL Fatal
    INTEGER Isnum, Kprint, Nout
    CHARACTER(6) :: Srnamt
    !     .. Scalars in Common ..
    INTEGER infot
    !     .. Local Scalars ..
    COMPLEX alpha, beta
    REAL ralpha, rbeta
    INTEGER kontrl
    !     .. Local Arrays ..
    COMPLEX a(2,1), b(2,1), c(2,1)
    !     .. External Subroutines ..
    EXTERNAL :: CGEMM, CHEMM, CHER2K, CHERK, CHKXER, CSYMM, CSYR2K, &
      CSYRK, CTRMM, CTRSM
    !* FIRST EXECUTABLE STATEMENT  CCHKE3
    CALL XGETF(kontrl)
    IF ( Kprint<=2 ) THEN
      CALL XSETF(0)
    ELSE
      CALL XSETF(1)
    ENDIF
    SELECT CASE (Isnum)
      CASE (2)
        infot = 1
        CALL XERCLR
        CALL CHEMM('/','U',0,0,alpha,a,1,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 2
        CALL XERCLR
        CALL CHEMM('L','/',0,0,alpha,a,1,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 3
        CALL XERCLR
        CALL CHEMM('L','U',-1,0,alpha,a,1,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 3
        CALL XERCLR
        CALL CHEMM('R','U',-1,0,alpha,a,1,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 3
        CALL XERCLR
        CALL CHEMM('L','L',-1,0,alpha,a,1,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 3
        CALL XERCLR
        CALL CHEMM('R','L',-1,0,alpha,a,1,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 4
        CALL XERCLR
        CALL CHEMM('L','U',0,-1,alpha,a,1,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 4
        CALL XERCLR
        CALL CHEMM('R','U',0,-1,alpha,a,1,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 4
        CALL XERCLR
        CALL CHEMM('L','L',0,-1,alpha,a,1,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 4
        CALL XERCLR
        CALL CHEMM('R','L',0,-1,alpha,a,1,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 7
        CALL XERCLR
        CALL CHEMM('L','U',2,0,alpha,a,1,b,2,beta,c,2)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 7
        CALL XERCLR
        CALL CHEMM('R','U',0,2,alpha,a,1,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 7
        CALL XERCLR
        CALL CHEMM('L','L',2,0,alpha,a,1,b,2,beta,c,2)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 7
        CALL XERCLR
        CALL CHEMM('R','L',0,2,alpha,a,1,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 9
        CALL XERCLR
        CALL CHEMM('L','U',2,0,alpha,a,2,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 9
        CALL XERCLR
        CALL CHEMM('R','U',2,0,alpha,a,1,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 9
        CALL XERCLR
        CALL CHEMM('L','L',2,0,alpha,a,2,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 9
        CALL XERCLR
        CALL CHEMM('R','L',2,0,alpha,a,1,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 12
        CALL XERCLR
        CALL CHEMM('L','U',2,0,alpha,a,2,b,2,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 12
        CALL XERCLR
        CALL CHEMM('R','U',2,0,alpha,a,1,b,2,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 12
        CALL XERCLR
        CALL CHEMM('L','L',2,0,alpha,a,2,b,2,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 12
        CALL XERCLR
        CALL CHEMM('R','L',2,0,alpha,a,1,b,2,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      CASE (3)
        infot = 1
        CALL XERCLR
        CALL CSYMM('/','U',0,0,alpha,a,1,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 2
        CALL XERCLR
        CALL CSYMM('L','/',0,0,alpha,a,1,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 3
        CALL XERCLR
        CALL CSYMM('L','U',-1,0,alpha,a,1,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 3
        CALL XERCLR
        CALL CSYMM('R','U',-1,0,alpha,a,1,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 3
        CALL XERCLR
        CALL CSYMM('L','L',-1,0,alpha,a,1,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 3
        CALL XERCLR
        CALL CSYMM('R','L',-1,0,alpha,a,1,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 4
        CALL XERCLR
        CALL CSYMM('L','U',0,-1,alpha,a,1,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 4
        CALL XERCLR
        CALL CSYMM('R','U',0,-1,alpha,a,1,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 4
        CALL XERCLR
        CALL CSYMM('L','L',0,-1,alpha,a,1,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 4
        CALL XERCLR
        CALL CSYMM('R','L',0,-1,alpha,a,1,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 7
        CALL XERCLR
        CALL CSYMM('L','U',2,0,alpha,a,1,b,2,beta,c,2)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 7
        CALL XERCLR
        CALL CSYMM('R','U',0,2,alpha,a,1,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 7
        CALL XERCLR
        CALL CSYMM('L','L',2,0,alpha,a,1,b,2,beta,c,2)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 7
        CALL XERCLR
        CALL CSYMM('R','L',0,2,alpha,a,1,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 9
        CALL XERCLR
        CALL CSYMM('L','U',2,0,alpha,a,2,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 9
        CALL XERCLR
        CALL CSYMM('R','U',2,0,alpha,a,1,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 9
        CALL XERCLR
        CALL CSYMM('L','L',2,0,alpha,a,2,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 9
        CALL XERCLR
        CALL CSYMM('R','L',2,0,alpha,a,1,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 12
        CALL XERCLR
        CALL CSYMM('L','U',2,0,alpha,a,2,b,2,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 12
        CALL XERCLR
        CALL CSYMM('R','U',2,0,alpha,a,1,b,2,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 12
        CALL XERCLR
        CALL CSYMM('L','L',2,0,alpha,a,2,b,2,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 12
        CALL XERCLR
        CALL CSYMM('R','L',2,0,alpha,a,1,b,2,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      CASE (4)
        infot = 1
        CALL XERCLR
        CALL CTRMM('/','U','N','N',0,0,alpha,a,1,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 2
        CALL XERCLR
        CALL CTRMM('L','/','N','N',0,0,alpha,a,1,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 3
        CALL XERCLR
        CALL CTRMM('L','U','/','N',0,0,alpha,a,1,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 4
        CALL XERCLR
        CALL CTRMM('L','U','N','/',0,0,alpha,a,1,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 5
        CALL XERCLR
        CALL CTRMM('L','U','N','N',-1,0,alpha,a,1,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 5
        CALL XERCLR
        CALL CTRMM('L','U','C','N',-1,0,alpha,a,1,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 5
        CALL XERCLR
        CALL CTRMM('L','U','T','N',-1,0,alpha,a,1,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 5
        CALL XERCLR
        CALL CTRMM('R','U','N','N',-1,0,alpha,a,1,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 5
        CALL XERCLR
        CALL CTRMM('R','U','C','N',-1,0,alpha,a,1,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 5
        CALL XERCLR
        CALL CTRMM('R','U','T','N',-1,0,alpha,a,1,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 5
        CALL XERCLR
        CALL CTRMM('L','L','N','N',-1,0,alpha,a,1,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 5
        CALL XERCLR
        CALL CTRMM('L','L','C','N',-1,0,alpha,a,1,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 5
        CALL XERCLR
        CALL CTRMM('L','L','T','N',-1,0,alpha,a,1,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 5
        CALL XERCLR
        CALL CTRMM('R','L','N','N',-1,0,alpha,a,1,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 5
        CALL XERCLR
        CALL CTRMM('R','L','C','N',-1,0,alpha,a,1,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 5
        CALL XERCLR
        CALL CTRMM('R','L','T','N',-1,0,alpha,a,1,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 6
        CALL XERCLR
        CALL CTRMM('L','U','N','N',0,-1,alpha,a,1,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 6
        CALL XERCLR
        CALL CTRMM('L','U','C','N',0,-1,alpha,a,1,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 6
        CALL XERCLR
        CALL CTRMM('L','U','T','N',0,-1,alpha,a,1,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 6
        CALL XERCLR
        CALL CTRMM('R','U','N','N',0,-1,alpha,a,1,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 6
        CALL XERCLR
        CALL CTRMM('R','U','C','N',0,-1,alpha,a,1,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 6
        CALL XERCLR
        CALL CTRMM('R','U','T','N',0,-1,alpha,a,1,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 6
        CALL XERCLR
        CALL CTRMM('L','L','N','N',0,-1,alpha,a,1,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 6
        CALL XERCLR
        CALL CTRMM('L','L','C','N',0,-1,alpha,a,1,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 6
        CALL XERCLR
        CALL CTRMM('L','L','T','N',0,-1,alpha,a,1,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 6
        CALL XERCLR
        CALL CTRMM('R','L','N','N',0,-1,alpha,a,1,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 6
        CALL XERCLR
        CALL CTRMM('R','L','C','N',0,-1,alpha,a,1,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 6
        CALL XERCLR
        CALL CTRMM('R','L','T','N',0,-1,alpha,a,1,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 9
        CALL XERCLR
        CALL CTRMM('L','U','N','N',2,0,alpha,a,1,b,2)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 9
        CALL XERCLR
        CALL CTRMM('L','U','C','N',2,0,alpha,a,1,b,2)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 9
        CALL XERCLR
        CALL CTRMM('L','U','T','N',2,0,alpha,a,1,b,2)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 9
        CALL XERCLR
        CALL CTRMM('R','U','N','N',0,2,alpha,a,1,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 9
        CALL XERCLR
        CALL CTRMM('R','U','C','N',0,2,alpha,a,1,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 9
        CALL XERCLR
        CALL CTRMM('R','U','T','N',0,2,alpha,a,1,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 9
        CALL XERCLR
        CALL CTRMM('L','L','N','N',2,0,alpha,a,1,b,2)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 9
        CALL XERCLR
        CALL CTRMM('L','L','C','N',2,0,alpha,a,1,b,2)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 9
        CALL XERCLR
        CALL CTRMM('L','L','T','N',2,0,alpha,a,1,b,2)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 9
        CALL XERCLR
        CALL CTRMM('R','L','N','N',0,2,alpha,a,1,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 9
        CALL XERCLR
        CALL CTRMM('R','L','C','N',0,2,alpha,a,1,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 9
        CALL XERCLR
        CALL CTRMM('R','L','T','N',0,2,alpha,a,1,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 11
        CALL XERCLR
        CALL CTRMM('L','U','N','N',2,0,alpha,a,2,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 11
        CALL XERCLR
        CALL CTRMM('L','U','C','N',2,0,alpha,a,2,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 11
        CALL XERCLR
        CALL CTRMM('L','U','T','N',2,0,alpha,a,2,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 11
        CALL XERCLR
        CALL CTRMM('R','U','N','N',2,0,alpha,a,1,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 11
        CALL XERCLR
        CALL CTRMM('R','U','C','N',2,0,alpha,a,1,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 11
        CALL XERCLR
        CALL CTRMM('R','U','T','N',2,0,alpha,a,1,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 11
        CALL XERCLR
        CALL CTRMM('L','L','N','N',2,0,alpha,a,2,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 11
        CALL XERCLR
        CALL CTRMM('L','L','C','N',2,0,alpha,a,2,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 11
        CALL XERCLR
        CALL CTRMM('L','L','T','N',2,0,alpha,a,2,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 11
        CALL XERCLR
        CALL CTRMM('R','L','N','N',2,0,alpha,a,1,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 11
        CALL XERCLR
        CALL CTRMM('R','L','C','N',2,0,alpha,a,1,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 11
        CALL XERCLR
        CALL CTRMM('R','L','T','N',2,0,alpha,a,1,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      CASE (5)
        infot = 1
        CALL XERCLR
        CALL CTRSM('/','U','N','N',0,0,alpha,a,1,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 2
        CALL XERCLR
        CALL CTRSM('L','/','N','N',0,0,alpha,a,1,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 3
        CALL XERCLR
        CALL CTRSM('L','U','/','N',0,0,alpha,a,1,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 4
        CALL XERCLR
        CALL CTRSM('L','U','N','/',0,0,alpha,a,1,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 5
        CALL XERCLR
        CALL CTRSM('L','U','N','N',-1,0,alpha,a,1,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 5
        CALL XERCLR
        CALL CTRSM('L','U','C','N',-1,0,alpha,a,1,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 5
        CALL XERCLR
        CALL CTRSM('L','U','T','N',-1,0,alpha,a,1,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 5
        CALL XERCLR
        CALL CTRSM('R','U','N','N',-1,0,alpha,a,1,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 5
        CALL XERCLR
        CALL CTRSM('R','U','C','N',-1,0,alpha,a,1,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 5
        CALL XERCLR
        CALL CTRSM('R','U','T','N',-1,0,alpha,a,1,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 5
        CALL XERCLR
        CALL CTRSM('L','L','N','N',-1,0,alpha,a,1,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 5
        CALL XERCLR
        CALL CTRSM('L','L','C','N',-1,0,alpha,a,1,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 5
        CALL XERCLR
        CALL CTRSM('L','L','T','N',-1,0,alpha,a,1,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 5
        CALL XERCLR
        CALL CTRSM('R','L','N','N',-1,0,alpha,a,1,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 5
        CALL XERCLR
        CALL CTRSM('R','L','C','N',-1,0,alpha,a,1,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 5
        CALL XERCLR
        CALL CTRSM('R','L','T','N',-1,0,alpha,a,1,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 6
        CALL XERCLR
        CALL CTRSM('L','U','N','N',0,-1,alpha,a,1,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 6
        CALL XERCLR
        CALL CTRSM('L','U','C','N',0,-1,alpha,a,1,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 6
        CALL XERCLR
        CALL CTRSM('L','U','T','N',0,-1,alpha,a,1,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 6
        CALL XERCLR
        CALL CTRSM('R','U','N','N',0,-1,alpha,a,1,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 6
        CALL XERCLR
        CALL CTRSM('R','U','C','N',0,-1,alpha,a,1,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 6
        CALL XERCLR
        CALL CTRSM('R','U','T','N',0,-1,alpha,a,1,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 6
        CALL XERCLR
        CALL CTRSM('L','L','N','N',0,-1,alpha,a,1,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 6
        CALL XERCLR
        CALL CTRSM('L','L','C','N',0,-1,alpha,a,1,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 6
        CALL XERCLR
        CALL CTRSM('L','L','T','N',0,-1,alpha,a,1,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 6
        CALL XERCLR
        CALL CTRSM('R','L','N','N',0,-1,alpha,a,1,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 6
        CALL XERCLR
        CALL CTRSM('R','L','C','N',0,-1,alpha,a,1,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 6
        CALL XERCLR
        CALL CTRSM('R','L','T','N',0,-1,alpha,a,1,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 9
        CALL XERCLR
        CALL CTRSM('L','U','N','N',2,0,alpha,a,1,b,2)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 9
        CALL XERCLR
        CALL CTRSM('L','U','C','N',2,0,alpha,a,1,b,2)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 9
        CALL XERCLR
        CALL CTRSM('L','U','T','N',2,0,alpha,a,1,b,2)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 9
        CALL XERCLR
        CALL CTRSM('R','U','N','N',0,2,alpha,a,1,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 9
        CALL XERCLR
        CALL CTRSM('R','U','C','N',0,2,alpha,a,1,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 9
        CALL XERCLR
        CALL CTRSM('R','U','T','N',0,2,alpha,a,1,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 9
        CALL XERCLR
        CALL CTRSM('L','L','N','N',2,0,alpha,a,1,b,2)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 9
        CALL XERCLR
        CALL CTRSM('L','L','C','N',2,0,alpha,a,1,b,2)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 9
        CALL XERCLR
        CALL CTRSM('L','L','T','N',2,0,alpha,a,1,b,2)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 9
        CALL XERCLR
        CALL CTRSM('R','L','N','N',0,2,alpha,a,1,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 9
        CALL XERCLR
        CALL CTRSM('R','L','C','N',0,2,alpha,a,1,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 9
        CALL XERCLR
        CALL CTRSM('R','L','T','N',0,2,alpha,a,1,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 11
        CALL XERCLR
        CALL CTRSM('L','U','N','N',2,0,alpha,a,2,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 11
        CALL XERCLR
        CALL CTRSM('L','U','C','N',2,0,alpha,a,2,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 11
        CALL XERCLR
        CALL CTRSM('L','U','T','N',2,0,alpha,a,2,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 11
        CALL XERCLR
        CALL CTRSM('R','U','N','N',2,0,alpha,a,1,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 11
        CALL XERCLR
        CALL CTRSM('R','U','C','N',2,0,alpha,a,1,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 11
        CALL XERCLR
        CALL CTRSM('R','U','T','N',2,0,alpha,a,1,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 11
        CALL XERCLR
        CALL CTRSM('L','L','N','N',2,0,alpha,a,2,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 11
        CALL XERCLR
        CALL CTRSM('L','L','C','N',2,0,alpha,a,2,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 11
        CALL XERCLR
        CALL CTRSM('L','L','T','N',2,0,alpha,a,2,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 11
        CALL XERCLR
        CALL CTRSM('R','L','N','N',2,0,alpha,a,1,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 11
        CALL XERCLR
        CALL CTRSM('R','L','C','N',2,0,alpha,a,1,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 11
        CALL XERCLR
        CALL CTRSM('R','L','T','N',2,0,alpha,a,1,b,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      CASE (6)
        infot = 1
        CALL XERCLR
        CALL CHERK('/','N',0,0,ralpha,a,1,rbeta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 2
        CALL XERCLR
        CALL CHERK('U','T',0,0,ralpha,a,1,rbeta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 3
        CALL XERCLR
        CALL CHERK('U','N',-1,0,ralpha,a,1,rbeta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 3
        CALL XERCLR
        CALL CHERK('U','C',-1,0,ralpha,a,1,rbeta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 3
        CALL XERCLR
        CALL CHERK('L','N',-1,0,ralpha,a,1,rbeta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 3
        CALL XERCLR
        CALL CHERK('L','C',-1,0,ralpha,a,1,rbeta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 4
        CALL XERCLR
        CALL CHERK('U','N',0,-1,ralpha,a,1,rbeta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 4
        CALL XERCLR
        CALL CHERK('U','C',0,-1,ralpha,a,1,rbeta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 4
        CALL XERCLR
        CALL CHERK('L','N',0,-1,ralpha,a,1,rbeta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 4
        CALL XERCLR
        CALL CHERK('L','C',0,-1,ralpha,a,1,rbeta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 7
        CALL XERCLR
        CALL CHERK('U','N',2,0,ralpha,a,1,rbeta,c,2)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 7
        CALL XERCLR
        CALL CHERK('U','C',0,2,ralpha,a,1,rbeta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 7
        CALL XERCLR
        CALL CHERK('L','N',2,0,ralpha,a,1,rbeta,c,2)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 7
        CALL XERCLR
        CALL CHERK('L','C',0,2,ralpha,a,1,rbeta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 10
        CALL XERCLR
        CALL CHERK('U','N',2,0,ralpha,a,2,rbeta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 10
        CALL XERCLR
        CALL CHERK('U','C',2,0,ralpha,a,1,rbeta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 10
        CALL XERCLR
        CALL CHERK('L','N',2,0,ralpha,a,2,rbeta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 10
        CALL XERCLR
        CALL CHERK('L','C',2,0,ralpha,a,1,rbeta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      CASE (7)
        infot = 1
        CALL XERCLR
        CALL CSYRK('/','N',0,0,alpha,a,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 2
        CALL XERCLR
        CALL CSYRK('U','C',0,0,alpha,a,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 3
        CALL XERCLR
        CALL CSYRK('U','N',-1,0,alpha,a,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 3
        CALL XERCLR
        CALL CSYRK('U','T',-1,0,alpha,a,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 3
        CALL XERCLR
        CALL CSYRK('L','N',-1,0,alpha,a,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 3
        CALL XERCLR
        CALL CSYRK('L','T',-1,0,alpha,a,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 4
        CALL XERCLR
        CALL CSYRK('U','N',0,-1,alpha,a,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 4
        CALL XERCLR
        CALL CSYRK('U','T',0,-1,alpha,a,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 4
        CALL XERCLR
        CALL CSYRK('L','N',0,-1,alpha,a,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 4
        CALL XERCLR
        CALL CSYRK('L','T',0,-1,alpha,a,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 7
        CALL XERCLR
        CALL CSYRK('U','N',2,0,alpha,a,1,beta,c,2)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 7
        CALL XERCLR
        CALL CSYRK('U','T',0,2,alpha,a,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 7
        CALL XERCLR
        CALL CSYRK('L','N',2,0,alpha,a,1,beta,c,2)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 7
        CALL XERCLR
        CALL CSYRK('L','T',0,2,alpha,a,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 10
        CALL XERCLR
        CALL CSYRK('U','N',2,0,alpha,a,2,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 10
        CALL XERCLR
        CALL CSYRK('U','T',2,0,alpha,a,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 10
        CALL XERCLR
        CALL CSYRK('L','N',2,0,alpha,a,2,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 10
        CALL XERCLR
        CALL CSYRK('L','T',2,0,alpha,a,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      CASE (8)
        infot = 1
        CALL XERCLR
        CALL CHER2K('/','N',0,0,alpha,a,1,b,1,rbeta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 2
        CALL XERCLR
        CALL CHER2K('U','T',0,0,alpha,a,1,b,1,rbeta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 3
        CALL XERCLR
        CALL CHER2K('U','N',-1,0,alpha,a,1,b,1,rbeta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 3
        CALL XERCLR
        CALL CHER2K('U','C',-1,0,alpha,a,1,b,1,rbeta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 3
        CALL XERCLR
        CALL CHER2K('L','N',-1,0,alpha,a,1,b,1,rbeta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 3
        CALL XERCLR
        CALL CHER2K('L','C',-1,0,alpha,a,1,b,1,rbeta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 4
        CALL XERCLR
        CALL CHER2K('U','N',0,-1,alpha,a,1,b,1,rbeta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 4
        CALL XERCLR
        CALL CHER2K('U','C',0,-1,alpha,a,1,b,1,rbeta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 4
        CALL XERCLR
        CALL CHER2K('L','N',0,-1,alpha,a,1,b,1,rbeta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 4
        CALL XERCLR
        CALL CHER2K('L','C',0,-1,alpha,a,1,b,1,rbeta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 7
        CALL XERCLR
        CALL CHER2K('U','N',2,0,alpha,a,1,b,1,rbeta,c,2)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 7
        CALL XERCLR
        CALL CHER2K('U','C',0,2,alpha,a,1,b,1,rbeta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 7
        CALL XERCLR
        CALL CHER2K('L','N',2,0,alpha,a,1,b,1,rbeta,c,2)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 7
        CALL XERCLR
        CALL CHER2K('L','C',0,2,alpha,a,1,b,1,rbeta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 9
        CALL XERCLR
        CALL CHER2K('U','N',2,0,alpha,a,2,b,1,rbeta,c,2)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 9
        CALL XERCLR
        CALL CHER2K('U','C',0,2,alpha,a,2,b,1,rbeta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 9
        CALL XERCLR
        CALL CHER2K('L','N',2,0,alpha,a,2,b,1,rbeta,c,2)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 9
        CALL XERCLR
        CALL CHER2K('L','C',0,2,alpha,a,2,b,1,rbeta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 12
        CALL XERCLR
        CALL CHER2K('U','N',2,0,alpha,a,2,b,2,rbeta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 12
        CALL XERCLR
        CALL CHER2K('U','C',2,0,alpha,a,1,b,1,rbeta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 12
        CALL XERCLR
        CALL CHER2K('L','N',2,0,alpha,a,2,b,2,rbeta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 12
        CALL XERCLR
        CALL CHER2K('L','C',2,0,alpha,a,1,b,1,rbeta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      CASE (9)
        infot = 1
        CALL XERCLR
        CALL CSYR2K('/','N',0,0,alpha,a,1,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 2
        CALL XERCLR
        CALL CSYR2K('U','C',0,0,alpha,a,1,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 3
        CALL XERCLR
        CALL CSYR2K('U','N',-1,0,alpha,a,1,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 3
        CALL XERCLR
        CALL CSYR2K('U','T',-1,0,alpha,a,1,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 3
        CALL XERCLR
        CALL CSYR2K('L','N',-1,0,alpha,a,1,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 3
        CALL XERCLR
        CALL CSYR2K('L','T',-1,0,alpha,a,1,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 4
        CALL XERCLR
        CALL CSYR2K('U','N',0,-1,alpha,a,1,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 4
        CALL XERCLR
        CALL CSYR2K('U','T',0,-1,alpha,a,1,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 4
        CALL XERCLR
        CALL CSYR2K('L','N',0,-1,alpha,a,1,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 4
        CALL XERCLR
        CALL CSYR2K('L','T',0,-1,alpha,a,1,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 7
        CALL XERCLR
        CALL CSYR2K('U','N',2,0,alpha,a,1,b,1,beta,c,2)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 7
        CALL XERCLR
        CALL CSYR2K('U','T',0,2,alpha,a,1,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 7
        CALL XERCLR
        CALL CSYR2K('L','N',2,0,alpha,a,1,b,1,beta,c,2)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 7
        CALL XERCLR
        CALL CSYR2K('L','T',0,2,alpha,a,1,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 9
        CALL XERCLR
        CALL CSYR2K('U','N',2,0,alpha,a,2,b,1,beta,c,2)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 9
        CALL XERCLR
        CALL CSYR2K('U','T',0,2,alpha,a,2,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 9
        CALL XERCLR
        CALL CSYR2K('L','N',2,0,alpha,a,2,b,1,beta,c,2)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 9
        CALL XERCLR
        CALL CSYR2K('L','T',0,2,alpha,a,2,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 12
        CALL XERCLR
        CALL CSYR2K('U','N',2,0,alpha,a,2,b,2,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 12
        CALL XERCLR
        CALL CSYR2K('U','T',2,0,alpha,a,1,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 12
        CALL XERCLR
        CALL CSYR2K('L','N',2,0,alpha,a,2,b,2,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 12
        CALL XERCLR
        CALL CSYR2K('L','T',2,0,alpha,a,1,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      CASE DEFAULT
        infot = 1
        CALL XERCLR
        CALL CGEMM('/','N',0,0,0,alpha,a,1,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 1
        CALL XERCLR
        CALL CGEMM('/','C',0,0,0,alpha,a,1,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 1
        CALL XERCLR
        CALL CGEMM('/','T',0,0,0,alpha,a,1,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 2
        CALL XERCLR
        CALL CGEMM('N','/',0,0,0,alpha,a,1,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 2
        CALL XERCLR
        CALL CGEMM('C','/',0,0,0,alpha,a,1,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 2
        CALL XERCLR
        CALL CGEMM('T','/',0,0,0,alpha,a,1,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 3
        CALL XERCLR
        CALL CGEMM('N','N',-1,0,0,alpha,a,1,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 3
        CALL XERCLR
        CALL CGEMM('N','C',-1,0,0,alpha,a,1,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 3
        CALL XERCLR
        CALL CGEMM('N','T',-1,0,0,alpha,a,1,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 3
        CALL XERCLR
        CALL CGEMM('C','N',-1,0,0,alpha,a,1,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 3
        CALL XERCLR
        CALL CGEMM('C','C',-1,0,0,alpha,a,1,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 3
        CALL XERCLR
        CALL CGEMM('C','T',-1,0,0,alpha,a,1,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 3
        CALL XERCLR
        CALL CGEMM('T','N',-1,0,0,alpha,a,1,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 3
        CALL XERCLR
        CALL CGEMM('T','C',-1,0,0,alpha,a,1,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 3
        CALL XERCLR
        CALL CGEMM('T','T',-1,0,0,alpha,a,1,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 4
        CALL XERCLR
        CALL CGEMM('N','N',0,-1,0,alpha,a,1,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 4
        CALL XERCLR
        CALL CGEMM('N','C',0,-1,0,alpha,a,1,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 4
        CALL XERCLR
        CALL CGEMM('N','T',0,-1,0,alpha,a,1,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 4
        CALL XERCLR
        CALL CGEMM('C','N',0,-1,0,alpha,a,1,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 4
        CALL XERCLR
        CALL CGEMM('C','C',0,-1,0,alpha,a,1,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 4
        CALL XERCLR
        CALL CGEMM('C','T',0,-1,0,alpha,a,1,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 4
        CALL XERCLR
        CALL CGEMM('T','N',0,-1,0,alpha,a,1,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 4
        CALL XERCLR
        CALL CGEMM('T','C',0,-1,0,alpha,a,1,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 4
        CALL XERCLR
        CALL CGEMM('T','T',0,-1,0,alpha,a,1,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 5
        CALL XERCLR
        CALL CGEMM('N','N',0,0,-1,alpha,a,1,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 5
        CALL XERCLR
        CALL CGEMM('N','C',0,0,-1,alpha,a,1,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 5
        CALL XERCLR
        CALL CGEMM('N','T',0,0,-1,alpha,a,1,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 5
        CALL XERCLR
        CALL CGEMM('C','N',0,0,-1,alpha,a,1,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 5
        CALL XERCLR
        CALL CGEMM('C','C',0,0,-1,alpha,a,1,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 5
        CALL XERCLR
        CALL CGEMM('C','T',0,0,-1,alpha,a,1,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 5
        CALL XERCLR
        CALL CGEMM('T','N',0,0,-1,alpha,a,1,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 5
        CALL XERCLR
        CALL CGEMM('T','C',0,0,-1,alpha,a,1,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 5
        CALL XERCLR
        CALL CGEMM('T','T',0,0,-1,alpha,a,1,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 8
        CALL XERCLR
        CALL CGEMM('N','N',2,0,0,alpha,a,1,b,1,beta,c,2)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 8
        CALL XERCLR
        CALL CGEMM('N','C',2,0,0,alpha,a,1,b,1,beta,c,2)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 8
        CALL XERCLR
        CALL CGEMM('N','T',2,0,0,alpha,a,1,b,1,beta,c,2)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 8
        CALL XERCLR
        CALL CGEMM('C','N',0,0,2,alpha,a,1,b,2,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 8
        CALL XERCLR
        CALL CGEMM('C','C',0,0,2,alpha,a,1,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 8
        CALL XERCLR
        CALL CGEMM('C','T',0,0,2,alpha,a,1,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 8
        CALL XERCLR
        CALL CGEMM('T','N',0,0,2,alpha,a,1,b,2,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 8
        CALL XERCLR
        CALL CGEMM('T','C',0,0,2,alpha,a,1,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 8
        CALL XERCLR
        CALL CGEMM('T','T',0,0,2,alpha,a,1,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 10
        CALL XERCLR
        CALL CGEMM('N','N',0,0,2,alpha,a,1,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 10
        CALL XERCLR
        CALL CGEMM('C','N',0,0,2,alpha,a,2,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 10
        CALL XERCLR
        CALL CGEMM('T','N',0,0,2,alpha,a,2,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 10
        CALL XERCLR
        CALL CGEMM('N','C',0,2,0,alpha,a,1,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 10
        CALL XERCLR
        CALL CGEMM('C','C',0,2,0,alpha,a,1,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 10
        CALL XERCLR
        CALL CGEMM('T','C',0,2,0,alpha,a,1,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 10
        CALL XERCLR
        CALL CGEMM('N','T',0,2,0,alpha,a,1,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 10
        CALL XERCLR
        CALL CGEMM('C','T',0,2,0,alpha,a,1,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 10
        CALL XERCLR
        CALL CGEMM('T','T',0,2,0,alpha,a,1,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 13
        CALL XERCLR
        CALL CGEMM('N','N',2,0,0,alpha,a,2,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 13
        CALL XERCLR
        CALL CGEMM('N','C',2,0,0,alpha,a,2,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 13
        CALL XERCLR
        CALL CGEMM('N','T',2,0,0,alpha,a,2,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 13
        CALL XERCLR
        CALL CGEMM('C','N',2,0,0,alpha,a,1,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 13
        CALL XERCLR
        CALL CGEMM('C','C',2,0,0,alpha,a,1,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 13
        CALL XERCLR
        CALL CGEMM('C','T',2,0,0,alpha,a,1,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 13
        CALL XERCLR
        CALL CGEMM('T','N',2,0,0,alpha,a,1,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 13
        CALL XERCLR
        CALL CGEMM('T','C',2,0,0,alpha,a,1,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 13
        CALL XERCLR
        CALL CGEMM('T','T',2,0,0,alpha,a,1,b,1,beta,c,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
    END SELECT
    !
    IF ( Kprint>=2 ) THEN
      CALL XERDMP
      IF ( .NOT.Fatal ) THEN
        WRITE (Nout,FMT=99001) Srnamt
      ELSE
        WRITE (Nout,FMT=99002) Srnamt
      ENDIF
    ENDIF
    CALL XSETF(kontrl)
    RETURN
    !
    99001 FORMAT (' ',A6,' PASSED THE TESTS OF ERROR-EXITS')
    99002 FORMAT (' ******* ',A6,' FAILED THE TESTS OF ERROR-EXITS *****','**')
    !
    !     End of CCHKE3.
    !
  END SUBROUTINE CCHKE3
END MODULE TEST20_MOD
!** TEST20
PROGRAM TEST20
  USE TEST20_MOD
  IMPLICIT NONE
  !>
  !***
  !  Driver for testing SLATEC subprograms
  !***
  ! **Library:**   SLATEC
  !***
  ! **Category:**  D1B
  !***
  ! **Keywords:**  QUICK CHECK DRIVER
  !***
  ! **Type:**      COMPLEX (TEST18-S, TEST19-D, TEST20-C)
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
  !        complex Levels 2 and 3 BLAS routines
  !
  !***
  ! **References:**  Kirby W. Fong,  Thomas H. Jefferson, Tokihiko Suyehiro
  !                 and Lee Walton, Guide to the SLATEC Common Mathema-
  !                 tical Library, April 10, 1990.
  !***
  ! **Routines called:**  I1MACH, CBLAT2, CBLAT3, XERMAX, XSETF, XSETUN

  !* REVISION HISTORY  (YYMMDD)
  !   920601  DATE WRITTEN

  INTEGER ipass, kprint, lin, lun, nfail
  !     .. External Functions ..
  INTEGER, EXTERNAL :: I1MACH
  !* FIRST EXECUTABLE STATEMENT  TEST20
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
  ENDIF
  !
  !     Test complex Level 2 BLAS routines
  !
  CALL CBLAT2(lun,kprint,ipass)
  IF ( ipass==0 ) nfail = nfail + 1
  !
  !     Test complex Level 3 BLAS routines
  !
  CALL CBLAT3(lun,kprint,ipass)
  IF ( ipass==0 ) nfail = nfail + 1
  !
  !     Write PASS or FAIL message
  !
  IF ( nfail==0 ) THEN
    WRITE (lun,99001)
    99001 FORMAT (/' --------------TEST20 PASSED ALL TESTS----------------')
  ELSE
    WRITE (lun,99002) nfail
    99002 FORMAT (/' ************* WARNING -- ',I5,&
      ' TEST(S) FAILED IN PROGRAM TEST20 *************')
  ENDIF
  STOP
END PROGRAM TEST20
