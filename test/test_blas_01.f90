MODULE TEST17_MOD
  IMPLICIT NONE

CONTAINS
  !** BLACHK
  SUBROUTINE BLACHK(Lun,Kprint,Ipass)
    IMPLICIT NONE
    !>
    !***
    !  Quick check for Basic Linear Algebra Subprograms.
    !***
    ! **Library:**   SLATEC
    !***
    ! **Keywords:**  QUICK CHECK
    !***
    ! **Author:**  Lawson, C. L., (JPL)
    !***
    ! **Description:**
    !
    !     ********************************* TBLA ***************************
    !     TEST DRIVER FOR BASIC LINEAR ALGEBRA SUBPROGRAMS.
    !     C. L. LAWSON, JPL, 1974 DEC 10, 1975 MAY 28
    !
    !     UPDATED BY K. HASKELL - JUNE 23,1980
    !
    !***
    ! **Routines called:**  CHECK0, CHECK1, CHECK2, HEADER
    !***
    ! COMMON BLOCKS    COMBLA

    !* REVISION HISTORY  (YYMMDD)
    !   751210  DATE WRITTEN
    !   890618  REVISION DATE from Version 3.2
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    
    INTEGER ICAse, INCx, INCy, Kprint, Lun, MODe, N, NPRint
    REAL sdfac, sfac
    INTEGER Ipass, jtest(38)
    REAL(8) :: dfac, dqfac
    LOGICAL PASs
    COMMON /COMBLA/ NPRint, ICAse, N, INCx, INCy, MODe, PASs
    DATA sfac, sdfac, dfac, dqfac/.625E-1, .50, .625D-1, 0.625D-1/
    DATA jtest/1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, &
      1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, &
      1, 1, 1, 1, 1/
    !* FIRST EXECUTABLE STATEMENT  BLACHK
    NPRint = Lun
    Ipass = 1
    !
    IF ( Kprint>=2 ) WRITE (NPRint,99001)
    99001 FORMAT ('1','QUICK CHECK OF 38 BASIC LINEAR ALGEBRA SUBROUTINES'/)
    DO ICAse = 1, 38
      IF ( jtest(ICAse)/=0 ) THEN
        CALL HEADER(Kprint)
        !
        !         INITIALIZE  PASS, INCX, INCY, AND MODE FOR A NEW CASE.
        !         THE VALUE 9999 FOR INCX, INCY OR MODE WILL APPEAR IN THE
        !         DETAILED  OUTPUT, IF ANY, FOR CASES THAT DO NOT INVOLVE
        !         THESE PARAMETERS.
        !
        PASs = .TRUE.
        INCx = 9999
        INCy = 9999
        MODe = 9999
        SELECT CASE (ICAse)
          CASE (1,2,3,4,5,6,7,8,9,10,11,14,15,18,19,20,21,22,23,24,25)
            ! ICASE =  1-11, 14-15, OR 18-25
            CALL CHECK2(sfac,sdfac,dfac,dqfac,Kprint)
          CASE (26,27,28,29,30,31,32,33,34,35,36,37,38)
            ! ICASE = 26-38
            CALL CHECK1(sfac,dfac,Kprint)
          CASE DEFAULT
            ! ICASE = 12-13 OR 16-17
            CALL CHECK0(sfac,dfac,Kprint)
        END SELECT
        !                                                  PRINT
        IF ( Kprint>=2.AND.PASs ) WRITE (NPRint,99002)
        99002 FORMAT ('+',39X,'PASS')
        IF ( .NOT.PASs ) Ipass = 0
      ENDIF
    ENDDO
    IF ( Kprint>=2.AND.Ipass==1 ) WRITE (NPRint,99003)
    99003 FORMAT (/' ****************BLAS PASSED ALL TESTS****************')
    IF ( Kprint>=1.AND.Ipass==0 ) WRITE (NPRint,99004)
    99004 FORMAT (/' ****************BLAS FAILED SOME TESTS***************')
    RETURN
  END SUBROUTINE BLACHK
  !** HEADER
  SUBROUTINE HEADER(Kprint)
    IMPLICIT NONE
    !>
    !***
    !  Print header for BLAS quick checks.
    !***
    ! **Library:**   SLATEC
    !***
    ! **Author:**  Lawson, C. L., (JPL)
    !***
    ! **Routines called:**  (NONE)
    !***
    ! COMMON BLOCKS    COMBLA

    !* REVISION HISTORY  (YYMMDD)
    !   741212  DATE WRITTEN
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    !   920210  Minor modifications to prologue and code.  (WRB)
    
    INTEGER ICAse, INCx, INCy, Kprint, MODe, N, NPRint
    COMMON /COMBLA/ NPRint, ICAse, N, INCx, INCy, MODe, PASs
    LOGICAL PASs
    CHARACTER(6) :: l(38)
    !
    DATA l(1)/'  SDOT'/
    DATA l(2)/' DSDOT'/
    DATA l(3)/'SDSDOT'/
    DATA l(4)/'  DDOT'/
    DATA l(5)/'DQDOTI'/
    DATA l(6)/'DQDOTA'/
    DATA l(7)/' CDOTC'/
    DATA l(8)/' CDOTU'/
    DATA l(9)/' SAXPY'/
    DATA l(10)/' DAXPY'/
    DATA l(11)/' CAXPY'/
    DATA l(12)/' SROTG'/
    DATA l(13)/' DROTG'/
    DATA l(14)/'  SROT'/
    DATA l(15)/'  DROT'/
    DATA l(16)/'SROTMG'/
    DATA l(17)/'DROTMG'/
    DATA l(18)/' SROTM'/
    DATA l(19)/' DROTM'/
    DATA l(20)/' SCOPY'/
    DATA l(21)/' DCOPY'/
    DATA l(22)/' CCOPY'/
    DATA l(23)/' SSWAP'/
    DATA l(24)/' DSWAP'/
    DATA l(25)/' CSWAP'/
    DATA l(26)/' SNRM2'/
    DATA l(27)/' DNRM2'/
    DATA l(28)/'SCNRM2'/
    DATA l(29)/' SASUM'/
    DATA l(30)/' DASUM'/
    DATA l(31)/'SCASUM'/
    DATA l(32)/' SSCAL'/
    DATA l(33)/' DSCAL'/
    DATA l(34)/' CSCAL'/
    DATA l(35)/'CSSCAL'/
    DATA l(36)/'ISAMAX'/
    DATA l(37)/'IDAMAX'/
    DATA l(38)/'ICAMAX'/
    !* FIRST EXECUTABLE STATEMENT  HEADER
    IF ( Kprint>=2 ) WRITE (NPRint,99001) ICAse, l(ICAse)
    !
    99001 FORMAT (' Test of subprogram number',I3,2X,A)
    RETURN
  END SUBROUTINE HEADER
  !** CHECK0
  SUBROUTINE CHECK0(Sfac,Dfac,Kprint)
    IMPLICIT NONE
    !>
    !***
    !  (UNKNOWN)
    !***
    ! **Library:**   SLATEC
    !***
    ! **Author:**  Lawson, C. L., (JPL)
    !           Hanson, R. J., (SNLA)
    !           Wisniewski, J. A., (SNLA)
    !***
    ! **Description:**
    !
    !     THIS SUBROUTINE TESTS SUBPROGRAMS 12-13 AND 16-17.
    !     THESE SUBPROGRAMS HAVE NO ARRAY ARGUMENTS.
    !
    !     C. L. LAWSON, JPL, 1975 MAR 07, MAY 28
    !     R. J. HANSON, J. A. WISNIEWSKI, SANDIA LABS, APRIL 25,1977.
    !
    !***
    ! **Routines called:**  DROTG, DROTMG, DTEST, SROTG, SROTMG, STEST
    !***
    ! COMMON BLOCKS    COMBLA

    !* REVISION HISTORY  (YYMMDD)
    !   750307  DATE WRITTEN
    !   890911  Removed unnecessary intrinsics.  (WRB)
    !   890911  REVISION DATE from Version 3.2
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    
    INTEGER i, ICAse, INCx, INCy, jump, k, Kprint, MODe, N, NPRint
    REAL sa, sb, sc, Sfac, ss, zero
    COMMON /COMBLA/ NPRint, ICAse, N, INCx, INCy, MODe, PASs
    LOGICAL PASs
    REAL strue(9), stemp(9)
    REAL(8) :: dc, ds, da1(8), db1(8), dc1(8), ds1(8)
    REAL(8) :: da, datrue(8), dbtrue(8), dzero, Dfac, db
    REAL(8) :: dab(4,9), dtemp(9), dtrue(9,9), d12
    DATA zero, dzero/0., 0.D0/
    DATA da1/.3D0, .4D0, -.3D0, -.4D0, -.3D0, 0.D0, 0.D0, 1.D0/
    DATA db1/.4D0, .3D0, .4D0, .3D0, -.4D0, 0.D0, 1.D0, 0.D0/
    DATA dc1/.6D0, .8D0, -.6D0, .8D0, .6D0, 1.D0, 0.D0, 1.D0/
    DATA ds1/.8D0, .6D0, .8D0, -.6D0, .8D0, 0.D0, 1.D0, 0.D0/
    DATA datrue/.5D0, .5D0, .5D0, -.5D0, -.5D0, 0.D0, 1.D0, 1.D0/
    DATA dbtrue/0.D0, .6D0, 0.D0, -.6D0, 0.D0, 0.D0, 1.D0, 0.D0/
    !                                              INPUT FOR MODIFIED GIVENS
    DATA dab/.1D0, .3D0, 1.2D0, .2D0, .7D0, .2D0, .6D0, 4.2D0, 0.D0, &
      0.D0, 0.D0, 0.D0, 4.D0, -1.D0, 2.D0, 4.D0, 6.D-10, 2.D-2, &
      1.D5, 10.D0, 4.D10, 2.D-2, 1.D-5, 10.D0, 2.D-10, 4.D-2, &
      1.D5, 10.D0, 2.D10, 4.D-2, 1.D-5, 10.D0, 4.D0, -2.D0, 8.D0, &
      4.D0/
    !                                       TRUE RESULTS FOR MODIFIED GIVENS
    DATA dtrue/0.D0, 0.D0, 1.3D0, .2D0, 0.D0, 0.D0, 0.D0, .5D0, 0.D0, &
      0.D0, 0.D0, 4.5D0, 4.2D0, 1.D0, .5D0, 0.D0, 0.D0, 0.D0, &
      0.D0, 0.D0, 0.D0, 0.D0, -2.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      0.D0, 0.D0, 0.D0, 4.D0, -1.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      0.D0, 15.D-3, 0.D0, 10.D0, -1.D0, 0.D0, -1.D-4, 0.D0, 1.D0, &
      0.D0, 0.D0, 6144.D-5, 10.D0, -1.D0, 4096.D0, -1.D6, 0.D0, &
      1.D0, 0.D0, 0.D0, 15.D0, 10.D0, -1.D0, 5.D-5, 0.D0, 1.D0, &
      0.D0, 0.D0, 0.D0, 15.D0, 10.D0, -1.D0, 5.D5, -4096.D0, 1.D0, &
      4096.D-6, 0.D0, 0.D0, 7.D0, 4.D0, 0.D0, 0.D0, -.5D0, -.25D0, &
      0.D0/
    !                   4096 = 2 ** 12
    DATA d12/4096.D0/
    !* FIRST EXECUTABLE STATEMENT  CHECK0
    !
    !                   COMPUTE TRUE VALUES WHICH CANNOT BE PRESTORED
    !                   IN DECIMAL NOTATION.
    dtrue(1,1) = 12.D0/130.D0
    dtrue(2,1) = 36.D0/130.D0
    dtrue(7,1) = -1.D0/6.D0
    dtrue(1,2) = 14.D0/75.D0
    dtrue(2,2) = 49.D0/75.D0
    dtrue(9,2) = 1.D0/7.D0
    dtrue(1,5) = 45.D-11*(d12*d12)
    dtrue(3,5) = 4.D5/(3.D0*d12)
    dtrue(6,5) = 1.D0/d12
    dtrue(8,5) = 1.D4/(3.D0*d12)
    dtrue(1,6) = 4.D10/(1.5D0*d12*d12)
    dtrue(2,6) = 2.D-2/1.5D0
    dtrue(8,6) = 5.D-7*d12
    dtrue(1,7) = 4.D0/150.D0
    dtrue(2,7) = (2.D-10/1.5D0)*(d12*d12)
    dtrue(7,7) = -dtrue(6,5)
    dtrue(9,7) = 1.D4/d12
    dtrue(1,8) = dtrue(1,7)
    dtrue(2,8) = 2.D10/(1.5D0*d12*d12)
    dtrue(1,9) = 32.D0/7.D0
    dtrue(2,9) = -16.D0/7.D0
    dbtrue(1) = 1.D0/.6D0
    dbtrue(3) = -1.D0/.6D0
    dbtrue(5) = 1.D0/.6D0
    !
    jump = ICAse - 11
    DO k = 1, 9
      !                        SET N=K FOR IDENTIFICATION IN OUTPUT IF ANY.
      N = k
      !                             BRANCH TO SELECT SUBPROGRAM TO BE TESTED.
      !
      SELECT CASE (jump)
        CASE (2)
          ! 13. DROTG
          IF ( k>8 ) EXIT
          da = da1(k)
          db = db1(k)
          CALL DROTG(da,db,dc,ds)
          CALL DTEST(1,[da],[datrue(k)],datrue(k),Dfac,Kprint)
          CALL DTEST(1,[db],[dbtrue(k)],dbtrue(k),Dfac,Kprint)
          CALL DTEST(1,[dc],[dc1(k)],dc1(k),Dfac,Kprint)
          CALL DTEST(1,[ds],[ds1(k)],ds1(k),Dfac,Kprint)
        CASE (3,4)
          GOTO 100
        CASE (5)
          ! 16. SROTMG
          DO i = 1, 4
            stemp(i) = REAL( dab(i,k), 4 )
            stemp(i+4) = zero
          ENDDO
          stemp(9) = zero
          CALL SROTMG(stemp(1),stemp(2),stemp(3),stemp(4),stemp(5))
          !
          DO i = 1, 9
            strue(i) = REAL( dtrue(i,k), 4 )
          ENDDO
          CALL STEST(9,stemp,strue,strue,Sfac,Kprint)
        CASE (6)
          ! 17. DROTMG
          DO i = 1, 4
            dtemp(i) = dab(i,k)
            dtemp(i+4) = dzero
          ENDDO
          dtemp(9) = dzero
          CALL DROTMG(dtemp(1),dtemp(2),dtemp(3),dtemp(4),dtemp(5))
          CALL DTEST(9,dtemp,dtrue(1,k),dtrue(1,k),Dfac,Kprint)
        CASE DEFAULT
          ! 12. SROTG
          IF ( k>8 ) EXIT
          sa = REAL( da1(k), 4 )
          sb = REAL( db1(k), 4 )
          CALL SROTG(sa,sb,sc,ss)
          CALL STEST(1,[sa],[REAL(datrue(k))],[REAL(datrue(k))],Sfac,Kprint)
          CALL STEST(1,[sb],[REAL(dbtrue(k))],[REAL(dbtrue(k))],Sfac,Kprint)
          CALL STEST(1,[sc],[REAL(dc1(k))],[REAL(dc1(k))],Sfac,Kprint)
          CALL STEST(1,[ss],[REAL(ds1(k))],[REAL(ds1(k))],Sfac,Kprint)
      END SELECT
    ENDDO
    RETURN
    !                     THE FOLLOWING STOP SHOULD NEVER BE REACHED.
    100  STOP
  END SUBROUTINE CHECK0
  !** CHECK1
  SUBROUTINE CHECK1(Sfac,Dfac,Kprint)
    IMPLICIT NONE
    !>
    !***
    !  (UNKNOWN)
    !***
    ! **Library:**   SLATEC
    !***
    ! **Author:**  Lawson, C. L., (JPL)
    !***
    ! **Description:**
    !
    !     THIS SUBPROGRAM TESTS THE INCREMENTING AND ACCURACY OF THE LINEAR
    !     ALGEBRA SUBPROGRAMS 26 - 38 (SNRM2 TO ICAMAX). STORED RESULTS ARE
    !     COMPARED WITH THE RESULT RETURNED BY THE SUBPROGRAM.
    !
    !     THESE SUBPROGRAMS REQUIRE A SINGLE VECTOR ARGUMENT.
    !
    !     ICASE            DESIGNATES WHICH SUBPROGRAM TO TEST.
    !                      26 .LE. ICASE .LE. 38
    !     C. L. LAWSON, JPL, 1974 DEC 10, MAY 28
    !
    !***
    ! **Routines called:**  CSCAL, CSSCAL, DASUM, DNRM2, DSCAL, DTEST, ICAMAX,
    !                    IDAMAX, ISAMAX, ITEST, SASUM, SCASUM, SCNRM2,
    !                    SNRM2, SSCAL, STEST
    !***
    ! COMMON BLOCKS    COMBLA

    !* REVISION HISTORY  (YYMMDD)
    !   741210  DATE WRITTEN
    !   890911  Removed unnecessary intrinsics.  (WRB)
    !   890911  REVISION DATE from Version 3.2
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    
    INTEGER i, ICAMAX, ICAse, IDAMAX, INCx, INCy, ISAMAX, jump, &
      Kprint, len, MODe, N, np1, NPRint
    REAL sa, SASUM, SCASUM, SCNRM2, Sfac, SNRM2, stemp
    COMMON /COMBLA/ NPRint, ICAse, N, INCx, INCy, MODe, PASs
    LOGICAL PASs
    INTEGER itrue2(5), itrue3(5)
    REAL(8) :: da, dx(8)
    REAL(8) :: dv(8,5,2)
    REAL(8) :: Dfac
    REAL(8) :: DNRM2, DASUM
    REAL(8) :: dtrue1(5), dtrue3(5), dtrue5(8,5,2)
    REAL strue2(5), strue4(5), strue(8), sx(8)
    COMPLEX ca, cv(8,5,2), ctrue5(8,5,2), ctrue6(8,5,2), cx(8)
    !
    DATA sa, da, ca/.3, .3D0, (.4,-.7)/
    DATA dv/.1D0, 2.D0, 2.D0, 2.D0, 2.D0, 2.D0, 2.D0, 2.D0, .3D0, &
      3.D0, 3.D0, 3.D0, 3.D0, 3.D0, 3.D0, 3.D0, .3D0, -.4D0, &
      4.D0, 4.D0, 4.D0, 4.D0, 4.D0, 4.D0, .2D0, -.6D0, .3D0, &
      5.D0, 5.D0, 5.D0, 5.D0, 5.D0, .1D0, -.3D0, .5D0, -.1D0, &
      6.D0, 6.D0, 6.D0, 6.D0, .1D0, 8.D0, 8.D0, 8.D0, 8.D0, 8.D0, &
      8.D0, 8.D0, .3D0, 9.D0, 9.D0, 9.D0, 9.D0, 9.D0, 9.D0, 9.D0, &
      .3D0, 2.D0, -.4D0, 2.D0, 2.D0, 2.D0, 2.D0, 2.D0, .2D0, &
      3.D0, -.6D0, 5.D0, .3D0, 2.D0, 2.D0, 2.D0, .1D0, 4.D0, &
      -.3D0, 6.D0, -.5D0, 7.D0, -.1D0, 3.D0/
    !     COMPLEX TEST VECTORS
    DATA cv/(.1,.1), (1.,2.), (1.,2.), (1.,2.), (1.,2.), (1.,2.), &
      (1.,2.), (1.,2.), (.3,-.4), (3.,4.), (3.,4.), (3.,4.), (3.,4.)&
      , (3.,4.), (3.,4.), (3.,4.), (.1,-.3), (.5,-.1), (5.,6.), &
      (5.,6.), (5.,6.), (5.,6.), (5.,6.), (5.,6.), (.1,.1), (-.6,.1)&
      , (.1,-.3), (7.,8.), (7.,8.), (7.,8.), (7.,8.), (7.,8.), &
      (.3,.1), (.1,.4), (.4,.1), (.1,.2), (2.,3.), (2.,3.), (2.,3.), &
      (2.,3.), (.1,.1), (4.,5.), (4.,5.), (4.,5.), (4.,5.), (4.,5.), &
      (4.,5.), (4.,5.), (.3,-.4), (6.,7.), (6.,7.), (6.,7.), (6.,7.)&
      , (6.,7.), (6.,7.), (6.,7.), (.1,-.3), (8.,9.), (.5,-.1), &
      (2.,5.), (2.,5.), (2.,5.), (2.,5.), (2.,5.), (.1,.1), (3.,6.), &
      (-.6,.1), (4.,7.), (.1,-.3), (7.,2.), (7.,2.), (7.,2.), (.3,.1)&
      , (5.,8.), (.1,.4), (6.,9.), (.4,.1), (8.,3.), (.1,.2), (9.,4.)&
      /
    !
    DATA strue2/.0, .5, .6, .7, .7/
    DATA strue4/.0, .7, 1., 1.3, 1.7/
    DATA dtrue1/.0D0, .3D0, .5D0, .7D0, .6D0/
    DATA dtrue3/.0D0, .3D0, .7D0, 1.1D0, 1.D0/
    DATA dtrue5/.10D0, 2.D0, 2.D0, 2.D0, 2.D0, 2.D0, 2.D0, 2.D0, &
      .09D0, 3.D0, 3.D0, 3.D0, 3.D0, 3.D0, 3.D0, 3.D0, .09D0, &
      -.12D0, 4.D0, 4.D0, 4.D0, 4.D0, 4.D0, 4.D0, .06D0, -.18D0, &
      .09D0, 5.D0, 5.D0, 5.D0, 5.D0, 5.D0, .03D0, -.09D0, .15D0, &
      -.03D0, 6.D0, 6.D0, 6.D0, 6.D0, .10D0, 8.D0, 8.D0, 8.D0, &
      8.D0, 8.D0, 8.D0, 8.D0, .09D0, 9.D0, 9.D0, 9.D0, 9.D0, &
      9.D0, 9.D0, 9.D0, .09D0, 2.D0, -.12D0, 2.D0, 2.D0, 2.D0, &
      2.D0, 2.D0, .06D0, 3.D0, -.18D0, 5.D0, .09D0, 2.D0, 2.D0, &
      2.D0, .03D0, 4.D0, -.09D0, 6.D0, -.15D0, 7.D0, -.03D0, 3.D0/
    !
    DATA ctrue5/(.1,.1), (1.,2.), (1.,2.), (1.,2.), (1.,2.), (1.,2.), &
      (1.,2.), (1.,2.), (-.16,-.37), (3.,4.), (3.,4.), (3.,4.), &
      (3.,4.), (3.,4.), (3.,4.), (3.,4.), (-.17,-.19), (.13,-.39), &
      (5.,6.), (5.,6.), (5.,6.), (5.,6.), (5.,6.), (5.,6.), &
      (.11,-.03), (-.17,.46), (-.17,-.19), (7.,8.), (7.,8.), (7.,8.), &
      (7.,8.), (7.,8.), (.19,-.17), (.32,.09), (.23,-.24), (.18,.01), &
      (2.,3.), (2.,3.), (2.,3.), (2.,3.), (.1,.1), (4.,5.), (4.,5.), &
      (4.,5.), (4.,5.), (4.,5.), (4.,5.), (4.,5.), (-.16,-.37), &
      (6.,7.), (6.,7.), (6.,7.), (6.,7.), (6.,7.), (6.,7.), (6.,7.), &
      (-.17,-.19), (8.,9.), (.13,-.39), (2.,5.), (2.,5.), (2.,5.), &
      (2.,5.), (2.,5.), (.11,-.03), (3.,6.), (-.17,.46), (4.,7.), &
      (-.17,-.19), (7.,2.), (7.,2.), (7.,2.), (.19,-.17), (5.,8.), &
      (.32,.09), (6.,9.), (.23,-.24), (8.,3.), (.18,.01), (9.,4.)/
    !
    DATA ctrue6/(.1,.1), (1.,2.), (1.,2.), (1.,2.), (1.,2.), (1.,2.), &
      (1.,2.), (1.,2.), (.09,-.12), (3.,4.), (3.,4.), (3.,4.), &
      (3.,4.), (3.,4.), (3.,4.), (3.,4.), (.03,-.09), (.15,-.03), &
      (5.,6.), (5.,6.), (5.,6.), (5.,6.), (5.,6.), (5.,6.), (.03,.03)&
      , (-.18,.03), (.03,-.09), (7.,8.), (7.,8.), (7.,8.), (7.,8.), &
      (7.,8.), (.09,.03), (.03,.12), (.12,.03), (.03,.06), (2.,3.), &
      (2.,3.), (2.,3.), (2.,3.), (.1,.1), (4.,5.), (4.,5.), (4.,5.), &
      (4.,5.), (4.,5.), (4.,5.), (4.,5.), (.09,-.12), (6.,7.), &
      (6.,7.), (6.,7.), (6.,7.), (6.,7.), (6.,7.), (6.,7.), &
      (.03,-.09), (8.,9.), (.15,-.03), (2.,5.), (2.,5.), (2.,5.), &
      (2.,5.), (2.,5.), (.03,.03), (3.,6.), (-.18,.03), (4.,7.), &
      (.03,-.09), (7.,2.), (7.,2.), (7.,2.), (.09,.03), (5.,8.), &
      (.03,.12), (6.,9.), (.12,.03), (8.,3.), (.03,.06), (9.,4.)/
    !
    !
    DATA itrue2/0, 1, 2, 2, 3/
    DATA itrue3/0, 1, 2, 2, 2/
    !* FIRST EXECUTABLE STATEMENT  CHECK1
    jump = ICAse - 25
    DO INCx = 1, 2
      DO np1 = 1, 5
        N = np1 - 1
        len = 2*MAX(N,1)
        !                                                  SET VECTOR ARGUMENTS.
        DO i = 1, len
          sx(i) = REAL( dv(i,np1,INCx), 4 )
          dx(i) = dv(i,np1,INCx)
          cx(i) = cv(i,np1,INCx)
        ENDDO
        !
        !                        BRANCH TO INVOKE SUBPROGRAM TO BE TESTED.
        !
        SELECT CASE (jump)
          CASE (2)
            ! 27. DNRM2
            CALL DTEST(1,[DNRM2(N,dx,INCx)],dtrue1(np1),dtrue1(np1),Dfac,Kprint)
          CASE (3)
            ! 28. SCNRM2
            CALL STEST(1,[SCNRM2(N,cx,INCx)],strue2(np1),strue2(np1),Sfac,Kprint)
          CASE (4)
            ! 29. SASUM
            stemp = REAL( dtrue3(np1), 4 )
            CALL STEST(1,[SASUM(N,sx,INCx)],[stemp],[stemp],Sfac,Kprint)
          CASE (5)
            ! 30. DASUM
            CALL DTEST(1,[DASUM(N,dx,INCx)],dtrue3(np1),dtrue3(np1),Dfac,Kprint)
          CASE (6)
            ! 31. SCASUM
            CALL STEST(1,[SCASUM(N,cx,INCx)],strue4(np1),strue4(np1),Sfac,Kprint)
          CASE (7)
            ! 32. SSCALE
            CALL SSCAL(N,sa,sx,INCx)
            DO i = 1, len
              strue(i) = REAL( dtrue5(i,np1,INCx), 4 )
            ENDDO
            CALL STEST(len,sx,strue,strue,Sfac,Kprint)
          CASE (8)
            ! 33. DSCALE
            CALL DSCAL(N,da,dx,INCx)
            CALL DTEST(len,dx,dtrue5(1,np1,INCx),dtrue5(1,np1,INCx),Dfac,Kprint)
          CASE (9)
            ! 34. CSCALE
            CALL CSCAL(N,ca,cx,INCx)
            CALL CTEST(len,cx,ctrue5(:,np1,INCx),ctrue5(:,np1,INCx),Sfac,Kprint)
          CASE (10)
            ! 35. CSSCAL
            CALL CSSCAL(N,sa,cx,INCx)
            CALL CTEST(len,cx,ctrue6(:,np1,INCx),ctrue6(:,np1,INCx),Sfac,Kprint)
          CASE (11)
            ! 36. ISAMAX
            CALL ITEST(1,[ISAMAX(N,sx,INCx)],itrue2(np1),Kprint)
          CASE (12)
            ! 37. IDAMAX
            CALL ITEST(1,[IDAMAX(N,dx,INCx)],itrue2(np1),Kprint)
          CASE (13)
            ! 38. ICAMAX
            CALL ITEST(1,[ICAMAX(N,cx,INCx)],itrue3(np1),Kprint)
          CASE DEFAULT
            ! 26. SNRM2
            stemp = REAL( dtrue1(np1), 4 )
            CALL STEST(1,[SNRM2(N,sx,INCx)],[stemp],[stemp],Sfac,Kprint)
        END SELECT
        !
      ENDDO
    ENDDO
  END SUBROUTINE CHECK1
  !** CHECK2
  SUBROUTINE CHECK2(Sfac,Sdfac,Dfac,Dqfac,Kprint)
    IMPLICIT NONE
    !>
    !***
    !  (UNKNOWN)
    !***
    ! **Library:**   SLATEC
    !***
    ! **Author:**  Lawson, C. L., (JPL)
    !***
    ! **Description:**
    !
    !     THIS SUBPROGRAM TESTS THE BASIC LINEAR ALGEBRA SUBPROGRAMS 1-11,
    !     14-15, AND 18-25. SUBPROGRAMS IN THIS SET EACH REQUIRE TWO ARRAYS
    !     IN THE PARAMETER LIST.
    !
    !     C. L. LAWSON, JPL, 1975 FEB 26, APR 29, MAY 8, MAY 28
    !
    !***
    ! **Routines called:**  CAXPY, CCOPY, CDOTC, CDOTU, CSWAP, DAXPY, DCOPY,
    !                    DDOT, DQDOTA, DQDOTI, DROT, DROTM, DSDOT, DSWAP,
    !                    DTEST, SAXPY, SCOPY, SDOT, SDSDOT, SROT, SROTM,
    !                    SSWAP, STEST
    !***
    ! COMMON BLOCKS    COMBLA

    !* REVISION HISTORY  (YYMMDD)
    !   750226  DATE WRITTEN
    !   890911  Removed unnecessary intrinsics.  (WRB)
    !   890911  REVISION DATE from Version 3.2
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    
    INTEGER i, ICAse, INCx, INCy, j, ki, kn, kni, kpar, Kprint, &
      ksize, lenx, leny, MODe, mx, my, N, NPRint
    REAL sa, sb, sc, Sdfac, SDOT, SDSDOT, Sfac, ss
    COMMON /COMBLA/ NPRint, ICAse, N, INCx, INCy, MODe, PASs
    !
    LOGICAL PASs
    INTEGER incxs(4), incys(4), lens(4,2), ns(4)
    REAL sx(7), sy(7), stx(7), sty(7), ssize1(4), ssize2(14,2)
    REAL ssize(7), qc(30), sparam(5), st7b(4,4), ssize3(4)
    REAL(8) :: dx(7), da, dx1(7), dy1(7), dy(7), dt7(4,4), &
      dt8(7,4,4)
    REAL(8) :: dx2(7), dy2(7), dt2(4,4,2), dparam(5), dpar(5,4)
    REAL(8) :: DSDOT, DDOT, DQDOTI, DQDOTA, Dfac, Dqfac
    REAL(8) :: dt10x(7,4,4), dt10y(7,4,4), db
    REAL(8) :: dsize1(4), dsize2(7,2), dsize(7)
    REAL(8) :: dc, ds, dt9x(7,4,4), dt9y(7,4,4), dtx(7), dty(7)
    REAL(8) :: dt19x(7,4,16), dt19xa(7,4,4), dt19xb(7,4,4)
    REAL(8) :: dt19xc(7,4,4), dt19xd(7,4,4), dt19y(7,4,16)
    REAL(8) :: dt19ya(7,4,4), dt19yb(7,4,4), dt19yc(7,4,4)
    REAL(8) :: dt19yd(7,4,4)
    !
    EQUIVALENCE (dt19x(1,1,1),dt19xa(1,1,1))
    EQUIVALENCE (dt19x(1,1,5),dt19xb(1,1,1))
    EQUIVALENCE (dt19x(1,1,9),dt19xc(1,1,1))
    EQUIVALENCE (dt19x(1,1,13),dt19xd(1,1,1))
    EQUIVALENCE (dt19y(1,1,1),dt19ya(1,1,1))
    EQUIVALENCE (dt19y(1,1,5),dt19yb(1,1,1))
    EQUIVALENCE (dt19y(1,1,9),dt19yc(1,1,1))
    EQUIVALENCE (dt19y(1,1,13),dt19yd(1,1,1))
    COMPLEX cx(7), ca, cx1(7), cy1(7), cy(7), ct6(4,4), ct7(4,4)
    COMPLEX ct8(7,4,4), csize1(4), csize2(7,2)
    COMPLEX ct10x(7,4,4), ct10y(7,4,4)
    COMPLEX CDOTC, CDOTU
    DATA sa, da, ca, db, sb/.3, .3D0, (.4,-.7), .25D0, .1/
    DATA incxs/1, 2, -2, -1/
    DATA incys/1, -2, 1, -2/
    DATA lens/1, 1, 2, 4, 1, 1, 3, 7/
    DATA ns/0, 1, 2, 4/
    DATA sc, ss, dc, ds/.8, .6, .8D0, .6D0/
    DATA dx1/.6D0, .1D0, -.5D0, .8D0, .9D0, -.3D0, -.4D0/
    DATA dy1/.5D0, -.9D0, .3D0, .7D0, -.6D0, .2D0, .8D0/
    DATA dx2/1.D0, .01D0, .02D0, 1.D0, .06D0, 2.D0, 1.D0/
    DATA dy2/1.D0, .04D0, -.03D0, -1.D0, .05D0, 3.D0, -1.D0/
    !            THE TERMS D11(3,2) AND D11(4,2) WILL BE SET BY
    !            COMPUTATION AT RUN TIME.
    DATA cx1/(.7,-.8), (-.4,-.7), (-.1,-.9), (.2,-.8), (-.9,-.4), (.1,.4), (-.6,.6)/
    DATA cy1/(.6,-.6), (-.9,.5), (.7,-.6), (.1,-.5), (-.1,-.2), (-.5,-.3), (.8,-.7)/
    !
    !                             FOR DQDOTI AND DQDOTA
    !
    DATA dt2/0.25D0, 1.25D0, 1.2504D0, 0.2498D0, 0.25D0, 1.25D0, &
      0.24D0, 0.2492D0, 0.25D0, 1.25D0, 0.31D0, 0.2518D0, 0.25D0, &
      1.25D0, 1.2497D0, 0.2507D0, 0.D0, 2.D0, 2.0008D0, -.0004D0, &
      0.D0, 2.D0, -.02D0, -.0016D0, 0.D0, 2.D0, .12D0, .0036D0, &
      0.D0, 2.D0, 1.9994D0, .0014D0/
    DATA dt7/0.D0, .30D0, .21D0, .62D0, 0.D0, .30D0, -.07D0, .85D0, &
      0.D0, .30D0, -.79D0, -.74D0, 0.D0, .30D0, .33D0, 1.27D0/
    DATA st7b/.1, .4, .31, .72, .1, .4, .03, .95, .1, .4, -.69, &
      -.64, .1, .4, .43, 1.37/
    !
    !                       FOR CDOTU
    !
    DATA ct7/(0.,0.), (-.06,-.90), (.65,-.47), (-.34,-1.22), (0.,0.), &
      (-.06,-.90), (-.59,-1.46), (-1.04,-.04), (0.,0.), (-.06,-.90), &
      (-.83,.59), (.07,-.37), (0.,0.), (-.06,-.90), (-.76,-1.15), &
      (-1.33,-1.82)/
    !
    !                       FOR CDOTC
    !
    DATA ct6/(0.,0.), (.90,0.06), (.91,-.77), (1.80,-.10), (0.,0.), &
      (.90,0.06), (1.45,.74), (.20,.90), (0.,0.), (.90,0.06), &
      (-.55,.23), (.83,-.39), (0.,0.), (.90,0.06), (1.04,0.79), &
      (1.95,1.22)/
    !
    DATA dt8/.5D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, .68D0, 0.D0, &
      0.D0, 0.D0, 0.D0, 0.D0, 0.D0, .68D0, -.87D0, 0.D0, 0.D0, &
      0.D0, 0.D0, 0.D0, .68D0, -.87D0, .15D0, .94D0, 0.D0, 0.D0, &
      0.D0, .5D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, .68D0, &
      0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, .35D0, -.9D0, .48D0, &
      0.D0, 0.D0, 0.D0, 0.D0, .38D0, -.9D0, .57D0, .7D0, -.75D0, &
      .2D0, .98D0, .5D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      .68D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, .35D0, -.72D0, &
      0.D0, 0.D0, 0.D0, 0.D0, 0.D0, .38D0, -.63D0, .15D0, .88D0, &
      0.D0, 0.D0, 0.D0, .5D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      .68D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, .68D0, -.9D0, &
      .33D0, 0.D0, 0.D0, 0.D0, 0.D0, .68D0, -.9D0, .33D0, .7D0, &
      -.75D0, .2D0, 1.04D0/
    !
    DATA ct8/(.6,-.6), (0.,0.), (0.,0.), (0.,0.), (0.,0.), (0.,0.), &
      (0.,0.), (.32,-1.41), (0.,0.), (0.,0.), (0.,0.), (0.,0.), &
      (0.,0.), (0.,0.), (.32,-1.41), (-1.55,.5), (0.,0.), (0.,0.), &
      (0.,0.), (0.,0.), (0.,0.), (.32,-1.41), (-1.55,.5), (.03,-.89), &
      (-.38,-.96), (0.,0.), (0.,0.), (0.,0.), (.6,-.6), (0.,0.), &
      (0.,0.), (0.,0.), (0.,0.), (0.,0.), (0.,0.), (.32,-1.41), &
      (0.,0.), (0.,0.), (0.,0.), (0.,0.), (0.,0.), (0.,0.), &
      (-.07,-.89), (-.9,.5), (.42,-1.41), (0.,0.), (0.,0.), (0.,0.), &
      (0.,0.), (.78,.06), (-.9,.5), (.06,-.13), (.1,-.5), (-.77,-.49)&
      , (-.5,-.3), (.52,-1.51), (.6,-.6), (0.,0.), (0.,0.), (0.,0.), &
      (0.,0.), (0.,0.), (0.,0.), (.32,-1.41), (0.,0.), (0.,0.), &
      (0.,0.), (0.,0.), (0.,0.), (0.,0.), (-.07,-.89), (-1.18,-.31), &
      (0.,0.), (0.,0.), (0.,0.), (0.,0.), (0.,0.), (.78,.06), &
      (-1.54,.97), (.03,-.89), (-.18,-1.31), (0.,0.), (0.,0.), (0.,0.)&
      , (.6,-.6), (0.,0.), (0.,0.), (0.,0.), (0.,0.), (0.,0.), &
      (0.,0.), (.32,-1.41), (0.,0.), (0.,0.), (0.,0.), (0.,0.), &
      (0.,0.), (0.,0.), (.32,-1.41), (-.9,.5), (.05,-.6), (0.,0.), &
      (0.,0.), (0.,0.), (0.,0.), (.32,-1.41), (-.9,.5), (.05,-.6), &
      (.1,-.5), (-.77,-.49), (-.5,-.3), (.32,-1.16)/
    !
    !
    !                TRUE X VALUES AFTER ROTATION USING SROT OR DROT.
    DATA dt9x/.6D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, .78D0, 0.D0, &
      0.D0, 0.D0, 0.D0, 0.D0, 0.D0, .78D0, -.46D0, 0.D0, 0.D0, &
      0.D0, 0.D0, 0.D0, .78D0, -.46D0, -.22D0, 1.06D0, 0.D0, 0.D0, &
      0.D0, .6D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, .78D0, &
      0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, .66D0, .1D0, -.1D0, &
      0.D0, 0.D0, 0.D0, 0.D0, .96D0, .1D0, -.76D0, .8D0, .90D0, &
      -.3D0, -.02D0, .6D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      .78D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, -.06D0, .1D0, &
      -.1D0, 0.D0, 0.D0, 0.D0, 0.D0, .90D0, .1D0, -.22D0, .8D0, &
      .18D0, -.3D0, -.02D0, .6D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      0.D0, .78D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, .78D0, &
      .26D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, .78D0, .26D0, -.76D0, &
      1.12D0, 0.D0, 0.D0, 0.D0/
    !
    !                TRUE Y VALUES AFTER ROTATION USING SROT OR DROT.
    !
    DATA dt9y/.5D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, .04D0, 0.D0, &
      0.D0, 0.D0, 0.D0, 0.D0, 0.D0, .04D0, -.78D0, 0.D0, 0.D0, &
      0.D0, 0.D0, 0.D0, .04D0, -.78D0, .54D0, .08D0, 0.D0, 0.D0, &
      0.D0, .5D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, .04D0, &
      0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, .7D0, -.9D0, -.12D0, &
      0.D0, 0.D0, 0.D0, 0.D0, .64D0, -.9D0, -.30D0, .7D0, -.18D0, &
      .2D0, .28D0, .5D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      .04D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, .7D0, -1.08D0, &
      0.D0, 0.D0, 0.D0, 0.D0, 0.D0, .64D0, -1.26D0, .54D0, .20D0, &
      0.D0, 0.D0, 0.D0, .5D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      .04D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, .04D0, -.9D0, &
      .18D0, 0.D0, 0.D0, 0.D0, 0.D0, .04D0, -.9D0, .18D0, .7D0, &
      -.18D0, .2D0, .16D0/
    !
    DATA dt10x/.6D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, .5D0, 0.D0, &
      0.D0, 0.D0, 0.D0, 0.D0, 0.D0, .5D0, -.9D0, 0.D0, 0.D0, &
      0.D0, 0.D0, 0.D0, .5D0, -.9D0, .3D0, .7D0, 0.D0, 0.D0, &
      0.D0, .6D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, .5D0, 0.D0, &
      0.D0, 0.D0, 0.D0, 0.D0, 0.D0, .3D0, .1D0, .5D0, 0.D0, 0.D0, &
      0.D0, 0.D0, .8D0, .1D0, -.6D0, .8D0, .3D0, -.3D0, .5D0, &
      .6D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, .5D0, 0.D0, 0.D0, &
      0.D0, 0.D0, 0.D0, 0.D0, -.9D0, .1D0, .5D0, 0.D0, 0.D0, &
      0.D0, 0.D0, .7D0, .1D0, .3D0, .8D0, -.9D0, -.3D0, .5D0, &
      .6D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, .5D0, 0.D0, 0.D0, &
      0.D0, 0.D0, 0.D0, 0.D0, .5D0, .3D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      0.D0, .5D0, .3D0, -.6D0, .8D0, 0.D0, 0.D0, 0.D0/
    !
    DATA dt10y/.5D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, .6D0, 0.D0, &
      0.D0, 0.D0, 0.D0, 0.D0, 0.D0, .6D0, .1D0, 0.D0, 0.D0, 0.D0, &
      0.D0, 0.D0, .6D0, .1D0, -.5D0, .8D0, 0.D0, 0.D0, 0.D0, &
      .5D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, .6D0, 0.D0, 0.D0, &
      0.D0, 0.D0, 0.D0, 0.D0, -.5D0, -.9D0, .6D0, 0.D0, 0.D0, &
      0.D0, 0.D0, -.4D0, -.9D0, .9D0, .7D0, -.5D0, .2D0, .6D0, &
      .5D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, .6D0, 0.D0, 0.D0, &
      0.D0, 0.D0, 0.D0, 0.D0, -.5D0, .6D0, 0.D0, 0.D0, 0.D0, &
      0.D0, 0.D0, -.4D0, .9D0, -.5D0, .6D0, 0.D0, 0.D0, 0.D0, &
      .5D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, .6D0, 0.D0, 0.D0, &
      0.D0, 0.D0, 0.D0, 0.D0, .6D0, -.9D0, .1D0, 0.D0, 0.D0, &
      0.D0, 0.D0, .6D0, -.9D0, .1D0, .7D0, -.5D0, .2D0, .8D0/
    !
    DATA ct10x/(.7,-.8), (0.,0.), (0.,0.), (0.,0.), (0.,0.), (0.,0.), &
      (0.,0.), (.6,-.6), (0.,0.), (0.,0.), (0.,0.), (0.,0.), (0.,0.)&
      , (0.,0.), (.6,-.6), (-.9,.5), (0.,0.), (0.,0.), (0.,0.), &
      (0.,0.), (0.,0.), (.6,-.6), (-.9,.5), (.7,-.6), (.1,-.5), &
      (0.,0.), (0.,0.), (0.,0.), (.7,-.8), (0.,0.), (0.,0.), (0.,0.)&
      , (0.,0.), (0.,0.), (0.,0.), (.6,-.6), (0.,0.), (0.,0.), &
      (0.,0.), (0.,0.), (0.,0.), (0.,0.), (.7,-.6), (-.4,-.7), &
      (.6,-.6), (0.,0.), (0.,0.), (0.,0.), (0.,0.), (.8,-.7), &
      (-.4,-.7), (-.1,-.2), (.2,-.8), (.7,-.6), (.1,.4), (.6,-.6), &
      (.7,-.8), (0.,0.), (0.,0.), (0.,0.), (0.,0.), (0.,0.), (0.,0.)&
      , (.6,-.6), (0.,0.), (0.,0.), (0.,0.), (0.,0.), (0.,0.), &
      (0.,0.), (-.9,.5), (-.4,-.7), (.6,-.6), (0.,0.), (0.,0.), &
      (0.,0.), (0.,0.), (.1,-.5), (-.4,-.7), (.7,-.6), (.2,-.8), &
      (-.9,.5), (.1,.4), (.6,-.6), (.7,-.8), (0.,0.), (0.,0.), &
      (0.,0.), (0.,0.), (0.,0.), (0.,0.), (.6,-.6), (0.,0.), (0.,0.)&
      , (0.,0.), (0.,0.), (0.,0.), (0.,0.), (.6,-.6), (.7,-.6), &
      (0.,0.), (0.,0.), (0.,0.), (0.,0.), (0.,0.), (.6,-.6), (.7,-.6)&
      , (-.1,-.2), (.8,-.7), (0.,0.), (0.,0.), (0.,0.)/
    !
    DATA ct10y/(.6,-.6), (0.,0.), (0.,0.), (0.,0.), (0.,0.), (0.,0.), &
      (0.,0.), (.7,-.8), (0.,0.), (0.,0.), (0.,0.), (0.,0.), (0.,0.)&
      , (0.,0.), (.7,-.8), (-.4,-.7), (0.,0.), (0.,0.), (0.,0.), &
      (0.,0.), (0.,0.), (.7,-.8), (-.4,-.7), (-.1,-.9), (.2,-.8), &
      (0.,0.), (0.,0.), (0.,0.), (.6,-.6), (0.,0.), (0.,0.), (0.,0.)&
      , (0.,0.), (0.,0.), (0.,0.), (.7,-.8), (0.,0.), (0.,0.), &
      (0.,0.), (0.,0.), (0.,0.), (0.,0.), (-.1,-.9), (-.9,.5), &
      (.7,-.8), (0.,0.), (0.,0.), (0.,0.), (0.,0.), (-.6,.6), &
      (-.9,.5), (-.9,-.4), (.1,-.5), (-.1,-.9), (-.5,-.3), (.7,-.8), &
      (.6,-.6), (0.,0.), (0.,0.), (0.,0.), (0.,0.), (0.,0.), (0.,0.)&
      , (.7,-.8), (0.,0.), (0.,0.), (0.,0.), (0.,0.), (0.,0.), &
      (0.,0.), (-.1,-.9), (.7,-.8), (0.,0.), (0.,0.), (0.,0.), &
      (0.,0.), (0.,0.), (-.6,.6), (-.9,-.4), (-.1,-.9), (.7,-.8), &
      (0.,0.), (0.,0.), (0.,0.), (.6,-.6), (0.,0.), (0.,0.), (0.,0.)&
      , (0.,0.), (0.,0.), (0.,0.), (.7,-.8), (0.,0.), (0.,0.), &
      (0.,0.), (0.,0.), (0.,0.), (0.,0.), (.7,-.8), (-.9,.5), &
      (-.4,-.7), (0.,0.), (0.,0.), (0.,0.), (0.,0.), (.7,-.8), &
      (-.9,.5), (-.4,-.7), (.1,-.5), (-.1,-.9), (-.5,-.3), (.2,-.8)/
    !                        TRUE X RESULTS F0R ROTATIONS SROTM AND DROTM
    DATA dt19xa/.6D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, .6D0, 0.D0, &
      0.D0, 0.D0, 0.D0, 0.D0, 0.D0, .6D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      0.D0, 0.D0, .6D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, .6D0, &
      0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, -.8D0, 0.D0, 0.D0, &
      0.D0, 0.D0, 0.D0, 0.D0, -.9D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      0.D0, 0.D0, 3.5D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      .6D0, .1D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, -.8D0, 3.8D0, &
      0.D0, 0.D0, 0.D0, 0.D0, 0.D0, -.9D0, 2.8D0, 0.D0, 0.D0, &
      0.D0, 0.D0, 0.D0, 3.5D0, -.4D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      0.D0, .6D0, .1D0, -.5D0, .8D0, 0.D0, 0.D0, 0.D0, -.8D0, &
      3.8D0, -2.2D0, -1.2D0, 0.D0, 0.D0, 0.D0, -.9D0, 2.8D0, &
      -1.4D0, -1.3D0, 0.D0, 0.D0, 0.D0, 3.5D0, -.4D0, -2.2D0, &
      4.7D0, 0.D0, 0.D0, 0.D0/
    !
    DATA dt19xb/.6D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, .6D0, 0.D0, &
      0.D0, 0.D0, 0.D0, 0.D0, 0.D0, .6D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      0.D0, 0.D0, .6D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, .6D0, &
      0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, -.8D0, 0.D0, 0.D0, &
      0.D0, 0.D0, 0.D0, 0.D0, -.9D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      0.D0, 0.D0, 3.5D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      .6D0, .1D0, -.5D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, .1D0, &
      -3.0D0, 0.D0, 0.D0, 0.D0, 0.D0, -.3D0, .1D0, -2.0D0, 0.D0, &
      0.D0, 0.D0, 0.D0, 3.3D0, .1D0, -2.0D0, 0.D0, 0.D0, 0.D0, &
      0.D0, .6D0, .1D0, -.5D0, .8D0, .9D0, -.3D0, -.4D0, -2.0D0, &
      .1D0, 1.4D0, .8D0, .6D0, -.3D0, -2.8D0, -1.8D0, .1D0, 1.3D0, &
      .8D0, 0.D0, -.3D0, -1.9D0, 3.8D0, .1D0, -3.1D0, .8D0, 4.8D0, &
      -.3D0, -1.5D0/
    !
    DATA dt19xc/.6D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, .6D0, 0.D0, &
      0.D0, 0.D0, 0.D0, 0.D0, 0.D0, .6D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      0.D0, 0.D0, .6D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, .6D0, &
      0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, -.8D0, 0.D0, 0.D0, &
      0.D0, 0.D0, 0.D0, 0.D0, -.9D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      0.D0, 0.D0, 3.5D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      .6D0, .1D0, -.5D0, 0.D0, 0.D0, 0.D0, 0.D0, 4.8D0, .1D0, &
      -3.0D0, 0.D0, 0.D0, 0.D0, 0.D0, 3.3D0, .1D0, -2.0D0, 0.D0, &
      0.D0, 0.D0, 0.D0, 2.1D0, .1D0, -2.0D0, 0.D0, 0.D0, 0.D0, &
      0.D0, .6D0, .1D0, -.5D0, .8D0, .9D0, -.3D0, -.4D0, -1.6D0, &
      .1D0, -2.2D0, .8D0, 5.4D0, -.3D0, -2.8D0, -1.5D0, .1D0, &
      -1.4D0, .8D0, 3.6D0, -.3D0, -1.9D0, 3.7D0, .1D0, -2.2D0, &
      .8D0, 3.6D0, -.3D0, -1.5D0/
    !
    DATA dt19xd/.6D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, .6D0, 0.D0, &
      0.D0, 0.D0, 0.D0, 0.D0, 0.D0, .6D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      0.D0, 0.D0, .6D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, .6D0, &
      0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, -.8D0, 0.D0, 0.D0, &
      0.D0, 0.D0, 0.D0, 0.D0, -.9D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      0.D0, 0.D0, 3.5D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      .6D0, .1D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, -.8D0, -1.0D0, &
      0.D0, 0.D0, 0.D0, 0.D0, 0.D0, -.9D0, -.8D0, 0.D0, 0.D0, &
      0.D0, 0.D0, 0.D0, 3.5D0, .8D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      0.D0, .6D0, .1D0, -.5D0, .8D0, 0.D0, 0.D0, 0.D0, -.8D0, &
      -1.0D0, 1.4D0, -1.6D0, 0.D0, 0.D0, 0.D0, -.9D0, -.8D0, &
      1.3D0, -1.6D0, 0.D0, 0.D0, 0.D0, 3.5D0, .8D0, -3.1D0, 4.8D0, &
      0.D0, 0.D0, 0.D0/
    !                        TRUE Y RESULTS FOR ROTATIONS SROTM AND DROTM
    DATA dt19ya/.5D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, .5D0, 0.D0, &
      0.D0, 0.D0, 0.D0, 0.D0, 0.D0, .5D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      0.D0, 0.D0, .5D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, .5D0, &
      0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, .7D0, 0.D0, 0.D0, 0.D0, &
      0.D0, 0.D0, 0.D0, 1.7D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      0.D0, -2.6D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, .5D0, &
      -.9D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, .7D0, -4.8D0, 0.D0, &
      0.D0, 0.D0, 0.D0, 0.D0, 1.7D0, -.7D0, 0.D0, 0.D0, 0.D0, &
      0.D0, 0.D0, -2.6D0, 3.5D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      .5D0, -.9D0, .3D0, .7D0, 0.D0, 0.D0, 0.D0, .7D0, -4.8D0, &
      3.0D0, 1.1D0, 0.D0, 0.D0, 0.D0, 1.7D0, -.7D0, -.7D0, 2.3D0, &
      0.D0, 0.D0, 0.D0, -2.6D0, 3.5D0, -.7D0, -3.6D0, 0.D0, 0.D0, &
      0.D0/
    !
    DATA dt19yb/.5D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, .5D0, 0.D0, &
      0.D0, 0.D0, 0.D0, 0.D0, 0.D0, .5D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      0.D0, 0.D0, .5D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, .5D0, &
      0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, .7D0, 0.D0, 0.D0, 0.D0, &
      0.D0, 0.D0, 0.D0, 1.7D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      0.D0, -2.6D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, .5D0, &
      -.9D0, .3D0, 0.D0, 0.D0, 0.D0, 0.D0, 4.0D0, -.9D0, -.3D0, &
      0.D0, 0.D0, 0.D0, 0.D0, -.5D0, -.9D0, 1.5D0, 0.D0, 0.D0, &
      0.D0, 0.D0, -1.5D0, -.9D0, -1.8D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      .5D0, -.9D0, .3D0, .7D0, -.6D0, .2D0, .8D0, 3.7D0, -.9D0, &
      -1.2D0, .7D0, -1.5D0, .2D0, 2.2D0, -.3D0, -.9D0, 2.1D0, &
      .7D0, -1.6D0, .2D0, 2.0D0, -1.6D0, -.9D0, -2.1D0, .7D0, &
      2.9D0, .2D0, -3.8D0/
    !
    DATA dt19yc/.5D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, .5D0, 0.D0, &
      0.D0, 0.D0, 0.D0, 0.D0, 0.D0, .5D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      0.D0, 0.D0, .5D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, .5D0, &
      0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, .7D0, 0.D0, 0.D0, 0.D0, &
      0.D0, 0.D0, 0.D0, 1.7D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      0.D0, -2.6D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, .5D0, &
      -.9D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 4.0D0, -6.3D0, 0.D0, &
      0.D0, 0.D0, 0.D0, 0.D0, -.5D0, .3D0, 0.D0, 0.D0, 0.D0, &
      0.D0, 0.D0, -1.5D0, 3.0D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      .5D0, -.9D0, .3D0, .7D0, 0.D0, 0.D0, 0.D0, 3.7D0, -7.2D0, &
      3.0D0, 1.7D0, 0.D0, 0.D0, 0.D0, -.3D0, .9D0, -.7D0, 1.9D0, &
      0.D0, 0.D0, 0.D0, -1.6D0, 2.7D0, -.7D0, -3.4D0, 0.D0, 0.D0, &
      0.D0/
    !
    DATA dt19yd/.5D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, .5D0, 0.D0, &
      0.D0, 0.D0, 0.D0, 0.D0, 0.D0, .5D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      0.D0, 0.D0, .5D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, .5D0, &
      0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, .7D0, 0.D0, 0.D0, 0.D0, &
      0.D0, 0.D0, 0.D0, 1.7D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      0.D0, -2.6D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, .5D0, &
      -.9D0, .3D0, 0.D0, 0.D0, 0.D0, 0.D0, .7D0, -.9D0, 1.2D0, &
      0.D0, 0.D0, 0.D0, 0.D0, 1.7D0, -.9D0, .5D0, 0.D0, 0.D0, &
      0.D0, 0.D0, -2.6D0, -.9D0, -1.3D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      .5D0, -.9D0, .3D0, .7D0, -.6D0, .2D0, .8D0, .7D0, -.9D0, &
      1.2D0, .7D0, -1.5D0, .2D0, 1.6D0, 1.7D0, -.9D0, .5D0, .7D0, &
      -1.6D0, .2D0, 2.4D0, -2.6D0, -.9D0, -1.3D0, .7D0, 2.9D0, &
      .2D0, -4.0D0/
    !
    DATA ssize1/0., .3, 1.6, 3.2/
    DATA dsize1/0.D0, .3D0, 1.6D0, 3.2D0/
    DATA ssize3/.1, .4, 1.7, 3.3/
    !
    !                         FOR CDOTC AND CDOTU
    !
    DATA csize1/(0.,0.), (.9,.9), (1.63,1.73), (2.90,2.78)/
    DATA ssize2/0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., &
      0., 0., 1.17, 1.17, 1.17, 1.17, 1.17, 1.17, 1.17, 1.17, &
      1.17, 1.17, 1.17, 1.17, 1.17, 1.17/
    DATA dsize2/0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 1.17D0, &
      1.17D0, 1.17D0, 1.17D0, 1.17D0, 1.17D0, 1.17D0/
    !
    !                         FOR CAXPY
    !
    DATA csize2/(0.,0.), (0.,0.), (0.,0.), (0.,0.), (0.,0.), (0.,0.), &
      (0.,0.), (1.54,1.54), (1.54,1.54), (1.54,1.54), (1.54,1.54), &
      (1.54,1.54), (1.54,1.54), (1.54,1.54)/
    !
    !                         FOR SROTM AND DROTM
    !
    DATA dpar/ - 2.D0, 0.D0, 0.D0, 0.D0, 0.D0, -1.D0, 2.D0, -3.D0, &
      -4.D0, 5.D0, 0.D0, 0.D0, 2.D0, -3.D0, 0.D0, 1.D0, 5.D0, &
      2.D0, 0.D0, -4.D0/
    !* FIRST EXECUTABLE STATEMENT  CHECK2
    DO ki = 1, 4
      INCx = incxs(ki)
      INCy = incys(ki)
      mx = ABS(INCx)
      my = ABS(INCy)
      !
      DO kn = 1, 4
        N = ns(kn)
        ksize = MIN(2,kn)
        lenx = lens(kn,mx)
        leny = lens(kn,my)
        ! INITIALIZE ALL ARGUMENT ARRAYS.
        DO i = 1, 7
          sx(i) = REAL( dx1(i), 4 )
          sy(i) = REAL( dy1(i), 4 )
          dx(i) = dx1(i)
          dy(i) = dy1(i)
          cx(i) = cx1(i)
          cy(i) = cy1(i)
        ENDDO
        !
        ! BRANCH TO SELECT SUBPROGRAM TO BE TESTED.
        !
        SELECT CASE (ICAse)
          CASE (2)
            ! 2. DSDOT
            CALL STEST(1,[REAL(DSDOT(N,sx,INCx,sy,INCy))],[REAL(dt7(kn,ki))],&
              ssize1(kn),Sfac,Kprint)
          CASE (3)
            ! 3. SDSDOT
            CALL STEST(1,[SDSDOT(N,sb,sx,INCx,sy,INCy)],st7b(kn,ki),ssize3(kn),&
              Sfac,Kprint)
          CASE (4)
            ! 4. DDOT
            CALL DTEST(1,[DDOT(N,dx,INCx,dy,INCy)],dt7(kn,ki),dsize1(kn),Dfac,&
              Kprint)
          CASE (5)
            ! 5. DQDOTI
            !     DQDOTI AND DQDOTA ARE SUPPOSED TO USE EXTENDED
            !     PRECISION ARITHMETIC INTERNALLY.
            !     SET MODE = 1 OR 2 TO DISTINGUISH TESTS OF DQDOTI OR DQDOTA
            !     IN THE DIAGNOSTIC OUTPUT.
            !
            MODe = 1
            CALL DTEST(1,[DQDOTI(N,db,qc,dx2,INCx,dy2,INCy)],dt2(kn,ki,1),&
              dt2(kn,ki,1),Dqfac,Kprint)
          CASE (6)
            ! 6. DQDOTA
            !     TO TEST DQDOTA WE ACTUALLY TEST BOTH DQDOTI AND DQDOTA.
            !     THE OUTPUT VALUE OF QX FROM DQDOTI WILL BE USED AS INPUT
            !     TO DQDOTA.  QX IS SUPPOSED TO BE IN A MACHINE-DEPENDENT
            !     EXTENDED PRECISION FORM.
            !     MODE IS SET TO 1 OR 2 TO DISTINGUISH TESTS OF
            !     DQDOTI OR DQDOTA IN THE DIAGNOSTIC OUTPUT.
            !
            MODe = 1
            CALL DTEST(1,[DQDOTI(N,db,qc,dx2,INCx,dy2,INCy)],dt2(kn,ki,1),&
              dt2(kn,ki,1),Dqfac,Kprint)
            MODe = 2
            CALL DTEST(1,[DQDOTA(N,-db,qc,dx2,INCx,dy2,INCy)],dt2(kn,ki,2),&
              dt2(kn,ki,2),Dqfac,Kprint)
          CASE (7)
            ! 7. CDOTC
            CALL CTEST(1,[CDOTC(N,cx,INCx,cy,INCy)],ct6(kn,ki),csize1(kn),Sfac,Kprint)
          CASE (8)
            ! 8. CDOTU
            CALL CTEST(1,[CDOTU(N,cx,INCx,cy,INCy)],ct7(kn,ki),csize1(kn),Sfac,Kprint)
          CASE (9)
            ! 9. SAXPY
            CALL SAXPY(N,sa,sx,INCx,sy,INCy)
            DO j = 1, leny
              sty(j) = REAL( dt8(j,kn,ki), 4 )
            ENDDO
            CALL STEST(leny,sy,sty,ssize2(1,ksize),Sfac,Kprint)
          CASE (10)
            ! 10. DAXPY
            CALL DAXPY(N,da,dx,INCx,dy,INCy)
            CALL DTEST(leny,dy,dt8(1,kn,ki),dsize2(1,ksize),Dfac,Kprint)
          CASE (11)
            ! 11. CAXPY
            CALL CAXPY(N,ca,cx,INCx,cy,INCy)
            CALL CTEST(leny,cy,ct8(:,kn,ki),csize2(:,ksize),Sfac,Kprint)
          CASE (12,13,16,17)
            GOTO 100
          CASE (14)
            ! 14. SROT
            DO i = 1, 7
              sx(i) = REAL( dx1(i), 4 )
              sy(i) = REAL( dy1(i), 4 )
              stx(i) = REAL( dt9x(i,kn,ki), 4 )
              sty(i) = REAL( dt9y(i,kn,ki), 4 )
            ENDDO
            CALL SROT(N,sx,INCx,sy,INCy,sc,ss)
            CALL STEST(lenx,sx,stx,ssize2(1,ksize),Sfac,Kprint)
            CALL STEST(leny,sy,sty,ssize2(1,ksize),Sfac,Kprint)
          CASE (15)
            ! 15. DROT
            DO i = 1, 7
              dx(i) = dx1(i)
              dy(i) = dy1(i)
            ENDDO
            CALL DROT(N,dx,INCx,dy,INCy,dc,ds)
            CALL DTEST(lenx,dx,dt9x(1,kn,ki),dsize2(1,ksize),Dfac,Kprint)
            CALL DTEST(leny,dy,dt9y(1,kn,ki),dsize2(1,ksize),Dfac,Kprint)
          CASE (18)
            ! 18. SROTM
            kni = kn + 4*(ki-1)
            DO kpar = 1, 4
              DO i = 1, 7
                sx(i) = REAL( dx1(i), 4 )
                sy(i) = REAL( dy1(i), 4 )
                stx(i) = REAL( dt19x(i,kpar,kni), 4 )
                sty(i) = REAL( dt19y(i,kpar,kni), 4 )
              ENDDO
              !
              DO i = 1, 5
                sparam(i) = REAL( dpar(i,kpar), 4 )
              ENDDO
              ! SET MODE TO IDENTIFY DIAGNOSTIC OUTPUT, IF ANY
              MODe = INT(sparam(1))
              !
              DO i = 1, lenx
                ssize(i) = stx(i)
              ENDDO
              !  THE TRUE RESULTS DT19X(1,2,7) AND
              !  DT19X(5,3,8) ARE ZERO DUE TO CANCELLATION.
              !  DT19X(1,2,7) = 2.*.6 - 4.*.3 = 0
              !  DT19X(5,3,8) = .9 - 3.*.3 = 0
              !  FOR THESE CASES RESPECTIVELY SET SIZE( )
              !  EQUAL TO 2.4 AND 1.8
              IF ( (kpar==2).AND.(kni==7) ) ssize(1) = 2.4E0
              IF ( (kpar==3).AND.(kni==8) ) ssize(5) = 1.8E0
              !
              CALL SROTM(N,sx,INCx,sy,INCy,sparam)
              CALL STEST(lenx,sx,stx,ssize,Sfac,Kprint)
              CALL STEST(leny,sy,sty,sty,Sfac,Kprint)
            ENDDO
          CASE (19)
            ! 19. DROTM
            kni = kn + 4*(ki-1)
            DO kpar = 1, 4
              DO i = 1, 7
                dx(i) = dx1(i)
                dy(i) = dy1(i)
                dtx(i) = dt19x(i,kpar,kni)
                dty(i) = dt19y(i,kpar,kni)
              ENDDO
              !
              DO i = 1, 5
                dparam(i) = dpar(i,kpar)
              ENDDO
              ! SET MODE TO IDENTIFY DIAGNOSTIC OUTPUT, IF ANY
              MODe = INT(dparam(1))
              !
              DO i = 1, lenx
                dsize(i) = dtx(i)
              ENDDO
              !  SEE REMARK ABOVE ABOUT DT11X(1,2,7) AND DT11X(5,3,8).
              IF ( (kpar==2).AND.(kni==7) ) dsize(1) = 2.4D0
              IF ( (kpar==3).AND.(kni==8) ) dsize(5) = 1.8D0
              !
              CALL DROTM(N,dx,INCx,dy,INCy,dparam)
              CALL DTEST(lenx,dx,dtx,dsize,Dfac,Kprint)
              CALL DTEST(leny,dy,dty,dty,Dfac,Kprint)
            ENDDO
          CASE (20)
            ! 20. SCOPY
            DO i = 1, 7
              sty(i) = REAL( dt10y(i,kn,ki), 4 )
            ENDDO
            CALL SCOPY(N,sx,INCx,sy,INCy)
            CALL STEST(leny,sy,sty,ssize2(1,1),1.,Kprint)
          CASE (21)
            ! 21. DCOPY
            CALL DCOPY(N,dx,INCx,dy,INCy)
            CALL DTEST(leny,dy,dt10y(1,kn,ki),dsize2(1,1),1.D0,Kprint)
          CASE (22)
            ! 22. CCOPY
            CALL CCOPY(N,cx,INCx,cy,INCy)
            CALL CTEST(leny,cy,ct10y(:,kn,ki),CMPLX(ssize2(1:8:2,1),ssize2(2:8:2,1)),1.,Kprint)
          CASE (23)
            ! 23. SSWAP
            CALL SSWAP(N,sx,INCx,sy,INCy)
            DO i = 1, 7
              stx(i) = REAL( dt10x(i,kn,ki), 4 )
              sty(i) = REAL( dt10y(i,kn,ki), 4 )
            ENDDO
            CALL STEST(lenx,sx,stx,ssize2(1,1),1.,Kprint)
            CALL STEST(leny,sy,sty,ssize2(1,1),1.,Kprint)
          CASE (24)
            ! 24. DSWAP
            CALL DSWAP(N,dx,INCx,dy,INCy)
            CALL DTEST(lenx,dx,dt10x(1,kn,ki),dsize2(1,1),1.D0,Kprint)
            CALL DTEST(leny,dy,dt10y(1,kn,ki),dsize2(1,1),1.D0,Kprint)
          CASE (25)
            ! 25. CSWAP
            CALL CSWAP(N,cx,INCx,cy,INCy)
            CALL CTEST(lenx,cx,ct10x(:,kn,ki),CMPLX(ssize2(1:8:2,1),ssize2(2:8:2,1)),1.,Kprint)
            CALL CTEST(leny,cy,ct10y(:,kn,ki),CMPLX(ssize2(1:8:2,1),ssize2(2:8:2,1)),1.,Kprint)
          CASE DEFAULT
            !                                                              1. SDOT
            CALL STEST(1,[SDOT(N,sx,INCx,sy,INCy)],[REAL(dt7(kn,ki))],ssize1(kn),&
              Sfac,Kprint)
        END SELECT
        !
        !
        !
      ENDDO
    ENDDO
    RETURN
    !                 THE FOLLOWING STOP SHOULD NEVER BE REACHED.
    100  STOP
  END SUBROUTINE CHECK2
  !** ITEST
  SUBROUTINE ITEST(Len,Icomp,Itrue,Kprint)
    IMPLICIT NONE
    !>
    !***
    !  Compare arrays ICOMP and ITRUE.
    !***
    ! **Library:**   SLATEC
    !***
    ! **Type:**      INTEGER (ITEST-I)
    !***
    ! **Keywords:**  QUICK CHECK
    !***
    ! **Author:**  Lawson, C. L., (JPL)
    !***
    ! **Description:**
    !
    !   This subroutine compares the arrays ICOMP and ITRUE of length LEN
    !   for equality.  In the case of an unequal compare, appropriate
    !   messages are written.
    !
    !***
    ! **Routines called:**  (NONE)
    !***
    ! COMMON BLOCKS    COMBLA

    !* REVISION HISTORY  (YYMMDD)
    !   741210  DATE WRITTEN
    !   890831  Modified array declarations.  (WRB)
    !   890831  REVISION DATE from Version 3.2
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    !   920211  Code restructured and information added to the DESCRIPTION
    !           section.  (WRB)
    
    INTEGER i, ICAse, id, INCx, INCy, Kprint, Len, MODe, N, NPRint
    INTEGER Icomp(*), Itrue(*)
    LOGICAL PASs
    COMMON /COMBLA/ NPRint, ICAse, N, INCx, INCy, MODe, PASs
    !* FIRST EXECUTABLE STATEMENT  ITEST
    DO i = 1, Len
      IF ( Icomp(i)/=Itrue(i) ) THEN
        !
        !         Here ICOMP(I) is not equal to ITRUE(I).
        !
        IF ( PASs ) THEN
          !
          !           Print FAIL message and header.
          !
          PASs = .FALSE.
          IF ( Kprint>=3 ) THEN
            WRITE (NPRint,99001)
            99001 FORMAT ('+',39X,'FAIL')
            WRITE (NPRint,99002)
            99002 FORMAT ('0CASE  N INCX INCY MODE  I',29X,'COMP(I)',29X,'TRUE(I)',&
              2X,'DIFFERENCE'/1X)
          ENDIF
        ENDIF
        IF ( Kprint>=3 ) THEN
          id = Icomp(i) - Itrue(i)
          WRITE (NPRint,99003) ICAse, N, INCx, INCy, MODe, i, Icomp(i), &
            Itrue(i), id
          99003 FORMAT (1X,I4,I3,3I5,I3,2I36,I12)
        ENDIF
      ENDIF
    ENDDO
    RETURN
  END SUBROUTINE ITEST
  !** STEST
  SUBROUTINE STEST(Len,Scomp,Strue,Ssize,Sfac,Kprint)
    IMPLICIT NONE
    !>
    !***
    !  Compare arrays SCOMP and STRUE.
    !***
    ! **Library:**   SLATEC
    !***
    ! **Type:**      SINGLE PRECISION (STEST-S, DTEST-D)
    !***
    ! **Keywords:**  QUICK CHECK
    !***
    ! **Author:**  Lawson, C. L., (JPL)
    !***
    ! **Description:**
    !
    !   This subroutine compares arrays SCOMP and STRUE of length LEN to
    !   see if the term by term differences, multiplied by SFAC, are
    !   negligible.  In the case of a significant difference, appropriate
    !   messages are written.
    !
    !***
    ! **Routines called:**  R1MACH
    !***
    ! COMMON BLOCKS    COMBLA

    !* REVISION HISTORY  (YYMMDD)
    !   741210  DATE WRITTEN
    !   890831  Modified array declarations.  (WRB)
    !   890831  REVISION DATE from Version 3.2
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    !   900820  Modified IF test to use function DIFF and made cosmetic
    !           changes to routine.  (WRB)
    !   901005  Removed usage of DIFF in favour of R1MACH.  (RWC)
    !   910501  Added TYPE record.  (WRB)
    !   920211  Code restructured and information added to the DESCRIPTION
    !           section.  (WRB)
    
    INTEGER i, ICAse, INCx, INCy, Kprint, Len, MODe, N, NPRint
    REAL Scomp(*), Strue(*), Ssize(*), Sfac, sd, releps, R1MACH
    LOGICAL PASs
    COMMON /COMBLA/ NPRint, ICAse, N, INCx, INCy, MODe, PASs
    SAVE releps
    DATA releps/0.0E0/
    !* FIRST EXECUTABLE STATEMENT  STEST
    IF ( releps==0.0E0 ) releps = R1MACH(4)
    DO i = 1, Len
      sd = ABS(Scomp(i)-Strue(i))
      IF ( Sfac*sd>ABS(Ssize(i))*releps ) THEN
        !
        !         Here SCOMP(I) is not close to STRUE(I).
        !
        IF ( PASs ) THEN
          !
          !           Print FAIL message and header.
          !
          PASs = .FALSE.
          IF ( Kprint>=3 ) THEN
            WRITE (NPRint,99001)
            99001 FORMAT ('+',39X,'FAIL')
            WRITE (NPRint,99002)
            99002 FORMAT ('0CASE  N INCX INCY MODE  I',29X,'COMP(I)',29X,'TRUE(I)',&
              2X,'DIFFERENCE',5X,'SIZE(I)'/1X)
          ENDIF
        ENDIF
        IF ( Kprint>=3 ) WRITE (NPRint,99003) ICAse, N, INCx, INCy, MODe, &
          i, Scomp(i), Strue(i), sd, Ssize(i)
        99003 FORMAT (1X,I4,I3,3I5,I3,2E36.8,2E12.4)
      ENDIF
    ENDDO
    RETURN
  END SUBROUTINE STEST
  !** DTEST
  SUBROUTINE DTEST(Len,Dcomp,Dtrue,Dsize,Dfac,Kprint)
    IMPLICIT NONE
    !>
    !***
    !  Compare arrays DCOMP and DTRUE.
    !***
    ! **Library:**   SLATEC
    !***
    ! **Type:**      DOUBLE PRECISION (STEST-S, DTEST-D)
    !***
    ! **Keywords:**  QUICK CHECK
    !***
    ! **Author:**  Lawson, C. L., (JPL)
    !***
    ! **Description:**
    !
    !   This subroutine compares arrays DCOMP and DTRUE of length LEN to
    !   see if the term by term differences, multiplied by DFAC, are
    !   negligible.  In the case of a significant difference, appropriate
    !   messages are written.
    !
    !***
    ! **Routines called:**  D1MACH
    !***
    ! COMMON BLOCKS    COMBLA

    !* REVISION HISTORY  (YYMMDD)
    !   741210  DATE WRITTEN
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    !   900820  Modified IF test to use function DDIFF and made cosmetic
    !           changes to routine.  (WRB)
    !   901005  Removed usage of DDIFF in favour of D1MACH.  (RWC)
    !   910501  Added TYPE record.  (WRB)
    !   920211  Code restructured and information added to the DESCRIPTION
    !           section.  (WRB)
    
    INTEGER i, ICAse, INCx, INCy, Kprint, Len, MODe, N, NPRint
    REAL(8) :: Dcomp(*), Dtrue(*), Dsize(*), Dfac, dd, releps, D1MACH
    LOGICAL PASs
    COMMON /COMBLA/ NPRint, ICAse, N, INCx, INCy, MODe, PASs
    SAVE releps
    DATA releps/0.0D0/
    !* FIRST EXECUTABLE STATEMENT  DTEST
    IF ( releps==0.0D0 ) releps = D1MACH(4)
    DO i = 1, Len
      dd = ABS(Dcomp(i)-Dtrue(i))
      IF ( Dfac*dd>ABS(Dsize(i))*releps ) THEN
        !
        !         Here DCOMP(I) is not close to DTRUE(I).
        !
        IF ( PASs ) THEN
          !
          !           Print FAIL message and header.
          !
          PASs = .FALSE.
          IF ( Kprint>=3 ) THEN
            WRITE (NPRint,99001)
            99001 FORMAT ('+',39X,'FAIL')
            WRITE (NPRint,99002)
            99002 FORMAT ('0CASE  N INCX INCY MODE  I',29X,'COMP(I)',29X,'TRUE(I)',&
              2X,'DIFFERENCE',5X,'SIZE(I)'/1X)
          ENDIF
        ENDIF
        IF ( Kprint>=3 ) WRITE (NPRint,99003) ICAse, N, INCx, INCy, MODe, &
          i, Dcomp(i), Dtrue(i), dd, Dsize(i)
        99003 FORMAT (1X,I4,I3,3I5,I3,2D36.18,2D12.4)
      ENDIF
    ENDDO
    RETURN
  END SUBROUTINE DTEST

  SUBROUTINE CTEST(Len,Ccomp,Ctrue,Csize,Cfac,Kprint)
    IMPLICIT NONE
    !>
    !***
    !  Compare arrays DCOMP and DTRUE.
    !***
    ! **Library:**   SLATEC
    !***
    ! **Type:**      DOUBLE PRECISION (STEST-S, DTEST-D)
    !***
    ! **Keywords:**  QUICK CHECK
    !***
    ! **Author:**  Lawson, C. L., (JPL)
    !***
    ! **Description:**
    !
    !   This subroutine compares arrays DCOMP and DTRUE of length LEN to
    !   see if the term by term differences, multiplied by DFAC, are
    !   negligible.  In the case of a significant difference, appropriate
    !   messages are written.
    !
    !***
    ! **Routines called:**  D1MACH
    !***
    ! COMMON BLOCKS    COMBLA

    !* REVISION HISTORY  (YYMMDD)
    !   741210  DATE WRITTEN
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    !   900820  Modified IF test to use function DDIFF and made cosmetic
    !           changes to routine.  (WRB)
    !   901005  Removed usage of DDIFF in favour of D1MACH.  (RWC)
    !   910501  Added TYPE record.  (WRB)
    !   920211  Code restructured and information added to the DESCRIPTION
    !           section.  (WRB)
    
    INTEGER i, ICAse, INCx, INCy, Kprint, Len, MODe, N, NPRint
    COMPLEX :: Ccomp(*), Ctrue(*), Csize(*)
    REAL :: Cfac, dd, releps, R1MACH, CABS1
    LOGICAL PASs
    COMMON /COMBLA/ NPRint, ICAse, N, INCx, INCy, MODe, PASs
    SAVE releps
    DATA releps/0.0/
    !* FIRST EXECUTABLE STATEMENT  DTEST
    IF ( releps==0.0 ) releps = R1MACH(4)
    DO i = 1, Len
      dd = CABS1(Ccomp(i)-Ctrue(i))
      IF ( Cfac*dd>ABS(Csize(i))*releps ) THEN
        !
        !         Here DCOMP(I) is not close to DTRUE(I).
        !
        IF ( PASs ) THEN
          !
          !           Print FAIL message and header.
          !
          PASs = .FALSE.
          IF ( Kprint>=3 ) THEN
            WRITE (NPRint,99001)
            99001 FORMAT ('+',39X,'FAIL')
            WRITE (NPRint,99002)
            99002 FORMAT ('0CASE  N INCX INCY MODE  I',29X,'COMP(I)',29X,'TRUE(I)',&
              2X,'DIFFERENCE',5X,'SIZE(I)'/1X)
          ENDIF
        ENDIF
        IF ( Kprint>=3 ) WRITE (NPRint,99003) ICAse, N, INCx, INCy, MODe, &
          i, Ccomp(i), Ctrue(i), dd, Csize(i)
        99003 FORMAT (1X,I4,I3,3I5,I3,2D36.18,2D12.4)
      ENDIF
    ENDDO
    RETURN
  END SUBROUTINE CTEST
END MODULE TEST17_MOD
!** TEST17
PROGRAM TEST17
  USE TEST17_MOD
  IMPLICIT NONE
  !>
  !***
  !  Driver for testing SLATEC subprograms
  !***
  ! **Library:**   SLATEC
  !***
  ! **Category:**  D1
  !***
  ! **Type:**      ALL (TEST17-A)
  !***
  ! **Keywords:**  QUICK CHECK DRIVER
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
  !        BLAS SUBPROGRAMS
  !
  !***
  ! **References:**  Kirby W. Fong, Thomas H. Jefferson, Tokihiko Suyehiro
  !                 and Lee Walton, Guide to the SLATEC Common Mathema-
  !                 tical Library, April 10, 1990.
  !***
  ! **Routines called:**  BLACHK, I1MACH, XERMAX, XSETF, XSETUN

  !* REVISION HISTORY  (YYMMDD)
  !   890618  DATE WRITTEN
  !   890618  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900524  Cosmetic changes to code.  (WRB)
  
  INTEGER I1MACH
  INTEGER ipass, kprint, lin, lun, nfail
  !* FIRST EXECUTABLE STATEMENT  TEST17
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
  !     Test BLAS
  !
  CALL BLACHK(lun,kprint,ipass)
  IF ( ipass==0 ) nfail = nfail + 1
  !
  !     Write PASS or FAIL message
  !
  IF ( nfail==0 ) THEN
    WRITE (lun,99001)
    99001 FORMAT (/' --------------TEST17 PASSED ALL TESTS----------------')
  ELSE
    WRITE (lun,99002) nfail
    99002 FORMAT (/' ************* WARNING -- ',I5,&
      ' TEST(S) FAILED IN PROGRAM TEST17 *************')
  ENDIF
  STOP
END PROGRAM TEST17
