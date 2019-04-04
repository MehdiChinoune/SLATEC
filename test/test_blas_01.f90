MODULE TEST17_MOD
  IMPLICIT NONE
  INTEGER :: NPRint, ICAse, N, INCx, INCy, MODe
  LOGICAL :: PASs

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

    INTEGER Lun, Ipass, Kprint
    REAL, PARAMETER :: sfac = .625E-1, sdfac = .50
    REAL(8), PARAMETER :: dfac = .625D-1, dqfac = 0.625D-1
    INTEGER, PARAMETER :: jtest(38) = [ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, &
      1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 ]
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

    INTEGER :: Kprint
    CHARACTER(6), PARAMETER :: l(38)  = [ '  SDOT', ' DSDOT', 'SDSDOT', &
      '  DDOT', 'DQDOTI', 'DQDOTA', ' CDOTC', ' CDOTU', ' SAXPY', ' DAXPY', &
      ' CAXPY', ' SROTG', ' DROTG', '  SROT', '  DROT', 'SROTMG', 'DROTMG', &
      ' SROTM', ' DROTM', ' SCOPY', ' DCOPY', ' CCOPY', ' SSWAP', ' DSWAP', &
      ' CSWAP', ' SNRM2', ' DNRM2', 'SCNRM2', ' SASUM', ' DASUM', 'SCASUM', &
      ' SSCAL', ' DSCAL', ' CSCAL', 'CSSCAL', 'ISAMAX', 'IDAMAX', 'ICAMAX' ]
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

    INTEGER i, jump, k, Kprint
    REAL sa(1), sb(1), sc(1), Sfac, ss(1)
    REAL strue(9), stemp(9), stmp(1)
    REAL(8) :: dc(1), ds(1), da(10), Dfac, db(1), dtemp(9)
    REAL, PARAMETER :: zero = 0.
    REAL(8), PARAMETER :: dzero = 0.D0
    REAL(8), PARAMETER :: da1(8) = [ .3D0, .4D0, -.3D0, -.4D0, -.3D0, 0.D0, 0.D0, 1.D0 ]
    REAL(8), PARAMETER :: db1(8) = [ .4D0, .3D0, .4D0, .3D0, -.4D0, 0.D0, 1.D0, 0.D0 ]
    REAL(8), PARAMETER :: dc1(8) = [ .6D0, .8D0, -.6D0, .8D0, .6D0, 1.D0, 0.D0, 1.D0 ]
    REAL(8), PARAMETER :: ds1(8) = [ .8D0, .6D0, .8D0, -.6D0, .8D0, 0.D0, 1.D0, 0.D0 ]
    REAL(8), PARAMETER :: datrue(8) = [ .5D0, .5D0, .5D0, -.5D0, -.5D0, 0.D0, 1.D0, 1.D0 ]
    REAL(8), PARAMETER :: dbtrue(8) = [ 1.D0/.6D0, .6D0, -1.D0/.6D0, -.6D0, &
      1.D0/.6D0, 0.D0, 1.D0, 0.D0 ]
    ! INPUT FOR MODIFIED GIVENS
    REAL(8), PARAMETER :: dab(4,9) = RESHAPE( [ &
      .1D0, .3D0, 1.2D0, .2D0, .7D0, .2D0, .6D0, 4.2D0, &
      0.D0, 0.D0, 0.D0, 0.D0, 4.D0, -1.D0, 2.D0, 4.D0, &
      6.D-10, 2.D-2, 1.D5, 10.D0, 4.D10, 2.D-2, 1.D-5, 10.D0, &
      2.D-10, 4.D-2, 1.D5, 10.D0, 2.D10, 4.D-2, 1.D-5, 10.D0, &
      4.D0, -2.D0, 8.D0, 4.D0 ], [4,9] )
    ! TRUE RESULTS FOR MODIFIED GIVENS
    REAL(8) :: dtrue(9,9) = RESHAPE( [ &
      0.D0, 0.D0, 1.3D0, .2D0, 0.D0, 0.D0, 0.D0, .5D0, 0.D0, &
      0.D0, 0.D0, 4.5D0, 4.2D0, 1.D0, .5D0, 0.D0, 0.D0, 0.D0, &
      0.D0, 0.D0, 0.D0, 0.D0, -2.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      0.D0, 0.D0, 0.D0, 4.D0, -1.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      0.D0, 15.D-3, 0.D0, 10.D0, -1.D0, 0.D0, -1.D-4, 0.D0, 1.D0, &
      0.D0, 0.D0, 6144.D-5, 10.D0, -1.D0, 4096.D0, -1.D6, 0.D0, &
      1.D0, 0.D0, 0.D0, 15.D0, 10.D0, -1.D0, 5.D-5, 0.D0, 1.D0, &
      0.D0, 0.D0, 0.D0, 15.D0, 10.D0, -1.D0, 5.D5, -4096.D0, 1.D0, &
      4096.D-6, 0.D0, 0.D0, 7.D0, 4.D0, 0.D0, 0.D0, -.5D0, -.25D0, 0.D0 ], [9,9] )
    !                   4096 = 2 ** 12
    REAL(8), PARAMETER :: d12 = 4096.D0
    !* FIRST EXECUTABLE STATEMENT  CHECK0
    !
    ! COMPUTE TRUE VALUES WHICH CANNOT BE PRESTORED
    ! IN DECIMAL NOTATION.
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
          CALL DTEST(1,da,datrue(k),datrue(k),Dfac,Kprint)
          CALL DTEST(1,db,dbtrue(k),dbtrue(k),Dfac,Kprint)
          CALL DTEST(1,dc,dc1(k),dc1(k),Dfac,Kprint)
          CALL DTEST(1,ds,ds1(k),ds1(k),Dfac,Kprint)
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
          stmp = REAL(datrue(k))
          CALL STEST(1,sa,stmp,stmp,Sfac,Kprint)
          stmp = REAL(dbtrue(k))
          CALL STEST(1,sb,stmp,stmp,Sfac,Kprint)
          stmp = REAL(dc1(k))
          CALL STEST(1,sc,stmp,stmp,Sfac,Kprint)
          stmp = REAL(ds1(k))
          CALL STEST(1,ss,stmp,stmp,Sfac,Kprint)
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

    INTEGER i, ICAMAX, IDAMAX, ISAMAX, jump, len, np1, Kprint, itmp(1)
    REAL SASUM, SCASUM, SCNRM2, Sfac, SNRM2
    REAL(8) :: dx(8), Dfac, dtmp(1)
    REAL(8) :: DNRM2, DASUM
    REAL strue(8), sx(8), stmp(1), stmp2(1)
    COMPLEX cx(8)
    !
    REAL, PARAMETER :: sa = .3
    REAL(8), PARAMETER :: da = .3D0
    COMPLEX, PARAMETER :: ca = (.4,-.7)
    REAL(8), PARAMETER :: dv(8,5,2) = RESHAPE( [ &
      .1D0, 2.D0, 2.D0, 2.D0, 2.D0, 2.D0, 2.D0, 2.D0, .3D0, &
      3.D0, 3.D0, 3.D0, 3.D0, 3.D0, 3.D0, 3.D0, .3D0, -.4D0, &
      4.D0, 4.D0, 4.D0, 4.D0, 4.D0, 4.D0, .2D0, -.6D0, .3D0, &
      5.D0, 5.D0, 5.D0, 5.D0, 5.D0, .1D0, -.3D0, .5D0, -.1D0, &
      6.D0, 6.D0, 6.D0, 6.D0, .1D0, 8.D0, 8.D0, 8.D0, 8.D0, 8.D0, &
      8.D0, 8.D0, .3D0, 9.D0, 9.D0, 9.D0, 9.D0, 9.D0, 9.D0, 9.D0, &
      .3D0, 2.D0, -.4D0, 2.D0, 2.D0, 2.D0, 2.D0, 2.D0, .2D0, &
      3.D0, -.6D0, 5.D0, .3D0, 2.D0, 2.D0, 2.D0, .1D0, 4.D0, &
      -.3D0, 6.D0, -.5D0, 7.D0, -.1D0, 3.D0 ], [8,5,2] )
    !     COMPLEX TEST VECTORS
    COMPLEX, PARAMETER :: cv(8,5,2) = RESHAPE( [ &
      (.1,.1), (1.,2.), (1.,2.), (1.,2.), (1.,2.), (1.,2.), &
      (1.,2.), (1.,2.), (.3,-.4), (3.,4.), (3.,4.), (3.,4.), (3.,4.), &
      (3.,4.), (3.,4.), (3.,4.), (.1,-.3), (.5,-.1), (5.,6.), &
      (5.,6.), (5.,6.), (5.,6.), (5.,6.), (5.,6.), (.1,.1), (-.6,.1), &
      (.1,-.3), (7.,8.), (7.,8.), (7.,8.), (7.,8.), (7.,8.), &
      (.3,.1), (.1,.4), (.4,.1), (.1,.2), (2.,3.), (2.,3.), (2.,3.), &
      (2.,3.), (.1,.1), (4.,5.), (4.,5.), (4.,5.), (4.,5.), (4.,5.), &
      (4.,5.), (4.,5.), (.3,-.4), (6.,7.), (6.,7.), (6.,7.), (6.,7.), &
      (6.,7.), (6.,7.), (6.,7.), (.1,-.3), (8.,9.), (.5,-.1), &
      (2.,5.), (2.,5.), (2.,5.), (2.,5.), (2.,5.), (.1,.1), (3.,6.), &
      (-.6,.1), (4.,7.), (.1,-.3), (7.,2.), (7.,2.), (7.,2.), (.3,.1), &
      (5.,8.), (.1,.4), (6.,9.), (.4,.1), (8.,3.), (.1,.2), (9.,4.) ], [8,5,2] )
    !
    REAL, PARAMETER :: strue2(5) = [ .0, .5, .6, .7, .7 ]
    REAL, PARAMETER :: strue4(5) = [ .0, .7, 1., 1.3, 1.7 ]
    REAL(8), PARAMETER :: dtrue1(5) = [ .0D0, .3D0, .5D0, .7D0, .6D0 ]
    REAL(8), PARAMETER :: dtrue3(5) = [ .0D0, .3D0, .7D0, 1.1D0, 1.D0 ]
    REAL(8), PARAMETER :: dtrue5(8,5,2) = RESHAPE( [ &
      .10D0, 2.D0, 2.D0, 2.D0, 2.D0, 2.D0, 2.D0, 2.D0, &
      .09D0, 3.D0, 3.D0, 3.D0, 3.D0, 3.D0, 3.D0, 3.D0, .09D0, &
      -.12D0, 4.D0, 4.D0, 4.D0, 4.D0, 4.D0, 4.D0, .06D0, -.18D0, &
      .09D0, 5.D0, 5.D0, 5.D0, 5.D0, 5.D0, .03D0, -.09D0, .15D0, &
      -.03D0, 6.D0, 6.D0, 6.D0, 6.D0, .10D0, 8.D0, 8.D0, 8.D0, &
      8.D0, 8.D0, 8.D0, 8.D0, .09D0, 9.D0, 9.D0, 9.D0, 9.D0, &
      9.D0, 9.D0, 9.D0, .09D0, 2.D0, -.12D0, 2.D0, 2.D0, 2.D0, &
      2.D0, 2.D0, .06D0, 3.D0, -.18D0, 5.D0, .09D0, 2.D0, 2.D0, &
      2.D0, .03D0, 4.D0, -.09D0, 6.D0, -.15D0, 7.D0, -.03D0, 3.D0 ], [8,5,2] )
    !
    COMPLEX, PARAMETER :: ctrue5(8,5,2) = RESHAPE( [ &
      (.1,.1), (1.,2.), (1.,2.), (1.,2.), (1.,2.), (1.,2.), &
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
      (.32,.09), (6.,9.), (.23,-.24), (8.,3.), (.18,.01), (9.,4.) ], [8,5,2] )
    !
    COMPLEX, PARAMETER :: ctrue6(8,5,2) = RESHAPE( [ &
      (.1,.1), (1.,2.), (1.,2.), (1.,2.), (1.,2.), (1.,2.), &
      (1.,2.), (1.,2.), (.09,-.12), (3.,4.), (3.,4.), (3.,4.), &
      (3.,4.), (3.,4.), (3.,4.), (3.,4.), (.03,-.09), (.15,-.03), &
      (5.,6.), (5.,6.), (5.,6.), (5.,6.), (5.,6.), (5.,6.), (.03,.03), &
      (-.18,.03), (.03,-.09), (7.,8.), (7.,8.), (7.,8.), (7.,8.), &
      (7.,8.), (.09,.03), (.03,.12), (.12,.03), (.03,.06), (2.,3.), &
      (2.,3.), (2.,3.), (2.,3.), (.1,.1), (4.,5.), (4.,5.), (4.,5.), &
      (4.,5.), (4.,5.), (4.,5.), (4.,5.), (.09,-.12), (6.,7.), &
      (6.,7.), (6.,7.), (6.,7.), (6.,7.), (6.,7.), (6.,7.), &
      (.03,-.09), (8.,9.), (.15,-.03), (2.,5.), (2.,5.), (2.,5.), &
      (2.,5.), (2.,5.), (.03,.03), (3.,6.), (-.18,.03), (4.,7.), &
      (.03,-.09), (7.,2.), (7.,2.), (7.,2.), (.09,.03), (5.,8.), &
      (.03,.12), (6.,9.), (.12,.03), (8.,3.), (.03,.06), (9.,4.) ], [8,5,2] )
    !
    !
    INTEGER, PARAMETER :: itrue2(5) = [ 0, 1, 2, 2, 3 ]
    INTEGER, PARAMETER :: itrue3(5) = [ 0, 1, 2, 2, 2 ]
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
            dtmp = DNRM2(N,dx,INCx)
            CALL DTEST(1,dtmp,dtrue1(np1),dtrue1(np1),Dfac,Kprint)
          CASE (3)
            ! 28. SCNRM2
            stmp = SCNRM2(N,cx,INCx)
            CALL STEST(1,stmp,strue2(np1),strue2(np1),Sfac,Kprint)
          CASE (4)
            ! 29. SASUM
            stmp = SASUM(N,sx,INCx)
            stmp2 = REAL( dtrue3(np1), 4 )
            CALL STEST(1,stmp,stmp2,stmp2,Sfac,Kprint)
          CASE (5)
            ! 30. DASUM
            dtmp = DASUM(N,dx,INCx)
            CALL DTEST(1,dtmp,dtrue3(np1),dtrue3(np1),Dfac,Kprint)
          CASE (6)
            ! 31. SCASUM
            stmp = SCASUM(N,cx,INCx)
            CALL STEST(1,stmp,strue4(np1),strue4(np1),Sfac,Kprint)
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
            itmp = ISAMAX(N,sx,INCx)
            CALL ITEST(1,itmp,itrue2(np1),Kprint)
          CASE (12)
            ! 37. IDAMAX
            itmp = IDAMAX(N,dx,INCx)
            CALL ITEST(1,itmp,itrue2(np1),Kprint)
          CASE (13)
            ! 38. ICAMAX
            itmp = ICAMAX(N,cx,INCx)
            CALL ITEST(1,itmp,itrue3(np1),Kprint)
          CASE DEFAULT
            ! 26. SNRM2
            stmp = SNRM2(N,sx,INCx)
            stmp2 = REAL( dtrue1(np1), 4 )
            CALL STEST(1,stmp,stmp2,stmp2,Sfac,Kprint)
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

    INTEGER i, j, ki, kn, kni, kpar, ksize, lenx, leny, mx, my, Kprint
    REAL Sdfac, SDOT, SDSDOT, Sfac
    REAL sx(7), sy(7), stx(7), sty(7), ssize(7), qc(30), sparam(5), stmp(1), stmp2(1)
    REAL(8) :: dx(7), dy(7), dparam(5), dsize(7), dtx(7), dty(7), dtmp(1)
    REAL(8) :: DSDOT, DDOT, DQDOTI, DQDOTA, Dfac, Dqfac
    !
    COMPLEX cx(7), cy(7), ctmp(1), ctmp2(4)
    COMPLEX CDOTC, CDOTU
    REAL, PARAMETER :: sa = .3, sb = .1
    REAL(8), PARAMETER :: da = .3D0, db = .25D0
    COMPLEX, PARAMETER :: ca = (.4,-.7)
    INTEGER, PARAMETER :: incxs(4) = [ 1, 2, -2, -1 ]
    INTEGER, PARAMETER :: incys(4) = [ 1, -2, 1, -2 ]
    INTEGER, PARAMETER :: lens(4,2) = RESHAPE( [1, 1, 2, 4, 1, 1, 3, 7], [4,2] )
    INTEGER, PARAMETER :: ns(4) = [ 0, 1, 2, 4 ]
    REAL, PARAMETER :: sc = .8, ss = .6
    REAL(8), PARAMETER :: dc = .8D0, ds = .6D0
    REAL(8), PARAMETER :: dx1(7) = [ .6D0, .1D0, -.5D0, .8D0, .9D0, -.3D0, -.4D0 ]
    REAL(8), PARAMETER :: dy1(7) = [ .5D0, -.9D0, .3D0, .7D0, -.6D0, .2D0, .8D0 ]
    REAL(8), PARAMETER :: dx2(7) = [ 1.D0, .01D0, .02D0, 1.D0, .06D0, 2.D0, 1.D0 ]
    REAL(8), PARAMETER :: dy2(7) = [ 1.D0, .04D0, -.03D0, -1.D0, .05D0, 3.D0, -1.D0 ]
    !            THE TERMS D11(3,2) AND D11(4,2) WILL BE SET BY
    !            COMPUTATION AT RUN TIME.
    COMPLEX, PARAMETER :: cx1(7) = [ (.7,-.8), (-.4,-.7), (-.1,-.9), (.2,-.8), &
      (-.9,-.4), (.1,.4), (-.6,.6) ]
    COMPLEX, PARAMETER :: cy1(7) = [ (.6,-.6), (-.9,.5), (.7,-.6), (.1,-.5), &
      (-.1,-.2), (-.5,-.3), (.8,-.7) ]
    !
    !                             FOR DQDOTI AND DQDOTA
    !
    REAL(8), PARAMETER :: dt2(4,4,2) = RESHAPE( [ 0.25D0, 1.25D0, 1.2504D0, 0.2498D0, &
      0.25D0, 1.25D0, 0.24D0, 0.2492D0, 0.25D0, 1.25D0, 0.31D0, 0.2518D0, &
      0.25D0, 1.25D0, 1.2497D0, 0.2507D0, 0.D0, 2.D0, 2.0008D0, -.0004D0, &
      0.D0, 2.D0, -.02D0, -.0016D0, 0.D0, 2.D0, .12D0, .0036D0, &
      0.D0, 2.D0, 1.9994D0, .0014D0 ], [4,4,2] )
    REAL(8), PARAMETER :: dt7(4,4) = RESHAPE( [ 0.D0, .30D0, .21D0, .62D0, &
      0.D0, .30D0, -.07D0, .85D0, 0.D0, .30D0, -.79D0, -.74D0, &
      0.D0, .30D0, .33D0, 1.27D0 ], [4,4] )
    REAL, PARAMETER :: st7b(4,4) = RESHAPE( [ .1, .4, .31, .72, .1, .4, .03, .95, &
      .1, .4, -.69, -.64, .1, .4, .43, 1.37 ], [4,4] )
    !
    !                       FOR CDOTU
    !
    COMPLEX, PARAMETER :: ct7(4,4) = RESHAPE( [ &
      (0.,0.), (-.06,-.90), (.65,-.47), (-.34,-1.22), &
      (0.,0.), (-.06,-.90), (-.59,-1.46), (-1.04,-.04), &
      (0.,0.), (-.06,-.90), (-.83,.59), (.07,-.37), &
      (0.,0.), (-.06,-.90), (-.76,-1.15), (-1.33,-1.82) ], [4,4] )
    !
    !                       FOR CDOTC
    !
    COMPLEX, PARAMETER :: ct6(4,4) = RESHAPE( [ &
      (0.,0.), (.90,0.06), (.91,-.77), (1.80,-.10), &
      (0.,0.), (.90,0.06), (1.45,.74), (.20,.90), &
      (0.,0.), (.90,0.06), (-.55,.23), (.83,-.39), &
      (0.,0.), (.90,0.06), (1.04,0.79), (1.95,1.22) ], [4,4] )
    !
    REAL(8), PARAMETER :: dt8(7,4,4) = RESHAPE( [ &
      .5D0,  0.D0,   0.D0,   0.D0,   0.D0,  0.D0,  0.D0, &
      .68D0, 0.D0,   0.D0,   0.D0,   0.D0,  0.D0,  0.D0, &
      .68D0, -.87D0, 0.D0,   0.D0,   0.D0,  0.D0,  0.D0, &
      .68D0, -.87D0,  .15D0,  .94D0, 0.D0,  0.D0,  0.D0, &
      .5D0,  0.D0,   0.D0,   0.D0,   0.D0,  0.D0,  0.D0, &
      .68D0, 0.D0,   0.D0,   0.D0,   0.D0,  0.D0,  0.D0, &
      .35D0, -.9D0,   .48D0, 0.D0,   0.D0,  0.D0,  0.D0, &
      .38D0, -.9D0,   .57D0,  .7D0,  -.75D0, .2D0,  .98D0, &
      .5D0,  0.D0,   0.D0,   0.D0,   0.D0,  0.D0,  0.D0, &
      .68D0, 0.D0,   0.D0,   0.D0,   0.D0,  0.D0,  0.D0, &
      .35D0, -.72D0, 0.D0,   0.D0,   0.D0,  0.D0,  0.D0, &
      .38D0, -.63D0,  .15D0,  .88D0, 0.D0,  0.D0,  0.D0, &
      .5D0,  0.D0,   0.D0,   0.D0,   0.D0,  0.D0,  0.D0, &
      .68D0, 0.D0,   0.D0,   0.D0,   0.D0,  0.D0,  0.D0, &
      .68D0, -.9D0,   .33D0, 0.D0,   0.D0,  0.D0,  0.D0, &
      .68D0, -.9D0,   .33D0,  .7D0,  -.75D0, .2D0, 1.04D0 ], [7,4,4] )
    !
    COMPLEX, PARAMETER :: ct8(7,4,4) = RESHAPE( [ &
      (.6,-.6), (0.,0.), (0.,0.), (0.,0.), (0.,0.), (0.,0.), (0.,0.), &
      (.32,-1.41), (0.,0.), (0.,0.), (0.,0.), (0.,0.), (0.,0.), (0.,0.), &
      (.32,-1.41), (-1.55,.5), (0.,0.), (0.,0.), (0.,0.), (0.,0.), (0.,0.), &
      (.32,-1.41), (-1.55,.5), (.03,-.89), (-.38,-.96), (0.,0.), (0.,0.), (0.,0.), &
      (.6,-.6), (0.,0.), (0.,0.), (0.,0.), (0.,0.), (0.,0.), (0.,0.), &
      (.32,-1.41), (0.,0.), (0.,0.), (0.,0.), (0.,0.), (0.,0.), (0.,0.), &
      (-.07,-.89), (-.9,.5), (.42,-1.41), (0.,0.), (0.,0.), (0.,0.), (0.,0.), &
      (.78,.06), (-.9,.5), (.06,-.13), (.1,-.5), (-.77,-.49), (-.5,-.3), (.52,-1.51), &
      (.6,-.6), (0.,0.), (0.,0.), (0.,0.), (0.,0.), (0.,0.), (0.,0.), &
      (.32,-1.41), (0.,0.), (0.,0.), (0.,0.), (0.,0.), (0.,0.), (0.,0.), &
      (-.07,-.89), (-1.18,-.31), (0.,0.), (0.,0.), (0.,0.), (0.,0.), (0.,0.), &
      (.78,.06), (-1.54,.97), (.03,-.89), (-.18,-1.31), (0.,0.), (0.,0.), (0.,0.), &
      (.6,-.6), (0.,0.), (0.,0.), (0.,0.), (0.,0.), (0.,0.), (0.,0.), &
      (.32,-1.41), (0.,0.), (0.,0.), (0.,0.), (0.,0.), (0.,0.), (0.,0.), &
      (.32,-1.41), (-.9,.5), (.05,-.6), (0.,0.), (0.,0.), (0.,0.), (0.,0.), &
      (.32,-1.41), (-.9,.5), (.05,-.6), (.1,-.5), (-.77,-.49), (-.5,-.3), (.32,-1.16) &
      ], [7,4,4] )
    !                TRUE X VALUES AFTER ROTATION USING SROT OR DROT.
    REAL(8), PARAMETER :: dt9x(7,4,4) = RESHAPE( [ &
      .6D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      .78D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      .78D0, -.46D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      .78D0, -.46D0, -.22D0, 1.06D0, 0.D0, 0.D0, 0.D0, &
      .6D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      .78D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      .66D0, .1D0, -.1D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      .96D0, .1D0, -.76D0, .8D0, .90D0, -.3D0, -.02D0, &
      .6D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      .78D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      -.06D0, .1D0, -.1D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      .90D0, .1D0, -.22D0, .8D0, .18D0, -.3D0, -.02D0, &
      .6D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      .78D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      .78D0, .26D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      .78D0, .26D0, -.76D0, 1.12D0, 0.D0, 0.D0, 0.D0 ] , [7,4,4] )
    !
    !                TRUE Y VALUES AFTER ROTATION USING SROT OR DROT.
    !
    REAL(8), PARAMETER :: dt9y(7,4,4) = RESHAPE( [ &
      .5D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      .04D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      .04D0, -.78D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      .04D0, -.78D0, .54D0, .08D0, 0.D0, 0.D0, 0.D0, &
      .5D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      .04D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      .7D0, -.9D0, -.12D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      .64D0, -.9D0, -.30D0, .7D0, -.18D0, .2D0, .28D0, &
      .5D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      .04D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      .7D0, -1.08D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      .64D0, -1.26D0, .54D0, .20D0, 0.D0, 0.D0, 0.D0, &
      .5D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      .04D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      .04D0, -.9D0, .18D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      .04D0, -.9D0, .18D0, .7D0, -.18D0, .2D0, .16D0 ], [7,4,4] )
    !
    REAL(8), PARAMETER :: dt10x(7,4,4) = RESHAPE( [ &
      .6D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      .5D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      .5D0, -.9D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      .5D0, -.9D0, .3D0, .7D0, 0.D0, 0.D0, 0.D0, &
      .6D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      .5D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      .3D0, .1D0, .5D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      .8D0, .1D0, -.6D0, .8D0, .3D0, -.3D0, .5D0, &
      .6D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      .5D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      -.9D0, .1D0, .5D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      .7D0, .1D0, .3D0, .8D0, -.9D0, -.3D0, .5D0, &
      .6D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      .5D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      .5D0, .3D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      .5D0, .3D0, -.6D0, .8D0, 0.D0, 0.D0, 0.D0 ], [7,4,4] )
    !
    REAL(8), PARAMETER :: dt10y(7,4,4) = RESHAPE( [ &
      .5D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      .6D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      .6D0, .1D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      .6D0, .1D0, -.5D0, .8D0, 0.D0, 0.D0, 0.D0, &
      .5D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      .6D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      -.5D0, -.9D0, .6D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      -.4D0, -.9D0, .9D0, .7D0, -.5D0, .2D0, .6D0, &
      .5D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      .6D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      -.5D0, .6D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      -.4D0, .9D0, -.5D0, .6D0, 0.D0, 0.D0, 0.D0, &
      .5D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      .6D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      .6D0, -.9D0, .1D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      .6D0, -.9D0, .1D0, .7D0, -.5D0, .2D0, .8D0 ], [7,4,4] )
    !
    COMPLEX, PARAMETER :: ct10x(7,4,4) = RESHAPE( [ &
      (.7,-.8), (0.,0.), (0.,0.), (0.,0.), (0.,0.), (0.,0.), (0.,0.), &
      (.6,-.6), (0.,0.), (0.,0.), (0.,0.), (0.,0.), (0.,0.),  (0.,0.), &
      (.6,-.6), (-.9,.5), (0.,0.), (0.,0.), (0.,0.), (0.,0.), (0.,0.), &
      (.6,-.6), (-.9,.5), (.7,-.6), (.1,-.5), (0.,0.), (0.,0.), (0.,0.), &
      (.7,-.8), (0.,0.), (0.,0.), (0.,0.), (0.,0.), (0.,0.), (0.,0.), &
      (.6,-.6), (0.,0.), (0.,0.), (0.,0.), (0.,0.), (0.,0.), (0.,0.), &
      (.7,-.6), (-.4,-.7), (.6,-.6), (0.,0.), (0.,0.), (0.,0.), (0.,0.), &
      (.8,-.7), (-.4,-.7), (-.1,-.2), (.2,-.8), (.7,-.6), (.1,.4), (.6,-.6), &
      (.7,-.8), (0.,0.), (0.,0.), (0.,0.), (0.,0.), (0.,0.), (0.,0.), &
      (.6,-.6), (0.,0.), (0.,0.), (0.,0.), (0.,0.), (0.,0.), (0.,0.), &
      (-.9,.5), (-.4,-.7), (.6,-.6), (0.,0.), (0.,0.), (0.,0.), (0.,0.), &
      (.1,-.5), (-.4,-.7), (.7,-.6), (.2,-.8), (-.9,.5), (.1,.4), (.6,-.6), &
      (.7,-.8), (0.,0.), (0.,0.), (0.,0.), (0.,0.), (0.,0.), (0.,0.), &
      (.6,-.6), (0.,0.), (0.,0.), (0.,0.), (0.,0.), (0.,0.), (0.,0.), &
      (.6,-.6), (.7,-.6), (0.,0.), (0.,0.), (0.,0.), (0.,0.), (0.,0.), &
      (.6,-.6), (.7,-.6), (-.1,-.2), (.8,-.7), (0.,0.), (0.,0.), (0.,0.) ], [7,4,4] )
    COMPLEX, PARAMETER :: ct10y(7,4,4) = RESHAPE( [ &
      (.6,-.6), (0.,0.), (0.,0.), (0.,0.), (0.,0.), (0.,0.), (0.,0.), &
      (.7,-.8), (0.,0.), (0.,0.), (0.,0.), (0.,0.), (0.,0.), (0.,0.), &
      (.7,-.8), (-.4,-.7), (0.,0.), (0.,0.), (0.,0.), (0.,0.), (0.,0.), &
      (.7,-.8), (-.4,-.7), (-.1,-.9), (.2,-.8), (0.,0.), (0.,0.), (0.,0.), &
      (.6,-.6), (0.,0.), (0.,0.), (0.,0.), (0.,0.), (0.,0.), (0.,0.), &
      (.7,-.8), (0.,0.), (0.,0.), (0.,0.), (0.,0.), (0.,0.), (0.,0.), &
      (-.1,-.9), (-.9,.5), (.7,-.8), (0.,0.), (0.,0.), (0.,0.), (0.,0.), &
      (-.6,.6), (-.9,.5), (-.9,-.4), (.1,-.5), (-.1,-.9), (-.5,-.3), (.7,-.8), &
      (.6,-.6), (0.,0.), (0.,0.), (0.,0.), (0.,0.), (0.,0.), (0.,0.), &
      (.7,-.8), (0.,0.), (0.,0.), (0.,0.), (0.,0.), (0.,0.), (0.,0.), &
      (-.1,-.9), (.7,-.8), (0.,0.), (0.,0.), (0.,0.), (0.,0.), (0.,0.), &
      (-.6,.6), (-.9,-.4), (-.1,-.9), (.7,-.8), (0.,0.), (0.,0.), (0.,0.), &
      (.6,-.6), (0.,0.), (0.,0.), (0.,0.), (0.,0.), (0.,0.), (0.,0.), &
      (.7,-.8), (0.,0.), (0.,0.), (0.,0.), (0.,0.), (0.,0.), (0.,0.), &
      (.7,-.8), (-.9,.5), (-.4,-.7), (0.,0.), (0.,0.), (0.,0.), (0.,0.), &
      (.7,-.8), (-.9,.5), (-.4,-.7), (.1,-.5), (-.1,-.9), (-.5,-.3), (.2,-.8) &
      ], [7,4,4] )
    !                        TRUE X RESULTS F0R ROTATIONS SROTM AND DROTM
    REAL(8), PARAMETER :: dt19xa(7,4,4) = RESHAPE( [ &
      .6D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      .6D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      .6D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      .6D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      .6D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      -.8D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      -.9D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      3.5D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      .6D0, .1D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      -.8D0, 3.8D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      -.9D0, 2.8D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      3.5D0, -.4D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      .6D0, .1D0, -.5D0, .8D0, 0.D0, 0.D0, 0.D0, &
      -.8D0, 3.8D0, -2.2D0, -1.2D0, 0.D0, 0.D0, 0.D0, &
      -.9D0, 2.8D0, -1.4D0, -1.3D0, 0.D0, 0.D0, 0.D0, &
      3.5D0, -.4D0, -2.2D0, 4.7D0, 0.D0, 0.D0, 0.D0 ], [7,4,4] )
    !
    REAL(8), PARAMETER :: dt19xb(7,4,4) = RESHAPE( [ &
      .6D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      .6D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      .6D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      .6D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      .6D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      -.8D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      -.9D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      3.5D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      .6D0, .1D0, -.5D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      0.D0, .1D0, -3.0D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      -.3D0, .1D0, -2.0D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      3.3D0, .1D0, -2.0D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      .6D0, .1D0, -.5D0, .8D0, .9D0, -.3D0, -.4D0, &
      -2.0D0, .1D0, 1.4D0, .8D0, .6D0, -.3D0, -2.8D0, &
      -1.8D0, .1D0, 1.3D0, .8D0, 0.D0, -.3D0, -1.9D0, &
      3.8D0, .1D0, -3.1D0, .8D0, 4.8D0, -.3D0, -1.5D0 ], [7,4,4] )
    !
    REAL(8), PARAMETER :: dt19xc(7,4,4) = RESHAPE( [ &
      .6D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      .6D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      .6D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      .6D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      .6D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      -.8D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      -.9D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      3.5D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      .6D0, .1D0, -.5D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      4.8D0, .1D0, -3.0D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      3.3D0, .1D0, -2.0D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      2.1D0, .1D0, -2.0D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      .6D0, .1D0, -.5D0, .8D0, .9D0, -.3D0, -.4D0, &
      -1.6D0, .1D0, -2.2D0, .8D0, 5.4D0, -.3D0, -2.8D0, &
      -1.5D0, .1D0, -1.4D0, .8D0, 3.6D0, -.3D0, -1.9D0, &
      3.7D0, .1D0, -2.2D0, .8D0, 3.6D0, -.3D0, -1.5D0 ], [7,4,4] )
    !
    REAL(8), PARAMETER :: dt19xd(7,4,4) = RESHAPE( [ &
      .6D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      .6D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      .6D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      .6D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      .6D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      -.8D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      -.9D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      3.5D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      .6D0, .1D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      -.8D0, -1.0D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      -.9D0, -.8D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      3.5D0, .8D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      .6D0, .1D0, -.5D0, .8D0, 0.D0, 0.D0, 0.D0, &
      -.8D0, -1.0D0, 1.4D0, -1.6D0, 0.D0, 0.D0, 0.D0, &
      -.9D0, -.8D0, 1.3D0, -1.6D0, 0.D0, 0.D0, 0.D0, &
      3.5D0, .8D0, -3.1D0, 4.8D0, 0.D0, 0.D0, 0.D0 ], [7,4,4] )
    !                        TRUE Y RESULTS FOR ROTATIONS SROTM AND DROTM
    REAL(8), PARAMETER :: dt19ya(7,4,4) = RESHAPE( [ &
      .5D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      .5D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      .5D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      .5D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      .5D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      .7D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      1.7D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      -2.6D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      .5D0, -.9D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      .7D0, -4.8D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      1.7D0, -.7D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      -2.6D0, 3.5D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      .5D0, -.9D0, .3D0, .7D0, 0.D0, 0.D0, 0.D0, &
      .7D0, -4.8D0, 3.0D0, 1.1D0, 0.D0, 0.D0, 0.D0, &
      1.7D0, -.7D0, -.7D0, 2.3D0, 0.D0, 0.D0, 0.D0, &
      -2.6D0, 3.5D0, -.7D0, -3.6D0, 0.D0, 0.D0, 0.D0 ], [7,4,4] )
    !
    REAL(8), PARAMETER :: dt19yb(7,4,4) = RESHAPE( [ &
      .5D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      .5D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      .5D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      .5D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      .5D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      .7D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      1.7D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      -2.6D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      .5D0, -.9D0, .3D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      4.0D0, -.9D0, -.3D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      -.5D0, -.9D0, 1.5D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      -1.5D0, -.9D0, -1.8D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      .5D0, -.9D0, .3D0, .7D0, -.6D0, .2D0, .8D0, &
      3.7D0, -.9D0, -1.2D0, .7D0, -1.5D0, .2D0, 2.2D0, &
      -.3D0, -.9D0, 2.1D0, .7D0, -1.6D0, .2D0, 2.0D0, &
      -1.6D0, -.9D0, -2.1D0, .7D0, 2.9D0, .2D0, -3.8D0 ], [7,4,4] )
    !
    REAL(8), PARAMETER :: dt19yc(7,4,4) = RESHAPE( [ &
      .5D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      .5D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      .5D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      .5D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      .5D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      .7D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      1.7D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      -2.6D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      .5D0, -.9D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      4.0D0, -6.3D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      -.5D0, .3D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      -1.5D0, 3.0D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      .5D0, -.9D0, .3D0, .7D0, 0.D0, 0.D0, 0.D0, &
      3.7D0, -7.2D0, 3.0D0, 1.7D0, 0.D0, 0.D0, 0.D0, &
      -.3D0, .9D0, -.7D0, 1.9D0, 0.D0, 0.D0, 0.D0, &
      -1.6D0, 2.7D0, -.7D0, -3.4D0, 0.D0, 0.D0, 0.D0 ], [7,4,4] )
    !
    REAL(8), PARAMETER :: dt19yd(7,4,4) = RESHAPE( [ &
      .5D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      .5D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      .5D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      .5D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      .5D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      .7D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      1.7D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      -2.6D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      .5D0, -.9D0, .3D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      .7D0, -.9D0, 1.2D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      1.7D0, -.9D0, .5D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      -2.6D0, -.9D0, -1.3D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      .5D0, -.9D0, .3D0, .7D0, -.6D0, .2D0, .8D0, &
      .7D0, -.9D0, 1.2D0, .7D0, -1.5D0, .2D0, 1.6D0, &
      1.7D0, -.9D0, .5D0, .7D0, -1.6D0, .2D0, 2.4D0, &
      -2.6D0, -.9D0, -1.3D0, .7D0, 2.9D0, .2D0, -4.0D0 ], [7,4,4] )
    !
    REAL(8), PARAMETER :: dt19x(7,4,16) = RESHAPE( [dt19xa, dt19xb, dt19xc, dt19xd], &
      [7,4,16] )
    REAL(8), PARAMETER :: dt19y(7,4,16) = RESHAPE( [dt19ya, dt19yb, dt19yc, dt19yd], &
      [7,4,16] )
    !
    REAL, PARAMETER :: ssize1(4) = [ 0., .3, 1.6, 3.2 ]
    REAL(8), PARAMETER :: dsize1(4) = [ 0.D0, .3D0, 1.6D0, 3.2D0 ]
    REAL, PARAMETER :: ssize3(4) = [ .1, .4, 1.7, 3.3 ]
    !
    !                         FOR CDOTC AND CDOTU
    !
    COMPLEX, PARAMETER :: csize1(4) = [ (0.,0.), (.9,.9), (1.63,1.73), (2.90,2.78) ]
    REAL, PARAMETER :: ssize2(14,2) = RESHAPE( [ &
      0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., &
      1.17, 1.17, 1.17, 1.17, 1.17, 1.17, 1.17, 1.17, 1.17, 1.17, 1.17, 1.17, &
      1.17, 1.17 ], [14,2] )
    REAL(8), PARAMETER :: dsize2(7,2) = RESHAPE( [ &
      0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      1.17D0, 1.17D0, 1.17D0, 1.17D0, 1.17D0, 1.17D0, 1.17D0 ], [7,2] )
    !
    !                         FOR CAXPY
    !
    COMPLEX, PARAMETER :: csize2(7,2) = RESHAPE( [ &
      (0.,0.), (0.,0.), (0.,0.), (0.,0.), (0.,0.), (0.,0.), (0.,0.), &
      (1.54,1.54), (1.54,1.54), (1.54,1.54), (1.54,1.54), (1.54,1.54), &
      (1.54,1.54), (1.54,1.54) ], [7,2] )
    !
    !                         FOR SROTM AND DROTM
    !
    REAL(8), PARAMETER :: dpar(5,4) = RESHAPE( [ -2.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      -1.D0, 2.D0, -3.D0, -4.D0, 5.D0, 0.D0, 0.D0, 2.D0, -3.D0, 0.D0, &
      1.D0, 5.D0, 2.D0, 0.D0, -4.D0 ], [5,4] )
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
            stmp = REAL(DSDOT(N,sx,INCx,sy,INCy))
            stmp2 = REAL(dt7(kn,ki))
            CALL STEST(1,stmp,stmp2,ssize1(kn),Sfac,Kprint)
          CASE (3)
            ! 3. SDSDOT
            stmp = SDSDOT(N,sb,sx,INCx,sy,INCy)
            CALL STEST(1,stmp,st7b(kn,ki),ssize3(kn),Sfac,Kprint)
          CASE (4)
            ! 4. DDOT
            dtmp = DDOT(N,dx,INCx,dy,INCy)
            CALL DTEST(1,dtmp,dt7(kn,ki),dsize1(kn),Dfac,Kprint)
          CASE (5)
            ! 5. DQDOTI
            !     DQDOTI AND DQDOTA ARE SUPPOSED TO USE EXTENDED
            !     PRECISION ARITHMETIC INTERNALLY.
            !     SET MODE = 1 OR 2 TO DISTINGUISH TESTS OF DQDOTI OR DQDOTA
            !     IN THE DIAGNOSTIC OUTPUT.
            !
            MODe = 1
            dtmp = DQDOTI(N,db,qc,dx2,INCx,dy2,INCy)
            CALL DTEST(1,dtmp,dt2(kn,ki,1),dt2(kn,ki,1),Dqfac,Kprint)
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
            dtmp = DQDOTI(N,db,qc,dx2,INCx,dy2,INCy)
            CALL DTEST(1,dtmp,dt2(kn,ki,1),dt2(kn,ki,1),Dqfac,Kprint)
            MODe = 2
            dtmp = DQDOTA(N,-db,qc,dx2,INCx,dy2,INCy)
            CALL DTEST(1,dtmp,dt2(kn,ki,2),dt2(kn,ki,2),Dqfac,Kprint)
          CASE (7)
            ! 7. CDOTC
            ctmp = CDOTC(N,cx,INCx,cy,INCy)
            CALL CTEST(1,ctmp,ct6(kn,ki),csize1(kn),Sfac,Kprint)
          CASE (8)
            ! 8. CDOTU
            ctmp = CDOTU(N,cx,INCx,cy,INCy)
            CALL CTEST(1,ctmp,ct7(kn,ki),csize1(kn),Sfac,Kprint)
          CASE (9)
            ! 9. SAXPY
            CALL SAXPY(N,sa,sx,INCx,sy,INCy)
            DO j = 1, leny
              sty(j) = REAL( dt8(j,kn,ki), 4 )
            ENDDO
            CALL STEST(leny,sy,sty,ssize2(:,ksize),Sfac,Kprint)
          CASE (10)
            ! 10. DAXPY
            CALL DAXPY(N,da,dx,INCx,dy,INCy)
            CALL DTEST(leny,dy,dt8(1,kn,ki),dsize2(:,ksize),Dfac,Kprint)
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
            CALL STEST(lenx,sx,stx,ssize2(:,ksize),Sfac,Kprint)
            CALL STEST(leny,sy,sty,ssize2(:,ksize),Sfac,Kprint)
          CASE (15)
            ! 15. DROT
            DO i = 1, 7
              dx(i) = dx1(i)
              dy(i) = dy1(i)
            ENDDO
            CALL DROT(N,dx,INCx,dy,INCy,dc,ds)
            CALL DTEST(lenx,dx,dt9x(1,kn,ki),dsize2(:,ksize),Dfac,Kprint)
            CALL DTEST(leny,dy,dt9y(1,kn,ki),dsize2(:,ksize),Dfac,Kprint)
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
            CALL STEST(leny,sy,sty,ssize2,1.,Kprint)
          CASE (21)
            ! 21. DCOPY
            CALL DCOPY(N,dx,INCx,dy,INCy)
            CALL DTEST(leny,dy,dt10y(1,kn,ki),dsize2,1.D0,Kprint)
          CASE (22)
            ! 22. CCOPY
            CALL CCOPY(N,cx,INCx,cy,INCy)
            ctmp2 = CMPLX(ssize2(1:8:2,1),ssize2(2:8:2,1))
            CALL CTEST(leny,cy,ct10y(:,kn,ki),ctmp2,1.,Kprint)
          CASE (23)
            ! 23. SSWAP
            CALL SSWAP(N,sx,INCx,sy,INCy)
            DO i = 1, 7
              stx(i) = REAL( dt10x(i,kn,ki), 4 )
              sty(i) = REAL( dt10y(i,kn,ki), 4 )
            ENDDO
            CALL STEST(lenx,sx,stx,ssize2,1.,Kprint)
            CALL STEST(leny,sy,sty,ssize2,1.,Kprint)
          CASE (24)
            ! 24. DSWAP
            CALL DSWAP(N,dx,INCx,dy,INCy)
            CALL DTEST(lenx,dx,dt10x(1,kn,ki),dsize2,1.D0,Kprint)
            CALL DTEST(leny,dy,dt10y(1,kn,ki),dsize2,1.D0,Kprint)
          CASE (25)
            ! 25. CSWAP
            CALL CSWAP(N,cx,INCx,cy,INCy)
            ctmp2 = CMPLX(ssize2(1:8:2,1),ssize2(2:8:2,1))
            CALL CTEST(lenx,cx,ct10x(:,kn,ki),ctmp2,1.,Kprint)
            CALL CTEST(leny,cy,ct10y(:,kn,ki),ctmp2,1.,Kprint)
          CASE DEFAULT
            !                                                              1. SDOT
            stmp = SDOT(N,sx,INCx,sy,INCy)
            stmp2 = REAL(dt7(kn,ki))
            CALL STEST(1,stmp,stmp2,ssize1(kn),Sfac,Kprint)
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

    INTEGER i, Len, Kprint, id
    INTEGER Icomp(*), Itrue(*)
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
          WRITE (NPRint,99003) ICAse, N, INCx, INCy, MODe, i, Icomp(i), Itrue(i), id
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

    INTEGER i, Len, Kprint
    REAL Scomp(*), Strue(*), Ssize(*), Sfac, sd, R1MACH
    REAL :: releps = 0.0E0
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

    INTEGER i, Len, Kprint
    REAL(8) :: Dcomp(*), Dtrue(*), Dsize(*), Dfac, dd, D1MACH
    REAL(8) :: releps = 0.0D0
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

    INTEGER i, Len, Kprint
    COMPLEX :: Ccomp(*), Ctrue(*), Csize(*)
    REAL :: Cfac, dd, R1MACH, CABS1
    REAL :: releps = 0.0
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
