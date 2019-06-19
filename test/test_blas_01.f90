MODULE TEST17_MOD
  USE service, ONLY : SP, DP
  IMPLICIT NONE
  INTEGER :: nprint_com, icase_com, n_com, incx_com, incy_com, mode_com
  LOGICAL :: PASs

CONTAINS
  !** BLACHK
  SUBROUTINE BLACHK(Lun,Kprint,Ipass)
    !> Quick check for Basic Linear Algebra Subprograms.
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

    INTEGER :: Lun, Ipass, Kprint
    REAL, PARAMETER :: sfac = .625E-1
    REAL(DP), PARAMETER :: dfac = .625D-1, dqfac = 0.625D-1
    INTEGER, PARAMETER :: jtest(38) = [ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, &
      1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 ]
    !* FIRST EXECUTABLE STATEMENT  BLACHK
    nprint_com = Lun
    Ipass = 1
    !
    IF( Kprint>=2 ) WRITE (nprint_com,99001)
    99001 FORMAT ('1','QUICK CHECK OF 38 BASIC LINEAR ALGEBRA SUBROUTINES'/)
    DO icase_com = 1, 38
      IF( jtest(icase_com)/=0 ) THEN
        CALL HEADER(Kprint)
        !
        !         INITIALIZE  PASS, INCX, INCY, AND MODE FOR A NEW CASE.
        !         THE VALUE 9999 FOR INCX, INCY OR MODE WILL APPEAR IN THE
        !         DETAILED  OUTPUT, IF ANY, FOR CASES THAT DO NOT INVOLVE
        !         THESE PARAMETERS.
        !
        PASs = .TRUE.
        incx_com = 9999
        incy_com = 9999
        mode_com = 9999
        SELECT CASE (icase_com)
          CASE (5,6)
            ! ICASE =  1-11, 14-15, OR 18-25
            CALL CHECK2(sfac,dfac,dqfac,Kprint)
        END SELECT
        !                                                  PRINT
        IF( Kprint>=2 .AND. PASs ) WRITE (nprint_com,99002)
        99002 FORMAT ('+',39X,'PASS')
        IF( .NOT. PASs ) Ipass = 0
      END IF
    END DO
    IF( Kprint>=2 .AND. Ipass==1 ) WRITE (nprint_com,99003)
    99003 FORMAT (/' ****************BLAS PASSED ALL TESTS****************')
    IF( Kprint>=1 .AND. Ipass==0 ) WRITE (nprint_com,99004)
    99004 FORMAT (/' ****************BLAS FAILED SOME TESTS***************')
    RETURN
  END SUBROUTINE BLACHK
  !** HEADER
  SUBROUTINE HEADER(Kprint)
    !> Print header for BLAS quick checks.
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
    IF( Kprint>=2 ) WRITE (nprint_com,99001) icase_com, l(icase_com)
    !
    99001 FORMAT (' Test of subprogram number',I3,2X,A)
    RETURN
  END SUBROUTINE HEADER
  !** CHECK2
  SUBROUTINE CHECK2(Sfac,Dfac,Dqfac,Kprint)
    !> (UNKNOWN)
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
    USE linear, ONLY : DQDOTA, DQDOTI
    INTEGER :: i, ki, kn, ksize, lenx, leny, mx, my, Kprint, qc_i(30)
    REAL(SP) :: Sfac
    REAL(SP) :: sx(7), sy(7)
    REAL(DP) :: dx(7), dy(7), dtmp(1)
    REAL(DP) :: Dfac, Dqfac
    !
    COMPLEX(SP) :: cx(7), cy(7)
    REAL(DP), PARAMETER :: db = .25D0
    INTEGER, PARAMETER :: incxs(4) = [ 1, 2, -2, -1 ]
    INTEGER, PARAMETER :: incys(4) = [ 1, -2, 1, -2 ]
    INTEGER, PARAMETER :: lens(4,2) = RESHAPE( [1, 1, 2, 4, 1, 1, 3, 7], [4,2] )
    INTEGER, PARAMETER :: ns(4) = [ 0, 1, 2, 4 ]
    REAL(DP), PARAMETER :: dx1(7) = [ .6D0, .1D0, -.5D0, .8D0, .9D0, -.3D0, -.4D0 ]
    REAL(DP), PARAMETER :: dy1(7) = [ .5D0, -.9D0, .3D0, .7D0, -.6D0, .2D0, .8D0 ]
    REAL(DP), PARAMETER :: dx2(7) = [ 1.D0, .01D0, .02D0, 1.D0, .06D0, 2.D0, 1.D0 ]
    REAL(DP), PARAMETER :: dy2(7) = [ 1.D0, .04D0, -.03D0, -1.D0, .05D0, 3.D0, -1.D0 ]
    !            THE TERMS D11(3,2) AND D11(4,2) WILL BE SET BY
    !            COMPUTATION AT RUN TIME.
    COMPLEX(SP), PARAMETER :: cx1(7) = [ (.7,-.8), (-.4,-.7), (-.1,-.9), (.2,-.8), &
      (-.9,-.4), (.1,.4), (-.6,.6) ]
    COMPLEX(SP), PARAMETER :: cy1(7) = [ (.6,-.6), (-.9,.5), (.7,-.6), (.1,-.5), &
      (-.1,-.2), (-.5,-.3), (.8,-.7) ]
    !
    !                             FOR DQDOTI AND DQDOTA
    !
    REAL(DP), PARAMETER :: dt2(4,4,2) = RESHAPE( [ 0.25D0, 1.25D0, 1.2504D0, 0.2498D0, &
      0.25D0, 1.25D0, 0.24D0, 0.2492D0, 0.25D0, 1.25D0, 0.31D0, 0.2518D0, &
      0.25D0, 1.25D0, 1.2497D0, 0.2507D0, 0.D0, 2.D0, 2.0008D0, -.0004D0, &
      0.D0, 2.D0, -.02D0, -.0016D0, 0.D0, 2.D0, .12D0, .0036D0, &
      0.D0, 2.D0, 1.9994D0, .0014D0 ], [4,4,2] )
    !* FIRST EXECUTABLE STATEMENT  CHECK2
    DO ki = 1, 4
      incx_com = incxs(ki)
      incy_com = incys(ki)
      mx = ABS(incx_com)
      my = ABS(incy_com)
      !
      DO kn = 1, 4
        n_com = ns(kn)
        ksize = MIN(2,kn)
        lenx = lens(kn,mx)
        leny = lens(kn,my)
        ! INITIALIZE ALL ARGUMENT ARRAYS.
        DO i = 1, 7
          sx(i) = REAL( dx1(i), SP )
          sy(i) = REAL( dy1(i), SP )
          dx(i) = dx1(i)
          dy(i) = dy1(i)
          cx(i) = cx1(i)
          cy(i) = cy1(i)
        END DO
        !
        ! BRANCH TO SELECT SUBPROGRAM TO BE TESTED.
        !
        SELECT CASE (icase_com)
          CASE (5)
            ! 5. DQDOTI
            !     DQDOTI AND DQDOTA ARE SUPPOSED TO USE EXTENDED
            !     PRECISION ARITHMETIC INTERNALLY.
            !     SET MODE = 1 OR 2 TO DISTINGUISH TESTS OF DQDOTI OR DQDOTA
            !     IN THE DIAGNOSTIC OUTPUT.
            !
            mode_com = 1
            dtmp = DQDOTI(n_com,db,qc_i,dx2,incx_com,dy2,incy_com)
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
            mode_com = 1
            dtmp = DQDOTI(n_com,db,qc_i,dx2,incx_com,dy2,incy_com)
            CALL DTEST(1,dtmp,dt2(kn,ki,1),dt2(kn,ki,1),Dqfac,Kprint)
            mode_com = 2
            dtmp = DQDOTA(n_com,-db,qc_i,dx2,incx_com,dy2,incy_com)
            CALL DTEST(1,dtmp,dt2(kn,ki,2),dt2(kn,ki,2),Dqfac,Kprint)
          CASE (12,13,16,17)
            GOTO 100
        END SELECT
        !
        !
        !
      END DO
    END DO
    RETURN
    !                 THE FOLLOWING STOP SHOULD NEVER BE REACHED.
    100  STOP
  END SUBROUTINE CHECK2
  !** ITEST
  SUBROUTINE ITEST(Leng,Icomp,Itrue,Kprint)
    !> Compare arrays ICOMP and ITRUE.
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

    INTEGER :: i, Leng, Kprint, id
    INTEGER :: Icomp(*), Itrue(*)
    !* FIRST EXECUTABLE STATEMENT  ITEST
    DO i = 1, Leng
      IF( Icomp(i)/=Itrue(i) ) THEN
        !
        !         Here ICOMP(I) is not equal to ITRUE(I).
        !
        IF( PASs ) THEN
          !
          !           Print FAIL message and header.
          !
          PASs = .FALSE.
          IF( Kprint>=3 ) THEN
            WRITE (nprint_com,99001)
            99001 FORMAT ('+',39X,'FAIL')
            WRITE (nprint_com,99002)
            99002 FORMAT ('0CASE  N INCX INCY MODE  I',29X,'COMP(I)',29X,'TRUE(I)',&
              2X,'DIFFERENCE'/1X)
          END IF
        END IF
        IF( Kprint>=3 ) THEN
          id = Icomp(i) - Itrue(i)
          WRITE (nprint_com,99003) icase_com, n_com, incx_com, incy_com, mode_com, &
            i, Icomp(i), Itrue(i), id
          99003 FORMAT (1X,I4,I3,3I5,I3,2I36,I12)
        END IF
      END IF
    END DO
    RETURN
  END SUBROUTINE ITEST
  !** STEST
  SUBROUTINE STEST(Leng,Scomp,Strue,Ssize,Sfac,Kprint)
    !> Compare arrays SCOMP and STRUE.
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
    USE slatec, ONLY : R1MACH
    INTEGER :: i, Leng, Kprint
    REAL(SP) :: Scomp(*), Strue(*), Ssize(*), Sfac, sd
    REAL(SP) :: releps = 0.0E0
    !* FIRST EXECUTABLE STATEMENT  STEST
    IF( releps==0.0E0 ) releps = R1MACH(4)
    DO i = 1, Leng
      sd = ABS(Scomp(i)-Strue(i))
      IF( Sfac*sd>ABS(Ssize(i))*releps ) THEN
        !
        !         Here SCOMP(I) is not close to STRUE(I).
        !
        IF( PASs ) THEN
          !
          !           Print FAIL message and header.
          !
          PASs = .FALSE.
          IF( Kprint>=3 ) THEN
            WRITE (nprint_com,99001)
            99001 FORMAT ('+',39X,'FAIL')
            WRITE (nprint_com,99002)
            99002 FORMAT ('0CASE  N INCX INCY MODE  I',29X,'COMP(I)',29X,'TRUE(I)',&
              2X,'DIFFERENCE',5X,'SIZE(I)'/1X)
          END IF
        END IF
        IF( Kprint>=3 ) WRITE (nprint_com,99003) icase_com, n_com, incx_com, &
          incy_com, mode_com, i, Scomp(i), Strue(i), sd, Ssize(i)
        99003 FORMAT (1X,I4,I3,3I5,I3,2E36.8,2E12.4)
      END IF
    END DO
    RETURN
  END SUBROUTINE STEST
  !** DTEST
  SUBROUTINE DTEST(Leng,Dcomp,Dtrue,Dsize,Dfac,Kprint)
    !> Compare arrays DCOMP and DTRUE.
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
    USE slatec, ONLY : D1MACH
    INTEGER :: i, Leng, Kprint
    REAL(DP) :: Dcomp(*), Dtrue(*), Dsize(*), Dfac, dd
    REAL(DP) :: releps = 0.0D0
    !* FIRST EXECUTABLE STATEMENT  DTEST
    IF( releps==0.0D0 ) releps = D1MACH(4)
    DO i = 1, Leng
      dd = ABS(Dcomp(i)-Dtrue(i))
      IF( Dfac*dd>ABS(Dsize(i))*releps ) THEN
        !
        !         Here DCOMP(I) is not close to DTRUE(I).
        !
        IF( PASs ) THEN
          !
          !           Print FAIL message and header.
          !
          PASs = .FALSE.
          IF( Kprint>=3 ) THEN
            WRITE (nprint_com,99001)
            99001 FORMAT ('+',39X,'FAIL')
            WRITE (nprint_com,99002)
            99002 FORMAT ('0CASE  N INCX INCY MODE  I',29X,'COMP(I)',29X,'TRUE(I)',&
              2X,'DIFFERENCE',5X,'SIZE(I)'/1X)
          END IF
        END IF
        IF( Kprint>=3 ) WRITE (nprint_com,99003) icase_com, n_com, incx_com, incy_com, mode_com, &
          i, Dcomp(i), Dtrue(i), dd, Dsize(i)
        99003 FORMAT (1X,I4,I3,3I5,I3,2D36.18,2D12.4)
      END IF
    END DO
    RETURN
  END SUBROUTINE DTEST

  SUBROUTINE CTEST(Leng,Ccomp,Ctrue,Csize,Cfac,Kprint)
    !> Compare arrays DCOMP and DTRUE.
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
    ! **Routines called:**  R1MACH
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
    USE slatec, ONLY : R1MACH
    USE blas, ONLY : SCABS1
    INTEGER :: i, Leng, Kprint
    COMPLEX(SP) :: Ccomp(*), Ctrue(*), Csize(*)
    REAL(SP) :: Cfac, dd
    REAL(SP) :: releps = 0.0
    !* FIRST EXECUTABLE STATEMENT  DTEST
    IF( releps==0.0 ) releps = R1MACH(4)
    DO i = 1, Leng
      dd = SCABS1(Ccomp(i)-Ctrue(i))
      IF( Cfac*dd>ABS(Csize(i))*releps ) THEN
        !
        !         Here DCOMP(I) is not close to DTRUE(I).
        !
        IF( PASs ) THEN
          !
          !           Print FAIL message and header.
          !
          PASs = .FALSE.
          IF( Kprint>=3 ) THEN
            WRITE (nprint_com,99001)
            99001 FORMAT ('+',39X,'FAIL')
            WRITE (nprint_com,99002)
            99002 FORMAT ('0CASE  N INCX INCY MODE  I',29X,'COMP(I)',29X,'TRUE(I)',&
              2X,'DIFFERENCE',5X,'SIZE(I)'/1X)
          END IF
        END IF
        IF( Kprint>=3 ) WRITE (nprint_com,99003) icase_com, n_com, incx_com, incy_com, mode_com, &
          i, Ccomp(i), Ctrue(i), dd, Csize(i)
        99003 FORMAT (1X,I4,I3,3I5,I3,2D36.18,2D12.4)
      END IF
    END DO
    RETURN
  END SUBROUTINE CTEST
END MODULE TEST17_MOD
!** TEST17
PROGRAM TEST17
  USE TEST17_MOD, ONLY : BLACHK
  USE slatec, ONLY : I1MACH, XSETF, XSETUN, XERMAX
  USE common_mod, ONLY : GET_ARGUMENT
  IMPLICIT NONE
  !> Driver for testing SLATEC subprograms
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
  INTEGER :: ipass, kprint, lin, lun, nfail
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
  IF( kprint<=1 ) THEN
    CALL XSETF(0)
  ELSE
    CALL XSETF(1)
  END IF
  !
  !     Test BLAS
  !
  CALL BLACHK(lun,kprint,ipass)
  IF( ipass==0 ) nfail = nfail + 1
  !
  !     Write PASS or FAIL message
  !
  IF( nfail==0 ) THEN
    WRITE (lun,99001)
    99001 FORMAT (/' --------------TEST17 PASSED ALL TESTS----------------')
  ELSE
    WRITE (lun,99002) nfail
    99002 FORMAT (/' ************* WARNING -- ',I5,&
      ' TEST(S) FAILED IN PROGRAM TEST17 *************')
  END IF
  STOP
END PROGRAM TEST17
