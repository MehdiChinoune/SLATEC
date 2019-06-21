MODULE TEST04_MOD
  USE service, ONLY : SP, DP
  IMPLICIT NONE

CONTAINS
  !** CFNCK
  SUBROUTINE CFNCK(Lun,Kprint,Ipass)
    !> Quick check for the complex Fullerton special functions.
    !***
    ! **Library:**   SLATEC
    !***
    ! **Keywords:**  QUICK CHECK
    !***
    ! **Author:**  Boland, W. Robert, (LANL)
    !           Chow, Jeff, (LANL)
    !           Rivera, Shawn, (LANL)
    !***
    ! **Description:**
    !
    !     This subroutine does a quick check for the complex
    !     routines in the Fullerton special function library.
    !
    !     Parameter list-
    !
    !     LUN      input integer value to designate the external
    !              device unit for message output
    !     KPRINT   input integer value to specify amount of
    !              printing to be done by quick check
    !     IPASS    output value indicating whether tests passed or
    !              failed
    !
    !***
    ! **Routines called:**  C0LGMC, CACOS, CACOSH, CASIN, CASINH, CATAN,
    !                    CATAN2, CATANH, CBETA, CCBRT, CCOSH, CCOT, CEXPRL,
    !                    CGAMMA, CGAMR, CLBETA, CLNGAM, CLNREL, CLOG10,
    !                    CPSI, CSINH, CTAN, CTANH, R1MACH

    !* REVISION HISTORY  (YYMMDD)
    !   800901  DATE WRITTEN
    !   891115  REVISION DATE from Version 3.2
    !   891120  Checks of remainder of FNLIB routines added and code
    !           reorganized.  (WRB)
    !   900330  Prologue converted to Version 4.0 format.  (BAB)
    !   900727  Added EXTERNAL statement.  (WRB)
    USE slatec, ONLY : C0LGMC, CATAN2, CBETA, CCBRT, CCOT, CEXPRL, CGAMMA, CGAMR, &
      CLBETA, CLNGAM, CLNREL, CLOG10, CPSI, R1MACH
    INTEGER :: i, Lun, Kprint, Ipass
    REAL(SP) :: errmax, errtol, abserr, relerr
    COMPLEX(SP) :: w(48)
    !
    !     Constants to be used
    !
    COMPLEX(SP), PARAMETER :: c1 = (1._SP,0._SP), ci = (0._SP,1._SP)
    REAL(SP), PARAMETER :: sqrt2  = .14142135623730950488E1_SP
    REAL(SP), PARAMETER :: sqrt3  = .17320508075688772935E1_SP
    REAL(SP), PARAMETER :: pi = 3.14159265358979323846_SP
    !
    !     Complex values through different calculations are stored in C(*)
    !
    COMPLEX(SP), PARAMETER :: c(28) = [ (.121699028117870E1_SP,.326091563038355_SP), &
      (.866025403784438_SP,.500000000000000_SP), &
      (-.785398163397449_SP,-.658478948462413_SP), &
      (.785398163397449_SP,-.658478948462413_SP), &
      (.313314753080534E-1_SP,.541264220944095E-1_SP), &
      (-.313314753080534E-1_SP,.541264220944095E-1_SP), &
      (.183048772171245E1_SP,.000000000000000_SP), &
      (-.757236713834364E-1_SP,-.961745759068982_SP), &
      (.150514997831990_SP,-.341094088460459_SP), &
      (.301029995663980_SP,.227396058973639_SP), &
      (.000000000000000_SP,.636619772367581_SP), &
      (.137802461354738E1_SP,.909330673631480_SP), &
      (.303123109082158E-1_SP,-.244978663126864_SP), &
      (.693147180559947_SP,.523598775598298_SP), &
      (.365427607174532E-1_SP,-.612881308922810E-1_SP), &
      (.896860330225849E-2_SP,.244804656578857E-1_SP), &
      (.177245385090552E1_SP,.000000000000000_SP), &
      (.300694617260656_SP,-.424967879433124_SP), &
      (.110951302025214E1_SP,-.156806064476794E1_SP), &
      (.183074439659052E1_SP,.569607641036682_SP), &
      (-.340863758923258E1_SP,.142127515954291E1_SP), &
      (-.156059525546301E1_SP,.152533527872833E1_SP), &
      (-.211272372936533_SP,-.765528316537801_SP), &
      (.380273164249058E-1_SP,-.286343074460341_SP), &
      (-.268079774264798E1_SP,.130151697855085E1_SP), &
      (-.164841998888369E1_SP,.785398163397448_SP), &
      (-.196351002602143E1_SP,.000000000000000_SP), &
      (.161278484461574E1_SP,.147079632679497E1_SP) ]
    !* FIRST EXECUTABLE STATEMENT  CFNCK
    !
    !     Compute functional values
    !
    !     Exercise routines in Category C2.
    !
    w(1) = CCBRT(sqrt2*(1._SP+ci))
    w(2) = CCBRT(ci)
    !
    !     Exercise routines in Category C4A.
    !
    w(3) = CATAN2(-.5_SP-.866025403784438_SP*ci,c1)
    w(4) = CATAN2(.5_SP-.866025403784438_SP*ci,c1)
    w(5) = CATAN2(.3125E-1_SP+.541265877365273E-1_SP*ci,c1)
    w(6) = CATAN2(-.3125E-1_SP+.541265877365273E-1_SP*ci,c1)
    w(7) = CCOT(.5_SP+0._SP*ci)
    w(8) = CCOT(-1._SP+.5_SP*pi*ci)
    !
    !     Exercise routines in Category C4B.
    !
    w(9) = CLOG10(1._SP-ci)
    w(10) = CLOG10(sqrt3+ci)
    w(11) = CEXPRL(pi*ci)
    w(12) = CEXPRL(1._SP+ci)
    w(13) = CLNREL(-.25_SP*ci)
    w(14) = CLNREL(sqrt3-1._SP+ci)
    !
    !     Exercise routines in Category C7A.
    !
    w(15) = C0LGMC(.5_SP+.5_SP*ci)
    w(16) = C0LGMC(1._SP-1._SP*ci)
    w(17) = CGAMMA(.5_SP+0._SP*ci)
    w(18) = CGAMMA(.5_SP+ci)
    w(19) = CGAMR(.5_SP-ci)
    w(20) = CGAMR(1._SP+ci)
    w(21) = CLNGAM(1.1_SP+3.2_SP*ci)
    w(22) = CLNGAM(1.9_SP+2.4_SP*ci)
    !
    !     Exercise routines in Category C7B.
    !
    w(23) = CBETA(1._SP+ci,1._SP+ci)
    w(24) = CBETA(2._SP-ci,.5_SP+ci)
    w(25) = CLBETA(2._SP+ci,1._SP-2._SP*ci)
    w(26) = CLBETA(1._SP-ci,2._SP+ci)
    !
    !     Exercise routines in Category C7C.
    !
    w(27) = CPSI(.5_SP+0._SP*ci)
    w(28) = CPSI(1._SP+5._SP*ci)
    !
    !     Check for possible errors
    !
    errmax = R1MACH(4)
    errtol = SQRT(errmax)
    DO i = 1, 28
      abserr = ABS(c(i)-w(i))
      relerr = abserr/ABS(c(i))
      errmax = MAX(relerr,errmax)
      IF( relerr>errtol .AND. Kprint>=2 ) WRITE (Lun,99001) i, relerr, abserr
      99001 FORMAT (' For I  = ',I3,'  test fails with RELERR  = ',E38.30,&
        '  and ABSERR  = ',E38.30)
    END DO
    Ipass = 0
    IF( errmax<=errtol ) Ipass = 1
    IF( Ipass/=0 .AND. Kprint>=2 ) WRITE (Lun,99002)
    99002 FORMAT (' Complex Fullerton special function routines o.k.')
    RETURN
  END SUBROUTINE CFNCK
END MODULE TEST04_MOD
!** TEST04
PROGRAM TEST04
  USE TEST04_MOD, ONLY : CFNCK
  USE slatec, ONLY : I1MACH, XSETF, XSETUN, XERMAX
  USE common_mod, ONLY : GET_ARGUMENT
  IMPLICIT NONE
  !> Driver for testing SLATEC subprograms
  !***
  ! **Library:**   SLATEC
  !***
  ! **Category:**  C
  !***
  ! **Type:**      COMPLEX (TEST02-S, TEST03-D, TEST04-C)
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
  !        complex Fullerton routines
  !
  !***
  ! **References:**  Kirby W. Fong, Thomas H. Jefferson, Tokihiko Suyehiro
  !                 and Lee Walton, Guide to the SLATEC Common Mathema-
  !                 tical Library, April 10, 1990.
  !***
  ! **Routines called:**  CFNCK, I1MACH, XERMAX, XSETF, XSETUN

  !* REVISION HISTORY  (YYMMDD)
  !   890618  DATE WRITTEN
  !   890618  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900524  Cosmetic changes to code.  (WRB)
  INTEGER :: ipass, kprint, lin, lun, nfail
  !* FIRST EXECUTABLE STATEMENT  TEST04
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
  !     Test complex Fullerton routines
  !
  CALL CFNCK(lun,kprint,ipass)
  IF( ipass==0 ) nfail = nfail + 1
  !
  !     Write PASS or FAIL message
  !
  IF( nfail==0 ) THEN
    WRITE (lun,99001)
    99001 FORMAT (/' --------------TEST04 PASSED ALL TESTS----------------')
  ELSE
    WRITE (lun,99002) nfail
    99002 FORMAT (/' ************* WARNING -- ',I5,&
      ' TEST(S) FAILED IN PROGRAM TEST04  *************')
  END IF
  STOP
END PROGRAM TEST04
