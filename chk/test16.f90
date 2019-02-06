!*==TEST16.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK TEST16
      PROGRAM TEST16
      IMPLICIT NONE
!*--TEST165
!*** Start of declarations inserted by SPAG
      INTEGER I1MACH
!*** End of declarations inserted by SPAG
!***BEGIN PROLOGUE  TEST16
!***PURPOSE  Driver for testing SLATEC subprograms
!            DRC3JJ   DRC3JM   DRC6J
!***LIBRARY   SLATEC
!***CATEGORY  C19
!***TYPE      DOUBLE PRECISION (TEST15-S, TEST16-D)
!***KEYWORDS  3J COEFFICIENTS, 3J SYMBOLS, 6J COEFFICIENTS, 6J SYMBOLS,
!             CLEBSCH-GORDAN COEFFICIENTS, QUICK CHECK DRIVER,
!             RACAH COEFFICIENTS, VECTOR ADDITION COEFFICIENTS,
!             WIGNER COEFFICIENTS
!***AUTHOR  SLATEC Common Mathematical Library Committee
!***DESCRIPTION
!
! *Usage:
!     One input data record is required
!         READ (LIN, '(I1)') KPRINT
!
! *Arguments:
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
! *Description:
!     Driver for testing SLATEC subprograms
!        DRC3JJ   DRC3JM   DRC6J
!
!***REFERENCES  Kirby W. Fong, Thomas H. Jefferson, Tokihiko Suyehiro
!                 and Lee Walton, Guide to the SLATEC Common Mathema-
!                 tical Library, April 10, 1990.
!***ROUTINES CALLED  DQC36J, I1MACH, XERMAX, XSETF, XSETUN
!***REVISION HISTORY  (YYMMDD)
!   891130  DATE WRITTEN
!***END PROLOGUE  TEST16
      INTEGER ipass , kprint , lin , lun , nfail
!***FIRST EXECUTABLE STATEMENT  TEST16
      lun = I1MACH(2)
      lin = I1MACH(1)
      nfail = 0
!
!     Read KPRINT parameter
!
      READ (lin,'(I1)') kprint
      CALL XERMAX(1000)
      CALL XSETUN(lun)
      IF ( kprint<=1 ) THEN
        CALL XSETF(0)
      ELSE
        CALL XSETF(1)
      ENDIF
!
!     Test double precision 3J6J routines
!
      CALL DQC36J(lun,kprint,ipass)
      IF ( ipass==0 ) nfail = nfail + 1
!
!     Write PASS or FAIL message
!
      IF ( nfail==0 ) THEN
        WRITE (lun,99001)
99001   FORMAT (/' --------------TEST16 PASSED ALL TESTS----------------')
      ELSE
        WRITE (lun,99002) nfail
99002   FORMAT (/' ************* WARNING -- ',I5,
     &          ' TEST(S) FAILED IN PROGRAM TEST16  *************')
      ENDIF
      STOP
      END PROGRAM TEST16
