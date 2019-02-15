!DECK TEST10
PROGRAM TEST10
  IMPLICIT NONE
  INTEGER I1MACH
  !***BEGIN PROLOGUE  TEST10
  !***PURPOSE  Driver for testing SLATEC subprograms.
  !***LIBRARY   SLATEC
  !***CATEGORY  C7A, C10A4, C10B4, C10D
  !***TYPE      DOUBLE PRECISION (TEST09-S, TEST10-D)
  !***KEYWORDS  QUICK CHECK DRIVER
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
  !               ZABS, ZAIRY, ZBESH, ZBESI, ZBESJ, ZBESK, ZBESY, ZBIRY,
  !               ZDIV, ZEXP, ZSQRT, ZUOIK, ZWRSK
  !
  !***REFERENCES  Fong, Kirby W., Jefferson, Thomas H., Suyehiro,
  !                 Tokihiko, Walton, Lee, Guidelines to the SLATEC Common
  !                 Mathematical Library, March 21, 1989.
  !***ROUTINES CALLED  I1MACH, XERMAX, XSETF, XSETUN, ZQCAI, ZQCBH, ZQCBI,
  !                    ZQCBJ, ZQCBK, ZQCBY
  !***REVISION HISTORY  (YYMMDD)
  !   910411  DATE WRITTEN
  !   920128  Category corrected.  (WRB)
  !***END PROLOGUE  TEST10
  INTEGER ipass, kprint, lin, lun, nfail
  !***FIRST EXECUTABLE STATEMENT  TEST10
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
  !     Test Double Precision Complex Bessel Functions.
  !
  CALL ZQCAI(lun,kprint,ipass)
  IF ( ipass==0 ) nfail = nfail + 1
  CALL ZQCBH(lun,kprint,ipass)
  IF ( ipass==0 ) nfail = nfail + 1
  CALL ZQCBI(lun,kprint,ipass)
  IF ( ipass==0 ) nfail = nfail + 1
  CALL ZQCBJ(lun,kprint,ipass)
  IF ( ipass==0 ) nfail = nfail + 1
  CALL ZQCBK(lun,kprint,ipass)
  IF ( ipass==0 ) nfail = nfail + 1
  CALL ZQCBY(lun,kprint,ipass)
  IF ( ipass==0 ) nfail = nfail + 1
  !
  !     Write PASS or FAIL message
  !
  IF ( nfail==0 ) THEN
    WRITE (lun,99001)
    99001   FORMAT (/' --------------TEST10 PASSED ALL TESTS----------------')
  ELSE
    WRITE (lun,99002) nfail
    99002   FORMAT (/' ************* WARNING -- ',I5,&
      ' TEST(S) FAILED IN PROGRAM TEST10  *************')
  ENDIF
  STOP
END PROGRAM TEST10
