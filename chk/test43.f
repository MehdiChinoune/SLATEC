*DECK TEST43
      PROGRAM TEST43
C***BEGIN PROLOGUE  TEST43
C***PURPOSE  Driver for testing SLATEC subprograms
C***LIBRARY   SLATEC
C***CATEGORY  I1
C***TYPE      SINGLE PRECISION (TEST43-S, TEST44-D)
C***KEYWORDS  QUICK CHECK DRIVER
C***AUTHOR  SLATEC Common Mathematical Library Committee
C***DESCRIPTION
C
C *Usage:
C     One input data record is required
C         READ (LIN, '(I1)') KPRINT
C
C *Arguments:
C     KPRINT = 0  Quick checks - No printing.
C                 Driver       - Short pass or fail message printed.
C              1  Quick checks - No message printed for passed tests,
C                                short message printed for failed tests.
C                 Driver       - Short pass or fail message printed.
C              2  Quick checks - Print short message for passed tests,
C                                fuller information for failed tests.
C                 Driver       - Pass or fail message printed.
C              3  Quick checks - Print complete quick check results.
C                 Driver       - Pass or fail message printed.
C
C *Description:
C     Driver for testing SLATEC subprograms
C        DEABM    DEBDF    DERKF    BVSUP
C
C***REFERENCES  Kirby W. Fong, Thomas H. Jefferson, Tokihiko Suyehiro
C                 and Lee Walton, Guide to the SLATEC Common Mathema-
C                 tical Library, April 10, 1990.
C***ROUTINES CALLED  I1MACH, QXABM, QXBDF, QXBVSP, QXRKF, XERMAX, XSETF,
C                    XSETUN
C***REVISION HISTORY  (YYMMDD)
C   890618  DATE WRITTEN
C   890618  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900524  Cosmetic changes to code.  (WRB)
C***END PROLOGUE  TEST43
      INTEGER IPASS, KPRINT, LIN, LUN, NFAIL
C***FIRST EXECUTABLE STATEMENT  TEST43
      LUN = I1MACH(2)
      LIN = I1MACH(1)
      NFAIL = 0
C
C     Read KPRINT parameter
C
      READ (LIN, '(I1)') KPRINT
      CALL XERMAX(1000)
      CALL XSETUN(LUN)
      IF (KPRINT .LE. 1) THEN
         CALL XSETF(0)
      ELSE
         CALL XSETF(1)
      ENDIF
C
C     Test DEABM
C
      CALL QXABM(LUN,KPRINT,IPASS)
      IF (IPASS .EQ. 0) NFAIL = NFAIL + 1
C
C     Test DEBDF
C
      CALL QXBDF(LUN,KPRINT,IPASS)
      IF (IPASS .EQ. 0) NFAIL = NFAIL + 1
C
C     Test DERKF
C
      CALL QXRKF(LUN,KPRINT,IPASS)
      IF (IPASS .EQ. 0) NFAIL = NFAIL + 1
C
C     Test BVSUP
C
      CALL QXBVSP(LUN,KPRINT,IPASS)
      IF (IPASS .EQ. 0) NFAIL = NFAIL + 1
C
C     Write PASS or FAIL message
C
      IF (NFAIL .EQ. 0) THEN
         WRITE (LUN, 9000)
      ELSE
         WRITE (LUN, 9010) NFAIL
      ENDIF
      STOP
 9000 FORMAT (/' --------------TEST43 PASSED ALL TESTS----------------')
 9010 FORMAT (/' ************* WARNING -- ', I5,
     1        ' TEST(S) FAILED IN PROGRAM TEST43 *************')
      END
