*DECK TEST26
      PROGRAM TEST26
C***BEGIN PROLOGUE  TEST26
C***PURPOSE  Driver for testing SLATEC subprograms.
C***LIBRARY   SLATEC
C***CATEGORY  D2A4, D2B4
C***TYPE      DOUBLE PRECISION (TEST25-S, TEST26-D)
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
C       Double precision SLAP subprograms
C
C***REFERENCES  Fong, Kirby W., Jefferson, Thomas H., Suyehiro,
C                 Tokihiko, Walton, Lee, Guidelines to the SLATEC Common
C                 Mathematical Library, March 21, 1989.
C***ROUTINES CALLED  DLAPQC, I1MACH, XERMAX, XSETF, XSETUN
C***REVISION HISTORY  (YYMMDD)
C   920401  DATE WRITTEN
C   920511  Added complete declaration section.  (WRB)
C***END PROLOGUE  TEST26
C     .. Local Scalars ..
      INTEGER IPASS, KPRINT, LIN, LUN, NFAIL
C     .. External Functions ..
      INTEGER I1MACH
      EXTERNAL I1MACH
C     .. External Subroutines ..
      EXTERNAL DLAPQC, XERMAX, XSETF, XSETUN
C***FIRST EXECUTABLE STATEMENT  TEST26
      LUN = I1MACH(2)
      LIN = I1MACH(1)
      NFAIL = 0
C
C     Read KPRINT parameter
C
      READ (LIN, '(I1)') KPRINT
      CALL XSETUN(LUN)
      IF (KPRINT .LE. 1) THEN
         CALL XSETF(0)
      ELSE
         CALL XSETF(1)
      ENDIF
      CALL XERMAX(1000)
C
C     Test SLAP (double precision)
C
      CALL DLAPQC(LUN,KPRINT,IPASS)
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
 9000 FORMAT (/' --------------TEST26 PASSED ALL TESTS----------------')
 9010 FORMAT (/' ************* WARNING -- ', I5,
     1        ' TEST(S) FAILED IN PROGRAM TEST26 *************')
      END
