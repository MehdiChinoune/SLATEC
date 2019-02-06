*DECK TEST09
      PROGRAM TEST09
C***BEGIN PROLOGUE  TEST09
C***PURPOSE  Driver for testing SLATEC subprograms.
C***LIBRARY   SLATEC
C***CATEGORY  C7A, C10A4, C10B4, C10D
C***TYPE      SINGLE PRECISION (TEST09-S, TEST10-D)
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
C         CAIRY, CBESH, CBESI, CBESJ, CBESK, CBESY, CBIRY, CUOIK, CWRSK
C
C***REFERENCES  Fong, Kirby W., Jefferson, Thomas H., Suyehiro,
C                 Tokihiko, Walton, Lee, Guidelines to the SLATEC Common
C                 Mathematical Library, March 21, 1989.
C***ROUTINES CALLED  CQCAI, CQCBH, CQCBI, CQCBJ, CQCBK, CQCBY, I1MACH,
C                    XERMAX, XSETF, XSETUN
C***REVISION HISTORY  (YYMMDD)
C   910411  DATE WRITTEN
C   920128  Category corrected.  (WRB)
C***END PROLOGUE  TEST09
      INTEGER IPASS, KPRINT, LIN, LUN, NFAIL
C***FIRST EXECUTABLE STATEMENT  TEST09
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
C     Test Single Precision Complex Bessel Functions.
C
      CALL CQCAI(LUN,KPRINT,IPASS)
      IF (IPASS .EQ. 0) NFAIL = NFAIL + 1
      CALL CQCBH(LUN,KPRINT,IPASS)
      IF (IPASS .EQ. 0) NFAIL = NFAIL + 1
      CALL CQCBI(LUN,KPRINT,IPASS)
      IF (IPASS .EQ. 0) NFAIL = NFAIL + 1
      CALL CQCBJ(LUN,KPRINT,IPASS)
      IF (IPASS .EQ. 0) NFAIL = NFAIL + 1
      CALL CQCBK(LUN,KPRINT,IPASS)
      IF (IPASS .EQ. 0) NFAIL = NFAIL + 1
      CALL CQCBY(LUN,KPRINT,IPASS)
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
 9000 FORMAT (/' --------------TEST09 PASSED ALL TESTS----------------')
 9010 FORMAT (/' ************* WARNING -- ', I5,
     1        ' TEST(S) FAILED IN PROGRAM TEST09  *************')
      END
