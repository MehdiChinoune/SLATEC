*DECK QC6A
      SUBROUTINE QC6A (LUN, KPRINT, IPASS)
C***BEGIN PROLOGUE  QC6A
C***PURPOSE  Test subroutine AAAAAA.
C***LIBRARY   SLATEC
C***TYPE      ALL (QC6A-A)
C***AUTHOR  Boland, W. Robert, (LANL)
C***DESCRIPTION
C
C *Usage:
C
C        INTEGER  LUN, KPRINT, IPASS
C
C        CALL  QC6A (LUN, KPRINT, IPASS)
C
C *Arguments:
C
C     LUN   :IN  is the unit number to which output is to be written.
C
C     KPRINT:IN  controls the amount of output, as specified in the
C                SLATEC Guidelines.
C
C     IPASS:OUT  indicates whether the test passed or failed.
C                A value of one is good, indicating no failures.
C
C *Description:
C
C   This routine tests the SLATEC routine AAAAAA to see if the version
C   number in the SLATEC library source is the same as the quick check
C   version number.
C
C***ROUTINES CALLED  AAAAAA
C***REVISION HISTORY  (YYMMDD)
C   890713  DATE WRITTEN
C   921215  Updated for Version 4.0.  (WRB)
C   930701  Updated for Version 4.1.  (WRB)
C***END PROLOGUE  QC6A
C
C*Internal Notes:
C
C     Data set-up is done via a PARAMETER statement.
C
C**End
C
C  Declare arguments.
C
      INTEGER  LUN, KPRINT, IPASS
C
C  DECLARE VARIABLES.
C
      CHARACTER * 16 VER, VERSN
      PARAMETER  (VERSN = ' 4.1')
C
C***FIRST EXECUTABLE STATEMENT  QC6A
      IF (KPRINT.GE.3) WRITE (LUN, 9000)
      CALL AAAAAA (VER)
      IF (VER .EQ. VERSN) THEN
         IPASS = 1
         IF (KPRINT .GE. 3) THEN
            WRITE (LUN, 9010)
            WRITE (LUN, 9020) VER
         ENDIF
      ELSE
         IPASS = 0
         IF (KPRINT .GE. 3) WRITE (LUN, 9010)
         IF (KPRINT .GE. 2) WRITE (LUN, 9030) VER, VERSN
      ENDIF
C
C     Terminate.
C
      IF (KPRINT.GE.2 .AND. IPASS.EQ.1) WRITE (LUN, 90000)
      IF (KPRINT.GE.1 .AND. IPASS.EQ.0) WRITE (LUN, 90010)
      RETURN
C
C     Formats.
C
 9000 FORMAT ('1' // ' CODE TO TEST SLATEC ROUTINE AAAAAA')
 9010 FORMAT (/ ' QC6A RESULTS')
 9020 FORMAT (' *** Passed -- version number = ', A16)
 9030 FORMAT (' *** Failed -- version number from AAAAAA = ', A16,
     +        ' but expected version number = ', A16)
90000 FORMAT(/' ************QC6A   PASSED  ALL TESTS ****************')
90010 FORMAT(/' ************QC6A   FAILED SOME TESTS ****************')
C------------- LAST LINE OF QC6A FOLLOWS -----------------------------
      END
