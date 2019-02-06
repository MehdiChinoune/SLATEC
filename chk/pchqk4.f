*DECK PCHQK4
      SUBROUTINE PCHQK4 (LUN, KPRINT, IPASS)
C***BEGIN PROLOGUE  PCHQK4
C***PURPOSE  Test the PCHIP monotonicity checker PCHCM.
C***LIBRARY   SLATEC (PCHIP)
C***TYPE      SINGLE PRECISION (PCHQK4-S, DPCHQ4-D)
C***KEYWORDS  PCHIP MONOTONICITY CHECKER QUICK CHECK
C***AUTHOR  Fritsch, F. N., (LLNL)
C***DESCRIPTION
C
C              PCHIP QUICK CHECK NUMBER 4
C
C     TESTS THE MONOTONICITY CHECKER:  PCHCM.
C *Usage:
C
C        INTEGER  LUN, KPRINT, IPASS
C
C        CALL PCHQK4 (LUN, KPRINT, IPASS)
C
C *Arguments:
C
C     LUN   :IN  is the unit number to which output is to be written.
C
C     KPRINT:IN  controls the amount of output, as specified in the
C                SLATEC Guidelines.
C
C     IPASS:OUT  will contain a pass/fail flag.  IPASS=1 is good.
C                IPASS=0 indicates one or more tests failed.
C
C *Description:
C
C   This routine tests a constructed data set with three different
C   INCFD settings and compares with the expected results.  It then
C   runs a special test to check for bug in overall monotonicity found
C   in PCHMC.  Finally, it reverses the data and repeats all tests.
C
C***ROUTINES CALLED  PCHCM
C***REVISION HISTORY  (YYMMDD)
C   890208  DATE WRITTEN
C   890306  Changed LOUT to LUN and added it to call list.  (FNF)
C   890316  Removed DATA statements to suit new quick check standards.
C   890410  Changed PCHMC to PCHCM.
C   890410  Added a SLATEC 4.0 format prologue.
C   900314  Changed name from PCHQK3 to PCHQK4 and improved some output
C           formats.
C   900315  Revised prologue and improved some output formats.  (FNF)
C   900321  Removed IFAIL from call sequence for SLATEC standards and
C           made miscellaneous cosmetic changes.  (FNF)
C   900322  Added declarations so all variables are declared.  (FNF)
C   910708  Minor modifications in use of KPRINT.  (WRB)
C   930317  Improved output formats.  (FNF)
C***END PROLOGUE  PCHQK4
C
C*Internal Notes:
C
C     Data set-up is done via assignment statements to avoid modifying
C     DATA-loaded arrays, as required by the 1989 SLATEC Guidelines.
C     Run with KPRINT=3 to display the data.
C**End
C
C  Declare arguments.
C
      INTEGER  LUN, KPRINT, IPASS
C
C  DECLARE VARIABLES.
C
      INTEGER  MAXN, MAXN2, MAXN3, NB
      PARAMETER  (MAXN = 16,  MAXN2 = 8,  MAXN3 = 6,  NB = 7)
      INTEGER  I, IERR, IFAIL, INCFD, ISMEX1(MAXN), ISMEX2(MAXN2),
     *         ISMEX3(MAXN3), ISMEXB(NB), ISMON(MAXN), K, N, NS(3)
      REAL  D(MAXN), DB(NB), F(MAXN), FB(NB), X(MAXN)
      LOGICAL  SKIP
C
C  DEFINE EXPECTED RESULTS.
C
      DATA  ISMEX1 / 1, 1,-1, 1, 1,-1, 1, 1,-1, 1, 1,-1, 1, 1,-1, 2/
      DATA  ISMEX2 / 1, 2, 2, 1, 2, 2, 1, 2/
      DATA  ISMEX3 / 1, 1, 1, 1, 1, 1/
      DATA  ISMEXB / 1, 3, 1, -1, -3, -1, 2/
C
C  DEFINE TEST DATA.
C
      DATA  NS /16, 8, 6/
C
C***FIRST EXECUTABLE STATEMENT  PCHQK4
      IF (KPRINT .GE. 3)  WRITE (LUN, 1000)
      IF (KPRINT .GE. 2)  WRITE (LUN, 1001)
C
C       Define X, F, D.
      DO 1  I = 1, MAXN
         X(I) = I
         D(I) = 0.E0
    1 CONTINUE
      DO 2  I = 2, MAXN, 3
         D(I) = 2.E0
    2 CONTINUE
      DO 3  I = 1, 3
         F(I) = X(I)
         F(I+ 3) = F(I  ) + 1.E0
         F(I+ 6) = F(I+3) + 1.E0
         F(I+ 9) = F(I+6) + 1.E0
         F(I+12) = F(I+9) + 1.E0
    3 CONTINUE
      F(16) = 6.E0
C       Define FB, DB.
      FB(1) = 0.E0
      FB(2) = 2.E0
      FB(3) = 3.E0
      FB(4) = 5.E0
      DB(1) = 1.E0
      DB(2) = 3.E0
      DB(3) = 3.E0
      DB(4) = 0.E0
      DO 4  I = 1, 3
         FB(NB-I+1) =  FB(I)
         DB(NB-I+1) = -DB(I)
    4 CONTINUE
C
C  INITIALIZE.
C
      IFAIL = 0
C
      IF (KPRINT .GE. 3)  THEN
         WRITE (LUN, 1002)
         DO 10  I = 1, NB
            WRITE (LUN, 1010)  I, X(I), F(I), D(I), FB(I), DB(I)
   10    CONTINUE
         DO 20  I = NB+1, MAXN
            WRITE (LUN, 1010)  I, X(I), F(I), D(I)
   20    CONTINUE
      ENDIF
C
C  TRANSFER POINT FOR SECOND SET OF TESTS.
C
   25 CONTINUE
C
C  Loop over a series of values of INCFD.
C
      DO 30  INCFD = 1, 3
         N = NS(INCFD)
         SKIP = .FALSE.
C        -------------------------------------------------
         CALL PCHCM (N, X, F, D, INCFD, SKIP, ISMON, IERR)
C        -------------------------------------------------
         IF (KPRINT.GE.3)
     *      WRITE (LUN, 2000)  INCFD, IERR, (ISMON(I), I=1,N)
         IF ( IERR.NE.0 )  THEN
            IFAIL = IFAIL + 1
            IF (KPRINT.GE.3)  WRITE (LUN,2001)
         ELSE
            DO 29  I = 1, N
               IF (INCFD.EQ.1)  THEN
                  IF ( ISMON(I).NE.ISMEX1(I) )  THEN
                     IFAIL = IFAIL + 1
                     IF (KPRINT.GE.3)
     *                  WRITE (LUN, 2002)  (ISMEX1(K),K=1,N)
                     GO TO 30
                  ENDIF
               ELSE IF (INCFD.EQ.2) THEN
                  IF ( ISMON(I).NE.ISMEX2(I) )  THEN
                     IFAIL = IFAIL + 1
                     IF (KPRINT.GE.3)
     *                  WRITE (LUN, 2002)  (ISMEX2(K),K=1,N)
                     GO TO 30
                  ENDIF
               ELSE
                  IF ( ISMON(I).NE.ISMEX3(I) )  THEN
                     IFAIL = IFAIL + 1
                     IF (KPRINT.GE.3)
     *                  WRITE (LUN, 2002)  (ISMEX3(K),K=1,N)
                     GO TO 30
                  ENDIF
               ENDIF
   29       CONTINUE
         ENDIF
   30 CONTINUE
C
C  Test for -1,3,1 bug.
C
      SKIP = .FALSE.
C     ------------------------------------------------
      CALL PCHCM (NB, X, FB, DB, 1, SKIP, ISMON, IERR)
C     ------------------------------------------------
      IF (KPRINT.GE.3)
     *   WRITE (LUN, 2030)  IERR, (ISMON(I), I=1,NB)
      IF ( IERR.NE.0 )  THEN
         IFAIL = IFAIL + 1
         IF (KPRINT.GE.3)  WRITE (LUN,2001)
      ELSE
         DO 34  I = 1, NB
            IF ( ISMON(I).NE.ISMEXB(I) )  THEN
               IFAIL = IFAIL + 1
               IF (KPRINT.GE.3)
     *            WRITE (LUN, 2002)  (ISMEXB(K),K=1,NB)
               GO TO 35
            ENDIF
   34    CONTINUE
      ENDIF
   35 CONTINUE
C
      IF (F(1).LT.0.)  GO TO 90
C
C  Change sign and do again.
C
      IF (KPRINT.GE.3)  WRITE (LUN, 2050)
      DO 40  I = 1, MAXN
         F(I) = -F(I)
         D(I) = -D(I)
         IF ( ISMEX1(I).NE.2 )  ISMEX1(I) = -ISMEX1(I)
   40 CONTINUE
      DO 42  I = 1, MAXN2
         IF ( ISMEX2(I).NE.2 )  ISMEX2(I) = -ISMEX2(I)
   42 CONTINUE
      DO 43  I = 1, MAXN3
         IF ( ISMEX3(I).NE.2 )  ISMEX3(I) = -ISMEX3(I)
   43 CONTINUE
      DO 50  I = 1, NB
         FB(I) = -FB(I)
         DB(I) = -DB(I)
         IF ( ISMEXB(I).NE.2 )  ISMEXB(I) = -ISMEXB(I)
   50 CONTINUE
      GO TO 25
C
C  PRINT SUMMARY AND TERMINATE.
C
   90 CONTINUE
      IF ((KPRINT.GE.2).AND.(IFAIL.NE.0))  WRITE (LUN, 3001)  IFAIL
C
      IF (IFAIL.EQ.0)  THEN
         IPASS = 1
         IF (KPRINT.GE.2) WRITE(LUN,99998)
      ELSE
         IPASS = 0
         IF (KPRINT.GE.1) WRITE(LUN,99999)
      ENDIF
C
      RETURN
C
C  FORMATS.
C
 1000 FORMAT ('1'//10X,'TEST PCHIP MONOTONICITY CHECKER')
 1001 FORMAT (//10X,'PCHQK4 RESULTS'/10X,'--------------')
 1002 FORMAT (// 5X,'DATA:'
     *        // 9X,'I',4X,'X',5X,'F',5X,'D',5X,'FB',4X,'DB')
 1010 FORMAT (5X,I5,5F6.1)
 2000 FORMAT (/4X,'INCFD =',I2,':  IERR =',I3/15X,'ISMON =',16I3)
 2001 FORMAT (' *** Failed -- bad IERR value.')
 2002 FORMAT (' *** Failed -- expect:',16I3)
 2030 FORMAT (/4X,' Bug test:  IERR =',I3/15X,'ISMON =',7I3)
 2050 FORMAT (/4X,'Changing sign of data.....')
 3001 FORMAT (/' *** TROUBLE ***',I5,' MONOTONICITY TESTS FAILED.')
99998 FORMAT (/' ------------  PCHIP PASSED  ALL MONOTONICITY TESTS',
     *        ' ------------')
99999 FORMAT (/' ************  PCHIP FAILED SOME MONOTONICITY TESTS',
     *        ' ************')
C------------- LAST LINE OF PCHQK4 FOLLOWS -----------------------------
      END
