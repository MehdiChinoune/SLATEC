*DECK DPCHQ3
      SUBROUTINE DPCHQ3 (LUN, KPRINT, IPASS)
C***BEGIN PROLOGUE  DPCHQ3
C***PURPOSE  Test the PCHIP interpolators DPCHIC, DPCHIM, DPCHSP.
C***LIBRARY   SLATEC (PCHIP)
C***TYPE      DOUBLE PRECISION (PCHQK3-S, DPCHQ3-D)
C***KEYWORDS  PCHIP INTERPOLATOR QUICK CHECK
C***AUTHOR  Fritsch, F. N., (LLNL)
C***DESCRIPTION
C
C             DPCHIP QUICK CHECK NUMBER 3
C
C     TESTS THE INTERPOLATORS:  DPCHIC, DPCHIM, DPCHSP.
C *Usage:
C
C        INTEGER  LUN, KPRINT, IPASS
C
C        CALL DPCHQ3 (LUN, KPRINT, IPASS)
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
C   This routine interpolates a constructed data set with all three
C   DPCHIP interpolators and compares the results with those obtained
C   on a Cray X/MP.  Two different values of the DPCHIC parameter SWITCH
C   are used.
C
C *Remarks:
C     1. The Cray results are given only to nine significant figures,
C        so don't expect them to match to more.
C     2. The results will depend to some extent on the accuracy of
C        the EXP function.
C
C***ROUTINES CALLED  COMP, D1MACH, DPCHIC, DPCHIM, DPCHSP
C***REVISION HISTORY  (YYMMDD)
C   900309  DATE WRITTEN
C   900314  Converted to a subroutine and added a SLATEC 4.0 prologue.
C   900315  Revised prologue and improved some output formats.  (FNF)
C   900316  Made TOLD machine-dependent and added extra output when
C           KPRINT=3.  (FNF)
C   900320  Added E0's to DATA statement for X to reduce single/double
C           differences, and other minor cosmetic changes.
C   900320  Converted to double precision.
C   900321  Removed IFAIL from call sequence for SLATEC standards and
C           made miscellaneous cosmetic changes.  (FNF)
C   900322  Minor changes to reduce single/double differences.  (FNF)
C   900530  Tolerance (TOLD) and argument to DPCHIC changed.  (WRB)
C   900802  Modified TOLD formula and constants in DPCHIC calls to
C           correct DPCHQ3 failures.  (FNF)
C   901130  Several significant changes:  (FNF)
C           1. Changed comparison between DPCHIM and DPCHIC to only
C              require agreement to machine precision.
C           2. Revised to print more output when KPRINT=3.
C           3. Added 1P's to formats.
C   910708  Minor modifications in use of KPRINT.  (WRB)
C   930317  Improved output formats.  (FNF)
C***END PROLOGUE  DPCHQ3
C
C*Internal Notes:
C
C     TOLD is used to compare with stored Cray results.  Its value
C          should be consistent with significance of stored values.
C     TOLZ is used for cases in which exact equality is expected.
C     TOL  is used for cases in which agreement to machine precision
C          is expected.
C**End
C
C  Declare arguments.
C
      INTEGER  LUN, KPRINT, IPASS
      LOGICAL  COMP
      DOUBLE PRECISION  D1MACH
C
C  Declare variables.
C
      INTEGER  I, IC(2), IERR, IFAIL, N, NBAD, NBADZ, NWK
      PARAMETER  (N = 9,  NWK = 2*N)
      DOUBLE PRECISION  D(N), DC(N), DC5, DC6, DM(N), DS(N), ERR, F(N),
     *                MONE, TOL, TOLD, TOLZ, VC(2), X(N), WK(NWK), ZERO
      PARAMETER  (ZERO = 0.0D0,  MONE = -1.0D0)
      CHARACTER*6  RESULT
C
C  Initialize.
C
C       Data.
      DATA  IC /0, 0/
      DATA  X /-2.2D0,-1.2D0,-1.0D0,-0.5D0,-0.01D0, 0.5D0, 1.0D0,
     *          2.0D0, 2.2D0/
C
C       Results generated on Cray X/MP (9 sign. figs.)
      DATA  DM / 0.            , 3.80027352D-01, 7.17253009D-01,
     *           5.82014161D-01, 0.            ,-5.68208031D-01,
     *          -5.13501618D-01,-7.77910977D-02,-2.45611117D-03/
      DATA  DC5,DC6 / 1.76950158D-02,-5.69579814D-01/
      DATA  DS /-5.16830792D-02, 5.71455855D-01, 7.40530225D-01,
     *           7.63864934D-01, 1.92614386D-02,-7.65324380D-01,
     *          -7.28209035D-01,-7.98445427D-02,-2.85983446D-02/
C
C***FIRST EXECUTABLE STATEMENT  DPCHQ3
      IFAIL = 0
C
C        Set tolerances.
      TOL  = 10*D1MACH(4)
      TOLD = MAX( 1.0D-7, 10*TOL )
      TOLZ = ZERO
C
      IF (KPRINT .GE. 3)  WRITE (LUN, 1000)
      IF (KPRINT .GE. 2)  WRITE (LUN, 1001)
C
C  Set up data.
C
      DO 10  I = 1, N
         F(I) = EXP(-X(I)**2)
   10 CONTINUE
C
      IF (KPRINT .GE. 3)  THEN
         WRITE (LUN, 1002)
         DO 12  I = 1, 4
            WRITE (LUN, 1010)  X(I), F(I), DM(I), DS(I)
   12    CONTINUE
            WRITE (LUN, 1011)  X(5), F(5), DM(5), DC5, DS(5)
            WRITE (LUN, 1011)  X(6), F(6), DM(6), DC6, DS(6)
         DO 15  I = 7, N
            WRITE (LUN, 1010)  X(I), F(I), DM(I), DS(I)
   15    CONTINUE
      ENDIF
C
C  Test DPCHIM.
C
      IF (KPRINT.GE.3)  WRITE (LUN, 2000) 'IM'
C     --------------------------------
      CALL DPCHIM (N, X, F, D, 1, IERR)
C     --------------------------------
C        Expect IERR=1 (one monotonicity switch).
      IF ( KPRINT.GE.3 )  WRITE (LUN, 2001) 1
      IF ( .NOT.COMP (IERR, 1, LUN, KPRINT) )  THEN
         IFAIL = IFAIL + 1
      ELSE
         IF ( KPRINT.GE.3 )  WRITE (LUN, 2002)
         NBAD = 0
         NBADZ = 0
         DO 20  I = 1, N
            RESULT = '  OK'
C             D-values should agree with stored values.
C               (Zero values should agree exactly.)
            IF ( DM(I).EQ.ZERO )  THEN
               ERR = ABS( D(I) )
               IF ( ERR.GT.TOLZ )  THEN
                  NBADZ = NBADZ + 1
                  RESULT = '**BADZ'
               ENDIF
            ELSE
               ERR = ABS( (D(I)-DM(I))/DM(I) )
               IF ( ERR.GT.TOLD )  THEN
                  NBAD = NBAD + 1
                  RESULT = '**BAD'
               ENDIF
            ENDIF
            IF (KPRINT.GE.3)
     *         WRITE (LUN, 2003)  I, X(I), D(I), ERR, RESULT
   20    CONTINUE
         IF ( (NBADZ.NE.0).OR.(NBAD.NE.0) )  THEN
            IFAIL = IFAIL + 1
            IF ((NBADZ.NE.0).AND.(KPRINT.GE.2))
     *         WRITE (LUN, 2004)  NBAD
            IF ((NBAD.NE.0).AND.(KPRINT.GE.2))
     *         WRITE (LUN, 2005)  NBAD, 'IM', TOLD
         ELSE
            IF (KPRINT.GE.2)  WRITE (LUN, 2006)  'IM'
         ENDIF
      ENDIF
C
C  Test DPCHIC -- options set to reproduce DPCHIM.
C
      IF (KPRINT.GE.3)  WRITE (LUN, 2000) 'IC'
C     --------------------------------------------------------
      CALL DPCHIC (IC, VC, ZERO, N, X, F, DC, 1, WK, NWK, IERR)
C     --------------------------------------------------------
C        Expect IERR=0 .
      IF ( KPRINT.GE.3 )  WRITE (LUN, 2001) 0
      IF ( .NOT.COMP (IERR, 0, LUN, KPRINT) )  THEN
         IFAIL = IFAIL + 1
      ELSE
         IF ( KPRINT.GE.3 )  WRITE (LUN, 2002)
         NBAD = 0
         DO 30  I = 1, N
            RESULT = '  OK'
C           D-values should agree exactly with those computed by DPCHIM.
C            (To be generous, will only test to machine precision.)
            ERR = ABS( D(I)-DC(I) )
            IF ( ERR.GT.TOL )  THEN
               NBAD = NBAD + 1
               RESULT = '**BAD'
            ENDIF
            IF (KPRINT.GE.3)
     *         WRITE (LUN, 2003)  I, X(I), DC(I), ERR, RESULT
   30    CONTINUE
         IF ( NBAD.NE.0 )  THEN
            IFAIL = IFAIL + 1
            IF (KPRINT.GE.2)  WRITE (LUN, 2005)  NBAD, 'IC', TOL
         ELSE
            IF (KPRINT.GE.2)  WRITE (LUN, 2006)  'IC'
         ENDIF
      ENDIF
C
C  Test DPCHIC -- default nonzero switch derivatives.
C
      IF (KPRINT.GE.3)  WRITE (LUN, 2000) 'IC'
C     -------------------------------------------------------
      CALL DPCHIC (IC, VC, MONE, N, X, F, D, 1, WK, NWK, IERR)
C     -------------------------------------------------------
C        Expect IERR=0 .
      IF ( KPRINT.GE.3 )  WRITE (LUN, 2001) 0
      IF ( .NOT.COMP (IERR, 0, LUN, KPRINT) )  THEN
         IFAIL = IFAIL + 1
      ELSE
         IF ( KPRINT.GE.3 )  WRITE (LUN, 2002)
         NBAD = 0
         NBADZ = 0
         DO 40  I = 1, N
            RESULT = '  OK'
C            D-values should agree exactly with those computed in
C            previous call, except at points 5 and 6.
            IF ( (I.LT.5).OR.(I.GT.6) )  THEN
               ERR = ABS( D(I)-DC(I) )
               IF ( ERR.GT.TOLZ )  THEN
                  NBADZ = NBADZ + 1
                  RESULT = '**BADA'
               ENDIF
            ELSE
               IF ( I.EQ.5 )  THEN
                  ERR = ABS( (D(I)-DC5)/DC5 )
               ELSE
                  ERR = ABS( (D(I)-DC6)/DC6 )
               ENDIF
               IF ( ERR.GT.TOLD )  THEN
                  NBAD = NBAD + 1
                  RESULT = '**BAD'
               ENDIF
            ENDIF
            IF (KPRINT.GE.3)
     *         WRITE (LUN, 2003)  I, X(I), D(I), ERR, RESULT
   40    CONTINUE
         IF ( (NBADZ.NE.0).OR.(NBAD.NE.0) )  THEN
            IFAIL = IFAIL + 1
            IF ((NBADZ.NE.0).AND.(KPRINT.GE.2))
     *         WRITE (LUN, 2007)  NBAD
            IF ((NBAD.NE.0).AND.(KPRINT.GE.2))
     *         WRITE (LUN, 2005)  NBAD, 'IC', TOLD
         ELSE
            IF (KPRINT.GE.2)  WRITE (LUN, 2006)  'IC'
         ENDIF
      ENDIF
C
C  Test DPCHSP.
C
      IF (KPRINT.GE.3)  WRITE (LUN, 2000) 'SP'
C     -------------------------------------------------
      CALL DPCHSP (IC, VC, N, X, F, D, 1, WK, NWK, IERR)
C     -------------------------------------------------
C        Expect IERR=0 .
      IF ( KPRINT.GE.3 )  WRITE (LUN, 2001) 0
      IF ( .NOT.COMP (IERR, 0, LUN, KPRINT) )  THEN
         IFAIL = IFAIL + 1
      ELSE
         IF ( KPRINT.GE.3 )  WRITE (LUN, 2002)
         NBAD = 0
         DO 50  I = 1, N
            RESULT = '  OK'
C             D-values should agree with stored values.
            ERR = ABS( (D(I)-DS(I))/DS(I) )
            IF ( ERR.GT.TOLD )  THEN
               NBAD = NBAD + 1
               RESULT = '**BAD'
            ENDIF
            IF (KPRINT.GE.3)
     *         WRITE (LUN, 2003)  I, X(I), D(I), ERR, RESULT
   50    CONTINUE
         IF ( NBAD.NE.0 )  THEN
            IFAIL = IFAIL + 1
            IF (KPRINT.GE.2)  WRITE (LUN, 2005)  NBAD, 'SP', TOLD
         ELSE
            IF (KPRINT.GE.2)  WRITE (LUN, 2006)  'SP'
         ENDIF
      ENDIF
C
C  PRINT SUMMARY AND TERMINATE.
C
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
 1000 FORMAT ('1'//10X,'TEST DPCHIP INTERPOLATORS')
 1001 FORMAT (//10X,'DPCHQ3 RESULTS'/10X,'--------------')
 1002 FORMAT (// 5X,'DATA:'
     *         /39X,'---------- EXPECTED D-VALUES ----------'
     *         /12X,'X',9X,'F',18X,'DM',13X,'DC',13X,'DS')
 1010 FORMAT (5X,F10.2,1P,D15.5,4X,D15.5,15X,D15.5)
 1011 FORMAT (5X,F10.2,1P,D15.5,4X,3D15.5)
 2000 FORMAT (/5X,'DPCH',A2,' TEST:')
 2001 FORMAT (15X,'EXPECT  IERR =',I5)
 2002 FORMAT (/9X,'I',7X,'X',9X,'D',13X,'ERR')
 2003 FORMAT (5X,I5,F10.2,1P,2D15.5,2X,A)
 2004 FORMAT (/'    **',I5,' DPCHIM RESULTS FAILED TO BE EXACTLY ZERO.')
 2005 FORMAT (/'    **',I5,' DPCH',A2,' RESULTS FAILED TOLERANCE TEST.',
     *                     '  TOL =',1P,D10.3)
 2006 FORMAT (/5X,'  ALL DPCH',A2,' RESULTS OK.')
 2007 FORMAT (/'    **',I5,' DPCHIC RESULTS FAILED TO AGREE WITH',
     *                      ' PREVIOUS CALL.')
 3001 FORMAT (/' *** TROUBLE ***',I5,' INTERPOLATION TESTS FAILED.')
99998 FORMAT (/' ------------ DPCHIP PASSED  ALL INTERPOLATION TESTS',
     *        ' ------------')
99999 FORMAT (/' ************ DPCHIP FAILED SOME INTERPOLATION TESTS',
     *        ' ************')
C------------- LAST LINE OF DPCHQ3 FOLLOWS -----------------------------
      END
