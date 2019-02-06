*DECK PCHQK3
      SUBROUTINE PCHQK3 (LUN, KPRINT, IPASS)
C***BEGIN PROLOGUE  PCHQK3
C***PURPOSE  Test the PCHIP interpolators PCHIC, PCHIM, PCHSP.
C***LIBRARY   SLATEC (PCHIP)
C***TYPE      SINGLE PRECISION (PCHQK3-S, DPCHQ3-D)
C***KEYWORDS  PCHIP INTERPOLATOR QUICK CHECK
C***AUTHOR  Fritsch, F. N., (LLNL)
C***DESCRIPTION
C
C              PCHIP QUICK CHECK NUMBER 3
C
C     TESTS THE INTERPOLATORS:  PCHIC, PCHIM, PCHSP.
C *Usage:
C
C        INTEGER  LUN, KPRINT, IPASS
C
C        CALL PCHQK3 (LUN, KPRINT, IPASS)
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
C    PCHIP interpolators and compares the results with those obtained
C   on a Cray X/MP.  Two different values of the  PCHIC parameter SWITCH
C   are used.
C
C *Remarks:
C     1. The Cray results are given only to nine significant figures,
C        so don't expect them to match to more.
C     2. The results will depend to some extent on the accuracy of
C        the EXP function.
C
C***ROUTINES CALLED  COMP, PCHIC, PCHIM, PCHSP, R1MACH
C***REVISION HISTORY  (YYMMDD)
C   900309  DATE WRITTEN
C   900314  Converted to a subroutine and added a SLATEC 4.0 prologue.
C   900315  Revised prologue and improved some output formats.  (FNF)
C   900316  Made TOLD machine-dependent and added extra output when
C           KPRINT=3.  (FNF)
C   900320  Added E0's to DATA statement for X to reduce single/double
C           differences, and other minor cosmetic changes.
C   900321  Removed IFAIL from call sequence for SLATEC standards and
C           made miscellaneous cosmetic changes.  (FNF)
C   900322  Minor changes to reduce single/double differences.  (FNF)
C   900530  Tolerance (TOLD) changed.  (WRB)
C   900802  Modified TOLD formula and constants in PCHIC calls to be
C           compatible with DPCHQ3.  (FNF)
C   901130  Several significant changes:  (FNF)
C           1. Changed comparison between PCHIM and PCHIC to only
C              require agreement to machine precision.
C           2. Revised to print more output when KPRINT=3.
C           3. Added 1P's to formats.
C   910708  Minor modifications in use of KPRINT.  (WRB)
C   930317  Improved output formats.  (FNF)
C***END PROLOGUE  PCHQK3
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
      REAL  R1MACH
C
C  Declare variables.
C
      INTEGER  I, IC(2), IERR, IFAIL, N, NBAD, NBADZ, NWK
      PARAMETER  (N = 9,  NWK = 2*N)
      REAL  D(N), DC(N), DC5, DC6, DM(N), DS(N), ERR, F(N), MONE, TOL,
     *      TOLD, TOLZ, VC(2), X(N), WK(NWK), ZERO
      PARAMETER  (ZERO = 0.0E0,  MONE = -1.0E0)
      CHARACTER*6  RESULT
C
C  Initialize.
C
C       Data.
      DATA  IC /0, 0/
      DATA  X /-2.2E0,-1.2E0,-1.0E0,-0.5E0,-0.01E0, 0.5E0, 1.0E0,
     *          2.0E0, 2.2E0/
C
C       Results generated on Cray X/MP (9 sign. figs.)
      DATA  DM / 0.            , 3.80027352E-01, 7.17253009E-01,
     *           5.82014161E-01, 0.            ,-5.68208031E-01,
     *          -5.13501618E-01,-7.77910977E-02,-2.45611117E-03/
      DATA  DC5,DC6 / 1.76950158E-02,-5.69579814E-01/
      DATA  DS /-5.16830792E-02, 5.71455855E-01, 7.40530225E-01,
     *           7.63864934E-01, 1.92614386E-02,-7.65324380E-01,
     *          -7.28209035E-01,-7.98445427E-02,-2.85983446E-02/
C
C***FIRST EXECUTABLE STATEMENT  PCHQK3
      IFAIL = 0
C
C        Set tolerances.
      TOL  = 10*R1MACH(4)
      TOLD = MAX( 1.0E-7, 10*TOL )
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
C  Test PCHIM.
C
      IF (KPRINT.GE.3)  WRITE (LUN, 2000) 'IM'
C     --------------------------------
      CALL PCHIM (N, X, F, D, 1, IERR)
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
C  Test PCHIC -- options set to reproduce PCHIM.
C
      IF (KPRINT.GE.3)  WRITE (LUN, 2000) 'IC'
C     --------------------------------------------------------
      CALL PCHIC (IC, VC, ZERO, N, X, F, DC, 1, WK, NWK, IERR)
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
C           D-values should agree exactly with those computed by PCHIM.
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
C  Test PCHIC -- default nonzero switch derivatives.
C
      IF (KPRINT.GE.3)  WRITE (LUN, 2000) 'IC'
C     -------------------------------------------------------
      CALL PCHIC (IC, VC, MONE, N, X, F, D, 1, WK, NWK, IERR)
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
C  Test PCHSP.
C
      IF (KPRINT.GE.3)  WRITE (LUN, 2000) 'SP'
C     -------------------------------------------------
      CALL PCHSP (IC, VC, N, X, F, D, 1, WK, NWK, IERR)
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
 1000 FORMAT ('1'//10X,'TEST PCHIP INTERPOLATORS')
 1001 FORMAT (//10X,'PCHQK3 RESULTS'/10X,'--------------')
 1002 FORMAT (// 5X,'DATA:'
     *         /39X,'---------- EXPECTED D-VALUES ----------'
     *         /12X,'X',9X,'F',18X,'DM',13X,'DC',13X,'DS')
 1010 FORMAT (5X,F10.2,1P,E15.5,4X,E15.5,15X,E15.5)
 1011 FORMAT (5X,F10.2,1P,E15.5,4X,3E15.5)
 2000 FORMAT (/5X,'PCH',A2,' TEST:')
 2001 FORMAT (15X,'EXPECT  IERR =',I5)
 2002 FORMAT (/9X,'I',7X,'X',9X,'D',13X,'ERR')
 2003 FORMAT (5X,I5,F10.2,1P,2E15.5,2X,A)
 2004 FORMAT (/'    **',I5,'  PCHIM RESULTS FAILED TO BE EXACTLY ZERO.')
 2005 FORMAT (/'    **',I5,'  PCH',A2,' RESULTS FAILED TOLERANCE TEST.',
     *                     '  TOL =',1P,E10.3)
 2006 FORMAT (/5X,'  ALL  PCH',A2,' RESULTS OK.')
 2007 FORMAT (/'    **',I5,'  PCHIC RESULTS FAILED TO AGREE WITH',
     *                      ' PREVIOUS CALL.')
 3001 FORMAT (/' *** TROUBLE ***',I5,' INTERPOLATION TESTS FAILED.')
99998 FORMAT (/' ------------  PCHIP PASSED  ALL INTERPOLATION TESTS',
     *        ' ------------')
99999 FORMAT (/' ************  PCHIP FAILED SOME INTERPOLATION TESTS',
     *        ' ************')
C------------- LAST LINE OF PCHQK3 FOLLOWS -----------------------------
      END
