*DECK DLAPQC
      SUBROUTINE DLAPQC (LUN, KPRINT, IPASS)
C***BEGIN PROLOGUE  DLAPQC
C***PURPOSE  Quick check for testing Sparse Linear Algebra Package
C            (SLAP) Version 2.0.2.
C***LIBRARY   SLATEC (SLAP)
C***CATEGORY  D2A4, D2B4
C***TYPE      DOUBLE PRECISION (SLAPQC-S, DLAPQC-D)
C***KEYWORDS  QUICK CHECK, SLAP
C***AUTHOR  Mark K. Seager (LLNL)
C             seager@llnl.gov
C             Lawrence Livermore National Laboratory
C             PO BOX 808, L-300
C             Livermore, CA 94550
C             (510) 423-3141
C***DESCRIPTION
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
C              4  Quick checks - Print complete quick check results.
C                                Prints matrices, etc.  Very verbose!!
C                                                       --------------
C                 Driver       - Pass or fail message printed.
C
C *Description:
C         This is a SLATEC Quick Check program to test the *SLAP*
C         Version 2.0.2 package.  It generates a "random" matrix (See
C         DRMGEN) and then runs all the various methods with all the
C         various preconditioners and all the various stop tests.
C
C         It is assumed that the test is being run interactively and
C         that STDIN (STANDARD INPUT) is Fortran I/O unit I1MACH(1)
C         and STDOUT (STANDARD OUTPUT) is unit I1MACH(2).
C
C         *************************************************************
C         **** WARNING !!! WARNING !!! WARNING !!! WARNING !!! WARNING
C         *************************************************************
C         **** THIS PROGRAM WILL NOT FUNCTION PROPERLY IF THE FORTRAN
C         **** I/O UNITS I1MACH(1) and I1MACH(2) are not connected
C         **** to the program for I/O.
C         *************************************************************
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  D1MACH, DCPPLT, DFILL, DRMGEN, DS2Y, DSDBCG, DSDCG,
C                    DSDCGN, DSDCGS, DSDGMR, DSDOMN, DSGS, DSICCG,
C                    DSILUR, DSJAC, DSLUBC, DSLUCN, DSLUCS, DSLUGM,
C                    DSLUOM, DUTERR, XERMAX, XSETF, XSETUN
C***COMMON BLOCKS    DSLBLK
C***REVISION HISTORY  (YYMMDD)
C   880601  DATE WRITTEN
C   881213  Revised to meet the new SLATEC prologue standards.
C   890920  Modified to reduce single/double differences and to meet
C           SLATEC standards, as requested at July 1989 CML Meeting.
C   891003  Reduced MAXN to a more reasonable size for quick check.
C   920401  Made routine a SUBROUTINE and made necessary changes to
C           interface with a SLATEC quick check driver.  (WRB)
C   920407  COMMON BLOCK renamed DSLBLK.  (WRB)
C   920511  Added complete declaration section.  (WRB)
C   920602  Eliminated unnecessary variables IOUT and ISTDO and made
C           various cosmetic changes.  (FNF)
C   920602  Reduced problem size for a shorter-running test and
C           corrected lower limit in "DO 80" statement.  (FNF)
C   921021  Changed E's to 1P,D's in output formats.  (FNF)
C***END PROLOGUE  DLAPQC
C
C     The problem size, MAXN, should be large enough that the
C     iterative methods do 10-15 iterations, just to be sure that
C     the truncated methods run to the end of their ropes and enter
C     their error recovery mode.  Thus, for a more thorough test change
C     the following PARAMETER statement to:
C     PARAMETER (MAXN=69, MXNELT=5000, MAXIW=5000, MAXRW=5000)
C
C     .. Parameters ..
      INTEGER MAXN, MXNELT, MAXIW, MAXRW
      PARAMETER (MAXN=25, MXNELT=500, MAXIW=1000, MAXRW=1000)
C     .. Scalar Arguments ..
      INTEGER IPASS, KPRINT, LUN
C     .. Arrays in Common ..
      DOUBLE PRECISION SOLN(MAXN)
C     .. Local Scalars ..
      DOUBLE PRECISION DENS, ERR, FACTOR, TOL
      INTEGER IERR, ISYM, ITER, ITMAX, ITOL, ITOLGM, IUNIT, K, KASE,
     +        LENIW, LENW, N, NELT, NELTMX, NFAIL, NMAX, NSAVE
C     .. Local Arrays ..
      DOUBLE PRECISION A(MXNELT), F(MAXN), RWORK(MAXRW), XITER(MAXN)
      INTEGER IA(MXNELT), IWORK(MAXIW), JA(MXNELT)
C     .. External Functions ..
      DOUBLE PRECISION D1MACH
      EXTERNAL D1MACH
C     .. External Subroutines ..
      EXTERNAL DCPPLT, DFILL, DRMGEN, DS2Y, DSDBCG, DSDCG, DSDCGN,
     +         DSDCGS, DSDGMR, DSDOMN, DSGS, DSICCG, DSILUR, DSJAC,
     +         DSLUBC, DSLUCN, DSLUCS, DSLUGM, DSLUOM, DUTERR
C     .. Intrinsic Functions ..
      INTRINSIC MAX, REAL
C     .. Common blocks ..
      COMMON /DSLBLK/ SOLN
C
C     The following lines are for the braindamaged Sun FPE handler.
C
C$$$      integer oldmode, fpmode
C***FIRST EXECUTABLE STATEMENT  DLAPQC
C$$$      oldmode = fpmode( 62464 )
C
C     Maximum problem sizes.
C
      NELTMX = MXNELT
      NMAX   = MAXN
      LENIW  = MAXIW
      LENW   = MAXRW
C
C     Set some input data.
C
      N      = NMAX
      ITMAX  = N
      FACTOR = 1.2D0
C
C     Set to print intermediate results if KPRINT.GE.3.
C
      IF( KPRINT.LT.3 ) THEN
         IUNIT = 0
      ELSE
         IUNIT = LUN
      ENDIF
C
C     Set the Error tolerance to depend on the machine epsilon.
C
      TOL = MAX(1.0D3*D1MACH(3),1.0D-6)
      NFAIL = 0
C
C     Test routines using various convergence criteria.
C
      DO 80 KASE = 1, 3
         IF(KASE .EQ. 1 .OR. KASE .EQ. 2) ITOL = KASE
         IF(KASE .EQ. 3) ITOL = 11
C
C         Test routines using nonsymmetric (ISYM=0) and symmetric
C         storage (ISYM=1).  For ISYM=0 a really non-symmetric matrix
C         is generated.  The amount of non-symmetry is controlled by
C         user.
C
         DO 70 ISYM = 0, 1
            IF( KPRINT.GE.2 )  WRITE (LUN, 1050)  N, KASE, ISYM
C
C         Set up a random matrix.
C
            CALL DRMGEN( NELTMX, FACTOR, IERR, N, NELT,
     $           ISYM, IA, JA, A, F, SOLN, RWORK, IWORK, IWORK(N+1) )
            IF( IERR.NE.0 ) THEN
               WRITE(LUN,990) IERR
               NFAIL = NFAIL + 1
               GO TO 70
            ENDIF
            IF( ISYM.EQ.0 ) THEN
               DENS = REAL(NELT)/(N*N)
            ELSE
               DENS = REAL(2*NELT)/(N*N)
            ENDIF
            IF( KPRINT.GE.2 ) THEN
              WRITE(LUN,1020) N, NELT, DENS
              WRITE(LUN,1030) TOL
            ENDIF
C
C         Convert to the SLAP-Column format and
C         write out matrix in SLAP-Column format, if desired.
C
            CALL DS2Y( N, NELT, IA, JA, A, ISYM )
            IF( KPRINT.GE.4 ) THEN
               WRITE(LUN,1040) (K,IA(K),JA(K),A(K),K=1,NELT)
               CALL DCPPLT( N, NELT, IA, JA, A, ISYM, LUN )
            ENDIF
C
C**********************************************************************
C                    BEGINNING OF SLAP QUICK TESTS
C**********************************************************************
C
C         * * * * * *   DSJAC   * * * * * *
C
            IF( KPRINT.GE.3 ) THEN
              WRITE(LUN,1000) 'DSJAC ', ITOL, ISYM
            ENDIF
            CALL DFILL( N, XITER, 0.0D0 )
C
            CALL DSJAC(N, F, XITER, NELT, IA, JA, A, ISYM,
     $           ITOL, TOL, 2*ITMAX, ITER, ERR, IERR, IUNIT,
     $           RWORK, LENW, IWORK, LENIW )
C
            CALL DUTERR( 'DSJAC ',IERR,KPRINT,NFAIL,LUN,ITER,ERR )
C
C         * * * * *  DSGS  * * * * *
C
            IF( KPRINT.GE.3 ) THEN
              WRITE(LUN,1000) 'DSGS  ',ITOL,ISYM
            ENDIF
            CALL DFILL( N, XITER, 0.0D0 )
C
            CALL DSGS(N, F, XITER, NELT, IA, JA, A, ISYM,
     $           ITOL, TOL, ITMAX, ITER, ERR, IERR, IUNIT,
     $           RWORK, LENW, IWORK, LENIW )
C
            CALL DUTERR( 'DSGS  ',IERR,KPRINT,NFAIL,LUN,ITER,ERR )
C
C         * * * * * *   DSILUR   * * * * * *
C
            IF( KPRINT.GE.3 ) THEN
              WRITE(LUN,1000) 'DSILUR',ITOL,ISYM
            ENDIF
            CALL DFILL( N, XITER, 0.0D0 )
C
            CALL DSILUR(N, F, XITER, NELT, IA, JA, A, ISYM,
     $           ITOL, TOL, ITMAX, ITER, ERR, IERR, IUNIT,
     $           RWORK, LENW, IWORK, LENIW )
C
            CALL DUTERR( 'DSILUR',IERR,KPRINT,NFAIL,LUN,ITER,ERR )
C
C         * * * * * *   DSDCG    * * * * * *
C
            IF( ISYM.EQ.1 ) THEN
               IF( KPRINT.GE.3 ) THEN
                  WRITE(LUN,1000) 'DSDCG',ITOL,ISYM
               ENDIF
               CALL DFILL( N, XITER, 0.0D0 )
C
               CALL DSDCG(N, F, XITER, NELT, IA, JA, A, ISYM,
     $              ITOL, TOL, ITMAX, ITER, ERR, IERR, IUNIT,
     $              RWORK, LENW, IWORK, LENIW )
C
               CALL DUTERR( 'DSDCG ',IERR,KPRINT,NFAIL,LUN,ITER,ERR )
            ENDIF
C
C         * * * * * *    DSICCG    * * * * * *
C
            IF( ISYM.EQ.1 ) THEN
               IF( KPRINT.GE.3 ) THEN
                  WRITE(LUN,1000) 'DSICCG',ITOL,ISYM
               ENDIF
               CALL DFILL( N, XITER, 0.0D0 )
C
               CALL DSICCG(N, F, XITER, NELT, IA, JA, A, ISYM,
     $              ITOL, TOL, ITMAX, ITER, ERR, IERR, IUNIT, RWORK,
     $              LENW, IWORK, LENIW )
C
               CALL DUTERR( 'DSICCG',IERR,KPRINT,NFAIL,LUN,ITER,ERR )
            ENDIF
C
C         * * * * * *    DSDCGN   * * * * * *
C
            IF( KPRINT.GE.3 ) THEN
               WRITE(LUN,1000) 'DSDCGN',ITOL,ISYM
            ENDIF
            CALL DFILL( N, XITER, 0.0D0 )
C
            CALL DSDCGN(N, F, XITER, NELT, IA, JA, A, ISYM, ITOL,
     $           TOL, ITMAX, ITER, ERR, IERR, IUNIT, RWORK, LENW,
     $           IWORK, LENIW )
C
            CALL DUTERR( 'DSDCGN',IERR,KPRINT,NFAIL,LUN,ITER,ERR )
C
C         * * * * * *   DSLUCN   * * * * * *
C
            IF( KPRINT.GE.3 ) THEN
               WRITE(LUN,1000) 'DSLUCN',ITOL,ISYM
            ENDIF
            CALL DFILL( N, XITER, 0.0D0 )
C
            CALL DSLUCN(N, F, XITER, NELT, IA, JA, A, ISYM, ITOL,
     $           TOL, ITMAX, ITER, ERR, IERR, IUNIT, RWORK, LENW,
     $           IWORK, LENIW )
C
            CALL DUTERR( 'DSLUCN',IERR,KPRINT,NFAIL,LUN,ITER,ERR )
C
C         * * * * * *    DSDBCG   * * * * * *
C
            IF( KPRINT.GE.3 ) THEN
               WRITE(LUN,1000) 'DSDBCG',ITOL,ISYM
            ENDIF
            CALL DFILL( N, XITER, 0.0D0 )
C
            CALL DSDBCG(N, F, XITER, NELT, IA, JA, A, ISYM, ITOL,
     $           TOL, ITMAX, ITER, ERR, IERR, IUNIT, RWORK, LENW,
     $           IWORK, LENIW )
C
            CALL DUTERR( 'DSDBCG',IERR,KPRINT,NFAIL,LUN,ITER,ERR )
C
C         * * * * * *   DSLUBC   * * * * * *
C
            IF( KPRINT.GE.3 ) THEN
               WRITE(LUN,1000) 'DSLUBC',ITOL,ISYM
            ENDIF
            CALL DFILL( N, XITER, 0.0D0 )
C
            CALL DSLUBC(N, F, XITER, NELT, IA, JA, A, ISYM,
     $           ITOL, TOL, ITMAX, ITER, ERR, IERR, IUNIT,
     $           RWORK, LENW, IWORK, LENIW )
C
            CALL DUTERR( 'DSLUBC',IERR,KPRINT,NFAIL,LUN,ITER,ERR )
C
C         * * * * * *    DSDCGS   * * * * * *
C
            IF( KPRINT.GE.3 ) THEN
               WRITE(LUN,1000) 'DSDCGS',ITOL,ISYM
            ENDIF
            CALL DFILL( N, XITER, 0.0D0 )
C
            CALL DSDCGS(N, F, XITER, NELT, IA, JA, A, ISYM, ITOL,
     $           TOL, ITMAX, ITER, ERR, IERR, IUNIT, RWORK, LENW,
     $           IWORK, LENIW )
C
            CALL DUTERR( 'DSDCGS',IERR,KPRINT,NFAIL,LUN,ITER,ERR )
C
C         * * * * * *   DSLUCS   * * * * * *
C
            IF( KPRINT.GE.3 ) THEN
               WRITE(LUN,1000) 'DSLUCS',ITOL,ISYM
            ENDIF
            CALL DFILL( N, XITER, 0.0D0 )
C
            CALL DSLUCS(N, F, XITER, NELT, IA, JA, A, ISYM,
     $           ITOL, TOL, ITMAX, ITER, ERR, IERR, IUNIT,
     $           RWORK, LENW, IWORK, LENIW )
C
            CALL DUTERR( 'DSLUCS',IERR,KPRINT,NFAIL,LUN,ITER,ERR )
C
C         * * * * * *    DSDOMN   * * * * * *
C
CVD$ NOVECTOR
            DO 30 NSAVE = 0, 3
               IF( KPRINT.GE.3 ) THEN
                  WRITE(LUN,1010) 'DSDOMN',ITOL, ISYM, NSAVE
               ENDIF
               CALL DFILL( N, XITER, 0.0D0 )
C
               CALL DSDOMN(N, F, XITER, NELT, IA, JA, A,
     $              ISYM, NSAVE, ITOL, TOL, ITMAX, ITER, ERR, IERR,
     $              IUNIT, RWORK, LENW, IWORK, LENIW )
C
               CALL DUTERR( 'DSDOMN',IERR,KPRINT,NFAIL,LUN,ITER,ERR )
 30         CONTINUE
C
C         * * * * * *   DSLUOM   * * * * * *
C
CVD$ NOVECTOR
            DO 40 NSAVE=0,3
               IF( KPRINT.GE.3 ) THEN
                  WRITE(LUN,1010) 'DSLUOM',ITOL, ISYM, NSAVE
               ENDIF
               CALL DFILL( N, XITER, 0.0D0 )
C
               CALL DSLUOM(N, F, XITER, NELT, IA, JA, A,
     $              ISYM, NSAVE, ITOL, TOL, ITMAX, ITER, ERR, IERR,
     $              IUNIT, RWORK, LENW, IWORK, LENIW )
C
               CALL DUTERR( 'DSLUOM',IERR,KPRINT,NFAIL,LUN,ITER,ERR )
 40         CONTINUE
C
C         * * * * * *   DSDGMR   * * * * * *
C
CVD$ NOVECTOR
            DO 50 NSAVE = 5, 12
               IF( KPRINT.GE.3 ) THEN
                  WRITE(LUN,1010) 'DSDGMR',ITOL, ISYM, NSAVE
               ENDIF
               CALL DFILL( N, XITER, 0.0D0 )
               ITOLGM = 0
C
               CALL DSDGMR(N, F, XITER, NELT, IA, JA, A,
     $              ISYM, NSAVE, ITOLGM, TOL, ITMAX, ITER, ERR, IERR,
     $              IUNIT, RWORK, LENW, IWORK, LENIW )
C
               CALL DUTERR( 'DSDGMR',IERR,KPRINT,NFAIL,LUN,ITER,ERR )
 50         CONTINUE
C
C         * * * * * *   DSLUGM   * * * * * *
C
CVD$ NOVECTOR
            DO 60 NSAVE = 5, 12
               IF( KPRINT.GE.3 ) THEN
                  WRITE(LUN,1010) 'DSLUGM',ITOL, ISYM, NSAVE
               ENDIF
               CALL DFILL( N, XITER, 0.0D0 )
C
               CALL DSLUGM(N, F, XITER, NELT, IA, JA, A,
     $              ISYM, NSAVE, ITOL, TOL, ITMAX, ITER, ERR, IERR,
     $              IUNIT, RWORK, LENW, IWORK, LENIW )
C
               CALL DUTERR( 'DSLUGM',IERR,KPRINT,NFAIL,LUN,ITER,ERR )
 60         CONTINUE
 70      CONTINUE
 80   CONTINUE
C
      IF (NFAIL .EQ. 0) THEN
         IPASS = 1
         IF( KPRINT .GE. 2 )  WRITE (LUN, 5001)
      ELSE
         IPASS = 0
         IF( KPRINT .GE. 2 )  WRITE (LUN, 5002) NFAIL
      ENDIF
C
      RETURN
C
  990 FORMAT(/1X, 'DLAPQC -- Fatal error ', I1, ' generating ',
     $       '*RANDOM* Matrix.')
 1000 FORMAT(/1X,A6,' : ITOL = ',I2,'   ISYM = ',I1)
 1010 FORMAT(/1X,A6,' : ITOL = ',I2,'   ISYM = ',I1,' NSAVE = ',I2)
 1020 FORMAT(/'                * RANDOM Matrix of size',I5,'*'
     $     /'                ',
     $     'Number of non-zeros & Density = ', I5,1P,D16.7)
 1030 FORMAT('                Error tolerance = ',1P,D16.7)
 1040 FORMAT(/'  ***** SLAP Column Matrix *****'/
     $        ' Indx   ia   ja     a'/(1X,I4,1X,I4,1X,I4,1X,1P,D16.7))
 1050 FORMAT('1'/' Running tests with  N =',I3,',  KASE =',I2,
     $                          ',  ISYM =',I2)
 5001 FORMAT('--------- All double precision SLAP tests passed ',
     $       '---------')
 5002 FORMAT('*********',I3,' double precision SLAP tests failed ',
     $       '*********')
      END
