*DECK SLAPQC
      SUBROUTINE SLAPQC (LUN, KPRINT, IPASS)
C***BEGIN PROLOGUE  SLAPQC
C***PURPOSE  Quick check for testing Sparse Linear Algebra Package
C            (SLAP) Version 2.0.2.
C***LIBRARY   SLATEC (SLAP)
C***CATEGORY  D2A4, D2B4
C***TYPE      SINGLE PRECISION (SLAPQC-S, DLAPQC-D)
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
C         SRMGEN) and then runs all the various methods with all the
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
C***ROUTINES CALLED  OUTERR, R1MACH, SCPPLT, SRMGEN, SS2Y,
C                    SSDBCG, SSDCG, SSDCGN, SSDCGS, SSDGMR, SSDOMN,
C                    SSGS, SSICCG, SSILUR, SSJAC, SSLUBC, SSLUCN,
C                    SSLUCS, SSLUGM, SSLUOM, VFILL, XERMAX, XSETF,
C                    XSETUN
C***COMMON BLOCKS    SSLBLK
C***REVISION HISTORY  (YYMMDD)
C   880601  DATE WRITTEN
C   881213  Revised to meet the new SLATEC prologue standards.
C   890920  Modified to reduce single/double differences and to meet
C           SLATEC standards, as requested at July 1989 CML Meeting.
C   891003  Reduced MAXN to a more reasonable size for quick check.
C   920401  Made routine a SUBROUTINE and made necessary changes to
C           interface with a SLATEC quick check driver.  (WRB)
C   920407  COMMON BLOCK renamed SSLBLK.  (WRB)
C   920511  Added complete declaration section.  (WRB)
C   920602  Eliminated unnecessary variables IOUT and ISTDO and made
C           various cosmetic changes.  (FNF)
C   920602  Reduced problem size for a shorter-running test and
C           corrected lower limit in "DO 80" statement.  (FNF)
C   921021  Added 1P's to output formats.  (FNF)
C***END PROLOGUE  SLAPQC
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
      REAL SOLN(MAXN)
C     .. Local Scalars ..
      REAL DENS, ERR, FACTOR, TOL
      INTEGER IERR, ISYM, ITER, ITMAX, ITOL, ITOLGM, IUNIT, K, KASE,
     +        LENIW, LENW, N, NELT, NELTMX, NFAIL, NMAX, NSAVE
C     .. Local Arrays ..
      REAL A(MXNELT), F(MAXN), RWORK(MAXRW), XITER(MAXN)
      INTEGER IA(MXNELT), IWORK(MAXIW), JA(MXNELT)
C     .. External Functions ..
      REAL R1MACH
      EXTERNAL R1MACH
C     .. External Subroutines ..
      EXTERNAL OUTERR, SCPPLT, SRMGEN, SS2Y, SSDBCG, SSDCG, SSDCGN,
     +         SSDCGS, SSDGMR, SSDOMN, SSGS, SSICCG, SSILUR, SSJAC,
     +         SSLUBC, SSLUCN, SSLUCS, SSLUGM, SSLUOM, VFILL
C     .. Intrinsic Functions ..
      INTRINSIC MAX, REAL
C     .. Common blocks ..
      COMMON /SSLBLK/ SOLN
C
C     The following lines are for the braindamaged Sun FPE handler.
C
C$$$      integer oldmode, fpmode
C***FIRST EXECUTABLE STATEMENT  SLAPQC
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
      FACTOR = 1.2E0
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
      TOL = MAX(1.0E3*R1MACH(3),1.0E-6)
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
            CALL SRMGEN( NELTMX, FACTOR, IERR, N, NELT,
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
            CALL SS2Y( N, NELT, IA, JA, A, ISYM )
            IF( KPRINT.GE.4 ) THEN
               WRITE(LUN,1040) (K,IA(K),JA(K),A(K),K=1,NELT)
               CALL SCPPLT( N, NELT, IA, JA, A, ISYM, LUN )
            ENDIF
C
C**********************************************************************
C                    BEGINNING OF SLAP QUICK TESTS
C**********************************************************************
C
C         * * * * * *   SSJAC   * * * * * *
C
            IF( KPRINT.GE.3 ) THEN
              WRITE(LUN,1000) 'SSJAC ', ITOL, ISYM
            ENDIF
            CALL VFILL( N, XITER, 0.0E0 )
C
            CALL SSJAC(N, F, XITER, NELT, IA, JA, A, ISYM,
     $           ITOL, TOL, 2*ITMAX, ITER, ERR, IERR, IUNIT,
     $           RWORK, LENW, IWORK, LENIW )
C
            CALL OUTERR( 'SSJAC ',IERR,KPRINT,NFAIL,LUN,ITER,ERR )
C
C         * * * * *  SSGS  * * * * *
C
            IF( KPRINT.GE.3 ) THEN
              WRITE(LUN,1000) 'SSGS  ',ITOL,ISYM
            ENDIF
            CALL VFILL( N, XITER, 0.0E0 )
C
            CALL SSGS(N, F, XITER, NELT, IA, JA, A, ISYM,
     $           ITOL, TOL, ITMAX, ITER, ERR, IERR, IUNIT,
     $           RWORK, LENW, IWORK, LENIW )
C
            CALL OUTERR( 'SSGS  ',IERR,KPRINT,NFAIL,LUN,ITER,ERR )
C
C         * * * * * *   SSILUR   * * * * * *
C
            IF( KPRINT.GE.3 ) THEN
              WRITE(LUN,1000) 'SSILUR',ITOL,ISYM
            ENDIF
            CALL VFILL( N, XITER, 0.0E0 )
C
            CALL SSILUR(N, F, XITER, NELT, IA, JA, A, ISYM,
     $           ITOL, TOL, ITMAX, ITER, ERR, IERR, IUNIT,
     $           RWORK, LENW, IWORK, LENIW )
C
            CALL OUTERR( 'SSILUR',IERR,KPRINT,NFAIL,LUN,ITER,ERR )
C
C         * * * * * *   SSDCG    * * * * * *
C
            IF( ISYM.EQ.1 ) THEN
               IF( KPRINT.GE.3 ) THEN
                  WRITE(LUN,1000) 'SSDCG',ITOL,ISYM
               ENDIF
               CALL VFILL( N, XITER, 0.0E0 )
C
               CALL SSDCG(N, F, XITER, NELT, IA, JA, A, ISYM,
     $              ITOL, TOL, ITMAX, ITER, ERR, IERR, IUNIT,
     $              RWORK, LENW, IWORK, LENIW )
C
               CALL OUTERR( 'SSDCG ',IERR,KPRINT,NFAIL,LUN,ITER,ERR )
            ENDIF
C
C         * * * * * *    SSICCG    * * * * * *
C
            IF( ISYM.EQ.1 ) THEN
               IF( KPRINT.GE.3 ) THEN
                  WRITE(LUN,1000) 'SSICCG',ITOL,ISYM
               ENDIF
               CALL VFILL( N, XITER, 0.0E0 )
C
               CALL SSICCG(N, F, XITER, NELT, IA, JA, A, ISYM,
     $              ITOL, TOL, ITMAX, ITER, ERR, IERR, IUNIT, RWORK,
     $              LENW, IWORK, LENIW )
C
               CALL OUTERR( 'SSICCG',IERR,KPRINT,NFAIL,LUN,ITER,ERR )
            ENDIF
C
C         * * * * * *    SSDCGN   * * * * * *
C
            IF( KPRINT.GE.3 ) THEN
               WRITE(LUN,1000) 'SSDCGN',ITOL,ISYM
            ENDIF
            CALL VFILL( N, XITER, 0.0E0 )
C
            CALL SSDCGN(N, F, XITER, NELT, IA, JA, A, ISYM, ITOL,
     $           TOL, ITMAX, ITER, ERR, IERR, IUNIT, RWORK, LENW,
     $           IWORK, LENIW )
C
            CALL OUTERR( 'SSDCGN',IERR,KPRINT,NFAIL,LUN,ITER,ERR )
C
C         * * * * * *   SSLUCN   * * * * * *
C
            IF( KPRINT.GE.3 ) THEN
               WRITE(LUN,1000) 'SSLUCN',ITOL,ISYM
            ENDIF
            CALL VFILL( N, XITER, 0.0E0 )
C
            CALL SSLUCN(N, F, XITER, NELT, IA, JA, A, ISYM, ITOL,
     $           TOL, ITMAX, ITER, ERR, IERR, IUNIT, RWORK, LENW,
     $           IWORK, LENIW )
C
            CALL OUTERR( 'SSLUCN',IERR,KPRINT,NFAIL,LUN,ITER,ERR )
C
C         * * * * * *    SSDBCG   * * * * * *
C
            IF( KPRINT.GE.3 ) THEN
               WRITE(LUN,1000) 'SSDBCG',ITOL,ISYM
            ENDIF
            CALL VFILL( N, XITER, 0.0E0 )
C
            CALL SSDBCG(N, F, XITER, NELT, IA, JA, A, ISYM, ITOL,
     $           TOL, ITMAX, ITER, ERR, IERR, IUNIT, RWORK, LENW,
     $           IWORK, LENIW )
C
            CALL OUTERR( 'SSDBCG',IERR,KPRINT,NFAIL,LUN,ITER,ERR )
C
C         * * * * * *   SSLUBC   * * * * * *
C
            IF( KPRINT.GE.3 ) THEN
               WRITE(LUN,1000) 'SSLUBC',ITOL,ISYM
            ENDIF
            CALL VFILL( N, XITER, 0.0E0 )
C
            CALL SSLUBC(N, F, XITER, NELT, IA, JA, A, ISYM,
     $           ITOL, TOL, ITMAX, ITER, ERR, IERR, IUNIT,
     $           RWORK, LENW, IWORK, LENIW )
C
            CALL OUTERR( 'SSLUBC',IERR,KPRINT,NFAIL,LUN,ITER,ERR )
C
C         * * * * * *    SSDCGS   * * * * * *
C
            IF( KPRINT.GE.3 ) THEN
               WRITE(LUN,1000) 'SSDCGS',ITOL,ISYM
            ENDIF
            CALL VFILL( N, XITER, 0.0E0 )
C
            CALL SSDCGS(N, F, XITER, NELT, IA, JA, A, ISYM, ITOL,
     $           TOL, ITMAX, ITER, ERR, IERR, IUNIT, RWORK, LENW,
     $           IWORK, LENIW )
C
            CALL OUTERR( 'SSDCGS',IERR,KPRINT,NFAIL,LUN,ITER,ERR )
C
C         * * * * * *   SSLUCS   * * * * * *
C
            IF( KPRINT.GE.3 ) THEN
               WRITE(LUN,1000) 'SSLUCS',ITOL,ISYM
            ENDIF
            CALL VFILL( N, XITER, 0.0E0 )
C
            CALL SSLUCS(N, F, XITER, NELT, IA, JA, A, ISYM,
     $           ITOL, TOL, ITMAX, ITER, ERR, IERR, IUNIT,
     $           RWORK, LENW, IWORK, LENIW )
C
            CALL OUTERR( 'SSLUCS',IERR,KPRINT,NFAIL,LUN,ITER,ERR )
C
C         * * * * * *    SSDOMN   * * * * * *
C
CVD$ NOVECTOR
            DO 30 NSAVE = 0, 3
               IF( KPRINT.GE.3 ) THEN
                  WRITE(LUN,1010) 'SSDOMN',ITOL, ISYM, NSAVE
               ENDIF
               CALL VFILL( N, XITER, 0.0E0 )
C
               CALL SSDOMN(N, F, XITER, NELT, IA, JA, A,
     $              ISYM, NSAVE, ITOL, TOL, ITMAX, ITER, ERR, IERR,
     $              IUNIT, RWORK, LENW, IWORK, LENIW )
C
               CALL OUTERR( 'SSDOMN',IERR,KPRINT,NFAIL,LUN,ITER,ERR )
 30         CONTINUE
C
C         * * * * * *   SSLUOM   * * * * * *
C
CVD$ NOVECTOR
            DO 40 NSAVE=0,3
               IF( KPRINT.GE.3 ) THEN
                  WRITE(LUN,1010) 'SSLUOM',ITOL, ISYM, NSAVE
               ENDIF
               CALL VFILL( N, XITER, 0.0E0 )
C
               CALL SSLUOM(N, F, XITER, NELT, IA, JA, A,
     $              ISYM, NSAVE, ITOL, TOL, ITMAX, ITER, ERR, IERR,
     $              IUNIT, RWORK, LENW, IWORK, LENIW )
C
               CALL OUTERR( 'SSLUOM',IERR,KPRINT,NFAIL,LUN,ITER,ERR )
 40         CONTINUE
C
C         * * * * * *   SSDGMR   * * * * * *
C
CVD$ NOVECTOR
            DO 50 NSAVE = 5, 12
               IF( KPRINT.GE.3 ) THEN
                  WRITE(LUN,1010) 'SSDGMR',ITOL, ISYM, NSAVE
               ENDIF
               CALL VFILL( N, XITER, 0.0E0 )
               ITOLGM = 0
C
               CALL SSDGMR(N, F, XITER, NELT, IA, JA, A,
     $              ISYM, NSAVE, ITOLGM, TOL, ITMAX, ITER, ERR, IERR,
     $              IUNIT, RWORK, LENW, IWORK, LENIW )
C
               CALL OUTERR( 'SSDGMR',IERR,KPRINT,NFAIL,LUN,ITER,ERR )
 50         CONTINUE
C
C         * * * * * *   SSLUGM   * * * * * *
C
CVD$ NOVECTOR
            DO 60 NSAVE = 5, 12
               IF( KPRINT.GE.3 ) THEN
                  WRITE(LUN,1010) 'SSLUGM',ITOL, ISYM, NSAVE
               ENDIF
               CALL VFILL( N, XITER, 0.0E0 )
C
               CALL SSLUGM(N, F, XITER, NELT, IA, JA, A,
     $              ISYM, NSAVE, ITOL, TOL, ITMAX, ITER, ERR, IERR,
     $              IUNIT, RWORK, LENW, IWORK, LENIW )
C
               CALL OUTERR( 'SSLUGM',IERR,KPRINT,NFAIL,LUN,ITER,ERR )
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
  990 FORMAT(/1X, 'SLAPQC -- Fatal error ', I1, ' generating ',
     $       '*RANDOM* Matrix.')
 1000 FORMAT(/1X,A6,' : ITOL = ',I2,'   ISYM = ',I1)
 1010 FORMAT(/1X,A6,' : ITOL = ',I2,'   ISYM = ',I1,' NSAVE = ',I2)
 1020 FORMAT(/'                * RANDOM Matrix of size',I5,'*'
     $     /'                ',
     $     'Number of non-zeros & Density = ', I5,1P,E16.7)
 1030 FORMAT('                Error tolerance = ',1P,E16.7)
 1040 FORMAT(/'  ***** SLAP Column Matrix *****'/
     $        ' Indx   ia   ja     a'/(1X,I4,1X,I4,1X,I4,1X,1P,E16.7))
 1050 FORMAT('1'/' Running tests with  N =',I3,',  KASE =',I2,
     $                          ',  ISYM =',I2)
 5001 FORMAT('--------- All single precision SLAP tests passed ',
     $       '---------')
 5002 FORMAT('*********',I3,' single precision SLAP tests failed ',
     $       '*********')
      END
