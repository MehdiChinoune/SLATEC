MODULE TEST25_MOD
  IMPLICIT NONE

CONTAINS
  !DECK SLAPQC
  SUBROUTINE SLAPQC(Lun,Kprint,Ipass)
    IMPLICIT NONE
    !***BEGIN PROLOGUE  SLAPQC
    !***PURPOSE  Quick check for testing Sparse Linear Algebra Package
    !            (SLAP) Version 2.0.2.
    !***LIBRARY   SLATEC (SLAP)
    !***CATEGORY  D2A4, D2B4
    !***TYPE      SINGLE PRECISION (SLAPQC-S, DLAPQC-D)
    !***KEYWORDS  QUICK CHECK, SLAP
    !***AUTHOR  Mark K. Seager (LLNL)
    !             seager@llnl.gov
    !             Lawrence Livermore National Laboratory
    !             PO BOX 808, L-300
    !             Livermore, CA 94550
    !             (510) 423-3141
    !***DESCRIPTION
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
    !              4  Quick checks - Print complete quick check results.
    !                                Prints matrices, etc.  Very verbose!!
    !                                                       --------------
    !                 Driver       - Pass or fail message printed.
    !
    ! *Description:
    !         This is a SLATEC Quick Check program to test the *SLAP*
    !         Version 2.0.2 package.  It generates a "random" matrix (See
    !         SRMGEN) and then runs all the various methods with all the
    !         various preconditioners and all the various stop tests.
    !
    !         It is assumed that the test is being run interactively and
    !         that STDIN (STANDARD INPUT) is Fortran I/O unit I1MACH(1)
    !         and STDOUT (STANDARD OUTPUT) is unit I1MACH(2).
    !
    !         *************************************************************
    !         **** WARNING !!! WARNING !!! WARNING !!! WARNING !!! WARNING
    !         *************************************************************
    !         **** THIS PROGRAM WILL NOT FUNCTION PROPERLY IF THE FORTRAN
    !         **** I/O UNITS I1MACH(1) and I1MACH(2) are not connected
    !         **** to the program for I/O.
    !         *************************************************************
    !
    !***REFERENCES  (NONE)
    !***ROUTINES CALLED  OUTERR, R1MACH, SCPPLT, SRMGEN, SS2Y,
    !                    SSDBCG, SSDCG, SSDCGN, SSDCGS, SSDGMR, SSDOMN,
    !                    SSGS, SSICCG, SSILUR, SSJAC, SSLUBC, SSLUCN,
    !                    SSLUCS, SSLUGM, SSLUOM, VFILL, XERMAX, XSETF,
    !                    XSETUN
    !***COMMON BLOCKS    SSLBLK
    !***REVISION HISTORY  (YYMMDD)
    !   880601  DATE WRITTEN
    !   881213  Revised to meet the new SLATEC prologue standards.
    !   890920  Modified to reduce single/double differences and to meet
    !           SLATEC standards, as requested at July 1989 CML Meeting.
    !   891003  Reduced MAXN to a more reasonable size for quick check.
    !   920401  Made routine a SUBROUTINE and made necessary changes to
    !           interface with a SLATEC quick check driver.  (WRB)
    !   920407  COMMON BLOCK renamed SSLBLK.  (WRB)
    !   920511  Added complete declaration section.  (WRB)
    !   920602  Eliminated unnecessary variables IOUT and ISTDO and made
    !           various cosmetic changes.  (FNF)
    !   920602  Reduced problem size for a shorter-running test and
    !           corrected lower limit in "DO 80" statement.  (FNF)
    !   921021  Added 1P's to output formats.  (FNF)
    !***END PROLOGUE  SLAPQC
    !
    !     The problem size, MAXN, should be large enough that the
    !     iterative methods do 10-15 iterations, just to be sure that
    !     the truncated methods run to the end of their ropes and enter
    !     their error recovery mode.  Thus, for a more thorough test change
    !     the following PARAMETER statement to:
    !     PARAMETER (MAXN=69, MXNELT=5000, MAXIW=5000, MAXRW=5000)
    !
    !     .. Parameters ..
    INTEGER MAXN, MXNELT, MAXIW, MAXRW
    PARAMETER (MAXN=25,MXNELT=500,MAXIW=1000,MAXRW=1000)
    !     .. Scalar Arguments ..
    INTEGER Ipass, Kprint, Lun
    !     .. Arrays in Common ..
    REAL SOLn(MAXN)
    !     .. Local Scalars ..
    REAL dens, err, factor, tol
    INTEGER ierr, isym, iter, itmax, itol, itolgm, iunit, k, kase, &
      leniw, lenw, n, nelt, neltmx, nfail, nmax, nsave
    !     .. Local Arrays ..
    REAL a(MXNELT), f(MAXN), rwork(MAXRW), xiter(MAXN)
    INTEGER ia(MXNELT), iwork(MAXIW), ja(MXNELT)
    !     .. External Functions ..
    REAL R1MACH
    EXTERNAL R1MACH
    !     .. External Subroutines ..
    EXTERNAL SCPPLT, SS2Y, SSDBCG, SSDCG, SSDCGN, &
      SSDCGS, SSDGMR, SSDOMN, SSGS, SSICCG, SSILUR, SSJAC, &
      SSLUBC, SSLUCN, SSLUCS, SSLUGM, SSLUOM
    !     .. Intrinsic Functions ..
    INTRINSIC MAX, REAL
    !     .. Common blocks ..
    COMMON /SSLBLK/ SOLn
    !
    !     The following lines are for the braindamaged Sun FPE handler.
    !
    !$$$      integer oldmode, fpmode
    !***FIRST EXECUTABLE STATEMENT  SLAPQC
    !$$$      oldmode = fpmode( 62464 )
    !
    !     Maximum problem sizes.
    !
    neltmx = MXNELT
    nmax = MAXN
    leniw = MAXIW
    lenw = MAXRW
    !
    !     Set some input data.
    !
    n = nmax
    itmax = n
    factor = 1.2E0
    !
    !     Set to print intermediate results if KPRINT.GE.3.
    !
    IF ( Kprint<3 ) THEN
      iunit = 0
    ELSE
      iunit = Lun
    ENDIF
    !
    !     Set the Error tolerance to depend on the machine epsilon.
    !
    tol = MAX(1.0E3*R1MACH(3),1.0E-6)
    nfail = 0
    !
    !     Test routines using various convergence criteria.
    !
    DO kase = 1, 3
      IF ( kase==1.OR.kase==2 ) itol = kase
      IF ( kase==3 ) itol = 11
      !
      !         Test routines using nonsymmetric (ISYM=0) and symmetric
      !         storage (ISYM=1).  For ISYM=0 a really non-symmetric matrix
      !         is generated.  The amount of non-symmetry is controlled by
      !         user.
      !
      DO isym = 0, 1
        IF ( Kprint>=2 ) WRITE (Lun,99001) n, kase, isym
        99001 FORMAT ('1'/' Running tests with  N =',I3,',  KASE =',I2,',  ISYM =',&
          I2)
        !
        !         Set up a random matrix.
        !
        CALL SRMGEN(neltmx,factor,ierr,n,nelt,isym,ia,ja,a,f,SOLn,rwork,iwork,&
          iwork(n+1))
        IF ( ierr/=0 ) THEN
          WRITE (Lun,99002) ierr
          !
          99002 FORMAT (/1X,'SLAPQC -- Fatal error ',I1,' generating ',&
            '*RANDOM* Matrix.')
          nfail = nfail + 1
          CYCLE
        ENDIF
        IF ( isym==0 ) THEN
          dens = REAL(nelt)/(n*n)
        ELSE
          dens = REAL(2*nelt)/(n*n)
        ENDIF
        IF ( Kprint>=2 ) THEN
          WRITE (Lun,99003) n, nelt, dens
          99003 FORMAT (/'                * RANDOM Matrix of size',I5,&
            '*'/'                ','Number of non-zeros & Density = ',&
            I5,1P,E16.7)
          WRITE (Lun,99004) tol
          99004 FORMAT ('                Error tolerance = ',1P,E16.7)
        ENDIF
        !
        !         Convert to the SLAP-Column format and
        !         write out matrix in SLAP-Column format, if desired.
        !
        CALL SS2Y(n,nelt,ia,ja,a,isym)
        IF ( Kprint>=4 ) THEN
          WRITE (Lun,99005) (k,ia(k),ja(k),a(k),k=1,nelt)
          99005 FORMAT (/'  ***** SLAP Column Matrix *****'/' Indx   ia   ja     a'/&
            (1X,I4,1X,I4,1X,I4,1X,1P,E16.7))
          CALL SCPPLT(n,nelt,ia,ja,a,isym,Lun)
        ENDIF
        !
        !**********************************************************************
        !                    BEGINNING OF SLAP QUICK TESTS
        !**********************************************************************
        !
        !         * * * * * *   SSJAC   * * * * * *
        !
        IF ( Kprint>=3 ) WRITE (Lun,99008) 'SSJAC ', itol, isym
        CALL VFILL(n,xiter,0.0E0)
        !
        CALL SSJAC(n,f,xiter,nelt,ia,ja,a,isym,itol,tol,2*itmax,iter,err,ierr,&
          iunit,rwork,lenw,iwork,leniw)
        !
        CALL OUTERR('SSJAC ',ierr,Kprint,nfail,Lun,iter,err)
        !
        !         * * * * *  SSGS  * * * * *
        !
        IF ( Kprint>=3 ) WRITE (Lun,99008) 'SSGS  ', itol, isym
        CALL VFILL(n,xiter,0.0E0)
        !
        CALL SSGS(n,f,xiter,nelt,ia,ja,a,isym,itol,tol,itmax,iter,err,ierr,&
          iunit,rwork,lenw,iwork,leniw)
        !
        CALL OUTERR('SSGS  ',ierr,Kprint,nfail,Lun,iter,err)
        !
        !         * * * * * *   SSILUR   * * * * * *
        !
        IF ( Kprint>=3 ) WRITE (Lun,99008) 'SSILUR', itol, isym
        CALL VFILL(n,xiter,0.0E0)
        !
        CALL SSILUR(n,f,xiter,nelt,ia,ja,a,isym,itol,tol,itmax,iter,err,ierr,&
          iunit,rwork,lenw,iwork,leniw)
        !
        CALL OUTERR('SSILUR',ierr,Kprint,nfail,Lun,iter,err)
        !
        !         * * * * * *   SSDCG    * * * * * *
        !
        IF ( isym==1 ) THEN
          IF ( Kprint>=3 ) WRITE (Lun,99008) 'SSDCG', itol, isym
          CALL VFILL(n,xiter,0.0E0)
          !
          CALL SSDCG(n,f,xiter,nelt,ia,ja,a,isym,itol,tol,itmax,iter,err,ierr,&
            iunit,rwork,lenw,iwork,leniw)
          !
          CALL OUTERR('SSDCG ',ierr,Kprint,nfail,Lun,iter,err)
        ENDIF
        !
        !         * * * * * *    SSICCG    * * * * * *
        !
        IF ( isym==1 ) THEN
          IF ( Kprint>=3 ) WRITE (Lun,99008) 'SSICCG', itol, isym
          CALL VFILL(n,xiter,0.0E0)
          !
          CALL SSICCG(n,f,xiter,nelt,ia,ja,a,isym,itol,tol,itmax,iter,err,&
            ierr,iunit,rwork,lenw,iwork,leniw)
          !
          CALL OUTERR('SSICCG',ierr,Kprint,nfail,Lun,iter,err)
        ENDIF
        !
        !         * * * * * *    SSDCGN   * * * * * *
        !
        IF ( Kprint>=3 ) WRITE (Lun,99008) 'SSDCGN', itol, isym
        CALL VFILL(n,xiter,0.0E0)
        !
        CALL SSDCGN(n,f,xiter,nelt,ia,ja,a,isym,itol,tol,itmax,iter,err,ierr,&
          iunit,rwork,lenw,iwork,leniw)
        !
        CALL OUTERR('SSDCGN',ierr,Kprint,nfail,Lun,iter,err)
        !
        !         * * * * * *   SSLUCN   * * * * * *
        !
        IF ( Kprint>=3 ) WRITE (Lun,99008) 'SSLUCN', itol, isym
        CALL VFILL(n,xiter,0.0E0)
        !
        CALL SSLUCN(n,f,xiter,nelt,ia,ja,a,isym,itol,tol,itmax,iter,err,ierr,&
          iunit,rwork,lenw,iwork,leniw)
        !
        CALL OUTERR('SSLUCN',ierr,Kprint,nfail,Lun,iter,err)
        !
        !         * * * * * *    SSDBCG   * * * * * *
        !
        IF ( Kprint>=3 ) WRITE (Lun,99008) 'SSDBCG', itol, isym
        CALL VFILL(n,xiter,0.0E0)
        !
        CALL SSDBCG(n,f,xiter,nelt,ia,ja,a,isym,itol,tol,itmax,iter,err,ierr,&
          iunit,rwork,lenw,iwork,leniw)
        !
        CALL OUTERR('SSDBCG',ierr,Kprint,nfail,Lun,iter,err)
        !
        !         * * * * * *   SSLUBC   * * * * * *
        !
        IF ( Kprint>=3 ) WRITE (Lun,99008) 'SSLUBC', itol, isym
        CALL VFILL(n,xiter,0.0E0)
        !
        CALL SSLUBC(n,f,xiter,nelt,ia,ja,a,isym,itol,tol,itmax,iter,err,ierr,&
          iunit,rwork,lenw,iwork,leniw)
        !
        CALL OUTERR('SSLUBC',ierr,Kprint,nfail,Lun,iter,err)
        !
        !         * * * * * *    SSDCGS   * * * * * *
        !
        IF ( Kprint>=3 ) WRITE (Lun,99008) 'SSDCGS', itol, isym
        CALL VFILL(n,xiter,0.0E0)
        !
        CALL SSDCGS(n,f,xiter,nelt,ia,ja,a,isym,itol,tol,itmax,iter,err,ierr,&
          iunit,rwork,lenw,iwork,leniw)
        !
        CALL OUTERR('SSDCGS',ierr,Kprint,nfail,Lun,iter,err)
        !
        !         * * * * * *   SSLUCS   * * * * * *
        !
        IF ( Kprint>=3 ) WRITE (Lun,99008) 'SSLUCS', itol, isym
        CALL VFILL(n,xiter,0.0E0)
        !
        CALL SSLUCS(n,f,xiter,nelt,ia,ja,a,isym,itol,tol,itmax,iter,err,ierr,&
          iunit,rwork,lenw,iwork,leniw)
        !
        CALL OUTERR('SSLUCS',ierr,Kprint,nfail,Lun,iter,err)
        !
        !         * * * * * *    SSDOMN   * * * * * *
        !
        !VD$ NOVECTOR
        DO nsave = 0, 3
          IF ( Kprint>=3 ) WRITE (Lun,99009) 'SSDOMN', itol, isym, nsave
          CALL VFILL(n,xiter,0.0E0)
          !
          CALL SSDOMN(n,f,xiter,nelt,ia,ja,a,isym,nsave,itol,tol,itmax,iter,&
            err,ierr,iunit,rwork,lenw,iwork,leniw)
          !
          CALL OUTERR('SSDOMN',ierr,Kprint,nfail,Lun,iter,err)
        ENDDO
        !
        !         * * * * * *   SSLUOM   * * * * * *
        !
        !VD$ NOVECTOR
        DO nsave = 0, 3
          IF ( Kprint>=3 ) WRITE (Lun,99009) 'SSLUOM', itol, isym, nsave
          CALL VFILL(n,xiter,0.0E0)
          !
          CALL SSLUOM(n,f,xiter,nelt,ia,ja,a,isym,nsave,itol,tol,itmax,iter,&
            err,ierr,iunit,rwork,lenw,iwork,leniw)
          !
          CALL OUTERR('SSLUOM',ierr,Kprint,nfail,Lun,iter,err)
        ENDDO
        !
        !         * * * * * *   SSDGMR   * * * * * *
        !
        !VD$ NOVECTOR
        DO nsave = 5, 12
          IF ( Kprint>=3 ) WRITE (Lun,99009) 'SSDGMR', itol, isym, nsave
          CALL VFILL(n,xiter,0.0E0)
          itolgm = 0
          !
          CALL SSDGMR(n,f,xiter,nelt,ia,ja,a,isym,nsave,itolgm,tol,itmax,iter,&
            err,ierr,iunit,rwork,lenw,iwork,leniw)
          !
          CALL OUTERR('SSDGMR',ierr,Kprint,nfail,Lun,iter,err)
        ENDDO
        !
        !         * * * * * *   SSLUGM   * * * * * *
        !
        !VD$ NOVECTOR
        DO nsave = 5, 12
          IF ( Kprint>=3 ) WRITE (Lun,99009) 'SSLUGM', itol, isym, nsave
          CALL VFILL(n,xiter,0.0E0)
          !
          CALL SSLUGM(n,f,xiter,nelt,ia,ja,a,isym,nsave,itol,tol,itmax,iter,&
            err,ierr,iunit,rwork,lenw,iwork,leniw)
          !
          CALL OUTERR('SSLUGM',ierr,Kprint,nfail,Lun,iter,err)
        ENDDO
      ENDDO
    ENDDO
    !
    IF ( nfail==0 ) THEN
      Ipass = 1
      IF ( Kprint>=2 ) WRITE (Lun,99006)
      99006 FORMAT ('--------- All single precision SLAP tests passed ','---------')
    ELSE
      Ipass = 0
      IF ( Kprint>=2 ) WRITE (Lun,99007) nfail
      99007 FORMAT ('*********',I3,' single precision SLAP tests failed ',&
        '*********')
    ENDIF
    !
    RETURN
    99008 FORMAT (/1X,A6,' : ITOL = ',I2,'   ISYM = ',I1)
    99009 FORMAT (/1X,A6,' : ITOL = ',I2,'   ISYM = ',I1,' NSAVE = ',I2)
  END SUBROUTINE SLAPQC
  !DECK SRMGEN
  SUBROUTINE SRMGEN(Neltmx,Factor,Ierr,N,Nelt,Isym,Ia,Ja,A,F,Soln,Dsum,Itmp,&
      Idiag)
    IMPLICIT NONE
    !***BEGIN PROLOGUE  SRMGEN
    !***SUBSIDIARY
    !***PURPOSE  This routine generates a "Random" symmetric or
    !            non-symmetric matrix of size N for use in the SLAP
    !            Quick Checks.
    !***LIBRARY   SLATEC (SLAP)
    !***TYPE      SINGLE PRECISION (SRMGEN-S, DRMGEN-D)
    !***AUTHOR  Seager, Mark K., (LLNL)
    !             seager@llnl.gov
    !             Lawrence Livermore National Laboratory
    !             PO BOX 808, L-300
    !             Livermore, CA 94550
    !             (510) 423-3141
    !***DESCRIPTION
    !
    ! *Usage:
    !       INTEGER NELTMX, IERR, N, NELT, ISYM,
    !       INTEGER IA(NELTMX), JA(NELTMX), ITMP(N), IDIAG(N)
    !       REAL    FACTOR, A(NELTMX), F(N), SOLN(N), DSUM(N)
    !
    !       CALL SRMGEN( NELTMX, FACTOR, IERR, N, NELT, ISYM,
    !      $     IA, JA, A, F, SOLN, DSUM, ITMP, IDIAG )
    !
    ! *Arguments:
    !
    ! NELTMX :IN       Integer.
    !         Maximum number of non-zeros that can be created by this
    !         routine for storage in the IA, JA, A arrays,  see below.
    ! FACTOR :IN       Real.
    !         Non-zeros in the upper triangle are set to FACTOR times
    !         the corresponding entry in the lower triangle when a non-
    !         symmetric matrix is requested (See ISYM, below).
    ! IERR   :OUT      Integer.
    !         Return error flag.
    !             IERR = 0 => everything went OK.
    !                  = 1 => Ran out of space trying to create matrix.
    !                         Set NELTMX to something larger and retry.
    ! N      :IN       Integer.
    !         Size of the linear system to generate (number of unknowns).
    ! NELT   :OUT      Integer.
    !         Number of non-zeros stored in the IA, JA, A arrays, see below.
    ! ISYM   :IN       Integer.
    !         Flag to indicate the type of matrix to generate:
    !             ISYM = 0 => Non-Symmetric Matrix (See FACTOR, above).
    !                  = 1 => Symmetric Matrix.
    ! IA     :OUT      Integer IA(NELTMX).
    !         Stores the row indices for the non-zeros.
    ! JA     :OUT      Integer JA(NELTMX).
    !         Stores the column indices for the non-zeros.
    ! A      :OUT      Real A(NELTMX).
    !         Stores the values of the non-zeros.
    ! F      :OUT      Real F(N).
    !         The right hand side of the linear system.  Obtained by
    !         multiplying the matrix times SOLN, see below.
    ! SOLN   :OUT      Real SOLN(N).
    !         The true solution to the linear system.  Each component is
    !         chosen at random (0.0<SOLN(I)<1.0, I=1,N)
    ! DSUM   :WORK     Real DSUM(N).
    ! ITMP   :WORK     Integer ITMP(N).
    ! IDIAG  :WORK     Integer IDIAG(N).
    !
    ! *Description
    !         The matrix is generated by choosing a random number of
    !         entries for each column and then chosing negative random
    !         numbers for each off diagonal.   The diagonal elements
    !         are chosen to be positive and large enough so the matrix
    !         is slightly diagonally dominant.  The lower triangle of
    !         the matrix is generated and if isym.eq.0 (all matrix elements
    !         stored) the upper triangle elements are chosen so that they
    !         are FACTOR times the corresponding lower triangular element.
    !
    !***ROUTINES CALLED  ISMPL, RAND
    !***REVISION HISTORY  (YYMMDD)
    !   881120  DATE WRITTEN
    !   890919  Replaced SMPL with ISMPL.  (MKS)
    !   890920  Minor changes to reduce single/double differences.  (FNF)
    !   920511  Added complete declaration section.  (WRB)
    !***END PROLOGUE  SRMGEN
    !     .. Scalar Arguments ..
    REAL Factor
    INTEGER Ierr, Isym, N, Nelt, Neltmx
    !     .. Array Arguments ..
    REAL A(Neltmx), Dsum(N), F(N), Soln(N)
    INTEGER Ia(Neltmx), Idiag(N), Itmp(N), Ja(Neltmx)
    !     .. Local Scalars ..
    REAL dummy
    INTEGER i, icol, inum, irow, iseed, k, nl
    !     .. External Functions ..
    REAL RAND
    EXTERNAL RAND
    !     .. External Subroutines ..
    EXTERNAL ISMPL
    !     .. Intrinsic Functions ..
    INTRINSIC INT
    !***FIRST EXECUTABLE STATEMENT  SRMGEN
    !
    !     Start by setting the random number generator seed.  This is done
    !     for reproducibility in debugging.
    !
    !     Remove the seed setting call for production testing.
    !
    !     Note:  Double precision version did not work properly with
    !            certain compilers with literal arguments to RAND.
    !
    dummy = 16381.0
    iseed = INT( RAND(dummy) )
    Ierr = 0
    DO i = 1, N
      Idiag(i) = 0
      Dsum(i) = -1.0E0
    ENDDO
    dummy = 0.0
    Nelt = 0
    !
    !     Set the matrix elements.
    !     Loop over the columns.
    !
    !VD$ NOCONCUR
    DO icol = 1, N
      nl = N + 1 - icol
      !
      !         To keep things sparse divide by two, three or four or ...
      !
      inum = (INT(RAND(dummy)*nl)+1)/3
      CALL ISMPL(nl,inum,Itmp)
      !
      !         Set up this column (and row, if non-symmetric structure).
      !VD$ NOVECTOR
      !VD$ NOCONCUR
      DO irow = 1, inum
        Nelt = Nelt + 1
        IF ( Nelt>Neltmx ) THEN
          Ierr = 1
          RETURN
        ENDIF
        Ia(Nelt) = N + 1 - Itmp(irow)
        Ja(Nelt) = icol
        IF ( Ia(Nelt)==icol ) THEN
          Idiag(icol) = Nelt
        ELSE
          A(Nelt) = -RAND(dummy)
          Dsum(icol) = Dsum(icol) + A(Nelt)
          IF ( Isym==0 ) THEN
            !
            !         Copy this element into upper triangle.
            !
            Nelt = Nelt + 1
            IF ( Nelt>Neltmx ) THEN
              Ierr = 1
              RETURN
            ENDIF
            Ia(Nelt) = icol
            Ja(Nelt) = Ia(Nelt-1)
            A(Nelt) = A(Nelt-1)*Factor
            Dsum(Ja(Nelt)) = Dsum(Ja(Nelt)) + A(Nelt)
          ELSE
            Dsum(Ia(Nelt)) = Dsum(Ia(Nelt)) + A(Nelt)
          ENDIF
        ENDIF
      ENDDO
      IF ( Idiag(icol)==0 ) THEN
        !
        !           Add a diagonal to the column.
        !
        Nelt = Nelt + 1
        IF ( Nelt>Neltmx ) THEN
          Ierr = 1
          RETURN
        ENDIF
        Idiag(icol) = Nelt
        A(Nelt) = 0.0E0
        Ia(Nelt) = icol
        Ja(Nelt) = icol
      ENDIF
    ENDDO
    !
    !         Clean up the diagonals.
    !
    !VD$ NODEPCHK
    !LLL. OPTION ASSERT (NOHAZARD)
    !DIR$ IVDEP
    DO i = 1, N
      A(Idiag(i)) = -1.0001E0*Dsum(i)
    ENDDO
    !
    !         Set a random solution and determine the right-hand side.
    !
    !VD$ NOVECTOR
    !VD$ NOCONCUR
    DO i = 1, N
      Soln(i) = RAND(dummy)
      F(i) = 0.0E0
    ENDDO
    !
    !VD$ NOVECTOR
    !VD$ NOCONCUR
    DO k = 1, Nelt
      F(Ia(k)) = F(Ia(k)) + A(k)*Soln(Ja(k))
      IF ( Isym/=0.AND.Ia(k)/=Ja(k) ) F(Ja(k)) = F(Ja(k)) + A(k)*Soln(Ia(k))
    ENDDO
    !------------- LAST LINE OF SRMGEN FOLLOWS ----------------------------
  END SUBROUTINE SRMGEN
  !DECK VFILL
  SUBROUTINE VFILL(N,V,Val)
    IMPLICIT NONE
    !***BEGIN PROLOGUE  VFILL
    !***SUBSIDIARY
    !***PURPOSE  Fill a vector with a value.
    !***LIBRARY   SLATEC (SLAP)
    !***TYPE      SINGLE PRECISION (VFILL-S, DFILL-D)
    !***AUTHOR  Seager, Mark K., (LLNL)
    !             Lawrence Livermore National Laboratory
    !             PO BOX 808, L-300
    !             Livermore, CA 94550 (510) 423-3141
    !             seager@llnl.gov
    !***DESCRIPTION
    !
    ! *Usage:
    !     INTEGER  N
    !     REAL     V(N), VAL
    !
    !     CALL VFILL( N, V, VAL )
    !
    ! *Arguments:
    ! N      :IN       Integer.
    !         Length of the vector
    ! V      :OUT      Real V(N).
    !         Vector to be set.
    ! VAL    :IN       Real.
    !         Value to seed the vector with.
    !***REFERENCES  (NONE)
    !***ROUTINES CALLED  (NONE)
    !***REVISION HISTORY  (YYMMDD)
    !   871119  DATE WRITTEN
    !   881213  Previous REVISION DATE
    !   890920  Converted prologue to SLATEC 4.0 format.  (FNF)
    !   920511  Added complete declaration section.  (WRB)
    !***END PROLOGUE  VFILL
    !     .. Scalar Arguments ..
    REAL Val
    INTEGER N
    !     .. Array Arguments ..
    REAL V(*)
    !     .. Local Scalars ..
    INTEGER i, is, nr
    !     .. Intrinsic Functions ..
    INTRINSIC MOD
    !***FIRST EXECUTABLE STATEMENT  VFILL
    IF ( N<=0 ) RETURN
    nr = MOD(N,4)
    !
    !         The following construct assumes a zero pass do loop.
    !
    is = 1
    SELECT CASE (nr+1)
      CASE (1)
      CASE (2)
        is = 2
        V(1) = Val
      CASE (3)
        is = 3
        V(1) = Val
        V(2) = Val
      CASE DEFAULT
        is = 4
        V(1) = Val
        V(2) = Val
        V(3) = Val
    END SELECT
    DO i = is, N, 4
      V(i) = Val
      V(i+1) = Val
      V(i+2) = Val
      V(i+3) = Val
    ENDDO
    !------------- LAST LINE OF VFILL FOLLOWS -----------------------------
  END SUBROUTINE VFILL
  !DECK OUTERR
  SUBROUTINE OUTERR(Method,Ierr,Iout,Nfail,Istdo,Iter,Err)
    IMPLICIT NONE
    !***BEGIN PROLOGUE  OUTERR
    !***SUBSIDIARY
    !***PURPOSE  Output error messages for the SLAP Quick Check.
    !***LIBRARY   SLATEC (SLAP)
    !***TYPE      SINGLE PRECISION (OUTERR-S, DUTERR-D)
    !***AUTHOR  Seager, Mark K., (LLNL)
    !             Lawrence Livermore National Laboratory
    !             PO BOX 808, L-300
    !             Livermore, CA 94550 (510) 423-3141
    !             seager@llnl.gov
    !***ROUTINES CALLED  (NONE)
    !***REVISION HISTORY  (YYMMDD)
    !   881010  DATE WRITTEN
    !   881213  Previous REVISION DATE
    !   890920  Converted prologue to SLATEC 4.0 format.  (FNF)
    !   920511  Added complete declaration section.  (WRB)
    !   921021  Added 1P's to output formats.  (FNF)
    !***END PROLOGUE  OUTERR
    !     .. Scalar Arguments ..
    REAL Err
    INTEGER Ierr, Iout, Istdo, Iter, Nfail
    CHARACTER Method*6
    !***FIRST EXECUTABLE STATEMENT  OUTERR
    IF ( Ierr/=0 ) Nfail = Nfail + 1
    IF ( Iout==1.AND.Ierr/=0 ) THEN
      WRITE (Istdo,99001) Method
      99001 FORMAT (1X,A6,' : **** FAILURE ****')
    ENDIF
    IF ( Iout==2 ) THEN
      IF ( Ierr==0 ) THEN
        WRITE (Istdo,99002) Method
        99002 FORMAT (1X,A6,' : **** PASSED  ****')
      ELSE
        WRITE (Istdo,99004) Method, Ierr, Iter, Err
      ENDIF
    ENDIF
    IF ( Iout>=3 ) THEN
      IF ( Ierr==0 ) THEN
        WRITE (Istdo,99003) Method, Ierr, Iter, Err
        99003 FORMAT (' ***************** PASSED ***********************'/' **** ',&
          A6,' Quick Test PASSED: IERR = ',I5,&
          ' ****'/' ***************** PASSED ***********************'/&
          ' Iteration Count = ',I3,' Stop Test = ',1P,E12.6)
      ELSE
        WRITE (Istdo,99004) Method, Ierr, Iter, Err
      ENDIF
    ENDIF
    RETURN
    99004 FORMAT (' **************** WARNING ***********************'/' **** ',A6,&
      ' Quick Test FAILED: IERR = ',I5,&
      ' ****'/' **************** WARNING ***********************'/&
      ' Iteration Count = ',I3,' Stop Test = ',1P,E12.6)
    !------------- LAST LINE OF OUTERR FOLLOWS ----------------------------
  END SUBROUTINE OUTERR
END MODULE TEST25_MOD
!DECK TEST25
PROGRAM TEST25
  USE TEST25_MOD
  IMPLICIT NONE
  !***BEGIN PROLOGUE  TEST25
  !***PURPOSE  Driver for testing SLATEC subprograms.
  !***LIBRARY   SLATEC
  !***CATEGORY  D2A4, D2B4
  !***TYPE      SINGLE PRECISION (TEST25-S, TEST26-D)
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
  !       Single precision SLAP subprograms
  !
  !***REFERENCES  Fong, Kirby W., Jefferson, Thomas H., Suyehiro,
  !                 Tokihiko, Walton, Lee, Guidelines to the SLATEC Common
  !                 Mathematical Library, March 21, 1989.
  !***ROUTINES CALLED  I1MACH, SLAPQC, XERMAX, XSETF, XSETUN
  !***REVISION HISTORY  (YYMMDD)
  !   920401  DATE WRITTEN
  !   920511  Added complete declaration section.  (WRB)
  !***END PROLOGUE  TEST25
  !     .. Local Scalars ..
  INTEGER ipass, kprint, lin, lun, nfail
  !     .. External Functions ..
  INTEGER I1MACH
  EXTERNAL I1MACH
  !     .. External Subroutines ..
  EXTERNAL XERMAX, XSETF, XSETUN
  !***FIRST EXECUTABLE STATEMENT  TEST25
  lun = I1MACH(2)
  lin = I1MACH(1)
  nfail = 0
  !
  !     Read KPRINT parameter
  !
  READ (lin,'(I1)') kprint
  CALL XSETUN(lun)
  IF ( kprint<=1 ) THEN
    CALL XSETF(0)
  ELSE
    CALL XSETF(1)
  ENDIF
  CALL XERMAX(1000)
  !
  !     Test SLAP (single precision)
  !
  CALL SLAPQC(lun,kprint,ipass)
  IF ( ipass==0 ) nfail = nfail + 1
  !
  !     Write PASS or FAIL message
  !
  IF ( nfail==0 ) THEN
    WRITE (lun,99001)
    99001 FORMAT (/' --------------TEST25 PASSED ALL TESTS----------------')
  ELSE
    WRITE (lun,99002) nfail
    99002 FORMAT (/' ************* WARNING -- ',I5,&
      ' TEST(S) FAILED IN PROGRAM TEST25 *************')
  ENDIF
  STOP
END PROGRAM TEST25