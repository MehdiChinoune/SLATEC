MODULE TEST26_MOD
  USE service, ONLY : SP, DP
  IMPLICIT NONE

CONTAINS
  !** DLAPQC
  SUBROUTINE DLAPQC(Lun,Kprint,Ipass)
    !>
    !  Quick check for testing Sparse Linear Algebra Package
    !            (SLAP) Version 2.0.2.
    !***
    ! **Library:**   SLATEC (SLAP)
    !***
    ! **Category:**  D2A4, D2B4
    !***
    ! **Type:**      DOUBLE PRECISION (SLAPQC-S, DLAPQC-D)
    !***
    ! **Keywords:**  QUICK CHECK, SLAP
    !***
    ! **Author:**  Mark K. Seager (LLNL)
    !             seager@llnl.gov
    !             Lawrence Livermore National Laboratory
    !             PO BOX 808, L-300
    !             Livermore, CA 94550
    !             (510) 423-3141
    !***
    ! **Description:**
    !
    !- Arguments:
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
    !- Description:
    !         This is a SLATEC Quick Check program to test the *SLAP*
    !         Version 2.0.2 package.  It generates a "random" matrix (See
    !         DRMGEN) and then runs all the various methods with all the
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
    !***
    ! **References:**  (NONE)
    !***
    ! **Routines called:**  D1MACH, DCPPLT, DFILL, DRMGEN, DS2Y, DSDBCG, DSDCG,
    !                    DSDCGN, DSDCGS, DSDGMR, DSDOMN, DSGS, DSICCG,
    !                    DSILUR, DSJAC, DSLUBC, DSLUCN, DSLUCS, DSLUGM,
    !                    DSLUOM, DUTERR, XERMAX, XSETF, XSETUN
    !***
    ! COMMON BLOCKS    DSLBLK

    !* REVISION HISTORY  (YYMMDD)
    !   880601  DATE WRITTEN
    !   881213  Revised to meet the new SLATEC prologue standards.
    !   890920  Modified to reduce single/double differences and to meet
    !           SLATEC standards, as requested at July 1989 CML Meeting.
    !   891003  Reduced MAXN to a more reasonable size for quick check.
    !   920401  Made routine a SUBROUTINE and made necessary changes to
    !           interface with a SLATEC quick check driver.  (WRB)
    !   920407  COMMON BLOCK renamed DSLBLK.  (WRB)
    !   920511  Added complete declaration section.  (WRB)
    !   920602  Eliminated unnecessary variables IOUT and ISTDO and made
    !           various cosmetic changes.  (FNF)
    !   920602  Reduced problem size for a shorter-running test and
    !           corrected lower limit in "DO 80" statement.  (FNF)
    !   921021  Changed E's to 1P,D's in output formats.  (FNF)
    USE DSLBLK, ONLY : soln_com
    USE slatec, ONLY : D1MACH, DCPPLT, DS2Y, DSDBCG, DSDCG, DSDCGN, DSDCGS, DSDGMR, &
      DSDOMN, DSGS, DSICCG, DSILUR, DSJAC, DSLUBC, DSLUCN, DSLUCS, DSLUGM, DSLUOM, &
      XERMAX, XSETF, XSETUN
    !
    !     The problem size, MAXN, should be large enough that the
    !     iterative methods do 10-15 iterations, just to be sure that
    !     the truncated methods run to the end of their ropes and enter
    !     their error recovery mode.  Thus, for a more thorough test change
    !     the following PARAMETER statement to:
    !     , PARAMETER :: MAXN=69, MXNELT=5000, MAXIW=5000, MAXRW=5000)
    !
    !     .. Parameters ..
    INTEGER, PARAMETER :: MAXN = 25, MXNELT = 500, MAXIW = 1000, MAXRW = 1000
    !     .. Scalar Arguments ..
    INTEGER Ipass, Kprint, Lun
    !     .. Local Scalars ..
    REAL(DP) :: dens, err, factor, tol
    INTEGER ierr, isym, iter, itmax, itol, itolgm, iunit, k, kase, &
      leniw, lenw, n, nelt, neltmx, nfail, nmax, nsave
    !     .. Local Arrays ..
    REAL(DP) :: a(MXNELT), f(MAXN), rwork(MAXRW), xiter(MAXN)
    INTEGER ia(MXNELT), iwork(MAXIW), ja(MXNELT)
    !     .. Intrinsic Functions ..
    INTRINSIC MAX, REAL
    !
    !     The following lines are for the braindamaged Sun FPE handler.
    !
    !$$$      integer oldmode, fpmode
    !* FIRST EXECUTABLE STATEMENT  DLAPQC
    !$$$      oldmode = fpmode( 62464 )
    !
    !     Maximum problem sizes.
    !
    neltmx = MXNELT
    nmax = MAXN
    leniw = MAXIW
    lenw = MAXRW
    ALLOCATE( soln_com(nmax) )
    !
    !     Set some input data.
    !
    n = nmax
    itmax = n
    factor = 1.2D0
    !
    !     Set to print intermediate results if KPRINT.GE.3.
    !
    IF ( Kprint<3 ) THEN
      iunit = 0
    ELSE
      iunit = Lun
    END IF
    !
    !     Set the Error tolerance to depend on the machine epsilon.
    !
    tol = MAX(1.0D3*D1MACH(3),1.0D-6)
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
        99001 FORMAT ('1'/' Running tests with  N =',I3,',  KASE =',I2,',  ISYM =',I2)
        !
        !         Set up a random matrix.
        !
        CALL DRMGEN(neltmx,factor,ierr,n,nelt,isym,ia,ja,a,f,soln_com,rwork,&
          iwork,iwork(n+1))
        IF ( ierr/=0 ) THEN
          WRITE (Lun,99002) ierr
          !
          99002 FORMAT (/1X,'DLAPQC -- Fatal error ',I1,' generating ',&
            '*RANDOM* Matrix.')
          nfail = nfail + 1
          CYCLE
        END IF
        IF ( isym==0 ) THEN
          dens = REAL(nelt)/(n*n)
        ELSE
          dens = REAL(2*nelt)/(n*n)
        END IF
        IF ( Kprint>=2 ) THEN
          WRITE (Lun,99003) n, nelt, dens
          99003 FORMAT (/'                * RANDOM Matrix of size',I5,&
            '*'/'                ','Number of non-zeros & Density = ',&
            I5,1P,D16.7)
          WRITE (Lun,99004) tol
          99004 FORMAT ('                Error tolerance = ',1P,D16.7)
        END IF
        !
        !         Convert to the SLAP-Column format and
        !         write out matrix in SLAP-Column format, if desired.
        !
        CALL DS2Y(n,nelt,ia,ja,a)
        IF ( Kprint>=4 ) THEN
          WRITE (Lun,99005) (k,ia(k),ja(k),a(k),k=1,nelt)
          99005 FORMAT (/'  ***** SLAP Column Matrix *****'/' Indx   ia   ja     a'/&
            (1X,I4,1X,I4,1X,I4,1X,1P,D16.7))
          CALL DCPPLT(n,nelt,ia,ja,a,isym,Lun)
        END IF
        !
        !**********************************************************************
        !                    BEGINNING OF SLAP QUICK TESTS
        !**********************************************************************
        !
        !         * * * * * *   DSJAC   * * * * * *
        !
        IF ( Kprint>=3 ) WRITE (Lun,99008) 'DSJAC ', itol, isym
        CALL DFILL(n,xiter,0.0D0)
        !
        CALL DSJAC(n,f,xiter,nelt,ia,ja,a,isym,itol,tol,2*itmax,iter,err,ierr,&
          iunit,rwork,lenw,iwork,leniw)
        !
        CALL DUTERR('DSJAC ',ierr,Kprint,nfail,Lun,iter,err)
        !
        !         * * * * *  DSGS  * * * * *
        !
        IF ( Kprint>=3 ) WRITE (Lun,99008) 'DSGS  ', itol, isym
        CALL DFILL(n,xiter,0.0D0)
        !
        CALL DSGS(n,f,xiter,nelt,ia,ja,a,isym,itol,tol,itmax,iter,err,ierr,&
          iunit,rwork,lenw,iwork,leniw)
        !
        CALL DUTERR('DSGS  ',ierr,Kprint,nfail,Lun,iter,err)
        !
        !         * * * * * *   DSILUR   * * * * * *
        !
        IF ( Kprint>=3 ) WRITE (Lun,99008) 'DSILUR', itol, isym
        CALL DFILL(n,xiter,0.0D0)
        !
        CALL DSILUR(n,f,xiter,nelt,ia,ja,a,isym,itol,tol,itmax,iter,err,ierr,&
          iunit,rwork,lenw,iwork,leniw)
        !
        CALL DUTERR('DSILUR',ierr,Kprint,nfail,Lun,iter,err)
        !
        !         * * * * * *   DSDCG    * * * * * *
        !
        IF ( isym==1 ) THEN
          IF ( Kprint>=3 ) WRITE (Lun,99008) 'DSDCG', itol, isym
          CALL DFILL(n,xiter,0.0D0)
          !
          CALL DSDCG(n,f,xiter,nelt,ia,ja,a,isym,itol,tol,itmax,iter,err,ierr,&
            iunit,rwork,lenw,iwork,leniw)
          !
          CALL DUTERR('DSDCG ',ierr,Kprint,nfail,Lun,iter,err)
        END IF
        !
        !         * * * * * *    DSICCG    * * * * * *
        !
        IF ( isym==1 ) THEN
          IF ( Kprint>=3 ) WRITE (Lun,99008) 'DSICCG', itol, isym
          CALL DFILL(n,xiter,0.0D0)
          !
          CALL DSICCG(n,f,xiter,nelt,ia,ja,a,isym,itol,tol,itmax,iter,err,&
            ierr,iunit,rwork,lenw,iwork,leniw)
          !
          CALL DUTERR('DSICCG',ierr,Kprint,nfail,Lun,iter,err)
        END IF
        !
        !         * * * * * *    DSDCGN   * * * * * *
        !
        IF ( Kprint>=3 ) WRITE (Lun,99008) 'DSDCGN', itol, isym
        CALL DFILL(n,xiter,0.0D0)
        !
        CALL DSDCGN(n,f,xiter,nelt,ia,ja,a,isym,itol,tol,itmax,iter,err,ierr,&
          iunit,rwork,lenw,iwork,leniw)
        !
        CALL DUTERR('DSDCGN',ierr,Kprint,nfail,Lun,iter,err)
        !
        !         * * * * * *   DSLUCN   * * * * * *
        !
        IF ( Kprint>=3 ) WRITE (Lun,99008) 'DSLUCN', itol, isym
        CALL DFILL(n,xiter,0.0D0)
        !
        CALL DSLUCN(n,f,xiter,nelt,ia,ja,a,isym,itol,tol,itmax,iter,err,ierr,&
          iunit,rwork,lenw,iwork,leniw)
        !
        CALL DUTERR('DSLUCN',ierr,Kprint,nfail,Lun,iter,err)
        !
        !         * * * * * *    DSDBCG   * * * * * *
        !
        IF ( Kprint>=3 ) WRITE (Lun,99008) 'DSDBCG', itol, isym
        CALL DFILL(n,xiter,0.0D0)
        !
        CALL DSDBCG(n,f,xiter,nelt,ia,ja,a,isym,itol,tol,itmax,iter,err,ierr,&
          iunit,rwork,lenw,iwork,leniw)
        !
        CALL DUTERR('DSDBCG',ierr,Kprint,nfail,Lun,iter,err)
        !
        !         * * * * * *   DSLUBC   * * * * * *
        !
        IF ( Kprint>=3 ) WRITE (Lun,99008) 'DSLUBC', itol, isym
        CALL DFILL(n,xiter,0.0D0)
        !
        CALL DSLUBC(n,f,xiter,nelt,ia,ja,a,isym,itol,tol,itmax,iter,err,ierr,&
          iunit,rwork,lenw,iwork,leniw)
        !
        CALL DUTERR('DSLUBC',ierr,Kprint,nfail,Lun,iter,err)
        !
        !         * * * * * *    DSDCGS   * * * * * *
        !
        IF ( Kprint>=3 ) WRITE (Lun,99008) 'DSDCGS', itol, isym
        CALL DFILL(n,xiter,0.0D0)
        !
        CALL DSDCGS(n,f,xiter,nelt,ia,ja,a,isym,itol,tol,itmax,iter,err,ierr,&
          iunit,rwork,lenw,iwork,leniw)
        !
        CALL DUTERR('DSDCGS',ierr,Kprint,nfail,Lun,iter,err)
        !
        !         * * * * * *   DSLUCS   * * * * * *
        !
        IF ( Kprint>=3 ) WRITE (Lun,99008) 'DSLUCS', itol, isym
        CALL DFILL(n,xiter,0.0D0)
        !
        CALL DSLUCS(n,f,xiter,nelt,ia,ja,a,isym,itol,tol,itmax,iter,err,ierr,&
          iunit,rwork,lenw,iwork,leniw)
        !
        CALL DUTERR('DSLUCS',ierr,Kprint,nfail,Lun,iter,err)
        !
        !         * * * * * *    DSDOMN   * * * * * *
        !
        DO nsave = 0, 3
          IF ( Kprint>=3 ) WRITE (Lun,99009) 'DSDOMN', itol, isym, nsave
          CALL DFILL(n,xiter,0.0D0)
          !
          CALL DSDOMN(n,f,xiter,nelt,ia,ja,a,isym,nsave,itol,tol,itmax,iter,&
            err,ierr,iunit,rwork,lenw,iwork,leniw)
          !
          CALL DUTERR('DSDOMN',ierr,Kprint,nfail,Lun,iter,err)
        END DO
        !
        !         * * * * * *   DSLUOM   * * * * * *
        !
        DO nsave = 0, 3
          IF ( Kprint>=3 ) WRITE (Lun,99009) 'DSLUOM', itol, isym, nsave
          CALL DFILL(n,xiter,0.0D0)
          !
          CALL DSLUOM(n,f,xiter,nelt,ia,ja,a,isym,nsave,itol,tol,itmax,iter,&
            err,ierr,iunit,rwork,lenw,iwork,leniw)
          !
          CALL DUTERR('DSLUOM',ierr,Kprint,nfail,Lun,iter,err)
        END DO
        !
        !         * * * * * *   DSDGMR   * * * * * *
        !
        DO nsave = 5, 12
          IF ( Kprint>=3 ) WRITE (Lun,99009) 'DSDGMR', itol, isym, nsave
          CALL DFILL(n,xiter,0.0D0)
          itolgm = 0
          !
          CALL DSDGMR(n,f,xiter,nelt,ia,ja,a,isym,nsave,tol,itmax,iter,&
            err,ierr,iunit,rwork,lenw,iwork,leniw)
          !
          CALL DUTERR('DSDGMR',ierr,Kprint,nfail,Lun,iter,err)
        END DO
        !
        !         * * * * * *   DSLUGM   * * * * * *
        !
        DO nsave = 5, 12
          IF ( Kprint>=3 ) WRITE (Lun,99009) 'DSLUGM', itol, isym, nsave
          CALL DFILL(n,xiter,0.0D0)
          !
          CALL DSLUGM(n,f,xiter,nelt,ia,ja,a,isym,nsave,tol,itmax,iter,&
            err,ierr,iunit,rwork,lenw,iwork,leniw)
          !
          CALL DUTERR('DSLUGM',ierr,Kprint,nfail,Lun,iter,err)
        END DO
      END DO
    END DO
    !
    IF ( nfail==0 ) THEN
      Ipass = 1
      IF ( Kprint>=2 ) WRITE (Lun,99006)
      99006 FORMAT ('--------- All double precision SLAP tests passed ','---------')
    ELSE
      Ipass = 0
      IF ( Kprint>=2 ) WRITE (Lun,99007) nfail
      99007 FORMAT ('*********',I3,' double precision SLAP tests failed ',&
        '*********')
    END IF
    !
    RETURN
    99008 FORMAT (/1X,A6,' : ITOL = ',I2,'   ISYM = ',I1)
    99009 FORMAT (/1X,A6,' : ITOL = ',I2,'   ISYM = ',I1,' NSAVE = ',I2)
  END SUBROUTINE DLAPQC
  !** DRMGEN
  SUBROUTINE DRMGEN(Neltmx,Factor,Ierr,N,Nelt,Isym,Ia,Ja,A,F,Soln,Dsum,Itmp,Idiag)
    !>
    !  This routine generates a "Random" symmetric or non-symmetric matrix
    !  of size N for use in the SLAP Quick Checks.
    !***
    ! **Library:**   SLATEC (SLAP)
    !***
    ! **Type:**      DOUBLE PRECISION (SRMGEN-S, DRMGEN-D)
    !***
    ! **Author:**  Seager, Mark K., (LLNL)
    !             seager@llnl.gov
    !             Lawrence Livermore National Laboratory
    !             PO BOX 808, L-300
    !             Livermore, CA 94550
    !             (510) 423-3141
    !***
    ! **Description:**
    !
    !- Usage:
    !       INTEGER NELTMX, IERR, N, NELT, ISYM,
    !       INTEGER IA(NELTMX), JA(NELTMX), ITMP(N), IDIAG(N)
    !       DOUBLE PRECISION FACTOR, A(NELTMX), F(N), SOLN(N), DSUM(N)
    !
    !       CALL DRMGEN( NELTMX, FACTOR, IERR, N, NELT, ISYM,
    !      $     IA, JA, A, F, SOLN, DSUM, ITMP, IDIAG )
    !
    !- Arguments:
    !
    ! NELTMX :IN       Integer.
    !         Maximum number of non-zeros that can be created by this
    !         routine for storage in the IA, JA, A arrays,  see below.
    ! FACTOR :IN       Double Precision.
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
    ! A      :OUT      Double Precision A(NELTMX).
    !         Stores the values of the non-zeros.
    ! F      :OUT      Double Precision F(N).
    !         The right hand side of the linear system.  Obtained by
    !         multiplying the matrix times SOLN, see below.
    ! SOLN   :OUT      Double Precision SOLN(N).
    !         The true solution to the linear system.  Each component is
    !         chosen at random (0.0<SOLN(I)<1.0, I=1,N)
    ! DSUM   :WORK     Double Precision DSUM(N).
    ! ITMP   :WORK     Integer ITMP(N).
    ! IDIAG  :WORK     Integer IDIAG(N).
    !
    !- Description
    !         The matrix is generated by choosing a random number of
    !         entries for each column and then chosing negative random
    !         numbers for each off diagonal.   The diagonal elements
    !         are chosen to be positive and large enough so the matrix
    !         is slightly diagonally dominant.  The lower triangle of
    !         the matrix is generated and if isym.eq.0 (all matrix elements
    !         stored) the upper triangle elements are chosen so that they
    !         are FACTOR times the corresponding lower triangular element.
    !
    !***
    ! **Routines called:**  ISMPL, RAND

    !* REVISION HISTORY  (YYMMDD)
    !   881120  DATE WRITTEN
    !   890919  Replaced DMPL with ISMPL.  (MKS)
    !   890920  Minor changes to reduce single/double differences.  (FNF)
    !   920511  Added complete declaration section.  (WRB)
    USE slatec, ONLY : RAND
    USE common_mod, ONLY : ISMPL
    !     .. Scalar Arguments ..
    REAL(DP) :: Factor
    INTEGER Ierr, Isym, N, Nelt, Neltmx
    !     .. Array Arguments ..
    REAL(DP) :: A(Neltmx), Dsum(N), F(N), SOLn(N)
    INTEGER Ia(Neltmx), Idiag(N), Itmp(N), Ja(Neltmx)
    !     .. Local Scalars ..
    REAL(SP) dummy
    INTEGER i, icol, inum, irow, iseed, k, nl
    !     .. Intrinsic Functions ..
    INTRINSIC INT
    !* FIRST EXECUTABLE STATEMENT  DRMGEN
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
      Dsum(i) = -1.0D0
    END DO
    dummy = 0.0
    Nelt = 0
    !
    !     Set the matrix elements.
    !     Loop over the columns.
    !
    DO icol = 1, N
      nl = N + 1 - icol
      !
      !         To keep things sparse divide by two, three or four or ...
      !
      inum = (INT(RAND(dummy)*nl)+1)/3
      CALL ISMPL(nl,inum,Itmp)
      !
      !         Set up this column (and row, if non-symmetric structure).
      DO irow = 1, inum
        Nelt = Nelt + 1
        IF ( Nelt>Neltmx ) THEN
          Ierr = 1
          RETURN
        END IF
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
            END IF
            Ia(Nelt) = icol
            Ja(Nelt) = Ia(Nelt-1)
            A(Nelt) = A(Nelt-1)*Factor
            Dsum(Ja(Nelt)) = Dsum(Ja(Nelt)) + A(Nelt)
          ELSE
            Dsum(Ia(Nelt)) = Dsum(Ia(Nelt)) + A(Nelt)
          END IF
        END IF
      END DO
      IF ( Idiag(icol)==0 ) THEN
        !
        !           Add a diagonal to the column.
        !
        Nelt = Nelt + 1
        IF ( Nelt>Neltmx ) THEN
          Ierr = 1
          RETURN
        END IF
        Idiag(icol) = Nelt
        A(Nelt) = 0.0D0
        Ia(Nelt) = icol
        Ja(Nelt) = icol
      END IF
    END DO
    !
    !         Clean up the diagonals.
    !
    DO i = 1, N
      A(Idiag(i)) = -1.0001D0*Dsum(i)
    END DO
    !
    !         Set a random solution and determine the right-hand side.
    !
    DO i = 1, N
      SOLn(i) = RAND(dummy)
      F(i) = 0.0D0
    END DO
    !
    DO k = 1, Nelt
      F(Ia(k)) = F(Ia(k)) + A(k)*SOLn(Ja(k))
      IF ( Isym/=0.AND.Ia(k)/=Ja(k) ) F(Ja(k)) = F(Ja(k)) + A(k)*SOLn(Ia(k))
    END DO
    !------------- LAST LINE OF DRMGEN FOLLOWS ----------------------------
  END SUBROUTINE DRMGEN
  !** DFILL
  SUBROUTINE DFILL(N,V,Val)
    !>
    !  Fill a vector with a value.
    !***
    ! **Library:**   SLATEC (SLAP)
    !***
    ! **Type:**      DOUBLE PRECISION (VFILL-S, DFILL-D)
    !***
    ! **Author:**  Seager, Mark K., (LLNL)
    !             Lawrence Livermore National Laboratory
    !             PO BOX 808, L-300
    !             Livermore, CA 94550 (510) 423-3141
    !             seager@llnl.gov
    !***
    ! **Description:**
    !
    !- Usage:
    !     INTEGER  N
    !     DOUBLE PRECISION V(N), VAL
    !
    !     CALL DFILL( N, V, VAL )
    !
    !- Arguments:
    ! N      :IN       Integer.
    !         Length of the vector
    ! V      :OUT      Double Precision V(N).
    !         Vector to be set.
    ! VAL    :IN       Double Precision.
    !         Value to seed the vector with.
    !***
    ! **References:**  (NONE)
    !***
    ! **Routines called:**  (NONE)

    !* REVISION HISTORY  (YYMMDD)
    !   890404  DATE WRITTEN
    !   890920  Converted prologue to SLATEC 4.0 format.  (FNF)
    !   920511  Added complete declaration section.  (WRB)

    !     .. Scalar Arguments ..
    REAL(DP) :: Val
    INTEGER N
    !     .. Array Arguments ..
    REAL(DP) :: V(N)
    !     .. Local Scalars ..
    INTEGER i, is, nr
    !     .. Intrinsic Functions ..
    INTRINSIC MOD
    !* FIRST EXECUTABLE STATEMENT  DFILL
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
    END DO
    !------------- LAST LINE OF DFILL FOLLOWS -----------------------------
  END SUBROUTINE DFILL
  !** DUTERR
  SUBROUTINE DUTERR(Method,Ierr,Iout,Nfail,Istdo,Iter,Err)
    !>
    !  Output error messages for the SLAP Quick Check.
    !***
    ! **Library:**   SLATEC (SLAP)
    !***
    ! **Type:**      DOUBLE PRECISION (OUTERR-S, DUTERR-D)
    !***
    ! **Author:**  Seager, Mark K., (LLNL)
    !             Lawrence Livermore National Laboratory
    !             PO BOX 808, L-300
    !             Livermore, CA 94550 (510) 423-3141
    !             seager@llnl.gov
    !***
    ! **Routines called:**  (NONE)

    !* REVISION HISTORY  (YYMMDD)
    !   890404  DATE WRITTEN
    !   890920  Converted prologue to SLATEC 4.0 format.  (FNF)
    !   920511  Added complete declaration section.  (WRB)
    !   921021  Changed E's to 1P,D's in output formats.  (FNF)

    !     .. Scalar Arguments ..
    REAL(DP) :: Err
    INTEGER Ierr, Iout, Istdo, Iter, Nfail
    CHARACTER Method*6
    !* FIRST EXECUTABLE STATEMENT  DUTERR
    IF ( Ierr/=0 ) Nfail = Nfail + 1
    IF ( Iout==1.AND.Ierr/=0 ) THEN
      WRITE (Istdo,99001) Method
      99001 FORMAT (1X,A6,' : **** FAILURE ****')
    END IF
    IF ( Iout==2 ) THEN
      IF ( Ierr==0 ) THEN
        WRITE (Istdo,99002) Method
        99002 FORMAT (1X,A6,' : **** PASSED  ****')
      ELSE
        WRITE (Istdo,99004) Method, Ierr, Iter, Err
      END IF
    END IF
    IF ( Iout>=3 ) THEN
      IF ( Ierr==0 ) THEN
        WRITE (Istdo,99003) Method, Ierr, Iter, Err
        99003 FORMAT (' ***************** PASSED ***********************'/' **** ',&
          A6,' Quick Test PASSED: IERR = ',I5,&
          ' ****'/' ***************** PASSED ***********************'/&
          ' Iteration Count = ',I3,' Stop Test = ',1P,D13.6)
      ELSE
        WRITE (Istdo,99004) Method, Ierr, Iter, Err
      END IF
    END IF
    RETURN
    99004 FORMAT (' **************** WARNING ***********************'/' **** ',A6,&
      ' Quick Test FAILED: IERR = ',I5,&
      ' ****'/' **************** WARNING ***********************'/&
      ' Iteration Count = ',I3,' Stop Test = ',1P,D13.6)
    !------------- LAST LINE OF DUTERR FOLLOWS ----------------------------
  END SUBROUTINE DUTERR
END MODULE TEST26_MOD
!** TEST26
PROGRAM TEST26
  USE TEST26_MOD, ONLY : DLAPQC
  USE slatec, ONLY : I1MACH, XSETF, XSETUN, XERMAX
  USE common_mod, ONLY : GET_ARGUMENT
  IMPLICIT NONE
  !>
  !  Driver for testing SLATEC subprograms.
  !***
  ! **Library:**   SLATEC
  !***
  ! **Category:**  D2A4, D2B4
  !***
  ! **Type:**      DOUBLE PRECISION (TEST25-S, TEST26-D)
  !***
  ! **Keywords:**  QUICK CHECK DRIVER
  !***
  ! **Author:**  SLATEC Common Mathematical Library Committee
  !***
  ! **Description:**
  !
  !- Usage:
  !     One input data record is required
  !         READ (LIN, '(I1)') KPRINT
  !
  !- Arguments:
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
  !- Description:
  !     Driver for testing SLATEC subprograms
  !       Double precision SLAP subprograms
  !
  !***
  ! **References:**  Fong, Kirby W., Jefferson, Thomas H., Suyehiro,
  !                 Tokihiko, Walton, Lee, Guidelines to the SLATEC Common
  !                 Mathematical Library, March 21, 1989.
  !***
  ! **Routines called:**  DLAPQC, I1MACH, XERMAX, XSETF, XSETUN

  !* REVISION HISTORY  (YYMMDD)
  !   920401  DATE WRITTEN
  !   920511  Added complete declaration section.  (WRB)

  !     .. Local Scalars ..
  INTEGER ipass, kprint, lin, lun, nfail
  !* FIRST EXECUTABLE STATEMENT  TEST26
  lun = I1MACH(2)
  lin = I1MACH(1)
  nfail = 0
  !
  !     Read KPRINT parameter
  !
  CALL GET_ARGUMENT(kprint)
  CALL XSETUN(lun)
  IF ( kprint<=1 ) THEN
    CALL XSETF(0)
  ELSE
    CALL XSETF(1)
  END IF
  CALL XERMAX(1000)
  !
  !     Test SLAP (double precision)
  !
  CALL DLAPQC(lun,kprint,ipass)
  IF ( ipass==0 ) nfail = nfail + 1
  !
  !     Write PASS or FAIL message
  !
  IF ( nfail==0 ) THEN
    WRITE (lun,99001)
    99001 FORMAT (/' --------------TEST26 PASSED ALL TESTS----------------')
  ELSE
    WRITE (lun,99002) nfail
    99002 FORMAT (/' ************* WARNING -- ',I5,&
      ' TEST(S) FAILED IN PROGRAM TEST26 *************')
  END IF
  STOP
END PROGRAM TEST26
