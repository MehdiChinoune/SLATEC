!*==DLAPQC.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK DLAPQC
SUBROUTINE DLAPQC(Lun,Kprint,Ipass)
  IMPLICIT NONE
  !*--DLAPQC5
  !***BEGIN PROLOGUE  DLAPQC
  !***PURPOSE  Quick check for testing Sparse Linear Algebra Package
  !            (SLAP) Version 2.0.2.
  !***LIBRARY   SLATEC (SLAP)
  !***CATEGORY  D2A4, D2B4
  !***TYPE      DOUBLE PRECISION (SLAPQC-S, DLAPQC-D)
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
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  D1MACH, DCPPLT, DFILL, DRMGEN, DS2Y, DSDBCG, DSDCG,
  !                    DSDCGN, DSDCGS, DSDGMR, DSDOMN, DSGS, DSICCG,
  !                    DSILUR, DSJAC, DSLUBC, DSLUCN, DSLUCS, DSLUGM,
  !                    DSLUOM, DUTERR, XERMAX, XSETF, XSETUN
  !***COMMON BLOCKS    DSLBLK
  !***REVISION HISTORY  (YYMMDD)
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
  !***END PROLOGUE  DLAPQC
  !
  !     The problem size, MAXN, should be large enough that the
  !     iterative methods do 10-15 iterations, just to be sure that
  !     the truncated methods run to the end of their ropes and enter
  !     their error recovery mode.  Thus, for a more thorough test change
  !     the following PARAMETER statement to:
  !     PARAMETER (MAXN=69, MXNELT=5000, MAXIW=5000, MAXRW=5000)
  !
  !     .. Parameters ..
  INTEGER MAXN , MXNELT , MAXIW , MAXRW
  PARAMETER (MAXN=25,MXNELT=500,MAXIW=1000,MAXRW=1000)
  !     .. Scalar Arguments ..
  INTEGER Ipass , Kprint , Lun
  !     .. Arrays in Common ..
  REAL(8) :: SOLn(MAXN)
  !     .. Local Scalars ..
  REAL(8) :: dens , err , factor , tol
  INTEGER ierr , isym , iter , itmax , itol , itolgm , iunit , k , kase , &
    leniw , lenw , n , nelt , neltmx , nfail , nmax , nsave
  !     .. Local Arrays ..
  REAL(8) :: a(MXNELT) , f(MAXN) , rwork(MAXRW) , xiter(MAXN)
  INTEGER ia(MXNELT) , iwork(MAXIW) , ja(MXNELT)
  !     .. External Functions ..
  REAL(8) :: D1MACH
  EXTERNAL D1MACH
  !     .. External Subroutines ..
  EXTERNAL DCPPLT , DFILL , DRMGEN , DS2Y , DSDBCG , DSDCG , DSDCGN , &
    DSDCGS , DSDGMR , DSDOMN , DSGS , DSICCG , DSILUR , DSJAC , &
    DSLUBC , DSLUCN , DSLUCS , DSLUGM , DSLUOM , DUTERR
  !     .. Intrinsic Functions ..
  INTRINSIC MAX , REAL
  !     .. Common blocks ..
  COMMON /DSLBLK/ SOLn
  !
  !     The following lines are for the braindamaged Sun FPE handler.
  !
  !$$$      integer oldmode, fpmode
  !***FIRST EXECUTABLE STATEMENT  DLAPQC
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
  factor = 1.2D0
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
  tol = MAX(1.0D3*D1MACH(3),1.0D-6)
  nfail = 0
  !
  !     Test routines using various convergence criteria.
  !
  DO kase = 1 , 3
    IF ( kase==1.OR.kase==2 ) itol = kase
    IF ( kase==3 ) itol = 11
    !
    !         Test routines using nonsymmetric (ISYM=0) and symmetric
    !         storage (ISYM=1).  For ISYM=0 a really non-symmetric matrix
    !         is generated.  The amount of non-symmetry is controlled by
    !         user.
    !
    DO isym = 0 , 1
      IF ( Kprint>=2 ) WRITE (Lun,99001) n , kase , isym
      99001     FORMAT ('1'/' Running tests with  N =',I3,',  KASE =',I2,',  ISYM =',&
        I2)
      !
      !         Set up a random matrix.
      !
      CALL DRMGEN(neltmx,factor,ierr,n,nelt,isym,ia,ja,a,f,SOLn,rwork,iwork,&
        iwork(n+1))
      IF ( ierr/=0 ) THEN
        WRITE (Lun,99002) ierr
        !
        99002       FORMAT (/1X,'DLAPQC -- Fatal error ',I1,' generating ',&
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
        WRITE (Lun,99003) n , nelt , dens
        99003       FORMAT (/'                * RANDOM Matrix of size',I5,&
          '*'/'                ','Number of non-zeros & Density = ',&
          I5,1P,D16.7)
        WRITE (Lun,99004) tol
        99004       FORMAT ('                Error tolerance = ',1P,D16.7)
      ENDIF
      !
      !         Convert to the SLAP-Column format and
      !         write out matrix in SLAP-Column format, if desired.
      !
      CALL DS2Y(n,nelt,ia,ja,a,isym)
      IF ( Kprint>=4 ) THEN
        WRITE (Lun,99005) (k,ia(k),ja(k),a(k),k=1,nelt)
        99005       FORMAT (/'  ***** SLAP Column Matrix *****'/' Indx   ia   ja     a'/&
          (1X,I4,1X,I4,1X,I4,1X,1P,D16.7))
        CALL DCPPLT(n,nelt,ia,ja,a,isym,Lun)
      ENDIF
      !
      !**********************************************************************
      !                    BEGINNING OF SLAP QUICK TESTS
      !**********************************************************************
      !
      !         * * * * * *   DSJAC   * * * * * *
      !
      IF ( Kprint>=3 ) WRITE (Lun,99008) 'DSJAC ' , itol , isym
      CALL DFILL(n,xiter,0.0D0)
      !
      CALL DSJAC(n,f,xiter,nelt,ia,ja,a,isym,itol,tol,2*itmax,iter,err,ierr,&
        iunit,rwork,lenw,iwork,leniw)
      !
      CALL DUTERR('DSJAC ',ierr,Kprint,nfail,Lun,iter,err)
      !
      !         * * * * *  DSGS  * * * * *
      !
      IF ( Kprint>=3 ) WRITE (Lun,99008) 'DSGS  ' , itol , isym
      CALL DFILL(n,xiter,0.0D0)
      !
      CALL DSGS(n,f,xiter,nelt,ia,ja,a,isym,itol,tol,itmax,iter,err,ierr,&
        iunit,rwork,lenw,iwork,leniw)
      !
      CALL DUTERR('DSGS  ',ierr,Kprint,nfail,Lun,iter,err)
      !
      !         * * * * * *   DSILUR   * * * * * *
      !
      IF ( Kprint>=3 ) WRITE (Lun,99008) 'DSILUR' , itol , isym
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
        IF ( Kprint>=3 ) WRITE (Lun,99008) 'DSDCG' , itol , isym
        CALL DFILL(n,xiter,0.0D0)
        !
        CALL DSDCG(n,f,xiter,nelt,ia,ja,a,isym,itol,tol,itmax,iter,err,ierr,&
          iunit,rwork,lenw,iwork,leniw)
        !
        CALL DUTERR('DSDCG ',ierr,Kprint,nfail,Lun,iter,err)
      ENDIF
      !
      !         * * * * * *    DSICCG    * * * * * *
      !
      IF ( isym==1 ) THEN
        IF ( Kprint>=3 ) WRITE (Lun,99008) 'DSICCG' , itol , isym
        CALL DFILL(n,xiter,0.0D0)
        !
        CALL DSICCG(n,f,xiter,nelt,ia,ja,a,isym,itol,tol,itmax,iter,err,&
          ierr,iunit,rwork,lenw,iwork,leniw)
        !
        CALL DUTERR('DSICCG',ierr,Kprint,nfail,Lun,iter,err)
      ENDIF
      !
      !         * * * * * *    DSDCGN   * * * * * *
      !
      IF ( Kprint>=3 ) WRITE (Lun,99008) 'DSDCGN' , itol , isym
      CALL DFILL(n,xiter,0.0D0)
      !
      CALL DSDCGN(n,f,xiter,nelt,ia,ja,a,isym,itol,tol,itmax,iter,err,ierr,&
        iunit,rwork,lenw,iwork,leniw)
      !
      CALL DUTERR('DSDCGN',ierr,Kprint,nfail,Lun,iter,err)
      !
      !         * * * * * *   DSLUCN   * * * * * *
      !
      IF ( Kprint>=3 ) WRITE (Lun,99008) 'DSLUCN' , itol , isym
      CALL DFILL(n,xiter,0.0D0)
      !
      CALL DSLUCN(n,f,xiter,nelt,ia,ja,a,isym,itol,tol,itmax,iter,err,ierr,&
        iunit,rwork,lenw,iwork,leniw)
      !
      CALL DUTERR('DSLUCN',ierr,Kprint,nfail,Lun,iter,err)
      !
      !         * * * * * *    DSDBCG   * * * * * *
      !
      IF ( Kprint>=3 ) WRITE (Lun,99008) 'DSDBCG' , itol , isym
      CALL DFILL(n,xiter,0.0D0)
      !
      CALL DSDBCG(n,f,xiter,nelt,ia,ja,a,isym,itol,tol,itmax,iter,err,ierr,&
        iunit,rwork,lenw,iwork,leniw)
      !
      CALL DUTERR('DSDBCG',ierr,Kprint,nfail,Lun,iter,err)
      !
      !         * * * * * *   DSLUBC   * * * * * *
      !
      IF ( Kprint>=3 ) WRITE (Lun,99008) 'DSLUBC' , itol , isym
      CALL DFILL(n,xiter,0.0D0)
      !
      CALL DSLUBC(n,f,xiter,nelt,ia,ja,a,isym,itol,tol,itmax,iter,err,ierr,&
        iunit,rwork,lenw,iwork,leniw)
      !
      CALL DUTERR('DSLUBC',ierr,Kprint,nfail,Lun,iter,err)
      !
      !         * * * * * *    DSDCGS   * * * * * *
      !
      IF ( Kprint>=3 ) WRITE (Lun,99008) 'DSDCGS' , itol , isym
      CALL DFILL(n,xiter,0.0D0)
      !
      CALL DSDCGS(n,f,xiter,nelt,ia,ja,a,isym,itol,tol,itmax,iter,err,ierr,&
        iunit,rwork,lenw,iwork,leniw)
      !
      CALL DUTERR('DSDCGS',ierr,Kprint,nfail,Lun,iter,err)
      !
      !         * * * * * *   DSLUCS   * * * * * *
      !
      IF ( Kprint>=3 ) WRITE (Lun,99008) 'DSLUCS' , itol , isym
      CALL DFILL(n,xiter,0.0D0)
      !
      CALL DSLUCS(n,f,xiter,nelt,ia,ja,a,isym,itol,tol,itmax,iter,err,ierr,&
        iunit,rwork,lenw,iwork,leniw)
      !
      CALL DUTERR('DSLUCS',ierr,Kprint,nfail,Lun,iter,err)
      !
      !         * * * * * *    DSDOMN   * * * * * *
      !
      !VD$ NOVECTOR
      DO nsave = 0 , 3
        IF ( Kprint>=3 ) WRITE (Lun,99009) 'DSDOMN' , itol , isym , nsave
        CALL DFILL(n,xiter,0.0D0)
        !
        CALL DSDOMN(n,f,xiter,nelt,ia,ja,a,isym,nsave,itol,tol,itmax,iter,&
          err,ierr,iunit,rwork,lenw,iwork,leniw)
        !
        CALL DUTERR('DSDOMN',ierr,Kprint,nfail,Lun,iter,err)
      ENDDO
      !
      !         * * * * * *   DSLUOM   * * * * * *
      !
      !VD$ NOVECTOR
      DO nsave = 0 , 3
        IF ( Kprint>=3 ) WRITE (Lun,99009) 'DSLUOM' , itol , isym , nsave
        CALL DFILL(n,xiter,0.0D0)
        !
        CALL DSLUOM(n,f,xiter,nelt,ia,ja,a,isym,nsave,itol,tol,itmax,iter,&
          err,ierr,iunit,rwork,lenw,iwork,leniw)
        !
        CALL DUTERR('DSLUOM',ierr,Kprint,nfail,Lun,iter,err)
      ENDDO
      !
      !         * * * * * *   DSDGMR   * * * * * *
      !
      !VD$ NOVECTOR
      DO nsave = 5 , 12
        IF ( Kprint>=3 ) WRITE (Lun,99009) 'DSDGMR' , itol , isym , nsave
        CALL DFILL(n,xiter,0.0D0)
        itolgm = 0
        !
        CALL DSDGMR(n,f,xiter,nelt,ia,ja,a,isym,nsave,itolgm,tol,itmax,iter,&
          err,ierr,iunit,rwork,lenw,iwork,leniw)
        !
        CALL DUTERR('DSDGMR',ierr,Kprint,nfail,Lun,iter,err)
      ENDDO
      !
      !         * * * * * *   DSLUGM   * * * * * *
      !
      !VD$ NOVECTOR
      DO nsave = 5 , 12
        IF ( Kprint>=3 ) WRITE (Lun,99009) 'DSLUGM' , itol , isym , nsave
        CALL DFILL(n,xiter,0.0D0)
        !
        CALL DSLUGM(n,f,xiter,nelt,ia,ja,a,isym,nsave,itol,tol,itmax,iter,&
          err,ierr,iunit,rwork,lenw,iwork,leniw)
        !
        CALL DUTERR('DSLUGM',ierr,Kprint,nfail,Lun,iter,err)
      ENDDO
    ENDDO
  ENDDO
  !
  IF ( nfail==0 ) THEN
    Ipass = 1
    IF ( Kprint>=2 ) WRITE (Lun,99006)
    99006   FORMAT ('--------- All double precision SLAP tests passed ','---------')
  ELSE
    Ipass = 0
    IF ( Kprint>=2 ) WRITE (Lun,99007) nfail
    99007   FORMAT ('*********',I3,' double precision SLAP tests failed ',&
      '*********')
  ENDIF
  !
  RETURN
  99008 FORMAT (/1X,A6,' : ITOL = ',I2,'   ISYM = ',I1)
  99009 FORMAT (/1X,A6,' : ITOL = ',I2,'   ISYM = ',I1,' NSAVE = ',I2)
END SUBROUTINE DLAPQC
