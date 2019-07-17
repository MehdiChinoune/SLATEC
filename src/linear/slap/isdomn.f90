!** ISDOMN
INTEGER PURE FUNCTION ISDOMN(N,Itol,Tol,R,Z,Dz,Bnrm,Solnrm)
  !> Preconditioned Orthomin Stop Test.
  !  This routine calculates the stop test for the Orthomin iteration scheme.
  !  It returns a non-zero if the error estimate (the type of which is determined
  !  by ITOL) is less than the user specified tolerance TOL.
  !***
  ! **Library:**   SLATEC (SLAP)
  !***
  ! **Category:**  D2A4, D2B4
  !***
  ! **Type:**      DOUBLE PRECISION (ISSOMN-S, ISDOMN-D)
  !***
  ! **Keywords:**  ITERATIVE PRECONDITION, NON-SYMMETRIC LINEAR SYSTEM,
  !             ORTHOMIN, SLAP, SPARSE, STOP TEST
  !***
  ! **Author:**  Greenbaum, Anne, (Courant Institute)
  !           Seager, Mark K., (LLNL)
  !             Lawrence Livermore National Laboratory
  !             PO BOX 808, L-60
  !             Livermore, CA 94550 (510) 423-3141
  !             seager@llnl.gov
  !***
  ! **Description:**
  !
  !- Usage:
  !     INTEGER  N, NELT, IA(NELT), JA(NELT), ISYM, NSAVE, ITOL, ITMAX
  !     INTEGER  ITER, IERR, IUNIT, IWORK(USER DEFINED)
  !     DOUBLE PRECISION B(N), X(N), A(NELT), TOL, ERR, R(N), Z(N)
  !     DOUBLE PRECISION P(N,0:NSAVE), AP(N,0:NSAVE), EMAP(N,0:NSAVE)
  !     DOUBLE PRECISION DZ(N), CSAV(NSAVE), RWORK(USER DEFINED), AK
  !     DOUBLE PRECISION BNRM, SOLNRM
  !     EXTERNAL MSOLVE
  !
  !     IF( ISDOMN(N, B, X, NELT, IA, JA, A, ISYM, MSOLVE, NSAVE,
  !    $     ITOL, TOL, ITMAX, ITER, ERR, IERR, IUNIT, R, Z, P, AP,
  !    $     EMAP, DZ, CSAV, RWORK, IWORK, AK, BNRM, SOLNRM)
  !    $     /=0 ) THEN ITERATION CONVERGED
  !
  !- Arguments:
  ! N      :IN       Integer.
  !         Order of the matrix.
  ! B      :IN       Double Precision B(N).
  !         Right-hand side vector.
  ! X      :IN       Double Precision X(N).
  !         On input X is your initial guess for solution vector.
  !         On output X is the final approximate solution.
  ! NELT   :IN       Integer.
  !         Number of Non-Zeros stored in A.
  ! IA     :IN       Integer IA(NELT).
  ! JA     :IN       Integer JA(NELT).
  ! A      :IN       Double Precision A(NELT).
  !         These arrays should hold the matrix A in either the SLAP
  !         Triad format or the SLAP Column format.  See "Description"
  !         in the DSDOMN or DSLUOM prologue.
  ! ISYM   :IN       Integer.
  !         Flag to indicate symmetric storage format.
  !         If ISYM=0, all non-zero entries of the matrix are stored.
  !         If ISYM=1, the matrix is symmetric, and only the upper
  !         or lower triangle of the matrix is stored.
  ! MSOLVE :EXT      External.
  !         Name of a routine which solves a linear system MZ = R for
  !         Z given R with the preconditioning matrix M (M is supplied via
  !         RWORK and IWORK arrays).  The name of the MSOLVE routine must
  !         be declared external in the calling program.  The calling
  !         sequence to MSOLVE is:
  !             CALL MSOLVE(N, R, Z, NELT, IA, JA, A, ISYM, RWORK, IWORK)
  !         Where N is the number of unknowns, R is the right-hand side
  !         vector and Z is the solution upon return.  NELT, IA, JA, A and
  !         ISYM are defined as above.  RWORK is a double precision array
  !         that can be used to pass necessary preconditioning information
  !         and/or workspace to MSOLVE.  IWORK is an integer work array
  !         for the same purpose as RWORK.
  ! NSAVE  :IN       Integer.
  !         Number of direction vectors to save and orthogonalize against.
  ! ITOL   :IN       Integer.
  !         Flag to indicate type of convergence criterion.
  !         If ITOL=1, iteration stops when the 2-norm of the residual
  !         divided by the 2-norm of the right-hand side is less than TOL.
  !         If ITOL=2, iteration stops when the 2-norm of M-inv times the
  !         residual divided by the 2-norm of M-inv times the right hand
  !         side is less than TOL, where M-inv is the inverse of the
  !         diagonal of A.
  !         ITOL=11 is often useful for checking and comparing different
  !         routines.  For this case, the user must supply the "exact"
  !         solution or a very accurate approximation (one with an error
  !         much less than TOL) through a common block,
  !             COMMON /DSLBLK/ SOLN( )
  !         If ITOL=11, iteration stops when the 2-norm of the difference
  !         between the iterative approximation and the user-supplied
  !         solution divided by the 2-norm of the user-supplied solution
  !         is less than TOL.  Note that this requires the user to set up
  !         the "COMMON /DSLBLK/ SOLN(LENGTH)" in the calling routine.
  !         The routine with this declaration should be loaded before the
  !         stop test so that the correct length is used by the loader.
  !         This procedure is not standard Fortran and may not work
  !         correctly on your system (although it has worked on every
  !         system the authors have tried).  If ITOL is not 11 then this
  !         common block is indeed standard Fortran.
  ! TOL    :IN       Double Precision.
  !         Convergence criterion, as described above.
  ! ITMAX  :IN       Integer.
  !         Maximum number of iterations.
  ! ITER   :IN       Integer.
  !         Current iteration count.  (Must be zero on first call.)
  ! ERR    :OUT      Double Precision.
  !         Error estimate of error in final approximate solution, as
  !         defined by ITOL.
  ! IERR   :OUT      Integer.
  !         Error flag.  IERR is set to 3 if ITOL is not one of the
  !         acceptable values, see above.
  ! IUNIT  :IN       Integer.
  !         Unit number on which to write the error at each iteration,
  !         if this is desired for monitoring convergence.  If unit
  !         number is 0, no writing will occur.
  ! R      :IN       Double Precision R(N).
  !         The residual R = B-AX.
  ! Z      :WORK     Double Precision Z(N).
  ! P      :IN       Double Precision P(N,0:NSAVE).
  !         Workspace used to hold the conjugate direction vector(s).
  ! AP     :IN       Double Precision AP(N,0:NSAVE).
  !         Workspace used to hold the matrix A times the P vector(s).
  ! EMAP   :IN       Double Precision EMAP(N,0:NSAVE).
  !         Workspace used to hold M-inv times the AP vector(s).
  ! DZ     :WORK     Double Precision DZ(N).
  !         Workspace.
  ! CSAV   :DUMMY    Double Precision CSAV(NSAVE)
  !         Reserved for future use.
  ! RWORK  :WORK     Double Precision RWORK(USER DEFINED).
  !         Double Precision array that can be used for workspace in
  !         MSOLVE.
  ! IWORK  :WORK     Integer IWORK(USER DEFINED).
  !         Integer array that can be used for workspace in MSOLVE.
  ! AK     :IN       Double Precision.
  !         Current iterate Orthomin iteration parameter.
  ! BNRM   :OUT      Double Precision.
  !         Current solution B-norm, if ITOL = 1 or 2.
  ! SOLNRM :OUT      Double Precision.
  !         True solution norm, if ITOL = 11.
  !
  !- Function Return Values:
  !       0 : Error estimate (determined by ITOL) is *NOT* less than the
  !           specified tolerance, TOL.  The iteration must continue.
  !       1 : Error estimate (determined by ITOL) is less than the
  !           specified tolerance, TOL.  The iteration can be considered
  !           complete.
  !
  !- Cautions:
  !     This routine will attempt to write to the Fortran logical output
  !     unit IUNIT, if IUNIT /= 0.  Thus, the user must make sure that
  !     this logical unit is attached to a file or terminal before calling
  !     this routine with a non-zero value for IUNIT.  This routine does
  !     not check for the validity of a non-zero IUNIT unit number.
  !
  !***
  ! **See also:**  DOMN, DSDOMN, DSLUOM
  !***
  ! **Routines called:**  D1MACH, DNRM2
  !***
  ! COMMON BLOCKS    DSLBLK

  !* REVISION HISTORY  (YYMMDD)
  !   890404  DATE WRITTEN
  !   890404  Previous REVISION DATE
  !   890915  Made changes requested at July 1989 CML Meeting.  (MKS)
  !   890922  Numerous changes to prologue to make closer to SLATEC standard.  (FNF)
  !   890929  Numerous changes to reduce SP/DP differences.  (FNF)
  !   891003  Removed C***REFER TO line, per MKS.
  !   910411  Prologue converted to Version 4.0 format.  (BAB)
  !   910502  Removed MSOLVE from ROUTINES CALLED list.  (FNF)
  !   910506  Made subsidiary to DOMN.  (FNF)
  !   920407  COMMON BLOCK renamed DSLBLK.  (WRB)
  !   920511  Added complete declaration section.  (WRB)
  !   920930  Corrected to not print AK when ITER=0.  (FNF)
  !   921026  Changed 1.0E10 to D1MACH(2) and corrected D to E in
  !           output format.  (FNF)
  !   921113  Corrected C***CATEGORY line.  (FNF)
  USE service, ONLY : D1MACH
  !     .. Scalar Arguments ..
  INTEGER, INTENT(IN) :: Itol, N
  REAL(DP), INTENT(IN) :: Bnrm, Solnrm, Tol
  !     .. Array Arguments ..
  REAL(DP), INTENT(IN) :: Dz(N), R(N), Z(N)
  !     .. Local Scalars ..
  INTEGER :: ierr
  REAL(DP) :: eror
  !* FIRST EXECUTABLE STATEMENT  ISDOMN
  ISDOMN = 0
  !
  IF( Itol==1 ) THEN
    !         err = ||Residual||/||RightHandSide|| (2-Norms).
    eror = NORM2(R)/Bnrm
  ELSEIF( Itol==2 ) THEN
    !                  -1              -1
    !         err = ||M  Residual||/||M  RightHandSide|| (2-Norms).
    eror = NORM2(Z)/Bnrm
  ELSEIF( Itol==11 ) THEN
    !         err = ||x-TrueSolution||/||TrueSolution|| (2-Norms).
    eror = NORM2(Dz)/Solnrm
  ELSE
    !         If we get here ITOL is not one of the acceptable values.
    eror = D1MACH(2)
    ierr = 3
    ERROR STOP "Itol is not one of the acceptable values {1,2,11}"
  END IF
  !
  IF( eror<=Tol ) ISDOMN = 1
  !
  RETURN
  !------------- LAST LINE OF ISDOMN FOLLOWS ----------------------------
END FUNCTION ISDOMN