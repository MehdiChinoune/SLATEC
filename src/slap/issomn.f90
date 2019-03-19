!** ISSOMN
INTEGER FUNCTION ISSOMN(N,B,X,Nelt,Ia,Ja,A,Isym,MSOLVE,Nsave,Itol,Tol,&
    Itmax,Iter,Err,Ierr,Iunit,R,Z,P,Ap,Emap,Dz,Csav,&
    Rwork,Iwork,Ak,Bnrm,Solnrm)
  IMPLICIT NONE
  !>
  !***
  !  Preconditioned Orthomin Stop Test.
  !            This routine calculates the stop test for the Orthomin
  !            iteration scheme.  It returns a non-zero if the error
  !            estimate (the type of which is determined by ITOL) is
  !            less than the user specified tolerance TOL.
  !***
  ! **Library:**   SLATEC (SLAP)
  !***
  ! **Category:**  D2A4, D2B4
  !***
  ! **Type:**      SINGLE PRECISION (ISSOMN-S, ISDOMN-D)
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
  !     REAL     B(N), X(N), A(NELT), TOL, ERR, R(N), Z(N)
  !     REAL     P(N,0:NSAVE), AP(N,0:NSAVE), EMAP(N,0:NSAVE)
  !     REAL     DZ(N), CSAV(NSAVE), RWORK(USER DEFINED), AK
  !     REAL     BNRM, SOLNRM
  !     EXTERNAL MSOLVE
  !
  !     IF( ISSOMN(N, B, X, NELT, IA, JA, A, ISYM, MSOLVE, NSAVE,
  !    $     ITOL, TOL, ITMAX, ITER, ERR, IERR, IUNIT, R, Z, P, AP,
  !    $     EMAP, DZ, CSAV, RWORK, IWORK, AK, BNRM, SOLNRM)
  !    $     .NE.0 ) THEN ITERATION CONVERGED
  !
  !- Arguments:
  ! N      :IN       Integer.
  !         Order of the matrix.
  ! B      :IN       Real B(N).
  !         Right-hand side vector.
  ! X      :IN       Real X(N).
  !         On input X is your initial guess for solution vector.
  !         On output X is the final approximate solution.
  ! NELT   :IN       Integer.
  !         Number of Non-Zeros stored in A.
  ! IA     :IN       Integer IA(NELT).
  ! JA     :IN       Integer JA(NELT).
  ! A      :IN       Real A(NELT).
  !         These arrays should hold the matrix A in either the SLAP
  !         Triad format or the SLAP Column format.  See "Description"
  !         in the SSDOMN or SSLUOM prologue.
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
  !         ISYM are defined as above.  RWORK is a real array that can
  !         be used to pass necessary preconditioning information and/or
  !         workspace to MSOLVE.  IWORK is an integer work array for
  !         the same purpose as RWORK.
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
  !             COMMON /SSLBLK/ SOLN( )
  !         If ITOL=11, iteration stops when the 2-norm of the difference
  !         between the iterative approximation and the user-supplied
  !         solution divided by the 2-norm of the user-supplied solution
  !         is less than TOL.  Note that this requires the user to set up
  !         the "COMMON /SSLBLK/ SOLN(LENGTH)" in the calling routine.
  !         The routine with this declaration should be loaded before the
  !         stop test so that the correct length is used by the loader.
  !         This procedure is not standard Fortran and may not work
  !         correctly on your system (although it has worked on every
  !         system the authors have tried).  If ITOL is not 11 then this
  !         common block is indeed standard Fortran.
  ! TOL    :IN       Real.
  !         Convergence criterion, as described above.
  ! ITMAX  :IN       Integer.
  !         Maximum number of iterations.
  ! ITER   :IN       Integer.
  !         Current iteration count.  (Must be zero on first call.)
  ! ERR    :OUT      Real.
  !         Error estimate of error in final approximate solution, as
  !         defined by ITOL.
  ! IERR   :OUT      Integer.
  !         Error flag.  IERR is set to 3 if ITOL is not one of the
  !         acceptable values, see above.
  ! IUNIT  :IN       Integer.
  !         Unit number on which to write the error at each iteration,
  !         if this is desired for monitoring convergence.  If unit
  !         number is 0, no writing will occur.
  ! R      :IN       Real R(N).
  !         The residual R = B-AX.
  ! Z      :WORK     Real Z(N).
  ! P      :IN       Real P(N,0:NSAVE).
  !         Workspace used to hold the conjugate direction vector(s).
  ! AP     :IN       Real AP(N,0:NSAVE).
  !         Workspace used to hold the matrix A times the P vector(s).
  ! EMAP   :IN       Real EMAP(N,0:NSAVE).
  !         Workspace used to hold M-inv times the AP vector(s).
  ! DZ     :WORK     Real DZ(N).
  !         Workspace.
  ! CSAV   :DUMMY    Real CSAV(NSAVE)
  !         Reserved for future use.
  ! RWORK  :WORK     Real RWORK(USER DEFINED).
  !         Real array that can be used for workspace in MSOLVE.
  ! IWORK  :WORK     Integer IWORK(USER DEFINED).
  !         Integer array that can be used for workspace in MSOLVE.
  ! AK     :IN       Real.
  !         Current iterate Orthomin iteration parameter.
  ! BNRM   :OUT      Real.
  !         Current solution B-norm, if ITOL = 1 or 2.
  ! SOLNRM :OUT      Real.
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
  !     unit IUNIT, if IUNIT .ne. 0.  Thus, the user must make sure that
  !     this logical unit is attached to a file or terminal before calling
  !     this routine with a non-zero value for IUNIT.  This routine does
  !     not check for the validity of a non-zero IUNIT unit number.
  !
  !***
  ! **See also:**  SOMN, SSDOMN, SSLUOM
  !***
  ! **Routines called:**  R1MACH, SNRM2
  !***
  ! COMMON BLOCKS    SSLBLK

  !* REVISION HISTORY  (YYMMDD)
  !   871119  DATE WRITTEN
  !   881213  Previous REVISION DATE
  !   890915  Made changes requested at July 1989 CML Meeting.  (MKS)
  !   890922  Numerous changes to prologue to make closer to SLATEC
  !           standard.  (FNF)
  !   890929  Numerous changes to reduce SP/DP differences.  (FNF)
  !   891003  Removed C***REFER TO line, per MKS.
  !   910411  Prologue converted to Version 4.0 format.  (BAB)
  !   910502  Removed MSOLVE from ROUTINES CALLED list.  (FNF)
  !   910506  Made subsidiary to SOMN.  (FNF)
  !   920407  COMMON BLOCK renamed SSLBLK.  (WRB)
  !   920511  Added complete declaration section.  (WRB)
  !   920930  Corrected to not print AK when ITER=0.  (FNF)
  !   921026  Changed 1.0E10 to R1MACH(2).  (FNF)
  !   921113  Corrected C***CATEGORY line.  (FNF)
  
  !     .. Scalar Arguments ..
  REAL Ak, Bnrm, Err, Solnrm, Tol
  INTEGER Ierr, Isym, Iter, Itmax, Itol, Iunit, N, Nelt, Nsave
  !     .. Array Arguments ..
  REAL A(Nelt), Ap(N,0:Nsave), B(N), Csav(Nsave), Dz(N), &
    Emap(N,0:Nsave), P(N,0:Nsave), R(N), Rwork(*), X(N), Z(N)
  INTEGER Ia(Nelt), Iwork(*), Ja(Nelt)
  !     .. Subroutine Arguments ..
  EXTERNAL MSOLVE
  !     .. Arrays in Common ..
  REAL SOLn(25)
  !     .. Local Scalars ..
  INTEGER i
  !     .. External Functions ..
  REAL R1MACH, SNRM2
  EXTERNAL R1MACH, SNRM2
  !     .. Common blocks ..
  COMMON /SSLBLK/ SOLn
  !* FIRST EXECUTABLE STATEMENT  ISSOMN
  ISSOMN = 0
  !
  IF ( Itol==1 ) THEN
    !         err = ||Residual||/||RightHandSide|| (2-Norms).
    IF ( Iter==0 ) Bnrm = SNRM2(N,B,1)
    Err = SNRM2(N,R,1)/Bnrm
  ELSEIF ( Itol==2 ) THEN
    !                  -1              -1
    !         err = ||M  Residual||/||M  RightHandSide|| (2-Norms).
    IF ( Iter==0 ) THEN
      CALL MSOLVE(N,B,Dz,Nelt,Ia,Ja,A,Isym,Rwork,Iwork)
      Bnrm = SNRM2(N,Dz,1)
    ENDIF
    Err = SNRM2(N,Z,1)/Bnrm
  ELSEIF ( Itol==11 ) THEN
    !         err = ||x-TrueSolution||/||TrueSolution|| (2-Norms).
    IF ( Iter==0 ) Solnrm = SNRM2(N,SOLn,1)
    DO i = 1, N
      Dz(i) = X(i) - SOLn(i)
    ENDDO
    Err = SNRM2(N,Dz,1)/Solnrm
  ELSE
    !
    !         If we get here ITOL is not one of the acceptable values.
    Err = R1MACH(2)
    Ierr = 3
  ENDIF
  !
  IF ( Iunit/=0 ) THEN
    IF ( Iter==0 ) THEN
      WRITE (Iunit,99001) Nsave, N, Itol
      99001 FORMAT (' Preconditioned Orthomin(',I3,') for ','N, ITOL = ',I5,I5,&
        /' ITER','   Error Estimate','            Alpha')
      WRITE (Iunit,99002) Iter, Err
    ELSE
      WRITE (Iunit,99002) Iter, Err, Ak
    ENDIF
  ENDIF
  IF ( Err<=Tol ) ISSOMN = 1
  !
  RETURN
  99002 FORMAT (1X,I4,1X,E16.7,1X,E16.7)
  !------------- LAST LINE OF ISSOMN FOLLOWS ----------------------------
END FUNCTION ISSOMN
