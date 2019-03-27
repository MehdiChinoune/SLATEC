!** ISDCGS
INTEGER FUNCTION ISDCGS(N,B,X,Nelt,Ia,Ja,A,Isym,MATVEC,MSOLVE,Itol,Tol,&
    Itmax,Iter,Err,Ierr,Iunit,R,R0,P,Q,U,V1,V2,Rwork,Iwork,Ak,Bk,Bnrm,Solnrm)
  IMPLICIT NONE
  !>
  !***
  !  Preconditioned BiConjugate Gradient Squared Stop Test.
  !            This routine calculates the stop test for the BiConjugate
  !            Gradient Squared iteration scheme.  It returns a non-zero
  !            if the error estimate (the type of which is determined by
  !            ITOL) is less than the user specified tolerance TOL.
  !***
  ! **Library:**   SLATEC (SLAP)
  !***
  ! **Category:**  D2A4, D2B4
  !***
  ! **Type:**      DOUBLE PRECISION (ISSCGS-S, ISDCGS-D)
  !***
  ! **Keywords:**  ITERATIVE PRECONDITION, NON-SYMMETRIC LINEAR SYSTEM, SLAP,
  !             SPARSE, STOP TEST
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
  !     INTEGER  N, NELT, IA(NELT), JA(NELT), ISYM, ITOL, ITMAX, ITER
  !     INTEGER  IERR, IUNIT, IWORK(USER DEFINED)
  !     DOUBLE PRECISION B(N), X(N), A(N), TOL, ERR, R(N), R0(N), P(N)
  !     DOUBLE PRECISION Q(N), U(N), V1(N), V2(N)
  !     DOUBLE PRECISION RWORK(USER DEFINED), AK, BK, BNRM, SOLNRM
  !     EXTERNAL MATVEC, MSOLVE
  !
  !     IF( ISDCGS(N, B, X, NELT, IA, JA, A, ISYM, MATVEC, MSOLVE, ITOL,
  !    $     TOL, ITMAX, ITER, ERR, IERR, IUNIT, R, R0, P, Q, U, V1,
  !    $     V2, RWORK, IWORK, AK, BK, BNRM, SOLNRM) .NE. 0 )
  !    $     THEN ITERATION DONE
  !
  !- Arguments:
  ! N      :IN       Integer
  !         Order of the Matrix.
  ! B      :IN       Double Precision B(N).
  !         Right-hand side vector.
  ! X      :INOUT    Double Precision X(N).
  !         On input X is your initial guess for solution vector.
  !         On output X is the final approximate solution.
  ! NELT   :IN       Integer.
  !         Number of Non-Zeros stored in A.
  ! IA     :IN       Integer IA(NELT).
  ! JA     :IN       Integer JA(NELT).
  ! A      :IN       Double Precision A(NELT).
  !         These arrays contain the matrix data structure for A.
  !         It could take any form.  See "Description" in SLAP routine
  !         DCGS for more details.
  ! ISYM   :IN       Integer.
  !         Flag to indicate symmetric storage format.
  !         If ISYM=0, all non-zero entries of the matrix are stored.
  !         If ISYM=1, the matrix is symmetric, and only the upper
  !         or lower triangle of the matrix is stored.
  ! MATVEC :EXT      External.
  !         Name of a routine which  performs the matrix vector multiply
  !         operation  Y = A*X  given A and X.  The  name of  the MATVEC
  !         routine must  be declared external  in the  calling program.
  !         The calling sequence of MATVEC is:
  !             CALL MATVEC( N, X, Y, NELT, IA, JA, A, ISYM )
  !         Where N is the number of unknowns, Y is the product A*X upon
  !         return,  X is an input  vector.  NELT, IA,  JA, A, and  ISYM
  !         define the SLAP matrix data structure.
  ! MSOLVE :EXT      External.
  !         Name of a routine which solves a linear system MZ = R  for Z
  !         given R with the preconditioning matrix M (M is supplied via
  !         RWORK  and IWORK arrays).   The name  of  the MSOLVE routine
  !         must be declared  external  in the  calling   program.   The
  !         calling sequence of MSOLVE is:
  !             CALL MSOLVE(N, R, Z, NELT, IA, JA, A, ISYM, RWORK, IWORK)
  !         Where N is the number of unknowns, R is  the right-hand side
  !         vector, and Z is the solution upon return.  NELT, IA, JA, A,
  !         and ISYM define the SLAP matrix data structure.
  !         RWORK is a double precision array that can be used to pass
  !         necessary preconditioning information and/or workspace to
  !         MSOLVE.
  !         IWORK is an integer work array for the same purpose as RWORK.
  ! ITOL   :IN       Integer.
  !         Flag to indicate type of convergence criterion.
  !         If ITOL=1, iteration stops when the 2-norm of the residual
  !         divided by the 2-norm of the right-hand side is less than TOL.
  !         This routine must calculate the residual from R = A*X - B.
  !         This is unnatural and hence expensive for this type of iter-
  !         ative method.  ITOL=2 is *STRONGLY* recommended.
  !         If ITOL=2, iteration stops when the 2-norm of M-inv times the
  !         residual divided by the 2-norm of M-inv times the right hand
  !         side is less than TOL, where M-inv time a vector is the pre-
  !         conditioning step.  This is the *NATURAL* stopping for this
  !         iterative method and is *STRONGLY* recommended.
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
  !         ITMAX iterations.
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
  !         The residual r = b - Ax.
  ! R0     :WORK     Double Precision R0(N).
  ! P      :DUMMY    Double Precision P(N).
  ! Q      :DUMMY    Double Precision Q(N).
  ! U      :DUMMY    Double Precision U(N).
  ! V1     :DUMMY    Double Precision V1(N).
  !         Double Precision arrays used for workspace.
  ! V2     :WORK     Double Precision V2(N).
  !         If ITOL.eq.1 then V2 is used to hold A * X - B on every call.
  !         If ITOL.eq.2 then V2 is used to hold M-inv * B on the first
  !         call.
  !         If ITOL.eq.11 then V2 is used to X - SOLN.
  ! RWORK  :WORK     Double Precision RWORK(USER DEFINED).
  !         Double Precision array that can be used for workspace in
  !         MSOLVE.
  ! IWORK  :WORK     Integer IWORK(USER DEFINED).
  !         Integer array that can be used for workspace in MSOLVE.
  ! AK     :IN       Double Precision.
  !         Current iterate BiConjugate Gradient iteration parameter.
  ! BK     :IN       Double Precision.
  !         Current iterate BiConjugate Gradient iteration parameter.
  ! BNRM   :INOUT    Double Precision.
  !         Norm of the right hand side.  Type of norm depends on ITOL.
  !         Calculated only on the first call.
  ! SOLNRM :INOUT    Double Precision.
  !         2-Norm of the true solution, SOLN.  Only computed and used
  !         if ITOL = 11.
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
  ! **See also:**  DCGS
  !***
  ! **Routines called:**  D1MACH, DNRM2
  !***
  ! COMMON BLOCKS    DSLBLK

  !* REVISION HISTORY  (YYMMDD)
  !   890404  DATE WRITTEN
  !   890404  Previous REVISION DATE
  !   890915  Made changes requested at July 1989 CML Meeting.  (MKS)
  !   890922  Numerous changes to prologue to make closer to SLATEC
  !           standard.  (FNF)
  !   890929  Numerous changes to reduce SP/DP differences.  (FNF)
  !   891003  Removed C***REFER TO line, per MKS.
  !   910411  Prologue converted to Version 4.0 format.  (BAB)
  !   910502  Removed MATVEC and MSOLVE from ROUTINES CALLED list.  (FNF)
  !   910506  Made subsidiary to DCGS.  (FNF)
  !   920407  COMMON BLOCK renamed DSLBLK.  (WRB)
  !   920511  Added complete declaration section.  (WRB)
  !   920930  Corrected to not print AK,BK when ITER=0.  (FNF)
  !   921026  Changed 1.0E10 to D1MACH(2) and corrected D to E in
  !           output format.  (FNF)
  !   921113  Corrected C***CATEGORY line.  (FNF)

  !     .. Scalar Arguments ..
  REAL(8) :: Ak, Bk, Bnrm, Err, Solnrm, Tol
  INTEGER Ierr, Isym, Iter, Itmax, Itol, Iunit, N, Nelt
  !     .. Array Arguments ..
  REAL(8) :: A(Nelt), B(N), P(N), Q(N), R(N), R0(N), Rwork(*), &
    U(N), V1(N), V2(N), X(N)
  INTEGER Ia(Nelt), Iwork(*), Ja(Nelt)
  !     .. Subroutine Arguments ..
  EXTERNAL :: MATVEC, MSOLVE
  !     .. Arrays in Common ..
  REAL(8) :: SOLn(25)
  !     .. Local Scalars ..
  INTEGER i
  !     .. External Functions ..
  REAL(8), EXTERNAL :: D1MACH, DNRM2
  !     .. Common blocks ..
  COMMON /DSLBLK/ SOLn
  !* FIRST EXECUTABLE STATEMENT  ISDCGS
  ISDCGS = 0
  !
  IF ( Itol==1 ) THEN
    !         err = ||Residual||/||RightHandSide|| (2-Norms).
    IF ( Iter==0 ) Bnrm = DNRM2(N,B,1)
    CALL MATVEC(N,X,V2,Nelt,Ia,Ja,A,Isym)
    DO i = 1, N
      V2(i) = V2(i) - B(i)
    ENDDO
    Err = DNRM2(N,V2,1)/Bnrm
  ELSEIF ( Itol==2 ) THEN
    !                  -1              -1
    !         err = ||M  Residual||/||M  RightHandSide|| (2-Norms).
    IF ( Iter==0 ) THEN
      CALL MSOLVE(N,B,V2,Nelt,Ia,Ja,A,Isym,Rwork,Iwork)
      Bnrm = DNRM2(N,V2,1)
    ENDIF
    Err = DNRM2(N,R,1)/Bnrm
  ELSEIF ( Itol==11 ) THEN
    !         err = ||x-TrueSolution||/||TrueSolution|| (2-Norms).
    IF ( Iter==0 ) Solnrm = DNRM2(N,SOLn,1)
    DO i = 1, N
      V2(i) = X(i) - SOLn(i)
    ENDDO
    Err = DNRM2(N,V2,1)/Solnrm
  ELSE
    !
    !         If we get here ITOL is not one of the acceptable values.
    Err = D1MACH(2)
    Ierr = 3
  ENDIF
  !
  !         Print the error and Coefficients AK, BK on each step,
  !         if desired.
  IF ( Iunit/=0 ) THEN
    IF ( Iter==0 ) THEN
      WRITE (Iunit,99001) N, Itol
      99001 FORMAT (' Preconditioned BiConjugate Gradient Squared for ',&
        'N, ITOL = ',I5,I5,/' ITER','   Error Estimate',&
        '            Alpha','             Beta')
      WRITE (Iunit,99002) Iter, Err
    ELSE
      WRITE (Iunit,99002) Iter, Err, Ak, Bk
    ENDIF
  ENDIF
  IF ( Err<=Tol ) ISDCGS = 1
  !
  RETURN
  99002 FORMAT (1X,I4,1X,D16.7,1X,D16.7,1X,D16.7)
  !------------- LAST LINE OF ISDCGS FOLLOWS ----------------------------
END FUNCTION ISDCGS
