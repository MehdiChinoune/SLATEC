!** DSDOMN
PURE SUBROUTINE DSDOMN(N,B,X,Nelt,Ia,Ja,A,Isym,Nsave,Itol,Tol,Itmax,Iter,Err,&
    Ierr,Rwork,Lenw,Iwork,Leniw)
  !> Diagonally Scaled Orthomin Sparse Iterative Ax=b Solver.
  !  Routine to solve a general linear system  Ax = b  using
  !  the Orthomin method with diagonal scaling.
  !***
  ! **Library:**   SLATEC (SLAP)
  !***
  ! **Category:**  D2A4, D2B4
  !***
  ! **Type:**      DOUBLE PRECISION (SSDOMN-S, DSDOMN-D)
  !***
  ! **Keywords:**  ITERATIVE PRECONDITION, NON-SYMMETRIC LINEAR SYSTEM SOLVE,
  !             SLAP, SPARSE
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
  !     INTEGER N, NELT, IA(NELT), JA(NELT), ISYM, NSAVE, ITOL, ITMAX
  !     INTEGER ITER, IERR, IUNIT, LENW, IWORK(10), LENIW
  !     DOUBLE PRECISION B(N), X(N), A(NELT), TOL, ERR
  !     DOUBLE PRECISION RWORK(7*N+3*N*NSAVE+NSAVE)
  !
  !     CALL DSDOMN(N, B, X, NELT, IA, JA, A, ISYM, NSAVE, ITOL, TOL,
  !    $     ITMAX, ITER, ERR, IERR, IUNIT, RWORK, LENW, IWORK, LENIW )
  !
  !- Arguments:
  ! N      :IN       Integer.
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
  !         These arrays should hold the matrix A in either the SLAP
  !         Triad format or the SLAP Column format.  See "Description",
  !         below.  If the SLAP Triad format is chosen, it is changed
  !         internally to the SLAP Column format.
  ! ISYM   :IN       Integer.
  !         Flag to indicate symmetric storage format.
  !         If ISYM=0, all non-zero entries of the matrix are stored.
  !         If ISYM=1, the matrix is symmetric, and only the upper
  !         or lower triangle of the matrix is stored.
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
  !         is less than TOL.
  ! TOL    :INOUT    Double Precision.
  !         Convergence criterion, as described above.  (Reset if IERR=4.)
  ! ITMAX  :IN       Integer.
  !         Maximum number of iterations.
  ! ITER   :OUT      Integer.
  !         Number of iterations required to reach convergence, or
  !         ITMAX+1 if convergence criterion could not be achieved in
  !         ITMAX iterations.
  ! ERR    :OUT      Double Precision.
  !         Error estimate of error in final approximate solution, as
  !         defined by ITOL.
  ! IERR   :OUT      Integer.
  !         Return error flag.
  !           IERR = 0 => All went well.
  !           IERR = 1 => Insufficient space allocated for WORK or IWORK.
  !           IERR = 2 => Method failed to converge in ITMAX steps.
  !           IERR = 3 => Error in user input.
  !                       Check input values of N, ITOL.
  !           IERR = 4 => User error tolerance set too tight.
  !                       Reset to 500*eps_2_dp.  Iteration proceeded.
  !           IERR = 5 => Preconditioning matrix, M, is not positive
  !                       definite.  (r,z) < 0.
  !           IERR = 6 => Breakdown of method detected.
  !                       (p,Ap) < epsilon**2.
  ! IUNIT  :IN       Integer.
  !         Unit number on which to write the error at each iteration,
  !         if this is desired for monitoring convergence.  If unit
  !         number is 0, no writing will occur.
  ! RWORK  :WORK     Double Precision RWORK(LENW).
  !         Double Precision array used for workspace.
  ! LENW   :IN       Integer.
  !         Length of the double precision workspace, RWORK.
  !         LENW >= 7*N+NSAVE*(3*N+1).
  ! IWORK  :WORK     Integer IWORK(LENIW).
  !         Used to hold pointers into the RWORK array.
  ! LENIW  :IN       Integer.
  !         Length of the integer workspace, IWORK.  LENIW >= 10.
  !
  !- Description:
  !       This routine  is simply a driver  for  the DOMN routine.  It
  !       calls the DSDS  routine  to set  up the  preconditioning and
  !       then   calls DOMN with the   appropriate   MATVEC and MSOLVE
  !       routines.
  !
  !       The Sparse Linear Algebra Package (SLAP) utilizes two matrix
  !       data structures: 1) the  SLAP Triad  format or  2)  the SLAP
  !       Column format.  The user can hand this routine either of the
  !       of these data structures and SLAP  will figure out  which on
  !       is being used and act accordingly.
  !
  !       =================== S L A P Triad format ===================
  !
  !       In  this   format only the  non-zeros are  stored.  They may
  !       appear  in *ANY* order.   The user  supplies three arrays of
  !       length NELT, where  NELT  is the number  of non-zeros in the
  !       matrix:  (IA(NELT), JA(NELT),  A(NELT)).  For each  non-zero
  !       the  user puts   the row  and  column index   of that matrix
  !       element in the IA and JA arrays.  The  value of the non-zero
  !       matrix  element is  placed in  the corresponding location of
  !       the A  array.  This is  an extremely easy data  structure to
  !       generate.  On  the other hand it  is  not too  efficient  on
  !       vector  computers   for the  iterative  solution  of  linear
  !       systems.  Hence, SLAP  changes this input  data structure to
  !       the SLAP   Column  format for the  iteration (but   does not
  !       change it back).
  !
  !       Here is an example of the  SLAP Triad   storage format for a
  !       5x5 Matrix.  Recall that the entries may appear in any order.
  !
  !           5x5 Matrix      SLAP Triad format for 5x5 matrix on left.
  !                              1  2  3  4  5  6  7  8  9 10 11
  !       |11 12  0  0 15|   A: 51 12 11 33 15 53 55 22 35 44 21
  !       |21 22  0  0  0|  IA:  5  1  1  3  1  5  5  2  3  4  2
  !       | 0  0 33  0 35|  JA:  1  2  1  3  5  3  5  2  5  4  1
  !       | 0  0  0 44  0|
  !       |51  0 53  0 55|
  !
  !       =================== S L A P Column format ==================
  !
  !       In  this format   the non-zeros are    stored counting  down
  !       columns (except  for the diagonal  entry, which must  appear
  !       first  in each "column") and are  stored in the  double pre-
  !       cision array  A. In  other  words,  for each  column  in the
  !       matrix  first put  the diagonal entry in A.  Then put in the
  !       other non-zero  elements going  down the column  (except the
  !       diagonal)  in order.  The IA array  holds the  row index for
  !       each non-zero.  The JA array  holds the offsets into the IA,
  !       A  arrays  for  the  beginning  of  each  column.  That  is,
  !       IA(JA(ICOL)),A(JA(ICOL)) are the first elements of the ICOL-
  !       th column in IA and A, and IA(JA(ICOL+1)-1), A(JA(ICOL+1)-1)
  !       are  the last elements of the ICOL-th column.   Note that we
  !       always have JA(N+1)=NELT+1, where N is the number of columns
  !       in the matrix  and NELT  is the number  of non-zeros  in the
  !       matrix.
  !
  !       Here is an example of the  SLAP Column  storage format for a
  !       5x5 Matrix (in the A and IA arrays '|'  denotes the end of a
  !       column):
  !
  !           5x5 Matrix      SLAP Column format for 5x5 matrix on left.
  !                              1  2  3    4  5    6  7    8    9 10 11
  !       |11 12  0  0 15|   A: 11 21 51 | 22 12 | 33 53 | 44 | 55 15 35
  !       |21 22  0  0  0|  IA:  1  2  5 |  2  1 |  3  5 |  4 |  5  1  3
  !       | 0  0 33  0 35|  JA:  1  4  6    8  9   12
  !       | 0  0  0 44  0|
  !       |51  0 53  0 55|
  !
  !- Side Effects:
  !       The SLAP Triad format (IA, JA, A)  is modified internally to
  !       be the SLAP Column format.  See above.
  !
  !- Cautions:
  !     This routine will attempt to write to the Fortran logical output
  !     unit IUNIT, if IUNIT /= 0.  Thus, the user must make sure that
  !     this logical unit is attached to a file or terminal before calling
  !     this routine with a non-zero value for IUNIT.  This routine does
  !     not check for the validity of a non-zero IUNIT unit number.
  !
  !***
  ! **See also:**  DOMN, DSLUOM
  !***
  ! **References:**  (NONE)
  !***
  ! **Routines called:**  DCHKW, DOMN, DS2Y, DSDI, DSDS, DSMV

  !* REVISION HISTORY  (YYMMDD)
  !   890404  DATE WRITTEN
  !   890404  Previous REVISION DATE
  !   890915  Made changes requested at July 1989 CML Meeting.  (MKS)
  !   890921  Removed TeX from comments.  (FNF)
  !   890922  Numerous changes to prologue to make closer to SLATEC standard.  (FNF)
  !   890929  Numerous changes to reduce SP/DP differences.  (FNF)
  !   910411  Prologue converted to Version 4.0 format.  (BAB)
  !   920407  COMMON BLOCK renamed DSLBLK.  (WRB)
  !   920511  Added complete declaration section.  (WRB)
  !   921113  Corrected C***CATEGORY line.  (FNF)

  !     .. Parameters ..
  INTEGER, PARAMETER :: LOCRB = 1, LOCIB = 11
  !     .. Scalar Arguments ..
  INTEGER, INTENT(IN) :: Isym, Itmax, Itol, Leniw, Lenw, N, Nelt, Nsave
  INTEGER, INTENT(OUT) :: Ierr, Iter
  REAL(DP), INTENT(INOUT) :: Tol
  REAL(DP), INTENT(OUT) :: Err
  !     .. Array Arguments ..
  INTEGER, INTENT(INOUT) :: Ia(Nelt), Ja(Nelt), Iwork(Leniw)
  REAL(DP), INTENT(IN) :: B(N)
  REAL(DP), INTENT(INOUT) :: A(N), X(N), Rwork(Lenw)
  !     .. Local Scalars ..
  INTEGER :: locap, loccsa, locdin, locdz, locema, lociw, locp, locr, locw, locz
  !* FIRST EXECUTABLE STATEMENT  DSDOMN
  !
  Ierr = 0
  IF( N<1 .OR. Nelt<1 ) THEN
    Ierr = 3
    RETURN
  END IF
  !
  !         Change the SLAP input matrix IA, JA, A to SLAP-Column format.
  CALL DS2Y(N,Nelt,Ia,Ja,A)
  !
  !         Set up the workspace.
  lociw = LOCIB
  !
  locdin = LOCRB
  locr = locdin + N
  locz = locr + N
  locp = locz + N
  locap = locp + N*(Nsave+1)
  locema = locap + N*(Nsave+1)
  locdz = locema + N*(Nsave+1)
  loccsa = locdz + N
  locw = loccsa + Nsave
  !
  !         Check the workspace allocations.
  CALL DCHKW('DSDOMN',lociw,Leniw,locw,Lenw,Ierr,Iter,Err)
  IF( Ierr/=0 ) RETURN
  !
  Iwork(4) = locdin
  Iwork(9) = lociw
  Iwork(10) = locw
  !
  !         Compute the inverse of the diagonal of the matrix.
  CALL DSDS(N,Nelt,Ja,A,Rwork(locdin))
  !
  !         Perform the Diagonally Scaled Orthomin iteration algorithm.
  CALL DOMN(N,B,X,Nelt,Ia,Ja,A,Isym,DSMV,DSDI,Nsave,Itol,Tol,Itmax,Iter,&
    Ierr,Rwork(locr),Rwork(locz),Rwork(locp),Rwork(locap),&
    Rwork(locema),Rwork(locdz),Rwork(loccsa),Rwork,Iwork)
  !------------- LAST LINE OF DSDOMN FOLLOWS ----------------------------
END SUBROUTINE DSDOMN