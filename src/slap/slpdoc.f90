!DECK SLPDOC
SUBROUTINE SLPDOC
  IMPLICIT NONE
  !***BEGIN PROLOGUE  SLPDOC
  !***PURPOSE  Sparse Linear Algebra Package Version 2.0.2 Documentation.
  !            Routines to solve large sparse symmetric and nonsymmetric
  !            positive definite linear systems, Ax = b, using precondi-
  !            tioned iterative methods.
  !***LIBRARY   SLATEC (SLAP)
  !***CATEGORY  D2A4, D2B4, Z
  !***TYPE      SINGLE PRECISION (SLPDOC-S, DLPDOC-D)
  !***KEYWORDS  BICONJUGATE GRADIENT SQUARED, DOCUMENTATION,
  !             GENERALIZED MINIMUM RESIDUAL, ITERATIVE IMPROVEMENT,
  !             NORMAL EQUATIONS, ORTHOMIN,
  !             PRECONDITIONED CONJUGATE GRADIENT, SLAP,
  !             SPARSE ITERATIVE METHODS
  !***AUTHOR  Seager, Mark. K., (LLNL)
  !             User Systems Division
  !             Lawrence Livermore National Laboratory
  !             PO BOX 808, L-60
  !             Livermore, CA 94550
  !             (FTS) 543-3141, (510) 423-3141
  !             seager@llnl.gov
  !***DESCRIPTION
  !                                 The
  !                    Sparse Linear Algebra Package
  !
  !                @@@@@@@  @            @@@    @@@@@@@@
  !               @       @ @           @   @   @       @
  !               @         @          @     @  @       @
  !                @@@@@@@  @         @       @ @@@@@@@@
  !                       @ @         @@@@@@@@@ @
  !               @       @ @         @       @ @
  !                @@@@@@@  @@@@@@@@@ @       @ @
  !
  !      @       @                            @@@@@@@        @@@@@
  !      @       @                           @       @      @    @@
  !      @       @  @@@@@@@  @ @@                    @     @    @  @
  !      @       @ @       @ @@  @             @@@@@@      @   @   @
  !       @     @  @@@@@@@@@ @                @            @  @    @
  !        @   @   @         @               @         @@@  @@    @
  !         @@@     @@@@@@@  @               @@@@@@@@@ @@@   @@@@@
  !
  !
  !    =================================================================
  !    ========================== Introduction =========================
  !    =================================================================
  !      This package was  originally derived from a set of  iterative
  !      routines written by Anne Greenbaum, as announced in "Routines
  !      for Solving Large Sparse Linear Systems",  Tentacle, Lawrence
  !      Livermore  National  Laboratory,  Livermore  Computing Center
  !      (January 1986), pp 15-21.
  !
  !    This document  contains the specifications for  the  SLAP Version
  !    2.0 package, a Fortran 77  package  for  the  solution  of  large
  !    sparse   linear systems, Ax  =  b,  via  preconditioned iterative
  !    methods.   Included in  this  package are "core"  routines  to do
  !    Iterative   Refinement  (Jacobi's  method),  Conjugate  Gradient,
  !    Conjugate Gradient on the normal equations, AA'y = b,  (where x =
  !    A'y and  A' denotes the  transpose of   A), BiConjugate Gradient,
  !    BiConjugate  Gradient  Squared, Orthomin and  Generalized Minimum
  !    Residual Iteration.    These "core" routines   do  not  require a
  !    "fixed"   data  structure   for storing  the   matrix  A  and the
  !    preconditioning   matrix  M.   The  user  is free  to  choose any
  !    structure that facilitates  efficient solution  of the problem at
  !    hand.  The drawback  to this approach  is that the user must also
  !    supply at least two routines  (MATVEC and MSOLVE,  say).   MATVEC
  !    must calculate, y = Ax, given x and the user's data structure for
  !    A.  MSOLVE must solve,  r = Mz, for z (*NOT*  r) given r  and the
  !    user's data  structure for  M (or its  inverse).  The user should
  !    choose M so that  inv(M)*A  is approximately the identity and the
  !    solution step r = Mz is "easy" to  solve.  For some of the "core"
  !    routines (Orthomin,  BiConjugate Gradient and  Conjugate Gradient
  !    on the  normal equations)   the user must  also  supply  a matrix
  !    transpose times   vector  routine  (MTTVEC,  say)  and (possibly,
  !    depending    on the "core"  method)   a  routine  that solves the
  !    transpose  of   the   preconditioning    step     (MTSOLV,  say).
  !    Specifically, MTTVEC is a routine which calculates y = A'x, given
  !    x and the user's data structure for A (A' is the transpose of A).
  !    MTSOLV is a routine which solves the system r = M'z for z given r
  !    and the user's data structure for M.
  !
  !    This process of writing the matrix vector operations  can be time
  !    consuming and error  prone.  To alleviate  these problems we have
  !    written drivers   for  the  "core" methods  that  assume the user
  !    supplies one of two specific data structures (SLAP Triad and SLAP
  !    Column format), see  below.  Utilizing these  data structures  we
  !    have augmented   each  "core" method  with   two preconditioners:
  !    Diagonal  Scaling and Incomplete Factorization.  Diagonal scaling
  !    is easy to implement, vectorizes very  well and for problems that
  !    are  not too  ill-conditioned  reduces the  number  of iterations
  !    enough   to warrant its use.  On   the other  hand, an Incomplete
  !    factorization  (Incomplete  Cholesky for  symmetric systems   and
  !    Incomplete LU for nonsymmetric  systems) may  take much longer to
  !    calculate, but it reduces the iteration count (for most problems)
  !    significantly.  Our implementations  of IC and ILU  vectorize for
  !    machines with hardware gather scatter, but the vector lengths can
  !    be quite short if  the  number  of non-zeros  in a column is  not
  !    large.
  !
  !    =================================================================
  !    ==================== Supplied Data Structures ===================
  !    =================================================================
  !    The following describes the data   structures supplied  with  the
  !    package: SLAP Triad and Column formats.
  !
  !    ====================== S L A P Triad format =====================
  !
  !    In the SLAP Triad format only the non-zeros are stored.  They may
  !    appear in *ANY* order.  The user supplies three  arrays of length
  !    NELT, where NELT  is the   number of  non-zeros  in the   matrix:
  !    (IA(NELT),  JA(NELT), A(NELT)).  If  the matrix is symmetric then
  !    one need only store the lower triangle (including  the  diagonal)
  !    and NELT would be the corresponding  number  of non-zeros stored.
  !    For each non-zero the user puts the row and column  index of that
  !    matrix  element   in the  IA  and JA  arrays.  The  value  of the
  !    non-zero matrix element is placed  in  the corresponding location
  !    of  the A array.   This  is an extremely  easy  data structure to
  !    generate.  On the other hand, it is not very  efficient on vector
  !    computers for the iterative  solution of  linear systems.  Hence,
  !    SLAP changes this input data structure to  the SLAP Column format
  !    for the iteration (but does not change it back).
  !
  !    Here  is an example   of  the  SLAP  Triad storage  format  for a
  !    nonsymmetric 5x5 Matrix.  NELT=11.   Recall that the  entries may
  !    appear in any order.
  !
  !     5x5 Matrix       SLAP Triad format for 5x5 matrix on left.
  !                           1  2  3  4  5  6  7  8  9 10 11
  !    |11 12  0  0 15|   A: 51 12 11 33 15 53 55 22 35 44 21
  !    |21 22  0  0  0|  IA:  5  1  1  3  1  5  5  2  3  4  2
  !    | 0  0 33  0 35|  JA:  1  2  1  3  5  3  5  2  5  4  1
  !    | 0  0  0 44  0|
  !    |51  0 53  0 55|
  !
  !    ====================== S L A P Column format ====================
  !
  !    In the SLAP Column format  the non-zeros are stored counting down
  !    columns (except for the  diagonal entry,  which must appear first
  !    in each "column") and are stored  in the real array A.   In other
  !    words, for each column in the matrix first put the diagonal entry
  !    in A.   Then put in  the other non-zero  elements going  down the
  !    column (except the  diagonal) in order.   The IA  array holds the
  !    row index for each non-zero.  The JA array holds the offsets into
  !    the  IA,  A arrays for the  beginning  of each column.   That is,
  !    IA(JA(ICOL)), A(JA(ICOL)) are the  first elements of the  ICOL-th
  !    column in IA   and A.  IA(JA(ICOL+1)-1),  A(JA(ICOL+1)-1) are the
  !    last elements  of  the ICOL-th column.   Note that we always have
  !    JA(N+1) = NELT+1, where N is the number of  columns in the matrix
  !    and NELT is the number of non-zeros in the matrix.  If the matrix
  !    is  symmetric one need  only store the  lower triangle (including
  !    the diagonal) and  NELT  would be the   corresponding   number of
  !    non-zeros stored.
  !
  !    Here is  an  example of the  SLAP   Column storage format  for  a
  !    nonsymmetric 5x5 Matrix (in the  A and  IA arrays '|' denotes the
  !    end of a column):
  !
  !       5x5 Matrix      SLAP Column format for 5x5 matrix on left.
  !                           1  2  3    4  5    6  7    8    9 10 11
  !    |11 12  0  0 15|   A: 11 21 51 | 22 12 | 33 53 | 44 | 55 15 35
  !    |21 22  0  0  0|  IA:  1  2  5 |  2  1 |  3  5 |  4 |  5  1  3
  !    | 0  0 33  0 35|  JA:  1  4  6    8  9   12
  !    | 0  0  0 44  0|
  !    |51  0 53  0 55|
  !
  !    =================================================================
  !    ====================== Which Method To Use ======================
  !    =================================================================
  !
  !                          BACKGROUND
  !    In solving a large sparse linear system Ax = b using an iterative
  !    method, it   is  not necessary to actually   store  the matrix A.
  !    Rather, what is needed is a procedure  for multiplying the matrix
  !    A times a given vector y to obtain the matrix-vector product, Ay.
  !    SLAP has been written to take advantage of this fact.  The higher
  !    level routines in the package require storage only of the non-zero
  !    elements of   A (and  their  positions), and  even this   can  be
  !    avoided, if the  user  writes his own subroutine for  multiplying
  !    the matrix times a vector  and   calls the lower-level  iterative
  !    routines in the package.
  !
  !    If  the matrix A is ill-conditioned,  then most iterative methods
  !    will be slow to converge (if they converge  at all!).  To improve
  !    the  convergence  rate,  one  may use  a "matrix  splitting," or,
  !    "preconditioning matrix," say, M.  It is then necessary to solve,
  !    at each iteration, a linear system  with coefficient matrix M.  A
  !    good preconditioner  M should have  two  properties: (1) M should
  !    "approximate" A, in the sense that the  matrix inv(M)*A  (or some
  !    variant  thereof) is better conditioned  than the original matrix
  !    A; and  (2) linear  systems with coefficient  matrix M should  be
  !    much easier  to solve  than  the original system with coefficient
  !    matrix   A.   Preconditioning routines  in the   SLAP package are
  !    separate from the  iterative   routines,  so   that any of    the
  !    preconditioners provided in the package,   or one that the   user
  !    codes himself, can be used with any of the iterative routines.
  !
  !                        CHOICE OF PRECONDITIONER
  !    If you  willing   to live with   either the SLAP Triad or  Column
  !    matrix data structure  you  can then  choose one  of two types of
  !    preconditioners   to   use:   diagonal  scaling    or  incomplete
  !    factorization.  To  choose   between these two   methods requires
  !    knowing  something  about the computer you're going  to run these
  !    codes on  and how well incomplete factorization  approximates the
  !    inverse of your matrix.
  !
  !    Let us  suppose you have   a scalar  machine.   Then,  unless the
  !    incomplete factorization is very,  very poor this  is *GENERALLY*
  !    the method to choose.  It  will reduce the  number of  iterations
  !    significantly and is not all  that expensive  to compute.  So  if
  !    you have just one  linear system to solve  and  "just want to get
  !    the job  done" then try  incomplete factorization first.   If you
  !    are thinking of integrating some SLAP  iterative method into your
  !    favorite   "production  code" then  try incomplete  factorization
  !    first,  but  also check  to see that  diagonal  scaling is indeed
  !    slower for a large sample of test problems.
  !
  !    Let us now suppose  you have  a  vector  computer  with  hardware
  !    gather/scatter support (Cray X-MP, Y-MP, SCS-40 or Cyber 205, ETA
  !    10,  ETA Piper,  Convex C-1,  etc.).   Then  it is much harder to
  !    choose  between the  two  methods.   The  versions  of incomplete
  !    factorization in SLAP do in fact vectorize, but have short vector
  !    lengths and the factorization step is relatively  more expensive.
  !    Hence,  for  most problems (i.e.,  unless  your  problem  is  ill
  !    conditioned,  sic!)  diagonal  scaling is  faster,  with its very
  !    fast    set up  time    and  vectorized  (with   long    vectors)
  !    preconditioning step (even though  it  may take more iterations).
  !    If you have several systems (or  right hand sides) to  solve that
  !    can  utilize  the  same  preconditioner  then the   cost   of the
  !    incomplete factorization can   be  amortized over these  several
  !    solutions.  This situation gives more advantage to the incomplete
  !    factorization methods.  If  you have  a  vector  machine  without
  !    hardware  gather/scatter (Cray  1,  Cray  2  &  Cray 3) then  the
  !    advantages for incomplete factorization are even less.
  !
  !    If you're trying to shoehorn SLAP into your  favorite "production
  !    code" and can not easily generate either the SLAP Triad or Column
  !    format  then  you are  left  to   your  own  devices in terms  of
  !    preconditioning.  Also,  you may  find that the   preconditioners
  !    supplied with SLAP are not sufficient  for your problem.  In this
  !    situation we would  recommend  that you   talk  with a  numerical
  !    analyst  versed in   iterative   methods   about   writing  other
  !    preconditioning  subroutines (e.g.,  polynomial  preconditioning,
  !    shifted incomplete factorization,  SOR  or SSOR  iteration).  You
  !    can always "roll your own"  by using the "core" iterative methods
  !    and supplying your own MSOLVE and MATVEC (and possibly MTSOLV and
  !    MTTVEC) routines.
  !
  !                          SYMMETRIC SYSTEMS
  !    If your matrix is symmetric then you would want to use one of the
  !    symmetric system  solvers.    If  your  system  is  also positive
  !    definite,   (Ax,x) (Ax dot  product  with x) is  positive for all
  !    non-zero  vectors x,  then use   Conjugate Gradient (SCG,  SSDCG,
  !    SSICSG).  If you're  not sure it's SPD   (symmetric and  Positive
  !    Definite)  then try SCG anyway and  if it works, fine.  If you're
  !    sure your matrix is not  positive definite  then you  may want to
  !    try the iterative refinement   methods  (SIR)  or the  GMRES code
  !    (SGMRES) if SIR converges too slowly.
  !
  !                         NONSYMMETRIC SYSTEMS
  !    This   is currently  an  area  of  active research  in  numerical
  !    analysis  and   there   are   new  strategies  being   developed.
  !    Consequently take the following advice with a grain of salt.   If
  !    you matrix is positive definite, (Ax,x)  (Ax  dot product  with x
  !    is positive for all non-zero  vectors x), then you can use any of
  !    the    methods   for   nonsymmetric   systems (Orthomin,   GMRES,
  !    BiConjugate Gradient, BiConjugate Gradient  Squared and Conjugate
  !    Gradient applied to the normal equations).  If your system is not
  !    too ill conditioned then try  BiConjugate Gradient Squared (BCGS)
  !    or GMRES (SGMRES).  Both  of  these methods converge very quickly
  !    and do  not require A'  or M' ('  denotes transpose) information.
  !    SGMRES  does require  some  additional storage,  though.  If  the
  !    system is very  ill conditioned  or   nearly positive  indefinite
  !    ((Ax,x) is positive,  but may be  very small),  then GMRES should
  !    be the first choice,  but try the  other  methods  if you have to
  !    fine tune  the solution process for a  "production code".  If you
  !    have a great preconditioner for the normal  equations (i.e., M is
  !    an approximation to the inverse of AA' rather than  just  A) then
  !    this is not a bad route to travel.  Old wisdom would say that the
  !    normal equations are a disaster  (since it squares the  condition
  !    number of the system and SCG convergence is linked to this number
  !    of    infamy), but   some     preconditioners    (like incomplete
  !    factorization) can reduce the condition number back below that of
  !    the original system.
  !
  !    =================================================================
  !    ======================= Naming Conventions ======================
  !    =================================================================
  !    SLAP  iterative  methods,    matrix vector    and  preconditioner
  !    calculation  routines   follow a naming   convention  which, when
  !    understood, allows one to determine the iterative method and data
  !    structure(s) used.  The  subroutine  naming convention  takes the
  !    following form:
  !                          P[S][M]DESC
  !    where
  !        P  stands for the precision (or data type) of the routine and
  !           is required in all names,
  !        S  denotes whether or not the routine requires the SLAP Triad
  !           or Column format (it does if the second letter of the name
  !           is S and does not otherwise),
  !        M  stands for the type of preconditioner used (only appears
  !           in drivers for "core" routines), and
  !     DESC  is some number of letters describing the method or purpose
  !           of the routine.  The following is a list of the "DESC"
  !           fields for iterative methods and their meaning:
  !             BCG,BC:       BiConjugate Gradient
  !             CG:           Conjugate Gradient
  !             CGN,CN:       Conjugate Gradient on the Normal equations
  !             CGS,CS:       biConjugate Gradient Squared
  !             GMRES,GMR,GM: Generalized Minimum RESidual
  !             IR,R:         Iterative Refinement
  !             JAC:          JACobi's method
  !             GS:           Gauss-Seidel
  !             OMN,OM:       OrthoMiN
  !
  !    In the single precision version of SLAP, all routine names start
  !    with an S. The brackets around the S and M designate that these
  !    fields are optional.
  !
  !    Here are some examples of the routines:
  !    1) SBCG: Single precision BiConjugate Gradient "core" routine.
  !       One can deduce that this is a "core" routine, because the S and
  !       M fields are missing and BiConjugate Gradient is an iterative
  !       method.
  !    2) SSDBCG: Single precision, SLAP data structure BCG with Diagonal
  !       scaling.
  !    3) SSLUBC: Single precision, SLAP data structure BCG with incom-
  !       plete LU factorization as the preconditioning.
  !    4) SCG: Single precision Conjugate Gradient "core" routine.
  !    5) SSDCG: Single precision, SLAP data structure Conjugate Gradient
  !       with Diagonal scaling.
  !    6) SSICCG: Single precision, SLAP data structure Conjugate Gra-
  !       dient with Incomplete Cholesky factorization preconditioning.
  !
  !
  !    =================================================================
  !    ===================== USER CALLABLE ROUTINES ====================
  !    =================================================================
  !    The following is a list of  the "user callable" SLAP routines and
  !    their one line descriptions.  The headers denote  the  file names
  !    where the routines can be found, as distributed for UNIX systems.
  !
  !    Note:  Each core routine, SXXX, has a corresponding stop routine,
  !         ISSXXX.  If the stop routine does not have the specific stop
  !         test the user requires (e.g., weighted infinity norm),  then
  !         the user should modify the source for ISSXXX accordingly.
  !
  !    ============================= sir.f =============================
  !    SIR: Preconditioned Iterative Refinement Sparse Ax = b Solver.
  !    SSJAC: Jacobi's Method Iterative Sparse Ax = b Solver.
  !    SSGS: Gauss-Seidel Method Iterative Sparse Ax = b Solver.
  !    SSILUR: Incomplete LU Iterative Refinement Sparse Ax = b Solver.
  !
  !    ============================= scg.f =============================
  !    SCG: Preconditioned Conjugate Gradient Sparse Ax=b Solver.
  !    SSDCG: Diagonally Scaled Conjugate Gradient Sparse Ax=b Solver.
  !    SSICCG: Incomplete Cholesky Conjugate Gradient Sparse Ax=b Solver.
  !
  !    ============================= scgn.f ============================
  !    SCGN: Preconditioned CG Sparse Ax=b Solver for Normal Equations.
  !    SSDCGN: Diagonally Scaled CG Sparse Ax=b Solver for Normal Eqn's.
  !    SSLUCN: Incomplete LU CG Sparse Ax=b Solver for Normal Equations.
  !
  !    ============================= sbcg.f ============================
  !    SBCG: Preconditioned BiConjugate Gradient Sparse Ax = b Solver.
  !    SSDBCG: Diagonally Scaled BiConjugate Gradient Sparse Ax=b Solver.
  !    SSLUBC: Incomplete LU BiConjugate Gradient Sparse Ax=b Solver.
  !
  !    ============================= scgs.f ============================
  !    SCGS: Preconditioned BiConjugate Gradient Squared Ax=b Solver.
  !    SSDCGS: Diagonally Scaled CGS Sparse Ax=b Solver.
  !    SSLUCS: Incomplete LU BiConjugate Gradient Squared Ax=b Solver.
  !
  !    ============================= somn.f ============================
  !    SOMN: Preconditioned Orthomin Sparse Iterative Ax=b Solver.
  !    SSDOMN: Diagonally Scaled Orthomin Sparse Iterative Ax=b Solver.
  !    SSLUOM: Incomplete LU Orthomin Sparse Iterative Ax=b Solver.
  !
  !    ============================ sgmres.f ===========================
  !    SGMRES: Preconditioned GMRES Iterative Sparse Ax=b Solver.
  !    SSDGMR: Diagonally Scaled GMRES Iterative Sparse Ax=b Solver.
  !    SSLUGM: Incomplete LU GMRES Iterative Sparse Ax=b Solver.
  !
  !    ============================ smset.f ============================
  !       The following routines are used to set up preconditioners.
  !
  !    SSDS: Diagonal Scaling Preconditioner SLAP Set Up.
  !    SSDSCL: Diagonally Scales/Unscales a SLAP Column Matrix.
  !    SSD2S: Diagonal Scaling Preconditioner SLAP Normal Eqns Set Up.
  !    SS2LT: Lower Triangle Preconditioner SLAP Set Up.
  !    SSICS: Incomplete Cholesky Decomp. Preconditioner SLAP Set Up.
  !    SSILUS: Incomplete LU Decomposition Preconditioner SLAP Set Up.
  !
  !    ============================ smvops.f ===========================
  !       Most of the incomplete  factorization  (LL' and LDU) solvers
  !       in this  file require an  intermediate routine  to translate
  !       from the SLAP MSOLVE(N, R, Z, NELT, IA,  JA, A, ISYM, RWORK,
  !       IWORK) calling  convention to the calling  sequence required
  !       by  the solve routine.   This generally  is  accomplished by
  !       fishing out pointers to the preconditioner (stored in RWORK)
  !       from the  IWORK  array and then making a call to the routine
  !       that actually does the backsolve.
  !
  !    SSMV: SLAP Column Format Sparse Matrix Vector Product.
  !    SSMTV: SLAP Column Format Sparse Matrix (transpose) Vector Prod.
  !    SSDI: Diagonal Matrix Vector Multiply.
  !    SSLI: SLAP MSOLVE for Lower Triangle Matrix (set up for SSLI2).
  !    SSLI2: Lower Triangle Matrix Backsolve.
  !    SSLLTI: SLAP MSOLVE for LDL' (IC) Fact. (set up for SLLTI2).
  !    SLLTI2: Backsolve routine for LDL' Factorization.
  !    SSLUI: SLAP MSOLVE for LDU Factorization (set up for SSLUI2).
  !    SSLUI2: SLAP Backsolve for LDU Factorization.
  !    SSLUTI: SLAP MTSOLV for LDU Factorization (set up for SSLUI4).
  !    SSLUI4: SLAP Backsolve for LDU Factorization.
  !    SSMMTI: SLAP MSOLVE for LDU Fact of Normal Eq (set up for SSMMI2).
  !    SSMMI2: SLAP Backsolve for LDU Factorization of Normal Equations.
  !
  !    =========================== slaputil.f ==========================
  !       The following utility routines are useful additions to SLAP.
  !
  !    SBHIN: Read Sparse Linear System in the Boeing/Harwell Format.
  !    SCHKW: SLAP WORK/IWORK Array Bounds Checker.
  !    SCPPLT: Printer Plot of SLAP Column Format Matrix.
  !    SS2Y: SLAP Triad to SLAP Column Format Converter.
  !    QS2I1R: Quick Sort Integer array, moving integer and real arrays.
  !            (Used by SS2Y.)
  !    STIN: Read in SLAP Triad Format Linear System.
  !    STOUT: Write out SLAP Triad Format Linear System.
  !
  !
  !***REFERENCES  1. Mark K. Seager, A SLAP for the Masses, in
  !                  G. F. Carey, Ed., Parallel Supercomputing: Methods,
  !                  Algorithms and Applications, Wiley, 1989, pp.135-155.
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   880715  DATE WRITTEN
  !   890404  Previous REVISION DATE
  !   890915  Made changes requested at July 1989 CML Meeting.  (MKS)
  !   890921  Removed TeX from comments.  (FNF)
  !   890922  Numerous changes to prologue to make closer to SLATEC
  !           standard.  (FNF)
  !   890929  Numerous changes to reduce SP/DP differences.  (FNF)
  !           -----( This produced Version 2.0.1. )-----
  !   891003  Rearranged list of user callable routines to agree with
  !           order in source deck.  (FNF)
  !   891004  Updated reference.
  !   910411  Prologue converted to Version 4.0 format.  (BAB)
  !           -----( This produced Version 2.0.2. )-----
  !   910506  Minor improvements to prologue.  (FNF)
  !   920511  Added complete declaration section.  (WRB)
  !   920929  Corrected format of reference.  (FNF)
  !   921019  Improved one-line descriptions, reordering some.  (FNF)
  !***END PROLOGUE  SLPDOC
  !***FIRST EXECUTABLE STATEMENT  SLPDOC
  !
  !     This is a *DUMMY* subroutine and should never be called.
  !
  !------------- LAST LINE OF SLPDOC FOLLOWS -----------------------------
END SUBROUTINE SLPDOC
