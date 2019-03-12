!DECK SDRIV3
SUBROUTINE SDRIV3(N,T,Y,F,Nstate,Tout,Ntask,Nroot,Eps,Ewt,Ierror,Mint,&
    Miter,Impl,Ml,Mu,Mxord,Hmax,Work,Lenw,Iwork,Leniw,&
    JACOBN,FA,Nde,Mxstep,G,USERS,Ierflg)
  IMPLICIT NONE
  !***BEGIN PROLOGUE  SDRIV3
  !***PURPOSE  The function of SDRIV3 is to solve N ordinary differential
  !            equations of the form dY(I)/dT = F(Y(I),T), given the
  !            initial conditions Y(I) = YI.  The program has options to
  !            allow the solution of both stiff and non-stiff differential
  !            equations.  Other important options are available.  SDRIV3
  !            uses single precision arithmetic.
  !***LIBRARY   SLATEC (SDRIVE)
  !***CATEGORY  I1A2, I1A1B
  !***TYPE      SINGLE PRECISION (SDRIV3-S, DDRIV3-D, CDRIV3-C)
  !***KEYWORDS  GEAR'S METHOD, INITIAL VALUE PROBLEMS, ODE,
  !             ORDINARY DIFFERENTIAL EQUATIONS, SDRIVE, SINGLE PRECISION,
  !             STIFF
  !***AUTHOR  Kahaner, D. K., (NIST)
  !             National Institute of Standards and Technology
  !             Gaithersburg, MD  20899
  !           Sutherland, C. D., (LANL)
  !             Mail Stop D466
  !             Los Alamos National Laboratory
  !             Los Alamos, NM  87545
  !***DESCRIPTION
  !
  !  I.  ABSTRACT  .......................................................
  !
  !    The primary function of SDRIV3 is to solve N ordinary differential
  !    equations of the form dY(I)/dT = F(Y(I),T), given the initial
  !    conditions Y(I) = YI.  The program has options to allow the
  !    solution of both stiff and non-stiff differential equations.  In
  !    addition, SDRIV3 may be used to solve:
  !      1. The initial value problem, A*dY(I)/dT = F(Y(I),T), where A is
  !         a non-singular matrix depending on Y and T.
  !      2. The hybrid differential/algebraic initial value problem,
  !         A*dY(I)/dT = F(Y(I),T), where A is a vector (whose values may
  !         depend upon Y and T) some of whose components will be zero
  !         corresponding to those equations which are algebraic rather
  !         than differential.
  !    SDRIV3 is to be called once for each output point of T.
  !
  !  II.  PARAMETERS  ....................................................
  !
  !    The user should use parameter names in the call sequence of SDRIV3
  !    for those quantities whose value may be altered by SDRIV3.  The
  !    parameters in the call sequence are:
  !
  !    N      = (Input) The number of dependent functions whose solution
  !             is desired.  N must not be altered during a problem.
  !
  !    T      = The independent variable.  On input for the first call, T
  !             is the initial point.  On output, T is the point at which
  !             the solution is given.
  !
  !    Y      = The vector of dependent variables.  Y is used as input on
  !             the first call, to set the initial values.  On output, Y
  !             is the computed solution vector.  This array Y is passed
  !             in the call sequence of the user-provided routines F,
  !             JACOBN, FA, USERS, and G.  Thus parameters required by
  !             those routines can be stored in this array in components
  !             N+1 and above.  (Note: Changes by the user to the first
  !             N components of this array will take effect only after a
  !             restart, i.e., after setting NSTATE to 1 .)
  !
  !    F      = A subroutine supplied by the user.  The name must be
  !             declared EXTERNAL in the user's calling program.  This
  !             subroutine is of the form:
  !                   SUBROUTINE F (N, T, Y, YDOT)
  !                   REAL Y(*), YDOT(*)
  !                     .
  !                     .
  !                   YDOT(1) = ...
  !                     .
  !                     .
  !                   YDOT(N) = ...
  !                   END (Sample)
  !             This computes YDOT = F(Y,T), the right hand side of the
  !             differential equations.  Here Y is a vector of length at
  !             least N.  The actual length of Y is determined by the
  !             user's declaration in the program which calls SDRIV3.
  !             Thus the dimensioning of Y in F, while required by FORTRAN
  !             convention, does not actually allocate any storage.  When
  !             this subroutine is called, the first N components of Y are
  !             intermediate approximations to the solution components.
  !             The user should not alter these values.  Here YDOT is a
  !             vector of length N.  The user should only compute YDOT(I)
  !             for I from 1 to N.  Normally a return from F passes
  !             control back to  SDRIV3.  However, if the user would like
  !             to abort the calculation, i.e., return control to the
  !             program which calls SDRIV3, he should set N to zero.
  !             SDRIV3 will signal this by returning a value of NSTATE
  !             equal to 6 .  Altering the value of N in F has no effect
  !             on the value of N in the call sequence of SDRIV3.
  !
  !    NSTATE = An integer describing the status of integration.  The
  !             meaning of NSTATE is as follows:
  !               1  (Input) Means the first call to the routine.  This
  !                  value must be set by the user.  On all subsequent
  !                  calls the value of NSTATE should be tested by the
  !                  user, but must not be altered.  (As a convenience to
  !                  the user who may wish to put out the initial
  !                  conditions, SDRIV3 can be called with NSTATE=1, and
  !                  TOUT=T.  In this case the program will return with
  !                  NSTATE unchanged, i.e., NSTATE=1.)
  !               2  (Output) Means a successful integration.  If a normal
  !                  continuation is desired (i.e., a further integration
  !                  in the same direction), simply advance TOUT and call
  !                  again.  All other parameters are automatically set.
  !               3  (Output)(Unsuccessful) Means the integrator has taken
  !                  MXSTEP steps without reaching TOUT.  The user can
  !                  continue the integration by simply calling SDRIV3
  !                  again.
  !               4  (Output)(Unsuccessful) Means too much accuracy has
  !                  been requested.  EPS has been increased to a value
  !                  the program estimates is appropriate.  The user can
  !                  continue the integration by simply calling SDRIV3
  !                  again.
  !               5  (Output) A root was found at a point less than TOUT.
  !                  The user can continue the integration toward TOUT by
  !                  simply calling SDRIV3 again.
  !               6  (Output)(Unsuccessful) N has been set to zero in
  !                  SUBROUTINE F.
  !               7  (Output)(Unsuccessful) N has been set to zero in
  !                  FUNCTION G.  See description of G below.
  !               8  (Output)(Unsuccessful) N has been set to zero in
  !                  SUBROUTINE JACOBN.  See description of JACOBN below.
  !               9  (Output)(Unsuccessful) N has been set to zero in
  !                  SUBROUTINE FA.  See description of FA below.
  !              10  (Output)(Unsuccessful) N has been set to zero in
  !                  SUBROUTINE USERS.  See description of USERS below.
  !              11  (Output)(Successful) For NTASK = 2 or 3, T is beyond
  !                  TOUT.  The solution was obtained by interpolation.
  !                  The user can continue the integration by simply
  !                  advancing TOUT and calling SDRIV3 again.
  !              12  (Output)(Unsuccessful) The solution could not be
  !                  obtained.  The value of IERFLG (see description
  !                  below) for a "Recoverable" situation indicates the
  !                  type of difficulty encountered: either an illegal
  !                  value for a parameter or an inability to continue the
  !                  solution.  For this condition the user should take
  !                  corrective action and reset NSTATE to 1 before
  !                  calling SDRIV3 again.  Otherwise the program will
  !                  terminate the run.
  !
  !    TOUT   = (Input) The point at which the solution is desired.  The
  !             position of TOUT relative to T on the first call
  !             determines the direction of integration.
  !
  !    NTASK  = (Input) An index specifying the manner of returning the
  !             solution, according to the following:
  !               NTASK = 1  Means SDRIV3 will integrate past TOUT and
  !                          interpolate the solution.  This is the most
  !                          efficient mode.
  !               NTASK = 2  Means SDRIV3 will return the solution after
  !                          each internal integration step, or at TOUT,
  !                          whichever comes first.  In the latter case,
  !                          the program integrates exactly to TOUT.
  !               NTASK = 3  Means SDRIV3 will adjust its internal step to
  !                          reach TOUT exactly (useful if a singularity
  !                          exists beyond TOUT.)
  !
  !    NROOT  = (Input) The number of equations whose roots are desired.
  !             If NROOT is zero, the root search is not active.  This
  !             option is useful for obtaining output at points which are
  !             not known in advance, but depend upon the solution, e.g.,
  !             when some solution component takes on a specified value.
  !             The root search is carried out using the user-written
  !             function G (see description of G below.)  SDRIV3 attempts
  !             to find the value of T at which one of the equations
  !             changes sign.  SDRIV3 can find at most one root per
  !             equation per internal integration step, and will then
  !             return the solution either at TOUT or at a root, whichever
  !             occurs first in the direction of integration.  The initial
  !             point is never reported as a root.  The index of the
  !             equation whose root is being reported is stored in the
  !             sixth element of IWORK.
  !             NOTE: NROOT is never altered by this program.
  !
  !    EPS    = On input, the requested relative accuracy in all solution
  !             components.  EPS = 0 is allowed.  On output, the adjusted
  !             relative accuracy if the input value was too small.  The
  !             value of EPS should be set as large as is reasonable,
  !             because the amount of work done by SDRIV3 increases as EPS
  !             decreases.
  !
  !    EWT    = (Input) Problem zero, i.e., the smallest, nonzero,
  !             physically meaningful value for the solution.  (Array,
  !             possibly of length one.  See following description of
  !             IERROR.)  Setting EWT smaller than necessary can adversely
  !             affect the running time.
  !
  !    IERROR = (Input) Error control indicator.  A value of 3 is
  !             suggested for most problems.  Other choices and detailed
  !             explanations of EWT and IERROR are given below for those
  !             who may need extra flexibility.
  !
  !             These last three input quantities EPS, EWT and IERROR
  !             control the accuracy of the computed solution.  EWT and
  !             IERROR are used internally to compute an array YWT.  One
  !             step error estimates divided by YWT(I) are kept less than
  !             EPS in root mean square norm.
  !                 IERROR (Set by the user) =
  !                 1  Means YWT(I) = 1. (Absolute error control)
  !                                   EWT is ignored.
  !                 2  Means YWT(I) = ABS(Y(I)),  (Relative error control)
  !                                   EWT is ignored.
  !                 3  Means YWT(I) = MAX(ABS(Y(I)), EWT(1)).
  !                 4  Means YWT(I) = MAX(ABS(Y(I)), EWT(I)).
  !                    This choice is useful when the solution components
  !                    have differing scales.
  !                 5  Means YWT(I) = EWT(I).
  !             If IERROR is 3, EWT need only be dimensioned one.
  !             If IERROR is 4 or 5, the user must dimension EWT at least
  !             N, and set its values.
  !
  !    MINT   = (Input) The integration method indicator.
  !               MINT = 1  Means the Adams methods, and is used for
  !                         non-stiff problems.
  !               MINT = 2  Means the stiff methods of Gear (i.e., the
  !                         backward differentiation formulas), and is
  !                         used for stiff problems.
  !               MINT = 3  Means the program dynamically selects the
  !                         Adams methods when the problem is non-stiff
  !                         and the Gear methods when the problem is
  !                         stiff.  When using the Adams methods, the
  !                         program uses a value of MITER=0; when using
  !                         the Gear methods, the program uses the value
  !                         of MITER provided by the user.  Only a value
  !                         of IMPL = 0 and a value of MITER = 1, 2, 4, or
  !                         5 is allowed for this option.  The user may
  !                         not alter the value of MINT or MITER without
  !                         restarting, i.e., setting NSTATE to 1.
  !
  !    MITER  = (Input) The iteration method indicator.
  !               MITER = 0  Means functional iteration.  This value is
  !                          suggested for non-stiff problems.
  !               MITER = 1  Means chord method with analytic Jacobian.
  !                          In this case, the user supplies subroutine
  !                          JACOBN (see description below).
  !               MITER = 2  Means chord method with Jacobian calculated
  !                          internally by finite differences.
  !               MITER = 3  Means chord method with corrections computed
  !                          by the user-written routine USERS (see
  !                          description of USERS below.)  This option
  !                          allows all matrix algebra and storage
  !                          decisions to be made by the user.  When using
  !                          a value of MITER = 3, the subroutine FA is
  !                          not required, even if IMPL is not 0.  For
  !                          further information on using this option, see
  !                          Section IV-E below.
  !               MITER = 4  Means the same as MITER = 1 but the A and
  !                          Jacobian matrices are assumed to be banded.
  !               MITER = 5  Means the same as MITER = 2 but the A and
  !                          Jacobian matrices are assumed to be banded.
  !
  !    IMPL   = (Input) The implicit method indicator.
  !               IMPL = 0    Means solving dY(I)/dT = F(Y(I),T).
  !               IMPL = 1    Means solving A*dY(I)/dT = F(Y(I),T), non-
  !                           singular A (see description of FA below.)
  !                           Only MINT = 1 or 2, and MITER = 1, 2, 3, 4,
  !                           or 5 are allowed for this option.
  !               IMPL = 2,3  Means solving certain systems of hybrid
  !                           differential/algebraic equations (see
  !                           description of FA below.)  Only MINT = 2 and
  !                           MITER = 1, 2, 3, 4, or 5, are allowed for
  !                           this option.
  !               The value of IMPL must not be changed during a problem.
  !
  !    ML     = (Input) The lower half-bandwidth in the case of a banded
  !             A or Jacobian matrix.  (I.e., maximum(R-C) for nonzero
  !             A(R,C).)
  !
  !    MU     = (Input) The upper half-bandwidth in the case of a banded
  !             A or Jacobian matrix.  (I.e., maximum(C-R).)
  !
  !    MXORD  = (Input) The maximum order desired. This is .LE. 12 for
  !             the Adams methods and .LE. 5 for the Gear methods.  Normal
  !             value is 12 and 5, respectively.  If MINT is 3, the
  !             maximum order used will be MIN(MXORD, 12) when using the
  !             Adams methods, and MIN(MXORD, 5) when using the Gear
  !             methods.  MXORD must not be altered during a problem.
  !
  !    HMAX   = (Input) The maximum magnitude of the step size that will
  !             be used for the problem.  This is useful for ensuring that
  !             important details are not missed.  If this is not the
  !             case, a large value, such as the interval length, is
  !             suggested.
  !
  !    WORK
  !    LENW   = (Input)
  !             WORK is an array of LENW real words used
  !             internally for temporary storage.  The user must allocate
  !             space for this array in the calling program by a statement
  !             such as
  !                       REAL WORK(...)
  !             The following table gives the required minimum value for
  !             the length of WORK, depending on the value of IMPL and
  !             MITER.  LENW should be set to the value used.  The
  !             contents of WORK should not be disturbed between calls to
  !             SDRIV3.
  !
  !      IMPL =   0            1               2             3
  !              ---------------------------------------------------------
  ! MITER =  0   (MXORD+4)*N   Not allowed     Not allowed   Not allowed
  !              + 2*NROOT
  !              + 250
  !
  !         1,2  N*N +         2*N*N +         N*N +         N*(N + NDE)
  !              (MXORD+5)*N   (MXORD+5)*N     (MXORD+6)*N   + (MXORD+5)*N
  !              + 2*NROOT     + 2*NROOT       + 2*NROOT     + 2*NROOT
  !              + 250         + 250           + 250         + 250
  !
  !          3   (MXORD+4)*N   (MXORD+4)*N     (MXORD+4)*N   (MXORD+4)*N
  !              + 2*NROOT     + 2*NROOT       + 2*NROOT     + 2*NROOT
  !              + 250         + 250           + 250         + 250
  !
  !         4,5  (2*ML+MU+1)   2*(2*ML+MU+1)   (2*ML+MU+1)   (2*ML+MU+1)*
  !              *N +          *N +            *N +          (N+NDE) +
  !              (MXORD+5)*N   (MXORD+5)*N     (MXORD+6)*N   + (MXORD+5)*N
  !              + 2*NROOT     + 2*NROOT       + 2*NROOT     + 2*NROOT
  !              + 250         + 250           + 250         + 250
  !              ---------------------------------------------------------
  !
  !    IWORK
  !    LENIW  = (Input)
  !             IWORK is an integer array of length LENIW used internally
  !             for temporary storage.  The user must allocate space for
  !             this array in the calling program by a statement such as
  !                       INTEGER IWORK(...)
  !             The length of IWORK should be at least
  !               50      if MITER is 0 or 3, or
  !               N+50    if MITER is 1, 2, 4, or 5, or MINT is 3,
  !             and LENIW should be set to the value used.  The contents
  !             of IWORK should not be disturbed between calls to SDRIV3.
  !
  !    JACOBN = A subroutine supplied by the user, if MITER is 1 or 4.
  !             If this is the case, the name must be declared EXTERNAL in
  !             the user's calling program.  Given a system of N
  !             differential equations, it is meaningful to speak about
  !             the partial derivative of the I-th right hand side with
  !             respect to the J-th dependent variable.  In general there
  !             are N*N such quantities.  Often however the equations can
  !             be ordered so that the I-th differential equation only
  !             involves dependent variables with index near I, e.g., I+1,
  !             I-2.  Such a system is called banded.  If, for all I, the
  !             I-th equation depends on at most the variables
  !               Y(I-ML), Y(I-ML+1), ..., Y(I), Y(I+1), ..., Y(I+MU)
  !             then we call ML+MU+1 the bandwidth of the system.  In a
  !             banded system many of the partial derivatives above are
  !             automatically zero.  For the cases MITER = 1, 2, 4, and 5,
  !             some of these partials are needed.  For the cases
  !             MITER = 2 and 5 the necessary derivatives are
  !             approximated numerically by SDRIV3, and we only ask the
  !             user to tell SDRIV3 the value of ML and MU if the system
  !             is banded.  For the cases MITER = 1 and 4 the user must
  !             derive these partials algebraically and encode them in
  !             subroutine JACOBN.  By computing these derivatives the
  !             user can often save 20-30 per cent of the computing time.
  !             Usually, however, the accuracy is not much affected and
  !             most users will probably forego this option.  The optional
  !             user-written subroutine JACOBN has the form:
  !                   SUBROUTINE JACOBN (N, T, Y, DFDY, MATDIM, ML, MU)
  !                   REAL Y(*), DFDY(MATDIM,*)
  !                     .
  !                     .
  !                     Calculate values of DFDY
  !                     .
  !                     .
  !                   END (Sample)
  !             Here Y is a vector of length at least N.  The actual
  !             length of Y is determined by the user's declaration in the
  !             program which calls SDRIV3.  Thus the dimensioning of Y in
  !             JACOBN, while required by FORTRAN convention, does not
  !             actually allocate any storage.  When this subroutine is
  !             called, the first N components of Y are intermediate
  !             approximations to the solution components.  The user
  !             should not alter these values.  If the system is not
  !             banded (MITER=1), the partials of the I-th equation with
  !             respect to the J-th dependent function are to be stored in
  !             DFDY(I,J).  Thus partials of the I-th equation are stored
  !             in the I-th row of DFDY.  If the system is banded
  !             (MITER=4), then the partials of the I-th equation with
  !             respect to Y(J) are to be stored in DFDY(K,J), where
  !             K=I-J+MU+1 .  Normally a return from JACOBN passes control
  !             back to SDRIV3.  However, if the user would like to abort
  !             the calculation, i.e., return control to the program which
  !             calls SDRIV3, he should set N to zero.  SDRIV3 will signal
  !             this by returning a value of NSTATE equal to +8(-8).
  !             Altering the value of N in JACOBN has no effect on the
  !             value of N in the call sequence of SDRIV3.
  !
  !    FA     = A subroutine supplied by the user if IMPL is not zero, and
  !             MITER is not 3.  If so, the name must be declared EXTERNAL
  !             in the user's calling program.  This subroutine computes
  !             the array A, where A*dY(I)/dT = F(Y(I),T).
  !             There are three cases:
  !
  !               IMPL=1.
  !               Subroutine FA is of the form:
  !                   SUBROUTINE FA (N, T, Y, A, MATDIM, ML, MU, NDE)
  !                   REAL Y(*), A(MATDIM,*)
  !                     .
  !                     .
  !                     Calculate ALL values of A
  !                     .
  !                     .
  !                   END (Sample)
  !               In this case A is assumed to be a nonsingular matrix,
  !               with the same structure as DFDY (see JACOBN description
  !               above).  Programming considerations prevent complete
  !               generality.  If MITER is 1 or 2, A is assumed to be full
  !               and the user must compute and store all values of
  !               A(I,J), I,J=1, ... ,N.  If MITER is 4 or 5, A is assumed
  !               to be banded with lower and upper half bandwidth ML and
  !               MU.  The left hand side of the I-th equation is a linear
  !               combination of dY(I-ML)/dT, dY(I-ML+1)/dT, ... ,
  !               dY(I)/dT, ..., dY(I+MU-1)/dT, dY(I+MU)/dT.  Thus in the
  !               I-th equation, the coefficient of dY(J)/dT is to be
  !               stored in A(K,J), where K=I-J+MU+1.
  !               NOTE: The array A will be altered between calls to FA.
  !
  !               IMPL=2.
  !               Subroutine FA is of the form:
  !                   SUBROUTINE FA (N, T, Y, A, MATDIM, ML, MU, NDE)
  !                   REAL Y(*), A(*)
  !                     .
  !                     .
  !                     Calculate non-zero values of A(1),...,A(NDE)
  !                     .
  !                     .
  !                   END (Sample)
  !               In this case it is assumed that the system is ordered by
  !               the user so that the differential equations appear
  !               first, and the algebraic equations appear last.  The
  !               algebraic equations must be written in the form:
  !               0 = F(Y(I),T).  When using this option it is up to the
  !               user to provide initial values for the Y(I) that satisfy
  !               the algebraic equations as well as possible.  It is
  !               further assumed that A is a vector of length NDE.  All
  !               of the components of A, which may depend on T, Y(I),
  !               etc., must be set by the user to non-zero values.
  !
  !               IMPL=3.
  !               Subroutine FA is of the form:
  !                   SUBROUTINE FA (N, T, Y, A, MATDIM, ML, MU, NDE)
  !                   REAL Y(*), A(MATDIM,*)
  !                     .
  !                     .
  !                     Calculate ALL values of A
  !                     .
  !                     .
  !                   END (Sample)
  !               In this case A is assumed to be a nonsingular NDE by NDE
  !               matrix with the same structure as DFDY (see JACOBN
  !               description above).  Programming considerations prevent
  !               complete generality.  If MITER is 1 or 2, A is assumed
  !               to be full and the user must compute and store all
  !               values of A(I,J), I,J=1, ... ,NDE.  If MITER is 4 or 5,
  !               A is assumed to be banded with lower and upper half
  !               bandwidths ML and MU.  The left hand side of the I-th
  !               equation is a linear combination of dY(I-ML)/dT,
  !               dY(I-ML+1)/dT, ..., dY(I)/dT, ..., dY(I+MU-1)/dT,
  !               dY(I+MU)/dT.  Thus in the I-th equation, the coefficient
  !               of dY(J)/dT is to be stored in A(K,J), where K=I-J+MU+1.
  !               It is assumed that the system is ordered by the user so
  !               that the differential equations appear first, and the
  !               algebraic equations appear last.  The algebraic
  !               equations must be written in the form 0 = F(Y(I),T).
  !               When using this option it is up to the user to provide
  !               initial values for the Y(I) that satisfy the algebraic
  !               equations as well as possible.
  !               NOTE: For IMPL = 3, the array A will be altered between
  !               calls to FA.
  !             Here Y is a vector of length at least N.  The actual
  !             length of Y is determined by the user's declaration in the
  !             program which calls SDRIV3.  Thus the dimensioning of Y in
  !             FA, while required by FORTRAN convention, does not
  !             actually allocate any storage.  When this subroutine is
  !             called, the first N components of Y are intermediate
  !             approximations to the solution components.  The user
  !             should not alter these values.  FA is always called
  !             immediately after calling F, with the same values of T
  !             and Y.  Normally a return from FA passes control back to
  !             SDRIV3.  However, if the user would like to abort the
  !             calculation, i.e., return control to the program which
  !             calls SDRIV3, he should set N to zero.  SDRIV3 will signal
  !             this by returning a value of NSTATE equal to +9(-9).
  !             Altering the value of N in FA has no effect on the value
  !             of N in the call sequence of SDRIV3.
  !
  !    NDE    = (Input) The number of differential equations.  This is
  !             required only for IMPL = 2 or 3, with NDE .LT. N.
  !
  !    MXSTEP = (Input) The maximum number of internal steps allowed on
  !             one call to SDRIV3.
  !
  !    G      = A real FORTRAN function supplied by the user
  !             if NROOT is not 0.  In this case, the name must be
  !             declared EXTERNAL in the user's calling program.  G is
  !             repeatedly called with different values of IROOT to obtain
  !             the value of each of the NROOT equations for which a root
  !             is desired.  G is of the form:
  !                   REAL FUNCTION G (N, T, Y, IROOT)
  !                   REAL Y(*)
  !                   GO TO (10, ...), IROOT
  !              10   G = ...
  !                     .
  !                     .
  !                   END (Sample)
  !             Here, Y is a vector of length at least N, whose first N
  !             components are the solution components at the point T.
  !             The user should not alter these values.  The actual length
  !             of Y is determined by the user's declaration in the
  !             program which calls SDRIV3.  Thus the dimensioning of Y in
  !             G, while required by FORTRAN convention, does not actually
  !             allocate any storage.  Normally a return from G passes
  !             control back to  SDRIV3.  However, if the user would like
  !             to abort the calculation, i.e., return control to the
  !             program which calls SDRIV3, he should set N to zero.
  !             SDRIV3 will signal this by returning a value of NSTATE
  !             equal to +7(-7).  In this case, the index of the equation
  !             being evaluated is stored in the sixth element of IWORK.
  !             Altering the value of N in G has no effect on the value of
  !             N in the call sequence of SDRIV3.
  !
  !    USERS  = A subroutine supplied by the user, if MITER is 3.
  !             If this is the case, the name must be declared EXTERNAL in
  !             the user's calling program.  The routine USERS is called
  !             by SDRIV3 when certain linear systems must be solved.  The
  !             user may choose any method to form, store and solve these
  !             systems in order to obtain the solution result that is
  !             returned to SDRIV3.  In particular, this allows sparse
  !             matrix methods to be used.  The call sequence for this
  !             routine is:
  !
  !                SUBROUTINE USERS (Y, YH, YWT, SAVE1, SAVE2, T, H, EL,
  !               8                  IMPL, N, NDE, IFLAG)
  !                REAL Y(*), YH(*), YWT(*), SAVE1(*),
  !               8     SAVE2(*), T, H, EL
  !
  !             The input variable IFLAG indicates what action is to be
  !             taken.  Subroutine USERS should perform the following
  !             operations, depending on the value of IFLAG and IMPL.
  !
  !               IFLAG = 0
  !                 IMPL = 0.  USERS is not called.
  !                 IMPL = 1, 2 or 3.  Solve the system A*X = SAVE2,
  !                   returning the result in SAVE2.  The array SAVE1 can
  !                   be used as a work array.  For IMPL = 1, there are N
  !                   components to the system, and for IMPL = 2 or 3,
  !                   there are NDE components to the system.
  !
  !               IFLAG = 1
  !                 IMPL = 0.  Compute, decompose and store the matrix
  !                   (I - H*EL*J), where I is the identity matrix and J
  !                   is the Jacobian matrix of the right hand side.  The
  !                   array SAVE1 can be used as a work array.
  !                 IMPL = 1, 2 or 3. Compute, decompose and store the
  !                   matrix (A - H*EL*J).  The array SAVE1 can be used as
  !                   a work array.
  !
  !               IFLAG = 2
  !                 IMPL = 0.   Solve the system
  !                     (I - H*EL*J)*X = H*SAVE2 - YH - SAVE1,
  !                   returning the result in SAVE2.
  !                 IMPL = 1, 2 or 3.  Solve the system
  !                   (A - H*EL*J)*X = H*SAVE2 - A*(YH + SAVE1)
  !                   returning the result in SAVE2.
  !                 The array SAVE1 should not be altered.
  !             If IFLAG is 0 and IMPL is 1 or 2 and the matrix A is
  !             singular, or if IFLAG is 1 and one of the matrices
  !             (I - H*EL*J), (A - H*EL*J) is singular, the INTEGER
  !             variable IFLAG is to be set to -1 before RETURNing.
  !             Normally a return from USERS passes control back to
  !             SDRIV3.  However, if the user would like to abort the
  !             calculation, i.e., return control to the program which
  !             calls SDRIV3, he should set N to zero.  SDRIV3 will signal
  !             this by returning a value of NSTATE equal to +10(-10).
  !             Altering the value of N in USERS has no effect on the
  !             value of N in the call sequence of SDRIV3.
  !
  !    IERFLG = An error flag.  The error number associated with a
  !             diagnostic message (see Section III-A below) is the same
  !             as the corresponding value of IERFLG.  The meaning of
  !             IERFLG:
  !               0  The routine completed successfully. (No message is
  !                  issued.)
  !               3  (Warning) The number of steps required to reach TOUT
  !                  exceeds MXSTEP.
  !               4  (Warning) The value of EPS is too small.
  !              11  (Warning) For NTASK = 2 or 3, T is beyond TOUT.
  !                  The solution was obtained by interpolation.
  !              15  (Warning) The integration step size is below the
  !                  roundoff level of T.  (The program issues this
  !                  message as a warning but does not return control to
  !                  the user.)
  !              22  (Recoverable) N is not positive.
  !              23  (Recoverable) MINT is less than 1 or greater than 3 .
  !              24  (Recoverable) MITER is less than 0 or greater than
  !                  5 .
  !              25  (Recoverable) IMPL is less than 0 or greater than 3 .
  !              26  (Recoverable) The value of NSTATE is less than 1 or
  !                  greater than 12 .
  !              27  (Recoverable) EPS is less than zero.
  !              28  (Recoverable) MXORD is not positive.
  !              29  (Recoverable) For MINT = 3, either MITER = 0 or 3, or
  !                  IMPL = 0 .
  !              30  (Recoverable) For MITER = 0, IMPL is not 0 .
  !              31  (Recoverable) For MINT = 1, IMPL is 2 or 3 .
  !              32  (Recoverable) Insufficient storage has been allocated
  !                  for the WORK array.
  !              33  (Recoverable) Insufficient storage has been allocated
  !                  for the IWORK array.
  !              41  (Recoverable) The integration step size has gone
  !                  to zero.
  !              42  (Recoverable) The integration step size has been
  !                  reduced about 50 times without advancing the
  !                  solution.  The problem setup may not be correct.
  !              43  (Recoverable)  For IMPL greater than 0, the matrix A
  !                  is singular.
  !             999  (Fatal) The value of NSTATE is 12 .
  !
  !  III.  OTHER COMMUNICATION TO THE USER  ..............................
  !
  !    A. The solver communicates to the user through the parameters
  !       above.  In addition it writes diagnostic messages through the
  !       standard error handling program XERMSG.  A complete description
  !       of XERMSG is given in "Guide to the SLATEC Common Mathematical
  !       Library" by Kirby W. Fong et al..  At installations which do not
  !       have this error handling package the short but serviceable
  !       routine, XERMSG, available with this package, can be used.  That
  !       program uses the file named OUTPUT to transmit messages.
  !
  !    B. The first three elements of WORK and the first five elements of
  !       IWORK will contain the following statistical data:
  !         AVGH     The average step size used.
  !         HUSED    The step size last used (successfully).
  !         AVGORD   The average order used.
  !         IMXERR   The index of the element of the solution vector that
  !                  contributed most to the last error test.
  !         NQUSED   The order last used (successfully).
  !         NSTEP    The number of steps taken since last initialization.
  !         NFE      The number of evaluations of the right hand side.
  !         NJE      The number of evaluations of the Jacobian matrix.
  !
  !  IV.  REMARKS  .......................................................
  !
  !    A. Other routines used:
  !         SDNTP, SDZRO, SDSTP, SDNTL, SDPST, SDCOR, SDCST,
  !         SDPSC, and SDSCL;
  !         SGEFA, SGESL, SGBFA, SGBSL, and SNRM2 (from LINPACK)
  !         R1MACH (from the Bell Laboratories Machine Constants Package)
  !         XERMSG (from the SLATEC Common Math Library)
  !       The last seven routines above, not having been written by the
  !       present authors, are not explicitly part of this package.
  !
  !    B. On any return from SDRIV3 all information necessary to continue
  !       the calculation is contained in the call sequence parameters,
  !       including the work arrays.  Thus it is possible to suspend one
  !       problem, integrate another, and then return to the first.
  !
  !    C. If this package is to be used in an overlay situation, the user
  !       must declare in the primary overlay the variables in the call
  !       sequence to SDRIV3.
  !
  !    D. Changing parameters during an integration.
  !       The value of NROOT, EPS, EWT, IERROR, MINT, MITER, or HMAX may
  !       be altered by the user between calls to SDRIV3.  For example, if
  !       too much accuracy has been requested (the program returns with
  !       NSTATE = 4 and an increased value of EPS) the user may wish to
  !       increase EPS further.  In general, prudence is necessary when
  !       making changes in parameters since such changes are not
  !       implemented until the next integration step, which is not
  !       necessarily the next call to SDRIV3.  This can happen if the
  !       program has already integrated to a point which is beyond the
  !       new point TOUT.
  !
  !    E. As the price for complete control of matrix algebra, the SDRIV3
  !       USERS option puts all responsibility for Jacobian matrix
  !       evaluation on the user.  It is often useful to approximate
  !       numerically all or part of the Jacobian matrix.  However this
  !       must be done carefully.  The FORTRAN sequence below illustrates
  !       the method we recommend.  It can be inserted directly into
  !       subroutine USERS to approximate Jacobian elements in rows I1
  !       to I2 and columns J1 to J2.
  !              REAL DFDY(N,N), EPSJ, H, R, R1MACH,
  !             8     SAVE1(N), SAVE2(N), T, UROUND, Y(N), YJ, YWT(N)
  !              UROUND = R1MACH(4)
  !              EPSJ = SQRT(UROUND)
  !              DO 30 J = J1,J2
  !                R = EPSJ*MAX(ABS(YWT(J)), ABS(Y(J)))
  !                IF (R .EQ. 0.E0) R = YWT(J)
  !                YJ = Y(J)
  !                Y(J) = Y(J) + R
  !                CALL F (N, T, Y, SAVE1)
  !                IF (N .EQ. 0) RETURN
  !                Y(J) = YJ
  !                DO 20 I = I1,I2
  !         20       DFDY(I,J) = (SAVE1(I) - SAVE2(I))/R
  !         30     CONTINUE
  !       Many problems give rise to structured sparse Jacobians, e.g.,
  !       block banded.  It is possible to approximate them with fewer
  !       function evaluations than the above procedure uses; see Curtis,
  !       Powell and Reid, J. Inst. Maths Applics, (1974), Vol. 13,
  !       pp. 117-119.
  !
  !    F. When any of the routines JACOBN, FA, G, or USERS, is not
  !       required, difficulties associated with unsatisfied externals can
  !       be avoided by using the name of the routine which calculates the
  !       right hand side of the differential equations in place of the
  !       corresponding name in the call sequence of SDRIV3.
  !
  !***REFERENCES  C. W. Gear, Numerical Initial Value Problems in
  !                 Ordinary Differential Equations, Prentice-Hall, 1971.
  !***ROUTINES CALLED  R1MACH, SDNTP, SDSTP, SDZRO, SGBFA, SGBSL, SGEFA,
  !                    SGESL, SNRM2, XERMSG
  !***REVISION HISTORY  (YYMMDD)
  !   790601  DATE WRITTEN
  !   900329  Initial submission to SLATEC.
  !***END PROLOGUE  SDRIV3
  EXTERNAL F, JACOBN, FA, G, USERS
  REAL ae, big, Eps, Ewt(*), G, glast, gnow, h, Hmax, hsign, &
    hused, NROUND, re, R1MACH, size, SNRM2, sum, T, tlast, &
    Tout, troot, uround, Work(*), Y(*)
  INTEGER i, ia, IAVGH, IAVGRD, ICNVRG, idfdy, IEL, Ierflg, Ierror, &
    ifac, iflag, ignow, IH, IHMAX, IHOLD, IHSIGN, IHUSED, &
    IJROOT, IJSTPL, IJTASK, IMNT, IMNTLD, Impl, IMTR, IMTRLD, &
    IMTRSV, imxerr, IMXORD, IMXRDS, INDMXR, INDPRT, INDPVT, &
    INDTRT, INFE, info, INJE, INQ, INQUSE, INROOT, INRTLD, &
    INSTEP, INWAIT, IRC, IRMAX, iroot, IMACH1, IMACH4, isave1, &
    isave2, IT, ITOUT, ITQ, ITREND, itroot, Iwork(*), IYH, &
    iywt, j, jstate, jtroot, lenchk, Leniw, Lenw, liwchk, &
    matdim, maxord, Mint, Miter, Ml, Mu, Mxord, Mxstep, N, &
    Nde, ndecom, npar, Nroot, Nstate, nstepl, Ntask
  LOGICAL convrg
  CHARACTER intgr1*8, intgr2*8, rl1*16, rl2*16
  PARAMETER (NROUND=20.E0)
  PARAMETER (IAVGH=1,IHUSED=2,IAVGRD=3,IEL=4,IH=160,IHMAX=161,IHOLD=162,&
    IHSIGN=163,IRC=164,IRMAX=165,IT=166,ITOUT=167,ITQ=168,&
    ITREND=204,IMACH1=205,IMACH4=206,IYH=251,INDMXR=1,INQUSE=2,&
    INSTEP=3,INFE=4,INJE=5,INROOT=6,ICNVRG=7,IJROOT=8,IJTASK=9,&
    IMNTLD=10,IMTRLD=11,INQ=12,INRTLD=13,INDTRT=14,INWAIT=15,&
    IMNT=16,IMTRSV=17,IMTR=18,IMXRDS=19,IMXORD=20,INDPRT=21,&
    IJSTPL=22,INDPVT=51)
  !***FIRST EXECUTABLE STATEMENT  SDRIV3
  IF ( Nstate==12 ) THEN
    Ierflg = 999
    CALL XERMSG('SLATEC','SDRIV3',&
      'Illegal input.  The value of NSTATE is 12 .',Ierflg,2)
    RETURN
  ELSEIF ( Nstate<1.OR.Nstate>12 ) THEN
    WRITE (intgr1,'(I8)') Nstate
    Ierflg = 26
    CALL XERMSG('SLATEC','SDRIV3',&
      'Illegal input.  Improper value for NSTATE(= '//intgr1//&
      ').',Ierflg,1)
    Nstate = 12
    RETURN
  ENDIF
  npar = N
  IF ( Eps<0.E0 ) THEN
    WRITE (rl1,'(E16.8)') Eps
    Ierflg = 27
    CALL XERMSG('SLATEC','SDRIV3','Illegal input.  EPS, '//rl1//&
      ', is negative.',Ierflg,1)
    Nstate = 12
    RETURN
  ENDIF
  IF ( N<=0 ) THEN
    WRITE (intgr1,'(I8)') N
    Ierflg = 22
    CALL XERMSG('SLATEC','SDRIV3','Illegal input.  Number of equations, '//&
      intgr1//', is not positive.',Ierflg,1)
    Nstate = 12
    RETURN
  ENDIF
  IF ( Mxord<=0 ) THEN
    WRITE (intgr1,'(I8)') Mxord
    Ierflg = 28
    CALL XERMSG('SLATEC','SDRIV3','Illegal input.  Maximum order, '//&
      intgr1//', is not positive.',Ierflg,1)
    Nstate = 12
    RETURN
  ENDIF
  IF ( Mint<1.OR.Mint>3 ) THEN
    WRITE (intgr1,'(I8)') Mint
    Ierflg = 23
    CALL XERMSG('SLATEC','SDRIV3',&
      'Illegal input.  Improper value for the integration method flag, '//intgr1//' .',Ierflg,1)
    Nstate = 12
    RETURN
  ELSEIF ( Miter<0.OR.Miter>5 ) THEN
    WRITE (intgr1,'(I8)') Miter
    Ierflg = 24
    CALL XERMSG('SLATEC','SDRIV3',&
      'Illegal input.  Improper value for MITER(= '//intgr1//').',&
      Ierflg,1)
    Nstate = 12
    RETURN
  ELSEIF ( Impl<0.OR.Impl>3 ) THEN
    WRITE (intgr1,'(I8)') Impl
    Ierflg = 25
    CALL XERMSG('SLATEC','SDRIV3',&
      'Illegal input.  Improper value for IMPL(= '//intgr1//').',&
      Ierflg,1)
    Nstate = 12
    RETURN
  ELSEIF ( Mint==3.AND.(Miter==0.OR.Miter==3.OR.Impl/=0) ) THEN
    WRITE (intgr1,'(I8)') Miter
    WRITE (intgr2,'(I8)') Impl
    Ierflg = 29
    CALL XERMSG('SLATEC','SDRIV3',&
      'Illegal input.  For MINT = 3, the value of MITER, '//&
      intgr1//', and/or IMPL, '//intgr2//', is not allowed.',&
      Ierflg,1)
    Nstate = 12
    RETURN
  ELSEIF ( (Impl>=1.AND.Impl<=3).AND.Miter==0 ) THEN
    WRITE (intgr1,'(I8)') Impl
    Ierflg = 30
    CALL XERMSG('SLATEC','SDRIV3',&
      'Illegal input.  For MITER = 0, the value of IMPL, '//&
      intgr1//', is not allowed.',Ierflg,1)
    Nstate = 12
    RETURN
  ELSEIF ( (Impl==2.OR.Impl==3).AND.Mint==1 ) THEN
    WRITE (intgr1,'(I8)') Impl
    Ierflg = 31
    CALL XERMSG('SLATEC','SDRIV3',&
      'Illegal input.  For MINT = 1, the value of IMPL, '//&
      intgr1//', is not allowed.',Ierflg,1)
    Nstate = 12
    RETURN
  ENDIF
  IF ( Miter==0.OR.Miter==3 ) THEN
    liwchk = INDPVT - 1
  ELSEIF ( Miter==1.OR.Miter==2.OR.Miter==4.OR.Miter==5 ) THEN
    liwchk = INDPVT + N - 1
  ENDIF
  IF ( Leniw<liwchk ) THEN
    WRITE (intgr1,'(I8)') liwchk
    Ierflg = 33
    CALL XERMSG('SLATEC','SDRIV3','Illegal input.  Insufficient storage allocated&
      & for the IWORK array.  Based on the value of the input parameters involved,&
      & the required storage is '//intgr1//' .',Ierflg,&
      1)
    Nstate = 12
    RETURN
  ENDIF
  !                                                Allocate the WORK array
  !                                         IYH is the index of YH in WORK
  IF ( Mint==1.OR.Mint==3 ) THEN
    maxord = MIN(Mxord,12)
  ELSEIF ( Mint==2 ) THEN
    maxord = MIN(Mxord,5)
  ENDIF
  idfdy = IYH + (maxord+1)*N
  !                                             IDFDY is the index of DFDY
  !
  IF ( Miter==0.OR.Miter==3 ) THEN
    iywt = idfdy
  ELSEIF ( Miter==1.OR.Miter==2 ) THEN
    iywt = idfdy + N*N
  ELSEIF ( Miter==4.OR.Miter==5 ) THEN
    iywt = idfdy + (2*Ml+Mu+1)*N
  ENDIF
  !                                               IYWT is the index of YWT
  isave1 = iywt + N
  !                                           ISAVE1 is the index of SAVE1
  isave2 = isave1 + N
  !                                           ISAVE2 is the index of SAVE2
  ignow = isave2 + N
  !                                             IGNOW is the index of GNOW
  itroot = ignow + Nroot
  !                                           ITROOT is the index of TROOT
  ifac = itroot + Nroot
  !                                               IFAC is the index of FAC
  IF ( Miter==2.OR.Miter==5.OR.Mint==3 ) THEN
    ia = ifac + N
  ELSE
    ia = ifac
  ENDIF
  !                                                   IA is the index of A
  IF ( Impl==0.OR.Miter==3 ) THEN
    lenchk = ia - 1
  ELSEIF ( Impl==1.AND.(Miter==1.OR.Miter==2) ) THEN
    lenchk = ia - 1 + N*N
  ELSEIF ( Impl==1.AND.(Miter==4.OR.Miter==5) ) THEN
    lenchk = ia - 1 + (2*Ml+Mu+1)*N
  ELSEIF ( Impl==2.AND.Miter/=3 ) THEN
    lenchk = ia - 1 + N
  ELSEIF ( Impl==3.AND.(Miter==1.OR.Miter==2) ) THEN
    lenchk = ia - 1 + N*Nde
  ELSEIF ( Impl==3.AND.(Miter==4.OR.Miter==5) ) THEN
    lenchk = ia - 1 + (2*Ml+Mu+1)*Nde
  ENDIF
  IF ( Lenw<lenchk ) THEN
    WRITE (intgr1,'(I8)') lenchk
    Ierflg = 32
    CALL XERMSG('SLATEC','SDRIV3','Illegal input.  Insufficient storage&
      & allocated for the WORK array.  Based on the value of the input&
      & parameters involved, the required storage is '//intgr1//' .',Ierflg,1)
    Nstate = 12
    RETURN
  ENDIF
  IF ( Miter==0.OR.Miter==3 ) THEN
    matdim = 1
  ELSEIF ( Miter==1.OR.Miter==2 ) THEN
    matdim = N
  ELSEIF ( Miter==4.OR.Miter==5 ) THEN
    matdim = 2*Ml + Mu + 1
  ENDIF
  IF ( Impl==0.OR.Impl==1 ) THEN
    ndecom = N
  ELSEIF ( Impl==2.OR.Impl==3 ) THEN
    ndecom = Nde
  ENDIF
  IF ( Nstate==1 ) THEN
    !                                                  Initialize parameters
    IF ( Mint==1.OR.Mint==3 ) THEN
      Iwork(IMXORD) = MIN(Mxord,12)
    ELSEIF ( Mint==2 ) THEN
      Iwork(IMXORD) = MIN(Mxord,5)
    ENDIF
    Iwork(IMXRDS) = Mxord
    IF ( Mint==1.OR.Mint==2 ) THEN
      Iwork(IMNT) = Mint
      Iwork(IMTR) = Miter
      Iwork(IMNTLD) = Mint
      Iwork(IMTRLD) = Miter
    ELSEIF ( Mint==3 ) THEN
      Iwork(IMNT) = 1
      Iwork(IMTR) = 0
      Iwork(IMNTLD) = Iwork(IMNT)
      Iwork(IMTRLD) = Iwork(IMTR)
      Iwork(IMTRSV) = Miter
    ENDIF
    Work(IHMAX) = Hmax
    uround = R1MACH(4)
    Work(IMACH4) = uround
    Work(IMACH1) = R1MACH(1)
    IF ( Nroot/=0 ) THEN
      re = uround
      ae = Work(IMACH1)
    ENDIF
    h = (Tout-T)*(1.E0-4.E0*uround)
    h = SIGN(MIN(ABS(h),Hmax),h)
    Work(IH) = h
    hsign = SIGN(1.E0,h)
    Work(IHSIGN) = hsign
    Iwork(IJTASK) = 0
    Work(IAVGH) = 0.E0
    Work(IHUSED) = 0.E0
    Work(IAVGRD) = 0.E0
    Iwork(INDMXR) = 0
    Iwork(INQUSE) = 0
    Iwork(INSTEP) = 0
    Iwork(IJSTPL) = 0
    Iwork(INFE) = 0
    Iwork(INJE) = 0
    Iwork(INROOT) = 0
    Work(IT) = T
    Iwork(ICNVRG) = 0
    Iwork(INDPRT) = 0
    !                                                 Set initial conditions
    DO i = 1, N
      Work(i+IYH-1) = Y(i)
    ENDDO
    IF ( T==Tout ) RETURN
    GOTO 100
  ELSE
    uround = Work(IMACH4)
    IF ( Nroot/=0 ) THEN
      re = uround
      ae = Work(IMACH1)
    ENDIF
  ENDIF
  !                                             On a continuation, check
  !                                             that output points have
  !                                             been or will be overtaken.
  IF ( Iwork(ICNVRG)==1 ) THEN
    convrg = .TRUE.
  ELSE
    convrg = .FALSE.
  ENDIF
  T = Work(IT)
  h = Work(IH)
  hsign = Work(IHSIGN)
  IF ( Iwork(IJTASK)/=0 ) THEN
    !
    !                                   IWORK(IJROOT) flags unreported
    !                                   roots, and is set to the value of
    !                                   NTASK when a root was last selected.
    !                                   It is set to zero when all roots
    !                                   have been reported.  IWORK(INROOT)
    !                                   contains the index and WORK(ITOUT)
    !                                   contains the value of the root last
    !                                   selected to be reported.
    !                                   IWORK(INRTLD) contains the value of
    !                                   NROOT and IWORK(INDTRT) contains
    !                                   the value of ITROOT when the array
    !                                   of roots was last calculated.
    IF ( Nroot/=0 ) THEN
      IF ( Iwork(IJROOT)>0 ) THEN
        !                                      TOUT has just been reported.
        !                                      If TROOT .LE. TOUT, report TROOT.
        IF ( Nstate==5 ) THEN
          troot = T
          iroot = 0
          DO i = 1, Iwork(INRTLD)
            jtroot = i + Iwork(INDTRT) - 1
            IF ( Work(jtroot)*hsign<=troot*hsign ) THEN
              !
              !                                              Check for multiple roots.
              !
              IF ( Work(jtroot)==Work(ITOUT).AND.i>Iwork(INROOT) ) THEN
                iroot = i
                troot = Work(jtroot)
                EXIT
              ENDIF
              IF ( Work(jtroot)*hsign>Work(ITOUT)*hsign ) THEN
                iroot = i
                troot = Work(jtroot)
              ENDIF
            ENDIF
          ENDDO
          Iwork(INROOT) = iroot
          Work(ITOUT) = troot
          Iwork(IJROOT) = Ntask
          IF ( Ntask==1 ) THEN
            IF ( iroot==0 ) THEN
              Iwork(IJROOT) = 0
            ELSEIF ( Tout*hsign>=troot*hsign ) THEN
              CALL SDNTP(h,0,N,Iwork(INQ),T,troot,Work(IYH),Y)
              Nstate = 5
              T = troot
              Ierflg = 0
              GOTO 500
            ENDIF
          ELSEIF ( Ntask==2.OR.Ntask==3 ) THEN
            !
            !                                     If there are no more roots, or the
            !                                     user has altered TOUT to be less
            !                                     than a root, set IJROOT to zero.
            !
            IF ( iroot==0.OR.(Tout*hsign<troot*hsign) ) THEN
              Iwork(IJROOT) = 0
            ELSE
              CALL SDNTP(h,0,N,Iwork(INQ),T,troot,Work(IYH),Y)
              Nstate = 5
              Ierflg = 0
              T = troot
              GOTO 500
            ENDIF
          ENDIF
        ELSEIF ( Tout*hsign>=Work(ITOUT)*hsign ) THEN
          troot = Work(ITOUT)
          CALL SDNTP(h,0,N,Iwork(INQ),T,troot,Work(IYH),Y)
          T = troot
          Nstate = 5
          Ierflg = 0
          GOTO 500
          !                                         A root has just been reported.
          !                                         Select the next root.
        ENDIF
      ENDIF
    ENDIF
    !
    IF ( Ntask==1 ) THEN
      Nstate = 2
      IF ( T*hsign>=Tout*hsign ) THEN
        CALL SDNTP(h,0,N,Iwork(INQ),T,Tout,Work(IYH),Y)
        T = Tout
        Ierflg = 0
        GOTO 500
      ENDIF
    ELSEIF ( Ntask==2 ) THEN
      !                                                      Check if TOUT has
      !                                                      been reset .LT. T
      IF ( T*hsign>Tout*hsign ) THEN
        WRITE (rl1,'(E16.8)') T
        WRITE (rl2,'(E16.8)') Tout
        Ierflg = 11
        CALL XERMSG('SLATEC','SDRIV3',&
          'While integrating exactly to TOUT, T, '//rl1//&
          ', was beyond TOUT, '//rl2//' .  Solution obtained by interpolation.',Ierflg,0)
        Nstate = 11
        CALL SDNTP(h,0,N,Iwork(INQ),T,Tout,Work(IYH),Y)
        T = Tout
        GOTO 500
      ENDIF
      !                                   Determine if TOUT has been overtaken
      !
      IF ( ABS(Tout-T)<=NROUND*uround*MAX(ABS(T),ABS(Tout)) ) THEN
        T = Tout
        Nstate = 2
        Ierflg = 0
        GOTO 400
      ENDIF
      !                                             If there are no more roots
      !                                             to report, report T.
      IF ( Nstate==5 ) THEN
        Nstate = 2
        Ierflg = 0
        GOTO 400
      ENDIF
      Nstate = 2
      !                                                       See if TOUT will
      !                                                       be overtaken.
      IF ( (T+h)*hsign>Tout*hsign ) THEN
        h = Tout - T
        IF ( (T+h)*hsign>Tout*hsign ) h = h*(1.E0-4.E0*uround)
        Work(IH) = h
        IF ( h==0.E0 ) GOTO 600
        Iwork(IJTASK) = -1
      ENDIF
    ELSEIF ( Ntask==3 ) THEN
      Nstate = 2
      IF ( T*hsign>Tout*hsign ) THEN
        WRITE (rl1,'(E16.8)') T
        WRITE (rl2,'(E16.8)') Tout
        Ierflg = 11
        CALL XERMSG('SLATEC','SDRIV3',&
          'While integrating exactly to TOUT, T, '//rl1//&
          ', was beyond TOUT, '//rl2//' .  Solution obtained by interpolation.',Ierflg,0)
        Nstate = 11
        CALL SDNTP(h,0,N,Iwork(INQ),T,Tout,Work(IYH),Y)
        T = Tout
        GOTO 500
      ENDIF
      IF ( ABS(Tout-T)<=NROUND*uround*MAX(ABS(T),ABS(Tout)) ) THEN
        T = Tout
        Ierflg = 0
        GOTO 400
      ENDIF
      IF ( (T+h)*hsign>Tout*hsign ) THEN
        h = Tout - T
        IF ( (T+h)*hsign>Tout*hsign ) h = h*(1.E0-4.E0*uround)
        Work(IH) = h
        IF ( h==0.E0 ) GOTO 600
        Iwork(IJTASK) = -1
      ENDIF
    ENDIF
    !                         Implement changes in MINT, MITER, and/or HMAX.
    !
    IF ( (Mint/=Iwork(IMNTLD).OR.Miter/=Iwork(IMTRLD)).AND.Mint/=3.AND.&
      Iwork(IMNTLD)/=3 ) Iwork(IJTASK) = -1
    IF ( Hmax/=Work(IHMAX) ) THEN
      h = SIGN(MIN(ABS(h),Hmax),h)
      IF ( h/=Work(IH) ) THEN
        Iwork(IJTASK) = -1
        Work(IH) = h
      ENDIF
      Work(IHMAX) = Hmax
    ENDIF
  ENDIF
  !
  100  nstepl = Iwork(INSTEP)
  DO i = 1, N
    Y(i) = Work(i+IYH-1)
  ENDDO
  IF ( Nroot/=0 ) THEN
    DO i = 1, Nroot
      Work(i+ignow-1) = G(npar,T,Y,i)
      IF ( npar==0 ) THEN
        Iwork(INROOT) = i
        Nstate = 7
        RETURN
      ENDIF
    ENDDO
  ENDIF
  IF ( Ierror==1 ) THEN
    DO i = 1, N
      Work(i+iywt-1) = 1.E0
    ENDDO
    GOTO 300
  ELSEIF ( Ierror==5 ) THEN
    DO i = 1, N
      Work(i+iywt-1) = Ewt(i)
    ENDDO
    GOTO 300
  ENDIF
  !                                       Reset YWT array.  Looping point.
  200 CONTINUE
  IF ( Ierror==2 ) THEN
    DO i = 1, N
      IF ( Y(i)==0.E0 ) GOTO 250
      Work(i+iywt-1) = ABS(Y(i))
    ENDDO
    GOTO 300
    250 CONTINUE
    IF ( Iwork(IJTASK)==0 ) THEN
      CALL F(npar,T,Y,Work(isave2))
      IF ( npar==0 ) THEN
        Nstate = 6
        RETURN
      ENDIF
      Iwork(INFE) = Iwork(INFE) + 1
      IF ( Miter==3.AND.Impl/=0 ) THEN
        iflag = 0
        CALL USERS(Y,Work(IYH),Work(iywt),Work(isave1),Work(isave2),T,h,&
          Work(IEL),Impl,npar,ndecom,iflag)
        IF ( iflag==-1 ) GOTO 700
        IF ( npar==0 ) THEN
          Nstate = 10
          RETURN
        ENDIF
      ELSEIF ( Impl==1 ) THEN
        IF ( Miter==1.OR.Miter==2 ) THEN
          CALL FA(npar,T,Y,Work(ia),matdim,Ml,Mu,ndecom)
          IF ( npar==0 ) THEN
            Nstate = 9
            RETURN
          ENDIF
          CALL SGEFA(Work(ia),matdim,N,Iwork(INDPVT),info)
          IF ( info/=0 ) GOTO 700
          CALL SGESL(Work(ia),matdim,N,Iwork(INDPVT),Work(isave2),0)
        ELSEIF ( Miter==4.OR.Miter==5 ) THEN
          CALL FA(npar,T,Y,Work(ia+Ml),matdim,Ml,Mu,ndecom)
          IF ( npar==0 ) THEN
            Nstate = 9
            RETURN
          ENDIF
          CALL SGBFA(Work(ia),matdim,N,Ml,Mu,Iwork(INDPVT),info)
          IF ( info/=0 ) GOTO 700
          CALL SGBSL(Work(ia),matdim,N,Ml,Mu,Iwork(INDPVT),Work(isave2),0)
        ENDIF
      ELSEIF ( Impl==2 ) THEN
        CALL FA(npar,T,Y,Work(ia),matdim,Ml,Mu,ndecom)
        IF ( npar==0 ) THEN
          Nstate = 9
          RETURN
        ENDIF
        DO i = 1, ndecom
          IF ( Work(i+ia-1)==0.E0 ) GOTO 700
          Work(i+isave2-1) = Work(i+isave2-1)/Work(i+ia-1)
        ENDDO
      ELSEIF ( Impl==3 ) THEN
        IF ( Miter==1.OR.Miter==2 ) THEN
          CALL FA(npar,T,Y,Work(ia),matdim,Ml,Mu,ndecom)
          IF ( npar==0 ) THEN
            Nstate = 9
            RETURN
          ENDIF
          CALL SGEFA(Work(ia),matdim,Nde,Iwork(INDPVT),info)
          IF ( info/=0 ) GOTO 700
          CALL SGESL(Work(ia),matdim,Nde,Iwork(INDPVT),Work(isave2),0)
        ELSEIF ( Miter==4.OR.Miter==5 ) THEN
          CALL FA(npar,T,Y,Work(ia+Ml),matdim,Ml,Mu,ndecom)
          IF ( npar==0 ) THEN
            Nstate = 9
            RETURN
          ENDIF
          CALL SGBFA(Work(ia),matdim,Nde,Ml,Mu,Iwork(INDPVT),info)
          IF ( info/=0 ) GOTO 700
          CALL SGBSL(Work(ia),matdim,Nde,Ml,Mu,Iwork(INDPVT),Work(isave2),0)
        ENDIF
      ENDIF
    ENDIF
    DO j = i, N
      IF ( Y(j)/=0.E0 ) THEN
        Work(j+iywt-1) = ABS(Y(j))
      ELSEIF ( Iwork(IJTASK)==0 ) THEN
        Work(j+iywt-1) = ABS(h*Work(j+isave2-1))
      ELSE
        Work(j+iywt-1) = ABS(Work(j+IYH+N-1))
      ENDIF
      IF ( Work(j+iywt-1)==0.E0 ) Work(j+iywt-1) = uround
    ENDDO
  ELSEIF ( Ierror==3 ) THEN
    DO i = 1, N
      Work(i+iywt-1) = MAX(Ewt(1),ABS(Y(i)))
    ENDDO
  ELSEIF ( Ierror==4 ) THEN
    DO i = 1, N
      Work(i+iywt-1) = MAX(Ewt(i),ABS(Y(i)))
    ENDDO
  ENDIF
  !
  300 CONTINUE
  DO i = 1, N
    Work(i+isave2-1) = Y(i)/Work(i+iywt-1)
  ENDDO
  sum = SNRM2(N,Work(isave2),1)/SQRT(REAL(N))
  sum = MAX(1.E0,sum)
  IF ( Eps<sum*uround ) THEN
    Eps = sum*uround*(1.E0+10.E0*uround)
    WRITE (rl1,'(E16.8)') T
    WRITE (rl2,'(E16.8)') Eps
    Ierflg = 4
    CALL XERMSG('SLATEC','SDRIV3','At T, '//rl1//', the requested accuracy,&
      & EPS, was not obtainable with the machine precision.&
      & EPS has been increased to '//rl2//' .',Ierflg,0)
    Nstate = 4
    GOTO 400
  ENDIF
  IF ( ABS(h)>=uround*ABS(T) ) THEN
    Iwork(INDPRT) = 0
  ELSEIF ( Iwork(INDPRT)==0 ) THEN
    WRITE (rl1,'(E16.8)') T
    WRITE (rl2,'(E16.8)') h
    Ierflg = 15
    CALL XERMSG('SLATEC','SDRIV3','At T, '//rl1//', the step size, '//rl2//&
      ', is smaller than the roundoff level of T. This may occur if there is&
      & an abrupt change in the right hand side of the differential equations.',Ierflg,0)
    Iwork(INDPRT) = 1
  ENDIF
  IF ( Ntask/=2 ) THEN
    IF ( (Iwork(INSTEP)-nstepl)==Mxstep ) THEN
      WRITE (rl1,'(E16.8)') T
      WRITE (intgr1,'(I8)') Mxstep
      WRITE (rl2,'(E16.8)') Tout
      Ierflg = 3
      CALL XERMSG('SLATEC','SDRIV3','At T, '//rl1//', '//intgr1//&
        ' steps have been taken without reaching TOUT, '//&
        rl2//' .',Ierflg,0)
      Nstate = 3
      GOTO 400
    ENDIF
  ENDIF
  !
  !     CALL SDSTP (EPS, F, FA, HMAX, IMPL, IERROR, JACOBN, MATDIM,
  !    8            MAXORD, MINT, MITER, ML, MU, N, NDE, YWT, UROUND,
  !    8            USERS,  AVGH, AVGORD, H, HUSED, JTASK, MNTOLD, MTROLD,
  !    8            NFE, NJE, NQUSED, NSTEP, T, Y, YH,  A, CONVRG,
  !    8            DFDY, EL, FAC, HOLD, IPVT, JSTATE, JSTEPL, NQ, NWAIT,
  !    8            RC, RMAX, SAVE1, SAVE2, TQ, TREND, ISWFLG, MTRSV,
  !    8            MXRDSV)
  !
  CALL SDSTP(Eps,F,FA,Work(IHMAX),Impl,Ierror,JACOBN,matdim,Iwork(IMXORD),&
    Iwork(IMNT),Iwork(IMTR),Ml,Mu,npar,ndecom,Work(iywt),uround,&
    USERS,Work(IAVGH),Work(IAVGRD),Work(IH),hused,Iwork(IJTASK),&
    Iwork(IMNTLD),Iwork(IMTRLD),Iwork(INFE),Iwork(INJE),&
    Iwork(INQUSE),Iwork(INSTEP),Work(IT),Y,Work(IYH),Work(ia),&
    convrg,Work(idfdy),Work(IEL),Work(ifac),Work(IHOLD),&
    Iwork(INDPVT),jstate,Iwork(IJSTPL),Iwork(INQ),Iwork(INWAIT),&
    Work(IRC),Work(IRMAX),Work(isave1),Work(isave2),Work(ITQ),&
    Work(ITREND),Mint,Iwork(IMTRSV),Iwork(IMXRDS))
  T = Work(IT)
  h = Work(IH)
  IF ( convrg ) THEN
    Iwork(ICNVRG) = 1
  ELSE
    Iwork(ICNVRG) = 0
  ENDIF
  SELECT CASE (jstate)
    CASE (2)
      GOTO 600
    CASE (3)
      !
      WRITE (rl1,'(E16.8)') T
      Ierflg = 42
      CALL XERMSG('SLATEC','SDRIV3',&
        'At T, '//rl1//', the step size has been reduced about 50 '&
        //&
        'times without advancing the solution.  Often this occurs if the problem setup is incorrect.',Ierflg,1)
      Nstate = 12
      RETURN
    CASE (4,5)
      GOTO 700
    CASE (6,7,8,9,10)
      !
      Nstate = jstate
      RETURN
    CASE DEFAULT
      Iwork(IJTASK) = 1
      !                                 Determine if a root has been overtaken
      IF ( Nroot/=0 ) THEN
        iroot = 0
        DO i = 1, Nroot
          glast = Work(i+ignow-1)
          gnow = G(npar,T,Y,i)
          IF ( npar==0 ) THEN
            Iwork(INROOT) = i
            Nstate = 7
            RETURN
          ENDIF
          Work(i+ignow-1) = gnow
          IF ( glast*gnow>0.E0 ) THEN
            Work(i+itroot-1) = T + h
          ELSEIF ( gnow==0.E0 ) THEN
            Work(i+itroot-1) = T
            iroot = i
          ELSEIF ( glast==0.E0 ) THEN
            Work(i+itroot-1) = T + h
          ELSEIF ( ABS(hused)>=uround*ABS(T) ) THEN
            tlast = T - hused
            iroot = i
            troot = T
            CALL SDZRO(ae,G,h,npar,Iwork(INQ),iroot,re,T,Work(IYH),uround,&
              troot,tlast,gnow,glast,Y)
            DO j = 1, N
              Y(j) = Work(IYH+j-1)
            ENDDO
            IF ( npar==0 ) THEN
              Iwork(INROOT) = i
              Nstate = 7
              RETURN
            ENDIF
            Work(i+itroot-1) = troot
          ELSE
            Work(i+itroot-1) = T
            iroot = i
          ENDIF
        ENDDO
        IF ( iroot==0 ) THEN
          Iwork(IJROOT) = 0
          !                                                  Select the first root
        ELSE
          Iwork(IJROOT) = Ntask
          Iwork(INRTLD) = Nroot
          Iwork(INDTRT) = itroot
          troot = T + h
          DO i = 1, Nroot
            IF ( Work(i+itroot-1)*hsign<troot*hsign ) THEN
              troot = Work(i+itroot-1)
              iroot = i
            ENDIF
          ENDDO
          Iwork(INROOT) = iroot
          Work(ITOUT) = troot
          IF ( troot*hsign<=Tout*hsign ) THEN
            CALL SDNTP(h,0,N,Iwork(INQ),T,troot,Work(IYH),Y)
            Nstate = 5
            T = troot
            Ierflg = 0
            GOTO 500
          ENDIF
        ENDIF
      ENDIF
      !                               Test for NTASK condition to be satisfied
      Nstate = 2
      IF ( Ntask==1 ) THEN
        IF ( T*hsign<Tout*hsign ) GOTO 200
        CALL SDNTP(h,0,N,Iwork(INQ),T,Tout,Work(IYH),Y)
        T = Tout
        Ierflg = 0
        GOTO 500
        !                               TOUT is assumed to have been attained
        !                               exactly if T is within twenty roundoff
        !                               units of TOUT, relative to MAX(TOUT, T).
        !
      ELSEIF ( Ntask==2 ) THEN
        IF ( ABS(Tout-T)<=NROUND*uround*MAX(ABS(T),ABS(Tout)) ) THEN
          T = Tout
        ELSEIF ( (T+h)*hsign>Tout*hsign ) THEN
          h = Tout - T
          IF ( (T+h)*hsign>Tout*hsign ) h = h*(1.E0-4.E0*uround)
          Work(IH) = h
          IF ( h==0.E0 ) GOTO 600
          Iwork(IJTASK) = -1
        ENDIF
      ELSEIF ( Ntask==3 ) THEN
        IF ( ABS(Tout-T)<=NROUND*uround*MAX(ABS(T),ABS(Tout)) ) THEN
          T = Tout
        ELSE
          IF ( (T+h)*hsign>Tout*hsign ) THEN
            h = Tout - T
            IF ( (T+h)*hsign>Tout*hsign ) h = h*(1.E0-4.E0*uround)
            Work(IH) = h
            IF ( h==0.E0 ) GOTO 600
            Iwork(IJTASK) = -1
          ENDIF
          GOTO 200
        ENDIF
      ENDIF
      Ierflg = 0
  END SELECT
  !                                      All returns are made through this
  !                                      section.  IMXERR is determined.
  400 CONTINUE
  DO i = 1, N
    Y(i) = Work(i+IYH-1)
  ENDDO
  500 CONTINUE
  IF ( Iwork(IJTASK)==0 ) RETURN
  big = 0.E0
  imxerr = 1
  DO i = 1, N
    !                                            SIZE = ABS(ERROR(I)/YWT(I))
    size = ABS(Work(i+isave1-1)/Work(i+iywt-1))
    IF ( big<size ) THEN
      big = size
      imxerr = i
    ENDIF
  ENDDO
  Iwork(INDMXR) = imxerr
  Work(IHUSED) = hused
  RETURN
  !                                        Fatal errors are processed here
  !
  600  WRITE (rl1,'(E16.8)') T
  Ierflg = 41
  CALL XERMSG('SLATEC','SDRIV3',&
    'At T, '//rl1//', the attempted step size has gone to zero.  Often this occurs if the problem setup is incorrect.',&
    Ierflg,1)
  Nstate = 12
  RETURN
  !
  700  WRITE (rl1,'(E16.8)') T
  Ierflg = 43
  CALL XERMSG('SLATEC','SDRIV3',&
    'At T, '//rl1//', while solving A*YDOT = F, A is singular.',&
    Ierflg,1)
  Nstate = 12
END SUBROUTINE SDRIV3
