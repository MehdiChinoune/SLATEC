!** SDASSL
SUBROUTINE SDASSL(RES,Neq,T,Y,Yprime,Tout,Info,Rtol,Atol,Idid,Rwork,Lrw,&
    Iwork,Liw,JAC)
  !>
  !  This code solves a system of differential/algebraic
  !            equations of the form G(T,Y,YPRIME) = 0.
  !***
  ! **Library:**   SLATEC (DASSL)
  !***
  ! **Category:**  I1A2
  !***
  ! **Type:**      SINGLE PRECISION (SDASSL-S, DDASSL-D)
  !***
  ! **Keywords:**  BACKWARD DIFFERENTIATION FORMULAS, DASSL,
  !             DIFFERENTIAL/ALGEBRAIC, IMPLICIT DIFFERENTIAL SYSTEMS
  !***
  ! **Author:**  Petzold, Linda R., (LLNL)
  !             Computing and Mathematics Research Division
  !             Lawrence Livermore National Laboratory
  !             L - 316, P.O. Box 808,
  !             Livermore, CA.    94550
  !***
  ! **Description:**
  !
  !- Usage:
  !
  !      EXTERNAL RES, JAC
  !      INTEGER NEQ, INFO(N), IDID, LRW, LIW, IWORK(LIW), IPAR
  !      REAL T, Y(NEQ), YPRIME(NEQ), TOUT, RTOL, ATOL,
  !     *   RWORK(LRW), RPAR
  !
  !      CALL SDASSL (RES, NEQ, T, Y, YPRIME, TOUT, INFO, RTOL, ATOL,
  !     *   IDID, RWORK, LRW, IWORK, LIW, RPAR, IPAR, JAC)
  !
  !
  !- Arguments:
  !
  !  RES:EXT     This is a subroutine which you provide to define the
  !              differential/algebraic system.
  !
  !  NEQ:IN      This is the number of equations to be solved.
  !
  !  T:INOUT     This is the current value of the independent variable.
  !
  !  Y(*):INOUT  This array contains the solution components at T.
  !
  !  YPRIME(*):INOUT  This array contains the derivatives of the solution
  !              components at T.
  !
  !  TOUT:IN     This is a point at which a solution is desired.
  !
  !  INFO(N):IN  The basic task of the code is to solve the system from T
  !              to TOUT and return an answer at TOUT.  INFO is an integer
  !              array which is used to communicate exactly how you want
  !              this task to be carried out.  (See below for details.)
  !              N must be greater than or equal to 15.
  !
  !  RTOL,ATOL:INOUT  These quantities represent relative and absolute
  !              error tolerances which you provide to indicate how
  !              accurately you wish the solution to be computed.  You
  !              may choose them to be both scalars or else both vectors.
  !              Caution:  In Fortran 77, a scalar is not the same as an
  !                        array of length 1.  Some compilers may object
  !                        to using scalars for RTOL,ATOL.
  !
  !  IDID:OUT    This scalar quantity is an indicator reporting what the
  !              code did.  You must monitor this integer variable to
  !              decide  what action to take next.
  !
  !  RWORK:WORK  A real work array of length LRW which provides the
  !              code with needed storage space.
  !
  !  LRW:IN      The length of RWORK.  (See below for required length.)
  !
  !  IWORK:WORK  An integer work array of length LIW which provides the
  !              code with needed storage space.
  !
  !  LIW:IN      The length of IWORK.  (See below for required length.)
  !
  !  RPAR,IPAR:IN  These are real and integer parameter arrays which
  !              you can use for communication between your calling
  !              program and the RES subroutine (and the JAC subroutine)
  !
  !  JAC:EXT     This is the name of a subroutine which you may choose
  !              to provide for defining a matrix of partial derivatives
  !              described below.
  !
  !  Quantities which may be altered by SDASSL are:
  !     T, Y(*), YPRIME(*), INFO(1), RTOL, ATOL,
  !     IDID, RWORK(*) AND IWORK(*)
  !
  !- Description
  !
  !  Subroutine SDASSL uses the backward differentiation formulas of
  !  orders one through five to solve a system of the above form for Y and
  !  YPRIME.  Values for Y and YPRIME at the initial time must be given as
  !  input.  These values must be consistent, (that is, if T,Y,YPRIME are
  !  the given initial values, they must satisfy G(T,Y,YPRIME) = 0.).  The
  !  subroutine solves the system from T to TOUT.  It is easy to continue
  !  the solution to get results at additional TOUT.  This is the interval
  !  mode of operation.  Intermediate results can also be obtained easily
  !  by using the intermediate-output capability.
  !
  !  The following detailed description is divided into subsections:
  !    1. Input required for the first call to SDASSL.
  !    2. Output after any return from SDASSL.
  !    3. What to do to continue the integration.
  !    4. Error messages.
  !
  !
  !  -------- INPUT -- WHAT TO DO ON THE FIRST CALL TO SDASSL ------------
  !
  !  The first call of the code is defined to be the start of each new
  !  problem. Read through the descriptions of all the following items,
  !  provide sufficient storage space for designated arrays, set
  !  appropriate variables for the initialization of the problem, and
  !  give information about how you want the problem to be solved.
  !
  !
  !  RES -- Provide a subroutine of the form
  !             SUBROUTINE RES(T,Y,YPRIME,DELTA,IRES,RPAR,IPAR)
  !         to define the system of differential/algebraic
  !         equations which is to be solved. For the given values
  !         of T,Y and YPRIME, the subroutine should
  !         return the residual of the differential/algebraic
  !         system
  !             DELTA = G(T,Y,YPRIME)
  !         (DELTA(*) is a vector of length NEQ which is
  !         output for RES.)
  !
  !         Subroutine RES must not alter T,Y or YPRIME.
  !         You must declare the name RES in an external
  !         statement in your program that calls SDASSL.
  !         You must dimension Y,YPRIME and DELTA in RES.
  !
  !         IRES is an integer flag which is always equal to
  !         zero on input. Subroutine RES should alter IRES
  !         only if it encounters an illegal value of Y or
  !         a stop condition. Set IRES = -1 if an input value
  !         is illegal, and SDASSL will try to solve the problem
  !         without getting IRES = -1. If IRES = -2, SDASSL
  !         will return control to the calling program
  !         with IDID = -11.
  !
  !         RPAR and IPAR are real and integer parameter arrays which
  !         you can use for communication between your calling program
  !         and subroutine RES. They are not altered by SDASSL. If you
  !         do not need RPAR or IPAR, ignore these parameters by treat-
  !         ing them as dummy arguments. If you do choose to use them,
  !         dimension them in your calling program and in RES as arrays
  !         of appropriate length.
  !
  !  NEQ -- Set it to the number of differential equations.
  !         (NEQ .GE. 1)
  !
  !  T -- Set it to the initial point of the integration.
  !         T must be defined as a variable.
  !
  !  Y(*) -- Set this vector to the initial values of the NEQ solution
  !         components at the initial point. You must dimension Y of
  !         length at least NEQ in your calling program.
  !
  !  YPRIME(*) -- Set this vector to the initial values of the NEQ
  !         first derivatives of the solution components at the initial
  !         point.  You must dimension YPRIME at least NEQ in your
  !         calling program. If you do not know initial values of some
  !         of the solution components, see the explanation of INFO(11).
  !
  !  TOUT -- Set it to the first point at which a solution
  !         is desired. You can not take TOUT = T.
  !         integration either forward in T (TOUT .GT. T) or
  !         backward in T (TOUT .LT. T) is permitted.
  !
  !         The code advances the solution from T to TOUT using
  !         step sizes which are automatically selected so as to
  !         achieve the desired accuracy. If you wish, the code will
  !         return with the solution and its derivative at
  !         intermediate steps (intermediate-output mode) so that
  !         you can monitor them, but you still must provide TOUT in
  !         accord with the basic aim of the code.
  !
  !         The first step taken by the code is a critical one
  !         because it must reflect how fast the solution changes near
  !         the initial point. The code automatically selects an
  !         initial step size which is practically always suitable for
  !         the problem. By using the fact that the code will not step
  !         past TOUT in the first step, you could, if necessary,
  !         restrict the length of the initial step size.
  !
  !         For some problems it may not be permissible to integrate
  !         past a point TSTOP because a discontinuity occurs there
  !         or the solution or its derivative is not defined beyond
  !         TSTOP. When you have declared a TSTOP point (SEE INFO(4)
  !         and RWORK(1)), you have told the code not to integrate
  !         past TSTOP. In this case any TOUT beyond TSTOP is invalid
  !         input.
  !
  !  INFO(*) -- Use the INFO array to give the code more details about
  !         how you want your problem solved.  This array should be
  !         dimensioned of length 15, though SDASSL uses only the first
  !         eleven entries.  You must respond to all of the following
  !         items, which are arranged as questions.  The simplest use
  !         of the code corresponds to answering all questions as yes,
  !         i.e. setting all entries of INFO to 0.
  !
  !       INFO(1) - This parameter enables the code to initialize
  !              itself. You must set it to indicate the start of every
  !              new problem.
  !
  !          **** Is this the first call for this problem ...
  !                Yes - Set INFO(1) = 0
  !                 No - Not applicable here.
  !                      See below for continuation calls.  ****
  !
  !       INFO(2) - How much accuracy you want of your solution
  !              is specified by the error tolerances RTOL and ATOL.
  !              The simplest use is to take them both to be scalars.
  !              To obtain more flexibility, they can both be vectors.
  !              The code must be told your choice.
  !
  !          **** Are both error tolerances RTOL, ATOL scalars ...
  !                Yes - Set INFO(2) = 0
  !                      and input scalars for both RTOL and ATOL
  !                 No - Set INFO(2) = 1
  !                      and input arrays for both RTOL and ATOL ****
  !
  !       INFO(3) - The code integrates from T in the direction
  !              of TOUT by steps. If you wish, it will return the
  !              computed solution and derivative at the next
  !              intermediate step (the intermediate-output mode) or
  !              TOUT, whichever comes first. This is a good way to
  !              proceed if you want to see the behavior of the solution.
  !              If you must have solutions at a great many specific
  !              TOUT points, this code will compute them efficiently.
  !
  !          **** Do you want the solution only at
  !                TOUT (and not at the next intermediate step) ...
  !                 Yes - Set INFO(3) = 0
  !                  No - Set INFO(3) = 1 ****
  !
  !       INFO(4) - To handle solutions at a great many specific
  !              values TOUT efficiently, this code may integrate past
  !              TOUT and interpolate to obtain the result at TOUT.
  !              Sometimes it is not possible to integrate beyond some
  !              point TSTOP because the equation changes there or it is
  !              not defined past TSTOP. Then you must tell the code
  !              not to go past.
  !
  !           **** Can the integration be carried out without any
  !                restrictions on the independent variable T ...
  !                 Yes - Set INFO(4)=0
  !                  No - Set INFO(4)=1
  !                       and define the stopping point TSTOP by
  !                       setting RWORK(1)=TSTOP ****
  !
  !       INFO(5) - To solve differential/algebraic problems it is
  !              necessary to use a matrix of partial derivatives of the
  !              system of differential equations. If you do not
  !              provide a subroutine to evaluate it analytically (see
  !              description of the item JAC in the call list), it will
  !              be approximated by numerical differencing in this code.
  !              although it is less trouble for you to have the code
  !              compute partial derivatives by numerical differencing,
  !              the solution will be more reliable if you provide the
  !              derivatives via JAC. Sometimes numerical differencing
  !              is cheaper than evaluating derivatives in JAC and
  !              sometimes it is not - this depends on your problem.
  !
  !           **** Do you want the code to evaluate the partial
  !                derivatives automatically by numerical differences ...
  !                   Yes - Set INFO(5)=0
  !                    No - Set INFO(5)=1
  !                  and provide subroutine JAC for evaluating the
  !                  matrix of partial derivatives ****
  !
  !       INFO(6) - SDASSL will perform much better if the matrix of
  !              partial derivatives, DG/DY + CJ*DG/DYPRIME,
  !              (here CJ is a scalar determined by SDASSL)
  !              is banded and the code is told this. In this
  !              case, the storage needed will be greatly reduced,
  !              numerical differencing will be performed much cheaper,
  !              and a number of important algorithms will execute much
  !              faster. The differential equation is said to have
  !              half-bandwidths ML (lower) and MU (upper) if equation i
  !              involves only unknowns Y(J) with
  !                             I-ML .LE. J .LE. I+MU
  !              for all I=1,2,...,NEQ. Thus, ML and MU are the widths
  !              of the lower and upper parts of the band, respectively,
  !              with the main diagonal being excluded. If you do not
  !              indicate that the equation has a banded matrix of partial
  !              derivatives, the code works with a full matrix of NEQ**2
  !              elements (stored in the conventional way). Computations
  !              with banded matrices cost less time and storage than with
  !              full matrices if 2*ML+MU .LT. NEQ. If you tell the
  !              code that the matrix of partial derivatives has a banded
  !              structure and you want to provide subroutine JAC to
  !              compute the partial derivatives, then you must be careful
  !              to store the elements of the matrix in the special form
  !              indicated in the description of JAC.
  !
  !          **** Do you want to solve the problem using a full
  !               (dense) matrix (and not a special banded
  !               structure) ...
  !                Yes - Set INFO(6)=0
  !                 No - Set INFO(6)=1
  !                       and provide the lower (ML) and upper (MU)
  !                       bandwidths by setting
  !                       IWORK(1)=ML
  !                       IWORK(2)=MU ****
  !
  !
  !        INFO(7) -- You can specify a maximum (absolute value of)
  !              stepsize, so that the code
  !              will avoid passing over very
  !              large regions.
  !
  !          ****  Do you want the code to decide
  !                on its own maximum stepsize?
  !                Yes - Set INFO(7)=0
  !                 No - Set INFO(7)=1
  !                      and define HMAX by setting
  !                      RWORK(2)=HMAX ****
  !
  !        INFO(8) -- Differential/algebraic problems
  !              may occasionally suffer from
  !              severe scaling difficulties on the
  !              first step. If you know a great deal
  !              about the scaling of your problem, you can
  !              help to alleviate this problem by
  !              specifying an initial stepsize HO.
  !
  !          ****  Do you want the code to define
  !                its own initial stepsize?
  !                Yes - Set INFO(8)=0
  !                 No - Set INFO(8)=1
  !                      and define HO by setting
  !                      RWORK(3)=HO ****
  !
  !        INFO(9) -- If storage is a severe problem,
  !              you can save some locations by
  !              restricting the maximum order MAXORD.
  !              the default value is 5. for each
  !              order decrease below 5, the code
  !              requires NEQ fewer locations, however
  !              it is likely to be slower. In any
  !              case, you must have 1 .LE. MAXORD .LE. 5
  !          ****  Do you want the maximum order to
  !                default to 5?
  !                Yes - Set INFO(9)=0
  !                 No - Set INFO(9)=1
  !                      and define MAXORD by setting
  !                      IWORK(3)=MAXORD ****
  !
  !        INFO(10) --If you know that the solutions to your equations
  !               will always be nonnegative, it may help to set this
  !               parameter. However, it is probably best to
  !               try the code without using this option first,
  !               and only to use this option if that doesn't
  !               work very well.
  !           ****  Do you want the code to solve the problem without
  !                 invoking any special nonnegativity constraints?
  !                  Yes - Set INFO(10)=0
  !                   No - Set INFO(10)=1
  !
  !        INFO(11) --SDASSL normally requires the initial T,
  !               Y, and YPRIME to be consistent. That is,
  !               you must have G(T,Y,YPRIME) = 0 at the initial
  !               time. If you do not know the initial
  !               derivative precisely, you can let SDASSL try
  !               to compute it.
  !          ****   Are the initial T, Y, YPRIME consistent?
  !                 Yes - Set INFO(11) = 0
  !                  No - Set INFO(11) = 1,
  !                       and set YPRIME to an initial approximation
  !                       to YPRIME.  (If you have no idea what
  !                       YPRIME should be, set it to zero. Note
  !                       that the initial Y should be such
  !                       that there must exist a YPRIME so that
  !                       G(T,Y,YPRIME) = 0.)
  !
  !  RTOL, ATOL -- You must assign relative (RTOL) and absolute (ATOL
  !         error tolerances to tell the code how accurately you
  !         want the solution to be computed.  They must be defined
  !         as variables because the code may change them.  You
  !         have two choices --
  !               Both RTOL and ATOL are scalars. (INFO(2)=0)
  !               Both RTOL and ATOL are vectors. (INFO(2)=1)
  !         in either case all components must be non-negative.
  !
  !         The tolerances are used by the code in a local error
  !         test at each step which requires roughly that
  !               ABS(LOCAL ERROR) .LE. RTOL*ABS(Y)+ATOL
  !         for each vector component.
  !         (More specifically, a root-mean-square norm is used to
  !         measure the size of vectors, and the error test uses the
  !         magnitude of the solution at the beginning of the step.)
  !
  !         The true (global) error is the difference between the
  !         true solution of the initial value problem and the
  !         computed approximation.  Practically all present day
  !         codes, including this one, control the local error at
  !         each step and do not even attempt to control the global
  !         error directly.
  !         Usually, but not always, the true accuracy of the
  !         computed Y is comparable to the error tolerances. This
  !         code will usually, but not always, deliver a more
  !         accurate solution if you reduce the tolerances and
  !         integrate again.  By comparing two such solutions you
  !         can get a fairly reliable idea of the true error in the
  !         solution at the bigger tolerances.
  !
  !         Setting ATOL=0. results in a pure relative error test on
  !         that component.  Setting RTOL=0. results in a pure
  !         absolute error test on that component.  A mixed test
  !         with non-zero RTOL and ATOL corresponds roughly to a
  !         relative error test when the solution component is much
  !         bigger than ATOL and to an absolute error test when the
  !         solution component is smaller than the threshhold ATOL.
  !
  !         The code will not attempt to compute a solution at an
  !         accuracy unreasonable for the machine being used.  It will
  !         advise you if you ask for too much accuracy and inform
  !         you as to the maximum accuracy it believes possible.
  !
  !  RWORK(*) --  Dimension this real work array of length LRW in your
  !         calling program.
  !
  !  LRW -- Set it to the declared length of the RWORK array.
  !               You must have
  !                    LRW .GE. 40+(MAXORD+4)*NEQ+NEQ**2
  !               for the full (dense) JACOBIAN case (when INFO(6)=0), or
  !                    LRW .GE. 40+(MAXORD+4)*NEQ+(2*ML+MU+1)*NEQ
  !               for the banded user-defined JACOBIAN case
  !               (when INFO(5)=1 and INFO(6)=1), or
  !                     LRW .GE. 40+(MAXORD+4)*NEQ+(2*ML+MU+1)*NEQ
  !                           +2*(NEQ/(ML+MU+1)+1)
  !               for the banded finite-difference-generated JACOBIAN case
  !               (when INFO(5)=0 and INFO(6)=1)
  !
  !  IWORK(*) --  Dimension this integer work array of length LIW in
  !         your calling program.
  !
  !  LIW -- Set it to the declared length of the IWORK array.
  !               You must have LIW .GE. 20+NEQ
  !
  !  RPAR, IPAR -- These are parameter arrays, of real and integer
  !         type, respectively.  You can use them for communication
  !         between your program that calls SDASSL and the
  !         RES subroutine (and the JAC subroutine).  They are not
  !         altered by SDASSL.  If you do not need RPAR or IPAR,
  !         ignore these parameters by treating them as dummy
  !         arguments.  If you do choose to use them, dimension
  !         them in your calling program and in RES (and in JAC)
  !         as arrays of appropriate length.
  !
  !  JAC -- If you have set INFO(5)=0, you can ignore this parameter
  !         by treating it as a dummy argument.  Otherwise, you must
  !         provide a subroutine of the form
  !             SUBROUTINE JAC(T,Y,YPRIME,PD,CJ,RPAR,IPAR)
  !         to define the matrix of partial derivatives
  !             PD=DG/DY+CJ*DG/DYPRIME
  !         CJ is a scalar which is input to JAC.
  !         For the given values of T,Y,YPRIME, the
  !         subroutine must evaluate the non-zero partial
  !         derivatives for each equation and each solution
  !         component, and store these values in the
  !         matrix PD.  The elements of PD are set to zero
  !         before each call to JAC so only non-zero elements
  !         need to be defined.
  !
  !         Subroutine JAC must not alter T,Y,(*),YPRIME(*), or CJ.
  !         You must declare the name JAC in an EXTERNAL statement in
  !         your program that calls SDASSL.  You must dimension Y,
  !         YPRIME and PD in JAC.
  !
  !         The way you must store the elements into the PD matrix
  !         depends on the structure of the matrix which you
  !         indicated by INFO(6).
  !               *** INFO(6)=0 -- Full (dense) matrix ***
  !                   Give PD a first dimension of NEQ.
  !                   When you evaluate the (non-zero) partial derivative
  !                   of equation I with respect to variable J, you must
  !                   store it in PD according to
  !                   PD(I,J) = "DG(I)/DY(J)+CJ*DG(I)/DYPRIME(J)"
  !               *** INFO(6)=1 -- Banded JACOBIAN with ML lower and MU
  !                   upper diagonal bands (refer to INFO(6) description
  !                   of ML and MU) ***
  !                   Give PD a first dimension of 2*ML+MU+1.
  !                   when you evaluate the (non-zero) partial derivative
  !                   of equation I with respect to variable J, you must
  !                   store it in PD according to
  !                   IROW = I - J + ML + MU + 1
  !                   PD(IROW,J) = "DG(I)/DY(J)+CJ*DG(I)/DYPRIME(J)"
  !
  !         RPAR and IPAR are real and integer parameter arrays
  !         which you can use for communication between your calling
  !         program and your JACOBIAN subroutine JAC. They are not
  !         altered by SDASSL. If you do not need RPAR or IPAR,
  !         ignore these parameters by treating them as dummy
  !         arguments. If you do choose to use them, dimension
  !         them in your calling program and in JAC as arrays of
  !         appropriate length.
  !
  !
  !  OPTIONALLY REPLACEABLE NORM ROUTINE:
  !
  !     SDASSL uses a weighted norm SDANRM to measure the size
  !     of vectors such as the estimated error in each step.
  !     A FUNCTION subprogram
  !       REAL FUNCTION SDANRM(NEQ,V,WT,RPAR,IPAR)
  !       DIMENSION V(NEQ),WT(NEQ)
  !     is used to define this norm. Here, V is the vector
  !     whose norm is to be computed, and WT is a vector of
  !     weights.  A SDANRM routine has been included with SDASSL
  !     which computes the weighted root-mean-square norm
  !     given by
  !       SDANRM=SQRT((1/NEQ)*SUM(V(I)/WT(I))**2)
  !     this norm is suitable for most problems. In some
  !     special cases, it may be more convenient and/or
  !     efficient to define your own norm by writing a function
  !     subprogram to be called instead of SDANRM. This should,
  !     however, be attempted only after careful thought and
  !     consideration.
  !
  !
  !  -------- OUTPUT -- AFTER ANY RETURN FROM SDASSL ---------------------
  !
  !  The principal aim of the code is to return a computed solution at
  !  TOUT, although it is also possible to obtain intermediate results
  !  along the way. To find out whether the code achieved its goal
  !  or if the integration process was interrupted before the task was
  !  completed, you must check the IDID parameter.
  !
  !
  !  T -- The solution was successfully advanced to the
  !               output value of T.
  !
  !  Y(*) -- Contains the computed solution approximation at T.
  !
  !  YPRIME(*) -- Contains the computed derivative
  !               approximation at T.
  !
  !  IDID -- Reports what the code did.
  !
  !                     *** Task completed ***
  !                Reported by positive values of IDID
  !
  !           IDID = 1 -- A step was successfully taken in the
  !                   intermediate-output mode. The code has not
  !                   yet reached TOUT.
  !
  !           IDID = 2 -- The integration to TSTOP was successfully
  !                   completed (T=TSTOP) by stepping exactly to TSTOP.
  !
  !           IDID = 3 -- The integration to TOUT was successfully
  !                   completed (T=TOUT) by stepping past TOUT.
  !                   Y(*) is obtained by interpolation.
  !                   YPRIME(*) is obtained by interpolation.
  !
  !                    *** Task interrupted ***
  !                Reported by negative values of IDID
  !
  !           IDID = -1 -- A large amount of work has been expended.
  !                   (About 500 steps)
  !
  !           IDID = -2 -- The error tolerances are too stringent.
  !
  !           IDID = -3 -- The local error test cannot be satisfied
  !                   because you specified a zero component in ATOL
  !                   and the corresponding computed solution
  !                   component is zero. Thus, a pure relative error
  !                   test is impossible for this component.
  !
  !           IDID = -6 -- SDASSL had repeated error test
  !                   failures on the last attempted step.
  !
  !           IDID = -7 -- The corrector could not converge.
  !
  !           IDID = -8 -- The matrix of partial derivatives
  !                   is singular.
  !
  !           IDID = -9 -- The corrector could not converge.
  !                   there were repeated error test failures
  !                   in this step.
  !
  !           IDID =-10 -- The corrector could not converge
  !                   because IRES was equal to minus one.
  !
  !           IDID =-11 -- IRES equal to -2 was encountered
  !                   and control is being returned to the
  !                   calling program.
  !
  !           IDID =-12 -- SDASSL failed to compute the initial
  !                   YPRIME.
  !
  !
  !
  !           IDID = -13,..,-32 -- Not applicable for this code
  !
  !                    *** Task terminated ***
  !                Reported by the value of IDID=-33
  !
  !           IDID = -33 -- The code has encountered trouble from which
  !                   it cannot recover. A message is printed
  !                   explaining the trouble and control is returned
  !                   to the calling program. For example, this occurs
  !                   when invalid input is detected.
  !
  !  RTOL, ATOL -- These quantities remain unchanged except when
  !               IDID = -2. In this case, the error tolerances have been
  !               increased by the code to values which are estimated to
  !               be appropriate for continuing the integration. However,
  !               the reported solution at T was obtained using the input
  !               values of RTOL and ATOL.
  !
  !  RWORK, IWORK -- Contain information which is usually of no
  !               interest to the user but necessary for subsequent calls.
  !               However, you may find use for
  !
  !               RWORK(3)--Which contains the step size H to be
  !                       attempted on the next step.
  !
  !               RWORK(4)--Which contains the current value of the
  !                       independent variable, i.e., the farthest point
  !                       integration has reached. This will be different
  !                       from T only when interpolation has been
  !                       performed (IDID=3).
  !
  !               RWORK(7)--Which contains the stepsize used
  !                       on the last successful step.
  !
  !               IWORK(7)--Which contains the order of the method to
  !                       be attempted on the next step.
  !
  !               IWORK(8)--Which contains the order of the method used
  !                       on the last step.
  !
  !               IWORK(11)--Which contains the number of steps taken so
  !                        far.
  !
  !               IWORK(12)--Which contains the number of calls to RES
  !                        so far.
  !
  !               IWORK(13)--Which contains the number of evaluations of
  !                        the matrix of partial derivatives needed so
  !                        far.
  !
  !               IWORK(14)--Which contains the total number
  !                        of error test failures so far.
  !
  !               IWORK(15)--Which contains the total number
  !                        of convergence test failures so far.
  !                        (includes singular iteration matrix
  !                        failures.)
  !
  !
  !  -------- INPUT -- WHAT TO DO TO CONTINUE THE INTEGRATION ------------
  !                    (CALLS AFTER THE FIRST)
  !
  !  This code is organized so that subsequent calls to continue the
  !  integration involve little (if any) additional effort on your
  !  part. You must monitor the IDID parameter in order to determine
  !  what to do next.
  !
  !  Recalling that the principal task of the code is to integrate
  !  from T to TOUT (the interval mode), usually all you will need
  !  to do is specify a new TOUT upon reaching the current TOUT.
  !
  !  Do not alter any quantity not specifically permitted below,
  !  in particular do not alter NEQ,T,Y(*),YPRIME(*),RWORK(*),IWORK(*)
  !  or the differential equation in subroutine RES. Any such
  !  alteration constitutes a new problem and must be treated as such,
  !  i.e., you must start afresh.
  !
  !  You cannot change from vector to scalar error control or vice
  !  versa (INFO(2)), but you can change the size of the entries of
  !  RTOL, ATOL. Increasing a tolerance makes the equation easier
  !  to integrate. Decreasing a tolerance will make the equation
  !  harder to integrate and should generally be avoided.
  !
  !  You can switch from the intermediate-output mode to the
  !  interval mode (INFO(3)) or vice versa at any time.
  !
  !  If it has been necessary to prevent the integration from going
  !  past a point TSTOP (INFO(4), RWORK(1)), keep in mind that the
  !  code will not integrate to any TOUT beyond the currently
  !  specified TSTOP. Once TSTOP has been reached you must change
  !  the value of TSTOP or set INFO(4)=0. You may change INFO(4)
  !  or TSTOP at any time but you must supply the value of TSTOP in
  !  RWORK(1) whenever you set INFO(4)=1.
  !
  !  Do not change INFO(5), INFO(6), IWORK(1), or IWORK(2)
  !  unless you are going to restart the code.
  !
  !                 *** Following a completed task ***
  !  If
  !     IDID = 1, call the code again to continue the integration
  !                  another step in the direction of TOUT.
  !
  !     IDID = 2 or 3, define a new TOUT and call the code again.
  !                  TOUT must be different from T. You cannot change
  !                  the direction of integration without restarting.
  !
  !                 *** Following an interrupted task ***
  !               To show the code that you realize the task was
  !               interrupted and that you want to continue, you
  !               must take appropriate action and set INFO(1) = 1
  !  If
  !    IDID = -1, The code has taken about 500 steps.
  !                  If you want to continue, set INFO(1) = 1 and
  !                  call the code again. An additional 500 steps
  !                  will be allowed.
  !
  !    IDID = -2, The error tolerances RTOL, ATOL have been
  !                  increased to values the code estimates appropriate
  !                  for continuing. You may want to change them
  !                  yourself. If you are sure you want to continue
  !                  with relaxed error tolerances, set INFO(1)=1 and
  !                  call the code again.
  !
  !    IDID = -3, A solution component is zero and you set the
  !                  corresponding component of ATOL to zero. If you
  !                  are sure you want to continue, you must first
  !                  alter the error criterion to use positive values
  !                  for those components of ATOL corresponding to zero
  !                  solution components, then set INFO(1)=1 and call
  !                  the code again.
  !
  !    IDID = -4,-5  --- Cannot occur with this code.
  !
  !    IDID = -6, Repeated error test failures occurred on the
  !                  last attempted step in SDASSL. A singularity in the
  !                  solution may be present. If you are absolutely
  !                  certain you want to continue, you should restart
  !                  the integration. (Provide initial values of Y and
  !                  YPRIME which are consistent)
  !
  !    IDID = -7, Repeated convergence test failures occurred
  !                  on the last attempted step in SDASSL. An inaccurate
  !                  or ill-conditioned JACOBIAN may be the problem. If
  !                  you are absolutely certain you want to continue, you
  !                  should restart the integration.
  !
  !    IDID = -8, The matrix of partial derivatives is singular.
  !                  Some of your equations may be redundant.
  !                  SDASSL cannot solve the problem as stated.
  !                  It is possible that the redundant equations
  !                  could be removed, and then SDASSL could
  !                  solve the problem. It is also possible
  !                  that a solution to your problem either
  !                  does not exist or is not unique.
  !
  !    IDID = -9, SDASSL had multiple convergence test
  !                  failures, preceded by multiple error
  !                  test failures, on the last attempted step.
  !                  It is possible that your problem
  !                  is ill-posed, and cannot be solved
  !                  using this code. Or, there may be a
  !                  discontinuity or a singularity in the
  !                  solution. If you are absolutely certain
  !                  you want to continue, you should restart
  !                  the integration.
  !
  !    IDID =-10, SDASSL had multiple convergence test failures
  !                  because IRES was equal to minus one.
  !                  If you are absolutely certain you want
  !                  to continue, you should restart the
  !                  integration.
  !
  !    IDID =-11, IRES=-2 was encountered, and control is being
  !                  returned to the calling program.
  !
  !    IDID =-12, SDASSL failed to compute the initial YPRIME.
  !                  This could happen because the initial
  !                  approximation to YPRIME was not very good, or
  !                  if a YPRIME consistent with the initial Y
  !                  does not exist. The problem could also be caused
  !                  by an inaccurate or singular iteration matrix.
  !
  !    IDID = -13,..,-32  --- Cannot occur with this code.
  !
  !
  !                 *** Following a terminated task ***
  !
  !  If IDID= -33, you cannot continue the solution of this problem.
  !                  An attempt to do so will result in your
  !                  run being terminated.
  !
  !
  !  -------- ERROR MESSAGES ---------------------------------------------
  !
  !      The SLATEC error print routine XERMSG is called in the event of
  !   unsuccessful completion of a task.  Most of these are treated as
  !   "recoverable errors", which means that (unless the user has directed
  !   otherwise) control will be returned to the calling program for
  !   possible action after the message has been printed.
  !
  !   In the event of a negative value of IDID other than -33, an appro-
  !   priate message is printed and the "error number" printed by XERMSG
  !   is the value of IDID.  There are quite a number of illegal input
  !   errors that can lead to a returned value IDID=-33.  The conditions
  !   and their printed "error numbers" are as follows:
  !
  !   Error number       Condition
  !
  !        1       Some element of INFO vector is not zero or one.
  !        2       NEQ .le. 0
  !        3       MAXORD not in range.
  !        4       LRW is less than the required length for RWORK.
  !        5       LIW is less than the required length for IWORK.
  !        6       Some element of RTOL is .lt. 0
  !        7       Some element of ATOL is .lt. 0
  !        8       All elements of RTOL and ATOL are zero.
  !        9       INFO(4)=1 and TSTOP is behind TOUT.
  !       10       HMAX .lt. 0.0
  !       11       TOUT is behind T.
  !       12       INFO(8)=1 and H0=0.0
  !       13       Some element of WT is .le. 0.0
  !       14       TOUT is too close to T to start integration.
  !       15       INFO(4)=1 and TSTOP is behind T.
  !       16       --( Not used in this version )--
  !       17       ML illegal.  Either .lt. 0 or .gt. NEQ
  !       18       MU illegal.  Either .lt. 0 or .gt. NEQ
  !       19       TOUT = T.
  !
  !   If SDASSL is called again without any action taken to remove the
  !   cause of an unsuccessful return, XERMSG will be called with a fatal
  !   error flag, which will cause unconditional termination of the
  !   program.  There are two such fatal errors:
  !
  !   Error number -998:  The last step was terminated with a negative
  !       value of IDID other than -33, and no appropriate action was
  !       taken.
  !
  !   Error number -999:  The previous call was terminated because of
  !       illegal input (IDID=-33) and there is illegal input in the
  !       present call, as well.  (Suspect infinite loop.)
  !
  !  ---------------------------------------------------------------------
  !
  !***
  ! **References:**  A DESCRIPTION OF DASSL: A DIFFERENTIAL/ALGEBRAIC
  !                 SYSTEM SOLVER, L. R. PETZOLD, SAND82-8637,
  !                 SANDIA NATIONAL LABORATORIES, SEPTEMBER 1982.
  !***
  ! **Routines called:**  R1MACH, SDAINI, SDANRM, SDASTP, SDATRP, SDAWTS,
  !                    XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   830315  DATE WRITTEN
  !   880387  Code changes made.  All common statements have been
  !           replaced by a DATA statement, which defines pointers into
  !           RWORK, and PARAMETER statements which define pointers
  !           into IWORK.  As well the documentation has gone through
  !           grammatical changes.
  !   881005  The prologue has been changed to mixed case.
  !           The subordinate routines had revision dates changed to
  !           this date, although the documentation for these routines
  !           is all upper case.  No code changes.
  !   890511  Code changes made.  The DATA statement in the declaration
  !           section of SDASSL was replaced with a PARAMETER
  !           statement.  Also the statement S = 100.E0 was removed
  !           from the top of the Newton iteration in SDASTP.
  !           The subordinate routines had revision dates changed to
  !           this date.
  !   890517  The revision date syntax was replaced with the revision
  !           history syntax.  Also the "DECK" comment was added to
  !           the top of all subroutines.  These changes are consistent
  !           with new SLATEC guidelines.
  !           The subordinate routines had revision dates changed to
  !           this date.  No code changes.
  !   891013  Code changes made.
  !           Removed all occurrences of FLOAT.  All operations
  !           are now performed with "mixed-mode" arithmetic.
  !           Also, specific function names were replaced with generic
  !           function names to be consistent with new SLATEC guidelines.
  !           In particular:
  !              Replaced AMIN1 with MIN everywhere.
  !              Replaced MIN0 with MIN everywhere.
  !              Replaced AMAX1 with MAX everywhere.
  !              Replaced MAX0 with MAX everywhere.
  !           Also replaced REVISION DATE with REVISION HISTORY in all
  !           subordinate routines.
  !   901004  Miscellaneous changes to prologue to complete conversion
  !           to SLATEC 4.0 format.  No code changes.  (F.N.Fritsch)
  !   901009  Corrected GAMS classification code and converted subsidiary
  !           routines to 4.0 format.  No code changes.  (F.N.Fritsch)
  !   901010  Converted XERRWV calls to XERMSG calls.  (R.Clemens, AFWL)
  !   901019  Code changes made.
  !           Merged SLATEC 4.0 changes with previous changes made
  !           by C. Ulrich.  Below is a history of the changes made by
  !           C. Ulrich. (Changes in subsidiary routines are implied
  !           by this history)
  !           891228  Bug was found and repaired inside the SDASSL
  !                   and SDAINI routines.  SDAINI was incorrectly
  !                   returning the initial T with Y and YPRIME
  !                   computed at T+H.  The routine now returns T+H
  !                   rather than the initial T.
  !                   Cosmetic changes made to SDASTP.
  !           900904  Three modifications were made to fix a bug (inside
  !                   SDASSL) re interpolation for continuation calls and
  !                   cases where TN is very close to TSTOP:
  !
  !                   1) In testing for whether H is too large, just
  !                      compare H to (TSTOP - TN), rather than
  !                      (TSTOP - TN) * (1-4*UROUND), and set H to
  !                      TSTOP - TN.  This will force SDASTP to step
  !                      exactly to TSTOP under certain situations
  !                      (i.e. when H returned from SDASTP would otherwise
  !                      take TN beyond TSTOP).
  !
  !                   2) Inside the SDASTP loop, interpolate exactly to
  !                      TSTOP if TN is very close to TSTOP (rather than
  !                      interpolating to within roundoff of TSTOP).
  !
  !                   3) Modified IDID description for IDID = 2 to say
  !                      that the solution is returned by stepping exactly
  !                      to TSTOP, rather than TOUT.  (In some cases the
  !                      solution is actually obtained by extrapolating
  !                      over a distance near unit roundoff to TSTOP,
  !                      but this small distance is deemed acceptable in
  !                      these circumstances.)
  !   901026  Added explicit declarations for all variables and minor
  !           cosmetic changes to prologue, removed unreferenced labels,
  !           and improved XERMSG calls.  (FNF)
  !   901030  Added ERROR MESSAGES section and reworked other sections to
  !           be of more uniform format.  (FNF)
  !   910624  Fixed minor bug related to HMAX (six lines after label
  !           525).  (LRP)
  USE service, ONLY : XERMSG, R1MACH
  !     Declare arguments.
  INTERFACE
    SUBROUTINE RES(T,Y,Yprime,Delta,Ires)
      IMPORT SP
      INTEGER :: Ires
      REAL(SP) :: T, Y(:), Yprime(:), Delta(:)
    END SUBROUTINE
    SUBROUTINE JAC(T,Y,Yprime,Pd,Cj)
      IMPORT SP
      REAL(SP) :: T, Cj, Pd(:,:), Y(:), Yprime(:)
    END SUBROUTINE
  END INTERFACE
  INTEGER :: Neq, Idid, Lrw, Liw
  INTEGER :: Info(15), Iwork(Liw)
  REAL(SP) :: T, Tout
  REAL(SP) :: Y(Neq), Yprime(Neq), Rtol(:), Atol(:), Rwork(Lrw)
  !
  !     Declare local variables.
  !
  INTEGER i, itemp, leniw, lenpd, lenrw, le, lpd, lphi, lwm, lwt, mband, msave, &
    mxord, ntemp, nzflg
  REAL(SP) atoli, h, hmax, hmin, ho, r, rh, rtoli, tdist, tn, tnext, tstop, uround, ypnorm
  LOGICAL done
  !       Auxiliary variables for conversion of values to be included in
  !       error messages.
  CHARACTER(8) :: xern1, xern2
  CHARACTER(16) :: xern3, xern4
  !
  !     SET POINTERS INTO IWORK
  INTEGER, PARAMETER :: LML = 1, LMU = 2, LMXORD = 3, LMTYPE = 4, LNST = 11, &
    LNRE = 12, LNJE = 13, LNPD = 16, LJCALC = 5, LPHASE = 6, LK = 7, LKOLD = 8, &
    LNS = 9, LNSTL = 10
  !
  !     SET RELATIVE OFFSET INTO RWORK
  INTEGER, PARAMETER :: NPD = 1
  !
  !     SET POINTERS INTO RWORK
  INTEGER, PARAMETER :: LTSTOP = 1, LHMAX = 2, LH = 3, LTN = 4, LCJ = 5, LCJOLD = 6, &
    LHOLD = 7, LS = 8, LROUND = 9, LALPHA = 11, LBETA = 17, LGAMMA = 23, LPSI = 29, &
    LSIGMA = 35, LDELTA = 41
  !
  !* FIRST EXECUTABLE STATEMENT  SDASSL
  IF ( Info(1)==0 ) THEN
    !
    !-----------------------------------------------------------------------
    !     THIS BLOCK IS EXECUTED FOR THE INITIAL CALL ONLY.
    !     IT CONTAINS CHECKING OF INPUTS AND INITIALIZATIONS.
    !-----------------------------------------------------------------------
    !
    !     FIRST CHECK INFO ARRAY TO MAKE SURE ALL ELEMENTS OF INFO
    !     ARE EITHER ZERO OR ONE.
    DO i = 2, 11
      IF ( Info(i)/=0.AND.Info(i)/=1 ) GOTO 400
    END DO
    !
    IF ( Neq<=0 ) THEN
      !
      WRITE (xern1,'(I8)') Neq
      CALL XERMSG('SDASSL','NEQ = '//xern1//' .LE. 0',2,1)
      GOTO 1200
    ELSE
      !
      !     CHECK AND COMPUTE MAXIMUM ORDER
      mxord = 5
      IF ( Info(9)/=0 ) THEN
        mxord = Iwork(LMXORD)
        IF ( mxord<1.OR.mxord>5 ) THEN
          !
          WRITE (xern1,'(I8)') mxord
          CALL XERMSG('SDASSL','MAXORD = '//xern1//' NOT IN RANGE',3,1)
          GOTO 1200
        END IF
      END IF
      Iwork(LMXORD) = mxord
      !
      !     COMPUTE MTYPE,LENPD,LENRW.CHECK ML AND MU.
      IF ( Info(6)==0 ) THEN
        lenpd = Neq**2
        lenrw = 40 + (Iwork(LMXORD)+4)*Neq + lenpd
        IF ( Info(5)/=0 ) THEN
          Iwork(LMTYPE) = 1
        ELSE
          Iwork(LMTYPE) = 2
        END IF
      ELSEIF ( Iwork(LML)<0.OR.Iwork(LML)>=Neq ) THEN
        !
        WRITE (xern1,'(I8)') Iwork(LML)
        CALL XERMSG('SDASSL','ML = '//xern1//&
          ' ILLEGAL.  EITHER .LT. 0 OR .GT. NEQ',17,1)
        GOTO 1200
      ELSEIF ( Iwork(LMU)<0.OR.Iwork(LMU)>=Neq ) THEN
        !
        WRITE (xern1,'(I8)') Iwork(LMU)
        CALL XERMSG('SDASSL','MU = '//xern1//&
          ' ILLEGAL.  EITHER .LT. 0 OR .GT. NEQ',18,1)
        GOTO 1200
      ELSE
        lenpd = (2*Iwork(LML)+Iwork(LMU)+1)*Neq
        IF ( Info(5)/=0 ) THEN
          Iwork(LMTYPE) = 4
          lenrw = 40 + (Iwork(LMXORD)+4)*Neq + lenpd
        ELSE
          Iwork(LMTYPE) = 5
          mband = Iwork(LML) + Iwork(LMU) + 1
          msave = (Neq/mband) + 1
          lenrw = 40 + (Iwork(LMXORD)+4)*Neq + lenpd + 2*msave
        END IF
      END IF
      !
      !     CHECK LENGTHS OF RWORK AND IWORK
      leniw = 20 + Neq
      Iwork(LNPD) = lenpd
      IF ( Lrw<lenrw ) THEN
        !
        WRITE (xern1,'(I8)') lenrw
        WRITE (xern2,'(I8)') Lrw
        CALL XERMSG('SDASSL','RWORK LENGTH NEEDED, LENRW = '//&
          xern1//', EXCEEDS LRW = '//xern2,4,1)
        GOTO 1200
      ELSEIF ( Liw<leniw ) THEN
        !
        WRITE (xern1,'(I8)') leniw
        WRITE (xern2,'(I8)') Liw
        CALL XERMSG('SDASSL','IWORK LENGTH NEEDED, LENIW = '//&
          xern1//', EXCEEDS LIW = '//xern2,5,1)
        GOTO 1200
      ELSE
        !
        !     CHECK TO SEE THAT TOUT IS DIFFERENT FROM T
        IF ( Tout==T ) GOTO 1100
        !
        !     CHECK HMAX
        IF ( Info(7)/=0 ) THEN
          hmax = Rwork(LHMAX)
          IF ( hmax<=0.0E0 ) THEN
            !
            WRITE (xern3,'(1P,E15.6)') hmax
            CALL XERMSG('SDASSL','HMAX = '//xern3//' .LT. 0.0',10,1)
            GOTO 1200
          END IF
        END IF
        !
        !     INITIALIZE COUNTERS
        Iwork(LNST) = 0
        Iwork(LNRE) = 0
        Iwork(LNJE) = 0
        !
        Iwork(LNSTL) = 0
        Idid = 1
      END IF
    END IF
    !
    !-----------------------------------------------------------------------
    !     THIS BLOCK IS FOR CONTINUATION CALLS
    !     ONLY. HERE WE CHECK INFO(1), AND IF THE
    !     LAST STEP WAS INTERRUPTED WE CHECK WHETHER
    !     APPROPRIATE ACTION WAS TAKEN.
    !-----------------------------------------------------------------------
    !
  ELSEIF ( Info(1)==1 ) THEN
    Iwork(LNSTL) = Iwork(LNST)
  ELSE
    IF ( Info(1)/=-1 ) GOTO 400
    !
    !     IF WE ARE HERE, THE LAST STEP WAS INTERRUPTED
    !     BY AN ERROR CONDITION FROM SDASTP, AND
    !     APPROPRIATE ACTION WAS NOT TAKEN. THIS
    !     IS A FATAL ERROR.
    WRITE (xern1,'(I8)') Idid
    CALL XERMSG('SDASSL',&
      'THE LAST STEP TERMINATED WITH A NEGATIVE VALUE OF IDID = '&
      //xern1//' AND NO APPROPRIATE ACTION WAS TAKEN.  RUN TERMINATED',-998,2)
    RETURN
  END IF
  !
  !-----------------------------------------------------------------------
  !     THIS BLOCK IS EXECUTED ON ALL CALLS.
  !     THE ERROR TOLERANCE PARAMETERS ARE
  !     CHECKED, AND THE WORK ARRAY POINTERS
  !     ARE SET.
  !-----------------------------------------------------------------------
  !
  !     CHECK RTOL,ATOL
  nzflg = 0
  rtoli = Rtol(1)
  atoli = Atol(1)
  DO i = 1, Neq
    IF ( Info(2)==1 ) rtoli = Rtol(i)
    IF ( Info(2)==1 ) atoli = Atol(i)
    IF ( rtoli>0.0E0.OR.atoli>0.0E0 ) nzflg = 1
    IF ( rtoli<0.0E0 ) GOTO 500
    IF ( atoli<0.0E0 ) GOTO 600
  END DO
  IF ( nzflg==0 ) THEN
    !
    CALL XERMSG('SDASSL','ALL ELEMENTS OF RTOL AND ATOL ARE ZERO',8,1)
    GOTO 1200
  ELSE
    !
    !     SET UP RWORK STORAGE.IWORK STORAGE IS FIXED
    !     IN DATA STATEMENT.
    le = LDELTA + Neq
    lwt = le + Neq
    lphi = lwt + Neq
    lpd = lphi + (Iwork(LMXORD)+1)*Neq
    lwm = lpd
    ntemp = NPD + Iwork(LNPD)
    IF ( Info(1)==1 ) THEN
      !
      !-------------------------------------------------------
      !     THIS BLOCK IS FOR CONTINUATION CALLS ONLY. ITS
      !     PURPOSE IS TO CHECK STOP CONDITIONS BEFORE
      !     TAKING A STEP.
      !     ADJUST H IF NECESSARY TO MEET HMAX BOUND
      !-------------------------------------------------------
      !
      uround = Rwork(LROUND)
      done = .FALSE.
      tn = Rwork(LTN)
      h = Rwork(LH)
      IF ( Info(7)/=0 ) THEN
        rh = ABS(h)/Rwork(LHMAX)
        IF ( rh>1.0E0 ) h = h/rh
      END IF
      IF ( T==Tout ) GOTO 1100
      IF ( (T-Tout)*h>0.0E0 ) GOTO 800
      IF ( Info(4)==1 ) THEN
        IF ( Info(3)==1 ) THEN
          tstop = Rwork(LTSTOP)
          IF ( (tn-tstop)*h>0.0E0 ) GOTO 1000
          IF ( (tstop-Tout)*h<0.0E0 ) GOTO 700
          IF ( (tn-T)*h>0.0E0 ) THEN
            IF ( (tn-Tout)*h>0.0E0 ) THEN
              CALL SDATRP(tn,Tout,Y,Yprime,Neq,Iwork(LKOLD),Rwork(lphi),Rwork(LPSI))
              T = Tout
              Idid = 3
              done = .TRUE.
            ELSE
              CALL SDATRP(tn,tn,Y,Yprime,Neq,Iwork(LKOLD),Rwork(lphi),Rwork(LPSI))
              T = tn
              Idid = 1
              done = .TRUE.
            END IF
            GOTO 20
          END IF
        ELSE
          tstop = Rwork(LTSTOP)
          IF ( (tn-tstop)*h>0.0E0 ) GOTO 1000
          IF ( (tstop-Tout)*h<0.0E0 ) GOTO 700
          IF ( (tn-Tout)*h>=0.0E0 ) THEN
            CALL SDATRP(tn,Tout,Y,Yprime,Neq,Iwork(LKOLD),Rwork(lphi),Rwork(LPSI))
            T = Tout
            Idid = 3
            done = .TRUE.
            GOTO 20
          END IF
        END IF
        !     CHECK WHETHER WE ARE WITHIN ROUNDOFF OF TSTOP
        IF ( ABS(tn-tstop)>100.0E0*uround*(ABS(tn)+ABS(h)) ) THEN
          tnext = tn + h
          IF ( (tnext-tstop)*h>0.0E0 ) THEN
            h = tstop - tn
            Rwork(LH) = h
          END IF
        ELSE
          CALL SDATRP(tn,tstop,Y,Yprime,Neq,Iwork(LKOLD),Rwork(lphi),Rwork(LPSI))
          Idid = 2
          T = tstop
          done = .TRUE.
        END IF
      ELSEIF ( Info(3)==1 ) THEN
        IF ( (tn-T)*h>0.0E0 ) THEN
          IF ( (tn-Tout)*h>0.0E0 ) THEN
            CALL SDATRP(tn,Tout,Y,Yprime,Neq,Iwork(LKOLD),Rwork(lphi),Rwork(LPSI))
            T = Tout
            Idid = 3
            done = .TRUE.
          ELSE
            CALL SDATRP(tn,tn,Y,Yprime,Neq,Iwork(LKOLD),Rwork(lphi),Rwork(LPSI))
            T = tn
            Idid = 1
            done = .TRUE.
          END IF
        END IF
      ELSEIF ( (tn-Tout)*h>=0.0E0 ) THEN
        CALL SDATRP(tn,Tout,Y,Yprime,Neq,Iwork(LKOLD),Rwork(lphi),Rwork(LPSI))
        T = Tout
        Idid = 3
        done = .TRUE.
      END IF
      !
      20  IF ( done ) GOTO 300
    ELSE
      !
      !-----------------------------------------------------------------------
      !     THIS BLOCK IS EXECUTED ON THE INITIAL CALL
      !     ONLY. SET THE INITIAL STEP SIZE, AND
      !     THE ERROR WEIGHT VECTOR, AND PHI.
      !     COMPUTE INITIAL YPRIME, IF NECESSARY.
      !-----------------------------------------------------------------------
      !
      tn = T
      Idid = 1
      !
      !     SET ERROR WEIGHT VECTOR WT
      CALL SDAWTS(Neq,Info(2),Rtol,Atol,Y,Rwork(lwt))
      DO i = 1, Neq
        IF ( Rwork(lwt+i-1)<=0.0E0 ) GOTO 900
      END DO
      !
      !     COMPUTE UNIT ROUNDOFF AND HMIN
      uround = R1MACH(4)
      Rwork(LROUND) = uround
      hmin = 4.0E0*uround*MAX(ABS(T),ABS(Tout))
      !
      !     CHECK INITIAL INTERVAL TO SEE THAT IT IS LONG ENOUGH
      tdist = ABS(Tout-T)
      IF ( tdist<hmin ) THEN
        !
        WRITE (xern3,'(1P,E15.6)') Tout
        WRITE (xern4,'(1P,E15.6)') T
        CALL XERMSG('SDASSL','TOUT = '//xern3//&
          ' TOO CLOSE TO T = '//xern4//' TO START INTEGRATION',14,1)
        GOTO 1200
      ELSE
        !
        !     CHECK HO, IF THIS WAS INPUT
        IF ( Info(8)==0 ) THEN
          !
          !     COMPUTE INITIAL STEPSIZE, TO BE USED BY EITHER
          !     SDASTP OR SDAINI, DEPENDING ON INFO(11)
          ho = 0.001E0*tdist
          ypnorm = SDANRM(Neq,Yprime,Rwork(lwt))
          IF ( ypnorm>0.5E0/ho ) ho = 0.5E0/ypnorm
          ho = SIGN(ho,Tout-T)
        ELSE
          ho = Rwork(LH)
          IF ( (Tout-T)*ho<0.0E0 ) GOTO 800
          IF ( ho==0.0E0 ) THEN
            !
            CALL XERMSG('SDASSL','INFO(8)=1 AND H0=0.0',12,1)
            GOTO 1200
          END IF
        END IF
        !     ADJUST HO IF NECESSARY TO MEET HMAX BOUND
        IF ( Info(7)/=0 ) THEN
          rh = ABS(ho)/Rwork(LHMAX)
          IF ( rh>1.0E0 ) ho = ho/rh
        END IF
        !     COMPUTE TSTOP, IF APPLICABLE
        IF ( Info(4)/=0 ) THEN
          tstop = Rwork(LTSTOP)
          IF ( (tstop-T)*ho<0.0E0 ) GOTO 1000
          IF ( (T+ho-tstop)*ho>0.0E0 ) ho = tstop - T
          IF ( (tstop-Tout)*ho<0.0E0 ) GOTO 700
        END IF
        !
        !     COMPUTE INITIAL DERIVATIVE, UPDATING TN AND Y, IF APPLICABLE
        IF ( Info(11)/=0 ) THEN
          CALL SDAINI(tn,Y,Yprime,Neq,RES,JAC,ho,Rwork(lwt:lphi-1),Idid,&
            Rwork(lphi:lpd-1),Rwork(LDELTA:lwt-1),Rwork(le:lwt-1),Rwork(lwm:),&
            Iwork,hmin,Rwork(LROUND),Info(10),ntemp)
          IF ( Idid<0 ) GOTO 100
        END IF
        !
        !     LOAD H WITH HO.  STORE H IN RWORK(LH)
        h = ho
        Rwork(LH) = h
        !
        !     LOAD Y AND H*YPRIME INTO PHI(*,1) AND PHI(*,2)
        itemp = lphi + Neq
        DO i = 1, Neq
          Rwork(lphi+i-1) = Y(i)
          Rwork(itemp+i-1) = h*Yprime(i)
          !
        END DO
      END IF
    END IF
  END IF
  !
  !-------------------------------------------------------
  !     THE NEXT BLOCK CONTAINS THE CALL TO THE
  !     ONE-STEP INTEGRATOR SDASTP.
  !     THIS IS A LOOPING POINT FOR THE INTEGRATION STEPS.
  !     CHECK FOR TOO MANY STEPS.
  !     UPDATE WT.
  !     CHECK FOR TOO MUCH ACCURACY REQUESTED.
  !     COMPUTE MINIMUM STEPSIZE.
  !-------------------------------------------------------
  !
  !     CHECK FOR FAILURE TO COMPUTE INITIAL YPRIME
  100 CONTINUE
  IF ( Idid/=-12 ) THEN
    !
    !     CHECK FOR TOO MANY STEPS
    IF ( (Iwork(LNST)-Iwork(LNSTL))<500 ) THEN
      !
      !     UPDATE WT
      CALL SDAWTS(Neq,Info(2),Rtol,Atol,Rwork(lphi),Rwork(lwt))
      DO i = 1, Neq
        IF ( Rwork(i+lwt-1)<=0.0E0 ) THEN
          Idid = -3
          GOTO 200
        END IF
      END DO
      !
      !     TEST FOR TOO MUCH ACCURACY REQUESTED.
      r = SDANRM(Neq,Rwork(lphi),Rwork(lwt))*100.0E0*uround
      IF ( r<=1.0E0 ) THEN
        !
        !     COMPUTE MINIMUM STEPSIZE
        hmin = 4.0E0*uround*MAX(ABS(tn),ABS(Tout))
        !
        !     TEST H VS. HMAX
        IF ( Info(7)/=0 ) THEN
          rh = ABS(h)/Rwork(LHMAX)
          IF ( rh>1.0E0 ) h = h/rh
        END IF
        !
        CALL SDASTP(tn,Y,Yprime,Neq,RES,JAC,h,Rwork(lwt:lphi-1),Info(1),Idid,&
          Rwork(lphi:lpd-1),Rwork(LDELTA:le-1),Rwork(le:lwt-1),Rwork(lwm:),&
          Iwork,Rwork(LALPHA:LBETA-1),Rwork(LBETA:LGAMMA-1),Rwork(LGAMMA:LPSI-1),&
          Rwork(LPSI:LSIGMA-1),Rwork(LSIGMA:LDELTA-1),Rwork(LCJ),Rwork(LCJOLD),&
          Rwork(LHOLD),Rwork(LS),hmin,Rwork(LROUND),Iwork(LPHASE),&
          Iwork(LJCALC),Iwork(LK),Iwork(LKOLD),Iwork(LNS),Info(10),ntemp)
        !     MULTIPLY RTOL AND ATOL BY R AND RETURN
      ELSEIF ( Info(2)==1 ) THEN
        DO i = 1, Neq
          Rtol(i) = r*Rtol(i)
          Atol(i) = r*Atol(i)
        END DO
        Idid = -2
      ELSE
        Rtol(1) = r*Rtol(1)
        Atol(1) = r*Atol(1)
        Idid = -2
      END IF
    ELSE
      Idid = -1
    END IF
  END IF
  200 CONTINUE
  IF ( Idid<0 ) THEN
    !
    !-----------------------------------------------------------------------
    !     THIS BLOCK HANDLES ALL UNSUCCESSFUL
    !     RETURNS OTHER THAN FOR ILLEGAL INPUT.
    !-----------------------------------------------------------------------
    !
    itemp = -Idid
    SELECT CASE (itemp)
      CASE (2)
        !
        !     TOO MUCH ACCURACY FOR MACHINE PRECISION
        WRITE (xern3,'(1P,E15.6)') tn
        CALL XERMSG('SDASSL','AT T = '//xern3//&
          ' TOO MUCH ACCURACY REQUESTED FOR PRECISION OF MACHINE. RTOL AND ATOL WERE INCREASED TO APPROPRIATE VALUES',Idid,1)
      CASE (3)
        !
        !     WT(I) .LE. 0.0 FOR SOME I (NOT AT START OF PROBLEM)
        WRITE (xern3,'(1P,E15.6)') tn
        CALL XERMSG('SDASSL','AT T = '//xern3//&
          ' SOME ELEMENT OF WT HAS BECOME .LE. 0.0',Idid,1)
      CASE (4,5)
      CASE (6)
        !
        !     ERROR TEST FAILED REPEATEDLY OR WITH H=HMIN
        WRITE (xern3,'(1P,E15.6)') tn
        WRITE (xern4,'(1P,E15.6)') h
        CALL XERMSG('SDASSL','AT T = '//xern3//' AND STEPSIZE H = '//&
          xern4//' THE ERROR TEST FAILED REPEATEDLY OR WITH ABS(H)=HMIN',Idid,1)
      CASE (7)
        !
        !     CORRECTOR CONVERGENCE FAILED REPEATEDLY OR WITH H=HMIN
        WRITE (xern3,'(1P,E15.6)') tn
        WRITE (xern4,'(1P,E15.6)') h
        CALL XERMSG('SDASSL','AT T = '//xern3//' AND STEPSIZE H = '//xern4//&
          ' THE CORRECTOR FAILED TO CONVERGE REPEATEDLY OR WITH ABS(H)=HMIN',Idid,1)
      CASE (8)
        !
        !     THE ITERATION MATRIX IS SINGULAR
        WRITE (xern3,'(1P,E15.6)') tn
        WRITE (xern4,'(1P,E15.6)') h
        CALL XERMSG('SDASSL','AT T = '//xern3//' AND STEPSIZE H = '//&
          xern4//' THE ITERATION MATRIX IS SINGULAR',Idid,1)
      CASE (9)
        !
        !     CORRECTOR FAILURE PRECEDED BY ERROR TEST FAILURES.
        WRITE (xern3,'(1P,E15.6)') tn
        WRITE (xern4,'(1P,E15.6)') h
        CALL XERMSG('SDASSL','AT T = '//xern3//' AND STEPSIZE H = '//xern4//&
          ' THE CORRECTOR COULD NOT CONVERGE.  ALSO, THE ERROR TEST FAILED REPEATEDLY.',Idid,1)
      CASE (10)
        !
        !     CORRECTOR FAILURE BECAUSE IRES = -1
        WRITE (xern3,'(1P,E15.6)') tn
        WRITE (xern4,'(1P,E15.6)') h
        CALL XERMSG('SDASSL','AT T = '//xern3//' AND STEPSIZE H = '//xern4//&
          ' THE CORRECTOR COULD NOT CONVERGE BECAUSE IRES WAS EQUAL TO MINUS ONE',Idid,1)
      CASE (11)
        !
        !     FAILURE BECAUSE IRES = -2
        WRITE (xern3,'(1P,E15.6)') tn
        WRITE (xern4,'(1P,E15.6)') h
        CALL XERMSG('SDASSL','AT T = '//xern3//' AND STEPSIZE H = '//&
          xern4//' IRES WAS EQUAL TO MINUS TWO',Idid,1)
      CASE (12)
        !
        !     FAILED TO COMPUTE INITIAL YPRIME
        WRITE (xern3,'(1P,E15.6)') tn
        WRITE (xern4,'(1P,E15.6)') ho
        CALL XERMSG('SDASSL','AT T = '//xern3//' AND STEPSIZE H = '//&
          xern4//' THE INITIAL YPRIME COULD NOT BE COMPUTED',Idid,1)
      CASE DEFAULT
        !
        !     THE MAXIMUM NUMBER OF STEPS WAS TAKEN BEFORE
        !     REACHING TOUT
        WRITE (xern3,'(1P,E15.6)') tn
        CALL XERMSG('SDASSL','AT CURRENT T = '//xern3//&
          ' 500 STEPS TAKEN ON THIS CALL BEFORE REACHING TOUT',Idid,1)
    END SELECT
    !
    Info(1) = -1
    T = tn
    Rwork(LTN) = tn
    Rwork(LH) = h
    RETURN
    !
    !--------------------------------------------------------
    !     THIS BLOCK HANDLES THE CASE OF A SUCCESSFUL RETURN
    !     FROM SDASTP (IDID=1).  TEST FOR STOP CONDITIONS.
    !--------------------------------------------------------
    !
  ELSEIF ( Info(4)/=0 ) THEN
    IF ( Info(3)/=0 ) THEN
      IF ( (tn-Tout)*h>=0.0E0 ) THEN
        CALL SDATRP(tn,Tout,Y,Yprime,Neq,Iwork(LKOLD),Rwork(lphi),Rwork(LPSI))
        T = Tout
        Idid = 3
      ELSEIF ( ABS(tn-tstop)<=100.0E0*uround*(ABS(tn)+ABS(h)) ) THEN
        CALL SDATRP(tn,tstop,Y,Yprime,Neq,Iwork(LKOLD),Rwork(lphi),Rwork(LPSI))
        Idid = 2
        T = tstop
      ELSE
        T = tn
        Idid = 1
      END IF
    ELSEIF ( (tn-Tout)*h<0.0E0 ) THEN
      IF ( ABS(tn-tstop)<=100.0E0*uround*(ABS(tn)+ABS(h)) ) THEN
        CALL SDATRP(tn,tstop,Y,Yprime,Neq,Iwork(LKOLD),Rwork(lphi),&
          Rwork(LPSI))
        Idid = 2
        T = tstop
      ELSE
        tnext = tn + h
        IF ( (tnext-tstop)*h>0.0E0 ) h = tstop - tn
        GOTO 100
      END IF
    ELSE
      CALL SDATRP(tn,Tout,Y,Yprime,Neq,Iwork(LKOLD),Rwork(lphi),Rwork(LPSI))
      T = Tout
      Idid = 3
    END IF
  ELSEIF ( Info(3)/=0 ) THEN
    IF ( (tn-Tout)*h>=0.0E0 ) THEN
      CALL SDATRP(tn,Tout,Y,Yprime,Neq,Iwork(LKOLD),Rwork(lphi),Rwork(LPSI))
      Idid = 3
      T = Tout
    ELSE
      T = tn
      Idid = 1
    END IF
  ELSE
    IF ( (tn-Tout)*h<0.0E0 ) GOTO 100
    CALL SDATRP(tn,Tout,Y,Yprime,Neq,Iwork(LKOLD),Rwork(lphi),Rwork(LPSI))
    Idid = 3
    T = Tout
  END IF
  !
  !--------------------------------------------------------
  !     ALL SUCCESSFUL RETURNS FROM SDASSL ARE MADE FROM
  !     THIS BLOCK.
  !--------------------------------------------------------
  !
  300  Rwork(LTN) = tn
  Rwork(LH) = h
  RETURN
  !
  !-----------------------------------------------------------------------
  !     THIS BLOCK HANDLES ALL ERROR RETURNS DUE
  !     TO ILLEGAL INPUT, AS DETECTED BEFORE CALLING
  !     SDASTP. FIRST THE ERROR MESSAGE ROUTINE IS
  !     CALLED. IF THIS HAPPENS TWICE IN
  !     SUCCESSION, EXECUTION IS TERMINATED
  !
  !-----------------------------------------------------------------------
  400  CALL XERMSG('SDASSL',&
    'SOME ELEMENT OF INFO VECTOR IS NOT ZERO OR ONE',1,1)
  GOTO 1200
  !
  500  CALL XERMSG('SDASSL','SOME ELEMENT OF RTOL IS .LT. 0',6,1)
  GOTO 1200
  !
  600  CALL XERMSG('SDASSL','SOME ELEMENT OF ATOL IS .LT. 0',7,1)
  GOTO 1200
  !
  700  WRITE (xern3,'(1P,E15.6)') tstop
  WRITE (xern4,'(1P,E15.6)') Tout
  CALL XERMSG('SDASSL','INFO(4) = 1 AND TSTOP = '//xern3//&
    ' BEHIND TOUT = '//xern4,9,1)
  GOTO 1200
  !
  800  WRITE (xern3,'(1P,E15.6)') Tout
  WRITE (xern4,'(1P,E15.6)') T
  CALL XERMSG('SDASSL','TOUT = '//xern3//' BEHIND T = '//xern4,11,1)
  GOTO 1200
  !
  900  CALL XERMSG('SDASSL','SOME ELEMENT OF WT IS .LE. 0.0',13,1)
  GOTO 1200
  !
  1000 WRITE (xern3,'(1P,E15.6)') tstop
  WRITE (xern4,'(1P,E15.6)') T
  CALL XERMSG('SDASSL','INFO(4)=1 AND TSTOP = '//xern3//&
    ' BEHIND T = '//xern4,15,1)
  GOTO 1200
  !
  1100 WRITE (xern3,'(1P,E15.6)') Tout
  CALL XERMSG('SDASSL','TOUT = T = '//xern3,19,1)
  !
  1200 Idid = -33
  IF ( Info(1)==-1 ) CALL XERMSG('SDASSL',&
    'REPEATED OCCURRENCES OF ILLEGAL INPUT$$RUN TERMINATED. APPARENT INFINITE LOOP',&
    -999,2)
  !
  Info(1) = -1
  !-----------END OF SUBROUTINE SDASSL------------------------------------
END SUBROUTINE SDASSL
