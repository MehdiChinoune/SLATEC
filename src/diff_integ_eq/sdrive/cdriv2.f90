!** CDRIV2
SUBROUTINE CDRIV2(N,T,Y,F,Tout,Mstate,Nroot,Eps,Ewt,Mint,Work,Lenw,Iwork,&
    Leniw,G,Ierflg)
  IMPLICIT NONE
  !>
  !***
  !  The function of CDRIV2 is to solve N ordinary differential
  !            equations of the form dY(I)/dT = F(Y(I),T), given the
  !            initial conditions Y(I) = YI.  The program has options to
  !            allow the solution of both stiff and non-stiff differential
  !            equations.  CDRIV2 allows complex-valued differential
  !            equations.
  !***
  ! **Library:**   SLATEC (SDRIVE)
  !***
  ! **Category:**  I1A2, I1A1B
  !***
  ! **Type:**      COMPLEX (SDRIV2-S, DDRIV2-D, CDRIV2-C)
  !***
  ! **Keywords:**  COMPLEX VALUED, GEAR'S METHOD, INITIAL VALUE PROBLEMS,
  !             ODE, ORDINARY DIFFERENTIAL EQUATIONS, SDRIVE, STIFF
  !***
  ! **Author:**  Kahaner, D. K., (NIST)
  !             National Institute of Standards and Technology
  !             Gaithersburg, MD  20899
  !           Sutherland, C. D., (LANL)
  !             Mail Stop D466
  !             Los Alamos National Laboratory
  !             Los Alamos, NM  87545
  !***
  ! **Description:**
  !
  !  I.  PARAMETERS  .....................................................
  !
  !    The user should use parameter names in the call sequence of CDRIV2
  !    for those quantities whose value may be altered by CDRIV2.  The
  !    parameters in the call sequence are:
  !
  !    N      = (Input) The number of differential equations.
  !
  !    T      = (Real) The independent variable.  On input for the first
  !             call, T is the initial point.  On output, T is the point
  !             at which the solution is given.
  !
  !    Y      = (Complex) The vector of dependent variables.  Y is used as
  !             input on the first call, to set the initial values.  On
  !             output, Y is the computed solution vector.  This array Y
  !             is passed in the call sequence of the user-provided
  !             routines F and G.  Thus parameters required by F and G can
  !             be stored in this array in components N+1 and above.
  !             (Note: Changes by the user to the first N components of
  !             this array will take effect only after a restart, i.e.,
  !             after setting MSTATE to +1(-1).)
  !
  !    F      = A subroutine supplied by the user.  The name must be
  !             declared EXTERNAL in the user's calling program.  This
  !             subroutine is of the form:
  !                   SUBROUTINE F (N, T, Y, YDOT)
  !                   COMPLEX Y(*), YDOT(*)
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
  !             user's declaration in the program which calls CDRIV2.
  !             Thus the dimensioning of Y in F, while required by FORTRAN
  !             convention, does not actually allocate any storage.  When
  !             this subroutine is called, the first N components of Y are
  !             intermediate approximations to the solution components.
  !             The user should not alter these values.  Here YDOT is a
  !             vector of length N.  The user should only compute YDOT(I)
  !             for I from 1 to N.  Normally a return from F passes
  !             control back to  CDRIV2.  However, if the user would like
  !             to abort the calculation, i.e., return control to the
  !             program which calls CDRIV2, he should set N to zero.
  !             CDRIV2 will signal this by returning a value of MSTATE
  !             equal to +6(-6).  Altering the value of N in F has no
  !             effect on the value of N in the call sequence of CDRIV2.
  !
  !    TOUT   = (Input, Real) The point at which the solution is desired.
  !
  !    MSTATE = An integer describing the status of integration.  The user
  !             must initialize MSTATE to +1 or -1.  If MSTATE is
  !             positive, the routine will integrate past TOUT and
  !             interpolate the solution.  This is the most efficient
  !             mode.  If MSTATE is negative, the routine will adjust its
  !             internal step to reach TOUT exactly (useful if a
  !             singularity exists beyond TOUT.)  The meaning of the
  !             magnitude of MSTATE:
  !               1  (Input) Means the first call to the routine.  This
  !                  value must be set by the user.  On all subsequent
  !                  calls the value of MSTATE should be tested by the
  !                  user.  Unless CDRIV2 is to be reinitialized, only the
  !                  sign of MSTATE may be changed by the user.  (As a
  !                  convenience to the user who may wish to put out the
  !                  initial conditions, CDRIV2 can be called with
  !                  MSTATE=+1(-1), and TOUT=T.  In this case the program
  !                  will return with MSTATE unchanged, i.e.,
  !                  MSTATE=+1(-1).)
  !               2  (Output) Means a successful integration.  If a normal
  !                  continuation is desired (i.e., a further integration
  !                  in the same direction), simply advance TOUT and call
  !                  again.  All other parameters are automatically set.
  !               3  (Output)(Unsuccessful) Means the integrator has taken
  !                  1000 steps without reaching TOUT.  The user can
  !                  continue the integration by simply calling CDRIV2
  !                  again.  Other than an error in problem setup, the
  !                  most likely cause for this condition is trying to
  !                  integrate a stiff set of equations with the non-stiff
  !                  integrator option. (See description of MINT below.)
  !               4  (Output)(Unsuccessful) Means too much accuracy has
  !                  been requested.  EPS has been increased to a value
  !                  the program estimates is appropriate.  The user can
  !                  continue the integration by simply calling CDRIV2
  !                  again.
  !               5  (Output) A root was found at a point less than TOUT.
  !                  The user can continue the integration toward TOUT by
  !                  simply calling CDRIV2 again.
  !               6  (Output)(Unsuccessful) N has been set to zero in
  !                  SUBROUTINE F.
  !               7  (Output)(Unsuccessful) N has been set to zero in
  !                  FUNCTION G.  See description of G below.
  !               8  (Output)(Successful) For MSTATE negative, T is beyond
  !                  TOUT.  The solution was obtained by interpolation.
  !                  The user can continue the integration by simply
  !                  advancing TOUT and calling CDRIV2 again.
  !               9  (Output)(Unsuccessful) The solution could not be
  !                  obtained.  The value of IERFLG (see description
  !                  below) for a "Recoverable" situation indicates the
  !                  type of difficulty encountered: either an illegal
  !                  value for a parameter or an inability to continue the
  !                  solution.  For this condition the user should take
  !                  corrective action and reset MSTATE to +1(-1) before
  !                  calling CDRIV2 again.  Otherwise the program will
  !                  terminate the run.
  !
  !    NROOT  = (Input) The number of equations whose roots are desired.
  !             If NROOT is zero, the root search is not active.  This
  !             option is useful for obtaining output at points which are
  !             not known in advance, but depend upon the solution, e.g.,
  !             when some solution component takes on a specified value.
  !             The root search is carried out using the user-written
  !             function G (see description of G below.)  CDRIV2 attempts
  !             to find the value of T at which one of the equations
  !             changes sign.  CDRIV2 can find at most one root per
  !             equation per internal integration step, and will then
  !             return the solution either at TOUT or at a root, whichever
  !             occurs first in the direction of integration.  The initial
  !             point is never reported as a root.  The index of the
  !             equation whose root is being reported is stored in the
  !             sixth element of IWORK.
  !             NOTE: NROOT is never altered by this program.
  !
  !    EPS    = (Real) On input, the requested relative accuracy in all
  !             solution components.  EPS = 0 is allowed.  On output, the
  !             adjusted relative accuracy if the input value was too
  !             small.  The value of EPS should be set as large as is
  !             reasonable, because the amount of work done by CDRIV2
  !             increases as EPS decreases.
  !
  !    EWT    = (Input, Real) Problem zero, i.e., the smallest physically
  !             meaningful value for the solution.  This is used inter-
  !             nally to compute an array YWT(I) = MAX(ABS(Y(I)), EWT).
  !             One step error estimates divided by YWT(I) are kept less
  !             than EPS.  Setting EWT to zero provides pure relative
  !             error control.  However, setting EWT smaller than
  !             necessary can adversely affect the running time.
  !
  !    MINT   = (Input) The integration method flag.
  !               MINT = 1  Means the Adams methods, and is used for
  !                         non-stiff problems.
  !               MINT = 2  Means the stiff methods of Gear (i.e., the
  !                         backward differentiation formulas), and is
  !                         used for stiff problems.
  !               MINT = 3  Means the program dynamically selects the
  !                         Adams methods when the problem is non-stiff
  !                         and the Gear methods when the problem is
  !                         stiff.
  !             MINT may not be changed without restarting, i.e., setting
  !             the magnitude of MSTATE to 1.
  !
  !    WORK
  !    LENW   = (Input)
  !             WORK is an array of LENW complex words used
  !             internally for temporary storage.  The user must allocate
  !             space for this array in the calling program by a statement
  !             such as
  !                       COMPLEX WORK(...)
  !             The length of WORK should be at least
  !               16*N + 2*NROOT + 250         if MINT is 1, or
  !               N*N + 10*N + 2*NROOT + 250   if MINT is 2, or
  !               N*N + 17*N + 2*NROOT + 250   if MINT is 3,
  !             and LENW should be set to the value used.  The contents of
  !             WORK should not be disturbed between calls to CDRIV2.
  !
  !    IWORK
  !    LENIW  = (Input)
  !             IWORK is an integer array of length LENIW used internally
  !             for temporary storage.  The user must allocate space for
  !             this array in the calling program by a statement such as
  !                       INTEGER IWORK(...)
  !             The length of IWORK should be at least
  !               50      if MINT is 1, or
  !               N+50    if MINT is 2 or 3,
  !             and LENIW should be set to the value used.  The contents
  !             of IWORK should not be disturbed between calls to CDRIV2.
  !
  !    G      = A real FORTRAN function supplied by the user
  !             if NROOT is not 0.  In this case, the name must be
  !             declared EXTERNAL in the user's calling program.  G is
  !             repeatedly called with different values of IROOT to
  !             obtain the value of each of the NROOT equations for which
  !             a root is desired.  G is of the form:
  !                   REAL FUNCTION G (N, T, Y, IROOT)
  !                   COMPLEX Y(*)
  !                   GO TO (10, ...), IROOT
  !              10   G = ...
  !                     .
  !                     .
  !                   END (Sample)
  !             Here, Y is a vector of length at least N, whose first N
  !             components are the solution components at the point T.
  !             The user should not alter these values.  The actual length
  !             of Y is determined by the user's declaration in the
  !             program which calls CDRIV2.  Thus the dimensioning of Y in
  !             G, while required by FORTRAN convention, does not actually
  !             allocate any storage.  Normally a return from G passes
  !             control back to  CDRIV2.  However, if the user would like
  !             to abort the calculation, i.e., return control to the
  !             program which calls CDRIV2, he should set N to zero.
  !             CDRIV2 will signal this by returning a value of MSTATE
  !             equal to +7(-7).  In this case, the index of the equation
  !             being evaluated is stored in the sixth element of IWORK.
  !             Altering the value of N in G has no effect on the value of
  !             N in the call sequence of CDRIV2.
  !
  !    IERFLG = An error flag.  The error number associated with a
  !             diagnostic message (see Section II-A below) is the same as
  !             the corresponding value of IERFLG.  The meaning of IERFLG:
  !               0  The routine completed successfully. (No message is
  !                  issued.)
  !               3  (Warning) The number of steps required to reach TOUT
  !                  exceeds MXSTEP.
  !               4  (Warning) The value of EPS is too small.
  !              11  (Warning) For MSTATE negative, T is beyond TOUT.
  !                  The solution was obtained by interpolation.
  !              15  (Warning) The integration step size is below the
  !                  roundoff level of T.  (The program issues this
  !                  message as a warning but does not return control to
  !                  the user.)
  !              22  (Recoverable) N is not positive.
  !              23  (Recoverable) MINT is less than 1 or greater than 3 .
  !              26  (Recoverable) The magnitude of MSTATE is either 0 or
  !                  greater than 9 .
  !              27  (Recoverable) EPS is less than zero.
  !              32  (Recoverable) Insufficient storage has been allocated
  !                  for the WORK array.
  !              33  (Recoverable) Insufficient storage has been allocated
  !                  for the IWORK array.
  !              41  (Recoverable) The integration step size has gone
  !                  to zero.
  !              42  (Recoverable) The integration step size has been
  !                  reduced about 50 times without advancing the
  !                  solution.  The problem setup may not be correct.
  !             999  (Fatal) The magnitude of MSTATE is 9 .
  !
  !  II.  OTHER COMMUNICATION TO THE USER  ...............................
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
  !  III.  REMARKS  ......................................................
  !
  !    A. On any return from CDRIV2 all information necessary to continue
  !       the calculation is contained in the call sequence parameters,
  !       including the work arrays.  Thus it is possible to suspend one
  !       problem, integrate another, and then return to the first.
  !
  !    B. If this package is to be used in an overlay situation, the user
  !       must declare in the primary overlay the variables in the call
  !       sequence to CDRIV2.
  !
  !    C. When the routine G is not required, difficulties associated with
  !       an unsatisfied external can be avoided by using the name of the
  !       routine which calculates the right hand side of the differential
  !       equations in place of G in the call sequence of CDRIV2.
  !
  !  IV.  USAGE  .........................................................
  !
  !               PROGRAM SAMPLE
  !               EXTERNAL F
  !               PARAMETER(MINT = 1, NROOT = 0, N = ...,
  !              8          LENW = 16*N + 2*NROOT + 250, LENIW = 50)
  !         C                                 N is the number of equations
  !               COMPLEX WORK(LENW), Y(N)
  !               REAL EPS, EWT, T, TOUT
  !               INTEGER IWORK(LENIW)
  !               OPEN(FILE='TAPE6', UNIT=6, STATUS='NEW')
  !         C                                                Initial point
  !               T = 0.
  !         C                                       Set initial conditions
  !               DO 10 I = 1,N
  !          10     Y(I) = ...
  !               TOUT = T
  !               EWT = ...
  !               MSTATE = 1
  !               EPS = ...
  !          20   CALL CDRIV2 (N, T, Y, F, TOUT, MSTATE, NROOT, EPS, EWT,
  !              8             MINT, WORK, LENW, IWORK, LENIW, F, IERFLG)
  !         C                                 Next to last argument is not
  !         C                                    F if rootfinding is used.
  !               IF (MSTATE .GT. 2) STOP
  !               WRITE(6, 100) TOUT, (Y(I), I=1,N)
  !               TOUT = TOUT + 1.
  !               IF (TOUT .LE. 10.) GO TO 20
  !          100  FORMAT(...)
  !               END (Sample)
  !
  !***
  ! **References:**  C. W. Gear, Numerical Initial Value Problems in
  !                 Ordinary Differential Equations, Prentice-Hall, 1971.
  !***
  ! **Routines called:**  CDRIV3, XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   790601  DATE WRITTEN
  !   900329  Initial submission to SLATEC.

  EXTERNAL :: F, G
  COMPLEX Work(*), Y(*)
  REAL Eps, Ewt, ewtcom(1), G, hmax, T, Tout
  INTEGER Iwork(*)
  INTEGER Ierflg, ierror, Leniw, Lenw, Mint, miter, ml, Mstate, mu, mxord, N, &
    nde, Nroot, nstate, ntask
  CHARACTER intgr1*8
  INTEGER, PARAMETER :: IMPL = 0, MXSTEP = 1000
  !* FIRST EXECUTABLE STATEMENT  CDRIV2
  IF ( ABS(Mstate)==9 ) THEN
    Ierflg = 999
    CALL XERMSG('SLATEC','CDRIV2',&
      'Illegal input.  The magnitude of MSTATE IS 9 .',Ierflg,2)
    RETURN
  ELSEIF ( ABS(Mstate)==0.OR.ABS(Mstate)>9 ) THEN
    WRITE (intgr1,'(I8)') Mstate
    Ierflg = 26
    CALL XERMSG('SLATEC','CDRIV2',&
      'Illegal input.  The magnitude of MSTATE, '//intgr1//&
      ' is not in the range 1 to 8 .',Ierflg,1)
    Mstate = SIGN(9,Mstate)
    RETURN
  END IF
  IF ( Mint<1.OR.Mint>3 ) THEN
    WRITE (intgr1,'(I8)') Mint
    Ierflg = 23
    CALL XERMSG('SLATEC','CDRIV2',&
      'Illegal input.  Improper value for the integration method flag, '//intgr1//' .',Ierflg,1)
    Mstate = SIGN(9,Mstate)
    RETURN
  END IF
  IF ( Mstate>=0 ) THEN
    nstate = Mstate
    ntask = 1
  ELSE
    nstate = -Mstate
    ntask = 3
  END IF
  ewtcom(1) = Ewt
  IF ( Ewt/=0.E0 ) THEN
    ierror = 3
  ELSE
    ierror = 2
  END IF
  IF ( Mint==1 ) THEN
    miter = 0
    mxord = 12
  ELSEIF ( Mint==2 ) THEN
    miter = 2
    mxord = 5
  ELSEIF ( Mint==3 ) THEN
    miter = 2
    mxord = 12
  END IF
  hmax = 2.E0*ABS(Tout-T)
  CALL CDRIV3(N,T,Y,F,nstate,Tout,ntask,Nroot,Eps,ewtcom,ierror,Mint,miter,&
    IMPL,ml,mu,mxord,hmax,Work,Lenw,Iwork,Leniw,F,F,nde,MXSTEP,G,F,Ierflg)
  IF ( nstate<=7 ) THEN
    Mstate = SIGN(nstate,Mstate)
  ELSEIF ( nstate==11 ) THEN
    Mstate = SIGN(8,Mstate)
  ELSEIF ( nstate>11 ) THEN
    Mstate = SIGN(9,Mstate)
  END IF
END SUBROUTINE CDRIV2
