!*==DDRIV1.f90  processed by SPAG 6.72Dc at 11:01 on  6 Feb 2019
!DECK DDRIV1
      SUBROUTINE DDRIV1(N,T,Y,F,Tout,Mstate,Eps,Work,Lenw,Ierflg)
      IMPLICIT NONE
!*--DDRIV15
!*** Start of declarations inserted by SPAG
      REAL F
!*** End of declarations inserted by SPAG
!***BEGIN PROLOGUE  DDRIV1
!***PURPOSE  The function of DDRIV1 is to solve N (200 or fewer)
!            ordinary differential equations of the form
!            dY(I)/dT = F(Y(I),T), given the initial conditions
!            Y(I) = YI.  DDRIV1 uses double precision arithmetic.
!***LIBRARY   SLATEC (SDRIVE)
!***CATEGORY  I1A2, I1A1B
!***TYPE      DOUBLE PRECISION (SDRIV1-S, DDRIV1-D, CDRIV1-C)
!***KEYWORDS  DOUBLE PRECISION, GEAR'S METHOD, INITIAL VALUE PROBLEMS,
!             ODE, ORDINARY DIFFERENTIAL EQUATIONS, SDRIVE, STIFF
!***AUTHOR  Kahaner, D. K., (NIST)
!             National Institute of Standards and Technology
!             Gaithersburg, MD  20899
!           Sutherland, C. D., (LANL)
!             Mail Stop D466
!             Los Alamos National Laboratory
!             Los Alamos, NM  87545
!***DESCRIPTION
!
!   Version 92.1
!
!  I.  CHOOSING THE CORRECT ROUTINE  ...................................
!
!     SDRIV
!     DDRIV
!     CDRIV
!           These are the generic names for three packages for solving
!           initial value problems for ordinary differential equations.
!           SDRIV uses single precision arithmetic.  DDRIV uses double
!           precision arithmetic.  CDRIV allows complex-valued
!           differential equations, integrated with respect to a single,
!           real, independent variable.
!
!    As an aid in selecting the proper program, the following is a
!    discussion of the important options or restrictions associated with
!    each program:
!
!      A. DDRIV1 should be tried first for those routine problems with
!         no more than 200 differential equations (DDRIV2 and DDRIV3
!         have no such restriction.)  Internally this routine has two
!         important technical defaults:
!           1. Numerical approximation of the Jacobian matrix of the
!              right hand side is used.
!           2. The stiff solver option is used.
!         Most users of DDRIV1 should not have to concern themselves
!         with these details.
!
!      B. DDRIV2 should be considered for those problems for which
!         DDRIV1 is inadequate.  For example, DDRIV1 may have difficulty
!         with problems having zero initial conditions and zero
!         derivatives.  In this case DDRIV2, with an appropriate value
!         of the parameter EWT, should perform more efficiently.  DDRIV2
!         provides three important additional options:
!           1. The nonstiff equation solver (as well as the stiff
!              solver) is available.
!           2. The root-finding option is available.
!           3. The program can dynamically select either the non-stiff
!              or the stiff methods.
!         Internally this routine also defaults to the numerical
!         approximation of the Jacobian matrix of the right hand side.
!
!      C. DDRIV3 is the most flexible, and hence the most complex, of
!         the programs.  Its important additional features include:
!           1. The ability to exploit band structure in the Jacobian
!              matrix.
!           2. The ability to solve some implicit differential
!              equations, i.e., those having the form:
!                   A(Y,T)*dY/dT = F(Y,T).
!           3. The option of integrating in the one step mode.
!           4. The option of allowing the user to provide a routine
!              which computes the analytic Jacobian matrix of the right
!              hand side.
!           5. The option of allowing the user to provide a routine
!              which does all the matrix algebra associated with
!              corrections to the solution components.
!
!  II.  PARAMETERS  ....................................................
!
!       (REMEMBER--To run DDRIV1 correctly in double precision, ALL
!       non-integer arguments in the call sequence, including
!       arrays, MUST be declared double precision.)
!
!    The user should use parameter names in the call sequence of DDRIV1
!    for those quantities whose value may be altered by DDRIV1.  The
!    parameters in the call sequence are:
!
!    N      = (Input) The number of differential equations, N .LE. 200
!
!    T      = The independent variable.  On input for the first call, T
!             is the initial point.  On output, T is the point at which
!             the solution is given.
!
!    Y      = The vector of dependent variables.  Y is used as input on
!             the first call, to set the initial values.  On output, Y
!             is the computed solution vector.  This array Y is passed
!             in the call sequence of the user-provided routine F.  Thus
!             parameters required by F can be stored in this array in
!             components N+1 and above.  (Note: Changes by the user to
!             the first N components of this array will take effect only
!             after a restart, i.e., after setting MSTATE to +1(-1).)
!
!    F      = A subroutine supplied by the user.  The name must be
!             declared EXTERNAL in the user's calling program.  This
!             subroutine is of the form:
!                   SUBROUTINE F (N, T, Y, YDOT)
!                   DOUBLE PRECISION Y(*), YDOT(*)
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
!             user's declaration in the program which calls DDRIV1.
!             Thus the dimensioning of Y in F, while required by FORTRAN
!             convention, does not actually allocate any storage.  When
!             this subroutine is called, the first N components of Y are
!             intermediate approximations to the solution components.
!             The user should not alter these values.  Here YDOT is a
!             vector of length N.  The user should only compute YDOT(I)
!             for I from 1 to N.  Normally a return from F passes
!             control back to  DDRIV1.  However, if the user would like
!             to abort the calculation, i.e., return control to the
!             program which calls DDRIV1, he should set N to zero.
!             DDRIV1 will signal this by returning a value of MSTATE
!             equal to +5(-5).  Altering the value of N in F has no
!             effect on the value of N in the call sequence of DDRIV1.
!
!    TOUT   = (Input) The point at which the solution is desired.
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
!                  user.  Unless DDRIV1 is to be reinitialized, only the
!                  sign of MSTATE may be changed by the user.  (As a
!                  convenience to the user who may wish to put out the
!                  initial conditions, DDRIV1 can be called with
!                  MSTATE=+1(-1), and TOUT=T.  In this case the program
!                  will return with MSTATE unchanged, i.e.,
!                  MSTATE=+1(-1).)
!               2  (Output) Means a successful integration.  If a normal
!                  continuation is desired (i.e., a further integration
!                  in the same direction), simply advance TOUT and call
!                  again.  All other parameters are automatically set.
!               3  (Output)(Unsuccessful) Means the integrator has taken
!                  1000 steps without reaching TOUT.  The user can
!                  continue the integration by simply calling DDRIV1
!                  again.
!               4  (Output)(Unsuccessful) Means too much accuracy has
!                  been requested.  EPS has been increased to a value
!                  the program estimates is appropriate.  The user can
!                  continue the integration by simply calling DDRIV1
!                  again.
!               5  (Output)(Unsuccessful) N has been set to zero in
!                  SUBROUTINE F.
!               6  (Output)(Successful) For MSTATE negative, T is beyond
!                  TOUT.  The solution was obtained by interpolation.
!                  The user can continue the integration by simply
!                  advancing TOUT and calling DDRIV1 again.
!               7  (Output)(Unsuccessful) The solution could not be
!                  obtained.  The value of IERFLG (see description
!                  below) for a "Recoverable" situation indicates the
!                  type of difficulty encountered: either an illegal
!                  value for a parameter or an inability to continue the
!                  solution.  For this condition the user should take
!                  corrective action and reset MSTATE to +1(-1) before
!                  calling DDRIV1 again.  Otherwise the program will
!                  terminate the run.
!
!    EPS    = On input, the requested relative accuracy in all solution
!             components.  On output, the adjusted relative accuracy if
!             the input value was too small.  The value of EPS should be
!             set as large as is reasonable, because the amount of work
!             done by DDRIV1 increases as EPS decreases.
!
!    WORK
!    LENW   = (Input)
!             WORK is an array of LENW double precision words used
!             internally for temporary storage.  The user must allocate
!             space for this array in the calling program by a statement
!             such as
!                       DOUBLE PRECISION WORK(...)
!             The length of WORK should be at least N*N + 11*N + 300
!             and LENW should be set to the value used.  The contents of
!             WORK should not be disturbed between calls to DDRIV1.
!
!    IERFLG = An error flag.  The error number associated with a
!             diagnostic message (see Section IV-A below) is the same as
!             the corresponding value of IERFLG.  The meaning of IERFLG:
!               0  The routine completed successfully. (No message is
!                  issued.)
!               3  (Warning) The number of steps required to reach TOUT
!                  exceeds 1000 .
!               4  (Warning) The value of EPS is too small.
!              11  (Warning) For MSTATE negative, T is beyond TOUT.
!                  The solution was obtained by interpolation.
!              15  (Warning) The integration step size is below the
!                  roundoff level of T.  (The program issues this
!                  message as a warning but does not return control to
!                  the user.)
!              21  (Recoverable) N is greater than 200 .
!              22  (Recoverable) N is not positive.
!              26  (Recoverable) The magnitude of MSTATE is either 0 or
!                  greater than 7 .
!              27  (Recoverable) EPS is less than zero.
!              32  (Recoverable) Insufficient storage has been allocated
!                  for the WORK array.
!              41  (Recoverable) The integration step size has gone
!                  to zero.
!              42  (Recoverable) The integration step size has been
!                  reduced about 50 times without advancing the
!                  solution.  The problem setup may not be correct.
!             999  (Fatal) The magnitude of MSTATE is 7 .
!
!  III.  USAGE  ........................................................
!
!                PROGRAM SAMPLE
!                EXTERNAL F
!                DOUBLE PRECISION ALFA, EPS, T, TOUT
!          C                                N is the number of equations
!                PARAMETER(ALFA = 1.D0, N = 3, LENW = N*N + 11*N + 300)
!                DOUBLE PRECISION WORK(LENW), Y(N+1)
!          C                                               Initial point
!                T = 0.00001D0
!          C                                      Set initial conditions
!                Y(1) = 10.D0
!                Y(2) = 0.D0
!                Y(3) = 10.D0
!          C                                              Pass parameter
!                Y(4) = ALFA
!                TOUT = T
!                MSTATE = 1
!                EPS = .001D0
!           10   CALL DDRIV1 (N, T, Y, F, TOUT, MSTATE, EPS, WORK, LENW,
!               8             IERFLG)
!                IF (MSTATE .GT. 2) STOP
!                WRITE(*, '(4E12.3)') TOUT, (Y(I), I=1,3)
!                TOUT = 10.D0*TOUT
!                IF (TOUT .LT. 50.D0) GO TO 10
!                END
!
!                SUBROUTINE F (N, T, Y, YDOT)
!                DOUBLE PRECISION ALFA, T, Y(*), YDOT(*)
!                ALFA = Y(N+1)
!                YDOT(1) = 1.D0 + ALFA*(Y(2) - Y(1)) - Y(1)*Y(3)
!                YDOT(2) = ALFA*(Y(1) - Y(2)) - Y(2)*Y(3)
!                YDOT(3) = 1.D0 - Y(3)*(Y(1) + Y(2))
!                END
!
!  IV.  OTHER COMMUNICATION TO THE USER  ...............................
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
!    B. The number of evaluations of the right hand side can be found
!       in the WORK array in the location determined by:
!            LENW - (N + 50) + 4
!
!  V.  REMARKS  ........................................................
!
!    For other information, see Section IV of the writeup for DDRIV3.
!
!***REFERENCES  C. W. Gear, Numerical Initial Value Problems in
!                 Ordinary Differential Equations, Prentice-Hall, 1971.
!***ROUTINES CALLED  DDRIV3, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   790601  DATE WRITTEN
!   900329  Initial submission to SLATEC.
!***END PROLOGUE  DDRIV1
      EXTERNAL F
      DOUBLE PRECISION Eps , ewtcom(1) , hmax , T , Tout , Work(*) , Y(*)
      INTEGER i , IDLIW , Ierflg , IERROR , IMPL , leniw , Lenw , lenwcm , 
     &        lnwchk , MINT , MITER , ml , Mstate , mu , MXN , MXORD , MXSTEP , 
     &        N , nde , NROOT , nstate , ntask
      PARAMETER (MXN=200,IDLIW=50)
      INTEGER iwork(IDLIW+MXN)
      CHARACTER intgr1*8
      PARAMETER (NROOT=0,IERROR=2,MINT=2,MITER=2,IMPL=0,MXORD=5,MXSTEP=1000)
      DATA ewtcom(1)/1.D0/
!***FIRST EXECUTABLE STATEMENT  DDRIV1
      IF ( ABS(Mstate)==0.OR.ABS(Mstate)>7 ) THEN
        WRITE (intgr1,'(I8)') Mstate
        Ierflg = 26
        CALL XERMSG('SLATEC','DDRIV1',
     &              'Illegal input.  The magnitude of MSTATE, '//intgr1//
     &              ', is not in the range 1 to 6 .',Ierflg,1)
        Mstate = SIGN(7,Mstate)
        RETURN
      ELSEIF ( ABS(Mstate)==7 ) THEN
        Ierflg = 999
        CALL XERMSG('SLATEC','DDRIV1',
     &              'Illegal input.  The magnitude of MSTATE is 7 .',Ierflg,2)
        RETURN
      ENDIF
      IF ( N>MXN ) THEN
        WRITE (intgr1,'(I8)') N
        Ierflg = 21
        CALL XERMSG('SLATEC','DDRIV1',
     &              'Illegal input.  The number of equations, '//intgr1//
     &              ', is greater than the maximum allowed: 200 .',Ierflg,1)
        Mstate = SIGN(7,Mstate)
        RETURN
      ENDIF
      IF ( Mstate>0 ) THEN
        nstate = Mstate
        ntask = 1
      ELSE
        nstate = -Mstate
        ntask = 3
      ENDIF
      hmax = 2.D0*ABS(Tout-T)
      leniw = N + IDLIW
      lenwcm = Lenw - leniw
      IF ( lenwcm<(N*N+10*N+250) ) THEN
        lnwchk = N*N + 10*N + 250 + leniw
        WRITE (intgr1,'(I8)') lnwchk
        Ierflg = 32
        CALL XERMSG('SLATEC','DDRIV1',
     &              'Insufficient storage allocated for the work array.  '//
     &              'The required storage is at least '//intgr1//' .',Ierflg,1)
        Mstate = SIGN(7,Mstate)
        RETURN
      ENDIF
      IF ( nstate/=1 ) THEN
        DO i = 1 , leniw
          iwork(i) = Work(i+lenwcm)
        ENDDO
      ENDIF
      CALL DDRIV3(N,T,Y,F,nstate,Tout,ntask,NROOT,Eps,ewtcom,IERROR,MINT,MITER,
     &            IMPL,ml,mu,MXORD,hmax,Work,lenwcm,iwork,leniw,F,F,nde,MXSTEP,
     &            F,F,Ierflg)
      DO i = 1 , leniw
        Work(i+lenwcm) = iwork(i)
      ENDDO
      IF ( nstate<=4 ) THEN
        Mstate = SIGN(nstate,Mstate)
      ELSEIF ( nstate==6 ) THEN
        Mstate = SIGN(5,Mstate)
      ELSEIF ( Ierflg==11 ) THEN
        Mstate = SIGN(6,Mstate)
      ELSEIF ( Ierflg>11 ) THEN
        Mstate = SIGN(7,Mstate)
      ENDIF
      END SUBROUTINE DDRIV1
