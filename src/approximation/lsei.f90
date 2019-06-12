!** LSEI
SUBROUTINE LSEI(W,Mdw,Me,Ma,Mg,N,Prgopt,X,Rnorme,Rnorml,Mode,Ws,Ip)
  !>
  !  Solve a linearly constrained least squares problem with
  !            equality and inequality constraints, and optionally compute
  !            a covariance matrix.
  !***
  ! **Library:**   SLATEC
  !***
  ! **Category:**  K1A2A, D9
  !***
  ! **Type:**      SINGLE PRECISION (LSEI-S, DLSEI-D)
  !***
  ! **Keywords:**  CONSTRAINED LEAST SQUARES, CURVE FITTING, DATA FITTING,
  !             EQUALITY CONSTRAINTS, INEQUALITY CONSTRAINTS,
  !             QUADRATIC PROGRAMMING
  !***
  ! **Author:**  Hanson, R. J., (SNLA)
  !           Haskell, K. H., (SNLA)
  !***
  ! **Description:**
  !
  !     Abstract
  !
  !     This subprogram solves a linearly constrained least squares
  !     problem with both equality and inequality constraints, and, if the
  !     user requests, obtains a covariance matrix of the solution
  !     parameters.
  !
  !     Suppose there are given matrices E, A and G of respective
  !     dimensions ME by N, MA by N and MG by N, and vectors F, B and H of
  !     respective lengths ME, MA and MG.  This subroutine solves the
  !     linearly constrained least squares problem
  !
  !                   EX = F, (E ME by N) (equations to be exactly
  !                                       satisfied)
  !                   AX = B, (A MA by N) (equations to be
  !                                       approximately satisfied,
  !                                       least squares sense)
  !                   GX .GE. H,(G MG by N) (inequality constraints)
  !
  !     The inequalities GX .GE. H mean that every component of the
  !     product GX must be .GE. the corresponding component of H.
  !
  !     In case the equality constraints cannot be satisfied, a
  !     generalized inverse solution residual vector length is obtained
  !     for F-EX.  This is the minimal length possible for F-EX.
  !
  !     Any values ME .GE. 0, MA .GE. 0, or MG .GE. 0 are permitted.  The
  !     rank of the matrix E is estimated during the computation.  We call
  !     this value KRANKE.  It is an output parameter in IP(1) defined
  !     below.  Using a generalized inverse solution of EX=F, a reduced
  !     least squares problem with inequality constraints is obtained.
  !     The tolerances used in these tests for determining the rank
  !     of E and the rank of the reduced least squares problem are
  !     given in Sandia Tech. Rept. SAND-78-1290.  They can be
  !     modified by the user if new values are provided in
  !     the option list of the array PRGOPT(*).
  !
  !     The user must dimension all arrays appearing in the call list..
  !     W(MDW,N+1),PRGOPT(*),X(N),WS(2*(ME+N)+K+(MG+2)*(N+7)),IP(MG+2*N+2)
  !     where K=MAX(MA+MG,N).  This allows for a solution of a range of
  !     problems in the given working space.  The dimension of WS(*)
  !     given is a necessary overestimate.  Once a particular problem
  !     has been run, the output parameter IP(3) gives the actual
  !     dimension required for that problem.
  !
  !     The parameters for LSEI( ) are
  !
  !     Input..
  !
  !     W(*,*),MDW,   The array W(*,*) is doubly subscripted with
  !     ME,MA,MG,N    first dimensioning parameter equal to MDW.
  !                   For this discussion let us call M = ME+MA+MG.  Then
  !                   MDW must satisfy MDW .GE. M.  The condition
  !                   MDW .LT. M is an error.
  !
  !                   The array W(*,*) contains the matrices and vectors
  !
  !                                  (E  F)
  !                                  (A  B)
  !                                  (G  H)
  !
  !                   in rows and columns 1,...,M and 1,...,N+1
  !                   respectively.
  !
  !                   The integers ME, MA, and MG are the
  !                   respective matrix row dimensions
  !                   of E, A and G.  Each matrix has N columns.
  !
  !     PRGOPT(*)    This real-valued array is the option vector.
  !                  If the user is satisfied with the nominal
  !                  subprogram features set
  !
  !                  PRGOPT(1)=1 (or PRGOPT(1)=1.0)
  !
  !                  Otherwise PRGOPT(*) is a linked list consisting of
  !                  groups of data of the following form
  !
  !                  LINK
  !                  KEY
  !                  DATA SET
  !
  !                  The parameters LINK and KEY are each one word.
  !                  The DATA SET can be comprised of several words.
  !                  The number of items depends on the value of KEY.
  !                  The value of LINK points to the first
  !                  entry of the next group of data within
  !                  PRGOPT(*).  The exception is when there are
  !                  no more options to change.  In that
  !                  case, LINK=1 and the values KEY and DATA SET
  !                  are not referenced.  The general layout of
  !                  PRGOPT(*) is as follows.
  !
  !               ...PRGOPT(1) = LINK1 (link to first entry of next group)
  !               .  PRGOPT(2) = KEY1 (key to the option change)
  !               .  PRGOPT(3) = data value (data value for this change)
  !               .       .
  !               .       .
  !               .       .
  !               ...PRGOPT(LINK1)   = LINK2 (link to the first entry of
  !               .                       next group)
  !               .  PRGOPT(LINK1+1) = KEY2 (key to the option change)
  !               .  PRGOPT(LINK1+2) = data value
  !               ...     .
  !               .       .
  !               .       .
  !               ...PRGOPT(LINK) = 1 (no more options to change)
  !
  !                  Values of LINK that are nonpositive are errors.
  !                  A value of LINK .GT. NLINK=100000 is also an error.
  !                  This helps prevent using invalid but positive
  !                  values of LINK that will probably extend
  !                  beyond the program limits of PRGOPT(*).
  !                  Unrecognized values of KEY are ignored.  The
  !                  order of the options is arbitrary and any number
  !                  of options can be changed with the following
  !                  restriction.  To prevent cycling in the
  !                  processing of the option array, a count of the
  !                  number of options changed is maintained.
  !                  Whenever this count exceeds NOPT=1000, an error
  !                  message is printed and the subprogram returns.
  !
  !                  Options..
  !
  !                  KEY=1
  !                         Compute in W(*,*) the N by N
  !                  covariance matrix of the solution variables
  !                  as an output parameter.  Nominally the
  !                  covariance matrix will not be computed.
  !                  (This requires no user input.)
  !                  The data set for this option is a single value.
  !                  It must be nonzero when the covariance matrix
  !                  is desired.  If it is zero, the covariance
  !                  matrix is not computed.  When the covariance matrix
  !                  is computed, the first dimensioning parameter
  !                  of the array W(*,*) must satisfy MDW .GE. MAX(M,N).
  !
  !                  KEY=10
  !                         Suppress scaling of the inverse of the
  !                  normal matrix by the scale factor RNORM**2/
  !                  MAX(1, no. of degrees of freedom).  This option
  !                  only applies when the option for computing the
  !                  covariance matrix (KEY=1) is used.  With KEY=1 and
  !                  KEY=10 used as options the unscaled inverse of the
  !                  normal matrix is returned in W(*,*).
  !                  The data set for this option is a single value.
  !                  When it is nonzero no scaling is done.  When it is
  !                  zero scaling is done.  The nominal case is to do
  !                  scaling so if option (KEY=1) is used alone, the
  !                  matrix will be scaled on output.
  !
  !                  KEY=2
  !                         Scale the nonzero columns of the
  !                         entire data matrix.
  !                  (E)
  !                  (A)
  !                  (G)
  !
  !                  to have length one.  The data set for this
  !                  option is a single value.  It must be
  !                  nonzero if unit length column scaling
  !                  is desired.
  !
  !                  KEY=3
  !                         Scale columns of the entire data matrix
  !                  (E)
  !                  (A)
  !                  (G)
  !
  !                  with a user-provided diagonal matrix.
  !                  The data set for this option consists
  !                  of the N diagonal scaling factors, one for
  !                  each matrix column.
  !
  !                  KEY=4
  !                         Change the rank determination tolerance for
  !                  the equality constraint equations from
  !                  the nominal value of SQRT(SRELPR).  This quantity can
  !                  be no smaller than SRELPR, the arithmetic-
  !                  storage precision.  The quantity SRELPR is the
  !                  largest positive number such that T=1.+SRELPR
  !                  satisfies T .EQ. 1.  The quantity used
  !                  here is internally restricted to be at
  !                  least SRELPR.  The data set for this option
  !                  is the new tolerance.
  !
  !                  KEY=5
  !                         Change the rank determination tolerance for
  !                  the reduced least squares equations from
  !                  the nominal value of SQRT(SRELPR).  This quantity can
  !                  be no smaller than SRELPR, the arithmetic-
  !                  storage precision.  The quantity used
  !                  here is internally restricted to be at
  !                  least SRELPR.  The data set for this option
  !                  is the new tolerance.
  !
  !                  For example, suppose we want to change
  !                  the tolerance for the reduced least squares
  !                  problem, compute the covariance matrix of
  !                  the solution parameters, and provide
  !                  column scaling for the data matrix.  For
  !                  these options the dimension of PRGOPT(*)
  !                  must be at least N+9.  The Fortran statements
  !                  defining these options would be as follows:
  !
  !                  PRGOPT(1)=4 (link to entry 4 in PRGOPT(*))
  !                  PRGOPT(2)=1 (covariance matrix key)
  !                  PRGOPT(3)=1 (covariance matrix wanted)
  !
  !                  PRGOPT(4)=7 (link to entry 7 in PRGOPT(*))
  !                  PRGOPT(5)=5 (least squares equas.  tolerance key)
  !                  PRGOPT(6)=... (new value of the tolerance)
  !
  !                  PRGOPT(7)=N+9 (link to entry N+9 in PRGOPT(*))
  !                  PRGOPT(8)=3 (user-provided column scaling key)
  !
  !                  CALL SCOPY (N, D, 1, PRGOPT(9), 1)  (Copy the N
  !                    scaling factors from the user array D(*)
  !                    to PRGOPT(9)-PRGOPT(N+8))
  !
  !                  PRGOPT(N+9)=1 (no more options to change)
  !
  !                  The contents of PRGOPT(*) are not modified
  !                  by the subprogram.
  !                  The options for WNNLS( ) can also be included
  !                  in this array.  The values of KEY recognized
  !                  by WNNLS( ) are 6, 7 and 8.  Their functions
  !                  are documented in the usage instructions for
  !                  subroutine WNNLS( ).  Normally these options
  !                  do not need to be modified when using LSEI( ).
  !
  !     IP(1),       The amounts of working storage actually
  !     IP(2)        allocated for the working arrays WS(*) and
  !                  IP(*), respectively.  These quantities are
  !                  compared with the actual amounts of storage
  !                  needed by LSEI( ).  Insufficient storage
  !                  allocated for either WS(*) or IP(*) is an
  !                  error.  This feature was included in LSEI( )
  !                  because miscalculating the storage formulas
  !                  for WS(*) and IP(*) might very well lead to
  !                  subtle and hard-to-find execution errors.
  !
  !                  The length of WS(*) must be at least
  !
  !                  LW = 2*(ME+N)+K+(MG+2)*(N+7)
  !
  !                  where K = max(MA+MG,N)
  !                  This test will not be made if IP(1).LE.0.
  !
  !                  The length of IP(*) must be at least
  !
  !                  LIP = MG+2*N+2
  !                  This test will not be made if IP(2).LE.0.
  !
  !     Output..
  !
  !     X(*),RNORME,  The array X(*) contains the solution parameters
  !     RNORML        if the integer output flag MODE = 0 or 1.
  !                   The definition of MODE is given directly below.
  !                   When MODE = 0 or 1, RNORME and RNORML
  !                   respectively contain the residual vector
  !                   Euclidean lengths of F - EX and B - AX.  When
  !                   MODE=1 the equality constraint equations EX=F
  !                   are contradictory, so RNORME .NE. 0.  The residual
  !                   vector F-EX has minimal Euclidean length.  For
  !                   MODE .GE. 2, none of these parameters is defined.
  !
  !     MODE          Integer flag that indicates the subprogram
  !                   status after completion.  If MODE .GE. 2, no
  !                   solution has been computed.
  !
  !                   MODE =
  !
  !                   0  Both equality and inequality constraints
  !                      are compatible and have been satisfied.
  !
  !                   1  Equality constraints are contradictory.
  !                      A generalized inverse solution of EX=F was used
  !                      to minimize the residual vector length F-EX.
  !                      In this sense, the solution is still meaningful.
  !
  !                   2  Inequality constraints are contradictory.
  !
  !                   3  Both equality and inequality constraints
  !                      are contradictory.
  !
  !                   The following interpretation of
  !                   MODE=1,2 or 3 must be made.  The
  !                   sets consisting of all solutions
  !                   of the equality constraints EX=F
  !                   and all vectors satisfying GX .GE. H
  !                   have no points in common.  (In
  !                   particular this does not say that
  !                   each individual set has no points
  !                   at all, although this could be the
  !                   case.)
  !
  !                   4  Usage error occurred.  The value
  !                      of MDW is .LT. ME+MA+MG, MDW is
  !                      .LT. N and a covariance matrix is
  !                      requested, or the option vector
  !                      PRGOPT(*) is not properly defined,
  !                      or the lengths of the working arrays
  !                      WS(*) and IP(*), when specified in
  !                      IP(1) and IP(2) respectively, are not
  !                      long enough.
  !
  !     W(*,*)        The array W(*,*) contains the N by N symmetric
  !                   covariance matrix of the solution parameters,
  !                   provided this was requested on input with
  !                   the option vector PRGOPT(*) and the output
  !                   flag is returned with MODE = 0 or 1.
  !
  !     IP(*)         The integer working array has three entries
  !                   that provide rank and working array length
  !                   information after completion.
  !
  !                      IP(1) = rank of equality constraint
  !                              matrix.  Define this quantity
  !                              as KRANKE.
  !
  !                      IP(2) = rank of reduced least squares
  !                              problem.
  !
  !                      IP(3) = the amount of storage in the
  !                              working array WS(*) that was
  !                              actually used by the subprogram.
  !                              The formula given above for the length
  !                              of WS(*) is a necessary overestimate.
  !                              If exactly the same problem matrices
  !                              are used in subsequent executions,
  !                              the declared dimension of WS(*) can
  !                              be reduced to this output value.
  !     User Designated
  !     Working Arrays..
  !
  !     WS(*),IP(*)              These are respectively type real
  !                              and type integer working arrays.
  !                              Their required minimal lengths are
  !                              given above.
  !
  !***
  ! **References:**  K. H. Haskell and R. J. Hanson, An algorithm for
  !                 linear least squares problems with equality and
  !                 nonnegativity constraints, Report SAND77-0552, Sandia
  !                 Laboratories, June 1978.
  !               K. H. Haskell and R. J. Hanson, Selected algorithms for
  !                 the linearly constrained least squares problem - a
  !                 users guide, Report SAND78-1290, Sandia Laboratories,
  !                 August 1979.
  !               K. H. Haskell and R. J. Hanson, An algorithm for
  !                 linear least squares problems with equality and
  !                 nonnegativity constraints, Mathematical Programming
  !                 21 (1981), pp. 98-118.
  !               R. J. Hanson and K. H. Haskell, Two algorithms for the
  !                 linearly constrained least squares problem, ACM
  !                 Transactions on Mathematical Software, September 1982.
  !***
  ! **Routines called:**  H12, LSI, R1MACH, SASUM, SAXPY, SCOPY, SDOT, SNRM2,
  !                    SSCAL, SSWAP, XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   790701  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890618  Completely restructured and extensively revised (WRB & RWC)
  !   890831  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900510  Convert XERRWV calls to XERMSG calls.  (RWC)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  USE service, ONLY : XERMSG, R1MACH
  USE blas, ONLY : SAXPY, SSWAP
  USE linear, ONLY : H12
  INTEGER :: Ma, Mdw, Me, Mg, Mode, N, Ip(Mg+2*N+2)
  REAL(SP) :: Prgopt(:), Rnorme, Rnorml, W(Mdw,N+1), Ws(*), X(N)
  !
  REAL(SP) :: enorm, fnorm, gam, rb, rn, rnmax, sizee, sn, snmax, t, tau, uj, up, vj, &
    xnorm, xnrme
  INTEGER :: i, imax, j, jp1, k, key, kranke, last, lchk, link, m,  mapke1, mdeqc, &
    mend, mep1, n1, n2, next, nlink, nopt, np1, ntimes
  LOGICAL :: cov
  CHARACTER(8) :: xern1, xern2, xern3, xern4
  !
  REAL(SP), PARAMETER :: srelpr = R1MACH(4)
  !* FIRST EXECUTABLE STATEMENT  LSEI
  !
  !     Set the nominal tolerance used in the code for the equality
  !     constraint equations.
  !
  tau = SQRT(srelpr)
  !
  !     Check that enough storage was allocated in WS(*) and IP(*).
  !
  Mode = 4
  IF ( MIN(N,Me,Ma,Mg)<0 ) THEN
    WRITE (xern1,'(I8)') N
    WRITE (xern2,'(I8)') Me
    WRITE (xern3,'(I8)') Ma
    WRITE (xern4,'(I8)') Mg
    CALL XERMSG('LSEI', &
      'ALL OF THE VARIABLES N, ME, MA, MG MUST BE .GE. 0$$ENTERED ROUTINE WITH$$N  = '//&
      xern1//'$$ME = '//xern2//'$$MA = '//xern3//'$$MG = '//xern4,2,1)
    RETURN
  END IF
  !
  IF ( Ip(1)>0 ) THEN
    lchk = 2*(Me+N) + MAX(Ma+Mg,N) + (Mg+2)*(N+7)
    IF ( Ip(1)<lchk ) THEN
      WRITE (xern1,'(I8)') lchk
      CALL XERMSG('LSEI','INSUFFICIENT STORAGE ALLOCATED FOR WS(*), NEED LW = '//xern1,2,1)
      RETURN
    END IF
  END IF
  !
  IF ( Ip(2)>0 ) THEN
    lchk = Mg + 2*N + 2
    IF ( Ip(2)<lchk ) THEN
      WRITE (xern1,'(I8)') lchk
      CALL XERMSG('LSEI','INSUFFICIENT STORAGE ALLOCATED FOR IP(*), NEED LIP = '//xern1,2,1)
      RETURN
    END IF
  END IF
  !
  !     Compute number of possible right multiplying Householder
  !     transformations.
  !
  m = Me + Ma + Mg
  IF ( N<=0.OR.m<=0 ) THEN
    Mode = 0
    Rnorme = 0
    Rnorml = 0
    RETURN
  END IF
  !
  IF ( Mdw<m ) THEN
    CALL XERMSG('LSEI','MDW.LT.ME+MA+MG IS AN ERROR',2,1)
    RETURN
  END IF
  !
  np1 = N + 1
  kranke = MIN(Me,N)
  n1 = 2*kranke + 1
  n2 = n1 + N
  !
  !     Set nominal values.
  !
  !     The nominal column scaling used in the code is
  !     the identity scaling.
  !
  Ws(n1:n1+N-1) = 1.E0
  !
  !     No covariance matrix is nominally computed.
  !
  cov = .FALSE.
  !
  !     Process option vector.
  !     Define bound for number of options to change.
  !
  nopt = 1000
  ntimes = 0
  !
  !     Define bound for positive values of LINK.
  !
  nlink = 100000
  last = 1
  link = INT( Prgopt(1) )
  IF ( link==0.OR.link>nlink ) THEN
    CALL XERMSG('LSEI','THE OPTION VECTOR IS UNDEFINED',2,1)
    RETURN
  END IF
  DO
    !
    IF ( link>1 ) THEN
      ntimes = ntimes + 1
      IF ( ntimes>nopt ) THEN
        CALL XERMSG('LSEI',&
          'THE LINKS IN THE OPTION VECTOR ARE CYCLING.',2,1)
        RETURN
      END IF
      !
      key = INT( Prgopt(last+1) )
      IF ( key==1 ) THEN
        cov = Prgopt(last+2)/=0.E0
      ELSEIF ( key==2.AND.Prgopt(last+2)/=0.E0 ) THEN
        DO j = 1, N
          t = NORM2(W(1:m,j))
          IF ( t/=0.E0 ) t = 1.E0/t
          Ws(j+n1-1) = t
        END DO
      ELSEIF ( key==3 ) THEN
        Ws(n1:n1+N-1) = Prgopt(last+2:last+N+1)
      ELSEIF ( key==4 ) THEN
        tau = MAX(srelpr,Prgopt(last+2))
      END IF
      !
      next = INT( Prgopt(link) )
      IF ( next<=0.OR.next>nlink ) THEN
        CALL XERMSG('LSEI','THE OPTION VECTOR IS UNDEFINED',2,1)
        RETURN
      END IF
      !
      last = link
      link = next
      CYCLE
    END IF
    !
    DO j = 1, N
      W(1:m,j) = W(1:m,j)*Ws(n1+j-1)
    END DO
    !
    IF ( cov.AND.Mdw<N ) THEN
      CALL XERMSG('LSEI',&
        'MDW .LT. N WHEN COV MATRIX NEEDED, IS AN ERROR',2,1)
      RETURN
    END IF
    !
    !     Problem definition and option vector OK.
    !
    Mode = 0
    !
    !     Compute norm of equality constraint matrix and right side.
    !
    enorm = 0.E0
    DO j = 1, N
      enorm = MAX(enorm,SUM(ABS(W(1:Me,j))))
    END DO
    !
    fnorm = SUM(ABS(W(1:Me,np1)))
    snmax = 0.E0
    rnmax = 0.E0
    DO i = 1, kranke
      !
      !        Compute maximum ratio of vector lengths. Partition is at
      !        column I.
      !
      DO k = i, Me
        sn = NORM2(W(k,i:N))**2
        rn = NORM2(W(k,1:i-1))**2
        IF ( rn==0.E0.AND.sn>snmax ) THEN
          snmax = sn
          imax = k
        ELSEIF ( k==i.OR.sn*rnmax>rn*snmax ) THEN
          snmax = sn
          rnmax = rn
          imax = k
        END IF
      END DO
      !
      !        Interchange rows if necessary.
      !
      IF ( i/=imax ) CALL SSWAP(np1,W(i,1),Mdw,W(imax,1),Mdw)
      IF ( snmax>rnmax*tau**2 ) THEN
        !
        !        Eliminate elements I+1,...,N in row I.
        !
        CALL H12(1,i,i+1,N,W(i,1),Mdw,Ws(i),W(i+1,1),Mdw,1,m-i)
      ELSE
        kranke = i - 1
        EXIT
      END IF
    END DO
    EXIT
  END DO
  !
  !     Save diagonal terms of lower trapezoidal matrix.
  !
  DO k = 1, kranke
    Ws(kranke+k) = W(k,k)
  END DO
  !
  !     Use Householder transformation from left to achieve
  !     KRANKE by KRANKE upper triangular form.
  !
  IF ( kranke<Me ) THEN
    DO k = kranke, 1, -1
      !
      !           Apply transformation to matrix cols. 1,...,K-1.
      !
      CALL H12(1,k,kranke+1,Me,W(1,k),1,up,W,1,Mdw,k-1)
      !
      !           Apply to rt side vector.
      !
      CALL H12(2,k,kranke+1,Me,W(1,k),1,up,W(1,np1),1,1,1)
    END DO
  END IF
  !
  !     Solve for variables 1,...,KRANKE in new coordinates.
  !
  X(1:kranke) = W(1:kranke,np1)
  DO i = 1, kranke
    X(i) = (X(i)-DOT_PRODUCT(W(i,1:i-1),X(1:i-1)))/W(i,i)
  END DO
  !
  !     Compute residuals for reduced problem.
  !
  mep1 = Me + 1
  Rnorml = 0.E0
  DO i = mep1, m
    W(i,np1) = W(i,np1) - DOT_PRODUCT(W(i,1:kranke),X(1:kranke))
    sn = NORM2(W(i,1:kranke))**2
    rn = NORM2(W(i,kranke+1:N))**2
    IF ( rn<=sn*tau**2.AND.kranke<N ) W(i,kranke+1:N) = 0.E0
  END DO
  !
  !     Compute equality constraint equations residual length.
  !
  Rnorme = NORM2(W(kranke+1:Me,np1))
  !
  !     Move reduced problem data upward if KRANKE.LT.ME.
  !
  IF ( kranke<Me ) THEN
    DO j = 1, np1
      W(kranke+1:kranke+m-Me,j) = W(Me+1:m,j)
    END DO
  END IF
  !
  !     Compute solution of reduced problem.
  !
  CALL LSI(W(kranke+1,kranke+1),Mdw,Ma,Mg,N-kranke,Prgopt,X(kranke+1),&
    Rnorml,Mode,Ws(n2),Ip(2))
  !
  !     Test for consistency of equality constraints.
  !
  IF ( Me>0 ) THEN
    mdeqc = 0
    xnrme = SUM(ABS(W(1:kranke,np1)))
    IF ( Rnorme>tau*(enorm*xnrme+fnorm) ) mdeqc = 1
    Mode = Mode + mdeqc
    !
    !        Check if solution to equality constraints satisfies inequality
    !        constraints when there are no degrees of freedom left.
    !
    IF ( kranke==N.AND.Mg>0 ) THEN
      xnorm = SUM(ABS(X(1:N)))
      mapke1 = Ma + kranke + 1
      mend = Ma + kranke + Mg
      DO i = mapke1, mend
        sizee = SUM(ABS(W(i,1:N)))*xnorm + ABS(W(i,np1))
        IF ( W(i,np1)>tau*sizee ) THEN
          Mode = Mode + 2
          GOTO 100
        END IF
      END DO
    END IF
  END IF
  !
  !     Replace diagonal terms of lower trapezoidal matrix.
  !
  IF ( kranke>0 ) THEN
    DO i = 1, kranke
      W(i,i) = Ws(kranke+i)
    END DO
    !
    !        Reapply transformation to put solution in original coordinates.
    !
    DO i = kranke, 1, -1
      CALL H12(2,i,i+1,N,W(i,1),Mdw,Ws(i),X,1,1,1)
    END DO
    !
    !        Compute covariance matrix of equality constrained problem.
    !
    IF ( cov ) THEN
      DO j = MIN(kranke,N-1), 1, -1
        rb = Ws(j)*W(j,j)
        IF ( rb/=0.E0 ) rb = 1.E0/rb
        jp1 = j + 1
        DO i = jp1, N
          W(i,j) = rb*DOT_PRODUCT(W(i,j+1:N),W(j,j+1:N))
        END DO
        !
        gam = 0.5E0*rb*DOT_PRODUCT(W(j+1:N,j),W(j,j+1:N))
        CALL SAXPY(N-j,gam,W(j,jp1),Mdw,W(jp1,j),1)
        DO i = jp1, N
          DO k = i, N
            W(i,k) = W(i,k) + W(j,i)*W(k,j) + W(i,j)*W(j,k)
            W(k,i) = W(i,k)
          END DO
        END DO
        uj = Ws(j)
        vj = gam*uj
        W(j,j) = uj*vj + uj*vj
        DO i = jp1, N
          W(j,i) = uj*W(i,j) + vj*W(j,i)
        END DO
        W(j+1:N,j) = W(j,j+1:N)
      END DO
    END IF
  END IF
  !
  !     Apply the scaling to the covariance matrix.
  !
  IF ( cov ) THEN
    DO i = 1, N
      W(i,1:N) = W(i,1:N)*Ws(i+n1-1)
      W(1:N,i) = W(1:N,i)*Ws(i+n1-1)
    END DO
  END IF
  !
  !     Rescale solution vector.
  !
  100 CONTINUE
  IF ( Mode<=1 ) THEN
    DO j = 1, N
      X(j) = X(j)*Ws(n1+j-1)
    END DO
  END IF
  !
  Ip(1) = kranke
  Ip(3) = Ip(3) + 2*kranke + N
END SUBROUTINE LSEI
