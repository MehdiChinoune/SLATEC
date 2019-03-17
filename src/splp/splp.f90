!DECK SPLP
SUBROUTINE SPLP(USRMAT,Mrelas,Nvars,Costs,Prgopt,Dattrv,Bl,Bu,Ind,Info,&
    Primal,Duals,Ibasis,Work,Lw,Iwork,Liw)
  IMPLICIT NONE
  !***BEGIN PROLOGUE  SPLP
  !***PURPOSE  Solve linear programming problems involving at
  !            most a few thousand constraints and variables.
  !            Takes advantage of sparsity in the constraint matrix.
  !***LIBRARY   SLATEC
  !***CATEGORY  G2A2
  !***TYPE      SINGLE PRECISION (SPLP-S, DSPLP-D)
  !***KEYWORDS  LINEAR CONSTRAINTS, LINEAR OPTIMIZATION,
  !             LINEAR PROGRAMMING, LP, SPARSE CONSTRAINTS
  !***AUTHOR  Hanson, R. J., (SNLA)
  !           Hiebert, K. L., (SNLA)
  !***DESCRIPTION
  !
  !     These are the short usage instructions; for details about
  !     other features, options and methods for defining the matrix
  !     A, see the extended usage instructions which are contained in
  !     the Long Description section below.
  !
  !   |------------|
  !   |Introduction|
  !   |------------|
  !     The subprogram SPLP( ) solves a linear optimization problem.
  !     The problem statement is as follows
  !
  !                         minimize (transpose of costs)*x
  !                         subject to A*x=w.
  !
  !     The entries of the unknowns x and w may have simple lower or
  !     upper bounds (or both), or be free to take on any value.  By
  !     setting the bounds for x and w, the user is imposing the con-
  !     straints of the problem.  The matrix A has MRELAS rows and
  !     NVARS columns.  The vectors costs, x, and w respectively
  !     have NVARS, NVARS, and MRELAS number of entries.
  !
  !     The input for the problem includes the problem dimensions,
  !     MRELAS and NVARS, the array COSTS(*), data for the matrix
  !     A, and the bound information for the unknowns x and w, BL(*),
  !     BU(*), and IND(*).  Only the nonzero entries of the matrix A
  !     are passed to SPLP( ).
  !
  !     The output from the problem (when output flag INFO=1) includes
  !     optimal values for x and w in PRIMAL(*), optimal values for
  !     dual variables of the equations A*x=w and the simple bounds
  !     on x in  DUALS(*), and the indices of the basic columns,
  !     IBASIS(*).
  !
  !  |------------------------------|
  !  |Fortran Declarations Required:|
  !  |------------------------------|
  !
  !     DIMENSION COSTS(NVARS),PRGOPT(*),DATTRV(*),
  !    *BL(NVARS+MRELAS),BU(NVARS+MRELAS),IND(NVARS+MRELAS),
  !    *PRIMAL(NVARS+MRELAS),DUALS(MRELAS+NVARS),IBASIS(NVARS+MRELAS),
  !    *WORK(LW),IWORK(LIW)
  !
  !     EXTERNAL USRMAT
  !
  !     The dimensions of PRGOPT(*) and DATTRV(*) must be at least 1.
  !     The exact lengths will be determined by user-required options and
  !     data transferred to the subprogram USRMAT( ).
  !
  !     The values of LW and LIW, the lengths of the arrays WORK(*)
  !     and IWORK(*), must satisfy the inequalities
  !
  !               LW .GE. 4*NVARS+ 8*MRELAS+LAMAT+  LBM
  !               LIW.GE.   NVARS+11*MRELAS+LAMAT+2*LBM
  !
  !     It is an error if they do not both satisfy these inequalities.
  !     (The subprogram will inform the user of the required lengths
  !     if either LW or LIW is wrong.)  The values of LAMAT and LBM
  !     nominally are
  !
  !               LAMAT=4*NVARS+7
  !     and       LBM  =8*MRELAS
  !
  !     LAMAT determines the length of the sparse matrix storage area.
  !     The value of LBM determines the amount of storage available
  !     to decompose and update the active basis matrix.
  !
  !  |------|
  !  |Input:|
  !  |------|
  !
  !     MRELAS,NVARS
  !     ------------
  !     These parameters are respectively the number of constraints (the
  !     linear relations A*x=w that the unknowns x and w are to satisfy)
  !     and the number of entries in the vector x.  Both must be .GE. 1.
  !     Other values are errors.
  !
  !     COSTS(*)
  !     --------
  !     The NVARS entries of this array are the coefficients of the
  !     linear objective function.  The value COSTS(J) is the
  !     multiplier for variable J of the unknown vector x.  Each
  !     entry of this array must be defined.
  !
  !     USRMAT
  !     ------
  !     This is the name of a specific subprogram in the SPLP( ) package
  !     used to define the matrix A.  In this usage mode of SPLP( )
  !     the user places the nonzero entries of A in the
  !     array DATTRV(*) as given in the description of that parameter.
  !     The name USRMAT must appear in a Fortran EXTERNAL statement.
  !
  !     DATTRV(*)
  !     ---------
  !     The array DATTRV(*) contains data for the matrix A as follows:
  !     Each column (numbered J) requires (floating point) data con-
  !     sisting of the value (-J) followed by pairs of values.  Each pair
  !     consists of the row index immediately followed by the value
  !     of the matrix at that entry.  A value of J=0 signals that there
  !     are no more columns.  The required length of
  !     DATTRV(*) is 2*no. of nonzeros + NVARS + 1.
  !
  !     BL(*),BU(*),IND(*)
  !     ------------------
  !     The values of IND(*) are input parameters that define
  !     the form of the bounds for the unknowns x and w.  The values for
  !     the bounds are found in the arrays BL(*) and BU(*) as follows.
  !
  !     For values of J between 1 and NVARS,
  !          if IND(J)=1, then X(J) .GE. BL(J); BU(J) is not used.
  !          if IND(J)=2, then X(J) .LE. BU(J); BL(J) is not used.
  !          if IND(J)=3, then BL(J) .LE. X(J) .LE. BU(J),(BL(J)=BU(J) ok)
  !          if IND(J)=4, then X(J) is free to have any value,
  !          and BL(J), BU(J) are not used.
  !
  !     For values of I between NVARS+1 and NVARS+MRELAS,
  !          if IND(I)=1, then W(I-NVARS) .GE. BL(I); BU(I) is not used.
  !          if IND(I)=2, then W(I-NVARS) .LE. BU(I); BL(I) is not used.
  !          if IND(I)=3, then BL(I) .LE. W(I-NVARS) .LE. BU(I),
  !          (BL(I)=BU(I) is ok).
  !          if IND(I)=4, then W(I-NVARS) is free to have any value,
  !          and BL(I), BU(I) are not used.
  !
  !     A value of IND(*) not equal to 1,2,3 or 4 is an error.  When
  !     IND(I)=3, BL(I) must be .LE. BU(I).  The condition BL(I).GT.
  !     BU(I) indicates infeasibility and is an error.
  !
  !     PRGOPT(*)
  !     ---------
  !     This array is used to redefine various parameters within SPLP( ).
  !     Frequently, perhaps most of the time, a user will be satisfied
  !     and obtain the solutions with no changes to any of these
  !     parameters.  To try this, simply set PRGOPT(1)=1.E0.
  !
  !     For users with more sophisticated needs, SPLP( ) provides several
  !     options that may be used to take advantage of more detailed
  !     knowledge of the problem or satisfy other utilitarian needs.
  !     The complete description of how to use this option array to
  !     utilize additional subprogram features is found under the
  !     heading  of SPLP( ) Subprogram Options in the Extended
  !     Usage Instructions.
  !
  !     Briefly, the user should note the following value of the parameter
  !     KEY and the corresponding task or feature desired before turning
  !     to that document.
  !
  !     Value     Brief Statement of Purpose for Option
  !     of KEY
  !     ------    -------------------------------------
  !     50        Change from a minimization problem to a
  !               maximization problem.
  !     51        Change the amount of printed output.
  !               Normally, no printed output is obtained.
  !     52        Redefine the line length and precision used
  !               for the printed output.
  !     53        Redefine the values of LAMAT and LBM that
  !               were discussed above under the heading
  !               Fortran Declarations Required.
  !     54        Redefine the unit number where pages of the sparse
  !               data matrix A are stored.  Normally, the unit
  !               number is 1.
  !     55        A computation, partially completed, is
  !               being continued.  Read the up-to-date
  !               partial results from unit number 2.
  !     56        Redefine the unit number where the partial results
  !               are stored.  Normally, the unit number is 2.
  !     57        Save partial results on unit 2 either after
  !               maximum iterations or at the optimum.
  !     58        Redefine the value for the maximum number of
  !               iterations.  Normally, the maximum number of
  !               iterations is 3*(NVARS+MRELAS).
  !     59        Provide SPLP( ) with a starting (feasible)
  !               nonsingular basis.  Normally, SPLP( ) starts
  !               with the identity matrix columns corresponding
  !               to the vector w.
  !     60        The user has provided scale factors for the
  !               columns of A.  Normally, SPLP( ) computes scale
  !               factors that are the reciprocals of the max. norm
  !               of each column.
  !     61        The user has provided a scale factor
  !               for the vector costs.  Normally, SPLP( ) computes
  !               a scale factor equal to the reciprocal of the
  !               max. norm of the vector costs after the column
  !               scaling for the data matrix has been applied.
  !     62        Size parameters, namely the smallest and
  !               largest magnitudes of nonzero entries in
  !               the matrix A, are provided.  Values noted
  !               outside this range are to be considered errors.
  !     63        Redefine the tolerance required in
  !               evaluating residuals for feasibility.
  !               Normally, this value is set to RELPR,
  !               where RELPR = relative precision of the arithmetic.
  !     64        Change the criterion for bringing new variables
  !               into the basis from the steepest edge (best
  !               local move) to the minimum reduced cost.
  !     65        Redefine the value for the number of iterations
  !               between recalculating the error in the primal
  !               solution.  Normally, this value is equal to ten.
  !     66        Perform "partial pricing" on variable selection.
  !               Redefine the value for the number of negative
  !               reduced costs to compute (at most) when finding
  !               a variable to enter the basis.  Normally this
  !               value is set to NVARS.  This implies that no
  !               "partial pricing" is used.
  !     67        Adjust the tuning factor (normally one) to apply
  !               to the primal and dual error estimates.
  !     68        Pass  information to the  subprogram  FULMAT(),
  !               provided with the SPLP() package, so that a Fortran
  !               two-dimensional array can be used as the argument
  !               DATTRV(*).
  !     69        Pass an absolute tolerance to use for the feasibility
  !               test when the usual relative error test indicates
  !               infeasibility.  The nominal value of this tolerance,
  !               TOLABS, is zero.
  !
  !
  !  |---------------|
  !  |Working Arrays:|
  !  |---------------|
  !
  !     WORK(*),LW,
  !     IWORK(*),LIW
  !     ------------
  !     The arrays WORK(*) and IWORK(*) are respectively floating point
  !     and type INTEGER working arrays for SPLP( ) and its
  !     subprograms.  The lengths of these arrays are respectively
  !     LW and LIW.  These parameters must satisfy the inequalities
  !     noted above under the heading "Fortran Declarations Required:"
  !     It is an error if either value is too small.
  !
  !  |----------------------------|
  !  |Input/Output files required:|
  !  |----------------------------|
  !
  !     Fortran unit 1 is used by SPLP( ) to store the sparse matrix A
  !     out of high-speed memory.  A crude
  !     upper bound for the amount of information written on unit 1
  !     is 6*nz, where nz is the number of nonzero entries in A.
  !
  !  |-------|
  !  |Output:|
  !  |-------|
  !
  !     INFO,PRIMAL(*),DUALS(*)
  !     -----------------------
  !     The integer flag INFO indicates why SPLP( ) has returned to the
  !     user.  If INFO=1 the solution has been computed.  In this case
  !     X(J)=PRIMAL(J) and W(I)=PRIMAL(I+NVARS).  The dual variables
  !     for the equations A*x=w are in the array DUALS(I)=dual for
  !     equation number I.  The dual value for the component X(J) that
  !     has an upper or lower bound (or both) is returned in
  !     DUALS(J+MRELAS).  The only other values for INFO are .LT. 0.
  !     The meaning of these values can be found by reading
  !     the diagnostic message in the output file, or by looking for
  !     error number = (-INFO) in the Extended Usage Instructions
  !     under the heading:
  !
  !          List of SPLP( ) Error and Diagnostic Messages.
  !
  !     BL(*),BU(*),IND(*)
  !     ------------------
  !     These arrays are output parameters only under the (unusual)
  !     circumstances where the stated problem is infeasible, has an
  !     unbounded optimum value, or both.  These respective conditions
  !     correspond to INFO=-1,-2 or -3.    See the Extended
  !     Usage Instructions for further details.
  !
  !     IBASIS(I),I=1,...,MRELAS
  !     ------------------------
  !     This array contains the indices of the variables that are
  !     in the active basis set at the solution (INFO=1).  A value
  !     of IBASIS(I) between 1 and NVARS corresponds to the variable
  !     X(IBASIS(I)).  A value of IBASIS(I) between NVARS+1 and NVARS+
  !     MRELAS corresponds to the variable W(IBASIS(I)-NVARS).
  !
  ! *Long Description:
  !
  !     SUBROUTINE SPLP(USRMAT,MRELAS,NVARS,COSTS,PRGOPT,DATTRV,
  !    *           BL,BU,IND,INFO,PRIMAL,DUALS,IBASIS,WORK,LW,IWORK,LIW)
  !
  !     |------------|
  !     |Introduction|
  !     |------------|
  !     The subprogram SPLP( ) solves a linear optimization problem.
  !     The problem statement is as follows
  !
  !                         minimize (transpose of costs)*x
  !                         subject to A*x=w.
  !
  !     The entries of the unknowns x and w may have simple lower or
  !     upper bounds (or both), or be free to take on any value.  By
  !     setting the bounds for x and w, the user is imposing the con-
  !     straints of the problem.
  !
  !     (The problem may also be stated as a maximization
  !     problem.  This is done by means of input in the option array
  !     PRGOPT(*).)  The matrix A has MRELAS rows and NVARS columns.  The
  !     vectors costs, x, and w respectively have NVARS, NVARS, and
  !     MRELAS number of entries.
  !
  !     The input for the problem includes the problem dimensions,
  !     MRELAS and NVARS, the array COSTS(*), data for the matrix
  !     A, and the bound information for the unknowns x and w, BL(*),
  !     BU(*), and IND(*).
  !
  !     The output from the problem (when output flag INFO=1) includes
  !     optimal values for x and w in PRIMAL(*), optimal values for
  !     dual variables of the equations A*x=w and the simple bounds
  !     on x in  DUALS(*), and the indices of the basic columns in
  !     IBASIS(*).
  !
  !  |------------------------------|
  !  |Fortran Declarations Required:|
  !  |------------------------------|
  !
  !     DIMENSION COSTS(NVARS),PRGOPT(*),DATTRV(*),
  !    *BL(NVARS+MRELAS),BU(NVARS+MRELAS),IND(NVARS+MRELAS),
  !    *PRIMAL(NVARS+MRELAS),DUALS(MRELAS+NVARS),IBASIS(NVARS+MRELAS),
  !    *WORK(LW),IWORK(LIW)
  !
  !     EXTERNAL USRMAT (or 'NAME', if user provides the subprogram)
  !
  !     The dimensions of PRGOPT(*) and DATTRV(*) must be at least 1.
  !     The exact lengths will be determined by user-required options and
  !     data transferred to the subprogram USRMAT( ) ( or 'NAME').
  !
  !     The values of LW and LIW, the lengths of the arrays WORK(*)
  !     and IWORK(*), must satisfy the inequalities
  !
  !               LW .GE. 4*NVARS+ 8*MRELAS+LAMAT+  LBM
  !               LIW.GE.   NVARS+11*MRELAS+LAMAT+2*LBM
  !
  !     It is an error if they do not both satisfy these inequalities.
  !     (The subprogram will inform the user of the required lengths
  !     if either LW or LIW is wrong.)  The values of LAMAT and LBM
  !     nominally are
  !
  !               LAMAT=4*NVARS+7
  !     and       LBM  =8*MRELAS
  !
  !     These values will be as shown unless the user changes them by
  !     means of input in the option array PRGOPT(*).  The value of LAMAT
  !     determines the length of the sparse matrix "staging" area.
  !     For reasons of efficiency the user may want to increase the value
  !     of LAMAT.  The value of LBM determines the amount of storage
  !     available to decompose and update the active basis matrix.
  !     Due to exhausting the working space because of fill-in,
  !     it may be necessary for the user to increase the value of LBM.
  !     (If this situation occurs an informative diagnostic is printed
  !     and a value of INFO=-28 is obtained as an output parameter.)
  !
  !  |------|
  !  |Input:|
  !  |------|
  !
  !     MRELAS,NVARS
  !     ------------
  !     These parameters are respectively the number of constraints (the
  !     linear relations A*x=w that the unknowns x and w are to satisfy)
  !     and the number of entries in the vector x.  Both must be .GE. 1.
  !     Other values are errors.
  !
  !     COSTS(*)
  !     --------
  !     The NVARS entries of this array are the coefficients of the
  !     linear objective function.  The value COSTS(J) is the
  !     multiplier for variable J of the unknown vector x.  Each
  !     entry of this array must be defined.  This array can be changed
  !     by the user between restarts.  See options with KEY=55,57 for
  !     details of checkpointing and restarting.
  !
  !     USRMAT
  !     ------
  !     This is the name of a specific subprogram in the SPLP( ) package
  !     that is used to define the matrix entries when this data is passed
  !     to SPLP( ) as a linear array.  In this usage mode of SPLP( )
  !     the user gives information about the nonzero entries of A
  !     in DATTRV(*) as given under the description of that parameter.
  !     The name USRMAT must appear in a Fortran EXTERNAL statement.
  !     Users who are passing the matrix data with USRMAT( ) can skip
  !     directly to the description of the input parameter DATTRV(*).
  !     Also see option 68 for passing the constraint matrix data using
  !     a standard Fortran two-dimensional array.
  !
  !     If the user chooses to provide a subprogram 'NAME'( ) to
  !     define the matrix A, then DATTRV(*) may be used to pass floating
  !     point data from the user's program unit to the subprogram
  !     'NAME'( ). The content of DATTRV(*) is not changed in any way.
  !
  !     The subprogram 'NAME'( ) can be of the user's choice
  !     but it must meet Fortran standards and it must appear in a
  !     Fortran EXTERNAL statement.  The first statement of the subprogram
  !     has the form
  !
  !          SUBROUTINE 'NAME'(I,J,AIJ, INDCAT, PRGOPT, DATTRV, IFLAG)
  !
  !     The variables I,J, INDCAT, IFLAG(10) are type INTEGER,
  !          while  AIJ, PRGOPT(*),DATTRV(*) are type REAL.
  !
  !     The user interacts with the contents of IFLAG(*) to
  !     direct the appropriate action.  The algorithmic steps are
  !     as follows.
  !
  !          Test IFLAG(1).
  !
  !             IF(IFLAG(1).EQ.1) THEN
  !
  !               Initialize the necessary pointers and data
  !               for defining the matrix A.  The contents
  !               of IFLAG(K), K=2,...,10, may be used for
  !               storage of the pointers.  This array remains intact
  !               between calls to 'NAME'( ) by SPLP( ).
  !               RETURN
  !
  !             END IF
  !
  !             IF(IFLAG(1).EQ.2) THEN
  !
  !               Define one set of values for I,J,AIJ, and INDCAT.
  !               Each nonzero entry of A must be defined this way.
  !               These values can be defined in any convenient order.
  !               (It is most efficient to define the data by
  !               columns in the order 1,...,NVARS; within each
  !               column define the entries in the order 1,...,MRELAS.)
  !               If this is the last matrix value to be
  !               defined or updated, then set IFLAG(1)=3.
  !               (When I and J are positive and respectively no larger
  !               than MRELAS and NVARS, the value of AIJ is used to
  !               define (or update) row I and column J of A.)
  !               RETURN
  !
  !             END IF
  !
  !               END
  !
  !     Remarks:  The values of I and J are the row and column
  !               indices for the nonzero entries of the matrix A.
  !               The value of this entry is AIJ.
  !               Set INDCAT=0 if this value defines that entry.
  !               Set INDCAT=1 if this entry is to be updated,
  !                            new entry=old entry+AIJ.
  !               A value of I not between 1 and MRELAS, a value of J
  !               not between 1 and NVARS, or a value of INDCAT
  !               not equal to 0 or 1 are each errors.
  !
  !               The contents of IFLAG(K), K=2,...,10, can be used to
  !               remember the status (of the process of defining the
  !               matrix entries) between calls to 'NAME'( ) by SPLP( ).
  !               On entry to 'NAME'( ), only the values 1 or 2 will be
  !               in IFLAG(1).  More than 2*NVARS*MRELAS definitions of
  !               the matrix elements is considered an error because
  !               it suggests an infinite loop in the user-written
  !               subprogram 'NAME'( ).  Any matrix element not
  !               provided by 'NAME'( ) is defined to be zero.
  !
  !               The REAL arrays PRGOPT(*) and DATTRV(*) are passed as
  !               arguments directly from SPLP( ) to 'NAME'( ).
  !               The array PRGOPT(*) contains any user-defined program
  !               options.  In this usage mode the array DATTRV(*) may
  !               now contain any (type REAL) data that the user needs
  !               to define the matrix A.  Both arrays PRGOPT(*) and
  !               DATTRV(*) remain intact between calls to 'NAME'( )
  !               by SPLP( ).
  !     Here is a subprogram that communicates the matrix values for A,
  !     as represented in DATTRV(*), to SPLP( ).  This subprogram,
  !     called USRMAT( ), is included as part of the SPLP( ) package.
  !     This subprogram 'decodes' the array DATTRV(*) and defines the
  !     nonzero entries of the matrix  A for SPLP( ) to store.  This
  !     listing is presented here as a guide and example
  !     for the users who find it necessary to write their own subroutine
  !     for this purpose.  The contents of DATTRV(*) are given below in
  !     the description of that parameter.
  !
  !     SUBROUTINE USRMAT(I,J,AIJ, INDCAT,PRGOPT,DATTRV,IFLAG)
  !     DIMENSION PRGOPT(*),DATTRV(*),IFLAG(10)
  !
  !     IF(IFLAG(1).EQ.1) THEN
  !
  !     THIS IS THE INITIALIZATION STEP.  THE VALUES OF IFLAG(K),K=2,3,4,
  !     ARE RESPECTIVELY THE COLUMN INDEX, THE ROW INDEX (OR THE NEXT COL.
  !     INDEX), AND THE POINTER TO THE MATRIX ENTRY'S VALUE WITHIN
  !     DATTRV(*).  ALSO CHECK (DATTRV(1)=0.) SIGNIFYING NO DATA.
  !          IF(DATTRV(1).EQ.0.) THEN
  !          I = 0
  !          J = 0
  !          IFLAG(1) = 3
  !          ELSE
  !          IFLAG(2)=-DATTRV(1)
  !          IFLAG(3)= DATTRV(2)
  !          IFLAG(4)= 3
  !          END IF
  !
  !          RETURN
  !     ELSE
  !          J=IFLAG(2)
  !          I=IFLAG(3)
  !          L=IFLAG(4)
  !          IF(I.EQ.0) THEN
  !
  !     SIGNAL THAT ALL OF THE NONZERO ENTRIES HAVE BEEN DEFINED.
  !               IFLAG(1)=3
  !               RETURN
  !          ELSE IF(I.LT.0) THEN
  !
  !     SIGNAL THAT A SWITCH IS MADE TO A NEW COLUMN.
  !               J=-I
  !               I=DATTRV(L)
  !               L=L+1
  !          END IF
  !
  !          AIJ=DATTRV(L)
  !
  !     UPDATE THE INDICES AND POINTERS FOR THE NEXT ENTRY.
  !          IFLAG(2)=J
  !          IFLAG(3)=DATTRV(L+1)
  !          IFLAG(4)=L+2
  !
  !     INDCAT=0 DENOTES THAT ENTRIES OF THE MATRIX ARE ASSIGNED THE
  !     VALUES FROM DATTRV(*).  NO ACCUMULATION IS PERFORMED.
  !          INDCAT=0
  !          RETURN
  !     END IF
  !     END
  !
  !     DATTRV(*)
  !     ---------
  !     If the user chooses to use the provided subprogram USRMAT( ) then
  !     the array DATTRV(*) contains data for the matrix A as follows:
  !     Each column (numbered J) requires (floating point) data con-
  !     sisting of the value (-J) followed by pairs of values.  Each pair
  !     consists of the row index immediately followed by the value
  !     of the matrix at that entry.  A value of J=0 signals that there
  !     are no more columns.  (See "Example of SPLP( ) Usage," below.)
  !     The dimension of DATTRV(*) must be 2*no. of nonzeros
  !     + NVARS + 1 in this usage.  No checking of the array
  !     length is done by the subprogram package.
  !
  !     If the Save/Restore feature is in use (see options with
  !     KEY=55,57 for details of checkpointing and restarting)
  !     USRMAT( ) can be used to redefine entries of the matrix.
  !     The matrix entries are redefined or overwritten.  No accum-
  !     ulation is performed.
  !     Any other nonzero entry of A, defined in a previous call to
  !     SPLP( ), remain intact.
  !
  !     BL(*),BU(*),IND(*)
  !     ------------------
  !     The values of IND(*) are input parameters that define
  !     the form of the bounds for the unknowns x and w.  The values for
  !     the bounds are found in the arrays BL(*) and BU(*) as follows.
  !
  !     For values of J between 1 and NVARS,
  !          if IND(J)=1, then X(J) .GE. BL(J); BU(J) is not used.
  !          if IND(J)=2, then X(J) .LE. BU(J); BL(J) is not used.
  !          if IND(J)=3, then BL(J) .LE. X(J) .LE. BU(J),(BL(J)=BU(J) ok)
  !          if IND(J)=4, then X(J) is free to have any value,
  !          and BL(J), BU(J) are not used.
  !
  !     For values of I between NVARS+1 and NVARS+MRELAS,
  !          if IND(I)=1, then W(I-NVARS) .GE. BL(I); BU(I) is not used.
  !          if IND(I)=2, then W(I-NVARS) .LE. BU(I); BL(I) is not used.
  !          if IND(I)=3, then BL(I) .LE. W(I-NVARS) .LE. BU(I),
  !          (BL(I)=BU(I) is ok).
  !          if IND(I)=4, then W(I-NVARS) is free to have any value,
  !          and BL(I), BU(I) are not used.
  !
  !     A value of IND(*) not equal to 1,2,3 or 4 is an error.  When
  !     IND(I)=3, BL(I) must be .LE. BU(I).  The condition BL(I).GT.
  !     BU(I) indicates infeasibility and is an error.  These
  !     arrays can be changed by the user between restarts.  See
  !     options with KEY=55,57 for details of checkpointing and
  !     restarting.
  !
  !     PRGOPT(*)
  !     ---------
  !     This array is used to redefine various parameters within SPLP( ).
  !     Frequently, perhaps most of the time, a user will be satisfied
  !     and obtain the solutions with no changes to any of these
  !     parameters.  To try this, simply set PRGOPT(1)=1.E0.
  !
  !     For users with more sophisticated needs, SPLP( ) provides several
  !     options that may be used to take advantage of more detailed
  !     knowledge of the problem or satisfy other utilitarian needs.
  !     The complete description of how to use this option array to
  !     utilize additional subprogram features is found under the
  !     heading "Usage of SPLP( ) Subprogram Options."
  !
  !     Briefly, the user should note the following value of the parameter
  !     KEY and the corresponding task or feature desired before turning
  !     to that section.
  !
  !     Value     Brief Statement of Purpose for Option
  !     of KEY
  !     ------    -------------------------------------
  !     50        Change from a minimization problem to a
  !               maximization problem.
  !     51        Change the amount of printed output.
  !               Normally, no printed output is obtained.
  !     52        Redefine the line length and precision used
  !               for the printed output.
  !     53        Redefine the values of LAMAT and LBM that
  !               were discussed above under the heading
  !               Fortran Declarations Required.
  !     54        Redefine the unit number where pages of the sparse
  !               data matrix A are stored.  Normally, the unit
  !               number is 1.
  !     55        A computation, partially completed, is
  !               being continued.  Read the up-to-date
  !               partial results from unit number 2.
  !     56        Redefine the unit number where the partial results
  !               are stored.  Normally, the unit number is 2.
  !     57        Save partial results on unit 2 either after
  !               maximum iterations or at the optimum.
  !     58        Redefine the value for the maximum number of
  !               iterations.  Normally, the maximum number of
  !               iterations is 3*(NVARS+MRELAS).
  !     59        Provide SPLP( ) with a starting (feasible)
  !               nonsingular basis.  Normally, SPLP( ) starts
  !               with the identity matrix columns corresponding
  !               to the vector w.
  !     60        The user has provided scale factors for the
  !               columns of A.  Normally, SPLP( ) computes scale
  !               factors that are the reciprocals of the max. norm
  !               of each column.
  !     61        The user has provided a scale factor
  !               for the vector costs.  Normally, SPLP( ) computes
  !               a scale factor equal to the reciprocal of the
  !               max. norm of the vector costs after the column
  !               scaling for the data matrix has been applied.
  !     62        Size parameters, namely the smallest and
  !               largest magnitudes of nonzero entries in
  !               the matrix A, are provided.  Values noted
  !               outside this range are to be considered errors.
  !     63        Redefine the tolerance required in
  !               evaluating residuals for feasibility.
  !               Normally, this value is set to the value RELPR,
  !               where RELPR = relative precision of the arithmetic.
  !     64        Change the criterion for bringing new variables
  !               into the basis from the steepest edge (best
  !               local move) to the minimum reduced cost.
  !     65        Redefine the value for the number of iterations
  !               between recalculating the error in the primal
  !               solution.  Normally, this value is equal to ten.
  !     66        Perform "partial pricing" on variable selection.
  !               Redefine the value for the number of negative
  !               reduced costs to compute (at most) when finding
  !               a variable to enter the basis.  Normally this
  !               value is set to NVARS.  This implies that no
  !               "partial pricing" is used.
  !     67        Adjust the tuning factor (normally one) to apply
  !               to the primal and dual error estimates.
  !     68        Pass  information to the  subprogram  FULMAT(),
  !               provided with the SPLP() package, so that a Fortran
  !               two-dimensional array can be used as the argument
  !               DATTRV(*).
  !     69        Pass an absolute tolerance to use for the feasibility
  !               test when the usual relative error test indicates
  !               infeasibility.  The nominal value of this tolerance,
  !               TOLABS, is zero.
  !
  !
  !  |---------------|
  !  |Working Arrays:|
  !  |---------------|
  !
  !     WORK(*),LW,
  !     IWORK(*),LIW
  !     ------------
  !     The arrays WORK(*) and IWORK(*) are respectively floating point
  !     and type INTEGER working arrays for SPLP( ) and its
  !     subprograms.  The lengths of these arrays are respectively
  !     LW and LIW.  These parameters must satisfy the inequalities
  !     noted above under the heading "Fortran Declarations Required."
  !     It is an error if either value is too small.
  !
  !  |----------------------------|
  !  |Input/Output files required:|
  !  |----------------------------|
  !
  !     Fortran unit 1 is used by SPLP( ) to store the sparse matrix A
  !     out of high-speed memory.  This direct access file is opened
  !     within the package under the following two conditions.
  !     1. When the Save/Restore feature is used.  2. When the
  !     constraint matrix is so large that storage out of high-speed
  !     memory is required.  The user may need to close unit 1
  !     (with deletion from the job step) in the main program unit
  !     when several calls are made to SPLP( ).  A crude
  !     upper bound for the amount of information written on unit 1
  !     is 6*nz, where nz is the number of nonzero entries in A.
  !     The unit number may be redefined to any other positive value
  !     by means of input in the option array PRGOPT(*).
  !
  !     Fortran unit 2 is used by SPLP( ) only when the Save/Restore
  !     feature is desired.  Normally this feature is not used.  It is
  !     activated by means of input in the option array PRGOPT(*).
  !     On some computer systems the user may need to open unit
  !     2 before executing a call to SPLP( ).  This file is type
  !     sequential and is unformatted.
  !
  !     Fortran unit=I1MACH(2) (check local setting) is used by SPLP( )
  !     when the printed output feature (KEY=51) is used.  Normally
  !     this feature is not used.  It is activated by input in the
  !     options array PRGOPT(*).  For many computer systems I1MACH(2)=6.
  !
  !  |-------|
  !  |Output:|
  !  |-------|
  !
  !     INFO,PRIMAL(*),DUALS(*)
  !     -----------------------
  !     The integer flag INFO indicates why SPLP( ) has returned to the
  !     user.  If INFO=1 the solution has been computed.  In this case
  !     X(J)=PRIMAL(J) and W(I)=PRIMAL(I+NVARS).  The dual variables
  !     for the equations A*x=w are in the array DUALS(I)=dual for
  !     equation number I.  The dual value for the component X(J) that
  !     has an upper or lower bound (or both) is returned in
  !     DUALS(J+MRELAS).  The only other values for INFO are .LT. 0.
  !     The meaning of these values can be found by reading
  !     the diagnostic message in the output file, or by looking for
  !     error number = (-INFO) under the heading "List of SPLP( ) Error
  !     and Diagnostic Messages."
  !     The diagnostic messages are printed using the error processing
  !     subprogram XERMSG( ) with error category LEVEL=1.
  !     See the document "Brief Instr. for Using the Sandia Math.
  !     Subroutine Library," SAND79-2382, Nov., 1980, for further inform-
  !     ation about resetting the usual response to a diagnostic message.
  !
  !     BL(*),BU(*),IND(*)
  !     ------------------
  !     These arrays are output parameters only under the (unusual)
  !     circumstances where the stated problem is infeasible, has an
  !     unbounded optimum value, or both.  These respective conditions
  !     correspond to INFO=-1,-2 or -3.  For INFO=-1 or -3 certain comp-
  !     onents of the vectors x or w will not satisfy the input bounds.
  !     If component J of X or component I of W does not satisfy its input
  !     bound because of infeasibility, then IND(J)=-4 or IND(I+NVARS)=-4,
  !     respectively.  For INFO=-2 or -3 certain
  !     components of the vector x could not be used as basic variables
  !     because the objective function would have become unbounded.
  !     In particular if component J of x corresponds to such a variable,
  !     then IND(J)=-3.  Further, if the input value of IND(J)
  !                      =1, then BU(J)=BL(J);
  !                      =2, then BL(J)=BU(J);
  !                      =4, then BL(J)=0.,BU(J)=0.
  !
  !     (The J-th variable in x has been restricted to an appropriate
  !     feasible value.)
  !     The negative output value for IND(*) allows the user to identify
  !     those constraints that are not satisfied or those variables that
  !     would cause unbounded values of the objective function.  Note
  !     that the absolute value of IND(*), together with BL(*) and BU(*),
  !     are valid input to SPLP( ).  In the case of infeasibility the
  !     sum of magnitudes of the infeasible values is minimized.  Thus
  !     one could reenter SPLP( ) with these components of x or w now
  !     fixed at their present values.  This involves setting
  !     the appropriate components of IND(*) = 3, and BL(*) = BU(*).
  !
  !     IBASIS(I),I=1,...,MRELAS
  !     ------------------------
  !     This array contains the indices of the variables that are
  !     in the active basis set at the solution (INFO=1).  A value
  !     of IBASIS(I) between 1 and NVARS corresponds to the variable
  !     X(IBASIS(I)).  A value of IBASIS(I) between NVARS+1 and NVARS+
  !     MRELAS corresponds to the variable W(IBASIS(I)-NVARS).
  !
  !     Computing with the Matrix A after Calling SPLP( )
  !     -------------------------------------------------
  !     Following the return from SPLP( ), nonzero entries of the MRELAS
  !     by NVARS matrix A are available for usage by the user.  The method
  !     for obtaining the next nonzero in column J with a row index
  !     strictly greater than I in value, is completed by executing
  !
  !         CALL PNNZRS(I,AIJ,IPLACE,WORK,IWORK,J)
  !
  !     The value of I is also an output parameter.  If I.LE.0 on output,
  !     then there are no more nonzeroes in column J.  If I.GT.0, the
  !     output value for component number I of column J is in AIJ.  The
  !     parameters WORK(*) and IWORK(*) are the same arguments as in the
  !     call to SPLP( ).  The parameter IPLACE is a single INTEGER
  !     working variable.
  !
  !     The data structure used for storage of the matrix A within SPLP( )
  !     corresponds to sequential storage by columns as defined in
  !     SAND78-0785.  Note that the names of the subprograms LNNZRS(),
  !     LCHNGS(),LINITM(),LLOC(),LRWPGE(), and LRWVIR() have been
  !     changed to PNNZRS(),PCHNGS(),PINITM(),IPLOC(),PRWPGE(), and
  !     PRWVIR() respectively.  The error processing subprogram LERROR()
  !     is no longer used; XERMSG() is used instead.
  !
  !  |-------------------------------|
  !  |Subprograms Required by SPLP( )|
  !  |-------------------------------|
  !     Called by SPLP() are SPLPMN(),SPLPUP(),SPINIT(),SPOPT(),
  !                          SPLPDM(),SPLPCE(),SPINCW(),SPLPFL(),
  !                          SPLPFE(),SPLPMU().
  !
  !     Error Processing Subprograms XERMSG(),I1MACH(),R1MACH()
  !
  !     Sparse Matrix Subprograms PNNZRS(),PCHNGS(),PRWPGE(),PRWVIR(),
  !                               PINITM(),IPLOC()
  !
  !     Mass Storage File Subprograms SOPENM(),SCLOSM(),SREADP(),SWRITP()
  !
  !     Basic Linear Algebra Subprograms SCOPY(),SASUM(),SDOT()
  !
  !     Sparse Matrix Basis Handling Subprograms LA05AS(),LA05BS(),
  !                                             LA05CS(),LA05ED(),MC20AS()
  !
  !     Vector Output Subprograms SVOUT(),IVOUT()
  !
  !     Machine-sensitive Subprograms I1MACH( ),R1MACH( ),
  !                                SOPENM(),SCLOSM(),SREADP(),SWRITP().
  !     COMMON Block Used
  !     -----------------
  !     /LA05DS/ SMALL,LP,LENL,LENU,NCP,LROW,LCOL
  !     See the document AERE-R8269 for further details.
  !    |------------------------|
  !    |Example of SPLP( ) Usage|
  !    |------------------------|
  !     PROGRAM LPEX
  !     THE OPTIMIZATION PROBLEM IS TO FIND X1, X2, X3 THAT
  !     MINIMIZE X1 + X2 + X3, X1.GE.0, X2.GE.0, X3 UNCONSTRAINED.
  !
  !     THE UNKNOWNS X1,X2,X3 ARE TO SATISFY CONSTRAINTS
  !
  !        X1 -3*X2 +4*X3 = 5
  !        X1 -2*X2     .LE.3
  !            2*X2 - X3.GE.4
  !
  !     WE FIRST DEFINE THE DEPENDENT VARIABLES
  !          W1=X1 -3*X2 +4*X3
  !          W2=X1- 2*X2
  !          W3=    2*X2 -X3
  !
  !     WE NOW SHOW HOW TO USE SPLP( ) TO SOLVE THIS LINEAR OPTIMIZATION
  !     PROBLEM.  EACH REQUIRED STEP WILL BE SHOWN IN THIS EXAMPLE.
  !     DIMENSION COSTS(03),PRGOPT(01),DATTRV(18),BL(06),BU(06),IND(06),
  !    *PRIMAL(06),DUALS(06),IBASIS(06),WORK(079),IWORK(103)
  !
  !     EXTERNAL USRMAT
  !     MRELAS=3
  !     NVARS=3
  !
  !     DEFINE THE ARRAY COSTS(*) FOR THE OBJECTIVE FUNCTION.
  !     COSTS(01)=1.
  !     COSTS(02)=1.
  !     COSTS(03)=1.
  !
  !     PLACE THE NONZERO INFORMATION ABOUT THE MATRIX IN DATTRV(*).
  !     DEFINE COL. 1:
  !     DATTRV(01)=-1
  !     DATTRV(02)=1
  !     DATTRV(03)=1.
  !     DATTRV(04)=2
  !     DATTRV(05)=1.
  !
  !     DEFINE COL. 2:
  !     DATTRV(06)=-2
  !     DATTRV(07)=1
  !     DATTRV(08)=-3.
  !     DATTRV(09)=2
  !     DATTRV(10)=-2.
  !     DATTRV(11)=3
  !     DATTRV(12)=2.
  !
  !     DEFINE COL. 3:
  !     DATTRV(13)=-3
  !     DATTRV(14)=1
  !     DATTRV(15)=4.
  !     DATTRV(16)=3
  !     DATTRV(17)=-1.
  !
  !     DATTRV(18)=0
  !
  !     CONSTRAIN X1,X2 TO BE NONNEGATIVE. LET X3 HAVE NO BOUNDS.
  !     BL(1)=0.
  !     IND(1)=1
  !     BL(2)=0.
  !     IND(2)=1
  !     IND(3)=4
  !
  !     CONSTRAIN W1=5,W2.LE.3, AND W3.GE.4.
  !     BL(4)=5.
  !     BU(4)=5.
  !     IND(4)=3
  !     BU(5)=3.
  !     IND(5)=2
  !     BL(6)=4.
  !     IND(6)=1
  !
  !     INDICATE THAT NO MODIFICATIONS TO OPTIONS ARE IN USE.
  !     PRGOPT(01)=1
  !
  !     DEFINE THE WORKING ARRAY LENGTHS.
  !     LW=079
  !     LIW=103
  !     CALL SPLP(USRMAT,MRELAS,NVARS,COSTS,PRGOPT,DATTRV,
  !    *BL,BU,IND,INFO,PRIMAL,DUALS,IBASIS,WORK,LW,IWORK,LIW)
  !
  !     CALCULATE VAL, THE MINIMAL VALUE OF THE OBJECTIVE FUNCTION.
  !     VAL=SDOT(NVARS,COSTS,1,PRIMAL,1)
  !
  !     STOP
  !     END
  !    |------------------------|
  !    |End of Example of Usage |
  !    |------------------------|
  !
  !    |------------------------------------|
  !    |Usage of SPLP( ) Subprogram Options.|
  !    |------------------------------------|
  !
  !     Users frequently have a large variety of requirements for linear
  !     optimization software.  Allowing for these varied requirements
  !     is at cross purposes with the desire to keep the usage of SPLP( )
  !     as simple as possible. One solution to this dilemma is as follows.
  !     (1) Provide a version of SPLP( ) that solves a wide class of
  !     problems and is easy to use. (2) Identify parameters within SPLP()
  !     that certain users may want to change.  (3) Provide a means
  !     of changing any selected number of these parameters that does
  !     not require changing all of them.
  !
  !     Changing selected parameters is done by requiring
  !     that the user provide an option array, PRGOPT(*), to SPLP( ).
  !     The contents of PRGOPT(*) inform SPLP( ) of just those options
  !     that are going to be modified within the total set of possible
  !     parameters that can be modified.  The array PRGOPT(*) is a linked
  !     list consisting of groups of data of the following form
  !
  !          LINK
  !          KEY
  !          SWITCH
  !          data set
  !
  !     that describe the desired options.  The parameters LINK, KEY and
  !     switch are each one word and are always required.  The data set
  !     can be comprised of several words or can be empty.  The number of
  !     words in the data set for each option depends on the value of
  !     the parameter KEY.
  !
  !     The value of LINK points to the first entry of the next group
  !     of data within PRGOPT(*).  The exception is when there are no more
  !     options to change.  In that case, LINK=1 and the values for KEY,
  !     SWITCH and data set are not referenced.  The general layout of
  !     PRGOPT(*) is as follows:
  !          ...PRGOPT(1)=LINK1 (link to first entry of next group)
  !          .  PRGOPT(2)=KEY1 (KEY to the option change)
  !          .  PRGOPT(3)=SWITCH1 (on/off switch for the option)
  !          .  PRGOPT(4)=data value
  !          .       .
  !          .       .
  !          .       .
  !          ...PRGOPT(LINK1)=LINK2 (link to first entry of next group)
  !          .  PRGOPT(LINK1+1)=KEY2 (KEY to option change)
  !          .  PRGOPT(LINK1+2)=SWITCH2 (on/off switch for the option)
  !          .  PRGOPT(LINK1+3)=data value
  !          ...     .
  !          .       .
  !          .       .
  !          ...PRGOPT(LINK)=1 (no more options to change)
  !
  !     A value of LINK that is .LE.0 or .GT. 10000 is an error.
  !     In this case SPLP( ) returns with an error message, INFO=-14.
  !     This helps prevent using invalid but positive values of LINK that
  !     will probably extend beyond the program limits of PRGOPT(*).
  !     Unrecognized values of KEY are ignored.  If the value of SWITCH is
  !     zero then the option is turned off.  For any other value of SWITCH
  !     the option is turned on.  This is used to allow easy changing of
  !     options without rewriting PRGOPT(*).  The order of the options is
  !     arbitrary and any number of options can be changed with the
  !     following restriction.  To prevent cycling in processing of the
  !     option array PRGOPT(*), a count of the number of options changed
  !     is maintained.  Whenever this count exceeds 1000 an error message
  !     (INFO=-15) is printed and the subprogram returns.
  !
  !     In the following description of the options, the value of
  !     LATP indicates the amount of additional storage that a particular
  !     option requires.  The sum of all of these values (plus one) is
  !     the minimum dimension for the array PRGOPT(*).
  !
  !     If a user is satisfied with the nominal form of SPLP( ),
  !     set PRGOPT(1)=1 (or PRGOPT(1)=1.E0).
  !
  !     Options:
  !
  ! -----KEY = 50.  Change from a minimization problem to a maximization
  !     problem.
  !     If SWITCH=0  option is off; solve minimization problem.
  !              =1  option is on; solve maximization problem.
  !     data set =empty
  !     LATP=3
  !
  ! -----KEY = 51.  Change the amount of printed output.  The nominal form
  !     of SPLP( ) has no printed output.
  !     The first level of output (SWITCH=1) includes
  !
  !     (1) Minimum dimensions for the arrays COSTS(*),BL(*),BU(*),IND(*),
  !         PRIMAL(*),DUALS(*),IBASIS(*), and PRGOPT(*).
  !     (2) Problem dimensions MRELAS,NVARS.
  !     (3) The types of and values for the bounds on x and w,
  !         and the values of the components of the vector costs.
  !     (4) Whether optimization problem is minimization or
  !         maximization.
  !     (5) Whether steepest edge or smallest reduced cost criteria used
  !         for exchanging variables in the revised simplex method.
  !
  !     Whenever a solution has been found, (INFO=1),
  !
  !     (6) the value of the objective function,
  !     (7) the values of the vectors x and w,
  !     (8) the dual variables for the constraints A*x=w and the
  !         bounded components of x,
  !     (9) the indices of the basic variables,
  !    (10) the number of revised simplex method iterations,
  !    (11) the number of full decompositions of the basis matrix.
  !
  !     The second level of output (SWITCH=2) includes all for SWITCH=1
  !     plus
  !
  !    (12) the iteration number,
  !    (13) the column number to enter the basis,
  !    (14) the column number to leave the basis,
  !    (15) the length of the step taken.
  !
  !     The third level of output (SWITCH=3) includes all for SWITCH=2
  !     plus
  !    (16) critical quantities required in the revised simplex method.
  !          This output is rather voluminous.  It is intended to be used
  !          as a diagnostic tool in case of a failure in SPLP( ).
  !
  !     If SWITCH=0  option is off; no printed output.
  !              =1  summary output.
  !              =2  lots of output.
  !              =3  even more output.
  !     data set =empty
  !     LATP=3
  !
  ! -----KEY = 52.  Redefine the parameter, IDIGIT, which determines the
  !     format and precision used for the printed output.  In the printed
  !     output, at least ABS(IDIGIT) decimal digits per number is printed.
  !     If IDIGIT.LT.0, 72 printing columns are used.  IF IDIGIT.GT.0, 133
  !     printing columns are used.
  !     If SWITCH=0  option is off; IDIGIT=-4.
  !              =1  option is on.
  !     data set =IDIGIT
  !     LATP=4
  !
  ! -----KEY = 53.  Redefine LAMAT and LBM, the lengths of the portions of
  !     WORK(*) and IWORK(*) that are allocated to the sparse matrix
  !     storage and the sparse linear equation solver, respectively.
  !     LAMAT must be .GE. NVARS+7 and LBM must be positive.
  !     If SWITCH=0  option is off; LAMAT=4*NVARS+7
  !                                 LBM  =8*MRELAS.
  !              =1  option is on.
  !     data set =LAMAT
  !               LBM
  !     LATP=5
  !
  ! -----KEY = 54. Redefine IPAGEF, the file number where the pages of the
  !     sparse data matrix are stored.  IPAGEF must be positive and
  !     different from ISAVE (see option 56).
  !     If SWITCH=0  option is off; IPAGEF=1.
  !              =1  option is on.
  !     data set =IPAGEF
  !     LATP=4
  !
  ! -----KEY = 55.  Partial results have been computed and stored on unit
  !     number ISAVE (see option 56), during a previous run of
  !     SPLP( ).  This is a continuation from these partial results.
  !     The arrays COSTS(*),BL(*),BU(*),IND(*) do not have to have
  !     the same values as they did when the checkpointing occurred.
  !     This feature makes it possible for the user to do certain
  !     types of parameter studies such as changing costs and varying
  !     the constraints of the problem.  This file is rewound both be-
  !     fore and after reading the partial results.
  !     If SWITCH=0  option is off; start a new problem.
  !              =1  option is on; continue from partial results
  !                  that are stored in file ISAVE.
  !     data set = empty
  !     LATP=3
  !
  ! -----KEY = 56.  Redefine ISAVE, the file number where the partial
  !     results are stored (see option 57).  ISAVE must be positive and
  !     different from IPAGEF (see option 54).
  !     If SWITCH=0  option is off; ISAVE=2.
  !              =1  option is on.
  !     data set =ISAVE
  !     LATP=4
  !
  ! -----KEY = 57.  Save the partial results after maximum number of
  !     iterations, MAXITR, or at the optimum.  When this option is on,
  !     data essential to continuing the calculation is saved on a file
  !     using a Fortran binary write operation.  The data saved includes
  !     all the information about the sparse data matrix A.  Also saved
  !     is information about the current basis.  Nominally the partial
  !     results are saved on Fortran unit 2.  This unit number can be
  !     redefined (see option 56).  If the save option is on,
  !     this file must be opened (or declared) by the user prior to the
  !     call to SPLP( ).  A crude upper bound for the number of words
  !     written to this file is 6*nz.  Here nz= number of nonzeros in A.
  !     If SWITCH=0  option is off; do not save partial results.
  !              =1  option is on; save partial results.
  !     data set = empty
  !     LATP=3
  !
  ! -----KEY = 58.  Redefine the maximum number of iterations, MAXITR, to
  !     be taken before returning to the user.
  !     If SWITCH=0  option is off; MAXITR=3*(NVARS+MRELAS).
  !              =1  option is on.
  !     data set =MAXITR
  !     LATP=4
  !
  ! -----KEY = 59.  Provide SPLP( ) with exactly MRELAS indices which
  !     comprise a feasible, nonsingular basis.  The basis must define a
  !     feasible point: values for x and w such that A*x=w and all the
  !     stated bounds on x and w are satisfied.  The basis must also be
  !     nonsingular.  The failure of either condition will cause an error
  !     message (INFO=-23 or =-24, respectively).  Normally, SPLP( ) uses
  !     identity matrix columns which correspond to the components of w.
  !     This option would normally not be used when restarting from
  !     a previously saved run (KEY=57).
  !     In numbering the unknowns,
  !     the components of x are numbered (1-NVARS) and the components
  !     of w are numbered (NVARS+1)-(NVARS+MRELAS).  A value for an
  !     index .LE. 0 or .GT. (NVARS+MRELAS) is an error (INFO=-16).
  !     If SWITCH=0  option is off; SPLP( ) chooses the initial basis.
  !              =1  option is on; user provides the initial basis.
  !     data set =MRELAS indices of basis; order is arbitrary.
  !     LATP=MRELAS+3
  !
  ! -----KEY = 60.  Provide the scale factors for the columns of the data
  !     matrix A.  Normally, SPLP( ) computes the scale factors as the
  !     reciprocals of the max. norm of each column.
  !     If SWITCH=0  option is off; SPLP( ) computes the scale factors.
  !              =1  option is on; user provides the scale factors.
  !     data set =scaling for column J, J=1,NVARS; order is sequential.
  !     LATP=NVARS+3
  !
  ! -----KEY = 61.  Provide a scale factor, COSTSC, for the vector of
  !     costs.  Normally, SPLP( ) computes this scale factor to be the
  !     reciprocal of the max. norm of the vector costs after the column
  !     scaling has been applied.
  !     If SWITCH=0  option is off; SPLP( ) computes COSTSC.
  !              =1  option is on; user provides COSTSC.
  !     data set =COSTSC
  !     LATP=4
  !
  ! -----KEY = 62.  Provide size parameters, ASMALL and ABIG, the smallest
  !     and largest magnitudes of nonzero entries in the data matrix A,
  !     respectively.  When this option is on, SPLP( ) will check the
  !     nonzero entries of A to see if they are in the range of ASMALL and
  !     ABIG.  If an entry of A is not within this range, SPLP( ) returns
  !     an error message, INFO=-22. Both ASMALL and ABIG must be positive
  !     with ASMALL .LE. ABIG.  Otherwise,  an error message is returned,
  !     INFO=-17.
  !     If SWITCH=0  option is off; no checking of the data matrix is done
  !              =1  option is on; checking is done.
  !     data set =ASMALL
  !               ABIG
  !     LATP=5
  !
  ! -----KEY = 63.  Redefine the relative tolerance, TOLLS, used in
  !     checking if the residuals are feasible.  Normally,
  !     TOLLS=RELPR, where RELPR is the machine precision.
  !     If SWITCH=0  option is off; TOLLS=RELPR.
  !              =1  option is on.
  !     data set =TOLLS
  !     LATP=4
  !
  ! -----KEY = 64. Use the minimum reduced cost pricing strategy to choose
  !     columns to enter the basis.  Normally, SPLP( ) uses the steepest
  !     edge pricing strategy which is the best local move.  The steepest
  !     edge pricing strategy generally uses fewer iterations than the
  !     minimum reduced cost pricing, but each iteration costs more in the
  !     number of calculations done.  The steepest edge pricing is
  !     considered to be more efficient.  However, this is very problem
  !     dependent.  That is why SPLP( ) provides the option of either
  !     pricing strategy.
  !     If SWITCH=0  option is off; steepest option edge pricing is used.
  !              =1  option is on; minimum reduced cost pricing is used.
  !     data set =empty
  !     LATP=3
  !
  ! -----KEY = 65.  Redefine MXITBR, the number of iterations between
  !     recalculating the error in the primal solution.  Normally, MXITBR
  !     is set to 10.  The error in the primal solution is used to monitor
  !     the error in solving the linear system.  This is an expensive
  !     calculation and every tenth iteration is generally often enough.
  !     If SWITCH=0  option is off; MXITBR=10.
  !              =1  option is on.
  !     data set =MXITBR
  !     LATP=4
  !
  ! -----KEY = 66.  Redefine NPP, the number of negative reduced costs
  !     (at most) to be found at each iteration of choosing
  !     a variable to enter the basis.  Normally NPP is set
  !     to NVARS which implies that all of the reduced costs
  !     are computed at each such step.  This "partial
  !     pricing" may very well increase the total number
  !     of iterations required.  However it decreases the
  !     number of calculations at each iteration.
  !     therefore the effect on overall efficiency is quite
  !     problem-dependent.
  !
  !     if SWITCH=0 option is off; NPP=NVARS
  !              =1 option is on.
  !     data set =NPP
  !     LATP=4
  !
  ! -----KEY =  67.  Redefine the tuning factor (PHI) used to scale the
  !     error estimates for the primal and dual linear algebraic systems
  !     of equations.  Normally, PHI = 1.E0, but in some environments it
  !     may be necessary to reset PHI to the range 0.001-0.01.  This is
  !     particularly important for machines with short word lengths.
  !
  !     if SWITCH = 0 option is off; PHI=1.E0.
  !               = 1 option is on.
  !     Data Set  = PHI
  !     LATP=4
  !
  ! -----KEY = 68.  Used together with the subprogram FULMAT(), provided
  !     with the SPLP() package, for passing a standard Fortran two-
  !     dimensional array containing the constraint matrix.  Thus the sub-
  !     program FULMAT must be declared in a Fortran EXTERNAL statement.
  !     The two-dimensional array is passed as the argument DATTRV.
  !     The information about the array and problem dimensions are passed
  !     in the option array PRGOPT(*).  It is an error if FULMAT() is
  !     used and this information is not passed in PRGOPT(*).
  !
  !     if SWITCH = 0 option is off; this is an error is FULMAT() is
  !                                  used.
  !               = 1 option is on.
  !     Data Set  = IA = row dimension of two-dimensional array.
  !                 MRELAS = number of constraint equations.
  !                 NVARS  = number of dependent variables.
  !     LATP = 6
  ! -----KEY = 69.  Normally a relative tolerance (TOLLS, see option 63)
  !     is used to decide if the problem is feasible.  If this test fails
  !     an absolute test will be applied using the value TOLABS.
  !     Nominally TOLABS = zero.
  !     If SWITCH = 0 option is off; TOLABS = zero.
  !               = 1 option is on.
  !     Data set  = TOLABS
  !     LATP = 4
  !
  !    |-----------------------------|
  !    |Example of Option array Usage|
  !    |-----------------------------|
  !     To illustrate the usage of the option array, let us suppose that
  !     the user has the following nonstandard requirements:
  !
  !          a) Wants to change from minimization to maximization problem.
  !          b) Wants to limit the number of simplex steps to 100.
  !          c) Wants to save the partial results after 100 steps on
  !             Fortran unit 2.
  !
  !     After these 100 steps are completed the user wants to continue the
  !     problem (until completed) using the partial results saved on
  !     Fortran unit 2.  Here are the entries of the array PRGOPT(*)
  !     that accomplish these tasks.  (The definitions of the other
  !     required input parameters are not shown.)
  !
  !     CHANGE TO A MAXIMIZATION PROBLEM; KEY=50.
  !     PRGOPT(01)=4
  !     PRGOPT(02)=50
  !     PRGOPT(03)=1
  !
  !     LIMIT THE NUMBER OF SIMPLEX STEPS TO 100; KEY=58.
  !     PRGOPT(04)=8
  !     PRGOPT(05)=58
  !     PRGOPT(06)=1
  !     PRGOPT(07)=100
  !
  !     SAVE THE PARTIAL RESULTS, AFTER 100 STEPS, ON FORTRAN
  !     UNIT 2; KEY=57.
  !     PRGOPT(08)=11
  !     PRGOPT(09)=57
  !     PRGOPT(10)=1
  !
  !     NO MORE OPTIONS TO CHANGE.
  !     PRGOPT(11)=1
  !     The user makes the CALL statement for SPLP( ) at this point.
  !     Now to restart, using the partial results after 100 steps, define
  !     new values for the array PRGOPT(*):
  !
  !     AGAIN INFORM SPLP( ) THAT THIS IS A MAXIMIZATION PROBLEM.
  !     PRGOPT(01)=4
  !     PRGOPT(02)=50
  !     PRGOPT(03)=1
  !
  !     RESTART, USING SAVED PARTIAL RESULTS; KEY=55.
  !     PRGOPT(04)=7
  !     PRGOPT(05)=55
  !     PRGOPT(06)=1
  !
  !     NO MORE OPTIONS TO CHANGE.  THE SUBPROGRAM SPLP( ) IS NO LONGER
  !     LIMITED TO 100 SIMPLEX STEPS BUT WILL RUN UNTIL COMPLETION OR
  !     MAX.=3*(MRELAS+NVARS) ITERATIONS.
  !     PRGOPT(07)=1
  !     The user now makes a CALL to subprogram SPLP( ) to compute the
  !     solution.
  !    |-------------------------------------------|
  !    |End of Usage of SPLP( ) Subprogram Options.|
  !    |-------------------------------------------|
  !
  !     |----------------------------------------------|
  !     |List of SPLP( ) Error and Diagnostic Messages.|
  !     |----------------------------------------------|
  !      This section may be required to understand the meanings of the
  !     error flag =-INFO  that may be returned from SPLP( ).
  !
  ! -----1. There is no set of values for x and w that satisfy A*x=w and
  !     the stated bounds.  The problem can be made feasible by ident-
  !     ifying components of w that are now infeasible and then rede-
  !     signating them as free variables.  Subprogram SPLP( ) only
  !     identifies an infeasible problem; it takes no other action to
  !     change this condition.  Message:
  !     SPLP( ). THE PROBLEM APPEARS TO BE INFEASIBLE.
  !     ERROR NUMBER =         1
  !
  !     2. One of the variables in either the vector x or w was con-
  !     strained at a bound.  Otherwise the objective function value,
  !     (transpose of costs)*x, would not have a finite optimum.
  !     Message:
  !     SPLP( ). THE PROBLEM APPEARS TO HAVE NO FINITE SOLN.
  !     ERROR NUMBER =         2
  !
  !     3.  Both of the conditions of 1. and 2. above have occurred.
  !     Message:
  !     SPLP( ). THE PROBLEM APPEARS TO BE INFEASIBLE AND TO
  !     HAVE NO FINITE SOLN.
  !     ERROR NUMBER =         3
  !
  ! -----4.  The REAL and INTEGER working arrays, WORK(*) and IWORK(*),
  !     are not long enough. The values (I1) and (I2) in the message
  !     below will give you the minimum length required.  Also redefine
  !     LW and LIW, the lengths of these arrays.  Message:
  !     SPLP( ). WORK OR IWORK IS NOT LONG ENOUGH. LW MUST BE (I1)
  !     AND LIW MUST BE (I2).
  !               IN ABOVE MESSAGE, I1=         0
  !               IN ABOVE MESSAGE, I2=         0
  !     ERROR NUMBER =        4
  !
  ! -----5. and 6.  These error messages often mean that one or more
  !     arguments were left out of the call statement to SPLP( ) or
  !     that the values of MRELAS and NVARS have been over-written
  !     by garbage.  Messages:
  !     SPLP( ). VALUE OF MRELAS MUST BE .GT.0. NOW=(I1).
  !               IN ABOVE MESSAGE, I1=         0
  !     ERROR NUMBER =         5
  !
  !     SPLP( ). VALUE OF NVARS MUST BE .GT.0. NOW=(I1).
  !               IN ABOVE MESSAGE, I1=         0
  !     ERROR NUMBER =         6
  !
  ! -----7.,8., and 9.  These error messages can occur as the data matrix
  !     is being defined by either USRMAT( ) or the user-supplied sub-
  !     program, 'NAME'( ). They would indicate a mistake in the contents
  !     of DATTRV(*), the user-written subprogram or that data has been
  !     over-written.
  !     Messages:
  !     SPLP( ). MORE THAN 2*NVARS*MRELAS ITERS. DEFINING OR UPDATING
  !     MATRIX DATA.
  !     ERROR NUMBER =        7
  !
  !     SPLP( ). ROW INDEX (I1) OR COLUMN INDEX (I2) IS OUT OF RANGE.
  !               IN ABOVE MESSAGE, I1=         1
  !               IN ABOVE MESSAGE, I2=        12
  !     ERROR NUMBER =        8
  !
  !     SPLP( ). INDICATION FLAG (I1) FOR MATRIX DATA MUST BE
  !     EITHER 0 OR 1.
  !               IN ABOVE MESSAGE, I1=        12
  !     ERROR NUMBER =        9
  !
  ! -----10. and 11.  The type of bound (even no bound) and the bounds
  !     must be specified for each independent variable. If an independent
  !     variable has both an upper and lower bound, the bounds must be
  !     consistent.  The lower bound must be .LE. the upper bound.
  !     Messages:
  !     SPLP( ). INDEPENDENT VARIABLE (I1) IS NOT DEFINED.
  !               IN ABOVE MESSAGE, I1=         1
  !     ERROR NUMBER =        10
  !
  !     SPLP( ).  LOWER BOUND (R1) AND UPPER BOUND (R2) FOR INDEP.
  !     VARIABLE (I1) ARE NOT CONSISTENT.
  !               IN ABOVE MESSAGE, I1=         1
  !               IN ABOVE MESSAGE, R1=    0.
  !               IN ABOVE MESSAGE, R2=    -.1000000000E+01
  !     ERROR NUMBER =        11
  !
  ! -----12. and 13.  The type of bound (even no bound) and the bounds
  !     must be specified for each dependent variable.  If a dependent
  !     variable has both an upper and lower bound, the bounds must be
  !     consistent. The lower bound must be .LE. the upper bound.
  !     Messages:
  !     SPLP( ). DEPENDENT VARIABLE (I1) IS NOT DEFINED.
  !               IN ABOVE MESSAGE, I1=         1
  !     ERROR NUMBER =        12
  !
  !     SPLP( ).  LOWER BOUND (R1) AND UPPER BOUND (R2) FOR DEP.
  !      VARIABLE (I1) ARE NOT CONSISTENT.
  !               IN ABOVE MESSAGE, I1=         1
  !               IN ABOVE MESSAGE, R1=    0.
  !               IN ABOVE MESSAGE, R2=    -.1000000000E+01
  !     ERROR NUMBER =        13
  !
  ! -----14. - 21.  These error messages can occur when processing the
  !     option array, PRGOPT(*), supplied by the user.  They would
  !     indicate a mistake in defining PRGOPT(*) or that data has been
  !     over-written.  See heading Usage of SPLP( )
  !     Subprogram Options, for details on how to define PRGOPT(*).
  !     Messages:
  !     SPLP( ). THE USER OPTION ARRAY HAS UNDEFINED DATA.
  !     ERROR NUMBER =        14
  !
  !     SPLP( ). OPTION ARRAY PROCESSING IS CYCLING.
  !     ERROR NUMBER =        15
  !
  !     SPLP( ). AN INDEX OF USER-SUPPLIED BASIS IS OUT OF RANGE.
  !     ERROR NUMBER =        16
  !
  !     SPLP( ). SIZE PARAMETERS FOR MATRIX MUST BE SMALLEST AND LARGEST
  !     MAGNITUDES OF NONZERO ENTRIES.
  !     ERROR NUMBER =        17
  !
  !     SPLP( ). THE NUMBER OF REVISED SIMPLEX STEPS BETWEEN CHECK-POINTS
  !     MUST BE POSITIVE.
  !     ERROR NUMBER =        18
  !
  !     SPLP( ). FILE NUMBERS FOR SAVED DATA AND MATRIX PAGES MUST BE
  !     POSITIVE AND NOT EQUAL.
  !     ERROR NUMBER =        19
  !
  !     SPLP( ). USER-DEFINED VALUE OF LAMAT (I1)
  !     MUST BE .GE. NVARS+7.
  !               IN ABOVE MESSAGE, I1=         1
  !     ERROR NUMBER =         20
  !
  !     SPLP( ). USER-DEFINED VALUE OF LBM MUST BE .GE. 0.
  !     ERROR NUMBER =         21
  !
  ! -----22.  The user-option, number 62, to check the size of the matrix
  !     data has been used.  An element of the matrix does not lie within
  !     the range of ASMALL and ABIG, parameters provided by the user.
  !     (See the heading: Usage of SPLP( ) Subprogram Options,
  !     for details about this feature.)  Message:
  !     SPLP( ). A MATRIX ELEMENT'S SIZE IS OUT OF THE SPECIFIED RANGE.
  !     ERROR NUMBER =        22
  !
  ! -----23.  The user has provided an initial basis that is singular.
  !     In this case, the user can remedy this problem by letting
  !     subprogram SPLP( ) choose its own initial basis.  Message:
  !     SPLP( ). A SINGULAR INITIAL BASIS WAS ENCOUNTERED.
  !     ERROR NUMBER =         23
  !
  ! -----24.  The user has provided an initial basis which is infeasible.
  !     The x and w values it defines do not satisfy A*x=w and the stated
  !     bounds.  In this case, the user can let subprogram SPLP( )
  !     choose its own initial basis.  Message:
  !     SPLP( ). AN INFEASIBLE INITIAL BASIS WAS ENCOUNTERED.
  !     ERROR NUMBER =        24
  !
  ! -----25. Subprogram SPLP( ) has completed the maximum specified number
  !     of iterations.  (The nominal maximum number is 3*(MRELAS+NVARS).)
  !     The results, necessary to continue on from
  !     this point, can be saved on Fortran unit 2 by activating option
  !     KEY=57.  If the user anticipates continuing the calculation, then
  !     the contents of Fortran unit 2 must be retained intact.  This
  !     is not done by subprogram SPLP( ), so the user needs to save unit
  !     2 by using the appropriate system commands.  Message:
  !     SPLP( ). MAX. ITERS. (I1) TAKEN. UP-TO-DATE RESULTS
  !     SAVED ON FILE (I2). IF(I2)=0, NO SAVE.
  !               IN ABOVE MESSAGE, I1=       500
  !               IN ABOVE MESSAGE, I2=         2
  !     ERROR NUMBER =        25
  !
  ! -----26.  This error should never happen.  Message:
  !     SPLP( ). MOVED TO A SINGULAR POINT. THIS SHOULD NOT HAPPEN.
  !     ERROR NUMBER =        26
  !
  ! -----27.  The subprogram LA05A( ), which decomposes the basis matrix,
  !     has returned with an error flag (R1).  (See the document,
  !     "Fortran subprograms for handling sparse linear programming
  !     bases", AERE-R8269, J.K. Reid, Jan., 1976, H.M. Stationery Office,
  !     for an explanation of this error.)  Message:
  !     SPLP( ). LA05A( ) RETURNED ERROR FLAG (R1) BELOW.
  !               IN ABOVE MESSAGE, R1=    -.5000000000E+01
  !     ERROR NUMBER =        27
  !
  ! -----28.  The sparse linear solver package, LA05*( ), requires more
  !     space.  The value of LBM must be increased.  See the companion
  !     document, Usage of SPLP( ) Subprogram Options, for details on how
  !     to increase the value of LBM.  Message:
  !     SPLP( ). SHORT ON STORAGE FOR LA05*( ) PACKAGE. USE PRGOPT(*)
  !     TO GIVE MORE.
  !     ERROR NUMBER =        28
  !
  ! -----29.  The row dimension of the two-dimensional Fortran array,
  !     the number of constraint equations (MRELAS), and the number
  !     of variables (NVARS), were not passed to the subprogram
  !     FULMAT().  See KEY = 68 for details.  Message:
  !     FULMAT() OF SPLP() PACKAGE. ROW DIM., MRELAS, NVARS ARE
  !     MISSING FROM PRGOPT(*).
  !     ERROR NUMBER =        29
  !
  !     |------------------------------------------------------|
  !     |End of List of SPLP( ) Error and Diagnostic Messages. |
  !     |------------------------------------------------------|
  !***REFERENCES  R. J. Hanson and K. L. Hiebert, A sparse linear
  !                 programming subprogram, Report SAND81-0297, Sandia
  !                 National Laboratories, 1981.
  !***ROUTINES CALLED  SPLPMN, XERMSG
  !***REVISION HISTORY  (YYMMDD)
  !   811215  DATE WRITTEN
  !   890605  Corrected references to XERRWV.  (WRB)
  !   890605  Removed unreferenced labels.  (WRB)
  !   890605  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900510  Convert XERRWV calls to XERMSG calls.  (RWC)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  !***END PROLOGUE  SPLP
  INTEGER iadbig, ictmax, ictopt, Info, iopt, key, lamat, last, &
    lbasma, lbm, lcolnr, lcsc, lerd, lerp, libb, librc, &
    limat, lipr, Liw, liwork
  INTEGER liwr, lmx, lrg, lrhs, lrprim, lrz, Lw, lwork, lwr, lww, &
    Mrelas, nerr, next, Nvars
  REAL Bl(*), Bu(*), Costs(*), Dattrv(*), Duals(*), Prgopt(*), &
    Primal(*), Work(*), zero
  !
  INTEGER Ibasis(*), Ind(*), Iwork(*)
  CHARACTER(8) :: xern1, xern2
  !
  EXTERNAL USRMAT
  !
  !***FIRST EXECUTABLE STATEMENT  SPLP
  zero = 0.E0
  iopt = 1
  !
  !     VERIFY THAT MRELAS, NVARS .GT. 0.
  !
  IF ( Mrelas<=0 ) THEN
    WRITE (xern1,'(I8)') Mrelas
    CALL XERMSG('SLATEC','SPLP','VALUE OF MRELAS MUST BE .GT. 0.  NOW = '//xern1,5,1)
    Info = -5
    RETURN
  ENDIF
  !
  IF ( Nvars<=0 ) THEN
    WRITE (xern1,'(I8)') Nvars
    CALL XERMSG('SLATEC','SPLP','VALUE OF NVARS MUST BE .GT. 0.  NOW = '//xern1,6,1)
    Info = -6
    RETURN
  ENDIF
  !
  lmx = 4*Nvars + 7
  lbm = 8*Mrelas
  last = 1
  iadbig = 10000
  ictmax = 1000
  ictopt = 0
  DO
    !
    !     LOOK IN OPTION ARRAY FOR CHANGES TO WORK ARRAY LENGTHS.
    next = INT( Prgopt(last) )
    IF ( next<=0.OR.next>iadbig ) THEN
      !
      !     THE CHECKS FOR SMALL OR LARGE VALUES OF NEXT ARE TO PREVENT
      !     WORKING WITH UNDEFINED DATA.
      nerr = 14
      CALL XERMSG('SLATEC','SPLP',&
        'THE USER OPTION ARRAY HAS UNDEFINED DATA.',nerr,iopt)
      Info = -nerr
      RETURN
    ELSEIF ( next==1 ) THEN
      !
      !     CHECK LENGTH VALIDITY OF SPARSE MATRIX STAGING AREA.
      !
      IF ( lmx<Nvars+7 ) THEN
        WRITE (xern1,'(I8)') lmx
        CALL XERMSG('SLATEC','SPLP','USER-DEFINED VALUE OF LAMAT = '//&
          xern1//' MUST BE .GE. NVARS+7.',20,1)
        Info = -20
        RETURN
      ENDIF
      !
      !     TRIVIAL CHECK ON LENGTH OF LA05*() MATRIX AREA.
      IF ( lbm<0 ) EXIT
      !
      !     DEFINE POINTERS FOR STARTS OF SUBARRAYS USED IN WORK(*)
      !     AND IWORK(*) IN OTHER SUBPROGRAMS OF THE PACKAGE.
      lamat = 1
      lcsc = lamat + lmx
      lcolnr = lcsc + Nvars
      lerd = lcolnr + Nvars
      lerp = lerd + Mrelas
      lbasma = lerp + Mrelas
      lwr = lbasma + lbm
      lrz = lwr + Mrelas
      lrg = lrz + Nvars + Mrelas
      lrprim = lrg + Nvars + Mrelas
      lrhs = lrprim + Mrelas
      lww = lrhs + Mrelas
      lwork = lww + Mrelas - 1
      limat = 1
      libb = limat + lmx
      librc = libb + Nvars + Mrelas
      lipr = librc + 2*lbm
      liwr = lipr + 2*Mrelas
      liwork = liwr + 8*Mrelas - 1
      !
      !     CHECK ARRAY LENGTH VALIDITY OF WORK(*), IWORK(*).
      !
      IF ( Lw<lwork.OR.Liw<liwork ) THEN
        WRITE (xern1,'(I8)') lwork
        WRITE (xern2,'(I8)') liwork
        CALL XERMSG('SLATEC','SPLP','WORK OR IWORK IS NOT LONG ENOUGH. LW MUST BE = '//xern1//' AND LIW MUST BE = '//&
          xern2,4,1)
        Info = -4
        RETURN
      ENDIF
      !
      CALL SPLPMN(USRMAT,Mrelas,Nvars,Costs,Prgopt,Dattrv,Bl,Bu,Ind,Info,&
        Primal,Duals,Work(lamat),Work(lcsc),Work(lcolnr),&
        Work(lerd),Work(lerp),Work(lbasma),Work(lwr),Work(lrz),&
        Work(lrg),Work(lrprim),Work(lrhs),Work(lww),lmx,lbm,&
        Ibasis,Iwork(libb),Iwork(limat),Iwork(librc),Iwork(lipr),&
        Iwork(liwr))
      RETURN
    ELSEIF ( ictopt<=ictmax ) THEN
      key = INT( Prgopt(last+1) )
      !
      !     IF KEY = 53, USER MAY SPECIFY LENGTHS OF PORTIONS
      !    OF WORK(*) AND IWORK(*) THAT ARE ALLOCATED TO THE
      !     SPARSE MATRIX STORAGE AND SPARSE LINEAR EQUATION
      !     SOLVING.
      IF ( key==53 ) THEN
        IF ( Prgopt(last+2)/=zero ) THEN
          lmx = INT( Prgopt(last+3) )
          lbm = INT( Prgopt(last+4) )
        ENDIF
      ENDIF
      ictopt = ictopt + 1
      last = next
    ELSE
      nerr = 15
      CALL XERMSG('SLATEC','SPLP','OPTION ARRAY PROCESSING IS CYCLING.',&
        nerr,iopt)
      Info = -nerr
      RETURN
    ENDIF
  ENDDO
  nerr = 21
  CALL XERMSG('SLATEC','SPLP','USER-DEFINED VALUE OF LBM MUST BE .GE. 0.',&
    nerr,iopt)
  Info = -nerr
  RETURN
  !
  !     CALL SPLPMN(USRMAT,MRELAS,NVARS,COSTS,PRGOPT,DATTRV,
  !    1 BL,BU,IND,INFO,PRIMAL,DUALS,AMAT,
  !    2 CSC,COLNRM,ERD,ERP,BASMAT,
  !    3 WR,RZ,RG,RPRIM,RHS,
  !    4 WW,LMX,LBM,IBASIS,IBB,IMAT,
  !    5 IBRC,IPR,IWR)
  !
  RETURN
END SUBROUTINE SPLP
