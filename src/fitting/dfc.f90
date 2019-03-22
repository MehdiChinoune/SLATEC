!** DFC
SUBROUTINE DFC(Ndata,Xdata,Ydata,Sddata,Nord,Nbkpt,Bkpt,Nconst,Xconst,&
    Yconst,Nderiv,Mode,Coeff,W,Iw)
  IMPLICIT NONE
  !>
  !***
  !  Fit a piecewise polynomial curve to discrete data.
  !            The piecewise polynomials are represented as B-splines.
  !            The fitting is done in a weighted least squares sense.
  !            Equality and inequality constraints can be imposed on the
  !            fitted curve.
  !***
  ! **Library:**   SLATEC
  !***
  ! **Category:**  K1A1A1, K1A2A, L8A3
  !***
  ! **Type:**      DOUBLE PRECISION (FC-S, DFC-D)
  !***
  ! **Keywords:**  B-SPLINE, CONSTRAINED LEAST SQUARES, CURVE FITTING,
  !             WEIGHTED LEAST SQUARES
  !***
  ! **Author:**  Hanson, R. J., (SNLA)
  !***
  ! **Description:**
  !
  !      This subprogram fits a piecewise polynomial curve
  !      to discrete data.  The piecewise polynomials are
  !      represented as B-splines.
  !      The fitting is done in a weighted least squares sense.
  !      Equality and inequality constraints can be imposed on the
  !      fitted curve.
  !
  !      For a description of the B-splines and usage instructions to
  !      evaluate them, see
  !
  !      C. W. de Boor, Package for Calculating with B-Splines.
  !                     SIAM J. Numer. Anal., p. 441, (June, 1977).
  !
  !      For further documentation and discussion of constrained
  !      curve fitting using B-splines, see
  !
  !      R. J. Hanson, Constrained Least Squares Curve Fitting
  !                   to Discrete Data Using B-Splines, a User's
  !                   Guide. Sandia Labs. Tech. Rept. SAND-78-1291,
  !                   December, (1978).
  !
  !  Input.. All TYPE REAL variables are DOUBLE PRECISION
  !      NDATA,XDATA(*),
  !      YDATA(*),
  !      SDDATA(*)
  !                         The NDATA discrete (X,Y) pairs and the Y value
  !                         standard deviation or uncertainty, SD, are in
  !                         the respective arrays XDATA(*), YDATA(*), and
  !                         SDDATA(*).  No sorting of XDATA(*) is
  !                         required.  Any non-negative value of NDATA is
  !                         allowed.  A negative value of NDATA is an
  !                         error.  A zero value for any entry of
  !                         SDDATA(*) will weight that data point as 1.
  !                         Otherwise the weight of that data point is
  !                         the reciprocal of this entry.
  !
  !      NORD,NBKPT,
  !      BKPT(*)
  !                         The NBKPT knots of the B-spline of order NORD
  !                         are in the array BKPT(*).  Normally the
  !                         problem data interval will be included between
  !                         the limits BKPT(NORD) and BKPT(NBKPT-NORD+1).
  !                         The additional end knots BKPT(I),I=1,...,
  !                         NORD-1 and I=NBKPT-NORD+2,...,NBKPT, are
  !                         required to compute the functions used to fit
  !                         the data.  No sorting of BKPT(*) is required.
  !                         Internal to DFC( ) the extreme end knots may
  !                         be reduced and increased respectively to
  !                         accommodate any data values that are exterior
  !                         to the given knot values.  The contents of
  !                         BKPT(*) is not changed.
  !
  !                         NORD must be in the range 1 .LE. NORD .LE. 20.
  !                         The value of NBKPT must satisfy the condition
  !                         NBKPT .GE. 2*NORD.
  !                         Other values are considered errors.
  !
  !                         (The order of the spline is one more than the
  !                         degree of the piecewise polynomial defined on
  !                         each interval.  This is consistent with the
  !                         B-spline package convention.  For example,
  !                         NORD=4 when we are using piecewise cubics.)
  !
  !      NCONST,XCONST(*),
  !      YCONST(*),NDERIV(*)
  !                         The number of conditions that constrain the
  !                         B-spline is NCONST.  A constraint is specified
  !                         by an (X,Y) pair in the arrays XCONST(*) and
  !                         YCONST(*), and by the type of constraint and
  !                         derivative value encoded in the array
  !                         NDERIV(*).  No sorting of XCONST(*) is
  !                         required.  The value of NDERIV(*) is
  !                         determined as follows.  Suppose the I-th
  !                         constraint applies to the J-th derivative
  !                         of the B-spline.  (Any non-negative value of
  !                         J < NORD is permitted.  In particular the
  !                         value J=0 refers to the B-spline itself.)
  !                         For this I-th constraint, set
  !                          XCONST(I)=X,
  !                          YCONST(I)=Y, and
  !                          NDERIV(I)=ITYPE+4*J, where
  !
  !                          ITYPE = 0,      if (J-th deriv. at X) .LE. Y.
  !                                = 1,      if (J-th deriv. at X) .GE. Y.
  !                                = 2,      if (J-th deriv. at X) .EQ. Y.
  !                                = 3,      if (J-th deriv. at X) .EQ.
  !                                             (J-th deriv. at Y).
  !                          (A value of NDERIV(I)=-1 will cause this
  !                          constraint to be ignored.  This subprogram
  !                          feature is often useful when temporarily
  !                          suppressing a constraint while still
  !                          retaining the source code of the calling
  !                          program.)
  !
  !        MODE
  !                         An input flag that directs the least squares
  !                         solution method used by DFC( ).
  !
  !                         The variance function, referred to below,
  !                         defines the square of the probable error of
  !                         the fitted curve at any point, XVAL.
  !                         This feature of DFC( ) allows one to use the
  !                         square root of this variance function to
  !                         determine a probable error band around the
  !                         fitted curve.
  !
  !                         =1  a new problem.  No variance function.
  !
  !                         =2  a new problem.  Want variance function.
  !
  !                         =3  an old problem.  No variance function.
  !
  !                         =4  an old problem.  Want variance function.
  !
  !                         Any value of MODE other than 1-4 is an error.
  !
  !                         The user with a new problem can skip directly
  !                         to the description of the input parameters
  !                         IW(1), IW(2).
  !
  !                         If the user correctly specifies the new or old
  !                         problem status, the subprogram DFC( ) will
  !                         perform more efficiently.
  !                         By an old problem it is meant that subprogram
  !                         DFC( ) was last called with this same set of
  !                         knots, data points and weights.
  !
  !                         Another often useful deployment of this old
  !                         problem designation can occur when one has
  !                         previously obtained a Q-R orthogonal
  !                         decomposition of the matrix resulting from
  !                         B-spline fitting of data (without constraints)
  !                         at the breakpoints BKPT(I), I=1,...,NBKPT.
  !                         For example, this matrix could be the result
  !                         of sequential accumulation of the least
  !                         squares equations for a very large data set.
  !                         The user writes this code in a manner
  !                         convenient for the application.  For the
  !                         discussion here let
  !
  !                                      N=NBKPT-NORD, and K=N+3
  !
  !                         Let us assume that an equivalent least squares
  !                         system
  !
  !                                      RC=D
  !
  !                         has been obtained.  Here R is an N+1 by N
  !                         matrix and D is a vector with N+1 components.
  !                         The last row of R is zero.  The matrix R is
  !                         upper triangular and banded.  At most NORD of
  !                         the diagonals are nonzero.
  !                         The contents of R and D can be copied to the
  !                         working array W(*) as follows.
  !
  !                         The I-th diagonal of R, which has N-I+1
  !                         elements, is copied to W(*) starting at
  !
  !                                      W((I-1)*K+1),
  !
  !                         for I=1,...,NORD.
  !                         The vector D is copied to W(*) starting at
  !
  !                                      W(NORD*K+1)
  !
  !                         The input value used for NDATA is arbitrary
  !                         when an old problem is designated.  Because
  !                         of the feature of DFC( ) that checks the
  !                         working storage array lengths, a value not
  !                         exceeding NBKPT should be used.  For example,
  !                         use NDATA=0.
  !
  !                         (The constraints or variance function request
  !                         can change in each call to DFC( ).)  A new
  !                         problem is anything other than an old problem.
  !
  !      IW(1),IW(2)
  !                         The amounts of working storage actually
  !                         allocated for the working arrays W(*) and
  !                         IW(*).  These quantities are compared with the
  !                         actual amounts of storage needed in DFC( ).
  !                         Insufficient storage allocated for either
  !                         W(*) or IW(*) is an error.  This feature was
  !                         included in DFC( ) because misreading the
  !                         storage formulas for W(*) and IW(*) might very
  !                         well lead to subtle and hard-to-find
  !                         programming bugs.
  !
  !                         The length of W(*) must be at least
  !
  !                           NB=(NBKPT-NORD+3)*(NORD+1)+
  !                               2*MAX(NDATA,NBKPT)+NBKPT+NORD**2
  !
  !                         Whenever possible the code uses banded matrix
  !                         processors DBNDAC( ) and DBNDSL( ).  These
  !                         are utilized if there are no constraints,
  !                         no variance function is required, and there
  !                         is sufficient data to uniquely determine the
  !                         B-spline coefficients.  If the band processors
  !                         cannot be used to determine the solution,
  !                         then the constrained least squares code DLSEI
  !                         is used.  In this case the subprogram requires
  !                         an additional block of storage in W(*).  For
  !                         the discussion here define the integers NEQCON
  !                         and NINCON respectively as the number of
  !                         equality (ITYPE=2,3) and inequality
  !                         (ITYPE=0,1) constraints imposed on the fitted
  !                         curve.  Define
  !
  !                           L=NBKPT-NORD+1
  !
  !                         and note that
  !
  !                           NCONST=NEQCON+NINCON.
  !
  !                         When the subprogram DFC( ) uses DLSEI( ) the
  !                         length of the working array W(*) must be at
  !                         least
  !
  !                           LW=NB+(L+NCONST)*L+
  !                              2*(NEQCON+L)+(NINCON+L)+(NINCON+2)*(L+6)
  !
  !                         The length of the array IW(*) must be at least
  !
  !                           IW1=NINCON+2*L
  !
  !                         in any case.
  !
  !  Output.. All TYPE REAL variables are DOUBLE PRECISION
  !      MODE
  !                         An output flag that indicates the status
  !                         of the constrained curve fit.
  !
  !                         =-1  a usage error of DFC( ) occurred.  The
  !                              offending condition is noted with the
  !                              SLATEC library error processor, XERMSG.
  !                              In case the working arrays W(*) or IW(*)
  !                              are not long enough, the minimal
  !                              acceptable length is printed.
  !
  !                         = 0  successful constrained curve fit.
  !
  !                         = 1  the requested equality constraints
  !                              are contradictory.
  !
  !                         = 2  the requested inequality constraints
  !                              are contradictory.
  !
  !                         = 3  both equality and inequality constraints
  !                              are contradictory.
  !
  !      COEFF(*)
  !                         If the output value of MODE=0 or 1, this array
  !                         contains the unknowns obtained from the least
  !                         squares fitting process.  These N=NBKPT-NORD
  !                         parameters are the B-spline coefficients.
  !                         For MODE=1, the equality constraints are
  !                         contradictory.  To make the fitting process
  !                         more robust, the equality constraints are
  !                         satisfied in a least squares sense.  In this
  !                         case the array COEFF(*) contains B-spline
  !                         coefficients for this extended concept of a
  !                         solution.  If MODE=-1,2 or 3 on output, the
  !                         array COEFF(*) is undefined.
  !
  !  Working Arrays.. All Type REAL variables are DOUBLE PRECISION
  !      W(*),IW(*)
  !                         These arrays are respectively typed DOUBLE
  !                         PRECISION and INTEGER.
  !                         Their required lengths are specified as input
  !                         parameters in IW(1), IW(2) noted above.  The
  !                         contents of W(*) must not be modified by the
  !                         user if the variance function is desired.
  !
  !  Evaluating the
  !  Variance Function..
  !                         To evaluate the variance function (assuming
  !                         that the uncertainties of the Y values were
  !                         provided to DFC( ) and an input value of
  !                         MODE=2 or 4 was used), use the function
  !                         subprogram DCV( )
  !
  !                           VAR=DCV(XVAL,NDATA,NCONST,NORD,NBKPT,
  !                                  BKPT,W)
  !
  !                         Here XVAL is the point where the variance is
  !                         desired.  The other arguments have the same
  !                         meaning as in the usage of DFC( ).
  !
  !                         For those users employing the old problem
  !                         designation, let MDATA be the number of data
  !                         points in the problem.  (This may be different
  !                         from NDATA if the old problem designation
  !                         feature was used.)  The value, VAR, should be
  !                         multiplied by the quantity
  !
  !                         REAL(MAX(NDATA-N,1), 8)/REAL(MAX(MDATA-N,1), 8)
  !
  !                         The output of this subprogram is not defined
  !                         if an input value of MODE=1 or 3 was used in
  !                         FC( ) or if an output value of MODE=-1, 2, or
  !                         3 was obtained.  The variance function, except
  !                         for the scaling factor noted above, is given
  !                         by
  !
  !                           VAR=(transpose of B(XVAL))*C*B(XVAL)
  !
  !                         The vector B(XVAL) is the B-spline basis
  !                         function values at X=XVAL.
  !                         The covariance matrix, C, of the solution
  !                         coefficients accounts only for the least
  !                         squares equations and the explicitly stated
  !                         equality constraints.  This fact must be
  !                         considered when interpreting the variance
  !                         function from a data fitting problem that has
  !                         inequality constraints on the fitted curve.
  !
  !  Evaluating the
  !  Fitted Curve..
  !                         To evaluate derivative number IDER at XVAL,
  !                         use the function subprogram DBVALU( )
  !
  !                           F = DBVALU(BKPT,COEFF,NBKPT-NORD,NORD,IDER,
  !                                      XVAL,INBV,WORKB)
  !
  !                         The output of this subprogram will not be
  !                         defined unless an output value of MODE=0 or 1
  !                         was obtained from DFC( ), XVAL is in the data
  !                         interval, and IDER is nonnegative and .LT.
  !                         NORD.
  !
  !                         The first time DBVALU( ) is called, INBV=1
  !                         must be specified.  This value of INBV is the
  !                         overwritten by DBVALU( ).  The array WORKB(*)
  !                         must be of length at least 3*NORD, and must
  !                         not be the same as the W(*) array used in
  !                         the call to DFC( ).
  !
  !                         DBVALU( ) expects the breakpoint array BKPT(*)
  !                         to be sorted.
  !
  !***
  ! **References:**  R. J. Hanson, Constrained least squares curve fitting
  !                 to discrete data using B-splines, a users guide,
  !                 Report SAND78-1291, Sandia Laboratories, December
  !                 1978.
  !***
  ! **Routines called:**  DFCMN

  !* REVISION HISTORY  (YYMMDD)
  !   780801  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   891006  Cosmetic changes to prologue.  (WRB)
  !   891006  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900510  Convert references to XERRWV to references to XERMSG.  (RWC)
  !   900607  Editorial changes to Prologue to make Prologues for EFC,
  !           DEFC, FC, and DFC look as much the same as possible.  (RWC)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  
  REAL(8) :: Bkpt(*), Coeff(*), Sddata(*), W(*), Xconst(*), &
    Xdata(*), Yconst(*), Ydata(*)
  INTEGER Iw(*), Mode, Nbkpt, Nconst, Ndata, Nderiv(*), Nord
  !
  EXTERNAL :: DFCMN
  !
  INTEGER i1, i2, i3, i4, i5, i6, i7, mdg, mdw
  !
  !* FIRST EXECUTABLE STATEMENT  DFC
  IF ( Nord<1.OR.Nord>20 ) THEN
    CALL XERMSG('SLATEC','FCMN',&
      'IN FC, THE ORDER OF THE B-SPLINE MUST BE 1 THRU 20.',2,1)
    Mode = -1
    RETURN
    !
  ELSEIF ( Nbkpt<2*Nord ) THEN
    CALL XERMSG('SLATEC','FCMN',&
      'IN FC, THE NUMBER OF KNOTS MUST BE AT LEAST TWICE THE B-SPLINE ORDER.',2,1)
    Mode = -1
    RETURN
  ENDIF

  mdg = Nbkpt - Nord + 3
  mdw = Nbkpt - Nord + 1 + Nconst
  !                        USAGE IN DFCMN( ) OF W(*)..
  !     I1,...,I2-1      G(*,*)
  !
  !     I2,...,I3-1      XTEMP(*)
  !
  !     I3,...,I4-1      PTEMP(*)
  !
  !     I4,...,I5-1      BKPT(*) (LOCAL TO DFCMN( ))
  !
  !     I5,...,I6-1      BF(*,*)
  !
  !     I6,...,I7-1      W(*,*)
  !
  !     I7,...           WORK(*) FOR DLSEI( )
  !
  i1 = 1
  i2 = i1 + mdg*(Nord+1)
  i3 = i2 + MAX(Ndata,Nbkpt)
  i4 = i3 + MAX(Ndata,Nbkpt)
  i5 = i4 + Nbkpt
  i6 = i5 + Nord*Nord
  i7 = i6 + mdw*(Nbkpt-Nord+1)
  CALL DFCMN(Ndata,Xdata,Ydata,Sddata,Nord,Nbkpt,Bkpt,Nconst,Xconst,Yconst,&
    Nderiv,Mode,Coeff,W(i5),W(i2),W(i3),W(i4),W(i1),mdg,W(i6),mdw,&
    W(i7),Iw)
END SUBROUTINE DFC
