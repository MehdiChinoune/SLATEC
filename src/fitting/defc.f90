!** DEFC
SUBROUTINE DEFC(Ndata,Xdata,Ydata,Sddata,Nord,Nbkpt,Bkpt,Mdein,Mdeout,&
    Coeff,Lw,W)
  IMPLICIT NONE
  !>
  !***
  !  Fit a piecewise polynomial curve to discrete data.
  !            The piecewise polynomials are represented as B-splines.
  !            The fitting is done in a weighted least squares sense.
  !***
  ! **Library:**   SLATEC
  !***
  ! **Category:**  K1A1A1, K1A2A, L8A3
  !***
  ! **Type:**      DOUBLE PRECISION (EFC-S, DEFC-D)
  !***
  ! **Keywords:**  B-SPLINE, CONSTRAINED LEAST SQUARES, CURVE FITTING
  !***
  ! **Author:**  Hanson, R. J., (SNLA)
  !***
  ! **Description:**
  !
  !      This subprogram fits a piecewise polynomial curve
  !      to discrete data.  The piecewise polynomials are
  !      represented as B-splines.
  !      The fitting is done in a weighted least squares sense.
  !
  !      The data can be processed in groups of modest size.
  !      The size of the group is chosen by the user.  This feature
  !      may be necessary for purposes of using constrained curve fitting
  !      with subprogram DFC( ) on a very large data set.
  !
  !      For a description of the B-splines and usage instructions to
  !      evaluate them, see
  !
  !      C. W. de Boor, Package for Calculating with B-Splines.
  !                     SIAM J. Numer. Anal., p. 441, (June, 1977).
  !
  !      For further discussion of (constrained) curve fitting using
  !      B-splines, see
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
  !                         Internal to DEFC( ) the extreme end knots may
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
  !        MDEIN
  !                         An integer flag, with one of two possible
  !                         values (1 or 2), that directs the subprogram
  !                         action with regard to new data points provided
  !                         by the user.
  !
  !                         =1  The first time that DEFC( ) has been
  !                         entered.  There are NDATA points to process.
  !
  !                         =2  This is another entry to DEFC().  The sub-
  !                         program DEFC( ) has been entered with MDEIN=1
  !                         exactly once before for this problem.  There
  !                         are NDATA new additional points to merge and
  !                         process with any previous points.
  !                         (When using DEFC( ) with MDEIN=2 it is import-
  !                         ant that the set of knots remain fixed at the
  !                         same values for all entries to DEFC( ).)
  !       LW
  !                         The amount of working storage actually
  !                         allocated for the working array W(*).
  !                         This quantity is compared with the
  !                         actual amount of storage needed in DEFC( ).
  !                         Insufficient storage allocated for W(*) is
  !                         an error.  This feature was included in DEFC
  !                         because misreading the storage formula
  !                         for W(*) might very well lead to subtle
  !                         and hard-to-find programming bugs.
  !
  !                         The length of the array W(*) must satisfy
  !
  !                         LW .GE. (NBKPT-NORD+3)*(NORD+1)+
  !                                 (NBKPT+1)*(NORD+1)+
  !                               2*MAX(NDATA,NBKPT)+NBKPT+NORD**2
  !
  !  Output.. All TYPE REAL variables are DOUBLE PRECISION
  !      MDEOUT
  !                         An output flag that indicates the status
  !                         of the curve fit.
  !
  !                         =-1  A usage error of DEFC( ) occurred.  The
  !                         offending condition is noted with the SLATEC
  !                         library error processor, XERMSG( ).  In case
  !                         the working array W(*) is not long enough, the
  !                         minimal acceptable length is printed.
  !
  !                         =1  The B-spline coefficients for the fitted
  !                         curve have been returned in array COEFF(*).
  !
  !                         =2  Not enough data has been processed to
  !                         determine the B-spline coefficients.
  !                         The user has one of two options.  Continue
  !                         to process more data until a unique set
  !                         of coefficients is obtained, or use the
  !                         subprogram DFC( ) to obtain a specific
  !                         set of coefficients.  The user should read
  !                         the usage instructions for DFC( ) for further
  !                         details if this second option is chosen.
  !      COEFF(*)
  !                         If the output value of MDEOUT=1, this array
  !                         contains the unknowns obtained from the least
  !                         squares fitting process.  These N=NBKPT-NORD
  !                         parameters are the B-spline coefficients.
  !                         For MDEOUT=2, not enough data was processed to
  !                         uniquely determine the B-spline coefficients.
  !                         In this case, and also when MDEOUT=-1, all
  !                         values of COEFF(*) are set to zero.
  !
  !                         If the user is not satisfied with the fitted
  !                         curve returned by DEFC( ), the constrained
  !                         least squares curve fitting subprogram DFC( )
  !                         may be required.  The work done within DEFC( )
  !                         to accumulate the data can be utilized by
  !                         the user, if so desired.  This involves
  !                         saving the first (NBKPT-NORD+3)*(NORD+1)
  !                         entries of W(*) and providing this data
  !                         to DFC( ) with the "old problem" designation.
  !                         The user should read the usage instructions
  !                         for subprogram DFC( ) for further details.
  !
  !  Working Array.. All TYPE REAL variables are DOUBLE PRECISION
  !      W(*)
  !                         This array is typed DOUBLE PRECISION.
  !                         Its length is  specified as an input parameter
  !                         in LW as noted above.  The contents of W(*)
  !                         must not be modified by the user between calls
  !                         to DEFC( ) with values of MDEIN=1,2,2,... .
  !                         The first (NBKPT-NORD+3)*(NORD+1) entries of
  !                         W(*) are acceptable as direct input to DFC( )
  !                         for an "old problem" only when MDEOUT=1 or 2.
  !
  !  Evaluating the
  !  Fitted Curve..
  !                         To evaluate derivative number IDER at XVAL,
  !                         use the function subprogram DBVALU( ).
  !
  !                         F = DBVALU(BKPT,COEFF,NBKPT-NORD,NORD,IDER,
  !                                      XVAL,INBV,WORKB)
  !
  !                         The output of this subprogram will not be
  !                         defined unless an output value of MDEOUT=1
  !                         was obtained from DEFC( ), XVAL is in the data
  !                         interval, and IDER is nonnegative and .LT.
  !                         NORD.
  !
  !                         The first time DBVALU( ) is called, INBV=1
  !                         must be specified.  This value of INBV is the
  !                         overwritten by DBVALU( ).  The array WORKB(*)
  !                         must be of length at least 3*NORD, and must
  !                         not be the same as the W(*) array used in the
  !                         call to DEFC( ).
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
  ! **Routines called:**  DEFCMN

  !* REVISION HISTORY  (YYMMDD)
  !   800801  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900510  Change Prologue comments to refer to XERMSG.  (RWC)
  !   900607  Editorial changes to Prologue to make Prologues for EFC,
  !           DEFC, FC, and DFC look as much the same as possible.  (RWC)
  !   920501  Reformatted the REFERENCES section.  (WRB)

  !
  !      SUBROUTINE           FUNCTION/REMARKS
  !
  !      DBSPVN( )          COMPUTE FUNCTION VALUES OF B-SPLINES.  FROM
  !                         THE B-SPLINE PACKAGE OF DE BOOR NOTED ABOVE.
  !
  !      DBNDAC( ),         BANDED LEAST SQUARES MATRIX PROCESSORS.
  !      DBNDSL( )          FROM LAWSON-HANSON, SOLVING LEAST
  !                         SQUARES PROBLEMS.
  !
  !      DSORT( )           DATA SORTING SUBROUTINE, FROM THE
  !                         SANDIA MATH. LIBRARY, SAND77-1441.
  !
  !      XERMSG( )          ERROR HANDLING ROUTINE
  !                         FOR THE SLATEC MATH. LIBRARY.
  !                         SEE SAND78-1189, BY R. E. JONES.
  !
  !      DCOPY( ),DSCAL( )  SUBROUTINES FROM THE BLAS PACKAGE.
  !
  !                         WRITTEN BY R. HANSON, SANDIA NATL. LABS.,
  !                         ALB., N. M., AUGUST-SEPTEMBER, 1980.
  !
  REAL(8) :: Bkpt(*), Coeff(*), W(*), Sddata(*), Xdata(*), Ydata(*)
  INTEGER Lw, Mdein, Mdeout, Nbkpt, Ndata, Nord
  !
  EXTERNAL :: DEFCMN
  !
  INTEGER lbf, lbkpt, lg, lptemp, lww, lxtemp, mdg, mdw
  !
  !* FIRST EXECUTABLE STATEMENT  DEFC
  !     LWW=1               USAGE IN DEFCMN( ) OF W(*)..
  !     LWW,...,LG-1        W(*,*)
  !
  !     LG,...,LXTEMP-1     G(*,*)
  !
  !     LXTEMP,...,LPTEMP-1 XTEMP(*)
  !
  !     LPTEMP,...,LBKPT-1  PTEMP(*)
  !
  !     LBKPT,...,LBF       BKPT(*) (LOCAL TO DEFCMN( ))
  !
  !     LBF,...,LBF+NORD**2 BF(*,*)
  !
  mdg = Nbkpt + 1
  mdw = Nbkpt - Nord + 3
  lww = 1
  lg = lww + mdw*(Nord+1)
  lxtemp = lg + mdg*(Nord+1)
  lptemp = lxtemp + MAX(Ndata,Nbkpt)
  lbkpt = lptemp + MAX(Ndata,Nbkpt)
  lbf = lbkpt + Nbkpt
  CALL DEFCMN(Ndata,Xdata,Ydata,Sddata,Nord,Nbkpt,Bkpt,Mdein,Mdeout,Coeff,&
    W(lbf),W(lxtemp),W(lptemp),W(lbkpt),W(lg),mdg,W(lww),mdw,Lw)
END SUBROUTINE DEFC
