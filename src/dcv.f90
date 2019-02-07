!*==DCV.f90  processed by SPAG 6.72Dc at 11:01 on  6 Feb 2019
!DECK DCV
REAL(8) FUNCTION DCV(Xval,Ndata,Nconst,Nord,Nbkpt,Bkpt,W)
  IMPLICIT NONE
  !*--DCV5
  !***BEGIN PROLOGUE  DCV
  !***PURPOSE  Evaluate the variance function of the curve obtained
  !            by the constrained B-spline fitting subprogram DFC.
  !***LIBRARY   SLATEC
  !***CATEGORY  L7A3
  !***TYPE      DOUBLE PRECISION (CV-S, DCV-D)
  !***KEYWORDS  ANALYSIS OF COVARIANCE, B-SPLINE,
  !             CONSTRAINED LEAST SQUARES, CURVE FITTING
  !***AUTHOR  Hanson, R. J., (SNLA)
  !***DESCRIPTION
  !
  !     DCV( ) is a companion function subprogram for DFC( ).  The
  !     documentation for DFC( ) has complete usage instructions.
  !
  !     DCV( ) is used to evaluate the variance function of the curve
  !     obtained by the constrained B-spline fitting subprogram, DFC( ).
  !     The variance function defines the square of the probable error
  !     of the fitted curve at any point, XVAL.  One can use the square
  !     root of this variance function to determine a probable error band
  !     around the fitted curve.
  !
  !     DCV( ) is used after a call to DFC( ).  MODE, an input variable to
  !     DFC( ), is used to indicate if the variance function is desired.
  !     In order to use DCV( ), MODE must equal 2 or 4 on input to DFC( ).
  !     MODE is also used as an output flag from DFC( ).  Check to make
  !     sure that MODE = 0 after calling DFC( ), indicating a successful
  !     constrained curve fit.  The array SDDATA, as input to DFC( ), must
  !     also be defined with the standard deviation or uncertainty of the
  !     Y values to use DCV( ).
  !
  !     To evaluate the variance function after calling DFC( ) as stated
  !     above, use DCV( ) as shown here
  !
  !          VAR=DCV(XVAL,NDATA,NCONST,NORD,NBKPT,BKPT,W)
  !
  !     The variance function is given by
  !
  !      VAR=(transpose of B(XVAL))*C*B(XVAL)/REAL(MAX(NDATA-N,1), 8)
  !
  !     where N = NBKPT - NORD.
  !
  !     The vector B(XVAL) is the B-spline basis function values at
  !     X=XVAL.  The covariance matrix, C, of the solution coefficients
  !     accounts only for the least squares equations and the explicitly
  !     stated equality constraints.  This fact must be considered when
  !     interpreting the variance function from a data fitting problem
  !     that has inequality constraints on the fitted curve.
  !
  !     All the variables in the calling sequence for DCV( ) are used in
  !     DFC( ) except the variable XVAL.  Do not change the values of
  !     these variables between the call to DFC( ) and the use of DCV( ).
  !
  !     The following is a brief description of the variables
  !
  !     XVAL    The point where the variance is desired, a double
  !             precision variable.
  !
  !     NDATA   The number of discrete (X,Y) pairs for which DFC( )
  !             calculated a piece-wise polynomial curve.
  !
  !     NCONST  The number of conditions that constrained the B-spline in
  !             DFC( ).
  !
  !     NORD    The order of the B-spline used in DFC( ).
  !             The value of NORD must satisfy 1 < NORD < 20 .
  !
  !             (The order of the spline is one more than the degree of
  !             the piece-wise polynomial defined on each interval.  This
  !             is consistent with the B-spline package convention.  For
  !             example, NORD=4 when we are using piece-wise cubics.)
  !
  !     NBKPT   The number of knots in the array BKPT(*).
  !             The value of NBKPT must satisfy NBKPT .GE. 2*NORD.
  !
  !     BKPT(*) The double precision array of knots.  Normally the problem
  !             data interval will be included between the limits
  !             BKPT(NORD) and BKPT(NBKPT-NORD+1).  The additional end
  !             knots BKPT(I),I=1,...,NORD-1 and I=NBKPT-NORD+2,...,NBKPT,
  !             are required by DFC( ) to compute the functions used to
  !             fit the data.
  !
  !     W(*)    Double precision work array as used in DFC( ).  See DFC( )
  !             for the required length of W(*).  The contents of W(*)
  !             must not be modified by the user if the variance function
  !             is desired.
  !
  !***REFERENCES  R. J. Hanson, Constrained least squares curve fitting
  !                 to discrete data using B-splines, a users guide,
  !                 Report SAND78-1291, Sandia Laboratories, December
  !                 1978.
  !***ROUTINES CALLED  DDOT, DFSPVN
  !***REVISION HISTORY  (YYMMDD)
  !   780801  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890831  Modified array declarations.  (WRB)
  !   890911  Removed unnecessary intrinsics.  (WRB)
  !   891006  Cosmetic changes to prologue.  (WRB)
  !   891006  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  !***END PROLOGUE  DCV
  INTEGER i , ileft , ip , is , last , mdg , mdw , n , Nbkpt , Nconst , &
    Ndata , Nord
  REAL(8) :: Bkpt , DDOT , v , W , Xval , zero
  DIMENSION Bkpt(*) , W(*) , v(40)
  !***FIRST EXECUTABLE STATEMENT  DCV
  zero = 0.0D0
  mdg = Nbkpt - Nord + 3
  mdw = Nbkpt - Nord + 1 + Nconst
  is = mdg*(Nord+1) + 2*MAX(Ndata,Nbkpt) + Nbkpt + Nord**2
  last = Nbkpt - Nord + 1
  ileft = Nord
  DO WHILE ( Xval>=Bkpt(ileft+1).AND.ileft<last-1 )
    ileft = ileft + 1
  ENDDO
  CALL DFSPVN(Bkpt,Nord,1,Xval,ileft,v(Nord+1))
  ileft = ileft - Nord + 1
  ip = mdw*(ileft-1) + ileft + is
  n = Nbkpt - Nord
  DO i = 1 , Nord
    v(i) = DDOT(Nord,W(ip),1,v(Nord+1),1)
    ip = ip + mdw
  ENDDO
  DCV = MAX(DDOT(Nord,v,1,v(Nord+1),1),zero)
  !
  !     SCALE THE VARIANCE SO IT IS AN UNBIASED ESTIMATE.
  DCV = DCV/MAX(Ndata-n,1)
END FUNCTION DCV
