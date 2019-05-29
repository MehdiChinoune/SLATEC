!** PCHBS
SUBROUTINE PCHBS(N,X,F,D,Incfd,Knotyp,Nknots,T,Bcoef,Ndim,Kord,Ierr)
  !>
  !  Piecewise Cubic Hermite to B-Spline converter.
  !***
  ! **Library:**   SLATEC (PCHIP)
  !***
  ! **Category:**  E3
  !***
  ! **Type:**      SINGLE PRECISION (PCHBS-S, DPCHBS-D)
  !***
  ! **Keywords:**  B-SPLINES, CONVERSION, CUBIC HERMITE INTERPOLATION,
  !             PIECEWISE CUBIC INTERPOLATION
  !***
  ! **Author:**  Fritsch, F. N., (LLNL)
  !             Computing and Mathematics Research Division
  !             Lawrence Livermore National Laboratory
  !             P.O. Box 808  (L-316)
  !             Livermore, CA  94550
  !             FTS 532-4275, (510) 422-4275
  !***
  ! **Description:**
  !
  !- Usage:
  !
  !        INTEGER  N, INCFD, KNOTYP, NKNOTS, NDIM, KORD, IERR
  !        PARAMETER  (INCFD = ...)
  !        REAL  X(nmax), F(INCFD,nmax), D(INCFD,nmax), T(2*nmax+4),
  !       *      BCOEF(2*nmax)
  !
  !        CALL  PCHBS (N, X, F, D, INCFD, KNOTYP, NKNOTS, T, BCOEF,
  !       *             NDIM, KORD, IERR)
  !
  !- Arguments:
  !
  !     N:IN  is the number of data points, N.ge.2 .  (not checked)
  !
  !     X:IN  is the real array of independent variable values.  The
  !           elements of X must be strictly increasing:
  !                X(I-1) .LT. X(I),  I = 2(1)N.   (not checked)
  !           nmax, the dimension of X, must be .ge.N.
  !
  !     F:IN  is the real array of dependent variable values.
  !           F(1+(I-1)*INCFD) is the value corresponding to X(I).
  !           nmax, the second dimension of F, must be .ge.N.
  !
  !     D:IN  is the real array of derivative values at the data points.
  !           D(1+(I-1)*INCFD) is the value corresponding to X(I).
  !           nmax, the second dimension of D, must be .ge.N.
  !
  !     INCFD:IN  is the increment between successive values in F and D.
  !           This argument is provided primarily for 2-D applications.
  !           It may have the value 1 for one-dimensional applications,
  !           in which case F and D may be singly-subscripted arrays.
  !
  !     KNOTYP:IN  is a flag to control the knot sequence.
  !           The knot sequence T is normally computed from X by putting
  !           a double knot at each X and setting the end knot pairs
  !           according to the value of KNOTYP:
  !              KNOTYP = 0:  Quadruple knots at X(1) and X(N).  (default)
  !              KNOTYP = 1:  Replicate lengths of extreme subintervals:
  !                           T( 1 ) = T( 2 ) = X(1) - (X(2)-X(1))  ;
  !                           T(M+4) = T(M+3) = X(N) + (X(N)-X(N-1)).
  !              KNOTYP = 2:  Periodic placement of boundary knots:
  !                           T( 1 ) = T( 2 ) = X(1) - (X(N)-X(N-1));
  !                           T(M+4) = T(M+3) = X(N) + (X(2)-X(1))  .
  !              Here M=NDIM=2*N.
  !           If the input value of KNOTYP is negative, however, it is
  !           assumed that NKNOTS and T were set in a previous call.
  !           This option is provided for improved efficiency when used
  !           in a parametric setting.
  !
  !     NKNOTS:INOUT  is the number of knots.
  !           If KNOTYP.GE.0, then NKNOTS will be set to NDIM+4.
  !           If KNOTYP.LT.0, then NKNOTS is an input variable, and an
  !              error return will be taken if it is not equal to NDIM+4.
  !
  !     T:INOUT  is the array of 2*N+4 knots for the B-representation.
  !           If KNOTYP.GE.0, T will be returned by PCHBS with the
  !              interior double knots equal to the X-values and the
  !              boundary knots set as indicated above.
  !           If KNOTYP.LT.0, it is assumed that T was set by a
  !              previous call to PCHBS.  (This routine does **not**
  !              verify that T forms a legitimate knot sequence.)
  !
  !     BCOEF:OUT  is the array of 2*N B-spline coefficients.
  !
  !     NDIM:OUT  is the dimension of the B-spline space.  (Set to 2*N.)
  !
  !     KORD:OUT  is the order of the B-spline.  (Set to 4.)
  !
  !     IERR:OUT  is an error flag.
  !           Normal return:
  !              IERR = 0  (no errors).
  !           "Recoverable" errors:
  !              IERR = -4  if KNOTYP.GT.2 .
  !              IERR = -5  if KNOTYP.LT.0 and NKNOTS.NE.(2*N+4).
  !
  !- Description:
  !     PCHBS computes the B-spline representation of the PCH function
  !     determined by N,X,F,D.  To be compatible with the rest of PCHIP,
  !     PCHBS includes INCFD, the increment between successive values of
  !     the F- and D-arrays.
  !
  !     The output is the B-representation for the function:  NKNOTS, T,
  !     BCOEF, NDIM, KORD.
  !
  !- Caution:
  !     Since it is assumed that the input PCH function has been
  !     computed by one of the other routines in the package PCHIP,
  !     input arguments N, X, INCFD are **not** checked for validity.
  !
  !- Restrictions/assumptions:
  !     1. N.GE.2 .  (not checked)
  !     2. X(i).LT.X(i+1), i=1,...,N .  (not checked)
  !     3. INCFD.GT.0 .  (not checked)
  !     4. KNOTYP.LE.2 .  (error return if not)
  !    *5. NKNOTS = NDIM+4 = 2*N+4 .  (error return if not)
  !    *6. T(2*k+1) = T(2*k) = X(k), k=1,...,N .  (not checked)
  !
  !       * Indicates this applies only if KNOTYP.LT.0 .
  !
  !- Portability:
  !     Argument INCFD is used only to cause the compiler to generate
  !     efficient code for the subscript expressions (1+(I-1)*INCFD) .
  !     The normal usage, in which PCHBS is called with one-dimensional
  !     arrays F and D, is probably non-Fortran 77, in the strict sense,
  !     but it works on all systems on which PCHBS has been tested.
  !
  !- See Also:
  !     PCHIC, PCHIM, or PCHSP can be used to determine an interpolating
  !        PCH function from a set of data.
  !     The B-spline routine BVALU can be used to evaluate the
  !        B-representation that is output by PCHBS.
  !        (See BSPDOC for more information.)
  !
  !***
  ! **References:**  F. N. Fritsch, "Representations for parametric cubic
  !                 splines," Computer Aided Geometric Design 6 (1989),
  !                 pp.79-82.
  !***
  ! **Routines called:**  PCHKT, XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   870701  DATE WRITTEN
  !   900405  Converted Fortran to upper case.
  !   900405  Removed requirement that X be dimensioned N+1.
  !   900406  Modified to make PCHKT a subsidiary routine to simplify
  !           usage.  In the process, added argument INCFD to be com-
  !           patible with the rest of PCHIP.
  !   900410  Converted prologue to SLATEC 4.0 format.
  !   900410  Added calls to XERMSG and changed constant 3. to 3 to
  !           reduce single/double differences.
  !   900411  Added reference.
  !   900501  Corrected declarations.
  !   930317  Minor cosmetic changes.  (FNF)
  !   930514  Corrected problems with dimensioning of arguments and
  !           clarified DESCRIPTION.  (FNF)
  !   930604  Removed  NKNOTS from PCHKT call list.  (FNF)
  USE service, ONLY : XERMSG
  !
  !*Internal Notes:
  !
  !**End
  !
  !  Declare arguments.
  !
  INTEGER N, Incfd, Knotyp, Nknots, Ndim, Kord, Ierr
  REAL X(N), F(Incfd,N), D(Incfd,N), T(2*N+4), Bcoef(2*N)
  !
  !  Declare local variables.
  !
  INTEGER k, kk
  REAL dov3, hnew, hold
  CHARACTER(8) :: libnam, subnam
  !* FIRST EXECUTABLE STATEMENT  PCHBS
  !
  !  Initialize.
  !
  Ndim = 2*N
  Kord = 4
  Ierr = 0
  libnam = 'SLATEC'
  subnam = 'PCHBS'
  !
  !  Check argument validity.  Set up knot sequence if OK.
  !
  IF ( Knotyp>2 ) THEN
    Ierr = -1
    CALL XERMSG(libnam,subnam,'KNOTYP GREATER THAN 2',Ierr,1)
    RETURN
  END IF
  IF ( Knotyp>=0 ) THEN
    !          Set up knot sequence.
    Nknots = Ndim + 4
    CALL PCHKT(N,X,Knotyp,T)
  ELSEIF ( Nknots/=Ndim+4 ) THEN
    Ierr = -2
    CALL XERMSG(libnam,subnam,'KNOTYP.LT.0 AND NKNOTS.NE.(2*N+4)',Ierr,1)
    RETURN
  END IF
  !
  !  Compute B-spline coefficients.
  !
  hnew = T(3) - T(1)
  DO k = 1, N
    kk = 2*k
    hold = hnew
    !          The following requires mixed mode arithmetic.
    dov3 = D(1,k)/3
    Bcoef(kk-1) = F(1,k) - hold*dov3
    !          The following assumes T(2*K+1) = X(K).
    hnew = T(kk+3) - T(kk+1)
    Bcoef(kk) = F(1,k) + hnew*dov3
  END DO
  !
  !  Terminate.
  !
  !------------- LAST LINE OF PCHBS FOLLOWS ------------------------------
END SUBROUTINE PCHBS
