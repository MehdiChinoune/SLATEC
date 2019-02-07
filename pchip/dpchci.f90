!*==DPCHCI.f90  processed by SPAG 6.72Dc at 11:00 on  6 Feb 2019
!DECK DPCHCI
SUBROUTINE DPCHCI(N,H,Slope,D,Incfd)
  IMPLICIT NONE
  !*--DPCHCI5
  !***BEGIN PROLOGUE  DPCHCI
  !***SUBSIDIARY
  !***PURPOSE  Set interior derivatives for DPCHIC
  !***LIBRARY   SLATEC (PCHIP)
  !***TYPE      DOUBLE PRECISION (PCHCI-S, DPCHCI-D)
  !***AUTHOR  Fritsch, F. N., (LLNL)
  !***DESCRIPTION
  !
  !          DPCHCI:  DPCHIC Initial Derivative Setter.
  !
  !    Called by DPCHIC to set derivatives needed to determine a monotone
  !    piecewise cubic Hermite interpolant to the data.
  !
  !    Default boundary conditions are provided which are compatible
  !    with monotonicity.  If the data are only piecewise monotonic, the
  !    interpolant will have an extremum at each point where monotonicity
  !    switches direction.
  !
  !    To facilitate two-dimensional applications, includes an increment
  !    between successive values of the D-array.
  !
  !    The resulting piecewise cubic Hermite function should be identical
  !    (within roundoff error) to that produced by DPCHIM.
  !
  ! ----------------------------------------------------------------------
  !
  !  Calling sequence:
  !
  !        PARAMETER  (INCFD = ...)
  !        INTEGER  N
  !        DOUBLE PRECISION  H(N), SLOPE(N), D(INCFD,N)
  !
  !        CALL  DPCHCI (N, H, SLOPE, D, INCFD)
  !
  !   Parameters:
  !
  !     N -- (input) number of data points.
  !           If N=2, simply does linear interpolation.
  !
  !     H -- (input) real*8 array of interval lengths.
  !     SLOPE -- (input) real*8 array of data slopes.
  !           If the data are (X(I),Y(I)), I=1(1)N, then these inputs are:
  !                  H(I) =  X(I+1)-X(I),
  !              SLOPE(I) = (Y(I+1)-Y(I))/H(I),  I=1(1)N-1.
  !
  !     D -- (output) real*8 array of derivative values at data points.
  !           If the data are monotonic, these values will determine a
  !           a monotone cubic Hermite function.
  !           The value corresponding to X(I) is stored in
  !                D(1+(I-1)*INCFD),  I=1(1)N.
  !           No other entries in D are changed.
  !
  !     INCFD -- (input) increment between successive values in D.
  !           This argument is provided primarily for 2-D applications.
  !
  !    -------
  !    WARNING:  This routine does no validity-checking of arguments.
  !    -------
  !
  !  Fortran intrinsics used:  ABS, MAX, MIN.
  !
  !***SEE ALSO  DPCHIC
  !***ROUTINES CALLED  DPCHST
  !***REVISION HISTORY  (YYMMDD)
  !   820218  DATE WRITTEN
  !   820601  Modified end conditions to be continuous functions of
  !           data when monotonicity switches in next interval.
  !   820602  1. Modified formulas so end conditions are less prone
  !             to over/underflow problems.
  !           2. Minor modification to HSUM calculation.
  !   820805  Converted to SLATEC library version.
  !   870813  Minor cosmetic changes.
  !   890411  Added SAVE statements (Vers. 3.2).
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890831  Modified array declarations.  (WRB)
  !   890831  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900328  Added TYPE section.  (WRB)
  !   910408  Updated AUTHOR section in prologue.  (WRB)
  !   930503  Improved purpose.  (FNF)
  !***END PROLOGUE  DPCHCI
  !
  !  Programming notes:
  !     1. The function  DPCHST(ARG1,ARG2)  is assumed to return zero if
  !        either argument is zero, +1 if they are of the same sign, and
  !        -1 if they are of opposite sign.
  !**End
  !
  !  DECLARE ARGUMENTS.
  !
  INTEGER N, Incfd
  REAL(8) :: H(*), Slope(*), D(Incfd,*)
  !
  !  DECLARE LOCAL VARIABLES.
  !
  INTEGER i, nless1
  REAL(8) :: del1, del2, dmax, dmin, drat1, drat2, hsum, &
    hsumt3, three, w1, w2, zero
  SAVE zero, three
  REAL(8) :: DPCHST
  !
  !  INITIALIZE.
  !
  DATA zero/0.D0/, three/3.D0/
  !***FIRST EXECUTABLE STATEMENT  DPCHCI
  nless1 = N - 1
  del1 = Slope(1)
  !
  !  SPECIAL CASE N=2 -- USE LINEAR INTERPOLATION.
  !
  IF ( nless1>1 ) THEN
    !
    !  NORMAL CASE  (N .GE. 3).
    !
    del2 = Slope(2)
    !
    !  SET D(1) VIA NON-CENTERED THREE-POINT FORMULA, ADJUSTED TO BE
    !     SHAPE-PRESERVING.
    !
    hsum = H(1) + H(2)
    w1 = (H(1)+hsum)/hsum
    w2 = -H(1)/hsum
    D(1,1) = w1*del1 + w2*del2
    IF ( DPCHST(D(1,1),del1)<=zero ) THEN
      D(1,1) = zero
    ELSEIF ( DPCHST(del1,del2)<zero ) THEN
      !        NEED DO THIS CHECK ONLY IF MONOTONICITY SWITCHES.
      dmax = three*del1
      IF ( ABS(D(1,1))>ABS(dmax) ) D(1,1) = dmax
    ENDIF
    !
    !  LOOP THROUGH INTERIOR POINTS.
    !
    DO i = 2, nless1
      IF ( i/=2 ) THEN
        !
        hsum = H(i-1) + H(i)
        del1 = del2
        del2 = Slope(i)
      ENDIF
      !
      !        SET D(I)=0 UNLESS DATA ARE STRICTLY MONOTONIC.
      !
      D(1,i) = zero
      IF ( DPCHST(del1,del2)>zero ) THEN
        !
        !        USE BRODLIE MODIFICATION OF BUTLAND FORMULA.
        !
        hsumt3 = hsum + hsum + hsum
        w1 = (hsum+H(i-1))/hsumt3
        w2 = (hsum+H(i))/hsumt3
        dmax = MAX(ABS(del1),ABS(del2))
        dmin = MIN(ABS(del1),ABS(del2))
        drat1 = del1/dmax
        drat2 = del2/dmax
        D(1,i) = dmin/(w1*drat1+w2*drat2)
      ENDIF
      !
    ENDDO
    !
    !  SET D(N) VIA NON-CENTERED THREE-POINT FORMULA, ADJUSTED TO BE
    !     SHAPE-PRESERVING.
    !
    w1 = -H(N-1)/hsum
    w2 = (H(N-1)+hsum)/hsum
    D(1,N) = w1*del1 + w2*del2
    IF ( DPCHST(D(1,N),del2)<=zero ) THEN
      D(1,N) = zero
    ELSEIF ( DPCHST(del1,del2)<zero ) THEN
      !        NEED DO THIS CHECK ONLY IF MONOTONICITY SWITCHES.
      dmax = three*del2
      IF ( ABS(D(1,N))>ABS(dmax) ) D(1,N) = dmax
    ENDIF
  ELSE
    D(1,1) = del1
    D(1,N) = del1
  ENDIF
  !
  !  NORMAL RETURN.
  !
  !------------- LAST LINE OF DPCHCI FOLLOWS -----------------------------
END SUBROUTINE DPCHCI
