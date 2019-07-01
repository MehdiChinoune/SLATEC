!** PCHCS
PURE SUBROUTINE PCHCS(Switch,N,H,Slope,D,Incfd,Ierr)
  !> Adjusts derivative values for PCHIC
  !***
  ! **Library:**   SLATEC (PCHIP)
  !***
  ! **Type:**      SINGLE PRECISION (PCHCS-S, DPCHCS-D)
  !***
  ! **Author:**  Fritsch, F. N., (LLNL)
  !***
  ! **Description:**
  !
  !         PCHCS:  PCHIC Monotonicity Switch Derivative Setter.
  !
  !     Called by  PCHIC  to adjust the values of D in the vicinity of a
  !     switch in direction of monotonicity, to produce a more "visually
  !     pleasing" curve than that given by  PCHIM .
  !
  ! ----------------------------------------------------------------------
  !
  !  Calling sequence:
  !
  !        PARAMETER  (INCFD = ...)
  !        INTEGER  N, IERR
  !        REAL  SWITCH, H(N), SLOPE(N), D(INCFD,N)
  !
  !        CALL  PCHCS (SWITCH, N, H, SLOPE, D, INCFD, IERR)
  !
  !   Parameters:
  !
  !     SWITCH -- (input) indicates the amount of control desired over
  !           local excursions from data.
  !
  !     N -- (input) number of data points.  (assumes N>2 .)
  !
  !     H -- (input) real array of interval lengths.
  !     SLOPE -- (input) real array of data slopes.
  !           If the data are (X(I),Y(I)), I=1(1)N, then these inputs are:
  !                  H(I) =  X(I+1)-X(I),
  !              SLOPE(I) = (Y(I+1)-Y(I))/H(I),  I=1(1)N-1.
  !
  !     D -- (input) real array of derivative values at the data points,
  !           as determined by PCHCI.
  !          (output) derivatives in the vicinity of switches in direction
  !           of monotonicity may be adjusted to produce a more "visually
  !           pleasing" curve.
  !           The value corresponding to X(I) is stored in
  !                D(1+(I-1)*INCFD),  I=1(1)N.
  !           No other entries in D are changed.
  !
  !     INCFD -- (input) increment between successive values in D.
  !           This argument is provided primarily for 2-D applications.
  !
  !     IERR -- (output) error flag.  should be zero.
  !           If negative, trouble in PCHSW.  (should never happen.)
  !
  !    -------
  !    WARNING:  This routine does no validity-checking of arguments.
  !    -------
  !
  !  Fortran intrinsics used:  ABS, MAX, MIN.
  !
  !***
  ! **See also:**  PCHIC
  !***
  ! **Routines called:**  PCHST, PCHSW

  !* REVISION HISTORY  (YYMMDD)
  !   820218  DATE WRITTEN
  !   820617  Redesigned to (1) fix  problem with lack of continuity
  !           approaching a flat-topped peak (2) be cleaner and
  !           easier to verify.
  !           Eliminated subroutines PCHSA and PCHSX in the process.
  !   820622  1. Limited fact to not exceed one, so computed D is a
  !             convex combination of PCHCI value and PCHSD value.
  !           2. Changed fudge from 1 to 4 (based on experiments).
  !   820623  Moved PCHSD to an inline function (eliminating MSWTYP).
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

  !
  !  Programming notes:
  !     1. The function  PCHST(ARG1,ARG2)  is assumed to return zero if
  !        either argument is zero, +1 if they are of the same sign, and
  !        -1 if they are of opposite sign.
  !**End
  !
  !  DECLARE ARGUMENTS.
  !
  INTEGER, INTENT(IN) :: N, Incfd
  INTEGER, INTENT(OUT) :: Ierr
  REAL(SP), INTENT(IN) :: Switch, H(N), Slope(N)
  REAL(SP), INTENT(INOUT) :: D(Incfd,N)
  !
  !  DECLARE LOCAL VARIABLES.
  !
  INTEGER :: i, indx, k, nless1
  REAL(SP) :: del(3), dext, dfloc, dfmx, fact, slmax, wtave(2)
  !
  !  INITIALIZE.
  !
  REAL(SP), PARAMETER :: zero = 0., one = 1.
  REAL(SP), PARAMETER :: fudge = 4.
  !* FIRST EXECUTABLE STATEMENT  PCHCS
  Ierr = 0
  nless1 = N - 1
  !
  !  LOOP OVER SEGMENTS.
  !
  DO i = 2, nless1
    IF( PCHST(Slope(i-1),Slope(i))<0 ) THEN
      !             --------------------------
      !
      !
      !....... SLOPE SWITCHES MONOTONICITY AT I-TH POINT .....................
      !
      !           DO NOT CHANGE D IF 'UP-DOWN-UP'.
      IF( i>2 ) THEN
        IF( PCHST(Slope(i-2),Slope(i))>zero ) CYCLE
        !                   --------------------------
      END IF
      IF( i<nless1 ) THEN
        IF( PCHST(Slope(i+1),Slope(i-1))>zero ) CYCLE
        !                   ----------------------------
      END IF
      !
      !   ....... COMPUTE PROVISIONAL VALUE FOR D(1,I).
      !
      dext = PCHSD(Slope(i-1),Slope(i),H(i-1),H(i))
      !
      !   ....... DETERMINE WHICH INTERVAL CONTAINS THE EXTREMUM.
      !
      IF( PCHST(dext,Slope(i-1))<0 ) THEN
        !                -----------------------
        !
        !              DEXT AND SLOPE(I-1) HAVE OPPOSITE SIGNS --
        !                        EXTREMUM IS IN (X(I-1),X(I)).
        k = i - 1
        !              SET UP TO COMPUTE NEW VALUES FOR D(1,I-1) AND D(1,I).
        wtave(2) = dext
        IF( k>1 ) wtave(1) = PCHSD(Slope(k-1),Slope(k),H(k-1),H(k))
      ELSEIF( PCHST(dext,Slope(i-1))==0 ) THEN
        CYCLE
      ELSE
        !
        !              DEXT AND SLOPE(I) HAVE OPPOSITE SIGNS --
        !                        EXTREMUM IS IN (X(I),X(I+1)).
        k = i
        !              SET UP TO COMPUTE NEW VALUES FOR D(1,I) AND D(1,I+1).
        wtave(1) = dext
        IF( k<nless1 ) wtave(2) = PCHSD(Slope(k),Slope(k+1),H(k),H(k+1))
      END IF
    ELSEIF( PCHST(Slope(i-1),Slope(i))==0 ) THEN
      !
      !
      !....... AT LEAST ONE OF SLOPE(I-1) AND SLOPE(I) IS ZERO --
      !                     CHECK FOR FLAT-TOPPED PEAK .......................
      !
      IF( i==nless1 ) CYCLE
      IF( PCHST(Slope(i-1),Slope(i+1))>=zero ) CYCLE
      !                -----------------------------
      !
      !           WE HAVE FLAT-TOPPED PEAK ON (X(I),X(I+1)).
      k = i
      !           SET UP TO COMPUTE NEW VALUES FOR D(1,I) AND D(1,I+1).
      wtave(1) = PCHSD(Slope(k-1),Slope(k),H(k-1),H(k))
      wtave(2) = PCHSD(Slope(k),Slope(k+1),H(k),H(k+1))
    ELSE
      CYCLE
    END IF
    !
    !
    !....... AT THIS POINT WE HAVE DETERMINED THAT THERE WILL BE AN EXTREMUM
    !        ON (X(K),X(K+1)), WHERE K=I OR I-1, AND HAVE SET ARRAY WTAVE--
    !           WTAVE(1) IS A WEIGHTED AVERAGE OF SLOPE(K-1) AND SLOPE(K),
    !                    IF K>1
    !           WTAVE(2) IS A WEIGHTED AVERAGE OF SLOPE(K) AND SLOPE(K+1),
    !                    IF K<N-1
    !
    slmax = ABS(Slope(k))
    IF( k>1 ) slmax = MAX(slmax,ABS(Slope(k-1)))
    IF( k<nless1 ) slmax = MAX(slmax,ABS(Slope(k+1)))
    !
    IF( k>1 ) del(1) = Slope(k-1)/slmax
    del(2) = Slope(k)/slmax
    IF( k<nless1 ) del(3) = Slope(k+1)/slmax
    !
    IF( (k>1) .AND. (k<nless1) ) THEN
      !           NORMAL CASE -- EXTREMUM IS NOT IN A BOUNDARY INTERVAL.
      fact = fudge*ABS(del(3)*(del(1)-del(2))*(wtave(2)/slmax))
      D(1,k) = D(1,k) + MIN(fact,one)*(wtave(1)-D(1,k))
      fact = fudge*ABS(del(1)*(del(3)-del(2))*(wtave(1)/slmax))
      D(1,k+1) = D(1,k+1) + MIN(fact,one)*(wtave(2)-D(1,k+1))
    ELSE
      !           SPECIAL CASE K=1 (WHICH CAN OCCUR ONLY IF I=2) OR
      !                        K=NLESS1 (WHICH CAN OCCUR ONLY IF I=NLESS1).
      fact = fudge*ABS(del(2))
      D(1,i) = MIN(fact,one)*wtave(i-k+1)
      !              NOTE THAT I-K+1 = 1 IF K=I  (=NLESS1),
      !                        I-K+1 = 2 IF K=I-1(=1).
    END IF
    !
    !
    !....... ADJUST IF NECESSARY TO LIMIT EXCURSIONS FROM DATA.
    !
    IF( Switch>zero ) THEN
      !
      dfloc = H(k)*ABS(Slope(k))
      IF( k>1 ) dfloc = MAX(dfloc,H(k-1)*ABS(Slope(k-1)))
      IF( k<nless1 ) dfloc = MAX(dfloc,H(k+1)*ABS(Slope(k+1)))
      dfmx = Switch*dfloc
      indx = i - k + 1
      !        INDX = 1 IF K=I, 2 IF K=I-1.
      !        ---------------------------------------------------------------
      CALL PCHSW(dfmx,indx,D(1,k),D(1,k+1),H(k),Slope(k),Ierr)
      !        ---------------------------------------------------------------
      IF( Ierr/=0 ) RETURN
    END IF
    !
    !....... END OF SEGMENT LOOP.
    !
  END DO
  !
  !------------- LAST LINE OF PCHCS FOLLOWS ------------------------------
CONTAINS
  !
  !  DEFINE INLINE FUNCTION FOR WEIGHTED AVERAGE OF SLOPES.
  !
  REAL(SP) ELEMENTAL FUNCTION PCHSD(s1,s2,h1,h2)
    REAL(SP), INTENT(IN) :: s1, s2, h1, h2

    PCHSD = (h2/(h1+h2))*s1 + (h1/(h1+h2))*s2

  END FUNCTION PCHSD
  !
END SUBROUTINE PCHCS