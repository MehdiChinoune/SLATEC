!*==MPCMD.f90  processed by SPAG 6.72Dc at 11:01 on  6 Feb 2019
!DECK MPCMD
SUBROUTINE MPCMD(X,Dz)
  IMPLICIT NONE
  !*--MPCMD5
  !*** Start of declarations inserted by SPAG
  INTEGER i , LUN , M , MXR
  !*** End of declarations inserted by SPAG
  !***BEGIN PROLOGUE  MPCMD
  !***SUBSIDIARY
  !***PURPOSE  Subsidiary to DQDOTA and DQDOTI
  !***LIBRARY   SLATEC
  !***TYPE      ALL (MPCMD-A)
  !***AUTHOR  (UNKNOWN)
  !***DESCRIPTION
  !
  !  Converts multiple-precision X to double-precision DZ. Assumes
  !  X is in allowable range for double-precision numbers. There is
  !  some loss of accuracy if the exponent is large.
  !
  !  The argument X(*) is INTEGER array of size 30.  See the comments in
  !  the routine MPBLAS for the reason for this choice.
  !
  !***SEE ALSO  DQDOTA, DQDOTI, MPBLAS
  !***ROUTINES CALLED  MPCHK, MPERR
  !***COMMON BLOCKS    MPCOM
  !***REVISION HISTORY  (YYMMDD)
  !   791001  DATE WRITTEN
  !   ??????  Modified for use with BLAS.  Blank COMMON changed to named
  !           COMMON.  R given dimension 12.
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900402  Added TYPE section.  (WRB)
  !   930124  Increased Array size in MPCON for SUN -r8.  (RWC)
  !***END PROLOGUE  MPCMD
  REAL(8) :: db , Dz , dz2
  COMMON /MPCOM / B , T , M , LUN , MXR , R(30)
  INTEGER B , T , R , X(*) , tm
  !***FIRST EXECUTABLE STATEMENT  MPCMD
  CALL MPCHK(1,4)
  Dz = 0D0
  IF ( X(1)==0 ) RETURN
  db = REAL(B, 8)
  DO i = 1 , T
    Dz = db*Dz + REAL(X(i+2), 8)
    tm = i
    ! CHECK IF FULL DOUBLE-PRECISION ACCURACY ATTAINED
    dz2 = Dz + 1D0
    ! TEST BELOW NOT ALWAYS EQUIVALENT TO - IF (DZ2.LE.DZ) GO TO 20,
    ! FOR EXAMPLE ON CYBER 76.
    IF ( (dz2-Dz)<=0D0 ) EXIT
  ENDDO
  ! NOW ALLOW FOR EXPONENT
  Dz = Dz*(db**(X(2)-tm))
  ! CHECK REASONABLENESS OF RESULT.
  IF ( Dz>0D0 ) THEN
    ! LHS SHOULD BE .LE. 0.5 BUT ALLOW FOR SOME ERROR IN LOG
    IF ( ABS(REAL(X(2), 8)-(LOG(Dz)/LOG(REAL(B, 8))+0.5D0))<=0.6D0 ) THEN
      IF ( X(1)<0 ) Dz = -Dz
      RETURN
    ENDIF
  ENDIF
  ! FOLLOWING MESSAGE INDICATES THAT X IS TOO LARGE OR SMALL -
  ! TRY USING MPCMDE INSTEAD.
  WRITE (LUN,99001)
  99001 FORMAT (' *** FLOATING-POINT OVER/UNDER-FLOW IN MPCMD ***')
  CALL MPERR
END SUBROUTINE MPCMD
