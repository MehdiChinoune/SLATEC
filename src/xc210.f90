!*==XC210.f90  processed by SPAG 6.72Dc at 11:02 on  6 Feb 2019
!DECK XC210
SUBROUTINE XC210(K,Z,J,Ierror)
  IMPLICIT NONE
  !*--XC2105
  !*** Start of declarations inserted by SPAG
  INTEGER i, ic, id, Ierror, ii, it, ja, ka, ka1, ka2, m, nm1, &
    np1
  !*** End of declarations inserted by SPAG
  !***BEGIN PROLOGUE  XC210
  !***PURPOSE  To provide single-precision floating-point arithmetic
  !            with an extended exponent range.
  !***LIBRARY   SLATEC
  !***CATEGORY  A3D
  !***TYPE      SINGLE PRECISION (XC210-S, DXC210-D)
  !***KEYWORDS  EXTENDED-RANGE SINGLE-PRECISION ARITHMETIC
  !***AUTHOR  Lozier, Daniel W., (National Bureau of Standards)
  !           Smith, John M., (NBS and George Mason University)
  !***DESCRIPTION
  !     INTEGER K, J
  !     REAL Z
  !
  !                  GIVEN K THIS SUBROUTINE COMPUTES J AND Z
  !                  SUCH THAT  RADIX**K = Z*10**J, WHERE Z IS IN
  !                  THE RANGE 1/10 .LE. Z .LT. 1.
  !                  THE VALUE OF Z WILL BE ACCURATE TO FULL
  !                  SINGLE-PRECISION PROVIDED THE NUMBER
  !                  OF DECIMAL PLACES IN THE LARGEST
  !                  INTEGER PLUS THE NUMBER OF DECIMAL
  !                  PLACES CARRIED IN SINGLE-PRECISION DOES NOT
  !                  EXCEED 60. XC210 IS CALLED BY SUBROUTINE
  !                  XCON WHEN NECESSARY. THE USER SHOULD
  !                  NEVER NEED TO CALL XC210 DIRECTLY.
  !
  !***SEE ALSO  XSET
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  XERMSG
  !***COMMON BLOCKS    XBLK3
  !***REVISION HISTORY  (YYMMDD)
  !   820712  DATE WRITTEN
  !   890126  Revised to meet SLATEC CML recommendations.  (DWL and JMS)
  !   901019  Revisions to prologue.  (DWL and WRB)
  !   901106  Changed all specific intrinsics to generic.  (WRB)
  !           Corrected order of sections in prologue and added TYPE
  !           section.  (WRB)
  !           CALLs to XERROR changed to CALLs to XERMSG.  (WRB)
  !   920127  Revised PURPOSE section of prologue.  (DWL)
  !***END PROLOGUE  XC210
  INTEGER K, J
  REAL Z
  INTEGER NLG102, MLG102, LG102
  COMMON /XBLK3 / NLG102, MLG102, LG102(21)
  SAVE /XBLK3 /
  !
  !   THE CONDITIONS IMPOSED ON NLG102, MLG102, AND LG102 BY
  ! THIS SUBROUTINE ARE
  !
  !     (1) NLG102 .GE. 2
  !
  !     (2) MLG102 .GE. 1
  !
  !     (3) 2*MLG102*(MLG102 - 1) .LE. 2**NBITS - 1
  !
  ! THESE CONDITIONS MUST BE MET BY APPROPRIATE CODING
  ! IN SUBROUTINE XSET.
  !
  !***FIRST EXECUTABLE STATEMENT  XC210
  Ierror = 0
  IF ( K==0 ) THEN
    J = 0
    Z = 1.0
    GOTO 99999
  ELSE
    m = MLG102
    ka = ABS(K)
    ka1 = ka/m
    ka2 = MOD(ka,m)
    IF ( ka1<m ) THEN
      nm1 = NLG102 - 1
      np1 = NLG102 + 1
      it = ka2*LG102(np1)
      ic = it/m
      id = MOD(it,m)
      Z = id
      IF ( ka1>0 ) THEN
        DO ii = 1, nm1
          i = np1 - ii
          it = ka2*LG102(i) + ka1*LG102(i+1) + ic
          ic = it/m
          id = MOD(it,m)
          Z = Z/m + id
        ENDDO
        ja = ka*LG102(1) + ka1*LG102(2) + ic
      ELSE
        DO ii = 1, nm1
          i = np1 - ii
          it = ka2*LG102(i) + ic
          ic = it/m
          id = MOD(it,m)
          Z = Z/m + id
        ENDDO
        ja = ka*LG102(1) + ic
      ENDIF
      Z = Z/m
      IF ( K>0 ) THEN
        J = ja + 1
        Z = 10.0**(Z-1.0)
      ELSE
        J = -ja
        Z = 10.0**(-Z)
      ENDIF
      GOTO 99999
    ENDIF
  ENDIF
  !   THIS ERROR OCCURS IF K EXCEEDS  MLG102**2 - 1  IN MAGNITUDE.
  !
  CALL XERMSG('SLATEC','XC210','K too large',108,1)
  Ierror = 108
  RETURN
  99999 CONTINUE
  END SUBROUTINE XC210
