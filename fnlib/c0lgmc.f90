!*==C0LGMC.f90  processed by SPAG 6.72Dc at 10:56 on  6 Feb 2019
!DECK C0LGMC
      COMPLEX FUNCTION C0LGMC(Z)
      IMPLICIT NONE
!*--C0LGMC5
!*** Start of declarations inserted by SPAG
      REAL cabsz , R1MACH , rbig
!*** End of declarations inserted by SPAG
!***BEGIN PROLOGUE  C0LGMC
!***PURPOSE  Evaluate (Z+0.5)*LOG((Z+1.)/Z) - 1.0 with relative
!            accuracy.
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C7A
!***TYPE      COMPLEX (C0LGMC-C)
!***KEYWORDS  FNLIB, GAMMA FUNCTION, SPECIAL FUNCTIONS
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
! Evaluate  (Z+0.5)*LOG((Z+1.0)/Z) - 1.0  with relative error accuracy
! Let Q = 1.0/Z so that
!     (Z+0.5)*LOG(1+1/Z) - 1 = (Z+0.5)*(LOG(1+Q) - Q + Q*Q/2) - Q*Q/4
!        = (Z+0.5)*Q**3*C9LN2R(Q) - Q**2/4,
! where  C9LN2R  is (LOG(1+Q) - Q + 0.5*Q**2) / Q**3.
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  C9LN2R, R1MACH
!***REVISION HISTORY  (YYMMDD)
!   780401  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!***END PROLOGUE  C0LGMC
      COMPLEX Z , q , C9LN2R
      SAVE rbig
      DATA rbig/0.0/
!***FIRST EXECUTABLE STATEMENT  C0LGMC
      IF ( rbig==0.0 ) rbig = 1.0/R1MACH(3)
!
      cabsz = ABS(Z)
      IF ( cabsz>rbig ) C0LGMC = -(Z+0.5)*LOG(Z) - Z
      IF ( cabsz>rbig ) RETURN
!
      q = 1.0/Z
      IF ( cabsz<=1.23 ) C0LGMC = (Z+0.5)*LOG(1.0+q) - 1.0
      IF ( cabsz>1.23 ) C0LGMC = ((1.+.5*q)*C9LN2R(q)-.25)*q**2
!
      END FUNCTION C0LGMC
