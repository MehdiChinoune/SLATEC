!** C0LGMC
COMPLEX(SP) FUNCTION C0LGMC(Z)
  !> Evaluate (Z+0.5)*LOG((Z+1.)/Z) - 1.0 with relative
  !            accuracy.
  !***
  ! **Library:**   SLATEC (FNLIB)
  !***
  ! **Category:**  C7A
  !***
  ! **Type:**      COMPLEX (C0LGMC-C)
  !***
  ! **Keywords:**  FNLIB, GAMMA FUNCTION, SPECIAL FUNCTIONS
  !***
  ! **Author:**  Fullerton, W., (LANL)
  !***
  ! **Description:**
  !
  ! Evaluate  (Z+0.5)*LOG((Z+1.0)/Z) - 1.0  with relative error accuracy
  ! Let Q = 1.0/Z so that
  !     (Z+0.5)*LOG(1+1/Z) - 1 = (Z+0.5)*(LOG(1+Q) - Q + Q*Q/2) - Q*Q/4
  !        = (Z+0.5)*Q**3*C9LN2R(Q) - Q**2/4,
  ! where  C9LN2R  is (LOG(1+Q) - Q + 0.5*Q**2) / Q**3.
  !
  !***
  ! **References:**  (NONE)
  !***
  ! **Routines called:**  C9LN2R, R1MACH

  !* REVISION HISTORY  (YYMMDD)
  !   780401  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  USE service, ONLY : R1MACH
  COMPLEX(SP) :: Z
  REAL(SP) :: cabsz
  COMPLEX(SP) :: q
  REAL(SP), PARAMETER :: rbig = 1._SP/R1MACH(3)
  !* FIRST EXECUTABLE STATEMENT  C0LGMC
  !
  cabsz = ABS(Z)
  IF( cabsz>rbig ) C0LGMC = -(Z+0.5_SP)*LOG(Z) - Z
  IF( cabsz>rbig ) RETURN
  !
  q = 1._SP/Z
  IF( cabsz<=1.23 ) C0LGMC = (Z+0.5_SP)*LOG(1._SP+q) - 1._SP
  IF( cabsz>1.23 ) C0LGMC = ((1._SP+.5_SP*q)*C9LN2R(q)-.25_SP)*q**2
  !
END FUNCTION C0LGMC
