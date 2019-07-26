!** DPCHDF
REAL(DP) PURE FUNCTION DPCHDF(K,X,S)
  !> Computes divided differences for DPCHCE and DPCHSP
  !***
  ! **Library:**   SLATEC (PCHIP)
  !***
  ! **Type:**      DOUBLE PRECISION (PCHDF-S, DPCHDF-D)
  !***
  ! **Author:**  Fritsch, F. N., (LLNL)
  !***
  ! **Description:**
  !
  !          DPCHDF:   DPCHIP Finite Difference Formula
  !
  !     Uses a divided difference formulation to compute a K-point approx-
  !     imation to the derivative at X(K) based on the data in X and S.
  !
  !     Called by  DPCHCE  and  DPCHSP  to compute 3- and 4-point boundary
  !     derivative approximations.
  !
  ! ----------------------------------------------------------------------
  !
  !     On input:
  !        K      is the order of the desired derivative approximation.
  !               K must be at least 3 (error return if not).
  !        X      contains the K values of the independent variable.
  !               X need not be ordered, but the values **MUST** be
  !               distinct.  (Not checked here.)
  !        S      contains the associated slope values:
  !                  S(I) = (F(I+1)-F(I))/(X(I+1)-X(I)), I=1(1)K-1.
  !               (Note that S need only be of length K-1.)
  !
  !     On return:
  !        S      will be destroyed.
  !        IERR   will be set to -1 if K<2 .
  !        DPCHDF  will be set to the desired derivative approximation if
  !               IERR=0 or to zero if IERR=-1.
  !
  ! ----------------------------------------------------------------------
  !
  !***
  ! **See also:**  DPCHCE, DPCHSP
  !***
  ! **References:**  Carl de Boor, A Practical Guide to Splines, Springer-
  !                 Verlag, New York, 1978, pp. 10-16.
  !***
  ! **Routines called:**  XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   820503  DATE WRITTEN
  !   820805  Converted to SLATEC library version.
  !   870707  Corrected XERROR calls for d.p. name(s).
  !   870813  Minor cosmetic changes.
  !   890206  Corrected XERROR calls.
  !   890411  Added SAVE statements (Vers. 3.2).
  !   890411  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900328  Added TYPE section.  (WRB)
  !   910408  Updated AUTHOR and DATE WRITTEN sections in prologue.  (WRB)
  !   920429  Revised format and order of references.  (WRB,FNF)
  !   930503  Improved purpose.  (FNF)

  !
  !**End
  !
  !  DECLARE ARGUMENTS.
  !
  INTEGER, INTENT(IN) :: K
  REAL(DP), INTENT(IN) :: X(K), S(K)
  !
  !  DECLARE LOCAL VARIABLES.
  !
  INTEGER :: i, j
  REAL(DP) :: value, s_tmp(K)
  !
  !  CHECK FOR LEGAL VALUE OF K.
  !
  !* FIRST EXECUTABLE STATEMENT  DPCHDF
  IF( K<3 ) THEN
    !
    !  ERROR RETURN.
    !
    !     K<3 RETURN.
    ERROR STOP 'DPCHDF : K LESS THAN THREE'
  END IF
  !
  !  COMPUTE COEFFICIENTS OF INTERPOLATING POLYNOMIAL.
  !
  DO j = 2, K - 1
    DO i = 1, K - j
      s_tmp(i) = (S(i+1)-S(i))/(X(i+j)-X(i))
    END DO
  END DO
  !
  !  EVALUATE DERIVATIVE AT X(K).
  !
  value = s_tmp(1)
  DO i = 2, K - 1
    value = s_tmp(i) + value*(X(K)-X(i))
  END DO
  !
  !  NORMAL RETURN.
  !
  DPCHDF = value
  RETURN
  !------------- LAST LINE OF DPCHDF FOLLOWS -----------------------------
END FUNCTION DPCHDF