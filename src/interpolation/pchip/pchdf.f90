!** PCHDF
REAL FUNCTION PCHDF(K,X,S,Ierr)
  !>
  !***
  !  Computes divided differences for PCHCE and PCHSP
  !***
  ! **Library:**   SLATEC (PCHIP)
  !***
  ! **Type:**      SINGLE PRECISION (PCHDF-S, DPCHDF-D)
  !***
  ! **Author:**  Fritsch, F. N., (LLNL)
  !***
  ! **Description:**
  !
  !          PCHDF:   PCHIP Finite Difference Formula
  !
  !     Uses a divided difference formulation to compute a K-point approx-
  !     imation to the derivative at X(K) based on the data in X and S.
  !
  !     Called by  PCHCE  and  PCHSP  to compute 3- and 4-point boundary
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
  !        IERR   will be set to -1 if K.LT.2 .
  !        PCHDF  will be set to the desired derivative approximation if
  !               IERR=0 or to zero if IERR=-1.
  !
  ! ----------------------------------------------------------------------
  !
  !***
  ! **See also:**  PCHCE, PCHSP
  !***
  ! **References:**  Carl de Boor, A Practical Guide to Splines, Springer-
  !                 Verlag, New York, 1978, pp. 10-16.
  !***
  ! **Routines called:**  XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   820503  DATE WRITTEN
  !   820805  Converted to SLATEC library version.
  !   870813  Minor cosmetic changes.
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
  INTEGER K, Ierr
  REAL X(K), S(K)
  !
  !  DECLARE LOCAL VARIABLES.
  !
  INTEGER i, j
  REAL value
  REAL, PARAMETER :: zero = 0.
  !
  !  CHECK FOR LEGAL VALUE OF K.
  !
  !* FIRST EXECUTABLE STATEMENT  PCHDF
  IF ( K<3 ) THEN
    !
    !  ERROR RETURN.
    !
    !     K.LT.3 RETURN.
    Ierr = -1
    CALL XERMSG('SLATEC','PCHDF','K LESS THAN THREE',Ierr,1)
    PCHDF = zero
    RETURN
  END IF
  !
  !  COMPUTE COEFFICIENTS OF INTERPOLATING POLYNOMIAL.
  !
  DO j = 2, K - 1
    DO i = 1, K - j
      S(i) = (S(i+1)-S(i))/(X(i+j)-X(i))
    END DO
  END DO
  !
  !  EVALUATE DERIVATIVE AT X(K).
  !
  value = S(1)
  DO i = 2, K - 1
    value = S(i) + value*(X(K)-X(i))
  END DO
  !
  !  NORMAL RETURN.
  !
  Ierr = 0
  PCHDF = value
  RETURN
  !------------- LAST LINE OF PCHDF FOLLOWS ------------------------------
  RETURN
END FUNCTION PCHDF
