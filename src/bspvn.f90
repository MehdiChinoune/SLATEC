!*==BSPVN.f90  processed by SPAG 6.72Dc at 11:01 on  6 Feb 2019
!DECK BSPVN
      SUBROUTINE BSPVN(T,Jhigh,K,Index,X,Ileft,Vnikx,Work,Iwork)
      IMPLICIT NONE
!*--BSPVN5
!***BEGIN PROLOGUE  BSPVN
!***PURPOSE  Calculate the value of all (possibly) nonzero basis
!            functions at X.
!***LIBRARY   SLATEC
!***CATEGORY  E3, K6
!***TYPE      SINGLE PRECISION (BSPVN-S, DBSPVN-D)
!***KEYWORDS  EVALUATION OF B-SPLINE
!***AUTHOR  Amos, D. E., (SNLA)
!***DESCRIPTION
!
!     Written by Carl de Boor and modified by D. E. Amos
!
!     Abstract
!         BSPVN is the BSPLVN routine of the reference.
!
!         BSPVN calculates the value of all (possibly) nonzero basis
!         functions at X of order MAX(JHIGH,(J+1)*(INDEX-1)), where
!         T(K) .LE. X .LE. T(N+1) and J=IWORK is set inside the routine
!         on the first call when INDEX=1.  ILEFT is such that T(ILEFT)
!         .LE. X .LT. T(ILEFT+1).  A call to INTRV(T,N+1,X,ILO,ILEFT,
!         MFLAG) produces the proper ILEFT.  BSPVN calculates using the
!         basic algorithm needed in BSPVD.  If only basis functions are
!         desired, setting JHIGH=K and INDEX=1 can be faster than
!         calling BSPVD, but extra coding is required for derivatives
!         (INDEX=2) and BSPVD is set up for this purpose.
!
!         Left limiting values are set up as described in BSPVD.
!
!     Description of Arguments
!         Input
!          T       - knot vector of length N+K, where
!                    N = number of B-spline basis functions
!                    N = sum of knot multiplicities-K
!          JHIGH   - order of B-spline, 1 .LE. JHIGH .LE. K
!          K       - highest possible order
!          INDEX   - INDEX = 1 gives basis functions of order JHIGH
!                          = 2 denotes previous entry with WORK, IWORK
!                              values saved for subsequent calls to
!                              BSPVN.
!          X       - argument of basis functions,
!                    T(K) .LE. X .LE. T(N+1)
!          ILEFT   - largest integer such that
!                    T(ILEFT) .LE. X .LT. T(ILEFT+1)
!
!         Output
!          VNIKX   - vector of length K for spline values.
!          WORK    - a work vector of length 2*K
!          IWORK   - a work parameter.  Both WORK and IWORK contain
!                    information necessary to continue for INDEX = 2.
!                    When INDEX = 1 exclusively, these are scratch
!                    variables and can be used for other purposes.
!
!     Error Conditions
!         Improper input is a fatal error.
!
!***REFERENCES  Carl de Boor, Package for calculating with B-splines,
!                 SIAM Journal on Numerical Analysis 14, 3 (June 1977),
!                 pp. 441-472.
!***ROUTINES CALLED  XERMSG
!***REVISION HISTORY  (YYMMDD)
!   800901  DATE WRITTEN
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  BSPVN
!
      INTEGER Ileft , imjp1 , Index , ipj , Iwork , Jhigh , jp1 , jp1ml , K , l
      REAL T , vm , vmprev , Vnikx , Work , X
!     DIMENSION T(ILEFT+JHIGH)
      DIMENSION T(*) , Vnikx(*) , Work(*)
!     CONTENT OF J, DELTAM, DELTAP IS EXPECTED UNCHANGED BETWEEN CALLS.
!     WORK(I) = DELTAP(I), WORK(K+I) = DELTAM(I), I = 1,K
!***FIRST EXECUTABLE STATEMENT  BSPVN
      IF ( K<1 ) THEN
!
!
        CALL XERMSG('SLATEC','BSPVN','K DOES NOT SATISFY K.GE.1',2,1)
        RETURN
      ELSEIF ( Jhigh>K.OR.Jhigh<1 ) THEN
        CALL XERMSG('SLATEC','BSPVN','JHIGH DOES NOT SATISFY 1.LE.JHIGH.LE.K',2,
     &              1)
        RETURN
      ELSEIF ( Index<1.OR.Index>2 ) THEN
        CALL XERMSG('SLATEC','BSPVN','INDEX IS NOT 1 OR 2',2,1)
        RETURN
      ELSEIF ( X<T(Ileft).OR.X>T(Ileft+1) ) THEN
        CALL XERMSG('SLATEC','BSPVN',
     &              'X DOES NOT SATISFY T(ILEFT).LE.X.LE.T(ILEFT+1)',2,1)
        GOTO 99999
      ELSE
        IF ( Index/=2 ) THEN
          Iwork = 1
          Vnikx(1) = 1.0E0
          IF ( Iwork>=Jhigh ) GOTO 100
        ENDIF
        DO
!
          ipj = Ileft + Iwork
          Work(Iwork) = T(ipj) - X
          imjp1 = Ileft - Iwork + 1
          Work(K+Iwork) = X - T(imjp1)
          vmprev = 0.0E0
          jp1 = Iwork + 1
          DO l = 1 , Iwork
            jp1ml = jp1 - l
            vm = Vnikx(l)/(Work(l)+Work(K+jp1ml))
            Vnikx(l) = vm*Work(l) + vmprev
            vmprev = vm*Work(K+jp1ml)
          ENDDO
          Vnikx(jp1) = vmprev
          Iwork = jp1
          IF ( Iwork>=Jhigh ) EXIT
        ENDDO
      ENDIF
!
 100  RETURN
99999 END SUBROUTINE BSPVN
