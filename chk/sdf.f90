!*==SDF.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK SDF
      SUBROUTINE SDF(N,T,Y,Yp)
      IMPLICIT NONE
!*--SDF5
!***BEGIN PROLOGUE  SDF
!***SUBSIDIARY
!***PURPOSE  Quick check for SLATEC routines SDRIV1, SDRIV2 and SDRIV3.
!***LIBRARY   SLATEC (SDRIVE)
!***CATEGORY  I1A2, I1A1B
!***TYPE      SINGLE PRECISION (SDF-S, DDF-D, CDF-C)
!***KEYWORDS  QUICK CHECK, SDRIV1, SDRIV2, SDRIV3
!***AUTHOR  Kahaner, D. K., (NIST)
!             National Institute of Standards and Technology
!             Gaithersburg, MD  20899
!           Sutherland, C. D., (LANL)
!             Mail Stop D466
!             Los Alamos National Laboratory
!             Los Alamos, NM  87545
!***SEE ALSO  SDQCK
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   890405  DATE WRITTEN
!   890405  Revised to meet SLATEC standards.
!***END PROLOGUE  SDF
      REAL alfa , T , Y(*) , Yp(*)
      INTEGER N
!***FIRST EXECUTABLE STATEMENT  SDF
      alfa = Y(N+1)
      Yp(1) = 1.E0 + alfa*(Y(2)-Y(1)) - Y(1)*Y(3)
      Yp(2) = alfa*(Y(1)-Y(2)) - Y(2)*Y(3)
      Yp(3) = 1.E0 - Y(3)*(Y(1)+Y(2))
      END SUBROUTINE SDF
