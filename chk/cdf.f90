!*==CDF.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK CDF
      SUBROUTINE CDF(N,T,Y,Yp)
      IMPLICIT NONE
!*--CDF5
!***BEGIN PROLOGUE  CDF
!***SUBSIDIARY
!***PURPOSE  Quick check for SLATEC routines CDRIV1, CDRIV2 and CDRIV3.
!***LIBRARY   SLATEC (SDRIVE)
!***CATEGORY  I1A2, I1A1B
!***TYPE      COMPLEX (SDF-S, DDF-D, CDF-C)
!***KEYWORDS  CDRIV1, CDRIV2, CDRIV3, QUICK CHECK, SDRIVE
!***AUTHOR  Kahaner, D. K., (NIST)
!             National Institute of Standards and Technology
!             Gaithersburg, MD  20899
!           Sutherland, C. D., (LANL)
!             Mail Stop D466
!             Los Alamos National Laboratory
!             Los Alamos, NM  87545
!***SEE ALSO  CDQCK
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   890405  DATE WRITTEN
!   890405  Revised to meet SLATEC standards.
!***END PROLOGUE  CDF
      REAL T
      COMPLEX alfa , Y(*) , Yp(*)
      INTEGER N
!***FIRST EXECUTABLE STATEMENT  CDF
      alfa = Y(N+1)
      Yp(1) = 1.E0 + alfa*(Y(2)-Y(1)) - Y(1)*Y(3)
      Yp(2) = alfa*(Y(1)-Y(2)) - Y(2)*Y(3)
      Yp(3) = 1.E0 - Y(3)*(Y(1)+Y(2))
      END SUBROUTINE CDF
