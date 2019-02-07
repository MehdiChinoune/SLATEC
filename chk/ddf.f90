!*==DDF.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK DDF
SUBROUTINE DDF(N,T,Y,Yp)
  IMPLICIT NONE
  !*--DDF5
  !***BEGIN PROLOGUE  DDF
  !***SUBSIDIARY
  !***PURPOSE  Quick check for SLATEC routines DDRIV1, DDRIV2 and DDRIV3.
  !***LIBRARY   SLATEC (SDRIVE)
  !***CATEGORY  I1A2, I1A1B
  !***TYPE      DOUBLE PRECISION (SDF-S, DDF-D, CDF-C)
  !***KEYWORDS  DDRIV1, DDRIV2, DDRIV3, QUICK CHECK, SDRIVE
  !***AUTHOR  Kahaner, D. K., (NIST)
  !             National Institute of Standards and Technology
  !             Gaithersburg, MD  20899
  !           Sutherland, C. D., (LANL)
  !             Mail Stop D466
  !             Los Alamos National Laboratory
  !             Los Alamos, NM  87545
  !***SEE ALSO  DDQCK
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   890405  DATE WRITTEN
  !   890405  Revised to meet SLATEC standards.
  !***END PROLOGUE  DDF
  REAL(8) :: alfa , T , Y(*) , Yp(*)
  INTEGER N
  !***FIRST EXECUTABLE STATEMENT  DDF
  alfa = Y(N+1)
  Yp(1) = 1.D0 + alfa*(Y(2)-Y(1)) - Y(1)*Y(3)
  Yp(2) = alfa*(Y(1)-Y(2)) - Y(2)*Y(3)
  Yp(3) = 1.D0 - Y(3)*(Y(1)+Y(2))
END SUBROUTINE DDF
