*DECK SDF
      SUBROUTINE SDF (N, T, Y, YP)
C***BEGIN PROLOGUE  SDF
C***SUBSIDIARY
C***PURPOSE  Quick check for SLATEC routines SDRIV1, SDRIV2 and SDRIV3.
C***LIBRARY   SLATEC (SDRIVE)
C***CATEGORY  I1A2, I1A1B
C***TYPE      SINGLE PRECISION (SDF-S, DDF-D, CDF-C)
C***KEYWORDS  QUICK CHECK, SDRIV1, SDRIV2, SDRIV3
C***AUTHOR  Kahaner, D. K., (NIST)
C             National Institute of Standards and Technology
C             Gaithersburg, MD  20899
C           Sutherland, C. D., (LANL)
C             Mail Stop D466
C             Los Alamos National Laboratory
C             Los Alamos, NM  87545
C***SEE ALSO  SDQCK
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   890405  DATE WRITTEN
C   890405  Revised to meet SLATEC standards.
C***END PROLOGUE  SDF
      REAL ALFA, T, Y(*), YP(*)
      INTEGER N
C***FIRST EXECUTABLE STATEMENT  SDF
      ALFA = Y(N+1)
      YP(1) = 1.E0 + ALFA*(Y(2) - Y(1)) - Y(1)*Y(3)
      YP(2) = ALFA*(Y(1) - Y(2)) - Y(2)*Y(3)
      YP(3) = 1.E0 - Y(3)*(Y(1) + Y(2))
      RETURN
      END
