*DECK DJAC
      SUBROUTINE DJAC (T, U, PD, NROWPD, RPAR, IPAR)
C***BEGIN PROLOGUE  DJAC
C***SUBSIDIARY
C***PURPOSE  Evaluate Jacobian for DDEBDF quick check.
C***LIBRARY   SLATEC
C***TYPE      DOUBLE PRECISION (JAC-S, DJAC-D)
C***AUTHOR  Chow, Jeff (LANL)
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   810801  DATE WRITTEN
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900415  Minor clean-up of prologue and code and name changed from
C           DDJAC to DJAC.  (WRB)
C***END PROLOGUE  DJAC
      INTEGER IPAR, NROWPD
      DOUBLE PRECISION PD, R, R5, RPAR, RSQ, T, U, U1SQ, U2SQ, U1U2
      DIMENSION U(*), PD(NROWPD,*), RPAR(*), IPAR(*)
C***FIRST EXECUTABLE STATEMENT  DJAC
      U1SQ = U(1)*U(1)
      U2SQ = U(2)*U(2)
      U1U2 = U(1)*U(2)
      RSQ = U1SQ + U2SQ
      R = SQRT(RSQ)
      R5 = RSQ*RSQ*R
      PD(3,1) = (3.D0*U1SQ - RSQ)/R5
      PD(4,1) = 3.D0*U1U2/R5
      PD(3,2) = PD(4,1)
      PD(4,2) = (3.D0*U2SQ - RSQ)/R5
      PD(1,3) = 1.D0
      PD(2,4) = 1.D0
      RETURN
      END
