!DECK CMPCSG
SUBROUTINE CMPCSG(N,Ijump,Fnum,Fden,A)
  IMPLICIT NONE
  REAL dum, Fden, Fnum, pi, pibyn, PIMACH, x, y
  INTEGER i, Ijump, k, k1, k2, k3, k4, k5, N, np1
  !***BEGIN PROLOGUE  CMPCSG
  !***SUBSIDIARY
  !***PURPOSE  Subsidiary to CMGNBN
  !***LIBRARY   SLATEC
  !***TYPE      COMPLEX (COSGEN-S, CMPCSG-C)
  !***AUTHOR  (UNKNOWN)
  !***DESCRIPTION
  !
  !     This subroutine computes required cosine values in ascending
  !     order.  When IJUMP .GT. 1 the routine computes values
  !
  !        2*COS(J*PI/L), J=1,2,...,L and J .NE. 0(MOD N/IJUMP+1)
  !
  !     where L = IJUMP*(N/IJUMP+1).
  !
  !
  !     when IJUMP = 1 it computes
  !
  !            2*COS((J-FNUM)*PI/(N+FDEN)),  J=1, 2, ... ,N
  !
  !     where
  !        FNUM = 0.5, FDEN = 0.0,  for regular reduction values.
  !        FNUM = 0.0, FDEN = 1.0, for B-R and C-R when ISTAG = 1
  !        FNUM = 0.0, FDEN = 0.5, for B-R and C-R when ISTAG = 2
  !        FNUM = 0.5, FDEN = 0.5, for B-R and C-R when ISTAG = 2
  !                                in CMPOSN only.
  !
  !***SEE ALSO  CMGNBN
  !***ROUTINES CALLED  PIMACH
  !***REVISION HISTORY  (YYMMDD)
  !   801001  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900402  Added TYPE section.  (WRB)
  !***END PROLOGUE  CMPCSG
  COMPLEX A
  DIMENSION A(*)
  !
  !
  !***FIRST EXECUTABLE STATEMENT  CMPCSG
  pi = PIMACH(dum)
  IF ( N/=0 ) THEN
    IF ( Ijump==1 ) THEN
      np1 = N + 1
      y = pi/(N+Fden)
      DO i = 1, N
        x = np1 - i - Fnum
        A(i) = CMPLX(2.*COS(x*y),0.)
      ENDDO
    ELSE
      k3 = N/Ijump + 1
      k4 = k3 - 1
      pibyn = pi/(N+Ijump)
      DO k = 1, Ijump
        k1 = (k-1)*k3
        k5 = (k-1)*k4
        DO i = 1, k4
          x = k1 + i
          k2 = k5 + i
          A(k2) = CMPLX(-2.*COS(x*pibyn),0.)
        ENDDO
      ENDDO
    ENDIF
  ENDIF
END SUBROUTINE CMPCSG
