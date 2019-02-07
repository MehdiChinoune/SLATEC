!*==DDRES2.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK DDRES2
SUBROUTINE DDRES2(T,Y,Yprime,Delta,Ires,Rpar,Ipar)
  IMPLICIT NONE
  !*--DDRES25
  !***BEGIN PROLOGUE  DDRES2
  !***SUBSIDIARY
  !***PURPOSE  Second residual evaluator for DDASQC.
  !***LIBRARY   SLATEC (DASSL)
  !***TYPE      DOUBLE PRECISION (SDRES2-S, DDRES2-D)
  !***AUTHOR  PETZOLD, LINDA R., (LLNL)
  !***SEE ALSO  DDASQC
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   891013  DATE WRITTEN
  !   901001  Converted prologue to 4.0 format and made all argument
  !           declarations explicit.  (FNF)
  !   901030  Made all local declarations explicit.  (FNF)
  !***END PROLOGUE  DDRES2
  INTEGER Ires , Ipar(*)
  DOUBLE PRECISION T , Y(*) , Yprime(*) , Delta(*) , Rpar(*)
  INTEGER i , j , k , ng
  DOUBLE PRECISION alph1 , alph2 , d
  DATA alph1/1.0D0/ , alph2/1.0D0/ , ng/5/
  !***FIRST EXECUTABLE STATEMENT  DDRES2
  DO j = 1 , ng
    DO i = 1 , ng
      k = i + (j-1)*ng
      d = -2.0D0*Y(k)
      IF ( i/=1 ) d = d + Y(k-1)*alph1
      IF ( j/=1 ) d = d + Y(k-ng)*alph2
      Delta(k) = d - Yprime(k)
    ENDDO
  ENDDO
END SUBROUTINE DDRES2
