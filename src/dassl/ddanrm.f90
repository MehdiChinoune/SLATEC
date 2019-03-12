!DECK DDANRM
REAL(8) FUNCTION DDANRM(Neq,V,Wt,Rpar,Ipar)
  IMPLICIT NONE
  !***BEGIN PROLOGUE  DDANRM
  !***SUBSIDIARY
  !***PURPOSE  Compute vector norm for DDASSL.
  !***LIBRARY   SLATEC (DASSL)
  !***TYPE      DOUBLE PRECISION (SDANRM-S, DDANRM-D)
  !***AUTHOR  Petzold, Linda R., (LLNL)
  !***DESCRIPTION
  !-----------------------------------------------------------------------
  !     THIS FUNCTION ROUTINE COMPUTES THE WEIGHTED
  !     ROOT-MEAN-SQUARE NORM OF THE VECTOR OF LENGTH
  !     NEQ CONTAINED IN THE ARRAY V,WITH WEIGHTS
  !     CONTAINED IN THE ARRAY WT OF LENGTH NEQ.
  !        DDANRM=SQRT((1/NEQ)*SUM(V(I)/WT(I))**2)
  !-----------------------------------------------------------------------
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   830315  DATE WRITTEN
  !   901009  Finished conversion to SLATEC 4.0 format (F.N.Fritsch)
  !   901019  Merged changes made by C. Ulrich with SLATEC 4.0 format.
  !   901026  Added explicit declarations for all variables and minor
  !           cosmetic changes to prologue.  (FNF)
  !***END PROLOGUE  DDANRM
  !
  INTEGER Neq, Ipar(*)
  REAL(8) :: V(Neq), Wt(Neq), Rpar(*)
  !
  INTEGER i
  REAL(8) :: sum, vmax
  !
  !***FIRST EXECUTABLE STATEMENT  DDANRM
  DDANRM = 0.0D0
  vmax = 0.0D0
  DO i = 1, Neq
    IF ( ABS(V(i)/Wt(i))>vmax ) vmax = ABS(V(i)/Wt(i))
  ENDDO
  IF ( vmax>0.0D0 ) THEN
    sum = 0.0D0
    DO i = 1, Neq
      sum = sum + ((V(i)/Wt(i))/vmax)**2
    ENDDO
    DDANRM = vmax*SQRT(sum/Neq)
  ENDIF
  !------END OF FUNCTION DDANRM------
END FUNCTION DDANRM
