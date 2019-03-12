!DECK SDANRM
REAL FUNCTION SDANRM(Neq,V,Wt,Rpar,Ipar)
  IMPLICIT NONE
  !***BEGIN PROLOGUE  SDANRM
  !***SUBSIDIARY
  !***PURPOSE  Compute vector norm for SDASSL.
  !***LIBRARY   SLATEC (DASSL)
  !***TYPE      SINGLE PRECISION (SDANRM-S, DDANRM-D)
  !***AUTHOR  Petzold, Linda R., (LLNL)
  !***DESCRIPTION
  !-----------------------------------------------------------------------
  !     THIS FUNCTION ROUTINE COMPUTES THE WEIGHTED
  !     ROOT-MEAN-SQUARE NORM OF THE VECTOR OF LENGTH
  !     NEQ CONTAINED IN THE ARRAY V,WITH WEIGHTS
  !     CONTAINED IN THE ARRAY WT OF LENGTH NEQ.
  !        SDANRM=SQRT((1/NEQ)*SUM(V(I)/WT(I))**2)
  !-----------------------------------------------------------------------
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   830315  DATE WRITTEN
  !   901009  Finished conversion to SLATEC 4.0 format (F.N.Fritsch)
  !   901019  Merged changes made by C. Ulrich with SLATEC 4.0 format.
  !   901026  Added explicit declarations for all variables and minor
  !           cosmetic changes to prologue.  (FNF)
  !***END PROLOGUE  SDANRM
  !
  INTEGER Neq, Ipar(*)
  REAL V(Neq), Wt(Neq), Rpar(*)
  !
  INTEGER i
  REAL sum, vmax
  !
  !***FIRST EXECUTABLE STATEMENT  SDANRM
  SDANRM = 0.0E0
  vmax = 0.0E0
  DO i = 1, Neq
    IF ( ABS(V(i)/Wt(i))>vmax ) vmax = ABS(V(i)/Wt(i))
  ENDDO
  IF ( vmax>0.0E0 ) THEN
    sum = 0.0E0
    DO i = 1, Neq
      sum = sum + ((V(i)/Wt(i))/vmax)**2
    ENDDO
    SDANRM = vmax*SQRT(sum/Neq)
  ENDIF
  !------END OF FUNCTION SDANRM------
END FUNCTION SDANRM
