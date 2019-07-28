!** DDANRM
REAL(DP) PURE FUNCTION DDANRM(Neq,V,Wt)
  !> Compute vector norm for DDASSL.
  !***
  ! **Library:**   SLATEC (DASSL)
  !***
  ! **Type:**      DOUBLE PRECISION (SDANRM-S, DDANRM-D)
  !***
  ! **Author:**  Petzold, Linda R., (LLNL)
  !***
  ! **Description:**
  !-----------------------------------------------------------------------
  !     THIS FUNCTION ROUTINE COMPUTES THE WEIGHTED
  !     ROOT-MEAN-SQUARE NORM OF THE VECTOR OF LENGTH
  !     NEQ CONTAINED IN THE ARRAY V,WITH WEIGHTS
  !     CONTAINED IN THE ARRAY WT OF LENGTH NEQ.
  !        DDANRM=SQRT((1/NEQ)*SUM(V(I)/WT(I))**2)
  !-----------------------------------------------------------------------
  !***
  ! **Routines called:**  (NONE)

  !* REVISION HISTORY  (YYMMDD)
  !   830315  DATE WRITTEN
  !   901009  Finished conversion to SLATEC 4.0 format (F.N.Fritsch)
  !   901019  Merged changes made by C. Ulrich with SLATEC 4.0 format.
  !   901026  Added explicit declarations for all variables and minor
  !           cosmetic changes to prologue.  (FNF)
  !
  INTEGER, INTENT(IN) :: Neq
  REAL(DP), INTENT(IN) :: V(Neq), Wt(Neq)
  !* FIRST EXECUTABLE STATEMENT  DDANRM
  DDANRM = NORM2( V/Wt ) / SQRT(1._DP*Neq)
  !------END OF FUNCTION DDANRM------
END FUNCTION DDANRM