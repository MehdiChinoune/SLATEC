!** DDANRM
REAL(8) FUNCTION DDANRM(Neq,V,Wt)
  !>
  !  Compute vector norm for DDASSL.
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
  INTEGER Neq
  REAL(8) :: V(Neq), Wt(Neq)
  !
  INTEGER i
  REAL(8) :: summ, vmax
  !
  !* FIRST EXECUTABLE STATEMENT  DDANRM
  DDANRM = 0.0D0
  vmax = 0.0D0
  DO i = 1, Neq
    IF ( ABS(V(i)/Wt(i))>vmax ) vmax = ABS(V(i)/Wt(i))
  END DO
  IF ( vmax>0.0D0 ) THEN
    summ = 0.0D0
    DO i = 1, Neq
      summ = summ + ((V(i)/Wt(i))/vmax)**2
    END DO
    DDANRM = vmax*SQRT(summ/Neq)
  END IF
  !------END OF FUNCTION DDANRM------
END FUNCTION DDANRM
