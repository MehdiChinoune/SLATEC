!** DDAWTS
PURE SUBROUTINE DDAWTS(Neq,Iwt,Rtol,Atol,Y,Wt)
  !> Set error weight vector for DDASSL.
  !***
  ! **Library:**   SLATEC (DASSL)
  !***
  ! **Type:**      DOUBLE PRECISION (SDAWTS-S, DDAWTS-D)
  !***
  ! **Author:**  Petzold, Linda R., (LLNL)
  !***
  ! **Description:**
  !-----------------------------------------------------------------------
  !     THIS SUBROUTINE SETS THE ERROR WEIGHT VECTOR
  !     WT ACCORDING TO WT(I)=RTOL(I)*ABS(Y(I))+ATOL(I),
  !     I=1,-,N.
  !     RTOL AND ATOL ARE SCALARS IF IWT = 0,
  !     AND VECTORS IF IWT = 1.
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
  INTEGER, INTENT(IN) :: Neq, Iwt
  REAL(DP), INTENT(IN) :: Rtol(Neq), Atol(Neq), Y(Neq)
  REAL(DP), INTENT(OUT) :: Wt(Neq)
  !
  !* FIRST EXECUTABLE STATEMENT  DDAWTS
  IF( Iwt==0 ) THEN
    Wt = Rtol(1)*ABS(Y) + Atol(1)
  ELSE
    Wt = Rtol*ABS(Y) + Atol
  END IF
  !-----------END OF SUBROUTINE DDAWTS------------------------------------
END SUBROUTINE DDAWTS