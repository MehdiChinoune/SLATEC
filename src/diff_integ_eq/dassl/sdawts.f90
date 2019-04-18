!** SDAWTS
SUBROUTINE SDAWTS(Neq,Iwt,Rtol,Atol,Y,Wt,Rpar,Ipar)
  !>
  !  Set error weight vector for SDASSL.
  !***
  ! **Library:**   SLATEC (DASSL)
  !***
  ! **Type:**      SINGLE PRECISION (SDAWTS-S, DDAWTS-D)
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
  INTEGER Neq, Iwt, Ipar(*)
  REAL Rtol(*), Atol(*), Y(*), Wt(*), Rpar(*)
  !
  INTEGER i
  REAL atoli, rtoli
  !
  !* FIRST EXECUTABLE STATEMENT  SDAWTS
  rtoli = Rtol(1)
  atoli = Atol(1)
  DO i = 1, Neq
    IF ( Iwt/=0 ) THEN
      rtoli = Rtol(i)
      atoli = Atol(i)
    END IF
    Wt(i) = rtoli*ABS(Y(i)) + atoli
  END DO
  !-----------END OF SUBROUTINE SDAWTS------------------------------------
END SUBROUTINE SDAWTS
