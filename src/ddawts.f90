!*==DDAWTS.f90  processed by SPAG 6.72Dc at 11:01 on  6 Feb 2019
!DECK DDAWTS
SUBROUTINE DDAWTS(Neq,Iwt,Rtol,Atol,Y,Wt,Rpar,Ipar)
  IMPLICIT NONE
  !*--DDAWTS5
  !***BEGIN PROLOGUE  DDAWTS
  !***SUBSIDIARY
  !***PURPOSE  Set error weight vector for DDASSL.
  !***LIBRARY   SLATEC (DASSL)
  !***TYPE      DOUBLE PRECISION (SDAWTS-S, DDAWTS-D)
  !***AUTHOR  Petzold, Linda R., (LLNL)
  !***DESCRIPTION
  !-----------------------------------------------------------------------
  !     THIS SUBROUTINE SETS THE ERROR WEIGHT VECTOR
  !     WT ACCORDING TO WT(I)=RTOL(I)*ABS(Y(I))+ATOL(I),
  !     I=1,-,N.
  !     RTOL AND ATOL ARE SCALARS IF IWT = 0,
  !     AND VECTORS IF IWT = 1.
  !-----------------------------------------------------------------------
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   830315  DATE WRITTEN
  !   901009  Finished conversion to SLATEC 4.0 format (F.N.Fritsch)
  !   901019  Merged changes made by C. Ulrich with SLATEC 4.0 format.
  !   901026  Added explicit declarations for all variables and minor
  !           cosmetic changes to prologue.  (FNF)
  !***END PROLOGUE  DDAWTS
  !
  INTEGER Neq , Iwt , Ipar(*)
  DOUBLE PRECISION Rtol(*) , Atol(*) , Y(*) , Wt(*) , Rpar(*)
  !
  INTEGER i
  DOUBLE PRECISION atoli , rtoli
  !
  !***FIRST EXECUTABLE STATEMENT  DDAWTS
  rtoli = Rtol(1)
  atoli = Atol(1)
  DO i = 1 , Neq
    IF ( Iwt/=0 ) THEN
      rtoli = Rtol(i)
      atoli = Atol(i)
    ENDIF
    Wt(i) = rtoli*ABS(Y(i)) + atoli
  ENDDO
  !-----------END OF SUBROUTINE DDAWTS------------------------------------
END SUBROUTINE DDAWTS
