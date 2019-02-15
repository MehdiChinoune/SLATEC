!DECK SDASLV
SUBROUTINE SDASLV(Neq,Delta,Wm,Iwm)
  IMPLICIT NONE
  !***BEGIN PROLOGUE  SDASLV
  !***SUBSIDIARY
  !***PURPOSE  Linear system solver for SDASSL.
  !***LIBRARY   SLATEC (DASSL)
  !***TYPE      SINGLE PRECISION (SDASLV-S, DDASLV-D)
  !***AUTHOR  Petzold, Linda R., (LLNL)
  !***DESCRIPTION
  !-----------------------------------------------------------------------
  !     THIS ROUTINE MANAGES THE SOLUTION OF THE LINEAR
  !     SYSTEM ARISING IN THE NEWTON ITERATION.
  !     MATRICES AND REAL TEMPORARY STORAGE AND
  !     REAL INFORMATION ARE STORED IN THE ARRAY WM.
  !     INTEGER MATRIX INFORMATION IS STORED IN
  !     THE ARRAY IWM.
  !     FOR A DENSE MATRIX, THE LINPACK ROUTINE
  !     SGESL IS CALLED.
  !     FOR A BANDED MATRIX,THE LINPACK ROUTINE
  !     SGBSL IS CALLED.
  !-----------------------------------------------------------------------
  !***ROUTINES CALLED  SGBSL, SGESL
  !***REVISION HISTORY  (YYMMDD)
  !   830315  DATE WRITTEN
  !   901009  Finished conversion to SLATEC 4.0 format (F.N.Fritsch)
  !   901019  Merged changes made by C. Ulrich with SLATEC 4.0 format.
  !   901026  Added explicit declarations for all variables and minor
  !           cosmetic changes to prologue.  (FNF)
  !***END PROLOGUE  SDASLV
  !
  INTEGER Neq, Iwm(*)
  REAL Delta(*), Wm(*)
  !
  EXTERNAL SGBSL, SGESL
  !
  INTEGER LIPVT, LML, LMU, LMTYPE, meband, mtype, NPD
  PARAMETER (NPD=1)
  PARAMETER (LML=1)
  PARAMETER (LMU=2)
  PARAMETER (LMTYPE=4)
  PARAMETER (LIPVT=21)
  !
  !***FIRST EXECUTABLE STATEMENT  SDASLV
  mtype = Iwm(LMTYPE)
  SELECT CASE (mtype)
    CASE (3)
      !
      !     DUMMY SECTION FOR MTYPE=3
      RETURN
    CASE (4,5)
      !
      !     BANDED MATRIX
      meband = 2*Iwm(LML) + Iwm(LMU) + 1
      CALL SGBSL(Wm(NPD),meband,Neq,Iwm(LML),Iwm(LMU),Iwm(LIPVT),Delta,0)
      GOTO 99999
    CASE DEFAULT
  END SELECT
  !
  !     DENSE MATRIX
  CALL SGESL(Wm(NPD),Neq,Neq,Iwm(LIPVT),Delta,0)
  RETURN
  !------END OF SUBROUTINE SDASLV------
  99999 CONTINUE
  END SUBROUTINE SDASLV
