!*==DDASLV.f90  processed by SPAG 6.72Dc at 11:01 on  6 Feb 2019
!DECK DDASLV
SUBROUTINE DDASLV(Neq,Delta,Wm,Iwm)
  IMPLICIT NONE
  !*--DDASLV5
  !***BEGIN PROLOGUE  DDASLV
  !***SUBSIDIARY
  !***PURPOSE  Linear system solver for DDASSL.
  !***LIBRARY   SLATEC (DASSL)
  !***TYPE      DOUBLE PRECISION (SDASLV-S, DDASLV-D)
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
  !     DGESL IS CALLED.
  !     FOR A BANDED MATRIX,THE LINPACK ROUTINE
  !     DGBSL IS CALLED.
  !-----------------------------------------------------------------------
  !***ROUTINES CALLED  DGBSL, DGESL
  !***REVISION HISTORY  (YYMMDD)
  !   830315  DATE WRITTEN
  !   901009  Finished conversion to SLATEC 4.0 format (F.N.Fritsch)
  !   901019  Merged changes made by C. Ulrich with SLATEC 4.0 format.
  !   901026  Added explicit declarations for all variables and minor
  !           cosmetic changes to prologue.  (FNF)
  !***END PROLOGUE  DDASLV
  !
  INTEGER Neq , Iwm(*)
  REAL(8) :: Delta(*) , Wm(*)
  !
  EXTERNAL DGBSL , DGESL
  !
  INTEGER LIPVT , LML , LMU , LMTYPE , meband , mtype , NPD
  PARAMETER (NPD=1)
  PARAMETER (LML=1)
  PARAMETER (LMU=2)
  PARAMETER (LMTYPE=4)
  PARAMETER (LIPVT=21)
  !
  !***FIRST EXECUTABLE STATEMENT  DDASLV
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
      CALL DGBSL(Wm(NPD),meband,Neq,Iwm(LML),Iwm(LMU),Iwm(LIPVT),Delta,0)
      GOTO 99999
    CASE DEFAULT
  END SELECT
  !
  !     DENSE MATRIX
  CALL DGESL(Wm(NPD),Neq,Neq,Iwm(LIPVT),Delta,0)
  RETURN
  !------END OF SUBROUTINE DDASLV------
  99999 END SUBROUTINE DDASLV
