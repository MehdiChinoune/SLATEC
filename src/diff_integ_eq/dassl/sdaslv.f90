!** SDASLV
SUBROUTINE SDASLV(Neq,Delta,Wm,Iwm)
  !>
  !  Linear system solver for SDASSL.
  !***
  ! **Library:**   SLATEC (DASSL)
  !***
  ! **Type:**      SINGLE PRECISION (SDASLV-S, DDASLV-D)
  !***
  ! **Author:**  Petzold, Linda R., (LLNL)
  !***
  ! **Description:**
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
  !***
  ! **Routines called:**  SGBSL, SGESL

  !* REVISION HISTORY  (YYMMDD)
  !   830315  DATE WRITTEN
  !   901009  Finished conversion to SLATEC 4.0 format (F.N.Fritsch)
  !   901019  Merged changes made by C. Ulrich with SLATEC 4.0 format.
  !   901026  Added explicit declarations for all variables and minor
  !           cosmetic changes to prologue.  (FNF)
  USE lapack, ONLY : SGBTRS, SGETRS
  !
  INTEGER :: Neq, Iwm(:)
  REAL(SP) :: Wm(:)
  REAL(SP), TARGET :: Delta(Neq)
  !
  INTEGER :: meband, mtype, info
  REAL(SP), POINTER :: delta2(:,:)
  INTEGER, PARAMETER :: LML = 1
  INTEGER, PARAMETER :: LMU = 2
  INTEGER, PARAMETER :: LMTYPE = 4
  INTEGER, PARAMETER :: LIPVT = 21
  !
  !* FIRST EXECUTABLE STATEMENT  SDASLV
  mtype = Iwm(LMTYPE)
  delta2(1:Neq,1:1) => Delta
  SELECT CASE (mtype)
    CASE (3)
      !
      !     DUMMY SECTION FOR MTYPE=3
      RETURN
    CASE (4,5)
      !
      !     BANDED MATRIX
      meband = 2*Iwm(LML) + Iwm(LMU) + 1
      CALL SGBTRS('N',Neq,Iwm(LML),Iwm(LMU),1,Wm,meband,Iwm(LIPVT:),delta2,Neq,info)
      RETURN
    CASE DEFAULT
  END SELECT
  !
  !     DENSE MATRIX
  CALL SGETRS('N',Neq,1,Wm,Neq,Iwm(LIPVT:LIPVT+Neq-1),delta2,Neq,info)
  RETURN
  !------END OF SUBROUTINE SDASLV------
  RETURN
END SUBROUTINE SDASLV
