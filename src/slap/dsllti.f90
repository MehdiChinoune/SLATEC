!** DSLLTI
SUBROUTINE DSLLTI(N,B,X,Nelt,Ia,Ja,A,Isym,Rwork,Iwork)
  IMPLICIT NONE
  !>
  !***
  !  SLAP MSOLVE for LDL' (IC) Factorization.
  !            This routine acts as an interface between the SLAP generic
  !            MSOLVE calling convention and the routine that actually
  !                           -1
  !            computes (LDL')  B = X.
  !***
  ! **Library:**   SLATEC (SLAP)
  !***
  ! **Category:**  D2E
  !***
  ! **Type:**      DOUBLE PRECISION (SSLLTI-S, DSLLTI-D)
  !***
  ! **Keywords:**  ITERATIVE PRECONDITION, LINEAR SYSTEM SOLVE, SLAP, SPARSE
  !***
  ! **Author:**  Greenbaum, Anne, (Courant Institute)
  !           Seager, Mark K., (LLNL)
  !             Lawrence Livermore National Laboratory
  !             PO BOX 808, L-60
  !             Livermore, CA 94550 (510) 423-3141
  !             seager@llnl.gov
  !***
  ! **Description:**
  !       It is assumed that RWORK and IWORK have initialized with
  !       the information required for DLLTI2:
  !          IWORK(1) = NEL
  !          IWORK(2) = Starting location of IEL in IWORK.
  !          IWORK(3) = Starting location of JEL in IWORK.
  !          IWORK(4) = Starting location of EL in RWORK.
  !          IWORK(5) = Starting location of DINV in RWORK.
  !       See the DESCRIPTION of DLLTI2 for details.
  !***
  ! **References:**  (NONE)
  !***
  ! **Routines called:**  DLLTI2

  !* REVISION HISTORY  (YYMMDD)
  !   871119  DATE WRITTEN
  !   881213  Previous REVISION DATE
  !   890915  Made changes requested at July 1989 CML Meeting.  (MKS)
  !   890922  Numerous changes to prologue to make closer to SLATEC
  !           standard.  (FNF)
  !   890929  Numerous changes to reduce SP/DP differences.  (FNF)
  !   910411  Prologue converted to Version 4.0 format.  (BAB)
  !   910502  Corrected conversion error.  (FNF)
  !   920511  Added complete declaration section.  (WRB)
  !   921113  Corrected C***CATEGORY line.  (FNF)
  !   930701  Updated CATEGORY section.  (FNF, WRB)
  
  !     .. Scalar Arguments ..
  INTEGER Isym, N, Nelt
  !     .. Array Arguments ..
  REAL(8) :: A(Nelt), B(*), Rwork(*), X(*)
  INTEGER Ia(Nelt), Iwork(*), Ja(Nelt)
  !     .. Local Scalars ..
  INTEGER locdin, locel, lociel, locjel, nel
  !     .. External Subroutines ..
  EXTERNAL :: DLLTI2
  !* FIRST EXECUTABLE STATEMENT  DSLLTI
  nel = Iwork(1)
  lociel = Iwork(3)
  locjel = Iwork(2)
  locel = Iwork(4)
  locdin = Iwork(5)
  CALL DLLTI2(N,B,X,nel,Iwork(lociel),Iwork(locjel),Rwork(locel),Rwork(locdin))
  !
  !------------- LAST LINE OF DSLLTI FOLLOWS ----------------------------
END SUBROUTINE DSLLTI
