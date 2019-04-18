!** SSLI
SUBROUTINE SSLI(N,B,X,Nelt,Ia,Ja,A,Isym,Rwork,Iwork)
  !>
  !  SLAP MSOLVE for Lower Triangle Matrix.
  !            This routine acts as an interface between the SLAP generic
  !            MSOLVE calling convention and the routine that actually
  !                      -1
  !            computes L  B = X.
  !***
  ! **Library:**   SLATEC (SLAP)
  !***
  ! **Category:**  D2A3
  !***
  ! **Type:**      SINGLE PRECISION (SSLI-S, DSLI-D)
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
  !       the information required for SSLI2:
  !          IWORK(1) = NEL
  !          IWORK(2) = Starting location of IEL in IWORK.
  !          IWORK(3) = Starting location of JEL in IWORK.
  !          IWORK(4) = Starting location of EL in RWORK.
  !       See the DESCRIPTION of SSLI2 for details.
  !***
  ! **References:**  (NONE)
  !***
  ! **Routines called:**  SSLI2

  !* REVISION HISTORY  (YYMMDD)
  !   871119  DATE WRITTEN
  !   881213  Previous REVISION DATE
  !   890915  Made changes requested at July 1989 CML Meeting.  (MKS)
  !   890922  Numerous changes to prologue to make closer to SLATEC
  !           standard.  (FNF)
  !   890929  Numerous changes to reduce SP/DP differences.  (FNF)
  !   910411  Prologue converted to Version 4.0 format.  (BAB)
  !   920511  Added complete declaration section.  (WRB)
  !   921113  Corrected C***CATEGORY line.  (FNF)
  !   930701  Updated CATEGORY section.  (FNF, WRB)

  !     .. Scalar Arguments ..
  INTEGER Isym, N, Nelt
  !     .. Array Arguments ..
  REAL A(Nelt), B(N), Rwork(*), X(N)
  INTEGER Ia(Nelt), Iwork(*), Ja(Nelt)
  !     .. Local Scalars ..
  INTEGER locel, lociel, locjel, nel
  !* FIRST EXECUTABLE STATEMENT  SSLI
  !
  nel = Iwork(1)
  lociel = Iwork(2)
  locjel = Iwork(3)
  locel = Iwork(4)
  CALL SSLI2(N,B,X,nel,Iwork(lociel),Iwork(locjel),Rwork(locel))
  !
  !------------- LAST LINE OF SSLI FOLLOWS ----------------------------
END SUBROUTINE SSLI
