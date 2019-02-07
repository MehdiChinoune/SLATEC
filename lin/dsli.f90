!*==DSLI.f90  processed by SPAG 6.72Dc at 10:58 on  6 Feb 2019
!DECK DSLI
SUBROUTINE DSLI(N,B,X,Nelt,Ia,Ja,A,Isym,Rwork,Iwork)
  IMPLICIT NONE
  !*--DSLI5
  !***BEGIN PROLOGUE  DSLI
  !***PURPOSE  SLAP MSOLVE for Lower Triangle Matrix.
  !            This routine acts as an interface between the SLAP generic
  !            MSOLVE calling convention and the routine that actually
  !                      -1
  !            computes L  B = X.
  !***LIBRARY   SLATEC (SLAP)
  !***CATEGORY  D2A3
  !***TYPE      DOUBLE PRECISION (SSLI-S, DSLI-D)
  !***KEYWORDS  ITERATIVE PRECONDITION, LINEAR SYSTEM SOLVE, SLAP, SPARSE
  !***AUTHOR  Greenbaum, Anne, (Courant Institute)
  !           Seager, Mark K., (LLNL)
  !             Lawrence Livermore National Laboratory
  !             PO BOX 808, L-60
  !             Livermore, CA 94550 (510) 423-3141
  !             seager@llnl.gov
  !***DESCRIPTION
  !       It is assumed that RWORK and IWORK have initialized with
  !       the information required for DSLI2:
  !          IWORK(1) = NEL
  !          IWORK(2) = Starting location of IEL in IWORK.
  !          IWORK(3) = Starting location of JEL in IWORK.
  !          IWORK(4) = Starting location of EL in RWORK.
  !       See the DESCRIPTION of DSLI2 for details.
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  DSLI2
  !***REVISION HISTORY  (YYMMDD)
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
  !***END PROLOGUE  DSLI
  !     .. Scalar Arguments ..
  INTEGER Isym , N , Nelt
  !     .. Array Arguments ..
  DOUBLE PRECISION A(Nelt) , B(N) , Rwork(*) , X(N)
  INTEGER Ia(Nelt) , Iwork(10) , Ja(Nelt)
  !     .. Local Scalars ..
  INTEGER locel , lociel , locjel , nel
  !     .. External Subroutines ..
  EXTERNAL DSLI2
  !***FIRST EXECUTABLE STATEMENT  DSLI
  !
  nel = Iwork(1)
  lociel = Iwork(2)
  locjel = Iwork(3)
  locel = Iwork(4)
  CALL DSLI2(N,B,X,nel,Iwork(lociel),Iwork(locjel),Rwork(locel))
  !
  !------------- LAST LINE OF DSLI FOLLOWS ----------------------------
END SUBROUTINE DSLI
