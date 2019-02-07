!*==SSLUTI.f90  processed by SPAG 6.72Dc at 10:58 on  6 Feb 2019
!DECK SSLUTI
SUBROUTINE SSLUTI(N,B,X,Nelt,Ia,Ja,A,Isym,Rwork,Iwork)
  IMPLICIT NONE
  !*--SSLUTI5
  !***BEGIN PROLOGUE  SSLUTI
  !***PURPOSE  SLAP MTSOLV for LDU Factorization.
  !            This routine acts as an interface between the SLAP generic
  !            MTSOLV calling convention and the routine that actually
  !                           -T
  !            computes  (LDU)  B = X.
  !***LIBRARY   SLATEC (SLAP)
  !***CATEGORY  D2E
  !***TYPE      SINGLE PRECISION (SSLUTI-S, DSLUTI-D)
  !***KEYWORDS  ITERATIVE PRECONDITION, LINEAR SYSTEM SOLVE, SLAP, SPARSE
  !***AUTHOR  Greenbaum, Anne, (Courant Institute)
  !           Seager, Mark K., (LLNL)
  !             Lawrence Livermore National Laboratory
  !             PO BOX 808, L-60
  !             Livermore, CA 94550 (510) 423-3141
  !             seager@llnl.gov
  !***DESCRIPTION
  !       It is assumed that RWORK and IWORK have initialized with
  !       the information required for SSLUI4:
  !          IWORK(1) = Starting location of IL in IWORK.
  !          IWORK(2) = Starting location of JL in IWORK.
  !          IWORK(3) = Starting location of IU in IWORK.
  !          IWORK(4) = Starting location of JU in IWORK.
  !          IWORK(5) = Starting location of L in RWORK.
  !          IWORK(6) = Starting location of DINV in RWORK.
  !          IWORK(7) = Starting location of U in RWORK.
  !       See the DESCRIPTION of SSLUI4 for details.
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  SSLUI4
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
  !***END PROLOGUE  SSLUTI
  !     .. Scalar Arguments ..
  INTEGER Isym , N , Nelt
  !     .. Array Arguments ..
  REAL A(N) , B(N) , Rwork(*) , X(N)
  INTEGER Ia(Nelt) , Iwork(10) , Ja(Nelt)
  !     .. Local Scalars ..
  INTEGER locdin , locil , lociu , locjl , locju , locl , locu
  !     .. External Subroutines ..
  EXTERNAL SSLUI4
  !***FIRST EXECUTABLE STATEMENT  SSLUTI
  !
  !         Pull out the pointers to the L, D and U matrices and call
  !         the workhorse routine.
  !
  locil = Iwork(1)
  locjl = Iwork(2)
  lociu = Iwork(3)
  locju = Iwork(4)
  locl = Iwork(5)
  locdin = Iwork(6)
  locu = Iwork(7)
  !
  CALL SSLUI4(N,B,X,Iwork(locil),Iwork(locjl),Rwork(locl),Rwork(locdin),&
    Iwork(lociu),Iwork(locju),Rwork(locu))
  !
  !------------- LAST LINE OF SSLUTI FOLLOWS ----------------------------
END SUBROUTINE SSLUTI
