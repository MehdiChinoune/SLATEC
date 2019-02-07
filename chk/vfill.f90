!*==VFILL.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK VFILL
SUBROUTINE VFILL(N,V,Val)
  IMPLICIT NONE
  !*--VFILL5
  !***BEGIN PROLOGUE  VFILL
  !***SUBSIDIARY
  !***PURPOSE  Fill a vector with a value.
  !***LIBRARY   SLATEC (SLAP)
  !***TYPE      SINGLE PRECISION (VFILL-S, DFILL-D)
  !***AUTHOR  Seager, Mark K., (LLNL)
  !             Lawrence Livermore National Laboratory
  !             PO BOX 808, L-300
  !             Livermore, CA 94550 (510) 423-3141
  !             seager@llnl.gov
  !***DESCRIPTION
  !
  ! *Usage:
  !     INTEGER  N
  !     REAL     V(N), VAL
  !
  !     CALL VFILL( N, V, VAL )
  !
  ! *Arguments:
  ! N      :IN       Integer.
  !         Length of the vector
  ! V      :OUT      Real V(N).
  !         Vector to be set.
  ! VAL    :IN       Real.
  !         Value to seed the vector with.
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   871119  DATE WRITTEN
  !   881213  Previous REVISION DATE
  !   890920  Converted prologue to SLATEC 4.0 format.  (FNF)
  !   920511  Added complete declaration section.  (WRB)
  !***END PROLOGUE  VFILL
  !     .. Scalar Arguments ..
  REAL Val
  INTEGER N
  !     .. Array Arguments ..
  REAL V(*)
  !     .. Local Scalars ..
  INTEGER i, is, nr
  !     .. Intrinsic Functions ..
  INTRINSIC MOD
  !***FIRST EXECUTABLE STATEMENT  VFILL
  IF ( N<=0 ) RETURN
  nr = MOD(N,4)
  !
  !         The following construct assumes a zero pass do loop.
  !
  is = 1
  SELECT CASE (nr+1)
    CASE (1)
    CASE (2)
      is = 2
      V(1) = Val
    CASE (3)
      is = 3
      V(1) = Val
      V(2) = Val
    CASE DEFAULT
      is = 4
      V(1) = Val
      V(2) = Val
      V(3) = Val
  END SELECT
  DO i = is, N, 4
    V(i) = Val
    V(i+1) = Val
    V(i+2) = Val
    V(i+3) = Val
  ENDDO
  !------------- LAST LINE OF VFILL FOLLOWS -----------------------------
END SUBROUTINE VFILL
