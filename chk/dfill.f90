!DECK DFILL
SUBROUTINE DFILL(N,V,Val)
  IMPLICIT NONE
  !***BEGIN PROLOGUE  DFILL
  !***SUBSIDIARY
  !***PURPOSE  Fill a vector with a value.
  !***LIBRARY   SLATEC (SLAP)
  !***TYPE      DOUBLE PRECISION (VFILL-S, DFILL-D)
  !***AUTHOR  Seager, Mark K., (LLNL)
  !             Lawrence Livermore National Laboratory
  !             PO BOX 808, L-300
  !             Livermore, CA 94550 (510) 423-3141
  !             seager@llnl.gov
  !***DESCRIPTION
  !
  ! *Usage:
  !     INTEGER  N
  !     DOUBLE PRECISION V(N), VAL
  !
  !     CALL DFILL( N, V, VAL )
  !
  ! *Arguments:
  ! N      :IN       Integer.
  !         Length of the vector
  ! V      :OUT      Double Precision V(N).
  !         Vector to be set.
  ! VAL    :IN       Double Precision.
  !         Value to seed the vector with.
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   890404  DATE WRITTEN
  !   890920  Converted prologue to SLATEC 4.0 format.  (FNF)
  !   920511  Added complete declaration section.  (WRB)
  !***END PROLOGUE  DFILL
  !     .. Scalar Arguments ..
  REAL(8) :: Val
  INTEGER N
  !     .. Array Arguments ..
  REAL(8) :: V(*)
  !     .. Local Scalars ..
  INTEGER i, is, nr
  !     .. Intrinsic Functions ..
  INTRINSIC MOD
  !***FIRST EXECUTABLE STATEMENT  DFILL
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
  !------------- LAST LINE OF DFILL FOLLOWS -----------------------------
END SUBROUTINE DFILL
