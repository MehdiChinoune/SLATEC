*DECK VFILL
      SUBROUTINE VFILL (N, V, VAL)
C***BEGIN PROLOGUE  VFILL
C***SUBSIDIARY
C***PURPOSE  Fill a vector with a value.
C***LIBRARY   SLATEC (SLAP)
C***TYPE      SINGLE PRECISION (VFILL-S, DFILL-D)
C***AUTHOR  Seager, Mark K., (LLNL)
C             Lawrence Livermore National Laboratory
C             PO BOX 808, L-300
C             Livermore, CA 94550 (510) 423-3141
C             seager@llnl.gov
C***DESCRIPTION
C
C *Usage:
C     INTEGER  N
C     REAL     V(N), VAL
C
C     CALL VFILL( N, V, VAL )
C
C *Arguments:
C N      :IN       Integer.
C         Length of the vector
C V      :OUT      Real V(N).
C         Vector to be set.
C VAL    :IN       Real.
C         Value to seed the vector with.
C***REFERENCES  (NONE)
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   871119  DATE WRITTEN
C   881213  Previous REVISION DATE
C   890920  Converted prologue to SLATEC 4.0 format.  (FNF)
C   920511  Added complete declaration section.  (WRB)
C***END PROLOGUE  VFILL
C     .. Scalar Arguments ..
      REAL VAL
      INTEGER N
C     .. Array Arguments ..
      REAL V(*)
C     .. Local Scalars ..
      INTEGER I, IS, NR
C     .. Intrinsic Functions ..
      INTRINSIC MOD
C***FIRST EXECUTABLE STATEMENT  VFILL
      IF (N .LE. 0) RETURN
      NR=MOD(N,4)
C
C         The following construct assumes a zero pass do loop.
C
      IS=1
      GOTO(1,2,3,4), NR+1
    4   IS=4
        V(1)=VAL
        V(2)=VAL
        V(3)=VAL
        GOTO 1
    3   IS=3
        V(1)=VAL
        V(2)=VAL
        GOTO 1
    2   IS=2
        V(1)=VAL
    1 DO 10 I=IS,N,4
        V(I)  =VAL
        V(I+1)=VAL
        V(I+2)=VAL
        V(I+3)=VAL
 10   CONTINUE
      RETURN
C------------- LAST LINE OF VFILL FOLLOWS -----------------------------
      END
