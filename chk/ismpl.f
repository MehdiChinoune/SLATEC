*DECK ISMPL
      SUBROUTINE ISMPL (N, M, INDX)
C***BEGIN PROLOGUE  ISMPL
C***SUBSIDIARY
C***PURPOSE  Generate integer sample.
C            This routine picks M "random" integers in the range 1 to
C            N without any repetitions.
C***LIBRARY   SLATEC (SLAP)
C***TYPE      INTEGER (ISMPL-I)
C***AUTHOR  Seager, Mark K., (LLNL)
C             Lawrence Livermore National Laboratory
C             PO BOX 808, L-300
C             Livermore, CA 94550 (510) 423-3141
C             seager@llnl.gov
C***ROUTINES CALLED  RAND
C***REVISION HISTORY  (YYMMDD)
C   871119  DATE WRITTEN
C   881213  Previous REVISION DATE
C   890919  Changed to integer name ISMPL.  (MKS)
C   890920  Converted prologue to SLATEC 4.0 format.  (FNF)
C   920511  Added complete declaration section.  (WRB)
C***END PROLOGUE  ISMPL
C     .. Scalar Arguments ..
      INTEGER M, N
C     .. Array Arguments ..
      INTEGER INDX(M)
C     .. Local Scalars ..
      REAL DUMMY
      INTEGER I, ID, J
C     .. External Functions ..
      REAL RAND
      EXTERNAL RAND
C     .. Intrinsic Functions ..
      INTRINSIC INT
C***FIRST EXECUTABLE STATEMENT  ISMPL
C
C     Check the input
C
      DUMMY = 0.0
      IF( N*M.LT.0 .OR. M.GT.N ) RETURN
C
C     Set the indices.
      INDX(1) = INT( RAND(DUMMY)*N ) + 1
CVD$ NOCONCUR
      DO 30 I = 2, M
 10      ID = INT( RAND(DUMMY)*N ) + 1
C
C        Check to see if ID has already been chosen.
CVD$ NOVECTOR
CVD$ NOCONCUR
         DO 20 J = 1, I-1
            IF( ID.EQ.INDX(J) ) GOTO 10
 20      CONTINUE
         INDX(I) = ID
 30   CONTINUE
      RETURN
C------------- LAST LINE OF ISMPL FOLLOWS ------------------------------
      END
