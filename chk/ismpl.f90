!*==ISMPL.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK ISMPL
SUBROUTINE ISMPL(N,M,Indx)
  IMPLICIT NONE
  !*--ISMPL5
  !***BEGIN PROLOGUE  ISMPL
  !***SUBSIDIARY
  !***PURPOSE  Generate integer sample.
  !            This routine picks M "random" integers in the range 1 to
  !            N without any repetitions.
  !***LIBRARY   SLATEC (SLAP)
  !***TYPE      INTEGER (ISMPL-I)
  !***AUTHOR  Seager, Mark K., (LLNL)
  !             Lawrence Livermore National Laboratory
  !             PO BOX 808, L-300
  !             Livermore, CA 94550 (510) 423-3141
  !             seager@llnl.gov
  !***ROUTINES CALLED  RAND
  !***REVISION HISTORY  (YYMMDD)
  !   871119  DATE WRITTEN
  !   881213  Previous REVISION DATE
  !   890919  Changed to integer name ISMPL.  (MKS)
  !   890920  Converted prologue to SLATEC 4.0 format.  (FNF)
  !   920511  Added complete declaration section.  (WRB)
  !***END PROLOGUE  ISMPL
  !     .. Scalar Arguments ..
  INTEGER M, N
  !     .. Array Arguments ..
  INTEGER Indx(M)
  !     .. Local Scalars ..
  REAL dummy
  INTEGER i, id, j
  !     .. External Functions ..
  REAL RAND
  EXTERNAL RAND
  !     .. Intrinsic Functions ..
  INTRINSIC INT
  !***FIRST EXECUTABLE STATEMENT  ISMPL
  !
  !     Check the input
  !
  dummy = 0.0
  IF ( N*M<0.OR.M>N ) RETURN
  !
  !     Set the indices.
  Indx(1) = INT(RAND(dummy)*N) + 1
  !VD$ NOCONCUR
  DO i = 2, M
    DO
      id = INT(RAND(dummy)*N) + 1
      !
      !        Check to see if ID has already been chosen.
      !VD$ NOVECTOR
      !VD$ NOCONCUR
      DO j = 1, i - 1
        IF ( id==Indx(j) ) GOTO 50
      ENDDO
      Indx(i) = id
      EXIT
      50     ENDDO
    ENDDO
    !------------- LAST LINE OF ISMPL FOLLOWS ------------------------------
END SUBROUTINE ISMPL
