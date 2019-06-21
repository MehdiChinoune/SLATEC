!** ISMPL
SUBROUTINE ISMPL(N,M,Indx)
  !> Generate integer sample.
  !            This routine picks M "random" integers in the range 1 to
  !            N without any repetitions.
  !***
  ! **Library:**   SLATEC (SLAP)
  !***
  ! **Type:**      INTEGER (ISMPL-I)
  !***
  ! **Author:**  Seager, Mark K., (LLNL)
  !             Lawrence Livermore National Laboratory
  !             PO BOX 808, L-300
  !             Livermore, CA 94550 (510) 423-3141
  !             seager@llnl.gov
  !***
  ! **Routines called:**  RAND

  !* REVISION HISTORY  (YYMMDD)
  !   871119  DATE WRITTEN
  !   881213  Previous REVISION DATE
  !   890919  Changed to integer name ISMPL.  (MKS)
  !   890920  Converted prologue to SLATEC 4.0 format.  (FNF)
  !   920511  Added complete declaration section.  (WRB)
  USE slatec, ONLY : RAND
  !     .. Scalar Arguments ..
  INTEGER :: M, N
  !     .. Array Arguments ..
  INTEGER :: Indx(M)
  !     .. Local Scalars ..
  REAL(SP) :: dummy
  INTEGER :: i, id, j
  !     .. Intrinsic Functions ..
  INTRINSIC INT
  !* FIRST EXECUTABLE STATEMENT  ISMPL
  !
  !     Check the input
  !
  dummy = 0._SP
  IF( N*M<=0 .OR. M>N ) RETURN
  !
  !     Set the indices.
  Indx(1) = INT(RAND(dummy)*N) + 1
  DO i = 2, M
    DO
      id = INT(RAND(dummy)*N) + 1
      !
      !        Check to see if ID has already been chosen.
      DO j = 1, i - 1
        IF( id==Indx(j) ) GOTO 50
      END DO
      Indx(i) = id
      EXIT
      50 CONTINUE
    END DO
  END DO
  !------------- LAST LINE OF ISMPL FOLLOWS ------------------------------
END SUBROUTINE ISMPL
