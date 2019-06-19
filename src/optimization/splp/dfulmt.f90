!** DFULMT
SUBROUTINE DFULMT(I,J,Aij,Indcat,Prgopt,Dattrv,Iflag)
  !> Subsidiary to DSPLP
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      DOUBLE PRECISION (FULMAT-S, DFULMT-D)
  !***
  ! **Author:**  (UNKNOWN)
  !***
  ! **Description:**
  !
  !     DECODES A STANDARD TWO-DIMENSIONAL FORTRAN ARRAY PASSED
  !     IN THE ARRAY DATTRV(IA,*).  THE ROW DIMENSION IA AND THE
  !     MATRIX DIMENSIONS MRELAS AND NVARS MUST SIMULTANEOUSLY BE
  !     PASSED USING THE OPTION ARRAY, PRGOPT(*).  IT IS AN ERROR
  !     IF THIS DATA IS NOT PASSED TO DFULMT( ).
  !     EXAMPLE-- (FOR USE TOGETHER WITH DSPLP().)
  !      EXTERNAL DUSRMT
  !      DIMENSION DATTRV(IA,*)
  !      PRGOPT(01)=7
  !      PRGOPT(02)=68
  !      PRGOPT(03)=1
  !      PRGOPT(04)=IA
  !      PRGOPT(05)=MRELAS
  !      PRGOPT(06)=NVARS
  !      PRGOPT(07)=1
  !     CALL DSPLP(  ... DFULMT INSTEAD OF DUSRMT...)
  !
  !***
  ! **See also:**  DSPLP
  !***
  ! **Routines called:**  XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   811215  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900328  Added TYPE section.  (WRB)
  USE service, ONLY : XERMSG
  INTEGER :: I, Indcat, J
  INTEGER :: Iflag(10)
  REAL(DP) :: Aij, Dattrv(:), Prgopt(:)
  INTEGER :: key, level, lp, nerr, next
  REAL(DP), PARAMETER :: zero = 0.D0
  !* FIRST EXECUTABLE STATEMENT  DFULMT
  IF( Iflag(1)==1 ) THEN
    !     INITIALIZE POINTERS TO PROCESS FULL TWO-DIMENSIONAL FORTRAN
    !     ARRAYS.
    lp = 1
    DO
      next = INT( Prgopt(lp) )
      IF( next>1 ) THEN
        key = INT( Prgopt(lp+1) )
        IF( key/=68 ) THEN
          lp = next
        ELSEIF( Prgopt(lp+2)/=zero ) THEN
          Iflag(2) = 1
          Iflag(3) = 1
          Iflag(4) = INT( Prgopt(lp+3) )
          Iflag(5) = INT( Prgopt(lp+4) )
          Iflag(6) = INT( Prgopt(lp+5) )
          EXIT
        ELSE
          lp = next
        END IF
      ELSE
        nerr = 29
        level = 1
        CALL XERMSG('DFULMT',&
          'IN DSPLP, ROW DIM., MRELAS, NVARS ARE MISSING FROM PRGOPT.',nerr,level)
        Iflag(1) = 3
        EXIT
      END IF
    END DO
  ELSEIF( Iflag(1)==2 ) THEN
    DO
      I = Iflag(2)
      J = Iflag(3)
      IF( J>Iflag(6) ) THEN
        Iflag(1) = 3
        EXIT
      ELSEIF( I<=Iflag(5) ) THEN
        Aij = Dattrv(Iflag(4)*(J-1)+I)
        Iflag(2) = I + 1
        IF( Aij/=zero ) THEN
          Indcat = 0
          EXIT
        END IF
      ELSE
        Iflag(2) = 1
        Iflag(3) = J + 1
      END IF
    END DO
  END IF
END SUBROUTINE DFULMT
