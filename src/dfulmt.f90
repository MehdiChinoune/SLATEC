!*==DFULMT.f90  processed by SPAG 6.72Dc at 11:01 on  6 Feb 2019
!DECK DFULMT
SUBROUTINE DFULMT(I,J,Aij,Indcat,Prgopt,Dattrv,Iflag)
  IMPLICIT NONE
  !*--DFULMT5
  !*** Start of declarations inserted by SPAG
  INTEGER I , Indcat , J , key , level , lp , nerr , next
  !*** End of declarations inserted by SPAG
  !***BEGIN PROLOGUE  DFULMT
  !***SUBSIDIARY
  !***PURPOSE  Subsidiary to DSPLP
  !***LIBRARY   SLATEC
  !***TYPE      DOUBLE PRECISION (FULMAT-S, DFULMT-D)
  !***AUTHOR  (UNKNOWN)
  !***DESCRIPTION
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
  !***SEE ALSO  DSPLP
  !***ROUTINES CALLED  XERMSG
  !***REVISION HISTORY  (YYMMDD)
  !   811215  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900328  Added TYPE section.  (WRB)
  !***END PROLOGUE  DFULMT
  DOUBLE PRECISION Aij , zero , Dattrv(*) , Prgopt(*)
  INTEGER Iflag(10)
  SAVE zero
  !***FIRST EXECUTABLE STATEMENT  DFULMT
  IF ( Iflag(1)==1 ) THEN
    !     INITIALIZE POINTERS TO PROCESS FULL TWO-DIMENSIONAL FORTRAN
    !     ARRAYS.
    zero = 0.D0
    lp = 1
    DO
      next = Prgopt(lp)
      IF ( next>1 ) THEN
        key = Prgopt(lp+1)
        IF ( key/=68 ) THEN
          lp = next
        ELSEIF ( Prgopt(lp+2)/=zero ) THEN
          Iflag(2) = 1
          Iflag(3) = 1
          Iflag(4) = Prgopt(lp+3)
          Iflag(5) = Prgopt(lp+4)
          Iflag(6) = Prgopt(lp+5)
          EXIT
        ELSE
          lp = next
        ENDIF
      ELSE
        nerr = 29
        level = 1
        CALL XERMSG('SLATEC','DFULMT',&
          'IN DSPLP, ROW DIM., MRELAS, NVARS ARE MISSING FROM PRGOPT.'&
          ,nerr,level)
        Iflag(1) = 3
        EXIT
      ENDIF
    ENDDO
  ELSEIF ( Iflag(1)==2 ) THEN
    DO
      I = Iflag(2)
      J = Iflag(3)
      IF ( J>Iflag(6) ) THEN
        Iflag(1) = 3
        EXIT
      ELSEIF ( I<=Iflag(5) ) THEN
        Aij = Dattrv(Iflag(4)*(J-1)+I)
        Iflag(2) = I + 1
        IF ( Aij/=zero ) THEN
          Indcat = 0
          EXIT
        ENDIF
      ELSE
        Iflag(2) = 1
        Iflag(3) = J + 1
      ENDIF
    ENDDO
  ENDIF
END SUBROUTINE DFULMT
