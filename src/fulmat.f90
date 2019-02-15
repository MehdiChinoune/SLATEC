!DECK FULMAT
SUBROUTINE FULMAT(I,J,Aij,Indcat,Prgopt,Dattrv,Iflag)
  IMPLICIT NONE
  INTEGER I, Indcat, J, key, level, lp, nerr, next
  !***BEGIN PROLOGUE  FULMAT
  !***SUBSIDIARY
  !***PURPOSE  Subsidiary to SPLP
  !***LIBRARY   SLATEC
  !***TYPE      SINGLE PRECISION (FULMAT-S, DFULMT-D)
  !***AUTHOR  (UNKNOWN)
  !***DESCRIPTION
  !
  !     DECODES A STANDARD TWO-DIMENSIONAL FORTRAN ARRAY PASSED
  !     IN THE ARRAY DATTRV(IA,*).  THE ROW DIMENSION IA AND THE
  !     MATRIX DIMENSIONS MRELAS AND NVARS MUST SIMULTANEOUSLY BE
  !     PASSED USING THE OPTION ARRAY, PRGOPT(*).  IT IS AN ERROR
  !     IF THIS DATA IS NOT PASSED TO FULMAT( ).
  !     EXAMPLE-- (FOR USE TOGETHER WITH SPLP().)
  !      EXTERNAL USRMAT
  !      DIMENSION DATTRV(IA,*)
  !      PRGOPT(01)=7
  !      PRGOPT(02)=68
  !      PRGOPT(03)=1
  !      PRGOPT(04)=IA
  !      PRGOPT(05)=MRELAS
  !      PRGOPT(06)=NVARS
  !      PRGOPT(07)=1
  !     CALL SPLP(  ... FULMAT INSTEAD OF USRMAT...)
  !
  !***SEE ALSO  SPLP
  !***ROUTINES CALLED  XERMSG
  !***REVISION HISTORY  (YYMMDD)
  !   811215  DATE WRITTEN
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900328  Added TYPE section.  (WRB)
  !***END PROLOGUE  FULMAT
  REAL Aij, zero, Dattrv(*), Prgopt(*)
  INTEGER Iflag(10)
  SAVE zero
  !***FIRST EXECUTABLE STATEMENT  FULMAT
  IF ( Iflag(1)==1 ) THEN
    !     INITIALIZE POINTERS TO PROCESS FULL TWO-DIMENSIONAL FORTRAN
    !     ARRAYS.
    zero = 0.
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
        CALL XERMSG('SLATEC','FULMAT',&
          'IN SPLP PACKAGE, ROW DIM., MRELAS, NVARS ARE MISSING FROM '&
          //'PRGOPT.',nerr,level)
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
END SUBROUTINE FULMAT
