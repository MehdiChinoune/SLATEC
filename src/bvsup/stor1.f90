!DECK STOR1
SUBROUTINE STOR1(U,Yh,V,Yp,Ntemp,Ndisk,Ntape)
  IMPLICIT NONE
  !***BEGIN PROLOGUE  STOR1
  !***SUBSIDIARY
  !***PURPOSE  Subsidiary to BVSUP
  !***LIBRARY   SLATEC
  !***TYPE      SINGLE PRECISION (STOR1-S, DSTOR1-D)
  !***AUTHOR  Watts, H. A., (SNLA)
  !***DESCRIPTION
  !
  ! **********************************************************************
  !             0 -- Storage at output points.
  !     NTEMP =
  !             1 -- Temporary storage
  ! **********************************************************************
  !
  !***SEE ALSO  BVSUP
  !***ROUTINES CALLED  (NONE)
  !***COMMON BLOCKS    ML8SZ
  !***REVISION HISTORY  (YYMMDD)
  !   750601  DATE WRITTEN
  !   890921  Realigned order of variables in certain COMMON blocks.
  !           (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900328  Added TYPE section.  (WRB)
  !   910722  Updated AUTHOR section.  (ALS)
  !***END PROLOGUE  STOR1
  REAL C, U, V, XSAv, Yh, Yp
  INTEGER IGOfx, INHomo, IVP, j, NCOmp, nctnf, Ndisk, NFC, Ntape, Ntemp
  DIMENSION U(*), Yh(*), V(*), Yp(*)
  !
  ! **********************************************************************
  !
  COMMON /ML8SZ / C, XSAv, IGOfx, INHomo, IVP, NCOmp, NFC
  !
  ! **********************************************************************
  !
  !***FIRST EXECUTABLE STATEMENT  STOR1
  nctnf = NCOmp*NFC
  DO j = 1, nctnf
    U(j) = Yh(j)
  ENDDO
  IF ( INHomo/=1 ) THEN
    !
    !   ZERO PARTICULAR SOLUTION
    !
    IF ( Ntemp==1 ) RETURN
    DO j = 1, NCOmp
      V(j) = 0.
    ENDDO
    IF ( Ndisk==1 ) WRITE (Ntape) (V(j),j=1,NCOmp), (U(j),j=1,nctnf)
    !
    !   NONZERO PARTICULAR SOLUTION
    !
  ELSEIF ( Ntemp==0 ) THEN
    !
    DO j = 1, NCOmp
      V(j) = C*Yp(j)
    ENDDO
    !
    !  IS OUTPUT INFORMATION TO BE WRITTEN TO DISK
    !
    IF ( Ndisk==1 ) WRITE (Ntape) (V(j),j=1,NCOmp), (U(j),j=1,nctnf)
  ELSE
    !
    DO j = 1, NCOmp
      V(j) = Yp(j)
    ENDDO
    RETURN
  ENDIF
  !
END SUBROUTINE STOR1
