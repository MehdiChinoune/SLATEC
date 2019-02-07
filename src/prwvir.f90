!*==PRWVIR.f90  processed by SPAG 6.72Dc at 11:01 on  6 Feb 2019
!DECK PRWVIR
SUBROUTINE PRWVIR(Key,Ipage,Lpg,Sx,Ix)
  IMPLICIT NONE
  !*--PRWVIR5
  !*** Start of declarations inserted by SPAG
  INTEGER iaddr , Ipage , ipagef , istart , Ix , Key , Lpg
  !*** End of declarations inserted by SPAG
  !***BEGIN PROLOGUE  PRWVIR
  !***SUBSIDIARY
  !***PURPOSE  Subsidiary to SPLP
  !***LIBRARY   SLATEC
  !***TYPE      SINGLE PRECISION (PRWVIR-S, DPRWVR-D)
  !***AUTHOR  Hanson, R. J., (SNLA)
  !           Wisniewski, J. A., (SNLA)
  !***DESCRIPTION
  !
  !     PRWVIR LIMITS THE TYPE OF STORAGE TO A SEQUENTIAL SPARSE MATRIX
  !     STORAGE SCHEME.  THE PAGE STORAGE IS ON RANDOM ACCESS DISK.
  !     PRWVIR IS PART OF THE SPARSE LP PACKAGE, SPLP.
  !
  !     KEY       IS A FLAG WHICH INDICATES WHETHER A READ OR WRITE
  !               OPERATION IS TO BE PERFORMED. A VALUE OF KEY=1 INDICATES
  !               A READ. A VALUE OF KEY=2 INDICATES A WRITE.
  !     IPAGE     IS THE PAGE OF MATRIX MN WE ARE ACCESSING.
  !     LPG       IS THE LENGTH OF THE PAGE.
  !   SX(*),IX(*) IS THE MATRIX DATA.
  !
  !     THIS SUBROUTINE IS A MODIFICATION OF THE SUBROUTINE LRWVIR,
  !     SANDIA LABS. REPT. SAND78-0785.
  !     MODIFICATIONS BY K.L. HIEBERT AND R.J. HANSON
  !
  !***SEE ALSO  SPLP
  !***ROUTINES CALLED  SOPENM, SREADP, SWRITP
  !***REVISION HISTORY  (YYMMDD)
  !   811215  DATE WRITTEN
  !   891009  Removed unreferenced variables.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900328  Added TYPE section.  (WRB)
  !   910403  Updated AUTHOR and DESCRIPTION sections.  (WRB)
  !***END PROLOGUE  PRWVIR
  DIMENSION Ix(*)
  REAL Sx(*) , zero , one
  LOGICAL first
  SAVE zero , one
  DATA zero , one/0.E0 , 1.E0/
  !***FIRST EXECUTABLE STATEMENT  PRWVIR
  !
  !     COMPUTE STARTING ADDRESS OF PAGE.
  !
  ipagef = Sx(3)
  istart = Ix(3) + 5
  !
  !     OPEN RANDOM ACCESS FILE NUMBER IPAGEF, IF FIRST PAGE WRITE.
  !
  first = Sx(4)==zero
  IF ( first ) THEN
    CALL SOPENM(ipagef,Lpg)
    Sx(4) = one
  ENDIF
  !
  !     PERFORM EITHER A READ OR A WRITE.
  !
  iaddr = 2*Ipage - 1
  IF ( Key==1 ) THEN
    CALL SREADP(ipagef,Ix(istart),Sx(istart),Lpg,iaddr)
  ELSEIF ( Key==2 ) THEN
    CALL SWRITP(ipagef,Ix(istart),Sx(istart),Lpg,iaddr)
  ENDIF
END SUBROUTINE PRWVIR
