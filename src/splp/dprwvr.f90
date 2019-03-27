!** DPRWVR
SUBROUTINE DPRWVR(Key,Ipage,Lpg,Sx,Ix)
  IMPLICIT NONE
  !>
  !***
  !  Subsidiary to DSPLP
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      DOUBLE PRECISION (PRWVIR-S, DPRWVR-D)
  !***
  ! **Author:**  Hanson, R. J., (SNLA)
  !           Wisniewski, J. A., (SNLA)
  !***
  ! **Description:**
  !
  !     DPRWVR LIMITS THE TYPE OF STORAGE TO A SEQUENTIAL SPARSE MATRIX
  !     STORAGE SCHEME.  THE PAGE STORAGE IS ON RANDOM ACCESS DISK.
  !     DPRWVR IS PART OF THE SPARSE LP PACKAGE, DSPLP.
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
  !***
  ! **See also:**  DSPLP
  !***
  ! **Routines called:**  DREADP, DWRITP, SOPENM

  !* REVISION HISTORY  (YYMMDD)
  !   811215  DATE WRITTEN
  !   891009  Removed unreferenced variables.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900328  Added TYPE section.  (WRB)
  !   910403  Updated AUTHOR and DESCRIPTION sections.  (WRB)

  INTEGER iaddr, Ipage, ipagef, istart, Ix(*), Key, Lpg
  LOGICAL :: first
  REAL(8) :: Sx(*), zero, one
  SAVE zero, one
  DATA zero, one/0.D0, 1.D0/
  !* FIRST EXECUTABLE STATEMENT  DPRWVR
  !
  !     COMPUTE STARTING ADDRESS OF PAGE.
  !
  ipagef = INT( Sx(3) )
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
    CALL DREADP(ipagef,Ix(istart),Sx(istart),Lpg,iaddr)
  ELSEIF ( Key==2 ) THEN
    CALL DWRITP(ipagef,Ix(istart),Sx(istart),Lpg,iaddr)
  ENDIF
END SUBROUTINE DPRWVR
