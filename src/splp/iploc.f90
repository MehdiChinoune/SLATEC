!** IPLOC
INTEGER FUNCTION IPLOC(Loc,Sx,Ix)
  IMPLICIT NONE
  !>
  !***
  !  Subsidiary to SPLP
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      SINGLE PRECISION (IPLOC-S, IDLOC-D)
  !***
  ! **Keywords:**  RELATIVE ADDRESS DETERMINATION FUNCTION, SLATEC
  !***
  ! **Author:**  Hanson, R. J., (SNLA)
  !           Wisniewski, J. A., (SNLA)
  !***
  ! **Description:**
  !
  !   Given a "virtual" location,  IPLOC returns the relative working
  !   address of the vector component stored in SX, IX.  Any necessary
  !   page swaps are performed automatically for the user in this
  !   function subprogram.
  !
  !   LOC       is the "virtual" address of the data to be retrieved.
  !   SX ,IX    represent the matrix where the data is stored.
  !
  !***
  ! **See also:**  SPLP
  !***
  ! **Routines called:**  PRWPGE, XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   810306  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890606  Restructured to match double precision version.  (WRB)
  !   890606  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   910731  Added code to set IPLOC to 0 if LOC is non-positive.  (WRB)
  
  INTEGER ipage, itemp, k, key, lmx, lmxm1, Loc, lpg, np
  REAL Sx(*)
  INTEGER Ix(*)
  !* FIRST EXECUTABLE STATEMENT  IPLOC
  IF ( Loc<=0 ) THEN
    CALL XERMSG('SLATEC','IPLOC',&
      'A value of LOC, the first argument, .LE. 0 was encountered',55,1)
    IPLOC = 0
    RETURN
  ENDIF
  !
  !     Two cases exist:  (1.LE.LOC.LE.K) .OR. (LOC.GT.K).
  !
  k = Ix(3) + 4
  lmx = Ix(1)
  lmxm1 = lmx - 1
  IF ( Loc<=k ) THEN
    IPLOC = Loc
    RETURN
  ENDIF
  !
  !     Compute length of the page, starting address of the page, page
  !     number and relative working address.
  !
  lpg = lmx - k
  itemp = Loc - k - 1
  ipage = itemp/lpg + 1
  IPLOC = MOD(itemp,lpg) + k + 1
  np = ABS(Ix(lmxm1))
  !
  !     Determine if a page fault has occurred.  If so, write page NP
  !     and read page IPAGE.  Write the page only if it has been
  !     modified.
  !
  IF ( ipage/=np ) THEN
    IF ( Sx(lmx)==1.0 ) THEN
      Sx(lmx) = 0.0
      key = 2
      CALL PRWPGE(key,np,lpg,Sx,Ix)
    ENDIF
    key = 1
    CALL PRWPGE(key,ipage,lpg,Sx,Ix)
  ENDIF
END FUNCTION IPLOC
