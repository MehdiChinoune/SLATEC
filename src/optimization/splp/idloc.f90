!** IDLOC
PURE INTEGER FUNCTION IDLOC(Locc,Ix)
  !> Subsidiary to DSPLP
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      DOUBLE PRECISION (IPLOC-S, IDLOC-D)
  !***
  ! **Keywords:**  RELATIVE ADDRESS DETERMINATION FUNCTION, SLATEC
  !***
  ! **Author:**  Boland, W. Robert, (LANL)
  !           Nicol, Tom, (University of British Columbia)
  !***
  ! **Description:**
  !
  !   Given a "virtual" location,  IDLOC returns the relative working
  !   address of the vector component stored in SX, IX.  Any necessary
  !   page swaps are performed automatically for the user in this
  !   function subprogram.
  !
  !   LOC       is the "virtual" address of the data to be retrieved.
  !   SX ,IX    represent the matrix where the data is stored.
  !
  !***
  ! **See also:**  DSPLP
  !***
  ! **Routines called:**  DPRWPG, XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   890606  DATE WRITTEN
  !   890606  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   910731  Added code to set IDLOC to 0 if LOC is non-positive.  (WRB)

  INTEGER, INTENT(IN) :: Locc
  INTEGER, INTENT(IN) :: Ix(:)
  INTEGER :: ipage, itemp, k, lmx, lmxm1, lpg
  !* FIRST EXECUTABLE STATEMENT  IDLOC
  IF( Locc<=0 ) THEN
    ERROR STOP 'IDLOC : A value of LOC, the first argument, <= 0 was encountered'
  END IF
  !
  !     Two cases exist:  (1<=LOC<=K) .OR. (LOC>K).
  !
  k = Ix(3) + 4
  lmx = Ix(1)
  lmxm1 = lmx - 1
  IF( Locc<=k ) THEN
    IDLOC = Locc
    RETURN
  END IF
  !
  !     Compute length of the page, starting address of the page, page
  !     number and relative working address.
  !
  lpg = lmx - k
  itemp = Locc - k - 1
  ipage = itemp/lpg + 1
  IDLOC = MOD(itemp,lpg) + k + 1

END FUNCTION IDLOC