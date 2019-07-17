!** IPLOC
PURE INTEGER FUNCTION IPLOC(Locc,Ix)
  !> Subsidiary to SPLP
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

  INTEGER, INTENT(IN) :: Locc
  INTEGER, INTENT(IN) :: Ix(:)
  INTEGER :: ipage, itemp, k, lmx, lmxm1, lpg
  !* FIRST EXECUTABLE STATEMENT  IPLOC
  IF( Locc<=0 ) THEN
    ERROR STOP 'IPLO : A value of LOC, the first argument, <= 0 was encountered'
    IPLOC = 0
    RETURN
  END IF
  !
  !     Two cases exist:  (1<=LOC<=K) .OR. (LOC>K).
  !
  k = Ix(3) + 4
  lmx = Ix(1)
  lmxm1 = lmx - 1
  IF( Locc<=k ) THEN
    IPLOC = Locc
    RETURN
  END IF
  !
  !     Compute length of the page, starting address of the page, page
  !     number and relative working address.
  !
  lpg = lmx - k
  itemp = Locc - k - 1
  ipage = itemp/lpg + 1
  IPLOC = MOD(itemp,lpg) + k + 1

END FUNCTION IPLOC