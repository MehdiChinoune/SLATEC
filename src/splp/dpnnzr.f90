!DECK DPNNZR
SUBROUTINE DPNNZR(I,Xval,Iplace,Sx,Ix,Ircx)
  IMPLICIT NONE
  !***BEGIN PROLOGUE  DPNNZR
  !***SUBSIDIARY
  !***PURPOSE  Subsidiary to DSPLP
  !***LIBRARY   SLATEC
  !***TYPE      DOUBLE PRECISION (PNNZRS-S, DPNNZR-D)
  !***AUTHOR  Hanson, R. J., (SNLA)
  !           Wisniewski, J. A., (SNLA)
  !***DESCRIPTION
  !
  !     DPNNZR LIMITS THE TYPE OF STORAGE TO A SEQUENTIAL SCHEME.
  !     SPARSE MATRIX NON ZERO RETRIEVAL SUBROUTINE.
  !
  !     SUBROUTINE DPNNZR() GETS THE NEXT NONZERO VALUE IN ROW OR COLUMN
  !     +/- IRCX WITH AN INDEX GREATER THAN THE VALUE OF I.
  !
  !             I ABSOLUTE VALUE OF THIS SUBSCRIPT IS TO BE EXCEEDED
  !               IN THE SEARCH FOR THE NEXT NONZERO VALUE. A NEGATIVE
  !               OR ZERO VALUE OF I CAUSES THE SEARCH TO START AT
  !               THE BEGINNING OF THE VECTOR.  A POSITIVE VALUE
  !               OF I CAUSES THE SEARCH TO CONTINUE FROM THE LAST PLACE
  !               ACCESSED. ON OUTPUT, THE ARGUMENT I
  !               CONTAINS THE VALUE OF THE SUBSCRIPT FOUND.  AN OUTPUT
  !               VALUE OF I EQUAL TO ZERO INDICATES THAT ALL COMPONENTS
  !               WITH AN INDEX GREATER THAN THE INPUT VALUE OF I ARE
  !               ZERO.
  !          XVAL VALUE OF THE NONZERO ELEMENT FOUND.  ON OUTPUT,
  !               XVAL=0. WHENEVER I=0.
  !     IPLACE POINTER INFORMATION WHICH IS MAINTAINED BY THE PACKAGE.
  !   SX(*),IX(*) THE WORK ARRAYS WHICH ARE USED TO STORE THE SPARSE
  !               MATRIX.  THESE ARRAY CONTENTS ARE AUTOMATICALLY
  !               MAINTAINED BY THE PACKAGE FOR THE USER.
  !          IRCX POINTS TO THE VECTOR OF THE MATRIX BEING SCANNED.  A
  !               NEGATIVE VALUE OF IRCX INDICATES THAT ROW -IRCX IS TO BE
  !               SCANNED.  A POSITIVE VALUE OF IRCX INDICATES THAT
  !               COLUMN IRCX IS TO BE SCANNED.  A ZERO VALUE OF IRCX IS
  !               AN ERROR.
  !
  !     THIS SUBROUTINE IS A MODIFICATION OF THE SUBROUTINE LNNZRS,
  !     SANDIA LABS. REPT. SAND78-0785.
  !     MODIFICATIONS BY K.L. HIEBERT AND R.J. HANSON
  !     REVISED 811130-1000
  !     REVISED YYMMDD-HHMM
  !
  !***SEE ALSO  DSPLP
  !***ROUTINES CALLED  IDLOC, XERMSG
  !***REVISION HISTORY  (YYMMDD)
  !   811215  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890605  Removed unreferenced labels.  (WRB)
  !   890606  Changed references from IPLOC to IDLOC.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900328  Added TYPE section.  (WRB)
  !   910403  Updated AUTHOR and DESCRIPTION sections.  (WRB)
  !***END PROLOGUE  DPNNZR
  INTEGER I, i1, idiff, IDLOC, iend, ii, il, ilast, iopt, ipl, &
    Iplace, ipploc, Ircx, istart, Ix, j, l, ll, lmx, lpg
  INTEGER n20046, nerr, np
  DIMENSION Ix(*)
  REAL(8) :: Xval, Sx(*), zero
  SAVE zero
  DATA zero/0.D0/
  !***FIRST EXECUTABLE STATEMENT  DPNNZR
  iopt = 1
  !
  !     CHECK VALIDITY OF ROW/COL. INDEX.
  !
  IF ( Ircx==0 ) THEN
    nerr = 55
    CALL XERMSG('SLATEC','DPNNZR','IRCX=0',nerr,iopt)
  ENDIF
  !
  !     LMX IS THE LENGTH OF THE IN-MEMORY STORAGE AREA.
  !
  lmx = Ix(1)
  IF ( Ircx>=0 ) THEN
    !
    !     CHECK SUBSCRIPTS OF THE COLUMN. THE COL. NUMBER MUST BE .LE. N AND
    !     THE INDEX MUST BE .LE. M.
    !
    IF ( Ircx>Ix(3).OR.ABS(I)>Ix(2) ) THEN
      nerr = 55
      CALL XERMSG('SLATEC','DPNNZR',&
        'SUBSCRIPTS FOR ARRAY ELEMENT TO BE ACCESSED WERE OUT OF BOUNDS',nerr,iopt)
    ENDIF
    l = Ix(2)
  ELSE
    !
    !     CHECK SUBSCRIPTS OF THE ROW. THE ROW NUMBER MUST BE .LE. M AND
    !     THE INDEX MUST BE .LE. N.
    !
    IF ( Ix(2)<-Ircx.OR.Ix(3)<ABS(I) ) THEN
      nerr = 55
      CALL XERMSG('SLATEC','DPNNZR',&
        'SUBSCRIPTS FOR ARRAY ELEMENT TO BE ACCESSED WERE OUT OF BOUNDS.',nerr,iopt)
    ENDIF
    l = Ix(3)
  ENDIF
  !
  !     HERE L IS THE LARGEST POSSIBLE SUBSCRIPT WITHIN THE VECTOR.
  !
  j = ABS(Ircx)
  ll = Ix(3) + 4
  lpg = lmx - ll
  IF ( Ircx<=0 ) THEN
    !
    !     SEARCH A ROW FOR THE NEXT NONZERO.
    !     FIND ELEMENT J=ABS(IRCX) IN ROWS ABS(I)+1,...,L.
    !
    I = ABS(I)
    !
    !     CHECK FOR END OF VECTOR.
    !
    IF ( I/=l ) THEN
      i1 = I + 1
      ii = i1
      n20046 = l
      GOTO 200
    ELSE
      I = 0
      Xval = zero
      RETURN
    ENDIF
  ELSE
    !
    !     SEARCHING FOR THE NEXT NONZERO IN A COLUMN.
    !
    !     INITIALIZE STARTING LOCATIONS..
    IF ( I<=0 ) THEN
      IF ( j/=1 ) THEN
        Iplace = Ix(j+3) + 1
      ELSE
        Iplace = ll + 1
      ENDIF
    ENDIF
    !
    !     THE CASE I.LE.0 SIGNALS THAT THE SCAN FOR THE ENTRY
    !     IS TO BEGIN AT THE START OF THE VECTOR.
    !
    I = ABS(I)
    IF ( j/=1 ) THEN
      istart = Ix(j+3) + 1
    ELSE
      istart = ll + 1
    ENDIF
    iend = Ix(j+4)
    !
    !     VALIDATE IPLACE. SET TO START OF VECTOR IF OUT OF RANGE.
    !
    IF ( istart>Iplace.OR.Iplace>iend ) THEN
      IF ( j/=1 ) THEN
        Iplace = Ix(j+3) + 1
      ELSE
        Iplace = ll + 1
      ENDIF
    ENDIF
    !
    !     SCAN THROUGH SEVERAL PAGES, IF NECESSARY, TO FIND MATRIX ENTRY.
    !
    ipl = IDLOC(Iplace,Sx,Ix)
    !
    !     FIX UP IPLACE AND IPL IF THEY POINT TO PAGING DATA.
    !     THIS IS NECESSARY BECAUSE THERE IS CONTROL INFORMATION AT THE
    !     END OF EACH PAGE.
    !
    idiff = lmx - ipl
    IF ( idiff<=1.AND.Ix(lmx-1)>0 ) THEN
      !
      !     UPDATE THE RELATIVE ADDRESS IN A NEW PAGE.
      !
      Iplace = Iplace + idiff + 1
      ipl = IDLOC(Iplace,Sx,Ix)
    ENDIF
    np = ABS(Ix(lmx-1))
  ENDIF
  100  ilast = MIN(iend,np*lpg+ll-2)
  !
  !     THE VIRTUAL END OF THE DATA FOR THIS PAGE IS ILAST.
  !
  il = IDLOC(ilast,Sx,Ix)
  il = MIN(il,lmx-2)
  !
  !     THE RELATIVE END OF DATA FOR THIS PAGE IS IL.
  !     SEARCH FOR A NONZERO VALUE WITH AN INDEX .GT. I ON THE PRESENT
  !     PAGE.
  !
  DO WHILE ( .NOT.(ipl>=il.OR.(Ix(ipl)>I.AND.Sx(ipl)/=zero)) )
    ipl = ipl + 1
  ENDDO
  !
  !     TEST IF WE HAVE FOUND THE NEXT NONZERO.
  !
  IF ( Ix(ipl)<=I.OR.Sx(ipl)==zero.OR.ipl>il ) THEN
    !
    !     UPDATE TO SCAN THE NEXT PAGE.
    ipl = ll + 1
    np = np + 1
    IF ( ilast/=iend ) GOTO 100
    !
    !     NO DATA WAS FOUND. END OF VECTOR ENCOUNTERED.
    !
    I = 0
    Xval = zero
    il = il + 1
    IF ( il==lmx-1 ) il = il + 2
    !
    !     IF A NEW ITEM WOULD BE INSERTED, IPLACE POINTS TO THE PLACE
    !     TO PUT IT.
    !
    Iplace = (np-1)*lpg + il
    RETURN
  ELSE
    I = Ix(ipl)
    Xval = Sx(ipl)
    Iplace = (np-1)*lpg + ipl
    RETURN
  ENDIF
  200 CONTINUE
  IF ( (n20046-ii)<0 ) THEN
    !
    !     ORTHOGONAL SCAN FAILED. THE VALUE J WAS NOT A SUBSCRIPT
    !     IN ANY ROW.
    !
    I = 0
    Xval = zero
    RETURN
  ELSE
    !
    !     INITIALIZE IPPLOC FOR ORTHOGONAL SCAN.
    !     LOOK FOR J AS A SUBSCRIPT IN ROWS II, II=I+1,...,L.
    !
    IF ( ii/=1 ) THEN
      ipploc = Ix(ii+3) + 1
    ELSE
      ipploc = ll + 1
    ENDIF
    iend = Ix(ii+4)
    !
    !     SCAN THROUGH SEVERAL PAGES, IF NECESSARY, TO FIND MATRIX ENTRY.
    !
    ipl = IDLOC(ipploc,Sx,Ix)
    !
    !     FIX UP IPPLOC AND IPL TO POINT TO MATRIX DATA.
    !
    idiff = lmx - ipl
    IF ( idiff<=1.AND.Ix(lmx-1)>0 ) THEN
      ipploc = ipploc + idiff + 1
      ipl = IDLOC(ipploc,Sx,Ix)
    ENDIF
    np = ABS(Ix(lmx-1))
  ENDIF
  300  ilast = MIN(iend,np*lpg+ll-2)
  il = IDLOC(ilast,Sx,Ix)
  il = MIN(il,lmx-2)
  DO WHILE ( .NOT.(ipl>=il.OR.Ix(ipl)>=j) )
    ipl = ipl + 1
  ENDDO
  !
  !     TEST IF WE HAVE FOUND THE NEXT NONZERO.
  !
  IF ( Ix(ipl)/=j.OR.Sx(ipl)==zero.OR.ipl>il ) THEN
    IF ( Ix(ipl)>=j ) ilast = iend
    ipl = ll + 1
    np = np + 1
    IF ( ilast/=iend ) GOTO 300
    ii = ii + 1
    GOTO 200
  ENDIF
  I = ii
  Xval = Sx(ipl)
  RETURN
END SUBROUTINE DPNNZR
