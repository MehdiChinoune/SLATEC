!** PCHNGS
SUBROUTINE PCHNGS(Ii,Xval,Iplace,Sx,Ix,Ircx)
  !>
  !  Subsidiary to SPLP
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      SINGLE PRECISION (PCHNGS-S, DPCHNG-D)
  !***
  ! **Author:**  Hanson, R. J., (SNLA)
  !           Wisniewski, J. A., (SNLA)
  !***
  ! **Description:**
  !
  !     PCHNGS LIMITS THE TYPE OF STORAGE TO A SEQUENTIAL SCHEME.
  !     SPARSE MATRIX ELEMENT ALTERATION SUBROUTINE.
  !
  !     SUBROUTINE PCHNGS() CHANGES ELEMENT II IN VECTOR +/- IRCX TO THE
  !     VALUE XVAL.
  !
  !            II THE ABSOLUTE VALUE OF THIS INTEGER IS THE SUBSCRIPT FOR
  !               THE ELEMENT TO BE CHANGED.
  !          XVAL NEW VALUE OF THE MATRIX ELEMENT BEING CHANGED.
  !     IPLACE POINTER INFORMATION WHICH IS MAINTAINED BY THE PACKAGE.
  !   SX(*),IX(*) THE WORK ARRAYS WHICH ARE USED TO STORE THE SPARSE
  !               MATRIX. THESE ARRAYS ARE AUTOMATICALLY MAINTAINED BY THE
  !               PACKAGE FOR THE USER.
  !          IRCX POINTS TO THE VECTOR OF THE MATRIX BEING UPDATED.
  !               A NEGATIVE VALUE OF IRCX INDICATES THAT ROW -IRCX IS
  !               BEING UPDATED.  A POSITIVE VALUE OF IRCX INDICATES THAT
  !               COLUMN IRCX IS BEING UPDATED.  A ZERO VALUE OF IRCX IS
  !               AN ERROR.
  !
  !     SINCE DATA ITEMS ARE KEPT SORTED IN THE SEQUENTIAL DATA STRUCTURE,
  !     CHANGING A MATRIX ELEMENT CAN REQUIRE THE MOVEMENT OF ALL THE DATA
  !     ITEMS IN THE MATRIX. FOR THIS REASON, IT IS SUGGESTED THAT DATA
  !     ITEMS BE ADDED A COL. AT A TIME, IN ASCENDING COL. SEQUENCE.
  !     FURTHERMORE, SINCE DELETING ITEMS FROM THE DATA STRUCTURE MAY ALSO
  !     REQUIRE MOVING LARGE AMOUNTS OF DATA, ZERO ELEMENTS ARE EXPLICITLY
  !     STORED IN THE MATRIX.
  !
  !     THIS SUBROUTINE IS A MODIFICATION OF THE SUBROUTINE LCHNGS,
  !     SANDIA LABS. REPT. SAND78-0785.
  !     MODIFICATIONS BY K.L. HIEBERT AND R.J. HANSON
  !     REVISED 811130-1000
  !     REVISED YYMMDD-HHMM
  !
  !***
  ! **See also:**  SPLP
  !***
  ! **Routines called:**  IPLOC, PRWPGE, XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   811215  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900328  Added TYPE section.  (WRB)
  !   910403  Updated AUTHOR and DESCRIPTION sections.  (WRB)
  USE service, ONLY : XERMSG
  INTEGER i, iend, Ii, il, ilast, iopt, ipl, Iplace, Ircx, istart, &
    Ix(*), ixlast, j, jj, jstart, k, key, ll, lmx, lpg, n20055, nerr, np
  REAL Sx(*), Xval, sxlast, sxval
  REAL, PARAMETER :: zero = 0.E0, one = 1.E0
  !* FIRST EXECUTABLE STATEMENT  PCHNGS
  iopt = 1
  !
  !     DETERMINE NULL-CASES..
  IF ( Ii==0 ) RETURN
  !
  !     CHECK VALIDITY OF ROW/COL. INDEX.
  !
  IF ( Ircx==0 ) THEN
    nerr = 55
    CALL XERMSG('SLATEC','PCHNGS','IRCX=0.',nerr,iopt)
  END IF
  lmx = Ix(1)
  !
  !     LMX IS THE LENGTH OF THE IN-MEMORY STORAGE AREA.
  !
  IF ( Ircx>=0 ) THEN
    !
    !     CHECK SUBSCRIPTS OF THE COLUMN. THE COL. NUMBER MUST BE .LE. N AND
    !     THE INDEX MUST BE .LE. M.
    !
    IF ( Ix(3)<Ircx.OR.Ix(2)<ABS(Ii) ) THEN
      nerr = 55
      CALL XERMSG('SLATEC','PCHNGS',&
        'SUBSCRIPTS FOR ARRAY ELEMENT TO BE ACCESSED WERE OUT OF BOUNDS.',nerr,iopt)
    END IF
    !
    !     CHECK SUBSCRIPTS OF THE ROW. THE ROW NUMBER MUST BE .LE. M AND
    !     THE INDEX MUST BE .LE. N.
    !
  ELSEIF ( Ix(2)<-Ircx.OR.Ix(3)<ABS(Ii) ) THEN
    nerr = 55
    CALL XERMSG('SLATEC','PCHNGS',&
      'SUBSCRIPTS FOR ARRAY ELEMENT TO BE ACCESSED WERE OUT OF BOUNDS.',nerr,iopt)
  END IF
  !
  !     SET I TO BE THE ELEMENT OF ROW/COLUMN J TO BE CHANGED.
  !
  IF ( Ircx<=0 ) THEN
    i = ABS(Ircx)
    j = ABS(Ii)
  ELSE
    i = ABS(Ii)
    j = ABS(Ircx)
  END IF
  !
  !     THE INTEGER LL POINTS TO THE START OF THE MATRIX ELEMENT DATA.
  !
  ll = Ix(3) + 4
  Ii = ABS(Ii)
  lpg = lmx - ll
  !
  !     SET IPLACE TO START OUR SCAN FOR THE ELEMENT AT THE BEGINNING
  !     OF THE VECTOR.
  !
  IF ( j/=1 ) THEN
    Iplace = Ix(j+3) + 1
  ELSE
    Iplace = ll + 1
  END IF
  !
  !     IEND POINTS TO THE LAST ELEMENT OF THE VECTOR TO BE SCANNED.
  !
  iend = Ix(j+4)
  !
  !     SCAN THROUGH SEVERAL PAGES, IF NECESSARY, TO FIND MATRIX ELEMENT.
  !
  ipl = IPLOC(Iplace,Sx,Ix)
  np = ABS(Ix(lmx-1))
  !
  !     THE VIRTUAL END OF DATA FOR THIS PAGE IS ILAST.
  !
  100  ilast = MIN(iend,np*lpg+ll-2)
  !
  !     THE RELATIVE END OF DATA FOR THIS PAGE IS IL.
  !     SEARCH FOR A MATRIX VALUE WITH AN INDEX .GE. I ON THE PRESENT
  !     PAGE.
  !
  il = IPLOC(ilast,Sx,Ix)
  il = MIN(il,lmx-2)
  DO WHILE ( .NOT.(ipl>=il.OR.Ix(ipl)>=i) )
    ipl = ipl + 1
  END DO
  !
  !     SET IPLACE AND STORE DATA ITEM IF FOUND.
  !
  IF ( Ix(ipl)/=i.OR.ipl>il ) THEN
    !
    !     EXIT FROM LOOP IF ITEM WAS FOUND.
    !
    IF ( Ix(ipl)>i.AND.ipl<=il ) ilast = iend
    IF ( ilast/=iend ) THEN
      ipl = ll + 1
      np = np + 1
    END IF
    IF ( ilast/=iend ) GOTO 100
    !
    !     INSERT NEW DATA ITEM INTO LOCATION AT IPLACE(IPL).
    !
    IF ( ipl>il.OR.(ipl==il.AND.i>Ix(ipl)) ) THEN
      ipl = il + 1
      IF ( ipl==lmx-1 ) ipl = ipl + 2
    END IF
    Iplace = (np-1)*lpg + ipl
    !
    !     GO TO A NEW PAGE, IF NECESSARY, TO INSERT THE ITEM.
    !
    IF ( ipl<=lmx.OR.Ix(lmx-1)>=0 ) ipl = IPLOC(Iplace,Sx,Ix)
    iend = Ix(ll)
    np = ABS(Ix(lmx-1))
    !
    !     LOOP THROUGH ALL SUBSEQUENT PAGES OF THE MATRIX MOVING DATA DOWN.
    !     THIS IS NECESSARY TO MAKE ROOM FOR THE NEW MATRIX ELEMENT AND
    !     KEEP THE ENTRIES SORTED.
    !
    sxval = Xval
  ELSE
    Sx(ipl) = Xval
    Sx(lmx) = one
    RETURN
  END IF
  200  ilast = MIN(iend,np*lpg+ll-2)
  il = IPLOC(ilast,Sx,Ix)
  il = MIN(il,lmx-2)
  sxlast = Sx(il)
  ixlast = Ix(il)
  istart = ipl + 1
  IF ( istart<=il ) THEN
    k = istart + il
    DO jj = istart, il
      Sx(k-jj) = Sx(k-jj-1)
      Ix(k-jj) = Ix(k-jj-1)
    END DO
    Sx(lmx) = one
  END IF
  IF ( ipl<=lmx ) THEN
    Sx(ipl) = sxval
    Ix(ipl) = i
    sxval = sxlast
    i = ixlast
    Sx(lmx) = one
    IF ( Ix(lmx-1)>0 ) THEN
      ipl = ll + 1
      np = np + 1
    END IF
  END IF
  IF ( Ix(lmx-1)>0 ) GOTO 200
  np = ABS(Ix(lmx-1))
  !
  !     DETERMINE IF A NEW PAGE IS TO BE CREATED FOR THE LAST ELEMENT
  !     MOVED DOWN.
  !
  il = il + 1
  IF ( il==lmx-1 ) THEN
    !
    !     CREATE A NEW PAGE.
    !
    Ix(lmx-1) = np
    !
    !     WRITE THE OLD PAGE.
    !
    Sx(lmx) = zero
    key = 2
    CALL PRWPGE(key,np,lpg,Sx,Ix)
    Sx(lmx) = one
    !
    !     STORE LAST ELEMENT MOVED DOWN IN A NEW PAGE.
    !
    ipl = ll + 1
    np = np + 1
    Ix(lmx-1) = -np
    Sx(ipl) = sxval
    Ix(ipl) = i
    !
    !     LAST ELEMENT MOVED REMAINED ON THE OLD PAGE.
    !
  ELSEIF ( ipl/=il ) THEN
    Sx(il) = sxval
    Ix(il) = i
    Sx(lmx) = one
  END IF
  !
  !     INCREMENT POINTERS TO LAST ELEMENT IN VECTORS J,J+1,... .
  !
  jstart = j + 4
  jj = jstart
  n20055 = ll
  DO WHILE ( (n20055-jj)>=0 )
    Ix(jj) = Ix(jj) + 1
    IF ( MOD(Ix(jj)-ll,lpg)==lpg-1 ) Ix(jj) = Ix(jj) + 2
    jj = jj + 1
  END DO
  !
  !     IPLACE POINTS TO THE INSERTED DATA ITEM.
  !
  ipl = IPLOC(Iplace,Sx,Ix)
END SUBROUTINE PCHNGS
