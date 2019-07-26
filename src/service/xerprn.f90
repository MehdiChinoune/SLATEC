!** XERPRN
SUBROUTINE XERPRN(Prefix,Npref,Messg,Nwrap)
  !> Print error messages processed by XERMSG.
  !***
  ! **Library:**   SLATEC (XERROR)
  !***
  ! **Category:**  R3C
  !***
  ! **Type:**      ALL (XERPRN-A)
  !***
  ! **Keywords:**  ERROR MESSAGES, PRINTING, XERROR
  !***
  ! **Author:**  Fong, Kirby, (NMFECC at LLNL)
  !***
  ! **Description:**
  !
  ! This routine sends one or more lines to each of the (up to five)
  ! logical units to which error messages are to be sent.  This routine
  ! is called several times by XERMSG, sometimes with a single line to
  ! print and sometimes with a (potentially very long) message that may
  ! wrap around into multiple lines.
  !
  ! PREFIX  Input argument of type CHARACTER.  This argument contains
  !         characters to be put at the beginning of each line before
  !         the body of the message.  No more than 16 characters of
  !         PREFIX will be used.
  !
  ! NPREF   Input argument of type INTEGER.  This argument is the number
  !         of characters to use from PREFIX.  If it is negative, the
  !         intrinsic function LEN is used to determine its length.  If
  !         it is zero, PREFIX is not used.  If it exceeds 16 or if
  !         LEN(PREFIX) exceeds 16, only the first 16 characters will be
  !         used.  If NPREF is positive and the length of PREFIX is less
  !         than NPREF, a copy of PREFIX extended with blanks to length
  !         NPREF will be used.
  !
  ! MESSG   Input argument of type CHARACTER.  This is the text of a
  !         message to be printed.  If it is a long message, it will be
  !         broken into pieces for printing on multiple lines.  Each line
  !         will start with the appropriate prefix and be followed by a
  !         piece of the message.  NWRAP is the number of characters per
  !         piece; that is, after each NWRAP characters, we break and
  !         start a new line.  In addition the characters '$$' embedded
  !         in MESSG are a sentinel for a new line.  The counting of
  !         characters up to NWRAP starts over for each new line.  The
  !         value of NWRAP typically used by XERMSG is 72 since many
  !         older error messages in the SLATEC Library are laid out to
  !         rely on wrap-around every 72 characters.
  !
  ! NWRAP   Input argument of type INTEGER.  This gives the maximum size
  !         piece into which to break MESSG for printing on multiple
  !         lines.  An embedded '$$' ends a line, and the count restarts
  !         at the following character.  If a line break does not occur
  !         on a blank (it would split a word) that word is moved to the
  !         next line.  Values of NWRAP less than 16 will be treated as
  !         16.  Values of NWRAP greater than 132 will be treated as 132.
  !         The actual line length will be NPREF + NWRAP after NPREF has
  !         been adjusted to fall between 0 and 16 and NWRAP has been
  !         adjusted to fall between 16 and 132.
  !
  !***
  ! **References:**  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
  !                 Error-handling Package, SAND82-0800, Sandia
  !                 Laboratories, 1982.
  !***
  ! **Routines called:**  I1MACH, XGETUA

  !* REVISION HISTORY  (YYMMDD)
  !   880621  DATE WRITTEN
  !   880708  REVISED AFTER THE SLATEC CML SUBCOMMITTEE MEETING OF
  !           JUNE 29 AND 30 TO CHANGE THE NAME TO XERPRN AND TO REWORK
  !           THE HANDLING OF THE NEW LINE SENTINEL TO BEHAVE LIKE THE
  !           SLASH CHARACTER IN FORMAT STATEMENTS.
  !   890706  REVISED WITH THE HELP OF FRED FRITSCH AND REG CLEMENS TO
  !           STREAMLINE THE CODING AND FIX A BUG THAT CAUSED EXTRA BLANK
  !           LINES TO BE PRINTED.
  !   890721  REVISED TO ADD A NEW FEATURE.  A NEGATIVE VALUE OF NPREF
  !           CAUSES LEN(PREFIX) TO BE USED AS THE LENGTH.
  !   891013  REVISED TO CORRECT ERROR IN CALCULATING PREFIX LENGTH.
  !   891214  Prologue converted to Version 4.0 format.  (WRB)
  !   900510  Added code to break messages between words.  (RWC)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  USE ISO_FORTRAN_ENV, ONLY : ERROR_UNIT
  !
  INTEGER :: i, idelta, lenmsg, lpiece, lpref, lwrap, n, nextc
  CHARACTER(*) Prefix, Messg
  INTEGER :: Npref, Nwrap
  CHARACTER(148) :: cbuff
  CHARACTER(2), PARAMETER :: NEWLIN = '$$'
  !* FIRST EXECUTABLE STATEMENT  XERPRN

  !
  !       LPREF IS THE LENGTH OF THE PREFIX.  THE PREFIX IS PLACED AT THE
  !       BEGINNING OF CBUFF, THE CHARACTER BUFFER, AND KEPT THERE DURING
  !       THE REST OF THIS ROUTINE.
  !
  IF( Npref<0 ) THEN
    lpref = LEN(Prefix)
  ELSE
    lpref = Npref
  END IF
  lpref = MIN(16,lpref)
  IF( lpref/=0 ) cbuff(1:lpref) = Prefix
  !
  !       LWRAP IS THE MAXIMUM NUMBER OF CHARACTERS WE WANT TO TAKE AT ONE
  !       TIME FROM MESSG TO PRINT ON ONE LINE.
  !
  lwrap = MAX(16,MIN(132,Nwrap))
  !
  !       SET LENMSG TO THE LENGTH OF MESSG, IGNORE ANY TRAILING BLANKS.
  !
  lenmsg = LEN(Messg)
  n = lenmsg
  DO i = 1, n
    IF( Messg(lenmsg:lenmsg)/=' ' ) EXIT
    lenmsg = lenmsg - 1
  END DO
  !
  !       IF THE MESSAGE IS ALL BLANKS, THEN PRINT ONE BLANK LINE.
  !
  IF( lenmsg==0 ) THEN
    cbuff(lpref+1:lpref+1) = ' '
    WRITE (ERROR_UNIT,'(A)') cbuff(1:lpref+1)
    RETURN
  END IF
  !
  !       SET NEXTC TO THE POSITION IN MESSG WHERE THE NEXT SUBSTRING
  !       STARTS.  FROM THIS POSITION WE SCAN FOR THE NEW LINE SENTINEL.
  !       WHEN NEXTC EXCEEDS LENMSG, THERE IS NO MORE TO PRINT.
  !       WE LOOP BACK TO LABEL 50 UNTIL ALL PIECES HAVE BEEN PRINTED.
  !
  !       WE LOOK FOR THE NEXT OCCURRENCE OF THE NEW LINE SENTINEL.  THE
  !       INDEX INTRINSIC FUNCTION RETURNS ZERO IF THERE IS NO OCCURRENCE
  !       OR IF THE LENGTH OF THE FIRST ARGUMENT IS LESS THAN THE LENGTH
  !       OF THE SECOND ARGUMENT.
  !
  !       THERE ARE SEVERAL CASES WHICH SHOULD BE CHECKED FOR IN THE
  !       FOLLOWING ORDER.  WE ARE ATTEMPTING TO SET LPIECE TO THE NUMBER
  !       OF CHARACTERS THAT SHOULD BE TAKEN FROM MESSG STARTING AT
  !       POSITION NEXTC.
  !
  !       LPIECE = 0   THE NEW LINE SENTINEL DOES NOT OCCUR IN THE
  !                       REMAINDER OF THE CHARACTER STRING.  LPIECE
  !                       SHOULD BE SET TO LWRAP OR LENMSG+1-NEXTC,
  !                       WHICHEVER IS LESS.
  !
  !       LPIECE = 1   THE NEW LINE SENTINEL STARTS AT MESSG(NEXTC:
  !                       NEXTC).  LPIECE IS EFFECTIVELY ZERO, AND WE
  !                       PRINT NOTHING TO AVOID PRODUCING UNNECESSARY
  !                       BLANK LINES.  THIS TAKES CARE OF THE SITUATION
  !                       WHERE THE LIBRARY ROUTINE HAS A MESSAGE OF
  !                       EXACTLY 72 CHARACTERS FOLLOWED BY A NEW LINE
  !                       SENTINEL FOLLOWED BY MORE CHARACTERS.  NEXTC
  !                       SHOULD BE INCREMENTED BY 2.
  !
  !       LPIECE > LWRAP+1  REDUCE LPIECE TO LWRAP.
  !
  !       ELSE            THIS LAST CASE MEANS 2 <= LPIECE <= LWRAP+1
  !                       RESET LPIECE = LPIECE-1.  NOTE THAT THIS
  !                       PROPERLY HANDLES THE END CASE WHERE LPIECE =
  !                       LWRAP+1.  THAT IS, THE SENTINEL FALLS EXACTLY
  !                       AT THE END OF A LINE.
  !
  nextc = 1
  DO
    lpiece = INDEX(Messg(nextc:lenmsg),NEWLIN)
    IF( lpiece==0 ) THEN
      !
      !       THERE WAS NO NEW LINE SENTINEL FOUND.
      !
      idelta = 0
      lpiece = MIN(lwrap,lenmsg+1-nextc)
      IF( lpiece<lenmsg+1-nextc ) THEN
        DO i = lpiece + 1, 2, -1
          IF( Messg(nextc+i-1:nextc+i-1)==' ' ) THEN
            lpiece = i - 1
            idelta = 1
            EXIT
          END IF
        END DO
      END IF
      cbuff(lpref+1:lpref+lpiece) = Messg(nextc:nextc+lpiece-1)
      nextc = nextc + lpiece + idelta
    ELSEIF( lpiece==1 ) THEN
      !
      !       WE HAVE A NEW LINE SENTINEL AT MESSG(NEXTC:NEXTC+1).
      !       DON'T PRINT A BLANK LINE.
      !
      nextc = nextc + 2
      CYCLE
    ELSEIF( lpiece>lwrap+1 ) THEN
      !
      !       LPIECE SHOULD BE SET DOWN TO LWRAP.
      !
      idelta = 0
      lpiece = lwrap
      DO i = lpiece + 1, 2, -1
        IF( Messg(nextc+i-1:nextc+i-1)==' ' ) THEN
          lpiece = i - 1
          idelta = 1
          EXIT
        END IF
      END DO
      cbuff(lpref+1:lpref+lpiece) = Messg(nextc:nextc+lpiece-1)
      nextc = nextc + lpiece + idelta
    ELSE
      !
      !       IF WE ARRIVE HERE, IT MEANS 2 <= LPIECE <= LWRAP+1.
      !       WE SHOULD DECREMENT LPIECE BY ONE.
      !
      lpiece = lpiece - 1
      cbuff(lpref+1:lpref+lpiece) = Messg(nextc:nextc+lpiece-1)
      nextc = nextc + lpiece + 2
    END IF
    !
    !       PRINT
    !
    WRITE (ERROR_UNIT,'(A)') cbuff(1:lpref+lpiece)
    !
    IF( nextc>lenmsg ) EXIT
  END DO
END SUBROUTINE XERPRN
