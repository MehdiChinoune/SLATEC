!** J4SAVE
INTEGER FUNCTION J4SAVE(Iwhich,Ivalue,Iset)
  !>
  !***
  !  Save or recall global variables needed by error
  !            handling routines.
  !***
  ! **Library:**   SLATEC (XERROR)
  !***
  ! **Type:**      INTEGER (J4SAVE-I)
  !***
  ! **Keywords:**  ERROR MESSAGES, ERROR NUMBER, RECALL, SAVE, XERROR
  !***
  ! **Author:**  Jones, R. E., (SNLA)
  !***
  ! **Description:**
  !
  !     Abstract
  !        J4SAVE saves and recalls several global variables needed
  !        by the library error handling routines.
  !
  !     Description of Parameters
  !      --Input--
  !        IWHICH - Index of item desired.
  !                = 1 Refers to current error number.
  !                = 2 Refers to current error control flag.
  !                = 3 Refers to current unit number to which error
  !                    messages are to be sent.  (0 means use standard.)
  !                = 4 Refers to the maximum number of times any
  !                     message is to be printed (as set by XERMAX).
  !                = 5 Refers to the total number of units to which
  !                     each error message is to be written.
  !                = 6 Refers to the 2nd unit for error messages
  !                = 7 Refers to the 3rd unit for error messages
  !                = 8 Refers to the 4th unit for error messages
  !                = 9 Refers to the 5th unit for error messages
  !        IVALUE - The value to be set for the IWHICH-th parameter,
  !                 if ISET is .TRUE. .
  !        ISET   - If ISET=.TRUE., the IWHICH-th parameter will BE
  !                 given the value, IVALUE.  If ISET=.FALSE., the
  !                 IWHICH-th parameter will be unchanged, and IVALUE
  !                 is a dummy parameter.
  !      --Output--
  !        The (old) value of the IWHICH-th parameter will be returned
  !        in the function value, J4SAVE.
  !
  !***
  ! **See also:**  XERMSG
  !***
  ! **References:**  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
  !                 Error-handling Package, SAND82-0800, Sandia
  !                 Laboratories, 1982.
  !***
  ! **Routines called:**  (NONE)

  !* REVISION HISTORY  (YYMMDD)
  !   790801  DATE WRITTEN
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900205  Minor modifications to prologue.  (WRB)
  !   900402  Added TYPE section.  (WRB)
  !   910411  Added KEYWORDS section.  (WRB)
  !   920501  Reformatted the REFERENCES section.  (WRB)

  INTEGER Ivalue, Iwhich
  LOGICAL Iset
  INTEGER, SAVE :: iparam(9) = [ 0, 2, 0, 10, 1, 0, 0, 0, 0 ]
  !* FIRST EXECUTABLE STATEMENT  J4SAVE
  J4SAVE = iparam(Iwhich)
  IF ( Iset ) iparam(Iwhich) = Ivalue
END FUNCTION J4SAVE
