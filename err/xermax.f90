!*==XERMAX.f90  processed by SPAG 6.72Dc at 10:54 on  6 Feb 2019
!DECK XERMAX
SUBROUTINE XERMAX(Max)
  IMPLICIT NONE
  !*--XERMAX5
  !*** Start of declarations inserted by SPAG
  INTEGER J4SAVE , junk , Max
  !*** End of declarations inserted by SPAG
  !***BEGIN PROLOGUE  XERMAX
  !***PURPOSE  Set maximum number of times any error message is to be
  !            printed.
  !***LIBRARY   SLATEC (XERROR)
  !***CATEGORY  R3C
  !***TYPE      ALL (XERMAX-A)
  !***KEYWORDS  ERROR, XERROR
  !***AUTHOR  Jones, R. E., (SNLA)
  !***DESCRIPTION
  !
  !     Abstract
  !        XERMAX sets the maximum number of times any message
  !        is to be printed.  That is, non-fatal messages are
  !        not to be printed after they have occurred MAX times.
  !        Such non-fatal messages may be printed less than
  !        MAX times even if they occur MAX times, if error
  !        suppression mode (KONTRL=0) is ever in effect.
  !
  !     Description of Parameter
  !      --Input--
  !        MAX - the maximum number of times any one message
  !              is to be printed.
  !
  !***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
  !                 Error-handling Package, SAND82-0800, Sandia
  !                 Laboratories, 1982.
  !***ROUTINES CALLED  J4SAVE
  !***REVISION HISTORY  (YYMMDD)
  !   790801  DATE WRITTEN
  !   861211  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  !***END PROLOGUE  XERMAX
  !***FIRST EXECUTABLE STATEMENT  XERMAX
  junk = J4SAVE(4,Max,.TRUE.)
END SUBROUTINE XERMAX
