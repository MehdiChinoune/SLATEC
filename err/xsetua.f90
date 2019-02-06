!*==XSETUA.f90  processed by SPAG 6.72Dc at 10:54 on  6 Feb 2019
!DECK XSETUA
      SUBROUTINE XSETUA(Iunita,N)
      IMPLICIT NONE
!*--XSETUA5
!*** Start of declarations inserted by SPAG
      INTEGER i , index , Iunita , J4SAVE , junk , N
!*** End of declarations inserted by SPAG
!***BEGIN PROLOGUE  XSETUA
!***PURPOSE  Set logical unit numbers (up to 5) to which error
!            messages are to be sent.
!***LIBRARY   SLATEC (XERROR)
!***CATEGORY  R3B
!***TYPE      ALL (XSETUA-A)
!***KEYWORDS  ERROR, XERROR
!***AUTHOR  Jones, R. E., (SNLA)
!***DESCRIPTION
!
!     Abstract
!        XSETUA may be called to declare a list of up to five
!        logical units, each of which is to receive a copy of
!        each error message processed by this package.
!        The purpose of XSETUA is to allow simultaneous printing
!        of each error message on, say, a main output file,
!        an interactive terminal, and other files such as graphics
!        communication files.
!
!     Description of Parameters
!      --Input--
!        IUNIT - an array of up to five unit numbers.
!                Normally these numbers should all be different
!                (but duplicates are not prohibited.)
!        N     - the number of unit numbers provided in IUNIT
!                must have 1 .LE. N .LE. 5.
!
!***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
!                 Error-handling Package, SAND82-0800, Sandia
!                 Laboratories, 1982.
!***ROUTINES CALLED  J4SAVE, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   790801  DATE WRITTEN
!   861211  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900510  Change call to XERRWV to XERMSG.  (RWC)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  XSETUA
      DIMENSION Iunita(5)
      CHARACTER*8 xern1
!***FIRST EXECUTABLE STATEMENT  XSETUA
!
      IF ( N<1.OR.N>5 ) THEN
        WRITE (xern1,'(I8)') N
        CALL XERMSG('SLATEC','XSETUA','INVALID NUMBER OF UNITS, N = '//xern1,1,
     &              2)
        RETURN
      ENDIF
!
      DO i = 1 , N
        index = i + 4
        IF ( i==1 ) index = 3
        junk = J4SAVE(index,Iunita(i),.TRUE.)
      ENDDO
      junk = J4SAVE(5,N,.TRUE.)
      END SUBROUTINE XSETUA
