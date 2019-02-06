!*==XERSVE.f90  processed by SPAG 6.72Dc at 11:02 on  6 Feb 2019
!DECK XERSVE
      SUBROUTINE XERSVE(Librar,Subrou,Messg,Kflag,Nerr,Level,Icount)
      IMPLICIT NONE
!*--XERSVE5
!*** Start of declarations inserted by SPAG
      INTEGER i , I1MACH , Icount , iunit , Kflag , kount , kountx , kunit , 
     &        LENTAB , Level , levtab , Nerr , nertab , nmsg , nunit
!*** End of declarations inserted by SPAG
!***BEGIN PROLOGUE  XERSVE
!***SUBSIDIARY
!***PURPOSE  Record that an error has occurred.
!***LIBRARY   SLATEC (XERROR)
!***CATEGORY  R3
!***TYPE      ALL (XERSVE-A)
!***KEYWORDS  ERROR, XERROR
!***AUTHOR  Jones, R. E., (SNLA)
!***DESCRIPTION
!
! *Usage:
!
!        INTEGER  KFLAG, NERR, LEVEL, ICOUNT
!        CHARACTER * (len) LIBRAR, SUBROU, MESSG
!
!        CALL XERSVE (LIBRAR, SUBROU, MESSG, KFLAG, NERR, LEVEL, ICOUNT)
!
! *Arguments:
!
!        LIBRAR :IN    is the library that the message is from.
!        SUBROU :IN    is the subroutine that the message is from.
!        MESSG  :IN    is the message to be saved.
!        KFLAG  :IN    indicates the action to be performed.
!                      when KFLAG > 0, the message in MESSG is saved.
!                      when KFLAG=0 the tables will be dumped and
!                      cleared.
!                      when KFLAG < 0, the tables will be dumped and
!                      not cleared.
!        NERR   :IN    is the error number.
!        LEVEL  :IN    is the error severity.
!        ICOUNT :OUT   the number of times this message has been seen,
!                      or zero if the table has overflowed and does not
!                      contain this message specifically.  When KFLAG=0,
!                      ICOUNT will not be altered.
!
! *Description:
!
!   Record that this error occurred and possibly dump and clear the
!   tables.
!
!***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
!                 Error-handling Package, SAND82-0800, Sandia
!                 Laboratories, 1982.
!***ROUTINES CALLED  I1MACH, XGETUA
!***REVISION HISTORY  (YYMMDD)
!   800319  DATE WRITTEN
!   861211  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900413  Routine modified to remove reference to KFLAG.  (WRB)
!   900510  Changed to add LIBRARY NAME and SUBROUTINE to calling
!           sequence, use IF-THEN-ELSE, make number of saved entries
!           easily changeable, changed routine name from XERSAV to
!           XERSVE.  (RWC)
!   910626  Added LIBTAB and SUBTAB to SAVE statement.  (BKS)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  XERSVE
      PARAMETER (LENTAB=10)
      INTEGER lun(5)
      CHARACTER*(*) Librar , Subrou , Messg
      CHARACTER*8 libtab(LENTAB) , subtab(LENTAB) , lib , sub
      CHARACTER*20 mestab(LENTAB) , mes
      DIMENSION nertab(LENTAB) , levtab(LENTAB) , kount(LENTAB)
      SAVE libtab , subtab , mestab , nertab , levtab , kount , kountx , nmsg
      DATA kountx/0/ , nmsg/0/
!***FIRST EXECUTABLE STATEMENT  XERSVE
!
      IF ( Kflag<=0 ) THEN
!
!        Dump the table.
!
        IF ( nmsg==0 ) RETURN
!
!        Print to each unit.
!
        CALL XGETUA(lun,nunit)
        DO kunit = 1 , nunit
          iunit = lun(kunit)
          IF ( iunit==0 ) iunit = I1MACH(4)
!
!           Print the table header.
!
          WRITE (iunit,99001)
!
!     Formats.
!
99001     FORMAT ('0          ERROR MESSAGE SUMMARY'/
     &            ' LIBRARY    SUBROUTINE MESSAGE START             NERR',
     &            '     LEVEL     COUNT')
!
!           Print body of table.
!
          DO i = 1 , nmsg
            WRITE (iunit,99002) libtab(i) , subtab(i) , mestab(i) , nertab(i) , 
     &                          levtab(i) , kount(i)
99002       FORMAT (1X,A,3X,A,3X,A,3I10)
          ENDDO
!
!           Print number of other errors.
!
          IF ( kountx/=0 ) WRITE (iunit,99003) kountx
99003     FORMAT ('0OTHER ERRORS NOT INDIVIDUALLY TABULATED = ',I10)
          WRITE (iunit,99004)
99004     FORMAT (1X)
        ENDDO
!
!        Clear the error tables.
!
        IF ( Kflag==0 ) THEN
          nmsg = 0
          kountx = 0
        ENDIF
      ELSE
!
!        PROCESS A MESSAGE...
!        SEARCH FOR THIS MESSG, OR ELSE AN EMPTY SLOT FOR THIS MESSG,
!        OR ELSE DETERMINE THAT THE ERROR TABLE IS FULL.
!
        lib = Librar
        sub = Subrou
        mes = Messg
        DO i = 1 , nmsg
          IF ( lib==libtab(i).AND.sub==subtab(i).AND.mes==mestab(i).AND.
     &         Nerr==nertab(i).AND.Level==levtab(i) ) THEN
            kount(i) = kount(i) + 1
            Icount = kount(i)
            RETURN
          ENDIF
        ENDDO
!
        IF ( nmsg<LENTAB ) THEN
!
!           Empty slot found for new message.
!
          nmsg = nmsg + 1
          libtab(i) = lib
          subtab(i) = sub
          mestab(i) = mes
          nertab(i) = Nerr
          levtab(i) = Level
          kount(i) = 1
          Icount = 1
        ELSE
!
!           Table is full.
!
          kountx = kountx + 1
          Icount = 0
        ENDIF
      ENDIF
      RETURN
      END SUBROUTINE XERSVE
