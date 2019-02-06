!*==DEVERK.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK DEVERK
      SUBROUTINE DEVERK(Lout,Kprint,Fail)
      IMPLICIT NONE
!*--DEVERK5
!***BEGIN PROLOGUE  DEVERK
!***SUBSIDIARY
!***PURPOSE  Test error returns from DPCHIP evaluators for DPCHQ1.
!***LIBRARY   SLATEC (PCHIP)
!***TYPE      DOUBLE PRECISION (EVERCK-S, DEVERK-D)
!***KEYWORDS  PCHIP EVALUATOR QUICK CHECK
!***AUTHOR  Fritsch, F. N., (LLNL)
!***DESCRIPTION
!
! --------- CODE TO TEST ERROR RETURNS FROM DPCHIP EVALUATORS. ---------
!
!
!     FORTRAN LIBRARY ROUTINES USED:  (WRITE).
!     SLATEC LIBRARY ROUTINES USED:  DCHFDV, DCHFEV, DPCHFD, DPCHFE,
!                                    XERDMP, XGETF, XSETF.
!     OTHER ROUTINES USED:  COMP.
!
!***ROUTINES CALLED  COMP, DCHFDV, DCHFEV, DPCHFD, DPCHFE, XERDMP,
!                    XGETF, XSETF
!***REVISION HISTORY  (YYMMDD)
!   820601  DATE WRITTEN
!   820715  CONVERTED TO QUICK CHECK FOR SLATEC LIBRARY.
!   890207  ADDED CALLS TO ERROR HANDLER.
!   890316  Added call to XERDMP if KPRINT.GT.2 (FNF).
!   890706  Cosmetic changes to prologue.  (WRB)
!   890911  Removed unnecessary intrinsics.  (WRB)
!   891009  Removed unreferenced statement label.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900309  Added COMP to list of routines called.  (FNF)
!   900315  Revised prologue and improved some output formats.  (FNF)
!   900316  Deleted INCFD tests because some compilers object to them,
!           and made additional minor cosmetic changes.  (FNF)
!   900322  Made miscellaneous cosmetic changes.  (FNF)
!   910708  Minor modifications in use of KPRINT.  (WRB)
!   930504  Removed parens from constants in WRITE statements.  (FNF)
!***END PROLOGUE  DEVERK
!
!  Declare arguments.
!
      INTEGER Lout , Kprint
      LOGICAL Fail
!
!  DECLARATIONS.
!
      INTEGER i , ierr , kontrl , N , nerr , next(2)
      DOUBLE PRECISION d(10) , dum , f(10) , temp , x(10)
      LOGICAL COMP , skip
!
!  INITIALIZE.
!
      PARAMETER (N=10)
!***FIRST EXECUTABLE STATEMENT  DEVERK
      nerr = 0
!
      CALL XGETF(kontrl)
      IF ( Kprint<=2 ) THEN
        CALL XSETF(0)
      ELSE
        CALL XSETF(1)
      ENDIF
!
      IF ( Kprint>=3 ) WRITE (Lout,99001)
!
!  FORMATS.
!
99001 FORMAT ('1'//10X,'TEST ERROR RETURNS')
      IF ( Kprint>=2 ) WRITE (Lout,99002)
99002 FORMAT (//10X,'DEVERK RESULTS'/10X,'--------------')
!
!  FIRST, TEST DCHFEV AND DCHFDV.
!
      IF ( Kprint>=3 ) WRITE (Lout,99005) -1
      CALL DCHFEV(0.D0,1.D0,3.D0,7.D0,3.D0,6.D0,0,dum,dum,next,ierr)
      IF ( .NOT.COMP(ierr,-1,Lout,Kprint) ) nerr = nerr + 1
!
      IF ( Kprint>=3 ) WRITE (Lout,99005) -2
      CALL DCHFEV(1.D0,1.D0,3.D0,7.D0,3.D0,6.D0,1,dum,dum,next,ierr)
      IF ( .NOT.COMP(ierr,-2,Lout,Kprint) ) nerr = nerr + 1
!
      IF ( Kprint>=3 ) WRITE (Lout,99005) -1
      CALL DCHFDV(0.D0,1.D0,3.D0,7.D0,3.D0,6.D0,0,dum,dum,dum,next,ierr)
      IF ( .NOT.COMP(ierr,-1,Lout,Kprint) ) nerr = nerr + 1
!
      IF ( Kprint>=3 ) WRITE (Lout,99005) -2
      CALL DCHFDV(1.D0,1.D0,3.D0,7.D0,3.D0,6.D0,1,dum,dum,dum,next,ierr)
      IF ( .NOT.COMP(ierr,-2,Lout,Kprint) ) nerr = nerr + 1
!
!  SET UP PCH DEFINITION.
!
      DO i = 1 , N
        x(i) = i
        f(i) = i + 2
        d(i) = 1.D0
      ENDDO
!
!  SWAP POINTS 4 AND 7, SO X-ARRAY IS OUT OF ORDER.
!
      temp = x(4)
      x(4) = x(7)
      x(7) = temp
!
!  NOW, TEST DPCHFE AND DPCHFD.
!
      IF ( Kprint>=3 ) WRITE (Lout,99005) -1
      skip = .FALSE.
      CALL DPCHFE(1,x,f,d,1,skip,0,dum,dum,ierr)
      IF ( .NOT.COMP(ierr,-1,Lout,Kprint) ) nerr = nerr + 1
!
      IF ( Kprint>=3 ) WRITE (Lout,99005) -3
      skip = .FALSE.
      CALL DPCHFE(N,x,f,d,1,skip,0,dum,dum,ierr)
      IF ( .NOT.COMP(ierr,-3,Lout,Kprint) ) nerr = nerr + 1
!
      IF ( Kprint>=3 ) WRITE (Lout,99005) -4
      skip = .TRUE.
      CALL DPCHFE(N,x,f,d,1,skip,0,dum,dum,ierr)
      IF ( .NOT.COMP(ierr,-4,Lout,Kprint) ) nerr = nerr + 1
!
      IF ( Kprint>=3 ) WRITE (Lout,99005) -1
      skip = .FALSE.
      CALL DPCHFD(1,x,f,d,1,skip,0,dum,dum,dum,ierr)
      IF ( .NOT.COMP(ierr,-1,Lout,Kprint) ) nerr = nerr + 1
!
      IF ( Kprint>=3 ) WRITE (Lout,99005) -3
      skip = .FALSE.
      CALL DPCHFD(N,x,f,d,1,skip,0,dum,dum,dum,ierr)
      IF ( .NOT.COMP(ierr,-3,Lout,Kprint) ) nerr = nerr + 1
!
      IF ( Kprint>=3 ) WRITE (Lout,99005) -4
      skip = .TRUE.
      CALL DPCHFD(N,x,f,d,1,skip,0,dum,dum,dum,ierr)
      IF ( .NOT.COMP(ierr,-4,Lout,Kprint) ) nerr = nerr + 1
!
!  SUMMARIZE RESULTS.
!
      IF ( Kprint>2 ) CALL XERDMP
      IF ( nerr==0 ) THEN
        Fail = .FALSE.
        IF ( Kprint>=2 ) WRITE (Lout,99003)
99003   FORMAT (/' ALL ERROR RETURNS OK.')
      ELSE
        Fail = .TRUE.
        IF ( Kprint>=2 ) WRITE (Lout,99004) nerr
99004   FORMAT (//' ***** TROUBLE IN DEVERK *****'//5X,I5,
     &          ' TESTS FAILED TO GIVE EXPECTED RESULTS.')
      ENDIF
!
!  TERMINATE.
!
      CALL XSETF(kontrl)
      RETURN
99005 FORMAT (/' THIS CALL SHOULD RETURN IERR =',I3)
!------------- LAST LINE OF DEVERK FOLLOWS -----------------------------
      END SUBROUTINE DEVERK
