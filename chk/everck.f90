!DECK EVERCK
SUBROUTINE EVERCK(Lout,Kprint,Fail)
  IMPLICIT NONE
  !***BEGIN PROLOGUE  EVERCK
  !***SUBSIDIARY
  !***PURPOSE  Test error returns from PCHIP evaluators for PCHQK1.
  !***LIBRARY   SLATEC (PCHIP)
  !***TYPE      SINGLE PRECISION (EVERCK-S, DEVERK-D)
  !***KEYWORDS  PCHIP EVALUATOR QUICK CHECK
  !***AUTHOR  Fritsch, F. N., (LLNL)
  !***DESCRIPTION
  !
  ! --------- CODE TO TEST ERROR RETURNS FROM PCHIP EVALUATORS. ---------
  !
  !
  !     FORTRAN LIBRARY ROUTINES USED:  (WRITE).
  !     SLATEC LIBRARY ROUTINES USED:  CHFDV, CHFEV, PCHFD, PCHFE,
  !                                    XERDMP, XGETF, XSETF.
  !     OTHER ROUTINES USED:  COMP.
  !
  !***ROUTINES CALLED  CHFDV, CHFEV, COMP, PCHFD, PCHFE, XERDMP, XGETF,
  !                    XSETF
  !***REVISION HISTORY  (YYMMDD)
  !   820601  DATE WRITTEN
  !   820715  CONVERTED TO QUICK CHECK FOR SLATEC LIBRARY.
  !   890207  ADDED CALLS TO ERROR HANDLER.
  !   890316  Added call to XERDMP if KPRINT.GT.2 (FNF).
  !   890629  Appended E0 to real constants to reduce S.P./D.P.
  !           differences.
  !   890706  Cosmetic changes to prologue.  (WRB)
  !   891009  Removed unreferenced statement label.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900309  Added COMP to list of routines called.  (FNF)
  !   900315  Revised prologue and improved some output formats.  (FNF)
  !   900316  Deleted INCFD tests because some compilers object to them,
  !           and made additional minor cosmetic changes.  (FNF)
  !   900322  Made miscellaneous cosmetic changes.  (FNF)
  !   910708  Minor modifications in use of KPRINT.  (WRB)
  !   930504  Removed parens from constants in WRITE statements.  (FNF)
  !***END PROLOGUE  EVERCK
  !
  !  Declare arguments.
  !
  INTEGER Lout, Kprint
  LOGICAL Fail
  !
  !  DECLARATIONS.
  !
  INTEGER i, ierr, kontrl, N, nerr, next(2)
  REAL d(10), dum, f(10), temp, x(10)
  LOGICAL COMP, skip
  !
  !  INITIALIZE.
  !
  PARAMETER (N=10)
  !***FIRST EXECUTABLE STATEMENT  EVERCK
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
  99002 FORMAT (//10X,'EVERCK RESULTS'/10X,'--------------')
  !
  !  FIRST, TEST CHFEV AND CHFDV.
  !
  IF ( Kprint>=3 ) WRITE (Lout,99005) -1
  CALL CHFEV(0.E0,1.E0,3.E0,7.E0,3.E0,6.E0,0,dum,dum,next,ierr)
  IF ( .NOT.COMP(ierr,-1,Lout,Kprint) ) nerr = nerr + 1
  !
  IF ( Kprint>=3 ) WRITE (Lout,99005) -2
  CALL CHFEV(1.E0,1.E0,3.E0,7.E0,3.E0,6.E0,1,dum,dum,next,ierr)
  IF ( .NOT.COMP(ierr,-2,Lout,Kprint) ) nerr = nerr + 1
  !
  IF ( Kprint>=3 ) WRITE (Lout,99005) -1
  CALL CHFDV(0.E0,1.E0,3.E0,7.E0,3.E0,6.E0,0,dum,dum,dum,next,ierr)
  IF ( .NOT.COMP(ierr,-1,Lout,Kprint) ) nerr = nerr + 1
  !
  IF ( Kprint>=3 ) WRITE (Lout,99005) -2
  CALL CHFDV(1.E0,1.E0,3.E0,7.E0,3.E0,6.E0,1,dum,dum,dum,next,ierr)
  IF ( .NOT.COMP(ierr,-2,Lout,Kprint) ) nerr = nerr + 1
  !
  !  SET UP PCH DEFINITION.
  !
  DO i = 1, N
    x(i) = i
    f(i) = i + 2
    d(i) = 1.E0
  ENDDO
  !
  !  SWAP POINTS 4 AND 7, SO X-ARRAY IS OUT OF ORDER.
  !
  temp = x(4)
  x(4) = x(7)
  x(7) = temp
  !
  !  NOW, TEST PCHFE AND PCHFD.
  !
  IF ( Kprint>=3 ) WRITE (Lout,99005) -1
  skip = .FALSE.
  CALL PCHFE(1,x,f,d,1,skip,0,dum,dum,ierr)
  IF ( .NOT.COMP(ierr,-1,Lout,Kprint) ) nerr = nerr + 1
  !
  IF ( Kprint>=3 ) WRITE (Lout,99005) -3
  skip = .FALSE.
  CALL PCHFE(N,x,f,d,1,skip,0,dum,dum,ierr)
  IF ( .NOT.COMP(ierr,-3,Lout,Kprint) ) nerr = nerr + 1
  !
  IF ( Kprint>=3 ) WRITE (Lout,99005) -4
  skip = .TRUE.
  CALL PCHFE(N,x,f,d,1,skip,0,dum,dum,ierr)
  IF ( .NOT.COMP(ierr,-4,Lout,Kprint) ) nerr = nerr + 1
  !
  IF ( Kprint>=3 ) WRITE (Lout,99005) -1
  skip = .FALSE.
  CALL PCHFD(1,x,f,d,1,skip,0,dum,dum,dum,ierr)
  IF ( .NOT.COMP(ierr,-1,Lout,Kprint) ) nerr = nerr + 1
  !
  IF ( Kprint>=3 ) WRITE (Lout,99005) -3
  skip = .FALSE.
  CALL PCHFD(N,x,f,d,1,skip,0,dum,dum,dum,ierr)
  IF ( .NOT.COMP(ierr,-3,Lout,Kprint) ) nerr = nerr + 1
  !
  IF ( Kprint>=3 ) WRITE (Lout,99005) -4
  skip = .TRUE.
  CALL PCHFD(N,x,f,d,1,skip,0,dum,dum,dum,ierr)
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
    99004   FORMAT (//' ***** TROUBLE IN EVERCK *****'//5X,I5,&
      ' TESTS FAILED TO GIVE EXPECTED RESULTS.')
  ENDIF
  !
  !  TERMINATE.
  !
  CALL XSETF(kontrl)
  RETURN
  99005 FORMAT (/' THIS CALL SHOULD RETURN IERR =',I3)
  !------------- LAST LINE OF EVERCK FOLLOWS -----------------------------
END SUBROUTINE EVERCK
