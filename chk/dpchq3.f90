!*==DPCHQ3.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK DPCHQ3
SUBROUTINE DPCHQ3(Lun,Kprint,Ipass)
  IMPLICIT NONE
  !*--DPCHQ35
  !***BEGIN PROLOGUE  DPCHQ3
  !***PURPOSE  Test the PCHIP interpolators DPCHIC, DPCHIM, DPCHSP.
  !***LIBRARY   SLATEC (PCHIP)
  !***TYPE      DOUBLE PRECISION (PCHQK3-S, DPCHQ3-D)
  !***KEYWORDS  PCHIP INTERPOLATOR QUICK CHECK
  !***AUTHOR  Fritsch, F. N., (LLNL)
  !***DESCRIPTION
  !
  !             DPCHIP QUICK CHECK NUMBER 3
  !
  !     TESTS THE INTERPOLATORS:  DPCHIC, DPCHIM, DPCHSP.
  ! *Usage:
  !
  !        INTEGER  LUN, KPRINT, IPASS
  !
  !        CALL DPCHQ3 (LUN, KPRINT, IPASS)
  !
  ! *Arguments:
  !
  !     LUN   :IN  is the unit number to which output is to be written.
  !
  !     KPRINT:IN  controls the amount of output, as specified in the
  !                SLATEC Guidelines.
  !
  !     IPASS:OUT  will contain a pass/fail flag.  IPASS=1 is good.
  !                IPASS=0 indicates one or more tests failed.
  !
  ! *Description:
  !
  !   This routine interpolates a constructed data set with all three
  !   DPCHIP interpolators and compares the results with those obtained
  !   on a Cray X/MP.  Two different values of the DPCHIC parameter SWITCH
  !   are used.
  !
  ! *Remarks:
  !     1. The Cray results are given only to nine significant figures,
  !        so don't expect them to match to more.
  !     2. The results will depend to some extent on the accuracy of
  !        the EXP function.
  !
  !***ROUTINES CALLED  COMP, D1MACH, DPCHIC, DPCHIM, DPCHSP
  !***REVISION HISTORY  (YYMMDD)
  !   900309  DATE WRITTEN
  !   900314  Converted to a subroutine and added a SLATEC 4.0 prologue.
  !   900315  Revised prologue and improved some output formats.  (FNF)
  !   900316  Made TOLD machine-dependent and added extra output when
  !           KPRINT=3.  (FNF)
  !   900320  Added E0's to DATA statement for X to reduce single/double
  !           differences, and other minor cosmetic changes.
  !   900320  Converted to double precision.
  !   900321  Removed IFAIL from call sequence for SLATEC standards and
  !           made miscellaneous cosmetic changes.  (FNF)
  !   900322  Minor changes to reduce single/double differences.  (FNF)
  !   900530  Tolerance (TOLD) and argument to DPCHIC changed.  (WRB)
  !   900802  Modified TOLD formula and constants in DPCHIC calls to
  !           correct DPCHQ3 failures.  (FNF)
  !   901130  Several significant changes:  (FNF)
  !           1. Changed comparison between DPCHIM and DPCHIC to only
  !              require agreement to machine precision.
  !           2. Revised to print more output when KPRINT=3.
  !           3. Added 1P's to formats.
  !   910708  Minor modifications in use of KPRINT.  (WRB)
  !   930317  Improved output formats.  (FNF)
  !***END PROLOGUE  DPCHQ3
  !
  !*Internal Notes:
  !
  !     TOLD is used to compare with stored Cray results.  Its value
  !          should be consistent with significance of stored values.
  !     TOLZ is used for cases in which exact equality is expected.
  !     TOL  is used for cases in which agreement to machine precision
  !          is expected.
  !**End
  !
  !  Declare arguments.
  !
  INTEGER Lun , Kprint , Ipass
  LOGICAL COMP
  DOUBLE PRECISION D1MACH
  !
  !  Declare variables.
  !
  INTEGER i , ic(2) , ierr , ifail , N , nbad , nbadz , NWK
  PARAMETER (N=9,NWK=2*N)
  DOUBLE PRECISION d(N) , dc(N) , dc5 , dc6 , dm(N) , ds(N) , err , f(N) ,&
    MONE , tol , told , tolz , vc(2) , x(N) , wk(NWK) , ZERO
  PARAMETER (ZERO=0.0D0,MONE=-1.0D0)
  CHARACTER(6) :: result
  !
  !  Initialize.
  !
  !       Data.
  DATA ic/0 , 0/
  DATA x/ - 2.2D0 , -1.2D0 , -1.0D0 , -0.5D0 , -0.01D0 , 0.5D0 , 1.0D0 ,&
    2.0D0 , 2.2D0/
  !
  !       Results generated on Cray X/MP (9 sign. figs.)
  DATA dm/0. , 3.80027352D-01 , 7.17253009D-01 , 5.82014161D-01 , 0. ,&
    -5.68208031D-01 , -5.13501618D-01 , -7.77910977D-02 ,&
    -2.45611117D-03/
  DATA dc5 , dc6/1.76950158D-02 , -5.69579814D-01/
  DATA ds/ - 5.16830792D-02 , 5.71455855D-01 , 7.40530225D-01 ,&
    7.63864934D-01 , 1.92614386D-02 , -7.65324380D-01 , -7.28209035D-01 ,&
    -7.98445427D-02 , -2.85983446D-02/
  !
  !***FIRST EXECUTABLE STATEMENT  DPCHQ3
  ifail = 0
  !
  !        Set tolerances.
  tol = 10*D1MACH(4)
  told = MAX(1.0D-7,10*tol)
  tolz = ZERO
  !
  IF ( Kprint>=3 ) WRITE (Lun,99001)
  !
  !  FORMATS.
  !
  99001 FORMAT ('1'//10X,'TEST DPCHIP INTERPOLATORS')
  IF ( Kprint>=2 ) WRITE (Lun,99002)
  99002 FORMAT (//10X,'DPCHQ3 RESULTS'/10X,'--------------')
  !
  !  Set up data.
  !
  DO i = 1 , N
    f(i) = EXP(-x(i)**2)
  ENDDO
  !
  IF ( Kprint>=3 ) THEN
    WRITE (Lun,99003)
    99003   FORMAT (//5X,'DATA:'/39X,'---------- EXPECTED D-VALUES ----------'/12X,&
      'X',9X,'F',18X,'DM',13X,'DC',13X,'DS')
    DO i = 1 , 4
      WRITE (Lun,99009) x(i) , f(i) , dm(i) , ds(i)
    ENDDO
    WRITE (Lun,99010) x(5) , f(5) , dm(5) , dc5 , ds(5)
    WRITE (Lun,99010) x(6) , f(6) , dm(6) , dc6 , ds(6)
    DO i = 7 , N
      WRITE (Lun,99009) x(i) , f(i) , dm(i) , ds(i)
    ENDDO
  ENDIF
  !
  !  Test DPCHIM.
  !
  IF ( Kprint>=3 ) WRITE (Lun,99011) 'IM'
  !     --------------------------------
  CALL DPCHIM(N,x,f,d,1,ierr)
  !     --------------------------------
  !        Expect IERR=1 (one monotonicity switch).
  IF ( Kprint>=3 ) WRITE (Lun,99012) 1
  IF ( .NOT.COMP(ierr,1,Lun,Kprint) ) THEN
    ifail = ifail + 1
  ELSE
    IF ( Kprint>=3 ) WRITE (Lun,99013)
    nbad = 0
    nbadz = 0
    DO i = 1 , N
      result = '  OK'
      !             D-values should agree with stored values.
      !               (Zero values should agree exactly.)
      IF ( dm(i)==ZERO ) THEN
        err = ABS(d(i))
        IF ( err>tolz ) THEN
          nbadz = nbadz + 1
          result = '**BADZ'
        ENDIF
      ELSE
        err = ABS((d(i)-dm(i))/dm(i))
        IF ( err>told ) THEN
          nbad = nbad + 1
          result = '**BAD'
        ENDIF
      ENDIF
      IF ( Kprint>=3 ) WRITE (Lun,99014) i , x(i) , d(i) , err , result
    ENDDO
    IF ( (nbadz/=0).OR.(nbad/=0) ) THEN
      ifail = ifail + 1
      IF ( (nbadz/=0).AND.(Kprint>=2) ) WRITE (Lun,99004) nbad
      99004     FORMAT (/'    **',I5,' DPCHIM RESULTS FAILED TO BE EXACTLY ZERO.')
      IF ( (nbad/=0).AND.(Kprint>=2) ) WRITE (Lun,99015) nbad , 'IM' , told
    ELSE
      IF ( Kprint>=2 ) WRITE (Lun,99016) 'IM'
    ENDIF
  ENDIF
  !
  !  Test DPCHIC -- options set to reproduce DPCHIM.
  !
  IF ( Kprint>=3 ) WRITE (Lun,99011) 'IC'
  !     --------------------------------------------------------
  CALL DPCHIC(ic,vc,ZERO,N,x,f,dc,1,wk,NWK,ierr)
  !     --------------------------------------------------------
  !        Expect IERR=0 .
  IF ( Kprint>=3 ) WRITE (Lun,99012) 0
  IF ( .NOT.COMP(ierr,0,Lun,Kprint) ) THEN
    ifail = ifail + 1
  ELSE
    IF ( Kprint>=3 ) WRITE (Lun,99013)
    nbad = 0
    DO i = 1 , N
      result = '  OK'
      !           D-values should agree exactly with those computed by DPCHIM.
      !            (To be generous, will only test to machine precision.)
      err = ABS(d(i)-dc(i))
      IF ( err>tol ) THEN
        nbad = nbad + 1
        result = '**BAD'
      ENDIF
      IF ( Kprint>=3 ) WRITE (Lun,99014) i , x(i) , dc(i) , err , result
    ENDDO
    IF ( nbad/=0 ) THEN
      ifail = ifail + 1
      IF ( Kprint>=2 ) WRITE (Lun,99015) nbad , 'IC' , tol
    ELSE
      IF ( Kprint>=2 ) WRITE (Lun,99016) 'IC'
    ENDIF
  ENDIF
  !
  !  Test DPCHIC -- default nonzero switch derivatives.
  !
  IF ( Kprint>=3 ) WRITE (Lun,99011) 'IC'
  !     -------------------------------------------------------
  CALL DPCHIC(ic,vc,MONE,N,x,f,d,1,wk,NWK,ierr)
  !     -------------------------------------------------------
  !        Expect IERR=0 .
  IF ( Kprint>=3 ) WRITE (Lun,99012) 0
  IF ( .NOT.COMP(ierr,0,Lun,Kprint) ) THEN
    ifail = ifail + 1
  ELSE
    IF ( Kprint>=3 ) WRITE (Lun,99013)
    nbad = 0
    nbadz = 0
    DO i = 1 , N
      result = '  OK'
      !            D-values should agree exactly with those computed in
      !            previous call, except at points 5 and 6.
      IF ( (i<5).OR.(i>6) ) THEN
        err = ABS(d(i)-dc(i))
        IF ( err>tolz ) THEN
          nbadz = nbadz + 1
          result = '**BADA'
        ENDIF
      ELSE
        IF ( i==5 ) THEN
          err = ABS((d(i)-dc5)/dc5)
        ELSE
          err = ABS((d(i)-dc6)/dc6)
        ENDIF
        IF ( err>told ) THEN
          nbad = nbad + 1
          result = '**BAD'
        ENDIF
      ENDIF
      IF ( Kprint>=3 ) WRITE (Lun,99014) i , x(i) , d(i) , err , result
    ENDDO
    IF ( (nbadz/=0).OR.(nbad/=0) ) THEN
      ifail = ifail + 1
      IF ( (nbadz/=0).AND.(Kprint>=2) ) WRITE (Lun,99005) nbad
      99005     FORMAT (/'    **',I5,' DPCHIC RESULTS FAILED TO AGREE WITH',&
        ' PREVIOUS CALL.')
      IF ( (nbad/=0).AND.(Kprint>=2) ) WRITE (Lun,99015) nbad , 'IC' , told
    ELSE
      IF ( Kprint>=2 ) WRITE (Lun,99016) 'IC'
    ENDIF
  ENDIF
  !
  !  Test DPCHSP.
  !
  IF ( Kprint>=3 ) WRITE (Lun,99011) 'SP'
  !     -------------------------------------------------
  CALL DPCHSP(ic,vc,N,x,f,d,1,wk,NWK,ierr)
  !     -------------------------------------------------
  !        Expect IERR=0 .
  IF ( Kprint>=3 ) WRITE (Lun,99012) 0
  IF ( .NOT.COMP(ierr,0,Lun,Kprint) ) THEN
    ifail = ifail + 1
  ELSE
    IF ( Kprint>=3 ) WRITE (Lun,99013)
    nbad = 0
    DO i = 1 , N
      result = '  OK'
      !             D-values should agree with stored values.
      err = ABS((d(i)-ds(i))/ds(i))
      IF ( err>told ) THEN
        nbad = nbad + 1
        result = '**BAD'
      ENDIF
      IF ( Kprint>=3 ) WRITE (Lun,99014) i , x(i) , d(i) , err , result
    ENDDO
    IF ( nbad/=0 ) THEN
      ifail = ifail + 1
      IF ( Kprint>=2 ) WRITE (Lun,99015) nbad , 'SP' , told
    ELSE
      IF ( Kprint>=2 ) WRITE (Lun,99016) 'SP'
    ENDIF
  ENDIF
  !
  !  PRINT SUMMARY AND TERMINATE.
  !
  IF ( (Kprint>=2).AND.(ifail/=0) ) WRITE (Lun,99006) ifail
  99006 FORMAT (/' *** TROUBLE ***',I5,' INTERPOLATION TESTS FAILED.')
  !
  IF ( ifail==0 ) THEN
    Ipass = 1
    IF ( Kprint>=2 ) WRITE (Lun,99007)
    99007   FORMAT (/' ------------ DPCHIP PASSED  ALL INTERPOLATION TESTS',&
      ' ------------')
  ELSE
    Ipass = 0
    IF ( Kprint>=1 ) WRITE (Lun,99008)
    99008   FORMAT (/' ************ DPCHIP FAILED SOME INTERPOLATION TESTS',&
      ' ************')
  ENDIF
  !
  RETURN
  99009 FORMAT (5X,F10.2,1P,D15.5,4X,D15.5,15X,D15.5)
  99010 FORMAT (5X,F10.2,1P,D15.5,4X,3D15.5)
  99011 FORMAT (/5X,'DPCH',A2,' TEST:')
  99012 FORMAT (15X,'EXPECT  IERR =',I5)
  99013 FORMAT (/9X,'I',7X,'X',9X,'D',13X,'ERR')
  99014 FORMAT (5X,I5,F10.2,1P,2D15.5,2X,A)
  99015 FORMAT (/'    **',I5,' DPCH',A2,' RESULTS FAILED TOLERANCE TEST.',&
    '  TOL =',1P,D10.3)
  99016 FORMAT (/5X,'  ALL DPCH',A2,' RESULTS OK.')
  !------------- LAST LINE OF DPCHQ3 FOLLOWS -----------------------------
END SUBROUTINE DPCHQ3
