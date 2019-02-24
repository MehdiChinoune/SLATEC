MODULE TEST16_MOD
  IMPLICIT NONE

CONTAINS
  !DECK DQC36J
  SUBROUTINE DQC36J(Lun,Kprint,Ipass)
    IMPLICIT NONE
    !***BEGIN PROLOGUE  DQC36J
    !***SUBSIDIARY
    !***PURPOSE  THIS IS A QUICK CHECK PROGRAM FOR THE SUBROUTINES DRC3JJ,
    !            DRC3JM, AND DRC6J, WHICH CALCULATE THE WIGNER COEFFICIENTS,
    !            3J AND 6J.
    !***LIBRARY   SLATEC
    !***CATEGORY  C19
    !***TYPE      DOUBLE PRECISION (QC36J-S, DQC36J-D)
    !***KEYWORDS  3J COEFFICIENTS, 3J SYMBOLS, 6J COEFFICIENTS, 6J SYMBOLS,
    !             CLEBSCH-GORDAN COEFFICIENTS, QUICK CHECK,
    !             RACAH COEFFICIENTS, VECTOR ADDITION COEFFICIENTS,
    !             WIGNER COEFFICIENTS
    !***AUTHOR  LOZIER, DANIEL W., (NIST)
    !           MCCLAIN, MARJORIE A., (NIST)
    !           SMITH, JOHN M., (NIST AND GEORGE MASON UNIVERSITY)
    !***REFERENCES  MESSIAH, ALBERT., QUANTUM MECHANICS, VOLUME II,
    !               NORTH-HOLLAND PUBLISHING COMPANY, 1963.
    !***ROUTINES CALLED  D1MACH, DRC3JJ, DRC3JM, DRC6J, NUMXER, XERCLR,
    !                     XSETF
    !***REVISION HISTORY  (YYMMDD)
    !   891129  DATE WRITTEN
    !   910415  Mixed type expressions eliminated; precision of output
    !           formats made uniform for all tests; detail added to output
    !           when KPRINT=2 and a test fails; name of quick check added
    !           to output when KPRINT=3 or KPRINT=2 and a test fails; some
    !           output formats modified for clarity or adherence to SLATEC
    !           guidelines. These changes were done by D. W. Lozier.
    !   930115  Replaced direct calculation of 3j-6j symbols in tests 1, 2,
    !           and 4 with values stored in data statements.  This involved
    !           removing all calls to subroutine DRACAH.  These changes were
    !           made by M. McClain.
    !***END PROLOGUE  DQC36J
    !
    INTEGER Lun, Kprint, Ipass
    !
    CHARACTER string*36, fmt*30, fmt2*13
    INTEGER ipass1, ipass2, ipass3, ipass4, ipass5, NDIM, ier, index, &
      i, first, last, nsig, NUMXER, nerr, ierjj, ierjm
    PARAMETER (NDIM=15)
    REAL(8) :: tol, l1, l2, l3, m1, m2, m3, l1min, l1max, &
      m2min, m2max, diff(NDIM), D1MACH, x, jjval, jmval, &
      thrcof(NDIM), sixcof(NDIM), r3jj(8), r3jm(14), &
      r6j(15)
    !
    DATA r3jj/2.7888667551135851599272400859506249646427D-1, &
      -9.5346258924559231544677592152721599861388D-2, &
      -6.7419986246324208624649067643642846008908D-2, &
      1.5331103516796664129653641995122585402552D-1, &
      -1.5644655469368596972508355755184909201031D-1, &
      1.0994504121565505107947893271429777505797D-1, &
      -5.5362356931317194333395729256559987745156D-2, &
      1.7998354511377858329814092962590761537262D-2/
    !
    DATA r3jm/2.0915897328861524261384476677886072045904D-2, &
      8.5375655532152472212727551895879778672762D-2, &
      9.0829537086869251694343772675676175068677D-2, &
      -3.8905437784649939169989036459327065796765D-2, &
      -6.6373497016568063569146501153397525003444D-2, &
      6.4952404052838939503061387831391216401903D-2, &
      2.1589431059540375939250708046202926313913D-2, &
      -7.7891271178523921999229618972588887261359D-2, &
      3.5976437105954340188005810512211794384411D-2, &
      5.4730150002126342307937096038252488407360D-2, &
      -7.5967866595676151462927617736745078548338D-2, &
      -2.1922444553989211377558215380002910257762D-2, &
      1.0116774428077220242411199686231560525497D-1, &
      7.3482572624471970469595137204530687381176D-2/
    !
    DATA r6j/3.4909051383732997774596981092927782159095D-2, &
      -3.7430250396597916085929064401358002747549D-2, &
      1.8908663909595601841537964135129184202064D-2, &
      7.3424482549286434570947151839589351100581D-3, &
      -2.3589351850817944585847816357296508528608D-2, &
      1.9134769552154365200026782557432864615918D-2, &
      1.2880173977241722084434864685591278730958D-3, &
      -1.9300183662905265397749119277519305417805D-2, &
      1.6773059493828887697413611251392749162229D-2, &
      5.5011472748509487167380502058890639729979D-3, &
      -2.1354397908968309742136976853078409839580D-2, &
      3.4603644514353873082775312319159137043869D-3, &
      2.5209500547955845860442730268272167527589D-2, &
      1.4839905612217133028540464232557124565509D-2, &
      2.7085776806331855972407001825016114677027D-3/
    !
    !***FIRST EXECUTABLE STATEMENT  DQC36J
    !
    ! --- INITIALIZATION OF TESTS
    tol = 100.0D0*D1MACH(3)
    IF ( Kprint>=2 ) THEN
      WRITE (Lun,*) ' THIS IS DQC36J, A TEST PROGRAM FOR THE '//&
        'DOUBLE PRECISION 3J6J PACKAGE.'
      WRITE (Lun,*) ' AN EXPLANATION OF THE VARIOUS '//&
        'TESTS CAN BE FOUND IN THE PROGRAM COMMENTS.'
      WRITE (Lun,*)
    ENDIF
    !
    ! --- FIND NUMBER OF SIGNIFICANT FIGURES FOR FORMATTING
    x = 1.D0/3.D0
    WRITE (string,99001) x
    99001 FORMAT (F35.25)
    DO i = 1, 35
      IF ( string(i:i)=='3' ) THEN
        first = i
        EXIT
      ENDIF
    ENDDO
    DO i = first, 35
      IF ( string(i:i)/='3' ) THEN
        last = i - 1
        GOTO 100
      ENDIF
    ENDDO
    last = 36
    100  nsig = last - first + 1
    fmt(1:16) = '(1X,F5.1,T8,G35.'
    WRITE (fmt(17:18),'(I2)') nsig
    fmt(19:27) = ',T45,G35.'
    WRITE (fmt(28:29),'(I2)') nsig
    fmt(30:30) = ')'
    fmt2(1:10) = '(1X,A,G35.'
    WRITE (fmt2(11:12),'(I2)') nsig
    fmt2(13:13) = ')'
    !
    ! --- TEST 1: COMPARE DRC3JJ VALUES WITH FORMULA
    ipass1 = 1
    l2 = 4.5D0
    l3 = 3.5D0
    m2 = -3.5D0
    m3 = 2.5D0
    CALL DRC3JJ(l2,l3,m2,m3,l1min,l1max,thrcof,NDIM,ier)
    IF ( ier/=0 ) THEN
      ipass1 = 0
    ELSE
      DO index = 1, INT(l1max-l1min) + 1
        m1 = 1.0D0
        diff(index) = ABS(thrcof(index)-r3jj(index))
        IF ( diff(index)>ABS(r3jj(index))*tol ) ipass1 = 0
      ENDDO
    ENDIF
    IF ( Kprint>=3.OR.(Kprint==2.AND.ipass1==0) ) THEN
      WRITE (Lun,*) ' TEST 1, RECURRENCE IN L1, COMPARE VALUES OF 3J ', &
        'CALCULATED BY DRC3JJ TO'
      WRITE (Lun,*) ' VALUES CALCULATED BY EXPLICIT FORMULA FROM ', &
        'MESSIAH''S QUANTUM MECHANICS'
      WRITE (Lun,99009) l2, l3
      WRITE (Lun,99010) m1, m2, m3
      IF ( ier/=0 ) THEN
        WRITE (Lun,*) ' ERROR RETURNED FROM SUBROUTINE ', 'DRC3JJ: IER =', &
          ier
      ELSE
        WRITE (Lun,99002)
        99002     FORMAT ('    L1',T31,'DRC3JJ VALUE',T67,'FORMULA VALUE')
        DO index = 1, INT(l1max-l1min) + 1
          l1 = index+ l1min - 1
          WRITE (Lun,fmt) l1, thrcof(index), r3jj(index)
          IF ( diff(index)>ABS(r3jj(index))*tol ) WRITE (Lun,'(1X,A,F5.1)')&
            'DIFFERENCE EXCEEDS ERROR '//'TOLERANCE FOR L1 =', l1
        ENDDO
      ENDIF
    ENDIF
    IF ( ipass1==0 ) THEN
      IF ( Kprint>=1 ) THEN
        WRITE (Lun,*) ' ***** ***** TEST 1 FAILED ***** *****'
        WRITE (Lun,*)
      ENDIF
    ELSEIF ( Kprint>=2 ) THEN
      WRITE (Lun,*) ' ***** ***** TEST 1 PASSED ***** *****'
      WRITE (Lun,*)
    ENDIF
    !
    ! --- TEST 2: COMPARE DRC3JM VALUES WITH FORMULA
    ipass2 = 1
    l1 = 8.0D0
    l2 = 7.5D0
    l3 = 6.5D0
    m1 = 1.0D0
    CALL DRC3JM(l1,l2,l3,m1,m2min,m2max,thrcof,NDIM,ier)
    IF ( ier/=0 ) THEN
      ipass2 = 0
    ELSE
      DO index = 1, INT(m2max-m2min) + 1
        m2 = index+ m2min - 1
        m3 = -m1 - m2
        diff(index) = ABS(thrcof(index)-r3jm(index))
        IF ( diff(index)>ABS(r3jm(index))*tol ) ipass2 = 0
      ENDDO
    ENDIF
    IF ( Kprint>=3.OR.(Kprint==2.AND.ipass2==0) ) THEN
      WRITE (Lun,*) ' TEST 2, RECURRENCE IN M2, COMPARE VALUES OF 3J ', &
        'CALCULATED BY DRC3JM TO'
      WRITE (Lun,*) ' VALUES CALCULATED BY EXPLICIT FORMULA FROM ', &
        'MESSIAH''S QUANTUM MECHANICS'
      WRITE (Lun,99003) l1, l2, l3
      99003   FORMAT (' L1 = ',F5.1,'   L2 = ',F5.1,'   L3 = ',F5.1)
      WRITE (Lun,99004) m1
      99004   FORMAT (' M1 = ',F5.1,'                M3 = -(M1+M2)')
      IF ( ier/=0 ) THEN
        WRITE (Lun,*) ' ERROR RETURNED FROM SUBROUTINE ', 'DRC3JM: IER =', &
          ier
      ELSE
        WRITE (Lun,99005)
        99005     FORMAT ('    M2',T31,'DRC3JM VALUE',T67,'FORMULA VALUE')
        DO index = 1, INT(m2max-m2min) + 1
          m2 = index+ m2min - 1
          WRITE (Lun,fmt) m2, thrcof(index), r3jm(index)
          IF ( diff(index)>ABS(r3jm(index))*tol ) WRITE (Lun,'(1X,A,F5.1)')&
            'DIFFERENCE EXCEEDS ERROR '//'TOLERANCE FOR M2 =', m2
        ENDDO
      ENDIF
    ENDIF
    IF ( ipass2==0 ) THEN
      IF ( Kprint>=1 ) THEN
        WRITE (Lun,*) ' ***** ***** TEST 2 FAILED ***** *****'
        WRITE (Lun,*)
      ENDIF
    ELSEIF ( Kprint>=2 ) THEN
      WRITE (Lun,*) ' ***** ***** TEST 2 PASSED ***** *****'
      WRITE (Lun,*)
    ENDIF
    !
    ! --- TEST3: COMPARE COMMON VALUE OF DRC3JJ AND DRC3JM
    ipass3 = 1
    l1 = 100.0D0
    l2 = 2.0D0
    l3 = 100.0D0
    m1 = -10.0D0
    m2 = 0.0D0
    m3 = 10.0D0
    CALL DRC3JJ(l2,l3,m2,m3,l1min,l1max,thrcof,NDIM,ierjj)
    jjval = thrcof(3)
    CALL DRC3JM(l1,l2,l3,m1,m2min,m2max,thrcof,NDIM,ierjm)
    jmval = thrcof(3)
    IF ( ierjj/=0.OR.ierjm/=0 ) THEN
      ipass3 = 0
    ELSE
      diff(1) = ABS(jjval-jmval)
      IF ( diff(1)>0.5*ABS(jjval+jmval)*tol ) ipass3 = 0
    ENDIF
    IF ( Kprint>=3.OR.(Kprint==2.AND.ipass3==0) ) THEN
      WRITE (Lun,*) ' TEST 3, COMPARE A COMMON VALUE CALCULATED BY ', &
        'BOTH DRC3JJ AND DRC3JM'
      WRITE (Lun,*) ' L1 = 100.0   L2 =   2.0   L3 = 100.0'
      WRITE (Lun,*) ' M1 = -10.0   M2 =   0.0   M3 =  10.0'
      IF ( ierjj/=0 ) THEN
        WRITE (Lun,*) ' ERROR RETURNED FROM SUBROUTINE ', 'DRC3JJ: IER =', &
          ierjj
      ELSEIF ( ierjm/=0 ) THEN
        WRITE (Lun,*) ' ERROR RETURNED FROM SUBROUTINE ', 'DRC3JM: IER =', &
          ierjm
      ELSE
        WRITE (Lun,fmt2) 'DRC3JJ VALUE =', jjval
        WRITE (Lun,fmt2) 'DRC3JM VALUE =', jmval
        IF ( diff(1)>0.5*ABS(jjval+jmval)*tol ) WRITE (Lun,'(1X,A)')&
          'DIFFERENCE EXCEEDS ERROR TOLERANCE'
      ENDIF
    ENDIF
    IF ( ipass3==0 ) THEN
      IF ( Kprint>=1 ) THEN
        WRITE (Lun,*) ' ***** ***** TEST 3 FAILED ***** *****'
        WRITE (Lun,*)
      ENDIF
    ELSEIF ( Kprint>=2 ) THEN
      WRITE (Lun,*) ' ***** ***** TEST 3 PASSED ***** *****'
      WRITE (Lun,*)
    ENDIF
    !
    ! --- TEST 4: COMPARE DRC6J VALUES WITH FORMULA
    ipass4 = 1
    l2 = 8.0D0
    l3 = 7.0D0
    m1 = 6.5D0
    m2 = 7.5D0
    m3 = 7.5D0
    CALL DRC6J(l2,l3,m1,m2,m3,l1min,l1max,sixcof,NDIM,ier)
    IF ( ier/=0 ) THEN
      ipass4 = 0
    ELSE
      DO index = 1, INT(l1max-l1min) + 1
        diff(index) = ABS(sixcof(index)-r6j(index))
        IF ( diff(index)>ABS(r6j(index))*tol ) ipass4 = 0
      ENDDO
    ENDIF
    IF ( Kprint>=3.OR.(Kprint==2.AND.ipass4==0) ) THEN
      WRITE (Lun,*) ' TEST 4, RECURRENCE IN L1, COMPARE VALUES OF 6J ', &
        'CALCULATED BY DRC6J TO'
      WRITE (Lun,*) ' VALUES CALCULATED BY EXPLICIT FORMULA FROM ', &
        'MESSIAH''S QUANTUM MECHANICS'
      WRITE (Lun,99009) l2, l3
      WRITE (Lun,99010) m1, m2, m3
      IF ( ier/=0 ) THEN
        WRITE (Lun,*) ' ERROR RETURNED FROM SUBROUTINE ', 'DRC6J: IER =', &
          ier
      ELSE
        WRITE (Lun,99006)
        99006     FORMAT ('    L1',T32,'DRC6J VALUE',T67,'FORMULA VALUE')
        DO index = 1, INT(l1max-l1min) + 1
          l1 = index+ l1min - 1
          WRITE (Lun,fmt) l1, sixcof(index), r6j(index)
          IF ( diff(index)>ABS(r6j(index))*tol ) WRITE (Lun,'(1X,A,F5.1)')&
            'DIFFERENCE EXCEEDS ERROR '//'TOLERANCE FOR L1 =', l1
        ENDDO
      ENDIF
    ENDIF
    IF ( ipass4==0 ) THEN
      IF ( Kprint>=1 ) THEN
        WRITE (Lun,*) ' ***** ***** TEST 4 FAILED ***** *****'
        WRITE (Lun,*)
      ENDIF
    ELSEIF ( Kprint>=2 ) THEN
      WRITE (Lun,*) ' ***** ***** TEST 4 PASSED ***** *****'
      WRITE (Lun,*)
    ENDIF
    !
    ! --- TEST 5: CHECK INVALID INPUT
    ipass5 = 1
    IF ( Kprint<=2 ) THEN
      CALL XSETF(0)
    ELSE
      CALL XSETF(-1)
    ENDIF
    IF ( Kprint>=3 ) WRITE (Lun,*) ' TEST 5, CHECK FOR PROPER HANDLING ', &
      'OF INVALID INPUT'
    ! --- DRC3JJ: L2-ABS(M2) OR L3-ABS(M3) LESS THAN ZERO (IER=1)
    l2 = 2.0D0
    l3 = 100.0D0
    m1 = -6.0D0
    m2 = -4.0D0
    m3 = 10.0D0
    IF ( Kprint>=3 ) WRITE (Lun,*)
    CALL XERCLR
    CALL DRC3JJ(l2,l3,m2,m3,l1min,l1max,thrcof,NDIM,ier)
    IF ( NUMXER(nerr)/=ier ) ipass5 = 0
    ! --- DRC3JJ: L2+ABS(M2) OR L3+ABS(M3) NOT INTEGER (IER=2)
    l2 = 2.0D0
    l3 = 99.5D0
    m1 = -10.0D0
    m2 = 0.0D0
    m3 = 10.0D0
    IF ( Kprint>=3 ) WRITE (Lun,*)
    CALL XERCLR
    CALL DRC3JJ(l2,l3,m2,m3,l1min,l1max,thrcof,NDIM,ier)
    IF ( NUMXER(nerr)/=ier ) ipass5 = 0
    ! --- DRC3JJ: L1MAX-L1MIN NOT INTEGER (IER=3)
    l2 = 3.2D0
    l3 = 4.5D0
    m1 = -1.3D0
    m2 = 0.8D0
    m3 = 0.5D0
    IF ( Kprint>=3 ) WRITE (Lun,*)
    CALL XERCLR
    CALL DRC3JJ(l2,l3,m2,m3,l1min,l1max,thrcof,NDIM,ier)
    IF ( NUMXER(nerr)/=ier ) ipass5 = 0
    ! --- DRC3JJ: L1MIN GREATER THAN L1MAX (IER=4)
    !             (NO TEST -- THIS ERROR SHOULD NEVER OCCUR)
    ! --- DRC3JJ: DIMENSION OF THRCOF TOO SMALL (IER=5)
    l2 = 10.0D0
    l3 = 150.0D0
    m1 = -10.0D0
    m2 = 0.0D0
    m3 = 10.0D0
    IF ( Kprint>=3 ) WRITE (Lun,*)
    CALL XERCLR
    CALL DRC3JJ(l2,l3,m2,m3,l1min,l1max,thrcof,NDIM,ier)
    IF ( NUMXER(nerr)/=ier ) ipass5 = 0
    ! --- DRC3JM: L1-ABS(M1) .LT. ZERO OR L1+ABS(M1) NOT INTEGER (IER=1)
    l1 = 100.0D0
    l2 = 2.0D0
    l3 = 100.0D0
    m1 = 150.0D0
    IF ( Kprint>=3 ) WRITE (Lun,*)
    CALL XERCLR
    CALL DRC3JM(l1,l2,l3,m1,m2min,m2max,thrcof,NDIM,ier)
    IF ( NUMXER(nerr)/=ier ) ipass5 = 0
    ! --- DRC3JM: L1, L2, L3 DO NOT SATISFY TRIANGULAR CONDITION (IER=2)
    l1 = 20.0D0
    l2 = 5.0D0
    l3 = 10.0D0
    m1 = -10.0D0
    IF ( Kprint>=3 ) WRITE (Lun,*)
    CALL XERCLR
    CALL DRC3JM(l1,l2,l3,m1,m2min,m2max,thrcof,NDIM,ier)
    IF ( NUMXER(nerr)/=ier ) ipass5 = 0
    ! --- DRC3JM: L1+L2+L3 NOT INTEGER (IER=3)
    l1 = 1.0D0
    l2 = 1.3D0
    l3 = 1.5D0
    m1 = 0.0D0
    IF ( Kprint>=3 ) WRITE (Lun,*)
    CALL XERCLR
    CALL DRC3JM(l1,l2,l3,m1,m2min,m2max,thrcof,NDIM,ier)
    IF ( NUMXER(nerr)/=ier ) ipass5 = 0
    ! --- DRC3JM: M2MAX-M2MIN NOT INTEGER (IER=4)
    l1 = 1.0D0
    l2 = 1.3D0
    l3 = 1.7D0
    m1 = 0.0D0
    IF ( Kprint>=3 ) WRITE (Lun,*)
    CALL XERCLR
    CALL DRC3JM(l1,l2,l3,m1,m2min,m2max,thrcof,NDIM,ier)
    IF ( NUMXER(nerr)/=ier ) ipass5 = 0
    ! --- DRC3JM: M2MIN GREATER THAN M2MAX (IER=5)
    !             (NO TEST -- THIS ERROR SHOULD NEVER OCCUR)
    ! --- DRC3JM: DIMENSION OF THRCOF TOO SMALL (IER=6)
    l1 = 100.0D0
    l2 = 10.0D0
    l3 = 110.0D0
    m1 = -10.0D0
    IF ( Kprint>=3 ) WRITE (Lun,*)
    CALL XERCLR
    CALL DRC3JM(l1,l2,l3,m1,m2min,m2max,thrcof,NDIM,ier)
    IF ( NUMXER(nerr)/=ier ) ipass5 = 0
    ! --- DRC6J: L2+L3+L5+L6 OR L4+L2+L6 NOT INTEGER (IER=1)
    l2 = 0.5D0
    l3 = 1.0D0
    m1 = 0.5D0
    m2 = 2.0D0
    m3 = 3.0D0
    IF ( Kprint>=3 ) WRITE (Lun,*)
    CALL XERCLR
    CALL DRC6J(l2,l3,m1,m2,m3,l1min,l1max,sixcof,NDIM,ier)
    IF ( NUMXER(nerr)/=ier ) ipass5 = 0
    ! --- DRC6J: L4, L2, L6 TRIANGULAR CONDITION NOT SATISFIED (IER=2)
    l2 = 1.0D0
    l3 = 3.0D0
    m1 = 5.0D0
    m2 = 6.0D0
    m3 = 2.0D0
    IF ( Kprint>=3 ) WRITE (Lun,*)
    CALL XERCLR
    CALL DRC6J(l2,l3,m1,m2,m3,l1min,l1max,sixcof,NDIM,ier)
    IF ( NUMXER(nerr)/=ier ) ipass5 = 0
    ! --- DRC6J: L4, L5, L3 TRIANGULAR CONDITION NOT SATISFIED (IER=3)
    l2 = 4.0D0
    l3 = 1.0D0
    m1 = 5.0D0
    m2 = 3.0D0
    m3 = 2.0D0
    IF ( Kprint>=3 ) WRITE (Lun,*)
    CALL XERCLR
    CALL DRC6J(l2,l3,m1,m2,m3,l1min,l1max,sixcof,NDIM,ier)
    IF ( NUMXER(nerr)/=ier ) ipass5 = 0
    ! --- DRC6J: L1MAX-L1MIN NOT INTEGER (IER=4)
    l2 = 0.9D0
    l3 = 0.5D0
    m1 = 0.9D0
    m2 = 0.4D0
    m3 = 0.2D0
    IF ( Kprint>=3 ) WRITE (Lun,*)
    CALL XERCLR
    CALL DRC6J(l2,l3,m1,m2,m3,l1min,l1max,sixcof,NDIM,ier)
    IF ( NUMXER(nerr)/=ier ) ipass5 = 0
    ! --- DRC6J: L1MIN GREATER THAN L1MAX (IER=5)
    !            (NO TEST -- THIS ERROR SHOULD NEVER OCCUR)
    ! --- DRC6J: DIMENSION OF SIXCOF TOO SMALL (IER=6)
    l2 = 50.0D0
    l3 = 25.0D0
    m1 = 15.0D0
    m2 = 30.0D0
    m3 = 40.0D0
    IF ( Kprint>=3 ) WRITE (Lun,*)
    CALL XERCLR
    CALL DRC6J(l2,l3,m1,m2,m3,l1min,l1max,sixcof,NDIM,ier)
    IF ( NUMXER(nerr)/=ier ) ipass5 = 0
    IF ( ipass5==0 ) THEN
      IF ( Kprint>=1 ) THEN
        WRITE (Lun,*) ' ***** ***** TEST 5 FAILED ***** *****'
        WRITE (Lun,*)
      ENDIF
    ELSEIF ( Kprint>=2 ) THEN
      WRITE (Lun,*) ' ***** ***** TEST 5 PASSED ***** *****'
      WRITE (Lun,*)
    ENDIF
    !
    ! --- END OF TESTS
    IF ( (ipass1==0).OR.(ipass2==0).OR.(ipass3==0).OR.(ipass4==0).OR.&
        (ipass5==0) ) THEN
      Ipass = 0
      IF ( Kprint>=1 ) WRITE (Lun,99007)
      99007   FORMAT (' ***** DQC36J  FAILED SOME TESTS *****')
    ELSE
      Ipass = 1
      IF ( Kprint>=2 ) WRITE (Lun,99008)
      99008   FORMAT (' ***** DQC36J  PASSED ALL TESTS  *****')
    ENDIF
    99009 FORMAT ('              L2 = ',F5.1,'   L3 = ',F5.1)
    99010 FORMAT (' M1 = ',F5.1,'   M2 = ',F5.1,'   M3 = ',F5.1)
    !
  END SUBROUTINE DQC36J
END MODULE TEST16_MOD
!DECK TEST16
PROGRAM TEST16
  USE TEST16_MOD
  IMPLICIT NONE
  INTEGER I1MACH
  !***BEGIN PROLOGUE  TEST16
  !***PURPOSE  Driver for testing SLATEC subprograms
  !            DRC3JJ   DRC3JM   DRC6J
  !***LIBRARY   SLATEC
  !***CATEGORY  C19
  !***TYPE      DOUBLE PRECISION (TEST15-S, TEST16-D)
  !***KEYWORDS  3J COEFFICIENTS, 3J SYMBOLS, 6J COEFFICIENTS, 6J SYMBOLS,
  !             CLEBSCH-GORDAN COEFFICIENTS, QUICK CHECK DRIVER,
  !             RACAH COEFFICIENTS, VECTOR ADDITION COEFFICIENTS,
  !             WIGNER COEFFICIENTS
  !***AUTHOR  SLATEC Common Mathematical Library Committee
  !***DESCRIPTION
  !
  ! *Usage:
  !     One input data record is required
  !         READ (LIN, '(I1)') KPRINT
  !
  ! *Arguments:
  !     KPRINT = 0  Quick checks - No printing.
  !                 Driver       - Short pass or fail message printed.
  !              1  Quick checks - No message printed for passed tests,
  !                                short message printed for failed tests.
  !                 Driver       - Short pass or fail message printed.
  !              2  Quick checks - Print short message for passed tests,
  !                                fuller information for failed tests.
  !                 Driver       - Pass or fail message printed.
  !              3  Quick checks - Print complete quick check results.
  !                 Driver       - Pass or fail message printed.
  !
  ! *Description:
  !     Driver for testing SLATEC subprograms
  !        DRC3JJ   DRC3JM   DRC6J
  !
  !***REFERENCES  Kirby W. Fong, Thomas H. Jefferson, Tokihiko Suyehiro
  !                 and Lee Walton, Guide to the SLATEC Common Mathema-
  !                 tical Library, April 10, 1990.
  !***ROUTINES CALLED  DQC36J, I1MACH, XERMAX, XSETF, XSETUN
  !***REVISION HISTORY  (YYMMDD)
  !   891130  DATE WRITTEN
  !***END PROLOGUE  TEST16
  INTEGER ipass, kprint, lin, lun, nfail
  !***FIRST EXECUTABLE STATEMENT  TEST16
  lun = I1MACH(2)
  lin = I1MACH(1)
  nfail = 0
  !
  !     Read KPRINT parameter
  !
  READ (lin,'(I1)') kprint
  CALL XERMAX(1000)
  CALL XSETUN(lun)
  IF ( kprint<=1 ) THEN
    CALL XSETF(0)
  ELSE
    CALL XSETF(1)
  ENDIF
  !
  !     Test double precision 3J6J routines
  !
  CALL DQC36J(lun,kprint,ipass)
  IF ( ipass==0 ) nfail = nfail + 1
  !
  !     Write PASS or FAIL message
  !
  IF ( nfail==0 ) THEN
    WRITE (lun,99001)
    99001   FORMAT (/' --------------TEST16 PASSED ALL TESTS----------------')
  ELSE
    WRITE (lun,99002) nfail
    99002   FORMAT (/' ************* WARNING -- ',I5,&
      ' TEST(S) FAILED IN PROGRAM TEST16  *************')
  ENDIF
  STOP
END PROGRAM TEST16
