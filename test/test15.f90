MODULE TEST15_MOD
  IMPLICIT NONE

CONTAINS
  !DECK QC36J
  SUBROUTINE QC36J(Lun,Kprint,Ipass)
    IMPLICIT NONE
    !***BEGIN PROLOGUE  QC36J
    !***SUBSIDIARY
    !***PURPOSE  THIS IS A QUICK CHECK PROGRAM FOR THE SUBROUTINES RC3JJ,
    !            RC3JM, AND RC6J, WHICH CALCULATE THE WIGNER COEFFICIENTS,
    !            3J AND 6J.
    !***LIBRARY   SLATEC
    !***CATEGORY  C19
    !***TYPE      SINGLE PRECISION (QC36J-S, DQC36J-D)
    !***KEYWORDS  3J COEFFICIENTS, 3J SYMBOLS, 6J COEFFICIENTS, 6J SYMBOLS,
    !             CLEBSCH-GORDAN COEFFICIENTS, QUICK CHECK,
    !             RACAH COEFFICIENTS, VECTOR ADDITION COEFFICIENTS,
    !             WIGNER COEFFICIENTS
    !***AUTHOR  LOZIER, DANIEL W., (NIST)
    !           MCCLAIN, MARJORIE A., (NIST)
    !           SMITH, JOHN M., (NIST AND GEORGE MASON UNIVERSITY)
    !***REFERENCES  MESSIAH, ALBERT., QUANTUM MECHANICS, VOLUME II,
    !               NORTH-HOLLAND PUBLISHING COMPANY, 1963.
    !***ROUTINES CALLED  NUMXER, R1MACH, RC3JJ, RC3JM, RC6J, XERCLR, XSETF
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
    !           removing all calls to subroutine RACAH.  These changes were
    !           made by M. McClain.
    !***END PROLOGUE  QC36J
    !
    INTEGER Lun, Kprint, Ipass
    !
    CHARACTER string*36, fmt*30, fmt2*13
    INTEGER ipass1, ipass2, ipass3, ipass4, ipass5, NDIM, ier, index, &
      i, first, last, nsig, NUMXER, nerr, ierjj, ierjm
    PARAMETER (NDIM=15)
    REAL tol, l1, l2, l3, m1, m2, m3, l1min, l1max, m2min, m2max, &
      diff(NDIM), R1MACH, x, jjval, jmval, thrcof(NDIM), sixcof(NDIM)&
      , r3jj(8), r3jm(14), r6j(15)
    !
    DATA r3jj/2.78886675511358515993E-1, -9.53462589245592315447E-2, &
      -6.74199862463242086246E-2, 1.53311035167966641297E-1, &
      -1.56446554693685969725E-1, 1.09945041215655051079E-1, &
      -5.53623569313171943334E-2, 1.79983545113778583298E-2/
    !
    DATA r3jm/2.09158973288615242614E-2, 8.53756555321524722127E-2, &
      9.08295370868692516943E-2, -3.89054377846499391700E-2, &
      -6.63734970165680635691E-2, 6.49524040528389395031E-2, &
      2.15894310595403759392E-2, -7.78912711785239219992E-2, &
      3.59764371059543401880E-2, 5.47301500021263423079E-2, &
      -7.59678665956761514629E-2, -2.19224445539892113776E-2, &
      1.01167744280772202424E-1, 7.34825726244719704696E-2/
    !
    DATA r6j/3.49090513837329977746E-2, -3.74302503965979160859E-2, &
      1.89086639095956018415E-2, 7.34244825492864345709E-3, &
      -2.35893518508179445858E-2, 1.91347695521543652000E-2, &
      1.28801739772417220844E-3, -1.93001836629052653977E-2, &
      1.67730594938288876974E-2, 5.50114727485094871674E-3, &
      -2.13543979089683097421E-2, 3.46036445143538730828E-3, &
      2.52095005479558458604E-2, 1.48399056122171330285E-2, &
      2.70857768063318559724E-3/
    !
    !***FIRST EXECUTABLE STATEMENT  QC36J
    !
    ! --- INITIALIZATION OF TESTS
    tol = 100.0*R1MACH(3)
    IF ( Kprint>=2 ) THEN
      WRITE (Lun,*) ' THIS IS QC36J, A TEST PROGRAM FOR THE '//&
        'SINGLE PRECISION 3J6J PACKAGE.'
      WRITE (Lun,*) ' AN EXPLANATION OF THE VARIOUS '//&
        'TESTS CAN BE FOUND IN THE PROGRAM COMMENTS.'
      WRITE (Lun,*)
    ENDIF
    !
    ! --- FIND NUMBER OF SIGNIFICANT FIGURES FOR FORMATTING
    x = 1.0/3.0
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
    ! --- TEST 1: COMPARE RC3JJ VALUES WITH FORMULA
    ipass1 = 1
    l2 = 4.5
    l3 = 3.5
    m2 = -3.5
    m3 = 2.5
    CALL RC3JJ(l2,l3,m2,m3,l1min,l1max,thrcof,NDIM,ier)
    IF ( ier/=0 ) THEN
      ipass1 = 0
    ELSE
      DO index = 1, INT(l1max-l1min) + 1
        m1 = 1.0
        diff(index) = ABS(thrcof(index)-r3jj(index))
        IF ( diff(index)>ABS(r3jj(index))*tol ) ipass1 = 0
      ENDDO
    ENDIF
    IF ( Kprint>=3.OR.(Kprint==2.AND.ipass1==0) ) THEN
      WRITE (Lun,*) ' TEST 1, RECURRENCE IN L1, COMPARE VALUES OF 3J ', &
        'CALCULATED BY RC3JJ TO'
      WRITE (Lun,*) ' VALUES CALCULATED BY EXPLICIT FORMULA FROM ', &
        'MESSIAH''S QUANTUM MECHANICS'
      WRITE (Lun,99009) l2, l3
      WRITE (Lun,99010) m1, m2, m3
      IF ( ier/=0 ) THEN
        WRITE (Lun,*) ' ERROR RETURNED FROM SUBROUTINE ', 'RC3JJ: IER =', &
          ier
      ELSE
        WRITE (Lun,99002)
        99002 FORMAT ('    L1',T31,' RC3JJ VALUE',T67,'FORMULA VALUE')
        DO index = 1, INT(l1max-l1min) + 1
          l1 = index+ l1min - 1
          WRITE (Lun,fmt) l1, thrcof(index), r3jj(index)
          IF ( diff(index)>ABS(r3jj(index))*tol ) WRITE (Lun,'(1X,A,F5.1)')&
            'DIFFERENCE EXCEEDS ERROR TOLERANCE FOR L1 =', l1
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
    ! --- TEST 2: COMPARE RC3JM VALUES WITH FORMULA
    ipass2 = 1
    l1 = 8.0
    l2 = 7.5
    l3 = 6.5
    m1 = 1.0
    CALL RC3JM(l1,l2,l3,m1,m2min,m2max,thrcof,NDIM,ier)
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
        'CALCULATED BY RC3JM TO'
      WRITE (Lun,*) ' VALUES CALCULATED BY EXPLICIT FORMULA FROM ', &
        'MESSIAH''S QUANTUM MECHANICS'
      WRITE (Lun,99003) l1, l2, l3
      99003 FORMAT (' L1 = ',F5.1,'   L2 = ',F5.1,'   L3 = ',F5.1)
      WRITE (Lun,99004) m1
      99004 FORMAT (' M1 = ',F5.1,'                M3 = -(M1+M2)')
      IF ( ier/=0 ) THEN
        WRITE (Lun,*) ' ERROR RETURNED FROM SUBROUTINE ', 'RC3JM: IER =', &
          ier
      ELSE
        WRITE (Lun,99005)
        99005 FORMAT ('    M2',T31,' RC3JM VALUE',T67,'FORMULA VALUE')
        DO index = 1, INT(m2max-m2min) + 1
          m2 = index+ m2min - 1
          WRITE (Lun,fmt) m2, thrcof(index), r3jm(index)
          IF ( diff(index)>ABS(r3jm(index))*tol ) WRITE (Lun,'(1X,A,F5.1)')&
            'DIFFERENCE EXCEEDS ERROR TOLERANCE FOR M2 =', m2
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
    ! --- TEST3: COMPARE COMMON VALUE OF RC3JJ AND RC3JM
    ipass3 = 1
    l1 = 100.0
    l2 = 2.0
    l3 = 100.0
    m1 = -10.0
    m2 = 0.0
    m3 = 10.0
    CALL RC3JJ(l2,l3,m2,m3,l1min,l1max,thrcof,NDIM,ierjj)
    jjval = thrcof(3)
    CALL RC3JM(l1,l2,l3,m1,m2min,m2max,thrcof,NDIM,ierjm)
    jmval = thrcof(3)
    IF ( ierjj/=0.OR.ierjm/=0 ) THEN
      ipass3 = 0
    ELSE
      diff(1) = ABS(jjval-jmval)
      IF ( diff(1)>0.5*ABS(jjval+jmval)*tol ) ipass3 = 0
    ENDIF
    IF ( Kprint>=3.OR.(Kprint==2.AND.ipass3==0) ) THEN
      WRITE (Lun,*) ' TEST 3, COMPARE A COMMON VALUE CALCULATED BY ', &
        'BOTH RC3JJ AND RC3JM'
      WRITE (Lun,*) ' L1 = 100.0   L2 =   2.0   L3 = 100.0'
      WRITE (Lun,*) ' M1 = -10.0   M2 =   0.0   M3 =  10.0'
      IF ( ierjj/=0 ) THEN
        WRITE (Lun,*) ' ERROR RETURNED FROM SUBROUTINE ', 'RC3JJ: IER =', &
          ierjj
      ELSEIF ( ierjm/=0 ) THEN
        WRITE (Lun,*) ' ERROR RETURNED FROM SUBROUTINE ', 'RC3JM: IER =', &
          ierjm
      ELSE
        WRITE (Lun,fmt2) 'RC3JJ VALUE =', jjval
        WRITE (Lun,fmt2) 'RC3JM VALUE =', jmval
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
    ! --- TEST 4: COMPARE RC6J VALUES WITH FORMULA
    ipass4 = 1
    l2 = 8.0
    l3 = 7.0
    m1 = 6.5
    m2 = 7.5
    m3 = 7.5
    CALL RC6J(l2,l3,m1,m2,m3,l1min,l1max,sixcof,NDIM,ier)
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
        'CALCULATED BY RC6J TO'
      WRITE (Lun,*) ' VALUES CALCULATED BY EXPLICIT FORMULA FROM ', &
        'MESSIAH''S QUANTUM MECHANICS'
      WRITE (Lun,99009) l2, l3
      WRITE (Lun,99010) m1, m2, m3
      IF ( ier/=0 ) THEN
        WRITE (Lun,*) ' ERROR RETURNED FROM SUBROUTINE ', 'RC6J: IER =', ier
      ELSE
        WRITE (Lun,99006)
        99006 FORMAT ('    L1',T32,' RC6J VALUE',T67,'FORMULA VALUE')
        DO index = 1, INT(l1max-l1min) + 1
          l1 = index+ l1min - 1
          WRITE (Lun,fmt) l1, sixcof(index), r6j(index)
          IF ( diff(index)>ABS(r6j(index))*tol ) WRITE (Lun,'(1X,A,F5.1)')&
            'DIFFERENCE EXCEEDS ERROR TOLERANCE FOR L1 =', l1
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
    ! --- RC3JJ: L2-ABS(M2) OR L3-ABS(M3) LESS THAN ZERO (IER=1)
    l2 = 2.0
    l3 = 100.0
    m1 = -6.0
    m2 = -4.0
    m3 = 10.0
    IF ( Kprint>=3 ) WRITE (Lun,*)
    CALL XERCLR
    CALL RC3JJ(l2,l3,m2,m3,l1min,l1max,thrcof,NDIM,ier)
    IF ( NUMXER(nerr)/=ier ) ipass5 = 0
    ! --- RC3JJ: L2+ABS(M2) OR L3+ABS(M3) NOT INTEGER (IER=2)
    l2 = 2.0
    l3 = 99.5
    m1 = -10.0
    m2 = 0.0
    m3 = 10.0
    IF ( Kprint>=3 ) WRITE (Lun,*)
    CALL XERCLR
    CALL RC3JJ(l2,l3,m2,m3,l1min,l1max,thrcof,NDIM,ier)
    IF ( NUMXER(nerr)/=ier ) ipass5 = 0
    ! --- RC3JJ: L1MAX-L1MIN NOT INTEGER (IER=3)
    l2 = 3.2
    l3 = 4.5
    m1 = -1.3
    m2 = 0.8
    m3 = 0.5
    IF ( Kprint>=3 ) WRITE (Lun,*)
    CALL XERCLR
    CALL RC3JJ(l2,l3,m2,m3,l1min,l1max,thrcof,NDIM,ier)
    IF ( NUMXER(nerr)/=ier ) ipass5 = 0
    ! --- RC3JJ: L1MIN GREATER THAN L1MAX (IER=4)
    !            (NO TEST -- THIS ERROR SHOULD NEVER OCCUR)
    ! --- RC3JJ: DIMENSION OF THRCOF TOO SMALL (IER=5)
    l2 = 10.0
    l3 = 150.0
    m1 = -10.0
    m2 = 0.0
    m3 = 10.0
    IF ( Kprint>=3 ) WRITE (Lun,*)
    CALL XERCLR
    CALL RC3JJ(l2,l3,m2,m3,l1min,l1max,thrcof,NDIM,ier)
    IF ( NUMXER(nerr)/=ier ) ipass5 = 0
    ! --- RC3JM: L1-ABS(M1) LESS THAN ZERO OR L1+ABS(M1) NOT INTEGER (IER=1)
    l1 = 100.0
    l2 = 2.0
    l3 = 100.0
    m1 = 150.0
    IF ( Kprint>=3 ) WRITE (Lun,*)
    CALL XERCLR
    CALL RC3JM(l1,l2,l3,m1,m2min,m2max,thrcof,NDIM,ier)
    IF ( NUMXER(nerr)/=ier ) ipass5 = 0
    ! --- RC3JM: L1, L2, L3 DO NOT SATISFY TRIANGULAR CONDITION (IER=2)
    l1 = 20.0
    l2 = 5.0
    l3 = 10.0
    m1 = -10.0
    IF ( Kprint>=3 ) WRITE (Lun,*)
    CALL XERCLR
    CALL RC3JM(l1,l2,l3,m1,m2min,m2max,thrcof,NDIM,ier)
    IF ( NUMXER(nerr)/=ier ) ipass5 = 0
    ! --- RC3JM: L1+L2+L3 NOT INTEGER (IER=3)
    l1 = 1.0
    l2 = 1.3
    l3 = 1.5
    m1 = 0.0
    IF ( Kprint>=3 ) WRITE (Lun,*)
    CALL XERCLR
    CALL RC3JM(l1,l2,l3,m1,m2min,m2max,thrcof,NDIM,ier)
    IF ( NUMXER(nerr)/=ier ) ipass5 = 0
    ! --- RC3JM: M2MAX-M2MIN NOT INTEGER (IER=4)
    l1 = 1.0
    l2 = 1.3
    l3 = 1.7
    m1 = 0.0
    IF ( Kprint>=3 ) WRITE (Lun,*)
    CALL XERCLR
    CALL RC3JM(l1,l2,l3,m1,m2min,m2max,thrcof,NDIM,ier)
    IF ( NUMXER(nerr)/=ier ) ipass5 = 0
    ! --- RC3JM: M2MIN GREATER THAN M2MAX (IER=5)
    !            (NO TEST -- THIS ERROR SHOULD NEVER OCCUR)
    ! --- RC3JM: DIMENSION OF THRCOF TOO SMALL (IER=6)
    l1 = 100.0
    l2 = 10.0
    l3 = 110.0
    m1 = -10.0
    IF ( Kprint>=3 ) WRITE (Lun,*)
    CALL XERCLR
    CALL RC3JM(l1,l2,l3,m1,m2min,m2max,thrcof,NDIM,ier)
    IF ( NUMXER(nerr)/=ier ) ipass5 = 0
    ! --- RC6J: L2+L3+L5+L6 OR L4+L2+L6 NOT INTEGER (IER=1)
    l2 = 0.5
    l3 = 1.0
    m1 = 0.5
    m2 = 2.0
    m3 = 3.0
    IF ( Kprint>=3 ) WRITE (Lun,*)
    CALL XERCLR
    CALL RC6J(l2,l3,m1,m2,m3,l1min,l1max,sixcof,NDIM,ier)
    IF ( NUMXER(nerr)/=ier ) ipass5 = 0
    ! --- RC6J: L4, L2, L6 TRIANGULAR CONDITION NOT SATISFIED (IER=2)
    l2 = 1.0
    l3 = 3.0
    m1 = 5.0
    m2 = 6.0
    m3 = 2.0
    IF ( Kprint>=3 ) WRITE (Lun,*)
    CALL XERCLR
    CALL RC6J(l2,l3,m1,m2,m3,l1min,l1max,sixcof,NDIM,ier)
    IF ( NUMXER(nerr)/=ier ) ipass5 = 0
    ! --- RC6J: L4, L5, L3 TRIANGULAR CONDITION NOT SATISFIED (IER=3)
    l2 = 4.0
    l3 = 1.0
    m1 = 5.0
    m2 = 3.0
    m3 = 2.0
    IF ( Kprint>=3 ) WRITE (Lun,*)
    CALL XERCLR
    CALL RC6J(l2,l3,m1,m2,m3,l1min,l1max,sixcof,NDIM,ier)
    IF ( NUMXER(nerr)/=ier ) ipass5 = 0
    ! --- RC6J: L1MAX-L1MIN NOT INTEGER (IER=4)
    l2 = 0.9
    l3 = 0.5
    m1 = 0.9
    m2 = 0.4
    m3 = 0.2
    IF ( Kprint>=3 ) WRITE (Lun,*)
    CALL XERCLR
    CALL RC6J(l2,l3,m1,m2,m3,l1min,l1max,sixcof,NDIM,ier)
    IF ( NUMXER(nerr)/=ier ) ipass5 = 0
    ! --- RC6J: L1MIN GREATER THAN L1MAX (IER=5)
    !           (NO TEST -- THIS ERROR SHOULD NEVER OCCUR)
    ! --- RC6J: DIMENSION OF SIXCOF TOO SMALL (IER=6)
    l2 = 50.0
    l3 = 25.0
    m1 = 15.0
    m2 = 30.0
    m3 = 40.0
    IF ( Kprint>=3 ) WRITE (Lun,*)
    CALL XERCLR
    CALL RC6J(l2,l3,m1,m2,m3,l1min,l1max,sixcof,NDIM,ier)
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
      99007 FORMAT (' *****  QC36J  FAILED SOME TESTS *****')
    ELSE
      Ipass = 1
      IF ( Kprint>=2 ) WRITE (Lun,99008)
      99008 FORMAT (' *****  QC36J  PASSED ALL TESTS  *****')
    ENDIF
    99009 FORMAT ('              L2 = ',F5.1,'   L3 = ',F5.1)
    99010 FORMAT (' M1 = ',F5.1,'   M2 = ',F5.1,'   M3 = ',F5.1)
    !
  END SUBROUTINE QC36J
END MODULE TEST15_MOD
!DECK TEST15
PROGRAM TEST15
  USE TEST15_MOD
  IMPLICIT NONE
  INTEGER I1MACH
  !***BEGIN PROLOGUE  TEST15
  !***PURPOSE  Driver for testing SLATEC subprograms
  !            RC3JJ    RC3JM    RC6J
  !***LIBRARY   SLATEC
  !***CATEGORY  C19
  !***TYPE      SINGLE PRECISION (TEST15-S, TEST16-D)
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
  !        RC3JJ    RC3JM    RC6J
  !
  !***REFERENCES  Kirby W. Fong, Thomas H. Jefferson, Tokihiko Suyehiro
  !                 and Lee Walton, Guide to the SLATEC Common Mathema-
  !                 tical Library, April 10, 1990.
  !***ROUTINES CALLED  I1MACH, QC36J, XERMAX, XSETF, XSETUN
  !***REVISION HISTORY  (YYMMDD)
  !   891130  DATE WRITTEN
  !***END PROLOGUE  TEST15
  INTEGER ipass, kprint, lin, lun, nfail
  !***FIRST EXECUTABLE STATEMENT  TEST15
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
  !     Test single precision 3J6J routines
  !
  CALL QC36J(lun,kprint,ipass)
  IF ( ipass==0 ) nfail = nfail + 1
  !
  !     Write PASS or FAIL message
  !
  IF ( nfail==0 ) THEN
    WRITE (lun,99001)
    99001 FORMAT (/' --------------TEST15 PASSED ALL TESTS----------------')
  ELSE
    WRITE (lun,99002) nfail
    99002 FORMAT (/' ************* WARNING -- ',I5,&
      ' TEST(S) FAILED IN PROGRAM TEST15  *************')
  ENDIF
  STOP
END PROGRAM TEST15