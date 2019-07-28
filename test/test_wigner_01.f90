MODULE TEST15_MOD
  USE service, ONLY : SP, DP
  IMPLICIT NONE

CONTAINS
  !** QC36J
  SUBROUTINE QC36J(Lun,Kprint,Ipass)
    !> THIS IS A QUICK CHECK PROGRAM FOR THE SUBROUTINES RC3JJ,
    !            RC3JM, AND RC6J, WHICH CALCULATE THE WIGNER COEFFICIENTS,
    !            3J AND 6J.
    !***
    ! **Library:**   SLATEC
    !***
    ! **Category:**  C19
    !***
    ! **Type:**      SINGLE PRECISION (QC36J-S, DQC36J-D)
    !***
    ! **Keywords:**  3J COEFFICIENTS, 3J SYMBOLS, 6J COEFFICIENTS, 6J SYMBOLS,
    !             CLEBSCH-GORDAN COEFFICIENTS, QUICK CHECK,
    !             RACAH COEFFICIENTS, VECTOR ADDITION COEFFICIENTS,
    !             WIGNER COEFFICIENTS
    !***
    ! **Author:**  LOZIER, DANIEL W., (NIST)
    !           MCCLAIN, MARJORIE A., (NIST)
    !           SMITH, JOHN M., (NIST AND GEORGE MASON UNIVERSITY)
    !***
    ! **References:**  MESSIAH, ALBERT., QUANTUM MECHANICS, VOLUME II,
    !               NORTH-HOLLAND PUBLISHING COMPANY, 1963.
    !***
    ! **Routines called:**  NUMXER, R1MACH, RC3JJ, RC3JM, RC6J, XERCLR, XSETF

    !* REVISION HISTORY  (YYMMDD)
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
    USE slatec, ONLY : eps_2_sp, RC3JJ, RC3JM, RC6J
    !
    INTEGER :: Lun, Kprint, Ipass
    !
    CHARACTER string*36, fmt*30, fmt2*13
    INTEGER :: ipass1, ipass2, ipass3, ipass4, ipass5, ier, indexx, &
      i, first, last, nsig, ierjj, ierjm
    INTEGER, PARAMETER :: NDIM=15
    REAL(SP) :: tol, l1, l2, l3, m1, m2, m3, l1min, l1max, m2min, m2max, &
      diff(NDIM), x, jjval, jmval, thrcof(NDIM), sixcof(NDIM)
    !
    REAL(SP), PARAMETER :: r3jj(8) = [ 2.78886675511358515993E-1_SP, &
      -9.53462589245592315447E-2_SP, -6.74199862463242086246E-2_SP, &
      1.53311035167966641297E-1_SP, -1.56446554693685969725E-1_SP, &
      1.09945041215655051079E-1_SP, -5.53623569313171943334E-2_SP, &
      1.79983545113778583298E-2_SP ]
    !
    REAL(SP), PARAMETER :: r3jm(14) = [ 2.09158973288615242614E-2_SP, &
      8.53756555321524722127E-2_SP,  9.08295370868692516943E-2_SP, &
      -3.89054377846499391700E-2_SP, -6.63734970165680635691E-2_SP, &
      6.49524040528389395031E-2_SP, 2.15894310595403759392E-2_SP, &
      -7.78912711785239219992E-2_SP, 3.59764371059543401880E-2_SP, &
      5.47301500021263423079E-2_SP, -7.59678665956761514629E-2_SP, &
      -2.19224445539892113776E-2_SP, 1.01167744280772202424E-1_SP, &
      7.34825726244719704696E-2_SP ]
    !
    REAL(SP), PARAMETER :: r6j(15) = [ 3.49090513837329977746E-2_SP, &
      -3.74302503965979160859E-2_SP, 1.89086639095956018415E-2_SP, &
      7.34244825492864345709E-3_SP, -2.35893518508179445858E-2_SP, &
      1.91347695521543652000E-2_SP, 1.28801739772417220844E-3_SP, &
      -1.93001836629052653977E-2_SP, 1.67730594938288876974E-2_SP, &
      5.50114727485094871674E-3_SP, -2.13543979089683097421E-2_SP, &
      3.46036445143538730828E-3_SP, 2.52095005479558458604E-2_SP, &
      1.48399056122171330285E-2_SP, 2.70857768063318559724E-3_SP ]
    !
    !* FIRST EXECUTABLE STATEMENT  QC36J
    !
    ! --- INITIALIZATION OF TESTS
    tol = 100._SP*eps_2_sp
    IF( Kprint>=2 ) THEN
      WRITE (Lun,*) ' THIS IS QC36J, A TEST PROGRAM FOR THE '//&
        'SINGLE PRECISION 3J6J PACKAGE.'
      WRITE (Lun,*) ' AN EXPLANATION OF THE VARIOUS '//&
        'TESTS CAN BE FOUND IN THE PROGRAM COMMENTS.'
      WRITE (Lun,*)
    END IF
    !
    ! --- FIND NUMBER OF SIGNIFICANT FIGURES FOR FORMATTING
    x = 1._SP/3._SP
    WRITE (string,99001) x
    99001 FORMAT (F35.25)
    DO i = 1, 35
      IF( string(i:i)=='3' ) THEN
        first = i
        EXIT
      END IF
    END DO
    DO i = first, 35
      IF( string(i:i)/='3' ) THEN
        last = i - 1
        GOTO 100
      END IF
    END DO
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
    l2 = 4.5_SP
    l3 = 3.5_SP
    m2 = -3.5_SP
    m3 = 2.5_SP
    CALL RC3JJ(l2,l3,m2,m3,l1min,l1max,thrcof,NDIM,ier)
    IF( ier/=0 ) THEN
      ipass1 = 0
    ELSE
      DO indexx = 1, INT(l1max-l1min) + 1
        m1 = 1._SP
        diff(indexx) = ABS(thrcof(indexx)-r3jj(indexx))
        IF( diff(indexx)>ABS(r3jj(indexx))*tol ) ipass1 = 0
      END DO
    END IF
    IF( Kprint>=3 .OR. (Kprint==2 .AND. ipass1==0) ) THEN
      WRITE (Lun,*) ' TEST 1, RECURRENCE IN L1, COMPARE VALUES OF 3J ', &
        'CALCULATED BY RC3JJ TO'
      WRITE (Lun,*) ' VALUES CALCULATED BY EXPLICIT FORMULA FROM ', &
        'MESSIAH''S QUANTUM MECHANICS'
      WRITE (Lun,99009) l2, l3
      WRITE (Lun,99010) m1, m2, m3
      IF( ier/=0 ) THEN
        WRITE (Lun,*) ' ERROR RETURNED FROM SUBROUTINE ', 'RC3JJ: IER =', ier
      ELSE
        WRITE (Lun,99002)
        99002 FORMAT ('    L1',T31,' RC3JJ VALUE',T67,'FORMULA VALUE')
        DO indexx = 1, INT(l1max-l1min) + 1
          l1 = indexx+ l1min - 1
          WRITE (Lun,fmt) l1, thrcof(indexx), r3jj(indexx)
          IF( diff(indexx)>ABS(r3jj(indexx))*tol ) WRITE (Lun,'(1X,A,F5.1)')&
            'DIFFERENCE EXCEEDS ERROR TOLERANCE FOR L1 =', l1
        END DO
      END IF
    END IF
    IF( ipass1==0 ) THEN
      IF( Kprint>=1 ) THEN
        WRITE (Lun,*) ' ***** ***** TEST 1 FAILED ***** *****'
        WRITE (Lun,*)
      END IF
    ELSEIF( Kprint>=2 ) THEN
      WRITE (Lun,*) ' ***** ***** TEST 1 PASSED ***** *****'
      WRITE (Lun,*)
    END IF
    !
    ! --- TEST 2: COMPARE RC3JM VALUES WITH FORMULA
    ipass2 = 1
    l1 = 8._SP
    l2 = 7.5_SP
    l3 = 6.5_SP
    m1 = 1._SP
    CALL RC3JM(l1,l2,l3,m1,m2min,m2max,thrcof,NDIM,ier)
    IF( ier/=0 ) THEN
      ipass2 = 0
    ELSE
      DO indexx = 1, INT(m2max-m2min) + 1
        m2 = indexx+ m2min - 1
        m3 = -m1 - m2
        diff(indexx) = ABS(thrcof(indexx)-r3jm(indexx))
        IF( diff(indexx)>ABS(r3jm(indexx))*tol ) ipass2 = 0
      END DO
    END IF
    IF( Kprint>=3 .OR. (Kprint==2 .AND. ipass2==0) ) THEN
      WRITE (Lun,*) ' TEST 2, RECURRENCE IN M2, COMPARE VALUES OF 3J ', &
        'CALCULATED BY RC3JM TO'
      WRITE (Lun,*) ' VALUES CALCULATED BY EXPLICIT FORMULA FROM ', &
        'MESSIAH''S QUANTUM MECHANICS'
      WRITE (Lun,99003) l1, l2, l3
      99003 FORMAT (' L1 = ',F5.1,'   L2 = ',F5.1,'   L3 = ',F5.1)
      WRITE (Lun,99004) m1
      99004 FORMAT (' M1 = ',F5.1,'                M3 = -(M1+M2)')
      IF( ier/=0 ) THEN
        WRITE (Lun,*) ' ERROR RETURNED FROM SUBROUTINE ', 'RC3JM: IER =', ier
      ELSE
        WRITE (Lun,99005)
        99005 FORMAT ('    M2',T31,' RC3JM VALUE',T67,'FORMULA VALUE')
        DO indexx = 1, INT(m2max-m2min) + 1
          m2 = indexx+ m2min - 1
          WRITE (Lun,fmt) m2, thrcof(indexx), r3jm(indexx)
          IF( diff(indexx)>ABS(r3jm(indexx))*tol ) WRITE (Lun,'(1X,A,F5.1)')&
            'DIFFERENCE EXCEEDS ERROR TOLERANCE FOR M2 =', m2
        END DO
      END IF
    END IF
    IF( ipass2==0 ) THEN
      IF( Kprint>=1 ) THEN
        WRITE (Lun,*) ' ***** ***** TEST 2 FAILED ***** *****'
        WRITE (Lun,*)
      END IF
    ELSEIF( Kprint>=2 ) THEN
      WRITE (Lun,*) ' ***** ***** TEST 2 PASSED ***** *****'
      WRITE (Lun,*)
    END IF
    !
    ! --- TEST3: COMPARE COMMON VALUE OF RC3JJ AND RC3JM
    ipass3 = 1
    l1 = 100._SP
    l2 = 2._SP
    l3 = 100._SP
    m1 = -10._SP
    m2 = 0._SP
    m3 = 10._SP
    CALL RC3JJ(l2,l3,m2,m3,l1min,l1max,thrcof,NDIM,ierjj)
    jjval = thrcof(3)
    CALL RC3JM(l1,l2,l3,m1,m2min,m2max,thrcof,NDIM,ierjm)
    jmval = thrcof(3)
    IF( ierjj/=0 .OR. ierjm/=0 ) THEN
      ipass3 = 0
    ELSE
      diff(1) = ABS(jjval-jmval)
      IF( diff(1)>0.5_SP*ABS(jjval+jmval)*tol ) ipass3 = 0
    END IF
    IF( Kprint>=3 .OR. (Kprint==2 .AND. ipass3==0) ) THEN
      WRITE (Lun,*) ' TEST 3, COMPARE A COMMON VALUE CALCULATED BY ', &
        'BOTH RC3JJ AND RC3JM'
      WRITE (Lun,*) ' L1 = 100.0   L2 =   2.0   L3 = 100.0'
      WRITE (Lun,*) ' M1 = -10.0   M2 =   0.0   M3 =  10.0'
      IF( ierjj/=0 ) THEN
        WRITE (Lun,*) ' ERROR RETURNED FROM SUBROUTINE ', 'RC3JJ: IER =', ierjj
      ELSEIF( ierjm/=0 ) THEN
        WRITE (Lun,*) ' ERROR RETURNED FROM SUBROUTINE ', 'RC3JM: IER =', ierjm
      ELSE
        WRITE (Lun,fmt2) 'RC3JJ VALUE =', jjval
        WRITE (Lun,fmt2) 'RC3JM VALUE =', jmval
        IF( diff(1)>0.5_SP*ABS(jjval+jmval)*tol ) WRITE (Lun,'(1X,A)')&
          'DIFFERENCE EXCEEDS ERROR TOLERANCE'
      END IF
    END IF
    IF( ipass3==0 ) THEN
      IF( Kprint>=1 ) THEN
        WRITE (Lun,*) ' ***** ***** TEST 3 FAILED ***** *****'
        WRITE (Lun,*)
      END IF
    ELSEIF( Kprint>=2 ) THEN
      WRITE (Lun,*) ' ***** ***** TEST 3 PASSED ***** *****'
      WRITE (Lun,*)
    END IF
    !
    ! --- TEST 4: COMPARE RC6J VALUES WITH FORMULA
    ipass4 = 1
    l2 = 8._SP
    l3 = 7._SP
    m1 = 6.5_SP
    m2 = 7.5_SP
    m3 = 7.5_SP
    CALL RC6J(l2,l3,m1,m2,m3,l1min,l1max,sixcof,NDIM,ier)
    IF( ier/=0 ) THEN
      ipass4 = 0
    ELSE
      DO indexx = 1, INT(l1max-l1min) + 1
        diff(indexx) = ABS(sixcof(indexx)-r6j(indexx))
        IF( diff(indexx)>ABS(r6j(indexx))*tol ) ipass4 = 0
      END DO
    END IF
    IF( Kprint>=3 .OR. (Kprint==2 .AND. ipass4==0) ) THEN
      WRITE (Lun,*) ' TEST 4, RECURRENCE IN L1, COMPARE VALUES OF 6J ', &
        'CALCULATED BY RC6J TO'
      WRITE (Lun,*) ' VALUES CALCULATED BY EXPLICIT FORMULA FROM ', &
        'MESSIAH''S QUANTUM MECHANICS'
      WRITE (Lun,99009) l2, l3
      WRITE (Lun,99010) m1, m2, m3
      IF( ier/=0 ) THEN
        WRITE (Lun,*) ' ERROR RETURNED FROM SUBROUTINE ', 'RC6J: IER =', ier
      ELSE
        WRITE (Lun,99006)
        99006 FORMAT ('    L1',T32,' RC6J VALUE',T67,'FORMULA VALUE')
        DO indexx = 1, INT(l1max-l1min) + 1
          l1 = indexx+ l1min - 1
          WRITE (Lun,fmt) l1, sixcof(indexx), r6j(indexx)
          IF( diff(indexx)>ABS(r6j(indexx))*tol ) WRITE (Lun,'(1X,A,F5.1)')&
            'DIFFERENCE EXCEEDS ERROR TOLERANCE FOR L1 =', l1
        END DO
      END IF
    END IF
    IF( ipass4==0 ) THEN
      IF( Kprint>=1 ) THEN
        WRITE (Lun,*) ' ***** ***** TEST 4 FAILED ***** *****'
        WRITE (Lun,*)
      END IF
    ELSEIF( Kprint>=2 ) THEN
      WRITE (Lun,*) ' ***** ***** TEST 4 PASSED ***** *****'
      WRITE (Lun,*)
    END IF
    !
    ! --- TEST 5: CHECK INVALID INPUT
    ipass5 = 1
!    IF( Kprint<=2 ) THEN
!      control_xer = 0
!    ELSE
!      control_xer = -1
!    END IF
!    IF( Kprint>=3 ) WRITE (Lun,*) ' TEST 5, CHECK FOR PROPER HANDLING ', &
!      'OF INVALID INPUT'
    ! --- RC3JJ: L2-ABS(M2) OR L3-ABS(M3) LESS THAN ZERO (IER=1)
!    l2 = 2._SP
!    l3 = 100._SP
!    m1 = -6._SP
!    m2 = -4._SP
!    m3 = 10._SP
!    IF( Kprint>=3 ) WRITE (Lun,*)
!    num_xer = 0
!    CALL RC3JJ(l2,l3,m2,m3,l1min,l1max,thrcof,NDIM,ier)
!    IF( num_xer/=ier ) ipass5 = 0
    ! --- RC3JJ: L2+ABS(M2) OR L3+ABS(M3) NOT INTEGER (IER=2)
!    l2 = 2._SP
!    l3 = 99.5_SP
!    m1 = -10._SP
!    m2 = 0._SP
!    m3 = 10._SP
!    IF( Kprint>=3 ) WRITE (Lun,*)
!    num_xer = 0
!    CALL RC3JJ(l2,l3,m2,m3,l1min,l1max,thrcof,NDIM,ier)
!    IF( num_xer/=ier ) ipass5 = 0
    ! --- RC3JJ: L1MAX-L1MIN NOT INTEGER (IER=3)
!    l2 = 3.2_SP
!    l3 = 4.5_SP
!    m1 = -1.3_SP
!    m2 = 0.8_SP
!    m3 = 0.5_SP
!    IF( Kprint>=3 ) WRITE (Lun,*)
!    num_xer = 0
!    CALL RC3JJ(l2,l3,m2,m3,l1min,l1max,thrcof,NDIM,ier)
!    IF( num_xer/=ier ) ipass5 = 0
    ! --- RC3JJ: L1MIN GREATER THAN L1MAX (IER=4)
    !            (NO TEST -- THIS ERROR SHOULD NEVER OCCUR)
    ! --- RC3JJ: DIMENSION OF THRCOF TOO SMALL (IER=5)
!    l2 = 10._SP
!    l3 = 150._SP
!    m1 = -10._SP
!    m2 = 0._SP
!    m3 = 10._SP
!    IF( Kprint>=3 ) WRITE (Lun,*)
!    num_xer = 0
!    CALL RC3JJ(l2,l3,m2,m3,l1min,l1max,thrcof,NDIM,ier)
!    IF( num_xer/=ier ) ipass5 = 0
    ! --- RC3JM: L1-ABS(M1) LESS THAN ZERO OR L1+ABS(M1) NOT INTEGER (IER=1)
!    l1 = 100._SP
!    l2 = 2._SP
!    l3 = 100._SP
!    m1 = 150._SP
!    IF( Kprint>=3 ) WRITE (Lun,*)
!    num_xer = 0
!    CALL RC3JM(l1,l2,l3,m1,m2min,m2max,thrcof,NDIM,ier)
!    IF( num_xer/=ier ) ipass5 = 0
    ! --- RC3JM: L1, L2, L3 DO NOT SATISFY TRIANGULAR CONDITION (IER=2)
!    l1 = 20._SP
!    l2 = 5._SP
!    l3 = 10._SP
!    m1 = -10._SP
!    IF( Kprint>=3 ) WRITE (Lun,*)
!    num_xer = 0
!    CALL RC3JM(l1,l2,l3,m1,m2min,m2max,thrcof,NDIM,ier)
!    IF( num_xer/=ier ) ipass5 = 0
    ! --- RC3JM: L1+L2+L3 NOT INTEGER (IER=3)
!    l1 = 1._SP
!    l2 = 1.3_SP
!    l3 = 1.5_SP
!    m1 = 0._SP
!    IF( Kprint>=3 ) WRITE (Lun,*)
!    num_xer = 0
!    CALL RC3JM(l1,l2,l3,m1,m2min,m2max,thrcof,NDIM,ier)
!    IF( num_xer/=ier ) ipass5 = 0
    ! --- RC3JM: M2MAX-M2MIN NOT INTEGER (IER=4)
!    l1 = 1._SP
!    l2 = 1.3_SP
!    l3 = 1.7_SP
!    m1 = 0._SP
!    IF( Kprint>=3 ) WRITE (Lun,*)
!    num_xer = 0
!    CALL RC3JM(l1,l2,l3,m1,m2min,m2max,thrcof,NDIM,ier)
!    IF( num_xer/=ier ) ipass5 = 0
    ! --- RC3JM: M2MIN GREATER THAN M2MAX (IER=5)
    !            (NO TEST -- THIS ERROR SHOULD NEVER OCCUR)
    ! --- RC3JM: DIMENSION OF THRCOF TOO SMALL (IER=6)
!    l1 = 100._SP
!    l2 = 10._SP
!    l3 = 110._SP
!    m1 = -10._SP
!    IF( Kprint>=3 ) WRITE (Lun,*)
!    num_xer = 0
!    CALL RC3JM(l1,l2,l3,m1,m2min,m2max,thrcof,NDIM,ier)
!    IF( num_xer/=ier ) ipass5 = 0
    ! --- RC6J: L2+L3+L5+L6 OR L4+L2+L6 NOT INTEGER (IER=1)
!    l2 = 0.5_SP
!    l3 = 1._SP
!    m1 = 0.5_SP
!    m2 = 2._SP
!    m3 = 3._SP
!    IF( Kprint>=3 ) WRITE (Lun,*)
!    num_xer = 0
!    CALL RC6J(l2,l3,m1,m2,m3,l1min,l1max,sixcof,NDIM,ier)
!    IF( num_xer/=ier ) ipass5 = 0
    ! --- RC6J: L4, L2, L6 TRIANGULAR CONDITION NOT SATISFIED (IER=2)
!    l2 = 1._SP
!    l3 = 3._SP
!    m1 = 5._SP
!    m2 = 6._SP
!    m3 = 2._SP
!    IF( Kprint>=3 ) WRITE (Lun,*)
!    num_xer = 0
!    CALL RC6J(l2,l3,m1,m2,m3,l1min,l1max,sixcof,NDIM,ier)
!    IF( num_xer/=ier ) ipass5 = 0
    ! --- RC6J: L4, L5, L3 TRIANGULAR CONDITION NOT SATISFIED (IER=3)
!    l2 = 4._SP
!    l3 = 1._SP
!    m1 = 5._SP
!    m2 = 3._SP
!    m3 = 2._SP
!    IF( Kprint>=3 ) WRITE (Lun,*)
!    num_xer = 0
!    CALL RC6J(l2,l3,m1,m2,m3,l1min,l1max,sixcof,NDIM,ier)
!    IF( num_xer/=ier ) ipass5 = 0
    ! --- RC6J: L1MAX-L1MIN NOT INTEGER (IER=4)
!    l2 = 0.9_SP
!    l3 = 0.5_SP
!    m1 = 0.9_SP
!    m2 = 0.4_SP
!    m3 = 0.2_SP
!    IF( Kprint>=3 ) WRITE (Lun,*)
!    num_xer = 0
!    CALL RC6J(l2,l3,m1,m2,m3,l1min,l1max,sixcof,NDIM,ier)
!    IF( num_xer/=ier ) ipass5 = 0
    ! --- RC6J: L1MIN GREATER THAN L1MAX (IER=5)
    !           (NO TEST -- THIS ERROR SHOULD NEVER OCCUR)
    ! --- RC6J: DIMENSION OF SIXCOF TOO SMALL (IER=6)
!    l2 = 50._SP
!    l3 = 25._SP
!    m1 = 15._SP
!    m2 = 30._SP
!    m3 = 40._SP
!    IF( Kprint>=3 ) WRITE (Lun,*)
!    num_xer = 0
!    CALL RC6J(l2,l3,m1,m2,m3,l1min,l1max,sixcof,NDIM,ier)
!    IF( num_xer/=ier ) ipass5 = 0
    IF( ipass5==0 ) THEN
      IF( Kprint>=1 ) THEN
        WRITE (Lun,*) ' ***** ***** TEST 5 FAILED ***** *****'
        WRITE (Lun,*)
      END IF
    ELSEIF( Kprint>=2 ) THEN
      WRITE (Lun,*) ' ***** ***** TEST 5 PASSED ***** *****'
      WRITE (Lun,*)
    END IF
    !
    ! --- END OF TESTS
    IF( (ipass1==0) .OR. (ipass2==0) .OR. (ipass3==0) .OR. (ipass4==0) .OR. &
        (ipass5==0) ) THEN
      Ipass = 0
      IF( Kprint>=1 ) WRITE (Lun,99007)
      99007 FORMAT (' *****  QC36J  FAILED SOME TESTS *****')
    ELSE
      Ipass = 1
      IF( Kprint>=2 ) WRITE (Lun,99008)
      99008 FORMAT (' *****  QC36J  PASSED ALL TESTS  *****')
    END IF
    99009 FORMAT ('              L2 = ',F5.1,'   L3 = ',F5.1)
    99010 FORMAT (' M1 = ',F5.1,'   M2 = ',F5.1,'   M3 = ',F5.1)
    !
  END SUBROUTINE QC36J
END MODULE TEST15_MOD
!** TEST15
PROGRAM TEST15
  USE TEST15_MOD, ONLY : QC36J
  USE ISO_FORTRAN_ENV, ONLY : INPUT_UNIT, OUTPUT_UNIT
  USE common_mod, ONLY : GET_ARGUMENT
  IMPLICIT NONE
  !> Driver for testing SLATEC subprograms
  !            RC3JJ    RC3JM    RC6J
  !***
  ! **Library:**   SLATEC
  !***
  ! **Category:**  C19
  !***
  ! **Type:**      SINGLE PRECISION (TEST15-S, TEST16-D)
  !***
  ! **Keywords:**  3J COEFFICIENTS, 3J SYMBOLS, 6J COEFFICIENTS, 6J SYMBOLS,
  !             CLEBSCH-GORDAN COEFFICIENTS, QUICK CHECK DRIVER,
  !             RACAH COEFFICIENTS, VECTOR ADDITION COEFFICIENTS,
  !             WIGNER COEFFICIENTS
  !***
  ! **Author:**  SLATEC Common Mathematical Library Committee
  !***
  ! **Description:**
  !
  !- Usage:
  !     One input data record is required
  !         READ (LIN, '(I1)') KPRINT
  !
  !- Arguments:
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
  !- Description:
  !     Driver for testing SLATEC subprograms
  !        RC3JJ    RC3JM    RC6J
  !
  !***
  ! **References:**  Kirby W. Fong, Thomas H. Jefferson, Tokihiko Suyehiro
  !                 and Lee Walton, Guide to the SLATEC Common Mathema-
  !                 tical Library, April 10, 1990.
  !***
  ! **Routines called:**  I1MACH, QC36J, XERMAX, XSETF

  !* REVISION HISTORY  (YYMMDD)
  !   891130  DATE WRITTEN
  INTEGER :: ipass, kprint, lin, lun, nfail
  !* FIRST EXECUTABLE STATEMENT  TEST15
  lun = OUTPUT_UNIT
  lin = INPUT_UNIT
  nfail = 0
  !
  !     Read KPRINT parameter
  !
  CALL GET_ARGUMENT(kprint)
  !
  !     Test single precision 3J6J routines
  !
  CALL QC36J(lun,kprint,ipass)
  IF( ipass==0 ) nfail = nfail + 1
  !
  !     Write PASS or FAIL message
  !
  IF( nfail==0 ) THEN
    WRITE (lun,99001)
    99001 FORMAT (/' --------------TEST15 PASSED ALL TESTS----------------')
  ELSE
    WRITE (lun,99002) nfail
    99002 FORMAT (/' ************* WARNING -- ',I5,&
      ' TEST(S) FAILED IN PROGRAM TEST15  *************')
  END IF
  STOP
END PROGRAM TEST15
