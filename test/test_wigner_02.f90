MODULE TEST16_MOD
  USE service, ONLY : SP, DP
  IMPLICIT NONE

CONTAINS
  !** DQC36J
  SUBROUTINE DQC36J(Lun,Kprint,Ipass)
    !> THIS IS A QUICK CHECK PROGRAM FOR THE SUBROUTINES DRC3JJ,
    !            DRC3JM, AND DRC6J, WHICH CALCULATE THE WIGNER COEFFICIENTS,
    !            3J AND 6J.
    !***
    ! **Library:**   SLATEC
    !***
    ! **Category:**  C19
    !***
    ! **Type:**      DOUBLE PRECISION (QC36J-S, DQC36J-D)
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
    ! **Routines called:**  D1MACH, DRC3JJ, DRC3JM, DRC6J, NUMXER, XERCLR,
    !                     XSETF

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
    !           removing all calls to subroutine DRACAH.  These changes were
    !           made by M. McClain.
    USE slatec, ONLY : eps_2_dp, DRC3JJ, DRC3JM, DRC6J, control_xer
    !
    INTEGER :: Lun, Kprint, Ipass
    !
    CHARACTER string*36, fmt*30, fmt2*13
    INTEGER :: ipass1, ipass2, ipass3, ipass4, ipass5, ier, indexx, &
      i, first, last, nsig, ierjj, ierjm
    INTEGER, PARAMETER :: NDIM = 15
    REAL(DP) :: tol, l1, l2, l3, m1, m2, m3, l1min, l1max, m2min, m2max, &
      diff(NDIM), x, jjval, jmval, thrcof(NDIM), sixcof(NDIM)
    !
    REAL(DP), PARAMETER :: r3jj(8) = [ 2.7888667551135851599272400859506249646427E-1_DP, &
      -9.5346258924559231544677592152721599861388E-2_DP, &
      -6.7419986246324208624649067643642846008908E-2_DP, &
      1.5331103516796664129653641995122585402552E-1_DP, &
      -1.5644655469368596972508355755184909201031E-1_DP, &
      1.0994504121565505107947893271429777505797E-1_DP, &
      -5.5362356931317194333395729256559987745156E-2_DP, &
      1.7998354511377858329814092962590761537262E-2_DP ]
    !
    REAL(DP), PARAMETER :: r3jm(14) = [ 2.0915897328861524261384476677886072045904E-2_DP, &
      8.5375655532152472212727551895879778672762E-2_DP, &
      9.0829537086869251694343772675676175068677E-2_DP, &
      -3.8905437784649939169989036459327065796765E-2_DP, &
      -6.6373497016568063569146501153397525003444E-2_DP, &
      6.4952404052838939503061387831391216401903E-2_DP, &
      2.1589431059540375939250708046202926313913E-2_DP, &
      -7.7891271178523921999229618972588887261359E-2_DP, &
      3.5976437105954340188005810512211794384411E-2_DP, &
      5.4730150002126342307937096038252488407360E-2_DP, &
      -7.5967866595676151462927617736745078548338E-2_DP, &
      -2.1922444553989211377558215380002910257762E-2_DP, &
      1.0116774428077220242411199686231560525497E-1_DP, &
      7.3482572624471970469595137204530687381176E-2_DP ]
    !
    REAL(DP), PARAMETER :: r6j(15) = [ 3.4909051383732997774596981092927782159095E-2_DP, &
      -3.7430250396597916085929064401358002747549E-2_DP, &
      1.8908663909595601841537964135129184202064E-2_DP, &
      7.3424482549286434570947151839589351100581E-3_DP, &
      -2.3589351850817944585847816357296508528608E-2_DP, &
      1.9134769552154365200026782557432864615918E-2_DP, &
      1.2880173977241722084434864685591278730958E-3_DP, &
      -1.9300183662905265397749119277519305417805E-2_DP, &
      1.6773059493828887697413611251392749162229E-2_DP, &
      5.5011472748509487167380502058890639729979E-3_DP, &
      -2.1354397908968309742136976853078409839580E-2_DP, &
      3.4603644514353873082775312319159137043869E-3_DP, &
      2.5209500547955845860442730268272167527589E-2_DP, &
      1.4839905612217133028540464232557124565509E-2_DP, &
      2.7085776806331855972407001825016114677027E-3_DP ]
    !
    !* FIRST EXECUTABLE STATEMENT  DQC36J
    !
    ! --- INITIALIZATION OF TESTS
    tol = 100._DP*eps_2_dp
    IF( Kprint>=2 ) THEN
      WRITE (Lun,*) ' THIS IS DQC36J, A TEST PROGRAM FOR THE '//&
        'DOUBLE PRECISION 3J6J PACKAGE.'
      WRITE (Lun,*) ' AN EXPLANATION OF THE VARIOUS '//&
        'TESTS CAN BE FOUND IN THE PROGRAM COMMENTS.'
      WRITE (Lun,*)
    END IF
    !
    ! --- FIND NUMBER OF SIGNIFICANT FIGURES FOR FORMATTING
    x = 1._DP/3._DP
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
    ! --- TEST 1: COMPARE DRC3JJ VALUES WITH FORMULA
    ipass1 = 1
    l2 = 4.5_DP
    l3 = 3.5_DP
    m2 = -3.5_DP
    m3 = 2.5_DP
    CALL DRC3JJ(l2,l3,m2,m3,l1min,l1max,thrcof,NDIM,ier)
    IF( ier/=0 ) THEN
      ipass1 = 0
    ELSE
      DO indexx = 1, INT(l1max-l1min) + 1
        m1 = 1._DP
        diff(indexx) = ABS(thrcof(indexx)-r3jj(indexx))
        IF( diff(indexx)>ABS(r3jj(indexx))*tol ) ipass1 = 0
      END DO
    END IF
    IF( Kprint>=3 .OR. (Kprint==2 .AND. ipass1==0) ) THEN
      WRITE (Lun,*) ' TEST 1, RECURRENCE IN L1, COMPARE VALUES OF 3J ', &
        'CALCULATED BY DRC3JJ TO'
      WRITE (Lun,*) ' VALUES CALCULATED BY EXPLICIT FORMULA FROM ', &
        'MESSIAH''S QUANTUM MECHANICS'
      WRITE (Lun,99009) l2, l3
      WRITE (Lun,99010) m1, m2, m3
      IF( ier/=0 ) THEN
        WRITE (Lun,*) ' ERROR RETURNED FROM SUBROUTINE ', 'DRC3JJ: IER =', ier
      ELSE
        WRITE (Lun,99002)
        99002 FORMAT ('    L1',T31,'DRC3JJ VALUE',T67,'FORMULA VALUE')
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
    ! --- TEST 2: COMPARE DRC3JM VALUES WITH FORMULA
    ipass2 = 1
    l1 = 8._DP
    l2 = 7.5_DP
    l3 = 6.5_DP
    m1 = 1._DP
    CALL DRC3JM(l1,l2,l3,m1,m2min,m2max,thrcof,NDIM,ier)
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
        'CALCULATED BY DRC3JM TO'
      WRITE (Lun,*) ' VALUES CALCULATED BY EXPLICIT FORMULA FROM ', &
        'MESSIAH''S QUANTUM MECHANICS'
      WRITE (Lun,99003) l1, l2, l3
      99003 FORMAT (' L1 = ',F5.1,'   L2 = ',F5.1,'   L3 = ',F5.1)
      WRITE (Lun,99004) m1
      99004 FORMAT (' M1 = ',F5.1,'                M3 = -(M1+M2)')
      IF( ier/=0 ) THEN
        WRITE (Lun,*) ' ERROR RETURNED FROM SUBROUTINE ', 'DRC3JM: IER =', ier
      ELSE
        WRITE (Lun,99005)
        99005 FORMAT ('    M2',T31,'DRC3JM VALUE',T67,'FORMULA VALUE')
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
    ! --- TEST3: COMPARE COMMON VALUE OF DRC3JJ AND DRC3JM
    ipass3 = 1
    l1 = 100._DP
    l2 = 2._DP
    l3 = 100._DP
    m1 = -10._DP
    m2 = 0._DP
    m3 = 10._DP
    CALL DRC3JJ(l2,l3,m2,m3,l1min,l1max,thrcof,NDIM,ierjj)
    jjval = thrcof(3)
    CALL DRC3JM(l1,l2,l3,m1,m2min,m2max,thrcof,NDIM,ierjm)
    jmval = thrcof(3)
    IF( ierjj/=0 .OR. ierjm/=0 ) THEN
      ipass3 = 0
    ELSE
      diff(1) = ABS(jjval-jmval)
      IF( diff(1)>0.5_SP*ABS(jjval+jmval)*tol ) ipass3 = 0
    END IF
    IF( Kprint>=3 .OR. (Kprint==2 .AND. ipass3==0) ) THEN
      WRITE (Lun,*) ' TEST 3, COMPARE A COMMON VALUE CALCULATED BY ', &
        'BOTH DRC3JJ AND DRC3JM'
      WRITE (Lun,*) ' L1 = 100.0   L2 =   2.0   L3 = 100.0'
      WRITE (Lun,*) ' M1 = -10.0   M2 =   0.0   M3 =  10.0'
      IF( ierjj/=0 ) THEN
        WRITE (Lun,*) ' ERROR RETURNED FROM SUBROUTINE ', 'DRC3JJ: IER =', ierjj
      ELSEIF( ierjm/=0 ) THEN
        WRITE (Lun,*) ' ERROR RETURNED FROM SUBROUTINE ', 'DRC3JM: IER =', ierjm
      ELSE
        WRITE (Lun,fmt2) 'DRC3JJ VALUE =', jjval
        WRITE (Lun,fmt2) 'DRC3JM VALUE =', jmval
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
    ! --- TEST 4: COMPARE DRC6J VALUES WITH FORMULA
    ipass4 = 1
    l2 = 8._DP
    l3 = 7._DP
    m1 = 6.5_DP
    m2 = 7.5_DP
    m3 = 7.5_DP
    CALL DRC6J(l2,l3,m1,m2,m3,l1min,l1max,sixcof,NDIM,ier)
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
        'CALCULATED BY DRC6J TO'
      WRITE (Lun,*) ' VALUES CALCULATED BY EXPLICIT FORMULA FROM ', &
        'MESSIAH''S QUANTUM MECHANICS'
      WRITE (Lun,99009) l2, l3
      WRITE (Lun,99010) m1, m2, m3
      IF( ier/=0 ) THEN
        WRITE (Lun,*) ' ERROR RETURNED FROM SUBROUTINE ', 'DRC6J: IER =', ier
      ELSE
        WRITE (Lun,99006)
        99006 FORMAT ('    L1',T32,'DRC6J VALUE',T67,'FORMULA VALUE')
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
    IF( Kprint<=2 ) THEN
      control_xer = 0
    ELSE
      control_xer = -1
    END IF
    IF( Kprint>=3 ) WRITE (Lun,*) ' TEST 5, CHECK FOR PROPER HANDLING ', &
      'OF INVALID INPUT'
    ! --- DRC3JJ: L2-ABS(M2) OR L3-ABS(M3) LESS THAN ZERO (IER=1)
!    l2 = 2._DP
!    l3 = 100._DP
!    m1 = -6._DP
!    m2 = -4._DP
!    m3 = 10._DP
!    IF( Kprint>=3 ) WRITE (Lun,*)
!    num_xer = 0
!    CALL DRC3JJ(l2,l3,m2,m3,l1min,l1max,thrcof,NDIM,ier)
!    IF( num_xer/=ier ) ipass5 = 0
    ! --- DRC3JJ: L2+ABS(M2) OR L3+ABS(M3) NOT INTEGER (IER=2)
!    l2 = 2._DP
!    l3 = 99.5_DP
!    m1 = -10._DP
!    m2 = 0._DP
!    m3 = 10._DP
!    IF( Kprint>=3 ) WRITE (Lun,*)
!    num_xer = 0
!    CALL DRC3JJ(l2,l3,m2,m3,l1min,l1max,thrcof,NDIM,ier)
!    IF( num_xer/=ier ) ipass5 = 0
    ! --- DRC3JJ: L1MAX-L1MIN NOT INTEGER (IER=3)
!    l2 = 3.2_DP
!    l3 = 4.5_DP
!    m1 = -1.3_DP
!    m2 = 0.8_DP
!    m3 = 0.5_DP
!    IF( Kprint>=3 ) WRITE (Lun,*)
!    num_xer = 0
!    CALL DRC3JJ(l2,l3,m2,m3,l1min,l1max,thrcof,NDIM,ier)
!    IF( num_xer/=ier ) ipass5 = 0
    ! --- DRC3JJ: L1MIN GREATER THAN L1MAX (IER=4)
    !             (NO TEST -- THIS ERROR SHOULD NEVER OCCUR)
    ! --- DRC3JJ: DIMENSION OF THRCOF TOO SMALL (IER=5)
!    l2 = 10._DP
!    l3 = 150._DP
!    m1 = -10._DP
!    m2 = 0._DP
!    m3 = 10._DP
!    IF( Kprint>=3 ) WRITE (Lun,*)
!    num_xer = 0
!    CALL DRC3JJ(l2,l3,m2,m3,l1min,l1max,thrcof,NDIM,ier)
!    IF( num_xer/=ier ) ipass5 = 0
    ! --- DRC3JM: L1-ABS(M1) < ZERO OR L1+ABS(M1) NOT INTEGER (IER=1)
!    l1 = 100._DP
!    l2 = 2._DP
!    l3 = 100._DP
!    m1 = 150._DP
!    IF( Kprint>=3 ) WRITE (Lun,*)
!    num_xer = 0
!    CALL DRC3JM(l1,l2,l3,m1,m2min,m2max,thrcof,NDIM,ier)
!    IF( num_xer/=ier ) ipass5 = 0
    ! --- DRC3JM: L1, L2, L3 DO NOT SATISFY TRIANGULAR CONDITION (IER=2)
!    l1 = 20._DP
!    l2 = 5._DP
!    l3 = 10._DP
!    m1 = -10._DP
!    IF( Kprint>=3 ) WRITE (Lun,*)
!    num_xer = 0
!    CALL DRC3JM(l1,l2,l3,m1,m2min,m2max,thrcof,NDIM,ier)
!    IF( num_xer/=ier ) ipass5 = 0
    ! --- DRC3JM: L1+L2+L3 NOT INTEGER (IER=3)
!    l1 = 1._DP
!    l2 = 1.3_DP
!    l3 = 1.5_DP
!    m1 = 0._DP
!    IF( Kprint>=3 ) WRITE (Lun,*)
!    num_xer = 0
!    CALL DRC3JM(l1,l2,l3,m1,m2min,m2max,thrcof,NDIM,ier)
!    IF( num_xer/=ier ) ipass5 = 0
    ! --- DRC3JM: M2MAX-M2MIN NOT INTEGER (IER=4)
!    l1 = 1._DP
!    l2 = 1.3_DP
!    l3 = 1.7_DP
!    m1 = 0._DP
!    IF( Kprint>=3 ) WRITE (Lun,*)
!    num_xer = 0
!    CALL DRC3JM(l1,l2,l3,m1,m2min,m2max,thrcof,NDIM,ier)
!    IF( num_xer/=ier ) ipass5 = 0
    ! --- DRC3JM: M2MIN GREATER THAN M2MAX (IER=5)
    !             (NO TEST -- THIS ERROR SHOULD NEVER OCCUR)
    ! --- DRC3JM: DIMENSION OF THRCOF TOO SMALL (IER=6)
!    l1 = 100._DP
!    l2 = 10._DP
!    l3 = 110._DP
!    m1 = -10._DP
!    IF( Kprint>=3 ) WRITE (Lun,*)
!    num_xer = 0
!    CALL DRC3JM(l1,l2,l3,m1,m2min,m2max,thrcof,NDIM,ier)
!    IF( num_xer/=ier ) ipass5 = 0
    ! --- DRC6J: L2+L3+L5+L6 OR L4+L2+L6 NOT INTEGER (IER=1)
!    l2 = 0.5_DP
!    l3 = 1._DP
!    m1 = 0.5_DP
!    m2 = 2._DP
!    m3 = 3._DP
!    IF( Kprint>=3 ) WRITE (Lun,*)
!    num_xer = 0
!    CALL DRC6J(l2,l3,m1,m2,m3,l1min,l1max,sixcof,NDIM,ier)
!    IF( num_xer/=ier ) ipass5 = 0
    ! --- DRC6J: L4, L2, L6 TRIANGULAR CONDITION NOT SATISFIED (IER=2)
!    l2 = 1._DP
!    l3 = 3._DP
!    m1 = 5._DP
!    m2 = 6._DP
!    m3 = 2._DP
!    IF( Kprint>=3 ) WRITE (Lun,*)
!    num_xer = 0
!    CALL DRC6J(l2,l3,m1,m2,m3,l1min,l1max,sixcof,NDIM,ier)
!    IF( num_xer/=ier ) ipass5 = 0
    ! --- DRC6J: L4, L5, L3 TRIANGULAR CONDITION NOT SATISFIED (IER=3)
!    l2 = 4._DP
!    l3 = 1._DP
!    m1 = 5._DP
!    m2 = 3._DP
!    m3 = 2._DP
!    IF( Kprint>=3 ) WRITE (Lun,*)
!    num_xer = 0
!    CALL DRC6J(l2,l3,m1,m2,m3,l1min,l1max,sixcof,NDIM,ier)
!    IF( num_xer/=ier ) ipass5 = 0
    ! --- DRC6J: L1MAX-L1MIN NOT INTEGER (IER=4)
!    l2 = 0.9_DP
!    l3 = 0.5_DP
!    m1 = 0.9_DP
!    m2 = 0.4_DP
!    m3 = 0.2_DP
!    IF( Kprint>=3 ) WRITE (Lun,*)
!    num_xer = 0
!    CALL DRC6J(l2,l3,m1,m2,m3,l1min,l1max,sixcof,NDIM,ier)
!    IF( num_xer/=ier ) ipass5 = 0
    ! --- DRC6J: L1MIN GREATER THAN L1MAX (IER=5)
    !            (NO TEST -- THIS ERROR SHOULD NEVER OCCUR)
    ! --- DRC6J: DIMENSION OF SIXCOF TOO SMALL (IER=6)
!    l2 = 50._DP
!    l3 = 25._DP
!    m1 = 15._DP
!    m2 = 30._DP
!    m3 = 40._DP
!    IF( Kprint>=3 ) WRITE (Lun,*)
!    num_xer = 0
!    CALL DRC6J(l2,l3,m1,m2,m3,l1min,l1max,sixcof,NDIM,ier)
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
      99007 FORMAT (' ***** DQC36J  FAILED SOME TESTS *****')
    ELSE
      Ipass = 1
      IF( Kprint>=2 ) WRITE (Lun,99008)
      99008 FORMAT (' ***** DQC36J  PASSED ALL TESTS  *****')
    END IF
    99009 FORMAT ('              L2 = ',F5.1,'   L3 = ',F5.1)
    99010 FORMAT (' M1 = ',F5.1,'   M2 = ',F5.1,'   M3 = ',F5.1)
    !
  END SUBROUTINE DQC36J
END MODULE TEST16_MOD
!** TEST16
PROGRAM TEST16
  USE TEST16_MOD, ONLY : DQC36J
  USE ISO_FORTRAN_ENV, ONLY : INPUT_UNIT, OUTPUT_UNIT
  USE slatec, ONLY : control_xer, max_xer
  USE common_mod, ONLY : GET_ARGUMENT
  IMPLICIT NONE
  !> Driver for testing SLATEC subprograms
  !            DRC3JJ   DRC3JM   DRC6J
  !***
  ! **Library:**   SLATEC
  !***
  ! **Category:**  C19
  !***
  ! **Type:**      DOUBLE PRECISION (TEST15-S, TEST16-D)
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
  !        DRC3JJ   DRC3JM   DRC6J
  !
  !***
  ! **References:**  Kirby W. Fong, Thomas H. Jefferson, Tokihiko Suyehiro
  !                 and Lee Walton, Guide to the SLATEC Common Mathema-
  !                 tical Library, April 10, 1990.
  !***
  ! **Routines called:**  DQC36J, I1MACH, XERMAX, XSETF

  !* REVISION HISTORY  (YYMMDD)
  !   891130  DATE WRITTEN
  INTEGER :: ipass, kprint, lin, lun, nfail
  !* FIRST EXECUTABLE STATEMENT  TEST16
  lun = OUTPUT_UNIT
  lin = INPUT_UNIT
  nfail = 0
  !
  !     Read KPRINT parameter
  !
  CALL GET_ARGUMENT(kprint)
  max_xer = 1000
  IF( kprint<=1 ) THEN
    control_xer = 0
  ELSE
    control_xer = 1
  END IF
  !
  !     Test double precision 3J6J routines
  !
  CALL DQC36J(lun,kprint,ipass)
  IF( ipass==0 ) nfail = nfail + 1
  !
  !     Write PASS or FAIL message
  !
  IF( nfail==0 ) THEN
    WRITE (lun,99001)
    99001 FORMAT (/' --------------TEST16 PASSED ALL TESTS----------------')
  ELSE
    WRITE (lun,99002) nfail
    99002 FORMAT (/' ************* WARNING -- ',I5,&
      ' TEST(S) FAILED IN PROGRAM TEST16  *************')
  END IF
  STOP
END PROGRAM TEST16
