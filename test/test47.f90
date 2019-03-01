MODULE TEST47_MOD
  IMPLICIT NONE

CONTAINS
  !DECK CDQCK
  SUBROUTINE CDQCK(Lun,Kprint,Ipass)
    IMPLICIT NONE
    !***BEGIN PROLOGUE  CDQCK
    !***PURPOSE  Quick check for SLATEC routines CDRIV1, CDRIV2 and CDRIV3.
    !***LIBRARY   SLATEC (SDRIVE)
    !***CATEGORY  I1A2, I1A1B
    !***TYPE      COMPLEX (SDQCK-S, DDQCK-D, CDQCK-C)
    !***KEYWORDS  CDRIV1, CDRIV2, CDRIV3, QUICK CHECK, SDRIVE
    !***AUTHOR  Kahaner, D. K., (NIST)
    !             National Institute of Standards and Technology
    !             Gaithersburg, MD  20899
    !           Sutherland, C. D., (LANL)
    !             Mail Stop D466
    !             Los Alamos National Laboratory
    !             Los Alamos, NM  87545
    !***DESCRIPTION
    !
    !  For assistance in determining the cause of a failure of these
    !  routines contact C. D. Sutherland at commercial telephone number
    !  (505)667-6949, FTS telephone number 8-843-6949, or electronic mail
    !  address CDS@LANL.GOV .
    !
    !***ROUTINES CALLED  CDF, CDRIV1, CDRIV2, CDRIV3, R1MACH, XERCLR
    !***REVISION HISTORY  (YYMMDD)
    !   890405  DATE WRITTEN
    !   890405  Revised to meet SLATEC standards.
    !***END PROLOGUE  CDQCK
    REAL eps, ewt(1), HMAX, R1MACH, t, tout
    INTEGER ierflg, IERROR, IMPL, Ipass, Kprint, leniw, leniwx, lenw, &
      LENWMX, lenwx, LIWMX, Lun, mint, MITER, ML, mstate, MU, &
      MXORD, MXSTEP, N, nde, nfe, nje, NROOT, nstate, nstep, &
      NTASK, nx
    PARAMETER (HMAX=15.E0,IERROR=3,IMPL=0,LENWMX=342,LIWMX=53,MITER=5,ML=2,&
      MU=2,MXORD=5,MXSTEP=1000,N=3,NROOT=0,NTASK=1)
    COMPLEX alfa, work(LENWMX), y(N+1)
    INTEGER iwork(LIWMX)
    DATA ewt(1)/.00001E0/
    !***FIRST EXECUTABLE STATEMENT  CDQCK
    alfa = (1.E0,1.E0)
    eps = R1MACH(4)**(1.E0/3.E0)
    Ipass = 1
    !                                            Exercise CDRIV1 for problem
    !                                            with known solution.
    y(4) = alfa
    t = 0.E0
    y(1) = 10.E0
    y(2) = 0.E0
    y(3) = 10.E0
    tout = 10.E0
    mstate = 1
    lenw = 342
    CALL CDRIV1(N,t,y,CDF,tout,mstate,eps,work,lenw,ierflg)
    nstep = INT( work(lenw-(N+50)+3) )
    nfe = INT( work(lenw-(N+50)+4) )
    nje = INT( work(lenw-(N+50)+5) )
    IF ( mstate/=2 ) THEN
      IF ( Kprint==1 ) THEN
        WRITE (Lun,&
          '('' While using CDRIV1, a solution was not obtained.'' //)'&
          )
      ELSEIF ( Kprint>=2 ) THEN
        WRITE (Lun,&
          '('' While using CDRIV1, a solution was not obtained.'')'&
          )
        WRITE (Lun,&
          '('' The values of parameters, results, and statistical quantities are:'')')
        WRITE (Lun,*) ' MSTATE = ', mstate, ', Error number = ', ierflg
        WRITE (Lun,*) ' N ', N, ', EPS ', eps, ', LENW ', lenw
        WRITE (Lun,*) ' T ', t
        WRITE (Lun,*) ' Y(1) ', y(1)
        WRITE (Lun,*) ' Y(2) ', y(2)
        WRITE (Lun,*) ' Y(3) ', y(3)
        WRITE (Lun,*) ' Number of steps taken is  ', nstep
        WRITE (Lun,*) ' Number of evaluations of the right hand side is  ', &
          nfe
        WRITE (Lun,*) ' Number of evaluations of the Jacobian matrix is  ', &
          nje
        WRITE (Lun,'(/)')
      ENDIF
      Ipass = 0
    ELSEIF ( ABS(0.620174E0-ABS(y(1)))>eps**(2.E0/3.E0).OR.&
        ABS(0.392232E0-ABS(y(2)))>eps**(2.E0/3.E0).OR.ABS(1.E0-ABS(y(3)))&
        >eps**(2.E0/3.E0) ) THEN
      IF ( Kprint==1 ) THEN
        WRITE (Lun,&
          '('' CDRIV1:The solution determined is not accurate enough.'&
          ' //)')
      ELSEIF ( Kprint==2 ) THEN
        WRITE (Lun,&
          '('' CDRIV1:The solution determined is not accurate enough.'&
          ')')
        WRITE (Lun,&
          '('' The values of parameters, results, and statistical quantities are:'')')
        WRITE (Lun,*) ' EPS = ', eps
        WRITE (Lun,*) ' T ', t
        WRITE (Lun,*) ' Y(1) ', y(1)
        WRITE (Lun,*) ' Y(2) ', y(2)
        WRITE (Lun,*) ' Y(3) ', y(3)
        WRITE (Lun,*) ' Number of steps taken is  ', nstep
        WRITE (Lun,*) ' Number of evaluations of the right hand side is  ', &
          nfe
        WRITE (Lun,*) ' Number of evaluations of the Jacobian matrix is  ', &
          nje
        WRITE (Lun,'(/)')
      ENDIF
      Ipass = 0
    ELSEIF ( Kprint==2 ) THEN
      WRITE (Lun,&
        '('' CDRIV1:The solution determined met the expected values.'' //)')
    ELSEIF ( Kprint==3 ) THEN
      WRITE (Lun,&
        '('' CDRIV1:The solution determined met the expected values.'')')
      WRITE (Lun,'('' The values of results are '')')
      WRITE (Lun,*) ' T ', t
      WRITE (Lun,*) ' Y(1) ', y(1)
      WRITE (Lun,*) ' Y(2) ', y(2)
      WRITE (Lun,*) ' Y(3) ', y(3)
      WRITE (Lun,'(/)')
    ENDIF
    CALL XERCLR
    !                                         Run CDRIV1 with invalid input.
    nx = 201
    t = 0.E0
    y(1) = 10.E0
    y(2) = 0.E0
    y(3) = 10.E0
    y(4) = alfa
    tout = 10.E0
    mstate = 1
    lenw = 342
    CALL CDRIV1(nx,t,y,CDF,tout,mstate,eps,work,lenw,ierflg)
    IF ( ierflg/=21 ) THEN
      IF ( Kprint==1 ) THEN
        WRITE (Lun,&
          '('' CDRIV1:An invalid parameter has not been correctly detected.'' //)')
      ELSEIF ( Kprint>=2 ) THEN
        WRITE (Lun,&
          '('' CDRIV1:An invalid parameter has not been correctly detected.'')')
        WRITE (Lun,*) ' The value of N was set to ', nx
        WRITE (Lun,*) ' MSTATE = ', mstate, ', Error number = ', ierflg
        WRITE (Lun,&
          '('' The values of parameters, results, and statistical quantities are:'')')
        WRITE (Lun,*) ' EPS ', eps, ', LENW ', lenw
        WRITE (Lun,*) ' T ', t
        WRITE (Lun,*) ' Y(1) ', y(1)
        WRITE (Lun,*) ' Y(2) ', y(2)
        WRITE (Lun,*) ' Y(3) ', y(3)
        WRITE (Lun,*) ' Number of steps taken is  ', nstep
        WRITE (Lun,*) ' Number of evaluations of the right hand side is  ', &
          nfe
        WRITE (Lun,*) ' Number of evaluations of the Jacobian matrix is  ', &
          nje
        WRITE (Lun,'(/)')
      ENDIF
      Ipass = 0
    ELSEIF ( Kprint==2 ) THEN
      WRITE (Lun,&
        '('' CDRIV1:An invalid parameter has been correctly detected.'' //)')
    ELSEIF ( Kprint==3 ) THEN
      WRITE (Lun,&
        '('' CDRIV1:An invalid parameter has been correctly detected.'')')
      WRITE (Lun,*) ' The value of N was set to ', nx
      WRITE (Lun,*) ' MSTATE = ', mstate, ', Error number = ', ierflg
      WRITE (Lun,'(/)')
    ENDIF
    CALL XERCLR
    !                                            Exercise CDRIV2 for problem
    !                                            with known solution.
    t = 0.E0
    y(1) = 10.E0
    y(2) = 0.E0
    y(3) = 10.E0
    y(4) = alfa
    mstate = 1
    tout = 10.E0
    mint = 1
    lenw = 298
    leniw = 50
    CALL CDRIV2(N,t,y,CDF,tout,mstate,NROOT,eps,ewt,mint,work,lenw,iwork,&
      leniw,CDF,ierflg)
    nstep = iwork(3)
    nfe = iwork(4)
    nje = iwork(5)
    IF ( mstate/=2 ) THEN
      IF ( Kprint==1 ) THEN
        WRITE (Lun,&
          '('' While using CDRIV2, a solution was not obtained.'' //)'&
          )
      ELSEIF ( Kprint>=2 ) THEN
        WRITE (Lun,&
          '('' While using CDRIV2, a solution was not obtained.'')'&
          )
        WRITE (Lun,*) ' MSTATE = ', mstate, ', Error number = ', ierflg
        WRITE (Lun,&
          '('' The values of parameters, results, and statistical quantities are:'')')
        WRITE (Lun,*) ' EPS = ', eps, ', EWT ', ewt
        WRITE (Lun,*) ' MINT = ', mint, ', LENW ', lenw, ', LENIW ', &
          leniw
        WRITE (Lun,*) ' T ', t
        WRITE (Lun,*) ' Y(1) ', y(1)
        WRITE (Lun,*) ' Y(2) ', y(2)
        WRITE (Lun,*) ' Y(3) ', y(3)
        WRITE (Lun,*) ' Number of steps taken is  ', nstep
        WRITE (Lun,*) ' Number of evaluations of the right hand side is  ', &
          nfe
        WRITE (Lun,*) ' Number of evaluations of the Jacobian matrix is  ', &
          nje
        WRITE (Lun,'(/)')
      ENDIF
      Ipass = 0
    ELSEIF ( ABS(0.620174E0-ABS(y(1)))>eps**(2.E0/3.E0).OR.&
        ABS(0.392232E0-ABS(y(2)))>eps**(2.E0/3.E0).OR.ABS(1.E0-ABS(y(3)))&
        >eps**(2.E0/3.E0) ) THEN
      IF ( Kprint==1 ) THEN
        WRITE (Lun,&
          '('' CDRIV2:The solution determined is not accurate enough. //'')')
      ELSEIF ( Kprint==2 ) THEN
        WRITE (Lun,&
          '('' CDRIV2:The solution determined is not accurate enough.'')')
        WRITE (Lun,&
          '('' The values of parameters, results, and statistical quantities are:'')')
        WRITE (Lun,*) ' EPS = ', eps, ', EWT = ', ewt
        WRITE (Lun,*) ' T ', t
        WRITE (Lun,*) ' Y(1) ', y(1)
        WRITE (Lun,*) ' Y(2) ', y(2)
        WRITE (Lun,*) ' Y(3) ', y(3)
        WRITE (Lun,*) ' Number of steps taken is  ', nstep
        WRITE (Lun,*) ' Number of evaluations of the right hand side is  ', &
          nfe
        WRITE (Lun,*) ' Number of evaluations of the Jacobian matrix is  ', &
          nje
        WRITE (Lun,'(/)')
      ENDIF
      Ipass = 0
    ELSEIF ( Kprint==2 ) THEN
      WRITE (Lun,&
        '('' CDRIV2:The solution determined met the expected values.'' //)')
    ELSEIF ( Kprint==3 ) THEN
      WRITE (Lun,&
        '('' CDRIV2:The solution determined met the expected values.'')')
      WRITE (Lun,'('' The values of results are '')')
      WRITE (Lun,*) ' T ', t
      WRITE (Lun,*) ' Y(1) ', y(1)
      WRITE (Lun,*) ' Y(2) ', y(2)
      WRITE (Lun,*) ' Y(3) ', y(3)
      WRITE (Lun,'(/)')
    ENDIF
    CALL XERCLR
    !                                         Run CDRIV2 with invalid input.
    t = 0.E0
    y(1) = 10.E0
    y(2) = 0.E0
    y(3) = 10.E0
    y(4) = alfa
    tout = 10.E0
    mstate = 1
    mint = 1
    lenwx = 1
    leniw = 50
    CALL CDRIV2(N,t,y,CDF,tout,mstate,NROOT,eps,ewt,mint,work,lenwx,iwork,&
      leniw,CDF,ierflg)
    IF ( ierflg/=32 ) THEN
      IF ( Kprint==1 ) THEN
        WRITE (Lun,&
          '('' CDRIV2:An invalid parameter has not been correctly detected.'' //)')
      ELSEIF ( Kprint>=2 ) THEN
        WRITE (Lun,&
          '('' CDRIV2:An invalid parameter has not been correctly detected.'')')
        WRITE (Lun,*) ' The value of LENW was set to ', lenwx
        WRITE (Lun,*) ' MSTATE = ', mstate, ', Error number = ', ierflg
        WRITE (Lun,&
          '('' The values of parameters, results, and statistical quantities are:'')')
        WRITE (Lun,*) ' EPS ', eps, ', MINT ', mint, ', LENW ', lenw, &
          ', LENIW ', leniw
        WRITE (Lun,*) ' T ', t
        WRITE (Lun,*) ' Y(1) ', y(1)
        WRITE (Lun,*) ' Y(2) ', y(2)
        WRITE (Lun,*) ' Y(3) ', y(3)
        WRITE (Lun,*) ' Number of steps taken is  ', nstep
        WRITE (Lun,*) ' Number of evaluations of the right hand side is  ', &
          nfe
        WRITE (Lun,*) ' Number of evaluations of the Jacobian matrix is  ', &
          nje
        WRITE (Lun,'(/)')
      ENDIF
      Ipass = 0
    ELSEIF ( Kprint==2 ) THEN
      WRITE (Lun,&
        '('' CDRIV2:An invalid parameter has been correctly detected.'' //)')
    ELSEIF ( Kprint==3 ) THEN
      WRITE (Lun,&
        '('' CDRIV2:An invalid parameter has been correctly detected.'')')
      WRITE (Lun,*) ' The value of LENW was set to ', lenwx
      WRITE (Lun,*) ' MSTATE = ', mstate, ', Error number = ', ierflg
      WRITE (Lun,'(/)')
    ENDIF
    CALL XERCLR
    !                                            Exercise CDRIV3 for problem
    !                                            with known solution.
    t = 0.E0
    y(1) = 10.E0
    y(2) = 0.E0
    y(3) = 10.E0
    y(4) = alfa
    nstate = 1
    tout = 10.E0
    mint = 2
    lenw = 301
    leniw = 53
    CALL CDRIV3(N,t,y,CDF,nstate,tout,NTASK,NROOT,eps,ewt,IERROR,mint,MITER,&
      IMPL,ML,MU,MXORD,HMAX,work,lenw,iwork,leniw,CDF,CDF,nde,&
      MXSTEP,CDF,CDF,ierflg)
    nstep = iwork(3)
    nfe = iwork(4)
    nje = iwork(5)
    IF ( nstate/=2 ) THEN
      IF ( Kprint==1 ) THEN
        WRITE (Lun,&
          '('' While using CDRIV3, a solution was not obtained.'' //)'&
          )
      ELSEIF ( Kprint>=2 ) THEN
        WRITE (Lun,&
          '('' While using CDRIV3, a solution was not obtained.'')'&
          )
        WRITE (Lun,*) ' MSTATE = ', mstate, ', Error number = ', ierflg
        WRITE (Lun,&
          '('' The values of parameters, results, and statistical quantities are:'')')
        WRITE (Lun,*) ' EPS = ', eps, ', EWT = ', ewt, ', IERROR = ', &
          IERROR
        WRITE (Lun,*) ' MINT = ', mint, ', MITER = ', MITER, ', IMPL = ', &
          IMPL
        WRITE (Lun,*) ' T ', t
        WRITE (Lun,*) ' Y(1) ', y(1)
        WRITE (Lun,*) ' Y(2) ', y(2)
        WRITE (Lun,*) ' Y(3) ', y(3)
        WRITE (Lun,*) ' Number of steps taken is  ', nstep
        WRITE (Lun,*) ' Number of evaluations of the right hand side is  ', &
          nfe
        WRITE (Lun,*) ' Number of evaluations of the Jacobian matrix is  ', &
          nje
        WRITE (Lun,'(/)')
      ENDIF
      Ipass = 0
    ELSEIF ( ABS(0.620174E0-ABS(y(1)))>eps**(2.E0/3.E0).OR.&
        ABS(0.392232E0-ABS(y(2)))>eps**(2.E0/3.E0).OR.ABS(1.E0-ABS(y(3)))&
        >eps**(2.E0/3.E0) ) THEN
      IF ( Kprint==1 ) THEN
        WRITE (Lun,&
          '('' CDRIV3:The solution determined is not accurate enough.'&
          ' //)')
      ELSEIF ( Kprint>=2 ) THEN
        WRITE (Lun,&
          '('' CDRIV3:The solution determined is not accurate enough.'&
          ')')
        WRITE (Lun,&
          '('' The values of parameters, results, and statistical quantities are:'')')
        WRITE (Lun,*) ' EPS = ', eps, ', EWT = ', ewt, ', IERROR = ', &
          IERROR
        WRITE (Lun,*) ' MINT = ', mint, ', MITER = ', MITER, ', IMPL = ', &
          IMPL
        WRITE (Lun,*) ' T ', t
        WRITE (Lun,*) ' Y(1) ', y(1)
        WRITE (Lun,*) ' Y(2) ', y(2)
        WRITE (Lun,*) ' Y(3) ', y(3)
        WRITE (Lun,*) ' Number of steps taken is  ', nstep
        WRITE (Lun,*) ' Number of evaluations of the right hand side is  ', &
          nfe
        WRITE (Lun,*) ' Number of evaluations of the Jacobian matrix is  ', &
          nje
        WRITE (Lun,'(/)')
      ENDIF
      Ipass = 0
    ELSEIF ( Kprint==2 ) THEN
      WRITE (Lun,&
        '('' CDRIV3:The solution determined met the expected values.'' //)')
    ELSEIF ( Kprint==3 ) THEN
      WRITE (Lun,&
        '('' CDRIV3:The solution determined met the expected values.'')')
      WRITE (Lun,'('' The values of results are '')')
      WRITE (Lun,*) ' T ', t
      WRITE (Lun,*) ' Y(1) ', y(1)
      WRITE (Lun,*) ' Y(2) ', y(2)
      WRITE (Lun,*) ' Y(3) ', y(3)
      WRITE (Lun,'(/)')
    ENDIF
    CALL XERCLR
    !                                         Run CDRIV3 with invalid input.
    t = 0.E0
    y(1) = 10.E0
    y(2) = 0.E0
    y(3) = 10.E0
    y(4) = alfa
    nstate = 1
    tout = 10.E0
    mint = 2
    lenw = 301
    leniwx = 1
    CALL CDRIV3(N,t,y,CDF,nstate,tout,NTASK,NROOT,eps,ewt,IERROR,mint,MITER,&
      IMPL,ML,MU,MXORD,HMAX,work,lenw,iwork,leniwx,CDF,CDF,nde,&
      MXSTEP,CDF,CDF,ierflg)
    IF ( ierflg/=33 ) THEN
      IF ( Kprint==1 ) THEN
        WRITE (Lun,&
          '('' CDRIV3:An invalid parameter has not been correctly detected.'' //)')
      ELSEIF ( Kprint>=2 ) THEN
        WRITE (Lun,&
          '('' CDRIV3:An invalid parameter has not been correctly detected.'')')
        WRITE (Lun,*) ' The value of LENIW was set to ', leniwx
        WRITE (Lun,*) ' NSTATE = ', nstate, ', Error number = ', ierflg
        WRITE (Lun,&
          '('' The values of parameters, results, and statistical quantities are:'')')
        WRITE (Lun,*) ' EPS = ', eps, ', EWT = ', ewt, ', IERROR = ', &
          IERROR
        WRITE (Lun,*) ' MINT = ', mint, ', MITER = ', MITER, ', IMPL = ', &
          IMPL
        WRITE (Lun,*) ' T ', t
        WRITE (Lun,*) ' Y(1) ', y(1)
        WRITE (Lun,*) ' Y(2) ', y(2)
        WRITE (Lun,*) ' Y(3) ', y(3)
        WRITE (Lun,*) ' Number of steps taken is  ', nstep
        WRITE (Lun,*) ' Number of evaluations of the right hand side is  ', &
          nfe
        WRITE (Lun,*) ' Number of evaluations of the Jacobian matrix is  ', &
          nje
        WRITE (Lun,'(/)')
      ENDIF
      Ipass = 0
    ELSEIF ( Kprint==2 ) THEN
      WRITE (Lun,&
        '('' CDRIV3:An invalid parameter has been correctly detected.'' //)')
    ELSEIF ( Kprint==3 ) THEN
      WRITE (Lun,&
        '('' CDRIV3:An invalid parameter has been correctly detected.'')')
      WRITE (Lun,*) ' The value of LENIW was set to ', leniwx
      WRITE (Lun,*) ' NSTATE = ', nstate, ', Error number = ', ierflg
      WRITE (Lun,'(/)')
    ENDIF
    CALL XERCLR
  END SUBROUTINE CDQCK
  !DECK CDF
  SUBROUTINE CDF(N,T,Y,Yp)
    IMPLICIT NONE
    !***BEGIN PROLOGUE  CDF
    !***SUBSIDIARY
    !***PURPOSE  Quick check for SLATEC routines CDRIV1, CDRIV2 and CDRIV3.
    !***LIBRARY   SLATEC (SDRIVE)
    !***CATEGORY  I1A2, I1A1B
    !***TYPE      COMPLEX (SDF-S, DDF-D, CDF-C)
    !***KEYWORDS  CDRIV1, CDRIV2, CDRIV3, QUICK CHECK, SDRIVE
    !***AUTHOR  Kahaner, D. K., (NIST)
    !             National Institute of Standards and Technology
    !             Gaithersburg, MD  20899
    !           Sutherland, C. D., (LANL)
    !             Mail Stop D466
    !             Los Alamos National Laboratory
    !             Los Alamos, NM  87545
    !***SEE ALSO  CDQCK
    !***ROUTINES CALLED  (NONE)
    !***REVISION HISTORY  (YYMMDD)
    !   890405  DATE WRITTEN
    !   890405  Revised to meet SLATEC standards.
    !***END PROLOGUE  CDF
    REAL T
    COMPLEX alfa, Y(*), Yp(*)
    INTEGER N
    !***FIRST EXECUTABLE STATEMENT  CDF
    alfa = Y(N+1)
    Yp(1) = 1.E0 + alfa*(Y(2)-Y(1)) - Y(1)*Y(3)
    Yp(2) = alfa*(Y(1)-Y(2)) - Y(2)*Y(3)
    Yp(3) = 1.E0 - Y(3)*(Y(1)+Y(2))
  END SUBROUTINE CDF
END MODULE TEST47_MOD
!DECK TEST47
PROGRAM TEST47
  USE TEST47_MOD
  IMPLICIT NONE
  INTEGER I1MACH
  !***BEGIN PROLOGUE  TEST47
  !***PURPOSE  Driver for testing SLATEC subprograms
  !            CDRIV1  CDRIV2  CDRIV3
  !***LIBRARY   SLATEC
  !***CATEGORY  I1A2, I1A1B
  !***TYPE      COMPLEX (TEST45-S, TEST46-D, TEST47-C)
  !***KEYWORDS  SDRIVE, QUICK CHECK DRIVER
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
  !        CDRIV1  CDRIV2  CDRIV3
  !
  !***REFERENCES  Kirby W. Fong, Thomas H. Jefferson, Tokihiko Suyehiro
  !                 and Lee Walton, Guide to the SLATEC Common Mathema-
  !                 tical Library, April 10, 1990.
  !***ROUTINES CALLED  CDQCK, I1MACH, XERMAX, XSETF, XSETUN
  !***REVISION HISTORY  (YYMMDD)
  !   920801  DATE WRITTEN
  !***END PROLOGUE  TEST47
  INTEGER ipass, kprint, lin, lun, nfail
  !***FIRST EXECUTABLE STATEMENT  TEST47
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
  !     Test complex SDRIVE
  !
  CALL CDQCK(lun,kprint,ipass)
  IF ( ipass==0 ) nfail = nfail + 1
  !
  !     Write PASS or FAIL message
  !
  IF ( nfail==0 ) THEN
    WRITE (lun,99001)
    99001 FORMAT (/' --------------TEST47 PASSED ALL TESTS----------------')
  ELSE
    WRITE (lun,99002) nfail
    99002 FORMAT (/' ************* WARNING -- ',I5,&
      ' TEST(S) FAILED IN PROGRAM TEST47 *************')
  ENDIF
  STOP
END PROGRAM TEST47
