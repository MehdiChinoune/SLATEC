!*==CDQCK.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK CDQCK
      SUBROUTINE CDQCK(Lun,Kprint,Ipass)
      IMPLICIT NONE
!*--CDQCK5
!*** Start of declarations inserted by SPAG
      REAL CDF
!*** End of declarations inserted by SPAG
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
      EXTERNAL CDF
      REAL eps , ewt(1) , HMAX , R1MACH , t , tout
      INTEGER ierflg , IERROR , IMPL , Ipass , Kprint , leniw , leniwx , lenw , 
     &        LENWMX , lenwx , LIWMX , Lun , mint , MITER , ML , mstate , MU , 
     &        MXORD , MXSTEP , N , nde , nfe , nje , NROOT , nstate , nstep , 
     &        NTASK , nx
      PARAMETER (HMAX=15.E0,IERROR=3,IMPL=0,LENWMX=342,LIWMX=53,MITER=5,ML=2,
     &           MU=2,MXORD=5,MXSTEP=1000,N=3,NROOT=0,NTASK=1)
      COMPLEX alfa , work(LENWMX) , y(N+1)
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
      nstep = work(lenw-(N+50)+3)
      nfe = work(lenw-(N+50)+4)
      nje = work(lenw-(N+50)+5)
      IF ( mstate/=2 ) THEN
        IF ( Kprint==1 ) THEN
          WRITE (Lun,
     & '('' While using CDRIV1, a solution was not '',        ''obtained.'' //)'
     & )
        ELSEIF ( Kprint>=2 ) THEN
          WRITE (Lun,
     &    '('' While using CDRIV1, a solution was not '',        ''obtained.'')'
     &    )
          WRITE (Lun,
     &'('' The values of parameters, results, and '',        ''statistical quant
     &ities are:'')')
          WRITE (Lun,*) ' MSTATE = ' , mstate , ', Error number = ' , ierflg
          WRITE (Lun,*) ' N ' , N , ', EPS ' , eps , ', LENW ' , lenw
          WRITE (Lun,*) ' T ' , t
          WRITE (Lun,*) ' Y(1) ' , y(1)
          WRITE (Lun,*) ' Y(2) ' , y(2)
          WRITE (Lun,*) ' Y(3) ' , y(3)
          WRITE (Lun,*) ' Number of steps taken is  ' , nstep
          WRITE (Lun,*) ' Number of evaluations of the right hand side is  ' , 
     &                  nfe
          WRITE (Lun,*) ' Number of evaluations of the Jacobian matrix is  ' , 
     &                  nje
          WRITE (Lun,'(/)')
        ENDIF
        Ipass = 0
      ELSEIF ( ABS(0.620174E0-ABS(y(1)))>eps**(2.E0/3.E0).OR.
     &         ABS(0.392232E0-ABS(y(2)))>eps**(2.E0/3.E0).OR.ABS(1.E0-ABS(y(3)))
     &         >eps**(2.E0/3.E0) ) THEN
        IF ( Kprint==1 ) THEN
          WRITE (Lun,
     &'('' CDRIV1:The solution determined is not '',         ''accurate enough.'
     &' //)')
        ELSEIF ( Kprint==2 ) THEN
          WRITE (Lun,
     &'('' CDRIV1:The solution determined is not '',         ''accurate enough.'
     &')')
          WRITE (Lun,
     &'('' The values of parameters, results, and '',        ''statistical quant
     &ities are:'')')
          WRITE (Lun,*) ' EPS = ' , eps
          WRITE (Lun,*) ' T ' , t
          WRITE (Lun,*) ' Y(1) ' , y(1)
          WRITE (Lun,*) ' Y(2) ' , y(2)
          WRITE (Lun,*) ' Y(3) ' , y(3)
          WRITE (Lun,*) ' Number of steps taken is  ' , nstep
          WRITE (Lun,*) ' Number of evaluations of the right hand side is  ' , 
     &                  nfe
          WRITE (Lun,*) ' Number of evaluations of the Jacobian matrix is  ' , 
     &                  nje
          WRITE (Lun,'(/)')
        ENDIF
        Ipass = 0
      ELSEIF ( Kprint==2 ) THEN
        WRITE (Lun,
     &'('' CDRIV1:The solution determined met '',            ''the expected valu
     &es.'' //)')
      ELSEIF ( Kprint==3 ) THEN
        WRITE (Lun,
     &'('' CDRIV1:The solution determined met '',            ''the expected valu
     &es.'')')
        WRITE (Lun,'('' The values of results are '')')
        WRITE (Lun,*) ' T ' , t
        WRITE (Lun,*) ' Y(1) ' , y(1)
        WRITE (Lun,*) ' Y(2) ' , y(2)
        WRITE (Lun,*) ' Y(3) ' , y(3)
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
          WRITE (Lun,
     &'('' CDRIV1:An invalid parameter has not '',           ''been correctly de
     &tected.'' //)')
        ELSEIF ( Kprint>=2 ) THEN
          WRITE (Lun,
     &'('' CDRIV1:An invalid parameter has not '',           ''been correctly de
     &tected.'')')
          WRITE (Lun,*) ' The value of N was set to ' , nx
          WRITE (Lun,*) ' MSTATE = ' , mstate , ', Error number = ' , ierflg
          WRITE (Lun,
     &'('' The values of parameters, results, and '',        ''statistical quant
     &ities are:'')')
          WRITE (Lun,*) ' EPS ' , eps , ', LENW ' , lenw
          WRITE (Lun,*) ' T ' , t
          WRITE (Lun,*) ' Y(1) ' , y(1)
          WRITE (Lun,*) ' Y(2) ' , y(2)
          WRITE (Lun,*) ' Y(3) ' , y(3)
          WRITE (Lun,*) ' Number of steps taken is  ' , nstep
          WRITE (Lun,*) ' Number of evaluations of the right hand side is  ' , 
     &                  nfe
          WRITE (Lun,*) ' Number of evaluations of the Jacobian matrix is  ' , 
     &                  nje
          WRITE (Lun,'(/)')
        ENDIF
        Ipass = 0
      ELSEIF ( Kprint==2 ) THEN
        WRITE (Lun,
     &'('' CDRIV1:An invalid parameter has been '',          ''correctly detecte
     &d.'' //)')
      ELSEIF ( Kprint==3 ) THEN
        WRITE (Lun,
     &'('' CDRIV1:An invalid parameter has been '',          ''correctly detecte
     &d.'')')
        WRITE (Lun,*) ' The value of N was set to ' , nx
        WRITE (Lun,*) ' MSTATE = ' , mstate , ', Error number = ' , ierflg
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
      CALL CDRIV2(N,t,y,CDF,tout,mstate,NROOT,eps,ewt,mint,work,lenw,iwork,
     &            leniw,CDF,ierflg)
      nstep = iwork(3)
      nfe = iwork(4)
      nje = iwork(5)
      IF ( mstate/=2 ) THEN
        IF ( Kprint==1 ) THEN
          WRITE (Lun,
     & '('' While using CDRIV2, a solution was not '',        ''obtained.'' //)'
     & )
        ELSEIF ( Kprint>=2 ) THEN
          WRITE (Lun,
     &    '('' While using CDRIV2, a solution was not '',        ''obtained.'')'
     &    )
          WRITE (Lun,*) ' MSTATE = ' , mstate , ', Error number = ' , ierflg
          WRITE (Lun,
     &'('' The values of parameters, results, and '',        ''statistical quant
     &ities are:'')')
          WRITE (Lun,*) ' EPS = ' , eps , ', EWT ' , ewt
          WRITE (Lun,*) ' MINT = ' , mint , ', LENW ' , lenw , ', LENIW ' , 
     &                  leniw
          WRITE (Lun,*) ' T ' , t
          WRITE (Lun,*) ' Y(1) ' , y(1)
          WRITE (Lun,*) ' Y(2) ' , y(2)
          WRITE (Lun,*) ' Y(3) ' , y(3)
          WRITE (Lun,*) ' Number of steps taken is  ' , nstep
          WRITE (Lun,*) ' Number of evaluations of the right hand side is  ' , 
     &                  nfe
          WRITE (Lun,*) ' Number of evaluations of the Jacobian matrix is  ' , 
     &                  nje
          WRITE (Lun,'(/)')
        ENDIF
        Ipass = 0
      ELSEIF ( ABS(0.620174E0-ABS(y(1)))>eps**(2.E0/3.E0).OR.
     &         ABS(0.392232E0-ABS(y(2)))>eps**(2.E0/3.E0).OR.ABS(1.E0-ABS(y(3)))
     &         >eps**(2.E0/3.E0) ) THEN
        IF ( Kprint==1 ) THEN
          WRITE (Lun,
     &'('' CDRIV2:The solution determined is not '',         ''accurate enough. 
     &//'')')
        ELSEIF ( Kprint==2 ) THEN
          WRITE (Lun,
     &'('' CDRIV2:The solution determined is not '',         ''accurate enough.'
     &')')
          WRITE (Lun,
     &'('' The values of parameters, results, and '',        ''statistical quant
     &ities are:'')')
          WRITE (Lun,*) ' EPS = ' , eps , ', EWT = ' , ewt
          WRITE (Lun,*) ' T ' , t
          WRITE (Lun,*) ' Y(1) ' , y(1)
          WRITE (Lun,*) ' Y(2) ' , y(2)
          WRITE (Lun,*) ' Y(3) ' , y(3)
          WRITE (Lun,*) ' Number of steps taken is  ' , nstep
          WRITE (Lun,*) ' Number of evaluations of the right hand side is  ' , 
     &                  nfe
          WRITE (Lun,*) ' Number of evaluations of the Jacobian matrix is  ' , 
     &                  nje
          WRITE (Lun,'(/)')
        ENDIF
        Ipass = 0
      ELSEIF ( Kprint==2 ) THEN
        WRITE (Lun,
     &'('' CDRIV2:The solution determined met '',            ''the expected valu
     &es.'' //)')
      ELSEIF ( Kprint==3 ) THEN
        WRITE (Lun,
     &'('' CDRIV2:The solution determined met '',            ''the expected valu
     &es.'')')
        WRITE (Lun,'('' The values of results are '')')
        WRITE (Lun,*) ' T ' , t
        WRITE (Lun,*) ' Y(1) ' , y(1)
        WRITE (Lun,*) ' Y(2) ' , y(2)
        WRITE (Lun,*) ' Y(3) ' , y(3)
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
      CALL CDRIV2(N,t,y,CDF,tout,mstate,NROOT,eps,ewt,mint,work,lenwx,iwork,
     &            leniw,CDF,ierflg)
      IF ( ierflg/=32 ) THEN
        IF ( Kprint==1 ) THEN
          WRITE (Lun,
     &'('' CDRIV2:An invalid parameter has not '',           ''been correctly de
     &tected.'' //)')
        ELSEIF ( Kprint>=2 ) THEN
          WRITE (Lun,
     &'('' CDRIV2:An invalid parameter has not '',           ''been correctly de
     &tected.'')')
          WRITE (Lun,*) ' The value of LENW was set to ' , lenwx
          WRITE (Lun,*) ' MSTATE = ' , mstate , ', Error number = ' , ierflg
          WRITE (Lun,
     &'('' The values of parameters, results, and '',        ''statistical quant
     &ities are:'')')
          WRITE (Lun,*) ' EPS ' , eps , ', MINT ' , mint , ', LENW ' , lenw , 
     &                  ', LENIW ' , leniw
          WRITE (Lun,*) ' T ' , t
          WRITE (Lun,*) ' Y(1) ' , y(1)
          WRITE (Lun,*) ' Y(2) ' , y(2)
          WRITE (Lun,*) ' Y(3) ' , y(3)
          WRITE (Lun,*) ' Number of steps taken is  ' , nstep
          WRITE (Lun,*) ' Number of evaluations of the right hand side is  ' , 
     &                  nfe
          WRITE (Lun,*) ' Number of evaluations of the Jacobian matrix is  ' , 
     &                  nje
          WRITE (Lun,'(/)')
        ENDIF
        Ipass = 0
      ELSEIF ( Kprint==2 ) THEN
        WRITE (Lun,
     &'('' CDRIV2:An invalid parameter has been '',          ''correctly detecte
     &d.'' //)')
      ELSEIF ( Kprint==3 ) THEN
        WRITE (Lun,
     &'('' CDRIV2:An invalid parameter has been '',          ''correctly detecte
     &d.'')')
        WRITE (Lun,*) ' The value of LENW was set to ' , lenwx
        WRITE (Lun,*) ' MSTATE = ' , mstate , ', Error number = ' , ierflg
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
      CALL CDRIV3(N,t,y,CDF,nstate,tout,NTASK,NROOT,eps,ewt,IERROR,mint,MITER,
     &            IMPL,ML,MU,MXORD,HMAX,work,lenw,iwork,leniw,CDF,CDF,nde,
     &            MXSTEP,CDF,CDF,ierflg)
      nstep = iwork(3)
      nfe = iwork(4)
      nje = iwork(5)
      IF ( nstate/=2 ) THEN
        IF ( Kprint==1 ) THEN
          WRITE (Lun,
     & '('' While using CDRIV3, a solution was not '',        ''obtained.'' //)'
     & )
        ELSEIF ( Kprint>=2 ) THEN
          WRITE (Lun,
     &    '('' While using CDRIV3, a solution was not '',        ''obtained.'')'
     &    )
          WRITE (Lun,*) ' MSTATE = ' , mstate , ', Error number = ' , ierflg
          WRITE (Lun,
     &'('' The values of parameters, results, and '',        ''statistical quant
     &ities are:'')')
          WRITE (Lun,*) ' EPS = ' , eps , ', EWT = ' , ewt , ', IERROR = ' , 
     &                  IERROR
          WRITE (Lun,*) ' MINT = ' , mint , ', MITER = ' , MITER , ', IMPL = ' , 
     &                  IMPL
          WRITE (Lun,*) ' T ' , t
          WRITE (Lun,*) ' Y(1) ' , y(1)
          WRITE (Lun,*) ' Y(2) ' , y(2)
          WRITE (Lun,*) ' Y(3) ' , y(3)
          WRITE (Lun,*) ' Number of steps taken is  ' , nstep
          WRITE (Lun,*) ' Number of evaluations of the right hand side is  ' , 
     &                  nfe
          WRITE (Lun,*) ' Number of evaluations of the Jacobian matrix is  ' , 
     &                  nje
          WRITE (Lun,'(/)')
        ENDIF
        Ipass = 0
      ELSEIF ( ABS(0.620174E0-ABS(y(1)))>eps**(2.E0/3.E0).OR.
     &         ABS(0.392232E0-ABS(y(2)))>eps**(2.E0/3.E0).OR.ABS(1.E0-ABS(y(3)))
     &         >eps**(2.E0/3.E0) ) THEN
        IF ( Kprint==1 ) THEN
          WRITE (Lun,
     &'('' CDRIV3:The solution determined is not '',         ''accurate enough.'
     &' //)')
        ELSEIF ( Kprint>=2 ) THEN
          WRITE (Lun,
     &'('' CDRIV3:The solution determined is not '',         ''accurate enough.'
     &')')
          WRITE (Lun,
     &'('' The values of parameters, results, and '',        ''statistical quant
     &ities are:'')')
          WRITE (Lun,*) ' EPS = ' , eps , ', EWT = ' , ewt , ', IERROR = ' , 
     &                  IERROR
          WRITE (Lun,*) ' MINT = ' , mint , ', MITER = ' , MITER , ', IMPL = ' , 
     &                  IMPL
          WRITE (Lun,*) ' T ' , t
          WRITE (Lun,*) ' Y(1) ' , y(1)
          WRITE (Lun,*) ' Y(2) ' , y(2)
          WRITE (Lun,*) ' Y(3) ' , y(3)
          WRITE (Lun,*) ' Number of steps taken is  ' , nstep
          WRITE (Lun,*) ' Number of evaluations of the right hand side is  ' , 
     &                  nfe
          WRITE (Lun,*) ' Number of evaluations of the Jacobian matrix is  ' , 
     &                  nje
          WRITE (Lun,'(/)')
        ENDIF
        Ipass = 0
      ELSEIF ( Kprint==2 ) THEN
        WRITE (Lun,
     &'('' CDRIV3:The solution determined met '',            ''the expected valu
     &es.'' //)')
      ELSEIF ( Kprint==3 ) THEN
        WRITE (Lun,
     &'('' CDRIV3:The solution determined met '',            ''the expected valu
     &es.'')')
        WRITE (Lun,'('' The values of results are '')')
        WRITE (Lun,*) ' T ' , t
        WRITE (Lun,*) ' Y(1) ' , y(1)
        WRITE (Lun,*) ' Y(2) ' , y(2)
        WRITE (Lun,*) ' Y(3) ' , y(3)
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
      CALL CDRIV3(N,t,y,CDF,nstate,tout,NTASK,NROOT,eps,ewt,IERROR,mint,MITER,
     &            IMPL,ML,MU,MXORD,HMAX,work,lenw,iwork,leniwx,CDF,CDF,nde,
     &            MXSTEP,CDF,CDF,ierflg)
      IF ( ierflg/=33 ) THEN
        IF ( Kprint==1 ) THEN
          WRITE (Lun,
     &'('' CDRIV3:An invalid parameter has not '',           ''been correctly de
     &tected.'' //)')
        ELSEIF ( Kprint>=2 ) THEN
          WRITE (Lun,
     &'('' CDRIV3:An invalid parameter has not '',           ''been correctly de
     &tected.'')')
          WRITE (Lun,*) ' The value of LENIW was set to ' , leniwx
          WRITE (Lun,*) ' NSTATE = ' , nstate , ', Error number = ' , ierflg
          WRITE (Lun,
     &'('' The values of parameters, results, and '',        ''statistical quant
     &ities are:'')')
          WRITE (Lun,*) ' EPS = ' , eps , ', EWT = ' , ewt , ', IERROR = ' , 
     &                  IERROR
          WRITE (Lun,*) ' MINT = ' , mint , ', MITER = ' , MITER , ', IMPL = ' , 
     &                  IMPL
          WRITE (Lun,*) ' T ' , t
          WRITE (Lun,*) ' Y(1) ' , y(1)
          WRITE (Lun,*) ' Y(2) ' , y(2)
          WRITE (Lun,*) ' Y(3) ' , y(3)
          WRITE (Lun,*) ' Number of steps taken is  ' , nstep
          WRITE (Lun,*) ' Number of evaluations of the right hand side is  ' , 
     &                  nfe
          WRITE (Lun,*) ' Number of evaluations of the Jacobian matrix is  ' , 
     &                  nje
          WRITE (Lun,'(/)')
        ENDIF
        Ipass = 0
      ELSEIF ( Kprint==2 ) THEN
        WRITE (Lun,
     &'('' CDRIV3:An invalid parameter has been '',          ''correctly detecte
     &d.'' //)')
      ELSEIF ( Kprint==3 ) THEN
        WRITE (Lun,
     &'('' CDRIV3:An invalid parameter has been '',          ''correctly detecte
     &d.'')')
        WRITE (Lun,*) ' The value of LENIW was set to ' , leniwx
        WRITE (Lun,*) ' NSTATE = ' , nstate , ', Error number = ' , ierflg
        WRITE (Lun,'(/)')
      ENDIF
      CALL XERCLR
      END SUBROUTINE CDQCK
