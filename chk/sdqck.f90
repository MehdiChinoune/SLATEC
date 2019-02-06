!*==SDQCK.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK SDQCK
      SUBROUTINE SDQCK(Lun,Kprint,Ipass)
      IMPLICIT NONE
!*--SDQCK5
!*** Start of declarations inserted by SPAG
      REAL SDF
!*** End of declarations inserted by SPAG
!***BEGIN PROLOGUE  SDQCK
!***PURPOSE  Quick check for SLATEC routines SDRIV1, SDRIV2 and SDRIV3.
!***LIBRARY   SLATEC (SDRIVE)
!***CATEGORY  I1A2, I1A1B
!***TYPE      SINGLE PRECISION (SDQCK-S, DDQCK-D, CDQCK-C)
!***KEYWORDS  QUICK CHECK, SDRIV1, SDRIV2, SDRIV3
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
!***ROUTINES CALLED  R1MACH, SDF, SDRIV1, SDRIV2, SDRIV3, XERCLR
!***REVISION HISTORY  (YYMMDD)
!   890405  DATE WRITTEN
!   890405  Revised to meet SLATEC standards.
!***END PROLOGUE  SDQCK
      EXTERNAL SDF
      REAL ALFA , eps , ewt(1) , HMAX , R1MACH , t , tout
      INTEGER ierflg , IERROR , IMPL , Ipass , Kprint , leniw , leniwx , lenw , 
     &        LENWMX , lenwx , LIWMX , Lun , mint , MITER , ML , mstate , MU , 
     &        MXORD , MXSTEP , N , nde , nfe , nje , NROOT , nstate , nstep , 
     &        NTASK , nx
      PARAMETER (ALFA=1.E0,HMAX=15.E0,IERROR=3,IMPL=0,LENWMX=342,LIWMX=53,
     &           MITER=5,ML=2,MU=2,MXORD=5,MXSTEP=1000,N=3,NROOT=0,NTASK=1)
      REAL work(LENWMX) , y(N+1)
      INTEGER iwork(LIWMX)
      DATA ewt(1)/.00001E0/
!***FIRST EXECUTABLE STATEMENT  SDQCK
      eps = R1MACH(4)**(1.E0/3.E0)
      Ipass = 1
!                                            Exercise SDRIV1 for problem
!                                            with known solution.
      y(4) = ALFA
      t = 0.E0
      y(1) = 10.E0
      y(2) = 0.E0
      y(3) = 10.E0
      tout = 10.E0
      mstate = 1
      lenw = 342
      CALL SDRIV1(N,t,y,SDF,tout,mstate,eps,work,lenw,ierflg)
      nstep = work(lenw-(N+50)+3)
      nfe = work(lenw-(N+50)+4)
      nje = work(lenw-(N+50)+5)
      IF ( mstate/=2 ) THEN
        IF ( Kprint==1 ) THEN
          WRITE (Lun,
     & '('' While using SDRIV1, a solution was not '',        ''obtained.'' //)'
     & )
        ELSEIF ( Kprint>=2 ) THEN
          WRITE (Lun,
     &    '('' While using SDRIV1, a solution was not '',        ''obtained.'')'
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
          WRITE (Lun,'(//)')
        ENDIF
        Ipass = 0
      ELSEIF ( ABS(1.E0-y(1)*1.5E0)>eps**(2.E0/3.E0).OR.ABS(1.E0-y(2)*3.E0)
     &         >eps**(2.E0/3.E0).OR.ABS(1.E0-y(3))>eps**(2.E0/3.E0) ) THEN
        IF ( Kprint==1 ) THEN
          WRITE (Lun,
     &'('' SDRIV1:The solution determined is not '',         ''accurate enough.'
     &' //)')
        ELSEIF ( Kprint==2 ) THEN
          WRITE (Lun,
     &'('' SDRIV1:The solution determined is not '',         ''accurate enough.'
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
          WRITE (Lun,'(//)')
        ENDIF
        Ipass = 0
      ELSEIF ( Kprint==2 ) THEN
        WRITE (Lun,
     &'('' SDRIV1:The solution determined met '',            ''the expected valu
     &es.'' //)')
      ELSEIF ( Kprint==3 ) THEN
        WRITE (Lun,
     &'('' SDRIV1:The solution determined met '',            ''the expected valu
     &es.'')')
        WRITE (Lun,'('' The values of results are '')')
        WRITE (Lun,*) ' T ' , t
        WRITE (Lun,*) ' Y(1) ' , y(1)
        WRITE (Lun,*) ' Y(2) ' , y(2)
        WRITE (Lun,*) ' Y(3) ' , y(3)
        WRITE (Lun,'(/)')
      ENDIF
      CALL XERCLR
!                                         Run SDRIV1 with invalid input.
      nx = 201
      t = 0.E0
      y(1) = 10.E0
      y(2) = 0.E0
      y(3) = 10.E0
      y(4) = ALFA
      tout = 10.E0
      mstate = 1
      lenw = 342
      CALL SDRIV1(nx,t,y,SDF,tout,mstate,eps,work,lenw,ierflg)
      IF ( ierflg/=21 ) THEN
        IF ( Kprint==1 ) THEN
          WRITE (Lun,
     &'('' SDRIV1:An invalid parameter has not '',           ''been correctly de
     &tected.'' //)')
        ELSEIF ( Kprint>=2 ) THEN
          WRITE (Lun,
     &'('' SDRIV1:An invalid parameter has not '',           ''been correctly de
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
          WRITE (Lun,'(//)')
        ENDIF
        Ipass = 0
      ELSEIF ( Kprint==2 ) THEN
        WRITE (Lun,
     &'('' SDRIV1:An invalid parameter has been '',          ''correctly detecte
     &d.'' //)')
      ELSEIF ( Kprint==3 ) THEN
        WRITE (Lun,
     &'('' SDRIV1:An invalid parameter has been '',          ''correctly detecte
     &d.'')')
        WRITE (Lun,*) ' The value of N was set to ' , nx
        WRITE (Lun,*) ' MSTATE = ' , mstate , ', Error number = ' , ierflg
        WRITE (Lun,'(/)')
      ENDIF
      CALL XERCLR
!                                            Exercise SDRIV2 for problem
!                                            with known solution.
      t = 0.E0
      y(1) = 10.E0
      y(2) = 0.E0
      y(3) = 10.E0
      y(4) = ALFA
      mstate = 1
      tout = 10.E0
      mint = 1
      lenw = 298
      leniw = 50
      CALL SDRIV2(N,t,y,SDF,tout,mstate,NROOT,eps,ewt,mint,work,lenw,iwork,
     &            leniw,SDF,ierflg)
      nstep = iwork(3)
      nfe = iwork(4)
      nje = iwork(5)
      IF ( mstate/=2 ) THEN
        IF ( Kprint==1 ) THEN
          WRITE (Lun,
     & '('' While using SDRIV2, a solution was not '',        ''obtained.'' //)'
     & )
        ELSEIF ( Kprint>=2 ) THEN
          WRITE (Lun,
     &    '('' While using SDRIV2, a solution was not '',        ''obtained.'')'
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
          WRITE (Lun,'(//)')
        ENDIF
        Ipass = 0
      ELSEIF ( ABS(1.E0-y(1)*1.5E0)>eps**(2.E0/3.E0).OR.ABS(1.E0-y(2)*3.E0)
     &         >eps**(2.E0/3.E0).OR.ABS(1.E0-y(3))>eps**(2.E0/3.E0) ) THEN
        IF ( Kprint==1 ) THEN
          WRITE (Lun,
     &'('' SDRIV2:The solution determined is not '',         ''accurate enough.'
     &' //)')
        ELSEIF ( Kprint==2 ) THEN
          WRITE (Lun,
     &'('' SDRIV2:The solution determined is not '',         ''accurate enough.'
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
          WRITE (Lun,'(//)')
        ENDIF
        Ipass = 0
      ELSEIF ( Kprint==2 ) THEN
        WRITE (Lun,
     &'('' SDRIV2:The solution determined met '',            ''the expected valu
     &es.'' //)')
      ELSEIF ( Kprint==3 ) THEN
        WRITE (Lun,
     &'('' SDRIV2:The solution determined met '',            ''the expected valu
     &es.'')')
        WRITE (Lun,'('' The values of results are '')')
        WRITE (Lun,*) ' T ' , t
        WRITE (Lun,*) ' Y(1) ' , y(1)
        WRITE (Lun,*) ' Y(2) ' , y(2)
        WRITE (Lun,*) ' Y(3) ' , y(3)
        WRITE (Lun,'(/)')
      ENDIF
      CALL XERCLR
!                                         Run SDRIV2 with invalid input.
      t = 0.E0
      y(1) = 10.E0
      y(2) = 0.E0
      y(3) = 10.E0
      y(4) = ALFA
      tout = 10.E0
      mstate = 1
      mint = 1
      lenwx = 1
      leniw = 50
      CALL SDRIV2(N,t,y,SDF,tout,mstate,NROOT,eps,ewt,mint,work,lenwx,iwork,
     &            leniw,SDF,ierflg)
      IF ( ierflg/=32 ) THEN
        IF ( Kprint==1 ) THEN
          WRITE (Lun,
     &'('' SDRIV2:An invalid parameter has not '',           ''been correctly de
     &tected.'' //)')
        ELSEIF ( Kprint>=2 ) THEN
          WRITE (Lun,
     &'('' SDRIV2:An invalid parameter has not '',           ''been correctly de
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
          WRITE (Lun,'(//)')
        ENDIF
        Ipass = 0
      ELSEIF ( Kprint==2 ) THEN
        WRITE (Lun,
     &'('' SDRIV2:An invalid parameter has been '',          ''correctly detecte
     &d.'' //)')
      ELSEIF ( Kprint==3 ) THEN
        WRITE (Lun,
     &'('' SDRIV2:An invalid parameter has been '',          ''correctly detecte
     &d.'')')
        WRITE (Lun,*) ' The value of LENW was set to ' , lenwx
        WRITE (Lun,*) ' MSTATE = ' , mstate , ', Error number = ' , ierflg
        WRITE (Lun,'(/)')
      ENDIF
      CALL XERCLR
!                                            Exercise SDRIV3 for problem
!                                            with known solution.
      t = 0.E0
      y(1) = 10.E0
      y(2) = 0.E0
      y(3) = 10.E0
      y(4) = ALFA
      nstate = 1
      tout = 10.E0
      mint = 2
      lenw = 301
      leniw = 53
      CALL SDRIV3(N,t,y,SDF,nstate,tout,NTASK,NROOT,eps,ewt,IERROR,mint,MITER,
     &            IMPL,ML,MU,MXORD,HMAX,work,lenw,iwork,leniw,SDF,SDF,nde,
     &            MXSTEP,SDF,SDF,ierflg)
      nstep = iwork(3)
      nfe = iwork(4)
      nje = iwork(5)
      IF ( nstate/=2 ) THEN
        IF ( Kprint==1 ) THEN
          WRITE (Lun,
     & '('' While using SDRIV3, a solution was not '',        ''obtained.'' //)'
     & )
        ELSEIF ( Kprint>=2 ) THEN
          WRITE (Lun,
     &    '('' While using SDRIV3, a solution was not '',        ''obtained.'')'
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
          WRITE (Lun,'(//)')
        ENDIF
        Ipass = 0
      ELSEIF ( ABS(1.E0-y(1)*1.5E0)>eps**(2.E0/3.E0).OR.ABS(1.E0-y(2)*3.E0)
     &         >eps**(2.E0/3.E0).OR.ABS(1.E0-y(3))>eps**(2.E0/3.E0) ) THEN
        IF ( Kprint==1 ) THEN
          WRITE (Lun,
     &'('' SDRIV3:The solution determined is not '',         ''accurate enough.'
     &' //)')
        ELSEIF ( Kprint>=2 ) THEN
          WRITE (Lun,
     &'('' SDRIV3:The solution determined is not '',         ''accurate enough.'
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
          WRITE (Lun,'(//)')
        ENDIF
        Ipass = 0
      ELSEIF ( Kprint==2 ) THEN
        WRITE (Lun,
     &'('' SDRIV3:The solution determined met '',            ''the expected valu
     &es.'' //)')
      ELSEIF ( Kprint==3 ) THEN
        WRITE (Lun,
     &'('' SDRIV3:The solution determined met '',            ''the expected valu
     &es.'')')
        WRITE (Lun,'('' The values of results are '')')
        WRITE (Lun,*) ' T ' , t
        WRITE (Lun,*) ' Y(1) ' , y(1)
        WRITE (Lun,*) ' Y(2) ' , y(2)
        WRITE (Lun,*) ' Y(3) ' , y(3)
        WRITE (Lun,'(/)')
      ENDIF
      CALL XERCLR
!                                         Run SDRIV3 with invalid input.
      t = 0.E0
      y(1) = 10.E0
      y(2) = 0.E0
      y(3) = 10.E0
      y(4) = ALFA
      nstate = 1
      tout = 10.E0
      mint = 2
      lenw = 301
      leniwx = 1
      CALL SDRIV3(N,t,y,SDF,nstate,tout,NTASK,NROOT,eps,ewt,IERROR,mint,MITER,
     &            IMPL,ML,MU,MXORD,HMAX,work,lenw,iwork,leniwx,SDF,SDF,nde,
     &            MXSTEP,SDF,SDF,ierflg)
      IF ( ierflg/=33 ) THEN
        IF ( Kprint==1 ) THEN
          WRITE (Lun,
     &'('' SDRIV3:An invalid parameter has not '',           ''been correctly de
     &tected.'' //)')
        ELSEIF ( Kprint>=2 ) THEN
          WRITE (Lun,
     &'('' SDRIV3:An invalid parameter has not '',           ''been correctly de
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
          WRITE (Lun,'(//)')
        ENDIF
        Ipass = 0
      ELSEIF ( Kprint==2 ) THEN
        WRITE (Lun,
     &'('' SDRIV3:An invalid parameter has been '',          ''correctly detecte
     &d.'' //)')
      ELSEIF ( Kprint==3 ) THEN
        WRITE (Lun,
     &'('' SDRIV3:An invalid parameter has been '',          ''correctly detecte
     &d.'')')
        WRITE (Lun,*) ' The value of LENIW was set to ' , leniwx
        WRITE (Lun,*) ' NSTATE = ' , nstate , ', Error number = ' , ierflg
        WRITE (Lun,'(/)')
      ENDIF
      CALL XERCLR
      END SUBROUTINE SDQCK
