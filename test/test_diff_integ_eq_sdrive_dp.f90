MODULE TEST46_MOD
  USE service, ONLY : SP, DP
  IMPLICIT NONE

CONTAINS
  !** DDQCK
  SUBROUTINE DDQCK(Lun,Kprint,Ipass)
    !> Quick check for SLATEC routines DDRIV1, DDRIV2 and DDRIV3.
    !***
    ! **Library:**   SLATEC (SDRIVE)
    !***
    ! **Category:**  I1A2, I1A1B
    !***
    ! **Type:**      DOUBLE PRECISION (SDQCK-S, DDQCK-D, CDQCK-C)
    !***
    ! **Keywords:**  DDRIV1, DDRIV2, DDRIV3, QUICK CHECK, SDRIVE
    !***
    ! **Author:**  Kahaner, D. K., (NIST)
    !             National Institute of Standards and Technology
    !             Gaithersburg, MD  20899
    !           Sutherland, C. D., (LANL)
    !             Mail Stop D466
    !             Los Alamos National Laboratory
    !             Los Alamos, NM  87545
    !***
    ! **Description:**
    !
    !  For assistance in determining the cause of a failure of these
    !  routines contact C. D. Sutherland at commercial telephone number
    !  (505)667-6949, FTS telephone number 8-843-6949, or electronic mail
    !  address CDS@LANL.GOV .
    !
    !***
    ! **Routines called:**  D1MACH, DDF, DDRIV1, DDRIV2, DDRIV3, XERCLR

    !* REVISION HISTORY  (YYMMDD)
    !   890405  DATE WRITTEN
    !   890405  Revised to meet SLATEC standards.
    USE service, ONLY : eps_dp
    USE diff_integ_eq, ONLY : DDRIV1, DDRIV2, DDRIV3
    !
    REAL(DP) :: eps, t, tout
    INTEGER :: ierflg, Ipass, Kprint, leniw, lenw, Lun, mint, &
      mstate, nde, nfe, nje, nstate, nstep
    REAL(DP), PARAMETER :: ALFA = 1._DP, HMAX = 15._DP
    INTEGER, PARAMETER :: IERROR = 3, IMPL = 0, LENWMX = 342, LIWMX = 53, &
      MITER = 5, ML = 2, MU = 2, MXORD = 5, MXSTEP = 1000, N = 3, NROOT = 0, NTASK = 1
    REAL(DP) :: work(LENWMX), y(N+1)
    INTEGER :: iwork(LIWMX)
    REAL(DP), PARAMETER :: ewt(1) = .00001_DP
    !* FIRST EXECUTABLE STATEMENT  DDQCK
    eps = eps_dp**(1._DP/3._DP)
    Ipass = 1
    !  Exercise DDRIV1 for problem with known solution.
    y(4) = ALFA
    t = 0._DP
    y(1) = 10._DP
    y(2) = 0._DP
    y(3) = 10._DP
    tout = 10._DP
    mstate = 1
    lenw = 342
    CALL DDRIV1(N,t,y,DDF,tout,mstate,eps,work,lenw,ierflg)
    nstep = INT( work(lenw-(N+50)+3) )
    nfe = INT( work(lenw-(N+50)+4) )
    nje = INT( work(lenw-(N+50)+5) )
    IF( mstate/=2 ) THEN
      IF( Kprint==1 ) THEN
        WRITE (Lun,&
          '('' While using DDRIV1, a solution was not obtained.'' //)')
      ELSEIF( Kprint>=2 ) THEN
        WRITE (Lun,&
          '('' While using DDRIV1, a solution was not obtained.'')')
        WRITE (Lun,&
          '('' The values of parameters, results, and statistical quantities are:'')')
        WRITE (Lun,*) ' MSTATE = ', mstate, ', Error number = ', ierflg
        WRITE (Lun,*) ' N ', N, ', EPS ', eps, ', LENW ', lenw
        WRITE (Lun,*) ' T ', t
        WRITE (Lun,*) ' Y(1) ', y(1)
        WRITE (Lun,*) ' Y(2) ', y(2)
        WRITE (Lun,*) ' Y(3) ', y(3)
        WRITE (Lun,*) ' Number of steps taken is  ', nstep
        WRITE (Lun,*) ' Number of evaluations of the right hand side is  ', nfe
        WRITE (Lun,*) ' Number of evaluations of the Jacobian matrix is  ', nje
        WRITE (Lun,'(//)')
      END IF
      Ipass = 0
    ELSEIF( ABS(1._DP-y(1)*1.5_DP)>eps**(2._DP/3._DP) .OR. ABS(1._DP-y(2)*3._DP)&
        >eps**(2._DP/3._DP) .OR. ABS(1._DP-y(3))>eps**(2._DP/3._DP) ) THEN
      IF( Kprint==1 ) THEN
        WRITE (Lun,&
          '('' DDRIV1:The solution determined is not  accurate enough.'' //)')
      ELSEIF( Kprint==2 ) THEN
        WRITE (Lun,&
          '('' DDRIV1:The solution determined is not accurate enough.'')')
        WRITE (Lun,&
          '('' The values of parameters, results, and statistical quantities are:'')')
        WRITE (Lun,*) ' EPS = ', eps
        WRITE (Lun,*) ' T ', t
        WRITE (Lun,*) ' Y(1) ', y(1)
        WRITE (Lun,*) ' Y(2) ', y(2)
        WRITE (Lun,*) ' Y(3) ', y(3)
        WRITE (Lun,*) ' Number of steps taken is  ', nstep
        WRITE (Lun,*) ' Number of evaluations of the right hand side is  ', nfe
        WRITE (Lun,*) ' Number of evaluations of the Jacobian matrix is  ', nje
        WRITE (Lun,'(//)')
      END IF
      Ipass = 0
    ELSEIF( Kprint==2 ) THEN
      WRITE (Lun,&
        '('' DDRIV1:The solution determined met the expected values.'' //)')
    ELSEIF( Kprint==3 ) THEN
      WRITE (Lun,&
        '('' DDRIV1:The solution determined met the expected values.'')')
      WRITE (Lun,'('' The values of results are '')')
      WRITE (Lun,*) ' T ', t
      WRITE (Lun,*) ' Y(1) ', y(1)
      WRITE (Lun,*) ' Y(2) ', y(2)
      WRITE (Lun,*) ' Y(3) ', y(3)
      WRITE (Lun,'(/)')
    END IF
    !  Run DDRIV1 with invalid input.
!    nx = 201
!    t = 0._DP
!    y(1) = 10._DP
!    y(2) = 0._DP
!    y(3) = 10._DP
!    y(4) = ALFA
!    tout = 10._DP
!    mstate = 1
!    lenw = 342
!    CALL DDRIV1(nx,t,y,DDF,tout,mstate,eps,work,lenw,ierflg)
!    IF( ierflg/=21 ) THEN
!      IF( Kprint==1 ) THEN
!        WRITE (Lun,&
!          '('' DDRIV1:An invalid parameter has not been correctly detected.'' //)')
!      ELSEIF( Kprint>=2 ) THEN
!        WRITE (Lun,&
!          '('' DDRIV1:An invalid parameter has not been correctly detected.'')')
!        WRITE (Lun,*) ' The value of N was set to ', nx
!        WRITE (Lun,*) ' MSTATE = ', mstate, ', Error number = ', ierflg
!        WRITE (Lun,&
!          '('' The values of parameters, results, and statistical quantities are:'')')
!        WRITE (Lun,*) ' EPS ', eps, ', LENW ', lenw
!        WRITE (Lun,*) ' T ', t
!        WRITE (Lun,*) ' Y(1) ', y(1)
!        WRITE (Lun,*) ' Y(2) ', y(2)
!        WRITE (Lun,*) ' Y(3) ', y(3)
!        WRITE (Lun,*) ' Number of steps taken is  ', nstep
!        WRITE (Lun,*) ' Number of evaluations of the right hand side is  ', nfe
!        WRITE (Lun,*) ' Number of evaluations of the Jacobian matrix is  ', nje
!        WRITE (Lun,'(//)')
!      END IF
!      Ipass = 0
!    ELSEIF( Kprint==2 ) THEN
!      WRITE (Lun,&
!        '('' DDRIV1:An invalid parameter has been correctly detected.'' //)')
!    ELSEIF( Kprint==3 ) THEN
!      WRITE (Lun,&
!        '('' DDRIV1:An invalid parameter has been correctly detected.'')')
!      WRITE (Lun,*) ' The value of N was set to ', nx
!      WRITE (Lun,*) ' MSTATE = ', mstate, ', Error number = ', ierflg
!      WRITE (Lun,'(/)')
!    END IF
!    num_xer = 0
    !  Exercise DDRIV2 for problem with known solution.
    t = 0._DP
    y(1) = 10._DP
    y(2) = 0._DP
    y(3) = 10._DP
    y(4) = ALFA
    mstate = 1
    tout = 10._DP
    mint = 1
    lenw = 298
    leniw = 50
    CALL DDRIV2(N,t,y,DDF,tout,mstate,NROOT,eps,ewt(1),mint,work,lenw,iwork,&
      leniw,dum_G,ierflg)
    nstep = iwork(3)
    nfe = iwork(4)
    nje = iwork(5)
    IF( mstate/=2 ) THEN
      IF( Kprint==1 ) THEN
        WRITE (Lun,&
          '('' While using DDRIV2, a solution was not obtained.'' //)' )
      ELSEIF( Kprint>=2 ) THEN
        WRITE (Lun,&
          '('' While using DDRIV2, a solution was not obtained.'')')
        WRITE (Lun,*) ' MSTATE = ', mstate, ', Error number = ', ierflg
        WRITE (Lun,&
          '('' The values of parameters, results, and statistical quant&
          &ities are:'')')
        WRITE (Lun,*) ' EPS = ', eps, ', EWT ', ewt
        WRITE (Lun,*) ' MINT = ', mint, ', LENW ', lenw, ', LENIW ', leniw
        WRITE (Lun,*) ' T ', t
        WRITE (Lun,*) ' Y(1) ', y(1)
        WRITE (Lun,*) ' Y(2) ', y(2)
        WRITE (Lun,*) ' Y(3) ', y(3)
        WRITE (Lun,*) ' Number of steps taken is  ', nstep
        WRITE (Lun,*) ' Number of evaluations of the right hand side is  ', nfe
        WRITE (Lun,*) ' Number of evaluations of the Jacobian matrix is  ', nje
        WRITE (Lun,'(//)')
      END IF
      Ipass = 0
    ELSEIF( ABS(1._DP-y(1)*1.5_DP)>eps**(2._DP/3._DP) .OR. ABS(1._DP-y(2)*3._DP)&
        >eps**(2._DP/3._DP) .OR. ABS(1._DP-y(3))>eps**(2._DP/3._DP) ) THEN
      IF( Kprint==1 ) THEN
        WRITE (Lun,&
          '('' DDRIV2:The solution determined is not accurate enough.'' //)')
      ELSEIF( Kprint==2 ) THEN
        WRITE (Lun,&
          '('' DDRIV2:The solution determined is not accurate enough.'')')
        WRITE (Lun,&
          '('' The values of parameters, results, and statistical quantities are:'')')
        WRITE (Lun,*) ' EPS = ', eps, ', EWT = ', ewt
        WRITE (Lun,*) ' T ', t
        WRITE (Lun,*) ' Y(1) ', y(1)
        WRITE (Lun,*) ' Y(2) ', y(2)
        WRITE (Lun,*) ' Y(3) ', y(3)
        WRITE (Lun,*) ' Number of steps taken is  ', nstep
        WRITE (Lun,*) ' Number of evaluations of the right hand side is  ', nfe
        WRITE (Lun,*) ' Number of evaluations of the Jacobian matrix is  ', nje
        WRITE (Lun,'(//)')
      END IF
      Ipass = 0
    ELSEIF( Kprint==2 ) THEN
      WRITE (Lun,&
        '('' DDRIV2:The solution determined met the expected values.'' //)')
    ELSEIF( Kprint==3 ) THEN
      WRITE (Lun,&
        '('' DDRIV2:The solution determined met the expected values.'')')
      WRITE (Lun,'('' The values of results are '')')
      WRITE (Lun,*) ' T ', t
      WRITE (Lun,*) ' Y(1) ', y(1)
      WRITE (Lun,*) ' Y(2) ', y(2)
      WRITE (Lun,*) ' Y(3) ', y(3)
      WRITE (Lun,'(/)')
    END IF
    !  Run DDRIV2 with invalid input.
!    t = 0._DP
!    y(1) = 10._DP
!    y(2) = 0._DP
!    y(3) = 10._DP
!    y(4) = ALFA
!    tout = 10._DP
!    mstate = 1
!    mint = 1
!    lenwx = 1
!    leniw = 50
!    CALL DDRIV2(N,t,y,DDF,tout,mstate,NROOT,eps,ewt(1),mint,work,lenwx,iwork,&
!      leniw,dum_G,ierflg)
!    IF( ierflg/=32 ) THEN
!      IF( Kprint==1 ) THEN
!        WRITE (Lun,&
!          '('' DDRIV2:An invalid parameter has not been correctly detected.'' //)')
!      ELSEIF( Kprint>=2 ) THEN
!        WRITE (Lun,&
!          '('' DDRIV2:An invalid parameter has not been correctly detected.'')')
!        WRITE (Lun,*) ' The value of LENW was set to ', lenwx
!        WRITE (Lun,*) ' MSTATE = ', mstate, ', Error number = ', ierflg
!        WRITE (Lun,&
!          '('' The values of parameters, results, and statistical quantities are:'')')
!        WRITE (Lun,*) ' EPS ', eps, ', MINT ', mint, ', LENW ', lenw, &
!          ', LENIW ', leniw
!        WRITE (Lun,*) ' T ', t
!        WRITE (Lun,*) ' Y(1) ', y(1)
!        WRITE (Lun,*) ' Y(2) ', y(2)
!        WRITE (Lun,*) ' Y(3) ', y(3)
!        WRITE (Lun,*) ' Number of steps taken is  ', nstep
!        WRITE (Lun,*) ' Number of evaluations of the right hand side is  ', nfe
!        WRITE (Lun,*) ' Number of evaluations of the Jacobian matrix is  ', nje
!        WRITE (Lun,'(//)')
!      END IF
!      Ipass = 0
!    ELSEIF( Kprint==2 ) THEN
!      WRITE (Lun,&
!        '('' DDRIV2:An invalid parameter has been correctly detected.'' //)')
!    ELSEIF( Kprint==3 ) THEN
!      WRITE (Lun,&
!        '('' DDRIV2:An invalid parameter has been correctly detected.'')')
!      WRITE (Lun,*) ' The value of LENW was set to ', lenwx
!      WRITE (Lun,*) ' MSTATE = ', mstate, ', Error number = ', ierflg
!      WRITE (Lun,'(/)')
!    END IF
!    num_xer = 0
    !  Exercise DDRIV3 for problem with known solution.
    t = 0._DP
    y(1) = 10._DP
    y(2) = 0._DP
    y(3) = 10._DP
    y(4) = ALFA
    nstate = 1
    tout = 10._DP
    mint = 2
    lenw = 301
    leniw = 53
    CALL DDRIV3(N,t,y,DDF,nstate,tout,NTASK,NROOT,eps,ewt,IERROR,mint,MITER,&
      IMPL,ML,MU,MXORD,HMAX,work,lenw,iwork,leniw,dum_JACOBN,dum_FA,nde,&
      MXSTEP,dum_G,dum_USERS,ierflg)
    nstep = iwork(3)
    nfe = iwork(4)
    nje = iwork(5)
    IF( nstate/=2 ) THEN
      IF( Kprint==1 ) THEN
        WRITE (Lun,&
          '('' While using DDRIV3, a solution was not obtained.'' //)' )
      ELSEIF( Kprint>=2 ) THEN
        WRITE (Lun,&
          '('' While using DDRIV3, a solution was not obtained.'')' )
        WRITE (Lun,*) ' MSTATE = ', mstate, ', Error number = ', ierflg
        WRITE (Lun,&
          '('' The values of parameters, results, and statistical quantities are:'')')
        WRITE (Lun,*) ' EPS = ', eps, ', EWT = ', ewt, ', IERROR = ', IERROR
        WRITE (Lun,*) ' MINT = ', mint, ', MITER = ', MITER, ', IMPL = ', IMPL
        WRITE (Lun,*) ' T ', t
        WRITE (Lun,*) ' Y(1) ', y(1)
        WRITE (Lun,*) ' Y(2) ', y(2)
        WRITE (Lun,*) ' Y(3) ', y(3)
        WRITE (Lun,*) ' Number of steps taken is  ', nstep
        WRITE (Lun,*) ' Number of evaluations of the right hand side is  ', nfe
        WRITE (Lun,*) ' Number of evaluations of the Jacobian matrix is  ', nje
        WRITE (Lun,'(//)')
      END IF
      Ipass = 0
    ELSEIF( ABS(1._DP-y(1)*1.5_DP)>eps**(2._DP/3._DP) .OR. ABS(1._DP-y(2)*3._DP)&
        >eps**(2._DP/3._DP) .OR. ABS(1._DP-y(3))>eps**(2._DP/3._DP) ) THEN
      IF( Kprint==1 ) THEN
        WRITE (Lun,&
          '('' DDRIV3:The solution determined is not accurate enough.'' //)')
      ELSEIF( Kprint>=2 ) THEN
        WRITE (Lun,&
          '('' DDRIV3:The solution determined is not accurate enough.'')')
        WRITE (Lun,&
          '('' The values of parameters, results, and statistical quantities are:'')')
        WRITE (Lun,*) ' EPS = ', eps, ', EWT = ', ewt, ', IERROR = ', IERROR
        WRITE (Lun,*) ' MINT = ', mint, ', MITER = ', MITER, ', IMPL = ', IMPL
        WRITE (Lun,*) ' T ', t
        WRITE (Lun,*) ' Y(1) ', y(1)
        WRITE (Lun,*) ' Y(2) ', y(2)
        WRITE (Lun,*) ' Y(3) ', y(3)
        WRITE (Lun,*) ' Number of steps taken is  ', nstep
        WRITE (Lun,*) ' Number of evaluations of the right hand side is  ', nfe
        WRITE (Lun,*) ' Number of evaluations of the Jacobian matrix is  ', nje
        WRITE (Lun,'(//)')
      END IF
      Ipass = 0
    ELSEIF( Kprint==2 ) THEN
      WRITE (Lun,&
        '('' DDRIV3:The solution determined met the expected values.'' //)')
    ELSEIF( Kprint==3 ) THEN
      WRITE (Lun,&
        '('' DDRIV3:The solution determined met the expected values.'')')
      WRITE (Lun,'('' The values of results are '')')
      WRITE (Lun,*) ' T ', t
      WRITE (Lun,*) ' Y(1) ', y(1)
      WRITE (Lun,*) ' Y(2) ', y(2)
      WRITE (Lun,*) ' Y(3) ', y(3)
      WRITE (Lun,'(/)')
    END IF
    !  Run DDRIV3 with invalid input.
!    t = 0._DP
!    y(1) = 10._DP
!    y(2) = 0._DP
!    y(3) = 10._DP
!    y(4) = ALFA
!    nstate = 1
!    tout = 10._DP
!    mint = 2
!    lenw = 301
!    leniwx = 1
!    CALL DDRIV3(N,t,y,DDF,nstate,tout,NTASK,NROOT,eps,ewt,IERROR,mint,MITER,&
!      IMPL,ML,MU,MXORD,HMAX,work,lenw,iwork,leniwx,dum_JACOBN,dum_FA,nde,&
!      MXSTEP,dum_G,dum_USERS,ierflg)
!    IF( ierflg/=33 ) THEN
!      IF( Kprint==1 ) THEN
!        WRITE (Lun,&
!          '('' DDRIV3:An invalid parameter has not been correctly detected.'' //)')
!      ELSEIF( Kprint>=2 ) THEN
!        WRITE (Lun,&
!          '('' DDRIV3:An invalid parameter has not been correctly detected.'')')
!        WRITE (Lun,*) ' The value of LENIW was set to ', leniwx
!        WRITE (Lun,*) ' NSTATE = ', nstate, ', Error number = ', ierflg
!        WRITE (Lun,&
!          '('' The values of parameters, results, and statistical quantities are:'')')
!        WRITE (Lun,*) ' EPS = ', eps, ', EWT = ', ewt, ', IERROR = ', IERROR
!        WRITE (Lun,*) ' MINT = ', mint, ', MITER = ', MITER, ', IMPL = ', IMPL
!        WRITE (Lun,*) ' T ', t
!        WRITE (Lun,*) ' Y(1) ', y(1)
!        WRITE (Lun,*) ' Y(2) ', y(2)
!        WRITE (Lun,*) ' Y(3) ', y(3)
!        WRITE (Lun,*) ' Number of steps taken is  ', nstep
!        WRITE (Lun,*) ' Number of evaluations of the right hand side is  ', nfe
!        WRITE (Lun,*) ' Number of evaluations of the Jacobian matrix is  ', nje
!        WRITE (Lun,'(//)')
!      END IF
!      Ipass = 0
!    ELSEIF( Kprint==2 ) THEN
!      WRITE (Lun,&
!        '('' DDRIV3:An invalid parameter has been correctly detected.'' //)')
!    ELSEIF( Kprint==3 ) THEN
!      WRITE (Lun,&
!        '('' DDRIV3:An invalid parameter has been correctly detected.'')')
!      WRITE (Lun,*) ' The value of LENIW was set to ', leniwx
!      WRITE (Lun,*) ' NSTATE = ', nstate, ', Error number = ', ierflg
!      WRITE (Lun,'(/)')
!    END IF
!    num_xer = 0

  CONTAINS
    REAL(DP) PURE FUNCTION dum_G(N,T,Y,Iroot)
      INTEGER, INTENT(IN) :: N, Iroot
      REAL(DP), INTENT(IN) :: T
      REAL(DP), INTENT(IN) :: Y(N)
      dum_G = SUM(REAL(Y))+ T + Iroot
    END FUNCTION dum_G
    PURE SUBROUTINE dum_JACOBN(N,T,Y,Dfdy,Matdim,Ml,Mu)
      INTEGER, INTENT(IN) :: N, Matdim, Ml, Mu
      REAL(DP), INTENT(IN) :: T
      REAL(DP), INTENT(IN) :: Y(N)
      REAL(DP), INTENT(OUT) :: Dfdy(Matdim,N)
      Dfdy = T + Y(1) + Ml + Mu
    END SUBROUTINE dum_JACOBN
    PURE SUBROUTINE dum_USERS(Y,Yh,Ywt,Save1,Save2,T,H,El,Impl,N,Nde,Iflag)
      INTEGER, INTENT(IN) :: Impl, N, Nde, Iflag
      REAL(DP), INTENT(IN) :: T, H, El
      REAL(DP), INTENT(IN) :: Y(N), Yh(N,13), Ywt(N)
      REAL(DP), INTENT(INOUT) :: Save1(N), Save2(N)
      Save1 = Yh(:,1) + Y + Ywt
      Save2 = T + H + El + Impl + Nde + Iflag
    END SUBROUTINE dum_USERS
    PURE SUBROUTINE dum_FA(N,T,Y,A,Matdim,Ml,Mu,Nde)
      INTEGER, INTENT(IN) :: N, Matdim, Ml, Mu, Nde
      REAL(DP), INTENT(IN) :: T
      REAL(DP), INTENT(IN) :: Y(N)
      REAL(DP), INTENT(INOUT) :: A(:,:)
      A = Y(1) + T + Matdim + Ml + Mu + Nde
    END SUBROUTINE dum_FA
  END SUBROUTINE DDQCK
  !** DDF
  PURE SUBROUTINE DDF(N,T,Y,Yp)
    !> Quick check for SLATEC routines DDRIV1, DDRIV2 and DDRIV3.
    !***
    ! **Library:**   SLATEC (SDRIVE)
    !***
    ! **Category:**  I1A2, I1A1B
    !***
    ! **Type:**      DOUBLE PRECISION (SDF-S, DDF-D, CDF-C)
    !***
    ! **Keywords:**  DDRIV1, DDRIV2, DDRIV3, QUICK CHECK, SDRIVE
    !***
    ! **Author:**  Kahaner, D. K., (NIST)
    !             National Institute of Standards and Technology
    !             Gaithersburg, MD  20899
    !           Sutherland, C. D., (LANL)
    !             Mail Stop D466
    !             Los Alamos National Laboratory
    !             Los Alamos, NM  87545
    !***
    ! **See also:**  DDQCK
    !***
    ! **Routines called:**  (NONE)

    !* REVISION HISTORY  (YYMMDD)
    !   890405  DATE WRITTEN
    !   890405  Revised to meet SLATEC standards.

    INTEGER, INTENT(IN) :: N
    REAL(DP), INTENT(IN) :: T, Y(:)
    REAL(DP), INTENT(OUT) :: Yp(:)
    !
    REAL(DP) :: alfa
    !* FIRST EXECUTABLE STATEMENT  DDF
    alfa = Y(N+1)
    Yp(1) = 1._DP + alfa*(Y(2)-Y(1)) - Y(1)*Y(3)
    Yp(2) = alfa*(Y(1)-Y(2)) - Y(2)*Y(3)
    Yp(3) = 1._DP - Y(3)*(Y(1)+Y(2))
    !
  END SUBROUTINE DDF
  !
END MODULE TEST46_MOD
!** TEST46
PROGRAM TEST46
  USE TEST46_MOD, ONLY : DDQCK
  USE ISO_FORTRAN_ENV, ONLY : INPUT_UNIT, OUTPUT_UNIT
  USE common_mod, ONLY : GET_ARGUMENT
  IMPLICIT NONE
  !> Driver for testing SLATEC subprograms
  !            DDRIV1  DDRIV2  DDRIV3
  !***
  ! **Library:**   SLATEC
  !***
  ! **Category:**  I1A2, I1A1B
  !***
  ! **Type:**      DOUBLE PRECISION (TEST45-S, TEST46-D, TEST47-C)
  !***
  ! **Keywords:**  SDRIVE, QUICK CHECK DRIVER
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
  !        DDRIV1  DDRIV2  DDRIV3
  !
  !***
  ! **References:**  Kirby W. Fong, Thomas H. Jefferson, Tokihiko Suyehiro
  !                 and Lee Walton, Guide to the SLATEC Common Mathema-
  !                 tical Library, April 10, 1990.
  !***
  ! **Routines called:**  DDQCK, I1MACH, XERMAX, XSETF

  !* REVISION HISTORY  (YYMMDD)
  !   920801  DATE WRITTEN
  INTEGER :: ipass, kprint, lin, lun, nfail
  !* FIRST EXECUTABLE STATEMENT  TEST46
  lun = OUTPUT_UNIT
  lin = INPUT_UNIT
  nfail = 0
  !
  !     Read KPRINT parameter
  !
  CALL GET_ARGUMENT(kprint)
  !
  !     Test double precision SDRIVE
  !
  CALL DDQCK(lun,kprint,ipass)
  IF( ipass==0 ) nfail = nfail + 1
  !
  !     Write PASS or FAIL message
  !
  IF( nfail==0 ) THEN
    WRITE (lun,99001)
    99001 FORMAT (/' --------------TEST46 PASSED ALL TESTS----------------')
  ELSE
    WRITE (lun,99002) nfail
    99002 FORMAT (/' ************* WARNING -- ',I5,&
      ' TEST(S) FAILED IN PROGRAM TEST46 *************')
  END IF
  STOP
END PROGRAM TEST46
