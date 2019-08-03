MODULE TEST45_MOD
  USE service, ONLY : SP, DP
  IMPLICIT NONE

CONTAINS
  !** SDQCK
  SUBROUTINE SDQCK(Lun,Kprint,Ipass)
    !> Quick check for SLATEC routines SDRIV1, SDRIV2 and SDRIV3.
    !***
    ! **Library:**   SLATEC (SDRIVE)
    !***
    ! **Category:**  I1A2, I1A1B
    !***
    ! **Type:**      SINGLE PRECISION (SDQCK-S, DDQCK-D, CDQCK-C)
    !***
    ! **Keywords:**  QUICK CHECK, SDRIV1, SDRIV2, SDRIV3
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
    ! **Routines called:**  R1MACH, SDF, SDRIV1, SDRIV2, SDRIV3, XERCLR

    !* REVISION HISTORY  (YYMMDD)
    !   890405  DATE WRITTEN
    !   890405  Revised to meet SLATEC standards.
    USE service, ONLY : eps_sp
    USE diff_integ_eq, ONLY : SDRIV1, SDRIV2, SDRIV3
    !
    REAL(SP) :: eps, t, tout
    INTEGER :: ierflg, Ipass, Kprint, leniw, lenw, Lun, mint, mstate, nde, nfe, &
      nje, nstate, nstep
    REAL(SP), PARAMETER :: ALFA = 1._SP, HMAX = 15._SP
    INTEGER, PARAMETER :: IERROR = 3, IMPL = 0, LENWMX = 342, LIWMX = 53, &
      MITER = 5, ML = 2, MU = 2, MXORD = 5, MXSTEP = 1000, N = 3, NROOT = 0, NTASK = 1
    REAL(SP) :: work(LENWMX), y(N+1)
    INTEGER :: iwork(LIWMX)
    REAL(SP), PARAMETER :: ewt(1) = .00001_SP
    !* FIRST EXECUTABLE STATEMENT  SDQCK
    eps = eps_sp**(1._SP/3._SP)
    Ipass = 1
    !  Exercise SDRIV1 for problem with known solution.
    y(4) = ALFA
    t = 0._SP
    y(1) = 10._SP
    y(2) = 0._SP
    y(3) = 10._SP
    tout = 10._SP
    mstate = 1
    lenw = 342
    CALL SDRIV1(N,t,y,SDF,tout,mstate,eps,work,lenw,ierflg)
    nstep = INT( work(lenw-(N+50)+3) )
    nfe = INT( work(lenw-(N+50)+4) )
    nje = INT( work(lenw-(N+50)+5) )
    IF( mstate/=2 ) THEN
      IF( Kprint==1 ) THEN
        WRITE (Lun,&
          '('' While using SDRIV1, a solution was not obtained.'' //)' )
      ELSEIF( Kprint>=2 ) THEN
        WRITE (Lun,&
          '('' While using SDRIV1, a solution was not obtained.'')' )
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
    ELSEIF( ABS(1._SP-y(1)*1.5_SP)>eps**(2._SP/3._SP) .OR. ABS(1._SP-y(2)*3._SP)&
        >eps**(2._SP/3._SP) .OR. ABS(1._SP-y(3))>eps**(2._SP/3._SP) ) THEN
      IF( Kprint==1 ) THEN
        WRITE (Lun,&
          '('' SDRIV1:The solution determined is not accurate enough.'' //)')
      ELSEIF( Kprint==2 ) THEN
        WRITE (Lun,&
          '('' SDRIV1:The solution determined is not accurate enough.'')')
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
        '('' SDRIV1:The solution determined met the expected values.'' //)')
    ELSEIF( Kprint==3 ) THEN
      WRITE (Lun,&
        '('' SDRIV1:The solution determined met the expected values.'')')
      WRITE (Lun,'('' The values of results are '')')
      WRITE (Lun,*) ' T ', t
      WRITE (Lun,*) ' Y(1) ', y(1)
      WRITE (Lun,*) ' Y(2) ', y(2)
      WRITE (Lun,*) ' Y(3) ', y(3)
      WRITE (Lun,'(/)')
    END IF
    !  Run SDRIV1 with invalid input.
!    nx = 201
!    t = 0._SP
!    y(1) = 10._SP
!    y(2) = 0._SP
!    y(3) = 10._SP
!    y(4) = ALFA
!    tout = 10._SP
!    mstate = 1
!    lenw = 342
!    CALL SDRIV1(nx,t,y,SDF,tout,mstate,eps,work,lenw,ierflg)
!    IF( ierflg/=21 ) THEN
!      IF( Kprint==1 ) THEN
!        WRITE (Lun,&
!          '('' SDRIV1:An invalid parameter has not been correctly detected.'' //)')
!      ELSEIF( Kprint>=2 ) THEN
!        WRITE (Lun,&
!          '('' SDRIV1:An invalid parameter has not been correctly detected.'')')
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
!        '('' SDRIV1:An invalid parameter has been correctly detected.'' //)')
!    ELSEIF( Kprint==3 ) THEN
!      WRITE (Lun,&
!        '('' SDRIV1:An invalid parameter has been correctly detected.'')')
!      WRITE (Lun,*) ' The value of N was set to ', nx
!      WRITE (Lun,*) ' MSTATE = ', mstate, ', Error number = ', ierflg
!      WRITE (Lun,'(/)')
!    END IF
!    num_xer = 0
    !  Exercise SDRIV2 for problem with known solution.
    t = 0._SP
    y(1) = 10._SP
    y(2) = 0._SP
    y(3) = 10._SP
    y(4) = ALFA
    mstate = 1
    tout = 10._SP
    mint = 1
    lenw = 298
    leniw = 50
    CALL SDRIV2(N,t,y,SDF,tout,mstate,NROOT,eps,ewt(1),mint,work,lenw,iwork,&
      leniw,dum_G,ierflg)
    nstep = iwork(3)
    nfe = iwork(4)
    nje = iwork(5)
    IF( mstate/=2 ) THEN
      IF( Kprint==1 ) THEN
        WRITE (Lun,&
          '('' While using SDRIV2, a solution was not obtained.'' //)' )
      ELSEIF( Kprint>=2 ) THEN
        WRITE (Lun,&
          '('' While using SDRIV2, a solution was not obtained.'')' )
        WRITE (Lun,*) ' MSTATE = ', mstate, ', Error number = ', ierflg
        WRITE (Lun,&
          '('' The values of parameters, results, and statistical quantities are:'')')
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
    ELSEIF( ABS(1._SP-y(1)*1.5_SP)>eps**(2._SP/3._SP) .OR. ABS(1._SP-y(2)*3._SP)&
        >eps**(2._SP/3._SP) .OR. ABS(1._SP-y(3))>eps**(2._SP/3._SP) ) THEN
      IF( Kprint==1 ) THEN
        WRITE (Lun,&
          '('' SDRIV2:The solution determined is not accurate enough.'' //)')
      ELSEIF( Kprint==2 ) THEN
        WRITE (Lun,&
          '('' SDRIV2:The solution determined is not accurate enough.'')')
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
        '('' SDRIV2:The solution determined met the expected values.'' //)')
    ELSEIF( Kprint==3 ) THEN
      WRITE (Lun,&
        '('' SDRIV2:The solution determined met the expected values.'')')
      WRITE (Lun,'('' The values of results are '')')
      WRITE (Lun,*) ' T ', t
      WRITE (Lun,*) ' Y(1) ', y(1)
      WRITE (Lun,*) ' Y(2) ', y(2)
      WRITE (Lun,*) ' Y(3) ', y(3)
      WRITE (Lun,'(/)')
    END IF
    !  Run SDRIV2 with invalid input.
!    t = 0._SP
!    y(1) = 10._SP
!    y(2) = 0._SP
!    y(3) = 10._SP
!    y(4) = ALFA
!    tout = 10._SP
!    mstate = 1
!    mint = 1
!    lenwx = 1
!    leniw = 50
!    CALL SDRIV2(N,t,y,SDF,tout,mstate,NROOT,eps,ewt(1),mint,work,lenwx,iwork,&
!      leniw,dum_G,ierflg)
!    IF( ierflg/=32 ) THEN
!      IF( Kprint==1 ) THEN
!        WRITE (Lun,&
!          '('' SDRIV2:An invalid parameter has not been correctly detected.'' //)')
!      ELSEIF( Kprint>=2 ) THEN
!        WRITE (Lun,&
!          '('' SDRIV2:An invalid parameter has not been correctly detected.'')')
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
!        '('' SDRIV2:An invalid parameter has been correctly detected.'' //)')
!    ELSEIF( Kprint==3 ) THEN
!      WRITE (Lun,&
!        '('' SDRIV2:An invalid parameter has been correctly detected.'')')
!      WRITE (Lun,*) ' The value of LENW was set to ', lenwx
!      WRITE (Lun,*) ' MSTATE = ', mstate, ', Error number = ', ierflg
!      WRITE (Lun,'(/)')
!    END IF
!    num_xer = 0
    !  Exercise SDRIV3 for problem with known solution.
    t = 0._SP
    y(1) = 10._SP
    y(2) = 0._SP
    y(3) = 10._SP
    y(4) = ALFA
    nstate = 1
    tout = 10._SP
    mint = 2
    lenw = 301
    leniw = 53
    CALL SDRIV3(N,t,y,SDF,nstate,tout,NTASK,NROOT,eps,ewt,IERROR,mint,MITER,&
      IMPL,ML,MU,MXORD,HMAX,work,lenw,iwork,leniw,dum_JACOBN,dum_FA,nde,&
      MXSTEP,dum_G,dum_USERS,ierflg)
    nstep = iwork(3)
    nfe = iwork(4)
    nje = iwork(5)
    IF( nstate/=2 ) THEN
      IF( Kprint==1 ) THEN
        WRITE (Lun,&
          '('' While using SDRIV3, a solution was not obtained.'' //)' )
      ELSEIF( Kprint>=2 ) THEN
        WRITE (Lun,&
          '('' While using SDRIV3, a solution was not obtained.'')' )
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
    ELSEIF( ABS(1._SP-y(1)*1.5_SP)>eps**(2._SP/3._SP) .OR. ABS(1._SP-y(2)*3._SP)&
        >eps**(2._SP/3._SP) .OR. ABS(1._SP-y(3))>eps**(2._SP/3._SP) ) THEN
      IF( Kprint==1 ) THEN
        WRITE (Lun,&
          '('' SDRIV3:The solution determined is not accurate enough.''//)')
      ELSEIF( Kprint>=2 ) THEN
        WRITE (Lun,&
          '('' SDRIV3:The solution determined is not accurate enough.'')')
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
        '('' SDRIV3:The solution determined met the expected values.'' //)')
    ELSEIF( Kprint==3 ) THEN
      WRITE (Lun,&
        '('' SDRIV3:The solution determined met the expected values.'')')
      WRITE (Lun,'('' The values of results are '')')
      WRITE (Lun,*) ' T ', t
      WRITE (Lun,*) ' Y(1) ', y(1)
      WRITE (Lun,*) ' Y(2) ', y(2)
      WRITE (Lun,*) ' Y(3) ', y(3)
      WRITE (Lun,'(/)')
    END IF
    !  Run SDRIV3 with invalid input.
!    t = 0._SP
!    y(1) = 10._SP
!    y(2) = 0._SP
!    y(3) = 10._SP
!    y(4) = ALFA
!    nstate = 1
!    tout = 10._SP
!    mint = 2
!    lenw = 301
!    leniwx = 1
!    CALL SDRIV3(N,t,y,SDF,nstate,tout,NTASK,NROOT,eps,ewt,IERROR,mint,MITER,&
!      IMPL,ML,MU,MXORD,HMAX,work,lenw,iwork,leniwx,dum_JACOBN,dum_FA,nde,&
!      MXSTEP,dum_G,dum_USERS,ierflg)
!    IF( ierflg/=33 ) THEN
!      IF( Kprint==1 ) THEN
!        WRITE (Lun,&
!          '('' SDRIV3:An invalid parameter has not been correctly detected.'' //)')
!      ELSEIF( Kprint>=2 ) THEN
!        WRITE (Lun,&
!          '('' SDRIV3:An invalid parameter has not been correctly detected.'')')
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
!      WRITE (Lun,'('' SDRIV3:An invalid parameter has been correctly detected.'' //)')
!    ELSEIF( Kprint==3 ) THEN
!      WRITE (Lun,'('' SDRIV3:An invalid parameter has been correctly detected.'')')
!      WRITE (Lun,*) ' The value of LENIW was set to ', leniwx
!      WRITE (Lun,*) ' NSTATE = ', nstate, ', Error number = ', ierflg
!      WRITE (Lun,'(/)')
!    END IF
!    num_xer = 0

  CONTAINS
    REAL(SP) PURE FUNCTION dum_G(N,T,Y,Iroot)
      INTEGER, INTENT(IN) :: N, Iroot
      REAL(SP), INTENT(IN) :: T
      REAL(SP), INTENT(IN) :: Y(N)
      dum_G = SUM(REAL(Y))+ T + Iroot
    END FUNCTION dum_G
    PURE SUBROUTINE dum_JACOBN(N,T,Y,Dfdy,Matdim,Ml,Mu)
      INTEGER, INTENT(IN) :: N, Matdim, Ml, Mu
      REAL(SP), INTENT(IN) :: T
      REAL(SP), INTENT(IN) :: Y(N)
      REAL(SP), INTENT(OUT) :: Dfdy(Matdim,N)
      Dfdy = T + Y(1) + Ml + Mu
    END SUBROUTINE dum_JACOBN
    PURE SUBROUTINE dum_USERS(Y,Yh,Ywt,Save1,Save2,T,H,El,Impl,N,Nde,Iflag)
      INTEGER, INTENT(IN) :: Impl, N, Nde, Iflag
      REAL(SP), INTENT(IN) :: T, H, El
      REAL(SP), INTENT(IN) :: Y(N), Yh(N,13), Ywt(N)
      REAL(SP), INTENT(INOUT) :: Save1(N), Save2(N)
      Save1 = Yh(:,1) + Y + Ywt
      Save2 = T + H + El + Impl + Nde + Iflag
    END SUBROUTINE dum_USERS
    PURE SUBROUTINE dum_FA(N,T,Y,A,Matdim,Ml,Mu,Nde)
      INTEGER, INTENT(IN) :: N, Matdim, Ml, Mu, Nde
      REAL(SP), INTENT(IN) :: T
      REAL(SP), INTENT(IN) :: Y(N)
      REAL(SP), INTENT(INOUT) :: A(:,:)
      A = Y(1) + T + Matdim + Ml + Mu + Nde
    END SUBROUTINE dum_FA
  END SUBROUTINE SDQCK
  !** SDF
  PURE SUBROUTINE SDF(N,T,Y,Yp)
    !> Quick check for SLATEC routines SDRIV1, SDRIV2 and SDRIV3.
    !***
    ! **Library:**   SLATEC (SDRIVE)
    !***
    ! **Category:**  I1A2, I1A1B
    !***
    ! **Type:**      SINGLE PRECISION (SDF-S, DDF-D, CDF-C)
    !***
    ! **Keywords:**  QUICK CHECK, SDRIV1, SDRIV2, SDRIV3
    !***
    ! **Author:**  Kahaner, D. K., (NIST)
    !             National Institute of Standards and Technology
    !             Gaithersburg, MD  20899
    !           Sutherland, C. D., (LANL)
    !             Mail Stop D466
    !             Los Alamos National Laboratory
    !             Los Alamos, NM  87545
    !***
    ! **See also:**  SDQCK
    !***
    ! **Routines called:**  (NONE)

    !* REVISION HISTORY  (YYMMDD)
    !   890405  DATE WRITTEN
    !   890405  Revised to meet SLATEC standards.

    INTEGER, INTENT(IN) :: N
    REAL(SP), INTENT(IN) :: T, Y(:)
    REAL(SP), INTENT(OUT) :: Yp(:)
    !
    REAL(SP) :: alfa
    !* FIRST EXECUTABLE STATEMENT  SDF
    alfa = Y(N+1)
    Yp(1) = 1._SP + alfa*(Y(2)-Y(1)) - Y(1)*Y(3)
    Yp(2) = alfa*(Y(1)-Y(2)) - Y(2)*Y(3)
    Yp(3) = 1._SP - Y(3)*(Y(1)+Y(2))
    !
  END SUBROUTINE SDF
  !
END MODULE TEST45_MOD
!** TEST45
PROGRAM TEST45
  USE TEST45_MOD, ONLY : SDQCK
  USE ISO_FORTRAN_ENV, ONLY : INPUT_UNIT, OUTPUT_UNIT
  USE common_mod, ONLY : GET_ARGUMENT
  IMPLICIT NONE
  !> Driver for testing SLATEC subprograms
  !            SDRIV1  SDRIV2  SDRIV3
  !***
  ! **Library:**   SLATEC
  !***
  ! **Category:**  I1A2, I1A1B
  !***
  ! **Type:**      SINGLE PRECISION (TEST45-S, TEST46-D, TEST47-C)
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
  !        SDRIV1  SDRIV2  SDRIV3
  !
  !***
  ! **References:**  Kirby W. Fong, Thomas H. Jefferson, Tokihiko Suyehiro
  !                 and Lee Walton, Guide to the SLATEC Common Mathema-
  !                 tical Library, April 10, 1990.
  !***
  ! **Routines called:**  I1MACH, SDQCK, XERMAX, XSETF

  !* REVISION HISTORY  (YYMMDD)
  !   920801  DATE WRITTEN
  INTEGER :: ipass, kprint, lin, lun, nfail
  !* FIRST EXECUTABLE STATEMENT  TEST45
  lun = OUTPUT_UNIT
  lin = INPUT_UNIT
  nfail = 0
  !
  !     Read KPRINT parameter
  !
  CALL GET_ARGUMENT(kprint)
  !
  !     Test single precision SDRIVE
  !
  CALL SDQCK(lun,kprint,ipass)
  IF( ipass==0 ) nfail = nfail + 1
  !
  !     Write PASS or FAIL message
  !
  IF( nfail==0 ) THEN
    WRITE (lun,99001)
    99001 FORMAT (/' --------------TEST45 PASSED ALL TESTS----------------')
  ELSE
    WRITE (lun,99002) nfail
    99002 FORMAT (/' ************* WARNING -- ',I5,&
      ' TEST(S) FAILED IN PROGRAM TEST45 *************')
  END IF
  STOP
END PROGRAM TEST45
