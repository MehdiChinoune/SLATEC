!** BVSUP
SUBROUTINE BVSUP(Y,Nrowy,Ncomp,Xpts,Nxpts,A,Nrowa,Alpha,Nic,B,Nrowb,Beta,&
    Nfc,Igofx,Re,Ae,Iflag,Work,Ndw,Iwork,Ndiw,Neqivp)
  !> Solve a linear two-point boundary value problem using superposition coupled
  !  with an orthonormalization procedure and a variable-step integration scheme.
  !***
  ! **Library:**   SLATEC
  !***
  ! **Category:**  I1B1
  !***
  ! **Type:**      SINGLE PRECISION (BVSUP-S, DBVSUP-D)
  !***
  ! **Keywords:**  ORTHONORMALIZATION, SHOOTING,
  !             TWO-POINT BOUNDARY VALUE PROBLEM
  !***
  ! **Author:**  Scott, M. R., (SNLA)
  !           Watts, H. A., (SNLA)
  !***
  ! **Description:**
  !
  !- *********************************************************************
  !     Subroutine BVSUP solves a LINEAR two-point boundary-value problem
  !     of the form
  !                        dY/dX = MATRIX(X,U)*Y(X) + G(X,U)
  !                A*Y(Xinitial) = ALPHA,  B*Y(Xfinal) = BETA
  !
  !     Coupled with the solution of the initial value problem
  !
  !                        dU/dX = F(X,U)
  !                      U(Xinitial) = ETA
  !
  !- *********************************************************************
  !     Abstract
  !        The method of solution uses superposition coupled with an
  !     orthonormalization procedure and a variable-step integration
  !     scheme.  Each time the superposition solutions start to
  !     lose their numerical linear independence, the vectors are
  !     reorthonormalized before integration proceeds.  The underlying
  !     principle of the algorithm is then to piece together the
  !     intermediate (orthogonalized) solutions, defined on the various
  !     subintervals, to obtain the desired solutions.
  !
  !- *********************************************************************
  !     INPUT to BVSUP
  !- *********************************************************************
  !
  !     NROWY = Actual row dimension of Y in calling program.
  !             NROWY must be >= NCOMP
  !
  !     NCOMP = Number of components per solution vector.
  !             NCOMP is equal to number of original differential
  !             equations.  NCOMP = NIC + NFC.
  !
  !     XPTS = Desired output points for solution. They must be monotonic.
  !            Xinitial = XPTS(1)
  !            Xfinal = XPTS(NXPTS)
  !
  !     NXPTS = Number of output points
  !
  !     A(NROWA,NCOMP) = Boundary condition matrix at Xinitial,
  !                      must be contained in (NIC,NCOMP) sub-matrix.
  !
  !     NROWA = Actual row dimension of A in calling program,
  !             NROWA must be >= NIC.
  !
  !     ALPHA(NIC+NEQIVP) = Boundary conditions at Xinitial.
  !                         If NEQIVP > 0 (see below), the boundary
  !                         conditions at Xinitial for the initial value
  !                         equations must be stored starting in
  !                         position (NIC + 1) of ALPHA.
  !                         Thus,  ALPHA(NIC+K) = ETA(K).
  !
  !     NIC = Number of boundary conditions at Xinitial.
  !
  !     B(NROWB,NCOMP) = Boundary condition matrix at Xfinal,
  !                      must be contained in (NFC,NCOMP) sub-matrix.
  !
  !     NROWB = Actual row dimension of B in calling program,
  !             NROWB must be >= NFC.
  !
  !     BETA(NFC) = Boundary conditions at Xfinal.
  !
  !     NFC = Number of boundary conditions at Xfinal
  !
  !     IGOFX =0 -- The inhomogeneous term G(X) is identically zero.
  !           =1 -- The inhomogeneous term G(X) is not identically zero.
  !                 (if IGOFX=1, then subroutine GVEC (or UVEC) must be
  !                  supplied).
  !
  !     RE = Relative error tolerance used by the integrator
  !          (see one of the integrators)
  !
  !     AE = Absolute error tolerance used by the integrator
  !          (see one of the integrators)
  !- *NOTE-  RE and AE should not both be zero.
  !
  !     IFLAG = A status parameter used principally for output.
  !             However, for efficient solution of problems which
  !             are originally defined as complex valued (but
  !             converted to real systems to use this code), the
  !             user must set IFLAG=13 on input. See the comment below
  !             for more information on solving such problems.
  !
  !     WORK(NDW) = Floating point array used for internal storage.
  !
  !     NDW = Actual dimension of WORK array allocated by user.
  !           An estimate for NDW can be computed from the following
  !            NDW = 130 + NCOMP**2 * (6 + NXPTS/2 + expected number of
  !                                                orthonormalizations/8)
  !             For the DISK or TAPE storage mode,
  !            NDW = 6 * NCOMP**2 + 10 * NCOMP + 130
  !  However, when the ADAMS integrator is to be used, the estimates are
  !            NDW = 130 + NCOMP**2 * (13 + NXPTS/2 + expected number of
  !                                                orthonormalizations/8)
  !    and     NDW = 13 * NCOMP**2 + 22 * NCOMP + 130  , respectively.
  !
  !     IWORK(NDIW) = Integer array used for internal storage.
  !
  !     NDIW = Actual dimension of IWORK array allocated by user.
  !            An estimate for NDIW can be computed from the following
  !            NDIW = 68 + NCOMP * (1 + expected number of
  !                                        orthonormalizations)
  !- *NOTE --  The amount of storage required is problem dependent and may
  !            be difficult to predict in advance. Experience has shown
  !            that for most problems 20 or fewer orthonormalizations
  !            should suffice. If the problem cannot be completed with the
  !            allotted storage, then a message will be printed which
  !            estimates the amount of storage necessary. In any case, the
  !            user can examine the IWORK array for the actual storage
  !            requirements, as described in the output information below.
  !
  !     NEQIVP = Number of auxiliary initial value equations being added
  !              to the boundary value problem.
  !- *NOTE -- Occasionally the coefficients  MATRIX  and/or  G  may be
  !           functions which depend on the independent variable  X  and
  !           on  U, the solution of an auxiliary initial value problem.
  !           In order to avoid the difficulties associated with
  !           interpolation, the auxiliary equations may be solved
  !           simultaneously with the given boundary value problem.
  !           This initial value problem may be LINEAR or NONLINEAR.
  !                 See SAND77-1328 for an example.
  !
  !
  !     The user must supply subroutines FMAT, GVEC, UIVP and UVEC, when
  !     needed (they MUST be so named), to evaluate the derivatives
  !     as follows
  !
  !        A. FMAT must be supplied.
  !
  !              SUBROUTINE FMAT(X,Y,YP)
  !              X = Independent variable (input to FMAT)
  !              Y = Dependent variable vector (input to FMAT)
  !              YP = dY/dX = Derivative vector (output from FMAT)
  !
  !            Compute the derivatives for the HOMOGENEOUS problem
  !              YP(I) = dY(I)/dX = MATRIX(X) * Y(I) , I = 1,...,NCOMP
  !
  !            When (NEQIVP > 0) and  MATRIX  is dependent on  U  as
  !            well as on  X, the following common statement must be
  !            included in FMAT
  !                    COMMON /MLIVP/ NOFST
  !            For convenience, the  U  vector is stored at the bottom
  !            of the  Y  array.  Thus, during any call to FMAT,
  !            U(I) is referenced by  Y(NOFST + I).
  !
  !
  !            Subroutine BVDER calls FMAT NFC times to evaluate the
  !            homogeneous equations and, if necessary, it calls FMAT once
  !            in evaluating the particular solution. Since X remains
  !            unchanged in this sequence of calls it is possible to
  !            realize considerable computational savings for complicated
  !            and expensive evaluations of the MATRIX entries. To do this
  !            the user merely passes a variable, say XS, via COMMON where
  !            XS is defined in the main program to be any value except
  !            the initial X. Then the non-constant elements of MATRIX(X)
  !            appearing in the differential equations need only be
  !            computed if X is unequal to XS, whereupon XS is reset to X.
  !
  !
  !        B. If  NEQIVP > 0,  UIVP must also be supplied.
  !
  !              SUBROUTINE UIVP(X,U,UP)
  !              X = Independent variable (input to UIVP)
  !              U = Dependent variable vector (input to UIVP)
  !              UP = dU/dX = Derivative vector (output from UIVP)
  !
  !            Compute the derivatives for the auxiliary initial value eqs
  !              UP(I) = dU(I)/dX, I = 1,...,NEQIVP.
  !
  !            Subroutine BVDER calls UIVP once to evaluate the
  !            derivatives for the auxiliary initial value equations.
  !
  !
  !        C. If  NEQIVP = 0  and  IGOFX = 1,  GVEC must be supplied.
  !
  !              SUBROUTINE GVEC(X,G)
  !              X = Independent variable (input to GVEC)
  !              G = Vector of inhomogeneous terms G(X) (output from GVEC)
  !
  !            Compute the inhomogeneous terms G(X)
  !                G(I) = G(X) values for I = 1,...,NCOMP.
  !
  !            Subroutine BVDER calls GVEC in evaluating the particular
  !            solution provided G(X) is NOT identically zero. Thus, when
  !            IGOFX=0, the user need NOT write a GVEC subroutine. Also,
  !            the user does not have to bother with the computational
  !            savings scheme for GVEC as this is automatically achieved
  !            via the BVDER subroutine.
  !
  !
  !        D. If  NEQIVP > 0  and  IGOFX = 1,  UVEC must be supplied.
  !
  !              SUBROUTINE UVEC(X,U,G)
  !              X = Independent variable (input to UVEC)
  !              U = Dependent variable vector from the auxiliary initial
  !                  value problem    (input to UVEC)
  !              G = Array of inhomogeneous terms G(X,U)(output from UVEC)
  !
  !            Compute the inhomogeneous terms G(X,U)
  !                G(I) = G(X,U) values for I = 1,...,NCOMP.
  !
  !            Subroutine BVDER calls UVEC in evaluating the particular
  !            solution provided G(X,U) is NOT identically zero.  Thus,
  !            when IGOFX=0, the user need NOT write a UVEC subroutine.
  !
  !
  !
  !     The following is optional input to BVSUP to give the user more
  !     flexibility in use of the code.  See SAND75-0198, SAND77-1328 ,
  !     SAND77-1690,SAND78-0522, and SAND78-1501 for more information.
  !
  !- ***CAUTION -- The user MUST zero out IWORK(1),...,IWORK(15)
  !                prior to calling BVSUP. These locations define optional
  !                input and MUST be zero UNLESS set to special values by
  !                the user as described below.
  !
  !     IWORK(1) -- Number of orthonormalization points.
  !                 A value need be set only if IWORK(11) = 1
  !
  !     IWORK(9) -- Integrator and orthonormalization parameter
  !                 (default value is 1)
  !                 1 = RUNGE-KUTTA-FEHLBERG code using GRAM-SCHMIDT test.
  !                 2 = ADAMS code using GRAM-SCHMIDT TEST.
  !
  !     IWORK(11) -- Orthonormalization points parameter
  !                  (default value is 0)
  !                  0 - Orthonormalization points not pre-assigned.
  !                  1 - Orthonormalization points pre-assigned in
  !                      the first IWORK(1) positions of WORK.
  !
  !     IWORK(12) -- Storage parameter
  !                  (default value is 0)
  !                  0 - All storage IN CORE
  !                LUN - Homogeneous and inhomogeneous solutions at
  !                     output points and orthonormalization information
  !                     are stored on DISK.  The logical unit number to be
  !                     used for DISK I/O (NTAPE) is set to IWORK(12).
  !
  !     WORK(1),... -- Pre-assigned orthonormalization points, stored
  !                    monotonically, corresponding to the direction
  !                    of integration.
  !
  !
  !
  !                 ******************************
  !                 *** COMPLEX VALUED PROBLEM ***
  !                 ******************************
  !- *NOTE***
  !       Suppose the original boundary value problem is NC equations
  !     of the form
  !                   dW/dX = MAT(X,U)*W(X) + H(X,U)
  !                 R*W(Xinitial)=GAMMA, S*W(Xfinal)=DELTA
  !
  !     where all variables are complex valued. The BVSUP code can be
  !     used by converting to a real system of size 2*NC. To solve the
  !     larger dimensioned problem efficiently,  the user must initialize
  !     IFLAG=13 on input and order the vector components according to
  !     Y(1)=real(W(1)),...,Y(NC)=real(W(NC)),Y(NC+1)=imag(W(1)),....,
  !     Y(2*NC)=imag(W(NC)). Then define
  !                        ...........................
  !                        . real(MAT)    -imag(MAT) .
  !            MATRIX  =   .                         .
  !                        . imag(MAT)     real(MAT) .
  !                        ...........................
  !
  !     The matrices A,B and vectors G,ALPHA,BETA must be defined
  !     similarly. Further details can be found in SAND78-1501.
  !
  !
  !- *********************************************************************
  !     OUTPUT from BVSUP
  !- *********************************************************************
  !
  !     Y(NROWY,NXPTS) = Solution at specified output points.
  !
  !     IFLAG output values
  !            =-5 Algorithm ,for obtaining starting vectors for the
  !                special complex problem structure, was unable to obtain
  !                the initial vectors satisfying the necessary
  !                independence criteria.
  !            =-4 Rank of boundary condition matrix A is less than NIC,
  !                as determined by LSSUDS.
  !            =-2 Invalid input parameters.
  !            =-1 Insufficient number of storage locations allocated for
  !                WORK or IWORK.
  !
  !            =0 Indicates successful solution
  !
  !            =1 A computed solution is returned but UNIQUENESS of the
  !               solution of the boundary-value problem is questionable.
  !               For an eigenvalue problem, this should be treated as a
  !               successful execution since this is the expected mode
  !               of return.
  !            =2 A computed solution is returned but the EXISTENCE of the
  !               solution to the boundary-value problem is questionable.
  !            =3 A nontrivial solution approximation is returned although
  !               the boundary condition matrix B*Y(Xfinal) is found to be
  !               nonsingular (to the desired accuracy level) while the
  !               right hand side vector is zero. To eliminate this type
  !               of return, the accuracy of the eigenvalue parameter
  !               must be improved.
  !           ***NOTE- We attempt to diagnose the correct problem behavior
  !               and report possible difficulties by the appropriate
  !               error flag.  However, the user should probably resolve
  !               the problem using smaller error tolerances and/or
  !               perturbations in the boundary conditions or other
  !               parameters. This will often reveal the correct
  !               interpretation for the problem posed.
  !
  !            =13 Maximum number of orthonormalizations attained before
  !                reaching Xfinal.
  !            =20-flag from integrator (DERKF or DEABM) values can range
  !                from 21 to 25.
  !            =30 Solution vectors form a dependent set.
  !
  !     WORK(1),...,WORK(IWORK(1)) = Orthonormalization points
  !                                  determined by BVPOR.
  !
  !     IWORK(1) = Number of orthonormalizations performed by BVPOR.
  !
  !     IWORK(2) = Maximum number of orthonormalizations allowed as
  !                calculated from storage allocated by user.
  !
  !     IWORK(3),IWORK(4),IWORK(5),IWORK(6)   Give information about
  !                actual storage requirements for WORK and IWORK
  !                arrays.  In particular,
  !                       required storage for  WORK array is
  !        IWORK(3) + IWORK(4)*(expected number of orthonormalizations)
  !
  !                       required storage for IWORK array is
  !        IWORK(5) + IWORK(6)*(expected number of orthonormalizations)
  !
  !     IWORK(8) = Final value of exponent parameter used in tolerance
  !                test for orthonormalization.
  !
  !     IWORK(16) = Number of independent vectors returned from MGSBV.
  !                 It is only of interest when IFLAG=30 is obtained.
  !
  !     IWORK(17) = Numerically estimated rank of the boundary
  !                 condition matrix defined from B*Y(Xfinal)
  !
  !- *********************************************************************
  !
  !     Necessary machine constants are defined in the function
  !     routine R1MACH. The user must make sure that the values
  !     set in R1MACH are relevant to the computer being used.
  !
  !- *********************************************************************
  !
  !***
  ! **References:**  M. R. Scott and H. A. Watts, SUPORT - a computer code
  !                 for two-point boundary-value problems via
  !                 orthonormalization, SIAM Journal of Numerical
  !                 Analysis 14, (1977), pp. 40-70.
  !               B. L. Darlow, M. R. Scott and H. A. Watts, Modifications
  !                 of SUPORT, a linear boundary value problem solver
  !                 Part I - pre-assigning orthonormalization points,
  !                 auxiliary initial value problem, disk or tape storage,
  !                 Report SAND77-1328, Sandia Laboratories, Albuquerque,
  !                 New Mexico, 1977.
  !               B. L. Darlow, M. R. Scott and H. A. Watts, Modifications
  !                 of SUPORT, a linear boundary value problem solver
  !                 Part II - inclusion of an Adams integrator, Report
  !                 SAND77-1690, Sandia Laboratories, Albuquerque,
  !                 New Mexico, 1977.
  !               M. E. Lord and H. A. Watts, Modifications of SUPORT,
  !                 a linear boundary value problem solver Part III -
  !                 orthonormalization improvements, Report SAND78-0522,
  !                 Sandia Laboratories, Albuquerque, New Mexico, 1978.
  !               H. A. Watts, M. R. Scott and M. E. Lord, Computational
  !                 solution of complex*16 valued boundary problems,
  !                 Report SAND78-1501, Sandia Laboratories,
  !                 Albuquerque, New Mexico, 1978.
  !***
  ! **Routines called:**  EXBVP, MACON, XERMSG
  !***
  ! COMMON BLOCKS    ML15TO, ML17BW, ML18JR, ML5MCO, ML8SZ

  !* REVISION HISTORY  (YYMMDD)
  !   750601  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890831  Modified array declarations.  (WRB)
  !   890921  Realigned order of variables in certain COMMON blocks.  (WRB)
  !   890921  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900510  Convert XERRWV calls to XERMSG calls.  (RWC)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  USE ML, ONLY : re_com , ae_com, indpvt_com, inhomo_com, integ_com, k1_com, &
    k2_com, k3_com, k4_com, k5_com, k6_com, k7_com, k8_com, k9_com, k10_com, &
    k11_com, kkkzpw_com, l1_com, l2_com, lllint_com, mnswot_com, icoco_com, &
    igofx_com, ndisk_com, mxnon_com, needw_com, nfcc_com, ncomp_com, neediw_com, &
    neqivp_com, nps_com, ntp_com, xbeg_com, xsav_com, xend_com, nfc_com, nic_com, &
    nopg_com, nxpts_com, ntape_com, kkkint_com
  !
  INTEGER, INTENT(IN) :: Igofx, Ncomp, Ndiw, Ndw, Neqivp, Nic, Nrowa, Nrowb, &
    Nrowy, Nxpts
  INTEGER, INTENT(INOUT) :: Iflag, Iwork(Ndiw), Nfc
  REAL(SP), INTENT(IN) :: Ae, Re
  REAL(SP), INTENT(IN) :: A(Nrowa,Ncomp), Alpha(:), B(Nrowb,Ncomp), Beta(Nxpts), Xpts(Nxpts)
  REAL(SP), INTENT(INOUT) :: Work(Ndw)
  REAL(SP), INTENT(OUT) :: Y(Nrowy,Nxpts)
  !
  INTEGER :: nitemp, non, nrtemp, is, j, k, kkkcoe, kkkcof, kkkg, kkks, kkksto, &
    kkksud, kkksvc, kkku, kkkv, kkkws, kkkyhp, kpts, lllcof, lllip, llliws, &
    lllsud, lllsvc, mxnoni, mxnonr, ndeq, nxptsm
  CHARACTER(8) :: xern1, xern2, xern3, xern4
  !* FIRST EXECUTABLE STATEMENT  BVSUP
  !
  !- *********************************************************************
  !     TEST FOR INVALID INPUT
  !
  IF( Nrowy>=Ncomp ) THEN
    IF( Ncomp==Nic+Nfc ) THEN
      IF( Nxpts>=2 ) THEN
        IF( Nic>0 ) THEN
          IF( Nrowa>=Nic ) THEN
            IF( Nfc>0 ) THEN
              IF( Nrowb>=Nfc ) THEN
                IF( Igofx>=0 .AND. Igofx<=1 ) THEN
                  IF( Re>=0._SP ) THEN
                    IF( Ae>=0._SP ) THEN
                      IF( Re/=0._SP .OR. Ae/=0._SP ) THEN
                        is = 1
                        IF( Xpts(Nxpts)<Xpts(1) ) is = 2
                        nxptsm = Nxpts - 1
                        DO k = 1, nxptsm
                          IF( is==2 ) THEN
                            IF( Xpts(k)<=Xpts(k+1) ) GOTO 100
                          ELSEIF( Xpts(k+1)<=Xpts(k) ) THEN
                            GOTO 100
                          END IF
                        END DO
                        !
                        !- *********************************************************************
                        !     CHECK FOR DISK STORAGE
                        !
                        kpts = Nxpts
                        ndisk_com = 0
                        IF( Iwork(12)/=0 ) THEN
                          ntape_com = Iwork(12)
                          kpts = 1
                          ndisk_com = 1
                        END IF
                        !
                        !- *********************************************************************
                        !     SET INTEG PARAMETER ACCORDING TO CHOICE OF INTEGRATOR.
                        !
                        integ_com = 1
                        IF( Iwork(9)==2 ) integ_com = 2
                        !
                        !- *********************************************************************
                        !     COMPUTE INHOMO
                        !
                        IF( Igofx==1 ) GOTO 300
                        DO j = 1, Nic
                          IF( Alpha(j)/=0._SP ) GOTO 300
                        END DO
                        DO j = 1, Nfc
                          IF( Beta(j)/=0._SP ) GOTO 200
                        END DO
                        inhomo_com = 3
                        GOTO 400
                      END IF
                    END IF
                  END IF
                END IF
              END IF
            END IF
          END IF
        END IF
      END IF
    END IF
  END IF
  100 CONTINUE
  IFlag = -2
  RETURN
  200  inhomo_com = 2
  GOTO 400
  300  inhomo_com = 1
  !
  !- *********************************************************************
  !     TO TAKE ADVANTAGE OF THE SPECIAL STRUCTURE WHEN SOLVING A
  !     COMPLEX VALUED PROBLEM,WE INTRODUCE NFCC=NFC WHILE CHANGING
  !     THE INTERNAL VALUE OF NFC
  !
  400  nfcc_com = Nfc
  IF( Iflag==13 ) Nfc = Nfc/2
  !
  !- *********************************************************************
  !     DETERMINE NECESSARY STORAGE REQUIREMENTS
  !
  ! FOR BASIC ARRAYS IN BVPOR
  kkkyhp = Ncomp*(Nfc+1) + Neqivp
  kkku = Ncomp*Nfc*kpts
  kkkv = Ncomp*kpts
  kkkcoe = nfcc_com
  kkks = Nfc + 1
  kkksto = Ncomp*(Nfc+1) + Neqivp + 1
  kkkg = Ncomp
  !
  ! FOR ORTHONORMALIZATION RELATED MATTERS
  ntp_com = (nfcc_com*(nfcc_com+1))/2
  kkkzpw_com = 1 + ntp_com + nfcc_com
  lllip = nfcc_com
  !
  ! FOR ADDITIONAL REQUIRED WORK SPACE
  !   (LSSUDS)
  kkksud = 4*Nic + (Nrowa+1)*Ncomp
  lllsud = Nic
  !   (SVECS)
  kkksvc = 1 + 4*nfcc_com + 2*nfcc_com**2
  lllsvc = 2*nfcc_com
  !
  ndeq = Ncomp*Nfc + Neqivp
  IF( inhomo_com==1 ) ndeq = ndeq + Ncomp
  IF( integ_com==2 ) THEN
    !   (DEABM)
    kkkint_com = 130 + 21*ndeq
    lllint_com = 51
  ELSE
    !   (DERKF)
    kkkint_com = 33 + 7*ndeq
    lllint_com = 34
  END IF
  !
  !   (COEF)
  kkkcof = 5*nfcc_com + nfcc_com**2
  lllcof = 3 + nfcc_com
  !
  kkkws = MAX(kkksud,kkksvc,kkkint_com,kkkcof)
  llliws = MAX(lllsud,lllsvc,lllint_com,lllcof)
  !
  needw_com = kkkyhp + kkku + kkkv + kkkcoe + kkks + kkksto + kkkg + kkkzpw_com + kkkws
  neediw_com = 17 + lllip + llliws
  !- *********************************************************************
  !     COMPUTE THE NUMBER OF POSSIBLE ORTHONORMALIZATIONS WITH THE
  !     ALLOTTED STORAGE
  !
  Iwork(3) = needw_com
  Iwork(4) = kkkzpw_com
  Iwork(5) = neediw_com
  Iwork(6) = lllip
  nrtemp = Ndw - needw_com
  nitemp = Ndiw - neediw_com
  IF( nrtemp>=0 ) THEN
    IF( nitemp>=0 ) THEN
      !
      IF( ndisk_com==0 ) THEN
        !
        mxnonr = nrtemp/kkkzpw_com
        mxnoni = nitemp/lllip
        mxnon_com = MIN(mxnonr,mxnoni)
        non = mxnon_com
      ELSE
        non = 0
        mxnon_com = nrtemp
      END IF
      !
      Iwork(2) = mxnon_com
      !
      !- *********************************************************************
      !     CHECK FOR PRE-ASSIGNED ORTHONORMALIZATION POINTS
      !
      nopg_com = 0
      IF( Iwork(11)/=1 ) GOTO 500
      IF( mxnon_com>=Iwork(1) ) THEN
        nopg_com = 1
        mxnon_com = Iwork(1)
        Work(mxnon_com+1) = 2._SP*Xpts(Nxpts) - Xpts(1)
        GOTO 500
      END IF
    END IF
  END IF
  !
  Iflag = -1
  IF( ndisk_com/=1 ) THEN
    WRITE (xern1,'(I8)') needw_com
    WRITE (xern2,'(I8)') kkkzpw_com
    WRITE (xern3,'(I8)') neediw_com
    WRITE (xern4,'(I8)') lllip
    ERROR STOP 'BVSUP : REQUIRED STORAGE FOR WORK NOT SATISFIED.&
      & REQUIRED STORAGE FOR IWORK ARRAY NOT SATISFIED.'
  ELSE
    WRITE (xern1,'(I8)') needw_com
    WRITE (xern2,'(I8)') neediw_com
    ERROR STOP 'BVSUP : REQUIRED STORAGE FOR WORK NOT SATISFIED.&
      & REQUIRED STORAGE FOR IWORK ARRAY NOT SATISFIED'
  END IF
  RETURN
  !
  !- *********************************************************************
  !     ALLOCATE STORAGE FROM WORK AND IWORK ARRAYS
  !
  !  (Z)
  500  k1_com = 1 + (mxnon_com+1)
  !  (P)
  k2_com = k1_com + ntp_com*(non+1)
  !  (W)
  k3_com = k2_com + nfcc_com*(non+1)
  !  (YHP)
  k4_com = k3_com + kkkyhp
  !  (U)
  k5_com = k4_com + kkku
  !  (V)
  k6_com = k5_com + kkkv
  !  (COEF)
  k7_com = k6_com + kkkcoe
  !  (S)
  k8_com = k7_com + kkks
  !  (STOWA)
  k9_com = k8_com + kkksto
  !  (G)
  k10_com = k9_com + kkkg
  k11_com = k10_com + kkkws
  !            REQUIRED ADDITIONAL REAL WORK SPACE STARTS AT WORK(K10)
  !            AND EXTENDS TO WORK(K11-1)
  !
  !     FIRST 17 LOCATIONS OF IWORK ARE USED FOR OPTIONAL
  !     INPUT AND OUTPUT ITEMS
  !  (IP)
  l1_com = 18 + nfcc_com*(non+1)
  l2_com = l1_com + llliws
  !            REQUIRED INTEGER WORK SPACE STARTS AT IWORK(L1)
  !            AND EXTENDS TO IWORK(L2-1)
  !
  !- *********************************************************************
  !     SET INDICATOR FOR NORMALIZATION OF PARTICULAR SOLUTION
  !
  nps_com = 0
  IF( Iwork(10)==1 ) nps_com = 1
  !
  !- *********************************************************************
  !     SET PIVOTING PARAMETER
  !
  indpvt_com = 0
  IF( Iwork(15)==1 ) indpvt_com = 1
  !
  !- *********************************************************************
  !     SET OTHER COMMON BLOCK PARAMETERS
  !
  nfc_com = Nfc
  ncomp_com = Ncomp
  igofx_com = Igofx
  nxpts_com = Nxpts
  nic_com = Nic
  re_com = Re
  ae_com = Ae
  neqivp_com = Neqivp
  mnswot_com = 20
  IF( Iwork(13)==-1 ) mnswot_com = MAX(1,Iwork(14))
  xbeg_com = Xpts(1)
  xend_com = Xpts(Nxpts)
  xsav_com = xend_com
  icoco_com = 1
  IF( inhomo_com==3 .AND. nopg_com==1 ) Work(mxnon_com+1) = xend_com
  !
  !- *********************************************************************
  !
  CALL EXBVP(Y,Nrowy,Xpts,A,Nrowa,Alpha,B,Nrowb,Beta,Iflag,Work,Iwork)
  Nfc = nfcc_com
  Iwork(17) = Iwork(l1_com)
  !
END SUBROUTINE BVSUP