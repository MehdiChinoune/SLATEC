!DECK BVSUP
SUBROUTINE BVSUP(Y,Nrowy,Ncomp,Xpts,Nxpts,A,Nrowa,Alpha,Nic,B,Nrowb,Beta,&
    Nfc,Igofx,Re,Ae,Iflag,Work,Ndw,Iwork,Ndiw,Neqivp)
  IMPLICIT NONE
  REAL A, Ae, AED, Alpha, B, Beta, C, EPS, FOUru, PWCnd, PX, Re ,&
    RED, SQOvfl, SRU, TND, TOL, TWOu, URO, Work
  REAL X, XBEg, XENd, XOP, XOT, Xpts, XSAv, Y
  INTEGER ICOco, Iflag, Igofx, IGOfxd, INDpvt, INFo, INHomo, INTeg ,&
    is, ISTkop, IVP, Iwork, j, k, K1, K10, K11, K2, K3, K4
  INTEGER K5, K6, K7, K8, K9, kkkcoe, kkkcof, kkkg, KKKint, kkks ,&
    kkksto, kkksud, kkksvc, kkku, kkkv, kkkws, kkkyhp, KKKzpw ,&
    KNSwot, KOP
  INTEGER kpts, L1, L2, lllcof, LLLint, lllip, llliws, lllsud ,&
    lllsvc, LOTjp, LPAr, MNSwot, MXNon, mxnoni, mxnonr, Ncomp ,&
    NCOmpd, ndeq, NDIsk, Ndiw
  INTEGER Ndw, NEEdiw, NEEdw, NEQ, NEQivd, Neqivp, Nfc, NFCc, NFCd ,&
    Nic, NICd, nitemp, non, NOPg, NPS, Nrowa, Nrowb, Nrowy ,&
    nrtemp, NSWot
  INTEGER NTApe, NTP, NUMort, Nxpts, NXPtsd, nxptsm
  !***BEGIN PROLOGUE  BVSUP
  !***PURPOSE  Solve a linear two-point boundary value problem using
  !            superposition coupled with an orthonormalization procedure
  !            and a variable-step integration scheme.
  !***LIBRARY   SLATEC
  !***CATEGORY  I1B1
  !***TYPE      SINGLE PRECISION (BVSUP-S, DBVSUP-D)
  !***KEYWORDS  ORTHONORMALIZATION, SHOOTING,
  !             TWO-POINT BOUNDARY VALUE PROBLEM
  !***AUTHOR  Scott, M. R., (SNLA)
  !           Watts, H. A., (SNLA)
  !***DESCRIPTION
  !
  ! **********************************************************************
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
  ! **********************************************************************
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
  ! **********************************************************************
  !     INPUT to BVSUP
  ! **********************************************************************
  !
  !     NROWY = Actual row dimension of Y in calling program.
  !             NROWY must be .GE. NCOMP
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
  !             NROWA must be .GE. NIC.
  !
  !     ALPHA(NIC+NEQIVP) = Boundary conditions at Xinitial.
  !                         If NEQIVP .GT. 0 (see below), the boundary
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
  !             NROWB must be .GE. NFC.
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
  ! **NOTE-  RE and AE should not both be zero.
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
  ! **NOTE --  The amount of storage required is problem dependent and may
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
  ! **NOTE -- Occasionally the coefficients  MATRIX  and/or  G  may be
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
  !            When (NEQIVP .GT. 0) and  MATRIX  is dependent on  U  as
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
  !        B. If  NEQIVP .GT. 0,  UIVP must also be supplied.
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
  !        D. If  NEQIVP .GT. 0  and  IGOFX = 1,  UVEC must be supplied.
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
  ! ****CAUTION -- The user MUST zero out IWORK(1),...,IWORK(15)
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
  ! **NOTE***
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
  ! **********************************************************************
  !     OUTPUT from BVSUP
  ! **********************************************************************
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
  ! **********************************************************************
  !
  !     Necessary machine constants are defined in the function
  !     routine R1MACH. The user must make sure that the values
  !     set in R1MACH are relevant to the computer being used.
  !
  ! **********************************************************************
  !
  !***REFERENCES  M. R. Scott and H. A. Watts, SUPORT - a computer code
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
  !***ROUTINES CALLED  EXBVP, MACON, XERMSG
  !***COMMON BLOCKS    ML15TO, ML17BW, ML18JR, ML5MCO, ML8SZ
  !***REVISION HISTORY  (YYMMDD)
  !   750601  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890831  Modified array declarations.  (WRB)
  !   890921  Realigned order of variables in certain COMMON blocks.
  !           (WRB)
  !   890921  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900510  Convert XERRWV calls to XERMSG calls.  (RWC)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  !***END PROLOGUE  BVSUP
  ! **********************************************************************
  !
  !
  DIMENSION Y(Nrowy,*), A(Nrowa,*), Alpha(*), B(Nrowb,*), Beta(*) ,&
    Work(*), Iwork(*), Xpts(*)
  CHARACTER(8) :: xern1, xern2, xern3, xern4
  !
  ! **********************************************************************
  !     THE COMMON BLOCK BELOW IS USED TO COMMUNICATE WITH SUBROUTINE
  !     BVDER.  THE USER SHOULD NOT ALTER OR USE THIS COMMON BLOCK IN THE
  !     CALLING PROGRAM.
  !
  COMMON /ML8SZ / C, XSAv, IGOfxd, INHomo, IVP, NCOmpd, NFCd
  !
  ! **********************************************************************
  !     THESE COMMON BLOCKS AID IN REDUCING THE NUMBER OF SUBROUTINE
  !     ARGUMENTS PREVALENT IN THIS MODULAR STRUCTURE
  !
  COMMON /ML18JR/ AED, RED, TOL, NXPtsd, NICd, NOPg, MXNon, NDIsk ,&
    NTApe, NEQ, INDpvt, INTeg, NPS, NTP, NEQivd ,&
    NUMort, NFCc, ICOco
  COMMON /ML17BW/ KKKzpw, NEEdw, NEEdiw, K1, K2, K3, K4, K5, K6 ,&
    K7, K8, K9, K10, K11, L1, L2, KKKint, LLLint
  !
  ! **********************************************************************
  !     THIS COMMON BLOCK IS USED IN SUBROUTINES BVSUP,BVPOR,RKFAB,
  !     REORT, AND STWAY. IT CONTAINS INFORMATION NECESSARY
  !     FOR THE ORTHONORMALIZATION TESTING PROCEDURE AND A BACKUP
  !     RESTARTING CAPABILITY.
  !
  COMMON /ML15TO/ PX, PWCnd, TND, X, XBEg, XENd, XOT, XOP, INFo(15)&
    , ISTkop, KNSwot, KOP, LOTjp, MNSwot, NSWot
  !
  ! **********************************************************************
  !     THIS COMMON BLOCK CONTAINS THE MACHINE DEPENDENT PARAMETERS
  !     USED BY THE CODE
  !
  COMMON /ML5MCO/ URO, SRU, EPS, SQOvfl, TWOu, FOUru, LPAr
  !
  ! **********************************************************************
  !     SET UP MACHINE DEPENDENT CONSTANTS.
  !
  !***FIRST EXECUTABLE STATEMENT  BVSUP
  CALL MACON
  !
  ! **********************************************************************
  !     TEST FOR INVALID INPUT
  !
  IF ( Nrowy>=Ncomp ) THEN
    IF ( Ncomp==Nic+Nfc ) THEN
      IF ( Nxpts>=2 ) THEN
        IF ( Nic>0 ) THEN
          IF ( Nrowa>=Nic ) THEN
            IF ( Nfc>0 ) THEN
              IF ( Nrowb>=Nfc ) THEN
                IF ( Igofx>=0.AND.Igofx<=1 ) THEN
                  IF ( Re>=0.0 ) THEN
                    IF ( Ae>=0.0 ) THEN
                      IF ( Re/=0.0.OR.Ae/=0.0 ) THEN
                        is = 1
                        IF ( Xpts(Nxpts)<Xpts(1) ) is = 2
                        nxptsm = Nxpts - 1
                        DO k = 1, nxptsm
                          IF ( is==2 ) THEN
                            IF ( Xpts(k)<=Xpts(k+1) ) GOTO 100
                          ELSEIF ( Xpts(k+1)<=Xpts(k) ) THEN
                            GOTO 100
                          ENDIF
                        ENDDO
                        !
                        ! **********************************************************************
                        !     CHECK FOR DISK STORAGE
                        !
                        kpts = Nxpts
                        NDIsk = 0
                        IF ( Iwork(12)/=0 ) THEN
                          NTApe = Iwork(12)
                          kpts = 1
                          NDIsk = 1
                        ENDIF
                        !
                        ! **********************************************************************
                        !     SET INTEG PARAMETER ACCORDING TO CHOICE OF INTEGRATOR.
                        !
                        INTeg = 1
                        IF ( Iwork(9)==2 ) INTeg = 2
                        !
                        ! **********************************************************************
                        !     COMPUTE INHOMO
                        !
                        IF ( Igofx==1 ) GOTO 300
                        DO j = 1, Nic
                          IF ( Alpha(j)/=0.0 ) GOTO 300
                        ENDDO
                        DO j = 1, Nfc
                          IF ( Beta(j)/=0.0 ) GOTO 200
                        ENDDO
                        INHomo = 3
                        GOTO 400
                      ENDIF
                    ENDIF
                  ENDIF
                ENDIF
              ENDIF
            ENDIF
          ENDIF
        ENDIF
      ENDIF
    ENDIF
  ENDIF
  100 CONTINUE
  IFlag = -2
  RETURN
  200  INHomo = 2
  GOTO 400
  300  INHomo = 1
  !
  ! **********************************************************************
  !     TO TAKE ADVANTAGE OF THE SPECIAL STRUCTURE WHEN SOLVING A
  !     COMPLEX VALUED PROBLEM,WE INTRODUCE NFCC=NFC WHILE CHANGING
  !     THE INTERNAL VALUE OF NFC
  !
  400  NFCc = Nfc
  IF ( Iflag==13 ) Nfc = Nfc/2
  !
  ! **********************************************************************
  !     DETERMINE NECESSARY STORAGE REQUIREMENTS
  !
  ! FOR BASIC ARRAYS IN BVPOR
  kkkyhp = Ncomp*(Nfc+1) + Neqivp
  kkku = Ncomp*Nfc*kpts
  kkkv = Ncomp*kpts
  kkkcoe = NFCc
  kkks = Nfc + 1
  kkksto = Ncomp*(Nfc+1) + Neqivp + 1
  kkkg = Ncomp
  !
  ! FOR ORTHONORMALIZATION RELATED MATTERS
  NTP = (NFCc*(NFCc+1))/2
  KKKzpw = 1 + NTP + NFCc
  lllip = NFCc
  !
  ! FOR ADDITIONAL REQUIRED WORK SPACE
  !   (LSSUDS)
  kkksud = 4*Nic + (Nrowa+1)*Ncomp
  lllsud = Nic
  !   (SVECS)
  kkksvc = 1 + 4*NFCc + 2*NFCc**2
  lllsvc = 2*NFCc
  !
  ndeq = Ncomp*Nfc + Neqivp
  IF ( INHomo==1 ) ndeq = ndeq + Ncomp
  IF ( INTeg==2 ) THEN
    !   (DEABM)
    KKKint = 130 + 21*ndeq
    LLLint = 51
  ELSE
    !   (DERKF)
    KKKint = 33 + 7*ndeq
    LLLint = 34
  ENDIF
  !
  !   (COEF)
  kkkcof = 5*NFCc + NFCc**2
  lllcof = 3 + NFCc
  !
  kkkws = MAX(kkksud,kkksvc,KKKint,kkkcof)
  llliws = MAX(lllsud,lllsvc,LLLint,lllcof)
  !
  NEEdw = kkkyhp + kkku + kkkv + kkkcoe + kkks + kkksto + kkkg + KKKzpw +&
    kkkws
  NEEdiw = 17 + lllip + llliws
  ! **********************************************************************
  !     COMPUTE THE NUMBER OF POSSIBLE ORTHONORMALIZATIONS WITH THE
  !     ALLOTTED STORAGE
  !
  Iwork(3) = NEEdw
  Iwork(4) = KKKzpw
  Iwork(5) = NEEdiw
  Iwork(6) = lllip
  nrtemp = Ndw - NEEdw
  nitemp = Ndiw - NEEdiw
  IF ( nrtemp>=0 ) THEN
    IF ( nitemp>=0 ) THEN
      !
      IF ( NDIsk==0 ) THEN
        !
        mxnonr = nrtemp/KKKzpw
        mxnoni = nitemp/lllip
        MXNon = MIN(mxnonr,mxnoni)
        non = MXNon
      ELSE
        non = 0
        MXNon = nrtemp
      ENDIF
      !
      Iwork(2) = MXNon
      !
      ! **********************************************************************
      !     CHECK FOR PRE-ASSIGNED ORTHONORMALIZATION POINTS
      !
      NOPg = 0
      IF ( Iwork(11)/=1 ) GOTO 500
      IF ( MXNon>=Iwork(1) ) THEN
        NOPg = 1
        MXNon = Iwork(1)
        Work(MXNon+1) = 2.*Xpts(Nxpts) - Xpts(1)
        GOTO 500
      ENDIF
    ENDIF
  ENDIF
  !
  Iflag = -1
  IF ( NDIsk/=1 ) THEN
    WRITE (xern1,'(I8)') NEEdw
    WRITE (xern2,'(I8)') KKKzpw
    WRITE (xern3,'(I8)') NEEdiw
    WRITE (xern4,'(I8)') lllip
    CALL XERMSG('SLATEC','BVSUP','REQUIRED STORAGE FOR WORK ARRAY IS '//&
      xern1//' + '//xern2//&
      '*(EXPECTED NUMBER OF ORTHONORMALIZATIONS) $$REQUIRED STORAGE FOR IWORK ARRAY IS '//xern3//' + '//&
      xern4//'*(EXPECTED NUMBER OF ORTHONORMALIZATIONS)',1,0)
  ELSE
    WRITE (xern1,'(I8)') NEEdw
    WRITE (xern2,'(I8)') NEEdiw
    CALL XERMSG('SLATEC','BVSUP','REQUIRED STORAGE FOR WORK ARRAY IS '//&
      xern1//' + NUMBER OF ORTHONOMALIZATIONS. $$REQUIRED STORAGE FOR IWORK ARRAY IS '//xern2,1,0)
  ENDIF
  RETURN
  !
  ! **********************************************************************
  !     ALLOCATE STORAGE FROM WORK AND IWORK ARRAYS
  !
  !  (Z)
  500  K1 = 1 + (MXNon+1)
  !  (P)
  K2 = K1 + NTP*(non+1)
  !  (W)
  K3 = K2 + NFCc*(non+1)
  !  (YHP)
  K4 = K3 + kkkyhp
  !  (U)
  K5 = K4 + kkku
  !  (V)
  K6 = K5 + kkkv
  !  (COEF)
  K7 = K6 + kkkcoe
  !  (S)
  K8 = K7 + kkks
  !  (STOWA)
  K9 = K8 + kkksto
  !  (G)
  K10 = K9 + kkkg
  K11 = K10 + kkkws
  !            REQUIRED ADDITIONAL REAL WORK SPACE STARTS AT WORK(K10)
  !            AND EXTENDS TO WORK(K11-1)
  !
  !     FIRST 17 LOCATIONS OF IWORK ARE USED FOR OPTIONAL
  !     INPUT AND OUTPUT ITEMS
  !  (IP)
  L1 = 18 + NFCc*(non+1)
  L2 = L1 + llliws
  !            REQUIRED INTEGER WORK SPACE STARTS AT IWORK(L1)
  !            AND EXTENDS TO IWORK(L2-1)
  !
  ! **********************************************************************
  !     SET INDICATOR FOR NORMALIZATION OF PARTICULAR SOLUTION
  !
  NPS = 0
  IF ( Iwork(10)==1 ) NPS = 1
  !
  ! **********************************************************************
  !     SET PIVOTING PARAMETER
  !
  INDpvt = 0
  IF ( Iwork(15)==1 ) INDpvt = 1
  !
  ! **********************************************************************
  !     SET OTHER COMMON BLOCK PARAMETERS
  !
  NFCd = Nfc
  NCOmpd = Ncomp
  IGOfxd = Igofx
  NXPtsd = Nxpts
  NICd = Nic
  RED = Re
  AED = Ae
  NEQivd = Neqivp
  MNSwot = 20
  IF ( Iwork(13)==-1 ) MNSwot = MAX(1,Iwork(14))
  XBEg = Xpts(1)
  XENd = Xpts(Nxpts)
  XSAv = XENd
  ICOco = 1
  IF ( INHomo==3.AND.NOPg==1 ) Work(MXNon+1) = XENd
  !
  ! **********************************************************************
  !
  CALL EXBVP(Y,Nrowy,Xpts,A,Nrowa,Alpha,B,Nrowb,Beta,Iflag,Work,Iwork)
  Nfc = NFCc
  Iwork(17) = Iwork(L1)
END SUBROUTINE BVSUP
