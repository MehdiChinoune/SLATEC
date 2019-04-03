!** ISSGMR
INTEGER FUNCTION ISSGMR(N,B,X,Xl,Nelt,Ia,Ja,A,Isym,MSOLVE,Nmsl,Itol,Tol,&
    Itmax,Iter,Err,Iunit,R,Z,Dz,Rwork,Iwork,Rnrm,Bnrm,&
    Sb,Sx,Jscal,Kmp,Lgmr,Maxl,Maxlp1,V,Q,Snormw,Prod,R0nrm,Hes,Jpre)
  USE SSLBLK
  IMPLICIT NONE
  !>
  !***
  !  Generalized Minimum Residual Stop Test.
  !            This routine calculates the stop test for the Generalized
  !            Minimum RESidual (GMRES) iteration scheme.  It returns a
  !            non-zero if the error estimate (the type of which is
  !            determined by ITOL) is less than the user specified
  !            tolerance TOL.
  !***
  ! **Library:**   SLATEC (SLAP)
  !***
  ! **Category:**  D2A4, D2B4
  !***
  ! **Type:**      SINGLE PRECISION (ISSGMR-S, ISDGMR-D)
  !***
  ! **Keywords:**  GMRES, LINEAR SYSTEM, SLAP, SPARSE, STOP TEST
  !***
  ! **Author:**  Brown, Peter, (LLNL), pnbrown@llnl.gov
  !           Hindmarsh, Alan, (LLNL), alanh@llnl.gov
  !           Seager, Mark K., (LLNL), seager@llnl.gov
  !             Lawrence Livermore National Laboratory
  !             PO Box 808, L-60
  !             Livermore, CA 94550 (510) 423-3141
  !***
  ! **Description:**
  !
  !- Usage:
  !      INTEGER N, NELT, IA(NELT), JA(NELT), ISYM, NMSL, ITOL
  !      INTEGER ITMAX, ITER, IUNIT, IWORK(USER DEFINED), JSCAL
  !      INTEGER KMP, LGMR, MAXL, MAXLP1, JPRE
  !      REAL B(N), X(N), XL(MAXL), A(NELT), TOL, ERR, R(N), Z(N),
  !     $     DZ(N), RWORK(USER DEFINED), RNRM, BNRM, SB(N), SX(N),
  !     $     V(N,MAXLP1), Q(2*MAXL), SNORMW, PROD, R0NRM,
  !     $     HES(MAXLP1,MAXL)
  !      EXTERNAL MSOLVE
  !
  !      IF (ISSGMR(N, B, X, XL, NELT, IA, JA, A, ISYM, MSOLVE,
  !     $     NMSL, ITOL, TOL, ITMAX, ITER, ERR, IUNIT, R, Z, DZ,
  !     $     RWORK, IWORK, RNRM, BNRM, SB, SX, JSCAL,
  !     $     KMP, LGMR, MAXL, MAXLP1, V, Q, SNORMW, PROD, R0NRM,
  !     $     HES, JPRE) .NE. 0) THEN ITERATION DONE
  !
  !- Arguments:
  ! N      :IN       Integer.
  !         Order of the Matrix.
  ! B      :IN       Real B(N).
  !         Right-hand-side vector.
  ! X      :IN       Real X(N).
  !         Approximate solution vector as of the last restart.
  ! XL     :OUT      Real XL(N)
  !         An array of length N used to hold the approximate
  !         solution as of the current iteration.  Only computed by
  !         this routine when ITOL=11.
  ! NELT   :IN       Integer.
  !         Number of Non-Zeros stored in A.
  ! IA     :IN       Integer IA(NELT).
  ! JA     :IN       Integer JA(NELT).
  ! A      :IN       Real A(NELT).
  !         These arrays contain the matrix data structure for A.
  !         It could take any form.  See "Description", in the SGMRES,
  !         SSLUGM and SSDGMR routines for more details.
  ! ISYM   :IN       Integer.
  !         Flag to indicate symmetric storage format.
  !         If ISYM=0, all non-zero entries of the matrix are stored.
  !         If ISYM=1, the matrix is symmetric, and only the upper
  !         or lower triangle of the matrix is stored.
  ! MSOLVE :EXT      External.
  !         Name of a routine which solves a linear system Mz = r for  z
  !         given r with the preconditioning matrix M (M is supplied via
  !         RWORK and IWORK arrays.  The name of the MSOLVE routine must
  !         be declared external in the calling program.  The calling
  !         sequence to MSOLVE is:
  !             CALL MSOLVE(N, R, Z, NELT, IA, JA, A, ISYM, RWORK, IWORK)
  !         Where N is the number of unknowns, R is the right-hand side
  !         vector and Z is the solution upon return.  NELT, IA, JA, A and
  !         ISYM are defined as above.  RWORK is a real array that can
  !         be used to pass necessary preconditioning information and/or
  !         workspace to MSOLVE.  IWORK is an integer work array for
  !         the same purpose as RWORK.
  ! NMSL   :INOUT    Integer.
  !         A counter for the number of calls to MSOLVE.
  ! ITOL   :IN       Integer.
  !         Flag to indicate the type of convergence criterion used.
  !         ITOL=0  Means the  iteration stops when the test described
  !                 below on  the  residual RL  is satisfied.  This is
  !                 the  "Natural Stopping Criteria" for this routine.
  !                 Other values  of   ITOL  cause  extra,   otherwise
  !                 unnecessary, computation per iteration and     are
  !                 therefore much less efficient.
  !         ITOL=1  Means   the  iteration stops   when the first test
  !                 described below on  the residual RL  is satisfied,
  !                 and there  is either right  or  no preconditioning
  !                 being used.
  !         ITOL=2  Implies     that   the  user    is   using    left
  !                 preconditioning, and the second stopping criterion
  !                 below is used.
  !         ITOL=3  Means the  iteration stops   when  the  third test
  !                 described below on Minv*Residual is satisfied, and
  !                 there is either left  or no  preconditioning begin
  !                 used.
  !         ITOL=11 is    often  useful  for   checking  and comparing
  !                 different routines.  For this case, the  user must
  !                 supply  the  "exact" solution or  a  very accurate
  !                 approximation (one with  an  error much less  than
  !                 TOL) through a common block,
  !                     COMMON /SSLBLK/ SOLN( )
  !                 If ITOL=11, iteration stops when the 2-norm of the
  !                 difference between the iterative approximation and
  !                 the user-supplied solution  divided by the  2-norm
  !                 of the  user-supplied solution  is  less than TOL.
  !                 Note that this requires  the  user to  set up  the
  !                 "COMMON     /SSLBLK/ SOLN(LENGTH)"  in the calling
  !                 routine.  The routine with this declaration should
  !                 be loaded before the stop test so that the correct
  !                 length is used by  the loader.  This procedure  is
  !                 not standard Fortran and may not work correctly on
  !                 your   system (although  it  has  worked  on every
  !                 system the authors have tried).  If ITOL is not 11
  !                 then this common block is indeed standard Fortran.
  ! TOL    :IN       Real.
  !         Convergence criterion, as described above.
  ! ITMAX  :IN       Integer.
  !         Maximum number of iterations.
  ! ITER   :IN       Integer.
  !         The iteration for which to check for convergence.
  ! ERR    :OUT      Real.
  !         Error estimate of error in final approximate solution, as
  !         defined by ITOL.  Letting norm() denote the Euclidean
  !         norm, ERR is defined as follows..
  !
  !         If ITOL=0, then ERR = norm(SB*(B-A*X(L)))/norm(SB*B),
  !                               for right or no preconditioning, and
  !                         ERR = norm(SB*(M-inverse)*(B-A*X(L)))/
  !                                norm(SB*(M-inverse)*B),
  !                               for left preconditioning.
  !         If ITOL=1, then ERR = norm(SB*(B-A*X(L)))/norm(SB*B),
  !                               since right or no preconditioning
  !                               being used.
  !         If ITOL=2, then ERR = norm(SB*(M-inverse)*(B-A*X(L)))/
  !                                norm(SB*(M-inverse)*B),
  !                               since left preconditioning is being
  !                               used.
  !         If ITOL=3, then ERR =  Max  |(Minv*(B-A*X(L)))(i)/x(i)|
  !                               i=1,n
  !         If ITOL=11, then ERR = norm(SB*(X(L)-SOLN))/norm(SB*SOLN).
  ! IUNIT  :IN       Integer.
  !         Unit number on which to write the error at each iteration,
  !         if this is desired for monitoring convergence.  If unit
  !         number is 0, no writing will occur.
  ! R      :INOUT    Real R(N).
  !         Work array used in calling routine.  It contains
  !         information necessary to compute the residual RL = B-A*XL.
  ! Z      :WORK     Real Z(N).
  !         Workspace used to hold the pseudo-residual M z = r.
  ! DZ     :WORK     Real DZ(N).
  !         Workspace used to hold temporary vector(s).
  ! RWORK  :WORK     Real RWORK(USER DEFINED).
  !         Real array that can be used by MSOLVE.
  ! IWORK  :WORK     Integer IWORK(USER DEFINED).
  !         Integer array that can be used by MSOLVE.
  ! RNRM   :IN       Real.
  !         Norm of the current residual.  Type of norm depends on ITOL.
  ! BNRM   :IN       Real.
  !         Norm of the right hand side.  Type of norm depends on ITOL.
  ! SB     :IN       Real SB(N).
  !         Scaling vector for B.
  ! SX     :IN       Real SX(N).
  !         Scaling vector for X.
  ! JSCAL  :IN       Integer.
  !         Flag indicating if scaling arrays SB and SX are being
  !         used in the calling routine SPIGMR.
  !         JSCAL=0 means SB and SX are not used and the
  !                 algorithm will perform as if all
  !                 SB(i) = 1 and SX(i) = 1.
  !         JSCAL=1 means only SX is used, and the algorithm
  !                 performs as if all SB(i) = 1.
  !         JSCAL=2 means only SB is used, and the algorithm
  !                 performs as if all SX(i) = 1.
  !         JSCAL=3 means both SB and SX are used.
  ! KMP    :IN       Integer
  !         The number of previous vectors the new vector VNEW
  !         must be made orthogonal to.  (KMP .le. MAXL)
  ! LGMR   :IN       Integer
  !         The number of GMRES iterations performed on the current call
  !         to SPIGMR (i.e., # iterations since the last restart) and
  !         the current order of the upper Hessenberg
  !         matrix HES.
  ! MAXL   :IN       Integer
  !         The maximum allowable order of the matrix H.
  ! MAXLP1 :IN       Integer
  !         MAXPL1 = MAXL + 1, used for dynamic dimensioning of HES.
  ! V      :IN       Real V(N,MAXLP1)
  !         The N by (LGMR+1) array containing the LGMR
  !         orthogonal vectors V(*,1) to V(*,LGMR).
  ! Q      :IN       Real Q(2*MAXL)
  !         A real array of length 2*MAXL containing the components
  !         of the Givens rotations used in the QR decomposition
  !         of HES.
  ! SNORMW :IN       Real
  !         A scalar containing the scaled norm of VNEW before it
  !         is renormalized in SPIGMR.
  ! PROD   :IN       Real
  !         The product s1*s2*...*sl = the product of the sines of the
  !         Givens rotations used in the QR factorization of the
  !         Hessenberg matrix HES.
  ! R0NRM  :IN       Real
  !         The scaled norm of initial residual R0.
  ! HES    :IN       Real HES(MAXLP1,MAXL)
  !         The upper triangular factor of the QR decomposition
  !         of the (LGMR+1) by LGMR upper Hessenberg matrix whose
  !         entries are the scaled inner-products of A*V(*,I)
  !         and V(*,K).
  ! JPRE   :IN       Integer
  !         Preconditioner type flag.
  !         (See description of IGWK(4) in SGMRES.)
  !
  !- Description
  !       When using the GMRES solver,  the preferred value  for ITOL
  !       is 0.  This is due to the fact that when ITOL=0 the norm of
  !       the residual required in the stopping test is  obtained for
  !       free, since this value is already  calculated  in the GMRES
  !       algorithm.   The  variable  RNRM contains the   appropriate
  !       norm, which is equal to norm(SB*(RL - A*XL))  when right or
  !       no   preconditioning is  being  performed,   and equal   to
  !       norm(SB*Minv*(RL - A*XL))  when using left preconditioning.
  !       Here, norm() is the Euclidean norm.  Nonzero values of ITOL
  !       require  additional work  to  calculate the  actual  scaled
  !       residual  or its scaled/preconditioned  form,  and/or   the
  !       approximate solution XL.  Hence, these values of  ITOL will
  !       not be as efficient as ITOL=0.
  !
  !- Cautions:
  !     This routine will attempt to write to the Fortran logical output
  !     unit IUNIT, if IUNIT .ne. 0.  Thus, the user must make sure that
  !     this logical unit is attached to a file or terminal before calling
  !     this routine with a non-zero value for IUNIT.  This routine does
  !     not check for the validity of a non-zero IUNIT unit number.
  !
  !     This routine does not verify that ITOL has a valid value.
  !     The calling routine should make such a test before calling
  !     ISSGMR, as is done in SGMRES.
  !
  !***
  ! **See also:**  SGMRES
  !***
  ! **Routines called:**  R1MACH, SCOPY, SNRM2, SRLCAL, SSCAL, SXLCAL
  !***
  ! COMMON BLOCKS    SSLBLK

  !* REVISION HISTORY  (YYMMDD)
  !   871211  DATE WRITTEN
  !   881213  Previous REVISION DATE
  !   890915  Made changes requested at July 1989 CML Meeting.  (MKS)
  !   890922  Numerous changes to prologue to make closer to SLATEC
  !           standard.  (FNF)
  !   890929  Numerous changes to reduce SP/DP differences.  (FNF)
  !   910411  Prologue converted to Version 4.0 format.  (BAB)
  !   910502  Corrected conversion errors, etc.  (FNF)
  !   910502  Removed MSOLVE from ROUTINES CALLED list.  (FNF)
  !   910506  Made subsidiary to SGMRES.  (FNF)
  !   920407  COMMON BLOCK renamed SSLBLK.  (WRB)
  !   920511  Added complete declaration section.  (WRB)
  !   921113  Corrected C***CATEGORY line.  (FNF)

  !     .. Scalar Arguments ..
  REAL Bnrm, Err, Prod, R0nrm, Rnrm, Snormw, Tol
  INTEGER Isym, Iter, Itmax, Itol, Iunit, Jpre, Jscal, Kmp, Lgmr, &
    Maxl, Maxlp1, N, Nelt, Nmsl
  !     .. Array Arguments ..
  REAL A(*), B(*), Dz(*), Hes(Maxlp1,Maxl), Q(*), R(*), Rwork(*), &
    Sb(*), Sx(*), V(N,*), X(*), Xl(*), Z(*)
  INTEGER Ia(*), Iwork(*), Ja(*)
  !     .. Subroutine Arguments ..
  EXTERNAL :: MSOLVE
  !     .. Local Scalars ..
  REAL dxnrm, fuzz, rat, ratmax, tem
  INTEGER i, ielmax
  !     .. External Functions ..
  REAL, EXTERNAL :: R1MACH, SNRM2
  !     .. External Subroutines ..
  EXTERNAL :: SCOPY, SRLCAL, SSCAL, SXLCAL
  !     .. Intrinsic Functions ..
  INTRINSIC ABS, MAX, SQRT
  !     .. Save statement ..
  REAL, SAVE :: solnrm
  !* FIRST EXECUTABLE STATEMENT  ISSGMR
  ISSGMR = 0
  !
  !       Use input from SPIGMR to determine if stop conditions are met.
  !
  IF ( Itol==0 ) Err = Rnrm/Bnrm
  IF ( (Itol>0).AND.(Itol<=3) ) THEN
    !
    !       Use SRLCAL to calculate the scaled residual vector.
    !       Store answer in R.
    !
    IF ( Lgmr/=0 ) CALL SRLCAL(N,Kmp,Lgmr,Maxl,V,Q,R,Snormw,Prod,R0nrm)
    IF ( Itol<=2 ) THEN
      !         err = ||Residual||/||RightHandSide||(2-Norms).
      Err = SNRM2(N,R,1)/Bnrm
      !
      !         Unscale R by R0NRM*PROD when KMP < MAXL.
      !
      IF ( (Kmp<Maxl).AND.(Lgmr/=0) ) THEN
        tem = 1.0E0/(R0nrm*Prod)
        CALL SSCAL(N,tem,R,1)
      ENDIF
    ELSEIF ( Itol==3 ) THEN
      !         err = Max |(Minv*Residual)(i)/x(i)|
      !         When JPRE .lt. 0, R already contains Minv*Residual.
      IF ( Jpre>0 ) THEN
        CALL MSOLVE(N,R,Dz,Nelt,Ia,Ja,A,Isym,Rwork,Iwork)
        Nmsl = Nmsl + 1
      ENDIF
      !
      !         Unscale R by R0NRM*PROD when KMP < MAXL.
      !
      IF ( (Kmp<Maxl).AND.(Lgmr/=0) ) THEN
        tem = 1.0E0/(R0nrm*Prod)
        CALL SSCAL(N,tem,R,1)
      ENDIF
      !
      fuzz = R1MACH(1)
      ielmax = 1
      ratmax = ABS(Dz(1))/MAX(ABS(X(1)),fuzz)
      DO i = 2, N
        rat = ABS(Dz(i))/MAX(ABS(X(i)),fuzz)
        IF ( rat>ratmax ) THEN
          ielmax = i
          ratmax = rat
        ENDIF
      ENDDO
      Err = ratmax
      IF ( ratmax<=Tol ) ISSGMR = 1
      IF ( Iunit>0 ) WRITE (Iunit,99001) Iter, ielmax, ratmax
      99001 FORMAT (1X,' ITER = ',I5,' IELMAX = ',I5,' |R(IELMAX)/X(IELMAX)| = ',&
        E12.5)
      RETURN
    ENDIF
  ENDIF
  IF ( Itol==11 ) THEN
    !
    !       Use SXLCAL to calculate the approximate solution XL.
    !
    IF ( (Lgmr/=0).AND.(Iter>0) ) THEN
      CALL SXLCAL(N,Lgmr,X,Xl,Xl,Hes,Maxlp1,Q,V,R0nrm,Dz,Sx,Jscal,Jpre,&
        MSOLVE,Nmsl,Rwork,Iwork,Nelt,Ia,Ja,A,Isym)
    ELSEIF ( Iter==0 ) THEN
      !         Copy X to XL to check if initial guess is good enough.
      CALL SCOPY(N,X,1,Xl,1)
    ELSE
      !         Return since this is the first call to SPIGMR on a restart.
      RETURN
    ENDIF
    !
    IF ( (Jscal==0).OR.(Jscal==2) ) THEN
      !         err = ||x-TrueSolution||/||TrueSolution||(2-Norms).
      IF ( Iter==0 ) solnrm = SNRM2(N,SOLn,1)
      DO i = 1, N
        Dz(i) = Xl(i) - SOLn(i)
      ENDDO
      Err = SNRM2(N,Dz,1)/solnrm
    ELSE
      IF ( Iter==0 ) THEN
        solnrm = 0
        DO i = 1, N
          solnrm = solnrm + (Sx(i)*SOLn(i))**2
        ENDDO
        solnrm = SQRT(solnrm)
      ENDIF
      dxnrm = 0
      DO i = 1, N
        dxnrm = dxnrm + (Sx(i)*(Xl(i)-SOLn(i)))**2
      ENDDO
      dxnrm = SQRT(dxnrm)
      !         err = ||SX*(x-TrueSolution)||/||SX*TrueSolution|| (2-Norms).
      Err = dxnrm/solnrm
    ENDIF
  ENDIF
  !
  IF ( Iunit/=0 ) THEN
    IF ( Iter==0 ) THEN
      WRITE (Iunit,99002) N, Itol, Maxl, Kmp
      99002 FORMAT (' Generalized Minimum Residual(',I3,I3,') for ','N, ITOL = ',&
        I5,I5,/' ITER','   Natural Err Est','   Error Estimate')
    ENDIF
    WRITE (Iunit,99003) Iter, Rnrm/Bnrm, Err
    99003 FORMAT (1X,I4,1X,E16.7,1X,E16.7)
  ENDIF
  IF ( Err<=Tol ) ISSGMR = 1
  !
  RETURN
  !------------- LAST LINE OF ISSGMR FOLLOWS ----------------------------
END FUNCTION ISSGMR
