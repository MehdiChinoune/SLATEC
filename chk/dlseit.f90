!DECK DLSEIT
SUBROUTINE DLSEIT(Lun,Kprint,Ipass)
  IMPLICIT NONE
  !***BEGIN PROLOGUE  DLSEIT
  !***PURPOSE  Quick check for DLSEI.
  !***LIBRARY   SLATEC
  !***TYPE      DOUBLE PRECISION (LSEIQX-S, DLSEIT-D)
  !***KEYWORDS  QUICK CHECK
  !***AUTHOR  Hanson, R. J., (SNLA)
  !           Haskell, Karen, (SNLA)
  !***DESCRIPTION
  !
  !   The sample problem solved is from a paper by J. Stoer, in
  !   SIAM Journal of Numerical Analysis, June 1971.
  !
  !***ROUTINES CALLED  D1MACH, DAXPY, DCOPY, DDOT, DLSEI, DNRM2, DVOUT
  !***REVISION HISTORY  (YYMMDD)
  !   790216  DATE WRITTEN
  !   890618  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   901010  Restructured using IF-THEN-ELSE-ENDIF, modified tolerances
  !           to use D1MACH(4) rather than D1MACH(3) and cleaned up
  !           FORMATs.  (RWC)
  !   920722  Initialized IP(1) and IP(2) for CALL to DLSEI.  (BKS, WRB)
  !   930214  Declarations sections added, code revised to test error
  !           returns for all values of KPRINT and code polished.  (WRB)
  !***END PROLOGUE  DDLSEIT
  !     .. Scalar Arguments ..
  INTEGER Ipass, Kprint, Lun
  !     .. Local Scalars ..
  REAL(8) :: cnorm, relerr, relnrm, resnrm, rnorme, rnorml, &
    tnorm
  INTEGER i, idigit, jdigit, kontrl, ma, mdd, me, meap1, mep1, mg, &
    mode, n, nerr, np1
  LOGICAL fatal
  !     .. Local Arrays ..
  REAL(8) :: a(6,5), d(11,6), err(5), f(6), g(5,5), h(5), &
    prgopt(4), sol(5), work(105), x(5)
  INTEGER ip(17)
  !     .. External Functions ..
  REAL(8) :: D1MACH, DDOT, DNRM2
  INTEGER NUMXER
  EXTERNAL NUMXER, D1MACH, DDOT, DNRM2
  !     .. External Subroutines ..
  EXTERNAL DAXPY, DCOPY, DLSEI, DVOUT, XGETF, XSETF
  !     .. Intrinsic Functions ..
  INTRINSIC SQRT
  !     .. Data statements ..
  !
  !     Define the data arrays for the example.  The array A contains
  !     the least squares equations.  (There are no equality constraints
  !     in this example).
  !
  DATA a(1,1), a(1,2), a(1,3), a(1,4), a(1,5)/ - 74., 80., 18., &
    -11., -4./
  DATA a(2,1), a(2,2), a(2,3), a(2,4), a(2,5)/14., -69., 21., 28., &
    0./
  DATA a(3,1), a(3,2), a(3,3), a(3,4), a(3,5)/66., -72., -5., 7., &
    1./
  DATA a(4,1), a(4,2), a(4,3), a(4,4), a(4,5)/ - 12., 66., -30., &
    -23., 3./
  DATA a(5,1), a(5,2), a(5,3), a(5,4), a(5,5)/3., 8., -7., -4., 1./
  DATA a(6,1), a(6,2), a(6,3), a(6,4), a(6,5)/4., -12., 4., 4., 0./
  !
  !     The array G contains the inequality constraint equations,
  !     written in the sense
  !     (row vector)*(solution vector) .GE. (given value).
  !
  DATA g(1,1), g(1,2), g(1,3), g(1,4), g(1,5)/ - 1., -1., -1., -1., &
    -1./
  DATA g(2,1), g(2,2), g(2,3), g(2,4), g(2,5)/10., 10., -3., 5., 4./
  DATA g(3,1), g(3,2), g(3,3), g(3,4), g(3,5)/ - 8., 1., -2., -5., &
    3./
  DATA g(4,1), g(4,2), g(4,3), g(4,4), g(4,5)/8., -1., 2., 5., -3./
  DATA g(5,1), g(5,2), g(5,3), g(5,4), g(5,5)/ - 4., -2., 3., -5., &
    1./
  !
  !     Define the least squares right-side vector.
  !
  DATA f(1), f(2), f(3), f(4), f(5), f(6)/ - 5., -9., 708., 4165., &
    -13266., 8409./
  !
  !     Define the inequality constraint right-side vector.
  !
  DATA h(1), h(2), h(3), h(4), h(5)/ - 5., 20., -40., 11., -30./
  !
  !     Define the vector that is the known solution.
  !
  DATA sol(1), sol(2), sol(3), sol(4), sol(5)/1., 2., -1., 3., -4./
  !***FIRST EXECUTABLE STATEMENT  DDLSEIT
  IF ( Kprint>=2 ) WRITE (Lun,99001)
  !
  99001 FORMAT ('1TEST OF SUBROUTINE DLSEI')
  !
  !     Define the matrix dimensions, number of least squares equations,
  !     number of equality constraints, total number of equations, and
  !     number of variables.  Set ME=0 to indicate there are no equality
  !     constraints.
  !
  mdd = 11
  ma = 6
  mg = 5
  n = 5
  me = 0
  !
  ip(1) = 105
  ip(2) = 17
  !
  np1 = n + 1
  mep1 = me + 1
  meap1 = me + ma + 1
  !
  !     Copy the problem matrices.
  !
  DO i = 1, n
    !
    !        Copy the i-th column of the inequality constraint matrix into
    !        the work array.
    !
    CALL DCOPY(mg,g(1,i),1,d(meap1,i),1)
    !
    !        Copy the i-th column of the least squares matrix into the work
    !        array.
    !
    CALL DCOPY(ma,a(1,i),1,d(mep1,i),1)
  ENDDO
  !
  !     Copy the right-side vectors into the work array in compatible
  !     order.
  !
  CALL DCOPY(mg,h,1,d(meap1,np1),1)
  CALL DCOPY(ma,f,1,d(mep1,np1),1)
  !
  !     Use default program options in DLSEI, and set matrix-vector
  !     printing accuracy parameters.
  !
  prgopt(1) = 1
  idigit = -4
  jdigit = -11
  !
  !     Compute residual norm of known least squares solution.
  !     (to be used to check computed residual norm = RNORML.)
  !
  DO i = 1, ma
    work(i) = DDOT(n,d(i,1),mdd,sol,1) - f(i)
  ENDDO
  resnrm = DNRM2(ma,work,1)
  !
  !     Call DLSEI to get solution in X(*), least squares residual in
  !     RNORML.
  !
  CALL DLSEI(d,mdd,me,ma,mg,n,prgopt,x,rnorme,rnorml,mode,work,ip)
  !
  !     Compute relative error in problem variable solution and residual
  !     norm computation.
  !
  tnorm = DNRM2(n,sol,1)
  CALL DCOPY(n,sol,1,err,1)
  CALL DAXPY(n,-1.0D0,x,1,err,1)
  cnorm = DNRM2(n,err,1)
  relerr = cnorm/tnorm
  relnrm = (resnrm-rnorml)/resnrm
  !
  IF ( relerr<=70.0D0*SQRT(D1MACH(4)).AND.relnrm<=5.0D0*D1MACH(4) ) THEN
    Ipass = 1
    IF ( Kprint>=3 ) WRITE (Lun,99002)
    99002   FORMAT (/' DLSEI PASSED TEST')
  ELSE
    Ipass = 0
    IF ( Kprint>=2 ) WRITE (Lun,99003) relerr, relnrm
    99003   FORMAT (/' DLSEI FAILED TEST'/' RELERR = ',1P,D20.6/' RELNRM = ',D20.6)
  ENDIF
  !
  !     Print out known and computed solutions.
  !
  IF ( Kprint>=3 ) THEN
    CALL DVOUT(n,err,'('' RESIDUALS FROM KNOWN LEAST SQUARES SOLUTION'')',&
      idigit)
    CALL DVOUT(n,x,'(/'' SOLUTION COMPUTED BY DLSEI'')',jdigit)
  ENDIF
  !
  IF ( Kprint>=2 ) THEN
    IF ( Kprint/=2.OR.Ipass==0 ) THEN
      !
      !           Print out the known and computed residual norms.
      !
      CALL DVOUT(1,resnrm,&
        '(/'' RESIDUAL NORM OF KNOWN LEAST SQUARES SOLUTION'')',&
        jdigit)
      CALL DVOUT(1,rnorml,'(/'' RESIDUAL NORM COMPUTED BY DLSEI'')',jdigit)
      !
      !           Print out the computed solution relative error.
      !
      CALL DVOUT(1,relerr,'(/'' COMPUTED SOLUTION RELATIVE ERROR'')',idigit)
      !
      !           Print out the computed relative error in residual norm.
      !
      CALL DVOUT(1,relnrm,'(/'' COMPUTED RELATIVE ERROR IN RESIDUAL NORM'')'&
        ,idigit)
    ENDIF
  ENDIF
  !
  !     Check calls to error processor.
  !
  CALL XGETF(kontrl)
  IF ( Kprint<=2 ) THEN
    CALL XSETF(0)
  ELSE
    CALL XSETF(1)
  ENDIF
  fatal = .FALSE.
  CALL XERCLR
  !
  IF ( Kprint>=3 ) WRITE (Lun,99004)
  99004 FORMAT (/' 2 ERROR MESSAGES EXPECTED')
  !
  CALL DLSEI(d,0,me,ma,mg,n,prgopt,x,rnorme,rnorml,mode,work,ip)
  IF ( NUMXER(nerr)/=2 ) THEN
    Ipass = 0
    fatal = .TRUE.
  ENDIF
  CALL XERCLR
  !
  prgopt(1) = -1
  CALL DLSEI(d,mdd,me,ma,mg,n,prgopt,x,rnorme,rnorml,mode,work,ip)
  IF ( NUMXER(nerr)/=2 ) THEN
    Ipass = 0
    fatal = .TRUE.
  ENDIF
  CALL XERCLR
  !
  !     Restore KONTRL and check to see if the tests of error detection
  !     passed.
  !
  CALL XSETF(kontrl)
  IF ( fatal ) THEN
    IF ( Kprint>=2 ) THEN
      WRITE (Lun,99005)
      99005     FORMAT (/' AT LEAST ONE INCORRECT ARGUMENT TEST FAILED')
    ENDIF
  ELSEIF ( Kprint>=3 ) THEN
    WRITE (Lun,99006)
    99006   FORMAT (/' ALL INCORRECT ARGUMENT TESTS PASSED')
  ENDIF
  !
  !     Print PASS/FAIL message.
  !
  IF ( Ipass==1.AND.Kprint>=2 ) WRITE (Lun,99007)
  99007 FORMAT (/' ****************DLSEI PASSED ALL TESTS***************')
  IF ( Ipass==0.AND.Kprint>=1 ) WRITE (Lun,99008)
  99008 FORMAT (/' ****************DLSEI FAILED SOME TESTS**************')
  RETURN
END SUBROUTINE DLSEIT
