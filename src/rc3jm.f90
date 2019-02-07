!*==RC3JM.f90  processed by SPAG 6.72Dc at 11:01 on  6 Feb 2019
!DECK RC3JM
SUBROUTINE RC3JM(L1,L2,L3,M1,M2min,M2max,Thrcof,Ndim,Ier)
  IMPLICIT NONE
  !*--RC3JM5
  !***BEGIN PROLOGUE  RC3JM
  !***PURPOSE  Evaluate the 3j symbol g(M2) = (L1 L2   L3  )
  !                                           (M1 M2 -M1-M2)
  !            for all allowed values of M2, the other parameters
  !            being held fixed.
  !***LIBRARY   SLATEC
  !***CATEGORY  C19
  !***TYPE      SINGLE PRECISION (RC3JM-S, DRC3JM-D)
  !***KEYWORDS  3J COEFFICIENTS, 3J SYMBOLS, CLEBSCH-GORDAN COEFFICIENTS,
  !             RACAH COEFFICIENTS, VECTOR ADDITION COEFFICIENTS,
  !             WIGNER COEFFICIENTS
  !***AUTHOR  Gordon, R. G., Harvard University
  !           Schulten, K., Max Planck Institute
  !***DESCRIPTION
  !
  ! *Usage:
  !
  !        REAL L1, L2, L3, M1, M2MIN, M2MAX, THRCOF(NDIM)
  !        INTEGER NDIM, IER
  !
  !        CALL RC3JM (L1, L2, L3, M1, M2MIN, M2MAX, THRCOF, NDIM, IER)
  !
  ! *Arguments:
  !
  !     L1 :IN      Parameter in 3j symbol.
  !
  !     L2 :IN      Parameter in 3j symbol.
  !
  !     L3 :IN      Parameter in 3j symbol.
  !
  !     M1 :IN      Parameter in 3j symbol.
  !
  !     M2MIN :OUT  Smallest allowable M2 in 3j symbol.
  !
  !     M2MAX :OUT  Largest allowable M2 in 3j symbol.
  !
  !     THRCOF :OUT Set of 3j coefficients generated by evaluating the
  !                 3j symbol for all allowed values of M2.  THRCOF(I)
  !                 will contain g(M2MIN+I-1), I=1,2,...,M2MAX-M2MIN+1.
  !
  !     NDIM :IN    Declared length of THRCOF in calling program.
  !
  !     IER :OUT    Error flag.
  !                 IER=0 No errors.
  !                 IER=1 Either L1.LT.ABS(M1) or L1+ABS(M1) non-integer.
  !                 IER=2 ABS(L1-L2).LE.L3.LE.L1+L2 not satisfied.
  !                 IER=3 L1+L2+L3 not an integer.
  !                 IER=4 M2MAX-M2MIN not an integer.
  !                 IER=5 M2MAX less than M2MIN.
  !                 IER=6 NDIM less than M2MAX-M2MIN+1.
  !
  ! *Description:
  !
  !     Although conventionally the parameters of the vector addition
  !  coefficients satisfy certain restrictions, such as being integers
  !  or integers plus 1/2, the restrictions imposed on input to this
  !  subroutine are somewhat weaker. See, for example, Section 27.9 of
  !  Abramowitz and Stegun or Appendix C of Volume II of A. Messiah.
  !  The restrictions imposed by this subroutine are
  !       1. L1.GE.ABS(M1) and L1+ABS(M1) must be an integer;
  !       2. ABS(L1-L2).LE.L3.LE.L1+L2;
  !       3. L1+L2+L3 must be an integer;
  !       4. M2MAX-M2MIN must be an integer, where
  !          M2MAX=MIN(L2,L3-M1) and M2MIN=MAX(-L2,-L3-M1).
  !  If the conventional restrictions are satisfied, then these
  !  restrictions are met.
  !
  !     The user should be cautious in using input parameters that do
  !  not satisfy the conventional restrictions. For example, the
  !  the subroutine produces values of
  !       g(M2) = (0.75 1.50   1.75  )
  !               (0.25  M2  -0.25-M2)
  !  for M2=-1.5,-0.5,0.5,1.5 but none of the symmetry properties of the
  !  3j symbol, set forth on page 1056 of Messiah, is satisfied.
  !
  !     The subroutine generates g(M2MIN), g(M2MIN+1), ..., g(M2MAX)
  !  where M2MIN and M2MAX are defined above. The sequence g(M2) is
  !  generated by a three-term recurrence algorithm with scaling to
  !  control overflow. Both backward and forward recurrence are used to
  !  maintain numerical stability. The two recurrence sequences are
  !  matched at an interior point and are normalized from the unitary
  !  property of 3j coefficients and Wigner's phase convention.
  !
  !    The algorithm is suited to applications in which large quantum
  !  numbers arise, such as in molecular dynamics.
  !
  !***REFERENCES  1. Abramowitz, M., and Stegun, I. A., Eds., Handbook
  !                  of Mathematical Functions with Formulas, Graphs
  !                  and Mathematical Tables, NBS Applied Mathematics
  !                  Series 55, June 1964 and subsequent printings.
  !               2. Messiah, Albert., Quantum Mechanics, Volume II,
  !                  North-Holland Publishing Company, 1963.
  !               3. Schulten, Klaus and Gordon, Roy G., Exact recursive
  !                  evaluation of 3j and 6j coefficients for quantum-
  !                  mechanical coupling of angular momenta, J Math
  !                  Phys, v 16, no. 10, October 1975, pp. 1961-1970.
  !               4. Schulten, Klaus and Gordon, Roy G., Semiclassical
  !                  approximations to 3j and 6j coefficients for
  !                  quantum-mechanical coupling of angular momenta,
  !                  J Math Phys, v 16, no. 10, October 1975,
  !                  pp. 1971-1988.
  !               5. Schulten, Klaus and Gordon, Roy G., Recursive
  !                  evaluation of 3j and 6j coefficients, Computer
  !                  Phys Comm, v 11, 1976, pp. 269-278.
  !***ROUTINES CALLED  R1MACH, XERMSG
  !***REVISION HISTORY  (YYMMDD)
  !   750101  DATE WRITTEN
  !   880515  SLATEC prologue added by G. C. Nielson, NBS; parameters
  !           HUGE and TINY revised to depend on R1MACH.
  !   891229  Prologue description rewritten; other prologue sections
  !           revised; MMATCH (location of match point for recurrences)
  !           removed from argument list; argument IER changed to serve
  !           only as an error flag (previously, in cases without error,
  !           it returned the number of scalings); number of error codes
  !           increased to provide more precise error information;
  !           program comments revised; SLATEC error handler calls
  !           introduced to enable printing of error messages to meet
  !           SLATEC standards. These changes were done by D. W. Lozier,
  !           M. A. McClain and J. M. Smith of the National Institute
  !           of Standards and Technology, formerly NBS.
  !   910415  Mixed type expressions eliminated; variable C1 initialized;
  !           description of THRCOF expanded. These changes were done by
  !           D. W. Lozier.
  !***END PROLOGUE  RC3JM
  !
  INTEGER Ndim , Ier
  REAL L1 , L2 , L3 , M1 , M2min , M2max , Thrcof(Ndim)
  !
  INTEGER i , index , lstep , n , nfin , nfinp1 , nfinp2 , nfinp3 , nlim , &
    nstep2
  REAL a1 , a1s , c1 , c1old , c2 , cnorm , R1MACH , dv , eps , huge , m2 , &
    m3 , newfac , oldfac , one , ratio , sign1 , sign2 , srhuge , &
    srtiny , sum1 , sum2 , sumbac , sumfor , sumuni , thresh , tiny , &
    two , x , x1 , x2 , x3 , y , y1 , y2 , y3 , zero
  !
  DATA zero , eps , one , two/0.0 , 0.01 , 1.0 , 2.0/
  !
  !***FIRST EXECUTABLE STATEMENT  RC3JM
  Ier = 0
  !  HUGE is the square root of one twentieth of the largest floating
  !  point number, approximately.
  huge = SQRT(R1MACH(2)/20.0)
  srhuge = SQRT(huge)
  tiny = 1.0/huge
  srtiny = 1.0/srhuge
  !
  !     MMATCH = ZERO
  !
  !
  !  Check error conditions 1, 2, and 3.
  IF ( (L1-ABS(M1)+eps<zero).OR.(MOD(L1+ABS(M1)+eps,one)>=eps+eps) ) THEN
    Ier = 1
    CALL XERMSG('SLATEC','RC3JM','L1-ABS(M1) less than zero or '//&
      'L1+ABS(M1) not integer.',Ier,1)
    RETURN
  ELSEIF ( (L1+L2-L3<-eps).OR.(L1-L2+L3<-eps).OR.(-L1+L2+L3<-eps) ) THEN
    Ier = 2
    CALL XERMSG('SLATEC','RC3JM','L1, L2, L3 do not satisfy '//&
      'triangular condition.',Ier,1)
    RETURN
  ELSEIF ( MOD(L1+L2+L3+eps,one)>=eps+eps ) THEN
    Ier = 3
    CALL XERMSG('SLATEC','RC3JM','L1+L2+L3 not integer.',Ier,1)
    RETURN
  ENDIF
  !
  !
  !  Limits for M2
  M2min = MAX(-L2,-L3-M1)
  M2max = MIN(L2,L3-M1)
  !
  !  Check error condition 4.
  IF ( MOD(M2max-M2min+eps,one)>=eps+eps ) THEN
    Ier = 4
    CALL XERMSG('SLATEC','RC3JM','M2MAX-M2MIN not integer.',Ier,1)
    RETURN
  ENDIF
  IF ( M2min<M2max-eps ) THEN
    !
    !  This is reached in case that M1 and M2 take more than one value.
    !     MSCALE = 0
    nfin = INT(M2max-M2min+one+eps)
    IF ( Ndim<nfin ) THEN
      !
      !  Check error condition 6.
      Ier = 6
      CALL XERMSG('SLATEC','RC3JM','Dimension of result array for 3j '//&
        'coefficients too small.',Ier,1)
      RETURN
    ELSE
      !
      !
      !
      !  Start of forward recursion from M2 = M2MIN
      !
      m2 = M2min
      Thrcof(1) = srtiny
      newfac = 0.0
      c1 = 0.0
      sum1 = tiny
      !
      !
      lstep = 1
    ENDIF
  ELSEIF ( M2min<M2max+eps ) THEN
    !
    !
    !  This is reached in case that M2 and M3 can take only one value.
    !     MSCALE = 0
    Thrcof(1) = (-one)**INT(ABS(L2-L3-M1)+eps)/SQRT(L1+L2+L3+one)
    RETURN
  ELSE
    !
    !  Check error condition 5.
    Ier = 5
    CALL XERMSG('SLATEC','RC3JM','M2MIN greater than M2MAX.',Ier,1)
    RETURN
  ENDIF
  DO
    lstep = lstep + 1
    m2 = m2 + one
    m3 = -M1 - m2
    !
    !
    oldfac = newfac
    a1 = (L2-m2+one)*(L2+m2)*(L3+m3+one)*(L3-m3)
    newfac = SQRT(a1)
    !
    !
    dv = (L1+L2+L3+one)*(L2+L3-L1) - (L2-m2+one)*(L3+m3+one) - (L2+m2-one)&
      *(L3-m3-one)
    !
    !
    IF ( lstep>2 ) c1old = ABS(c1)
    c1 = -dv/newfac
    !
    IF ( lstep>2 ) THEN
      !
      !
      c2 = -oldfac/newfac
      !
      !  Recursion to the next 3j coefficient
      x = c1*Thrcof(lstep-1) + c2*Thrcof(lstep-2)
      Thrcof(lstep) = x
      sumfor = sum1
      sum1 = sum1 + x*x
      IF ( lstep/=nfin ) THEN
        !
        !  See if last unnormalized 3j coefficient exceeds SRHUGE
        !
        IF ( ABS(x)>=srhuge ) THEN
          !
          !  This is reached if last 3j coefficient larger than SRHUGE,
          !  so that the recursion series THRCOF(1), ... , THRCOF(LSTEP)
          !  has to be rescaled to prevent overflow
          !
          !     MSCALE = MSCALE + 1
          DO i = 1 , lstep
            IF ( ABS(Thrcof(i))<srtiny ) Thrcof(i) = zero
            Thrcof(i) = Thrcof(i)/srhuge
          ENDDO
          sum1 = sum1/huge
          sumfor = sumfor/huge
          x = x/srhuge
        ENDIF
        !
        !
        !  As long as ABS(C1) is decreasing, the recursion proceeds towards
        !  increasing 3j values and, hence, is numerically stable.  Once
        !  an increase of ABS(C1) is detected, the recursion direction is
        !  reversed.
        !
        IF ( c1old>ABS(c1) ) CYCLE
      ENDIF
      !
      !
      !  Keep three 3j coefficients around MMATCH for comparison later
      !  with backward recursion values.
      !
      !     MMATCH = M2 - 1
      nstep2 = nfin - lstep + 3
      x1 = x
      x2 = Thrcof(lstep-1)
      x3 = Thrcof(lstep-2)
      !
      !  Starting backward recursion from M2MAX taking NSTEP2 steps, so
      !  that forwards and backwards recursion overlap at the three points
      !  M2 = MMATCH+1, MMATCH, MMATCH-1.
      !
      nfinp1 = nfin + 1
      nfinp2 = nfin + 2
      nfinp3 = nfin + 3
      Thrcof(nfin) = srtiny
      sum2 = tiny
      !
      !
      !
      m2 = M2max + two
      lstep = 1
      DO
        lstep = lstep + 1
        m2 = m2 - one
        m3 = -M1 - m2
        oldfac = newfac
        a1s = (L2-m2+two)*(L2+m2-one)*(L3+m3+two)*(L3-m3-one)
        newfac = SQRT(a1s)
        dv = (L1+L2+L3+one)*(L2+L3-L1) - (L2-m2+one)*(L3+m3+one)&
          - (L2+m2-one)*(L3-m3-one)
        c1 = -dv/newfac
        IF ( lstep>2 ) THEN
          !
          c2 = -oldfac/newfac
          !
          !  Recursion to the next 3j coefficient
          !
          y = c1*Thrcof(nfinp2-lstep) + c2*Thrcof(nfinp3-lstep)
          !
          IF ( lstep==nstep2 ) EXIT
          !
          Thrcof(nfinp1-lstep) = y
          sumbac = sum2
          sum2 = sum2 + y*y
          !
          !
          !  See if last 3j coefficient exceeds SRHUGE
          !
          IF ( ABS(y)>=srhuge ) THEN
            !
            !  This is reached if last 3j coefficient larger than SRHUGE,
            !  so that the recursion series THRCOF(NFIN), ... , THRCOF(NFIN-LSTEP+1)
            !  has to be rescaled to prevent overflow.
            !
            !     MSCALE = MSCALE + 1
            DO i = 1 , lstep
              index = nfin - i + 1
              IF ( ABS(Thrcof(index))<srtiny ) Thrcof(index) = zero
              Thrcof(index) = Thrcof(index)/srhuge
            ENDDO
            sum2 = sum2/huge
            !
            sumbac = sumbac/huge
          ENDIF
        ELSE
          !
          !  If M2 = M2MAX + 1 the third term in the recursion equation vanishes
          !
          y = srtiny*c1
          Thrcof(nfin-1) = y
          IF ( lstep==nstep2 ) EXIT
          sumbac = sum2
          sum2 = sum2 + y*y
        ENDIF
      ENDDO
      !
      !
      !
      !  The forward recursion 3j coefficients X1, X2, X3 are to be matched
      !  with the corresponding backward recursion values Y1, Y2, Y3.
      !
      y3 = y
      y2 = Thrcof(nfinp2-lstep)
      y1 = Thrcof(nfinp3-lstep)
      !
      !
      !  Determine now RATIO such that YI = RATIO * XI  (I=1,2,3) holds
      !  with minimal error.
      !
      ratio = (x1*y1+x2*y2+x3*y3)/(x1*x1+x2*x2+x3*x3)
      nlim = nfin - nstep2 + 1
      !
      IF ( ABS(ratio)<one ) THEN
        !
        nlim = nlim + 1
        ratio = one/ratio
        DO n = nlim , nfin
          Thrcof(n) = ratio*Thrcof(n)
        ENDDO
        sumuni = sumfor + ratio*ratio*sumbac
      ELSE
        !
        DO n = 1 , nlim
          Thrcof(n) = ratio*Thrcof(n)
        ENDDO
        sumuni = ratio*ratio*sumfor + sumbac
      ENDIF
      EXIT
    ELSE
      !
      !
      !  If M2 = M2MIN + 1, the third term in the recursion equation vanishes,
      !  hence
      !
      x = srtiny*c1
      Thrcof(2) = x
      sum1 = sum1 + tiny*c1*c1
      IF ( lstep==nfin ) THEN
        !
        sumuni = sum1
        EXIT
      ENDIF
    ENDIF
  ENDDO
  !
  !
  !  Normalize 3j coefficients
  !
  cnorm = one/SQRT((L1+L1+one)*sumuni)
  !
  !  Sign convention for last 3j coefficient determines overall phase
  !
  sign1 = SIGN(one,Thrcof(nfin))
  sign2 = (-one)**INT(ABS(L2-L3-M1)+eps)
  IF ( sign1*sign2<=0 ) cnorm = -cnorm
  !
  IF ( ABS(cnorm)<one ) THEN
    !
    thresh = tiny/ABS(cnorm)
    DO n = 1 , nfin
      IF ( ABS(Thrcof(n))<thresh ) Thrcof(n) = zero
      Thrcof(n) = cnorm*Thrcof(n)
    ENDDO
    GOTO 99999
  ENDIF
  !
  DO n = 1 , nfin
    Thrcof(n) = cnorm*Thrcof(n)
  ENDDO
  RETURN
  !
  !
  !
  99999 END SUBROUTINE RC3JM
