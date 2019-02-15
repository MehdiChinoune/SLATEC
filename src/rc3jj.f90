!DECK RC3JJ
SUBROUTINE RC3JJ(L2,L3,M2,M3,L1min,L1max,Thrcof,Ndim,Ier)
  IMPLICIT NONE
  !***BEGIN PROLOGUE  RC3JJ
  !***PURPOSE  Evaluate the 3j symbol f(L1) = (  L1   L2 L3)
  !                                           (-M2-M3 M2 M3)
  !            for all allowed values of L1, the other parameters
  !            being held fixed.
  !***LIBRARY   SLATEC
  !***CATEGORY  C19
  !***TYPE      SINGLE PRECISION (RC3JJ-S, DRC3JJ-D)
  !***KEYWORDS  3J COEFFICIENTS, 3J SYMBOLS, CLEBSCH-GORDAN COEFFICIENTS,
  !             RACAH COEFFICIENTS, VECTOR ADDITION COEFFICIENTS,
  !             WIGNER COEFFICIENTS
  !***AUTHOR  Gordon, R. G., Harvard University
  !           Schulten, K., Max Planck Institute
  !***DESCRIPTION
  !
  ! *Usage:
  !
  !        REAL L2, L3, M2, M3, L1MIN, L1MAX, THRCOF(NDIM)
  !        INTEGER NDIM, IER
  !
  !        CALL RC3JJ (L2, L3, M2, M3, L1MIN, L1MAX, THRCOF, NDIM, IER)
  !
  ! *Arguments:
  !
  !     L2 :IN      Parameter in 3j symbol.
  !
  !     L3 :IN      Parameter in 3j symbol.
  !
  !     M2 :IN      Parameter in 3j symbol.
  !
  !     M3 :IN      Parameter in 3j symbol.
  !
  !     L1MIN :OUT  Smallest allowable L1 in 3j symbol.
  !
  !     L1MAX :OUT  Largest allowable L1 in 3j symbol.
  !
  !     THRCOF :OUT Set of 3j coefficients generated by evaluating the
  !                 3j symbol for all allowed values of L1.  THRCOF(I)
  !                 will contain f(L1MIN+I-1), I=1,2,...,L1MAX+L1MIN+1.
  !
  !     NDIM :IN    Declared length of THRCOF in calling program.
  !
  !     IER :OUT    Error flag.
  !                 IER=0 No errors.
  !                 IER=1 Either L2.LT.ABS(M2) or L3.LT.ABS(M3).
  !                 IER=2 Either L2+ABS(M2) or L3+ABS(M3) non-integer.
  !                 IER=3 L1MAX-L1MIN not an integer.
  !                 IER=4 L1MAX less than L1MIN.
  !                 IER=5 NDIM less than L1MAX-L1MIN+1.
  !
  ! *Description:
  !
  !     Although conventionally the parameters of the vector addition
  !  coefficients satisfy certain restrictions, such as being integers
  !  or integers plus 1/2, the restrictions imposed on input to this
  !  subroutine are somewhat weaker. See, for example, Section 27.9 of
  !  Abramowitz and Stegun or Appendix C of Volume II of A. Messiah.
  !  The restrictions imposed by this subroutine are
  !       1. L2 .GE. ABS(M2) and L3 .GE. ABS(M3);
  !       2. L2+ABS(M2) and L3+ABS(M3) must be integers;
  !       3. L1MAX-L1MIN must be a non-negative integer, where
  !          L1MAX=L2+L3 and L1MIN=MAX(ABS(L2-L3),ABS(M2+M3)).
  !  If the conventional restrictions are satisfied, then these
  !  restrictions are met.
  !
  !     The user should be cautious in using input parameters that do
  !  not satisfy the conventional restrictions. For example, the
  !  the subroutine produces values of
  !       f(L1) = ( L1  2.5  5.8)
  !               (-0.3 1.5 -1.2)
  !  for L1=3.3,4.3,...,8.3 but none of the symmetry properties of the 3j
  !  symbol, set forth on page 1056 of Messiah, is satisfied.
  !
  !     The subroutine generates f(L1MIN), f(L1MIN+1), ..., f(L1MAX)
  !  where L1MIN and L1MAX are defined above. The sequence f(L1) is
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
  !                  approximations to 3j  and 6j coefficients for
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
  !           revised; LMATCH (location of match point for recurrences)
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
  !***END PROLOGUE  RC3JJ
  !
  INTEGER Ndim, Ier
  REAL L2, L3, M2, M3, L1min, L1max, Thrcof(Ndim)
  !
  INTEGER i, index, lstep, n, nfin, nfinp1, nfinp2, nfinp3, nlim, &
    nstep2
  REAL a1, a1s, a2, a2s, c1, c1old, c2, cnorm, R1MACH, denom, dv, &
    eps, huge, l1, m1, newfac, oldfac, one, ratio, sign1, &
    sign2, srhuge, srtiny, sum1, sum2, sumbac, sumfor, sumuni, &
    three, thresh, tiny, two, x, x1, x2, x3, y, y1, y2, y3, &
    zero
  !
  DATA zero, eps, one, two, three/0.0, 0.01, 1.0, 2.0, 3.0/
  !
  !***FIRST EXECUTABLE STATEMENT  RC3JJ
  Ier = 0
  !  HUGE is the square root of one twentieth of the largest floating
  !  point number, approximately.
  huge = SQRT(R1MACH(2)/20.0)
  srhuge = SQRT(huge)
  tiny = 1.0/huge
  srtiny = 1.0/srhuge
  !
  !     LMATCH = ZERO
  m1 = -M2 - M3
  !
  !  Check error conditions 1 and 2.
  IF ( (L2-ABS(M2)+eps<zero).OR.(L3-ABS(M3)+eps<zero) ) THEN
    Ier = 1
    CALL XERMSG('SLATEC','RC3JJ','L2-ABS(M2) or L3-ABS(M3) '//&
      'less than zero.',Ier,1)
    RETURN
  ELSEIF ( (MOD(L2+ABS(M2)+eps,one)>=eps+eps).OR.&
      (MOD(L3+ABS(M3)+eps,one)>=eps+eps) ) THEN
    Ier = 2
    CALL XERMSG('SLATEC','RC3JJ','L2+ABS(M2) or L3+ABS(M3) '//&
      'not integer.',Ier,1)
    RETURN
  ENDIF
  !
  !
  !
  !  Limits for L1
  !
  L1min = MAX(ABS(L2-L3),ABS(m1))
  L1max = L2 + L3
  !
  !  Check error condition 3.
  IF ( MOD(L1max-L1min+eps,one)>=eps+eps ) THEN
    Ier = 3
    CALL XERMSG('SLATEC','RC3JJ','L1MAX-L1MIN not integer.',Ier,1)
    RETURN
  ENDIF
  IF ( L1min<L1max-eps ) THEN
    !
    !  This is reached in case that L1 takes more than one value,
    !  i.e. L1MIN < L1MAX.
    !
    !     LSCALE = 0
    nfin = INT(L1max-L1min+one+eps)
    IF ( Ndim<nfin ) THEN
      !
      !  Check error condition 5.
      Ier = 5
      CALL XERMSG('SLATEC','RC3JJ','Dimension of result array for 3j '//&
        'coefficients too small.',Ier,1)
      RETURN
    ELSE
      !
      !
      !  Starting forward recursion from L1MIN taking NSTEP1 steps
      !
      l1 = L1min
      newfac = 0.0
      c1 = 0.0
      Thrcof(1) = srtiny
      sum1 = (l1+l1+one)*tiny
      !
      !
      lstep = 1
    ENDIF
  ELSEIF ( L1min<L1max+eps ) THEN
    !
    !  This is reached in case that L1 can take only one value,
    !  i.e. L1MIN = L1MAX
    !
    !     LSCALE = 0
    Thrcof(1) = (-one)**INT(ABS(L2+M2-L3+M3)+eps)/SQRT(L1min+L2+L3+one)
    RETURN
  ELSE
    !
    !  Check error condition 4.
    Ier = 4
    CALL XERMSG('SLATEC','RC3JJ','L1MIN greater than L1MAX.',Ier,1)
    RETURN
  ENDIF
  100  lstep = lstep + 1
  l1 = l1 + one
  !
  !
  oldfac = newfac
  a1 = (l1+L2+L3+one)*(l1-L2+L3)*(l1+L2-L3)*(-l1+L2+L3+one)
  a2 = (l1+m1)*(l1-m1)
  newfac = SQRT(a1*a2)
  IF ( l1<one+eps ) THEN
    !
    !  If L1 = 1, (L1-1) has to be factored out of DV, hence
    !
    c1 = -(l1+l1-one)*l1*(M3-M2)/newfac
  ELSE
    !
    !
    dv = -L2*(L2+one)*m1 + L3*(L3+one)*m1 + l1*(l1-one)*(M3-M2)
    denom = (l1-one)*newfac
    !
    !
    IF ( lstep>2 ) c1old = ABS(c1)
    c1 = -(l1+l1-one)*dv/denom
  ENDIF
  !
  IF ( lstep>2 ) THEN
    !
    !
    c2 = -l1*oldfac/denom
    !
    !  Recursion to the next 3j coefficient X
    !
    x = c1*Thrcof(lstep-1) + c2*Thrcof(lstep-2)
    Thrcof(lstep) = x
    sumfor = sum1
    sum1 = sum1 + (l1+l1+one)*x*x
    IF ( lstep/=nfin ) THEN
      !
      !  See if last unnormalized 3j coefficient exceeds SRHUGE
      !
      IF ( ABS(x)>=srhuge ) THEN
        !
        !  This is reached if last 3j coefficient larger than SRHUGE,
        !  so that the recursion series THRCOF(1), ..., THRCOF(LSTEP)
        !  has to be rescaled to prevent overflow
        !
        !     LSCALE = LSCALE + 1
        DO i = 1, lstep
          IF ( ABS(Thrcof(i))<srtiny ) Thrcof(i) = zero
          Thrcof(i) = Thrcof(i)/srhuge
        ENDDO
        sum1 = sum1/huge
        sumfor = sumfor/huge
        x = x/srhuge
      ENDIF
      !
      !  As long as ABS(C1) is decreasing, the recursion proceeds towards
      !  increasing 3j values and, hence, is numerically stable.  Once
      !  an increase of ABS(C1) is detected, the recursion direction is
      !  reversed.
      !
      IF ( c1old>ABS(c1) ) GOTO 100
    ENDIF
    !
    !
    !  Keep three 3j coefficients around LMATCH for comparison with
    !  backward recursion.
    !
    !     LMATCH = L1 - 1
    x1 = x
    x2 = Thrcof(lstep-1)
    x3 = Thrcof(lstep-2)
    nstep2 = nfin - lstep + 3
    !
    !
    !
    !
    !  Starting backward recursion from L1MAX taking NSTEP2 steps, so
    !  that forward and backward recursion overlap at three points
    !  L1 = LMATCH+1, LMATCH, LMATCH-1.
    !
    nfinp1 = nfin + 1
    nfinp2 = nfin + 2
    nfinp3 = nfin + 3
    l1 = L1max
    Thrcof(nfin) = srtiny
    sum2 = tiny*(l1+l1+one)
    !
    l1 = l1 + two
    lstep = 1
    DO
      lstep = lstep + 1
      l1 = l1 - one
      !
      oldfac = newfac
      a1s = (l1+L2+L3)*(l1-L2+L3-one)*(l1+L2-L3-one)*(-l1+L2+L3+two)
      a2s = (l1+m1-one)*(l1-m1-one)
      newfac = SQRT(a1s*a2s)
      !
      dv = -L2*(L2+one)*m1 + L3*(L3+one)*m1 + l1*(l1-one)*(M3-M2)
      !
      denom = l1*newfac
      c1 = -(l1+l1-one)*dv/denom
      IF ( lstep>2 ) THEN
        !
        !
        c2 = -(l1-one)*oldfac/denom
        !
        !  Recursion to the next 3j coefficient Y
        !
        y = c1*Thrcof(nfinp2-lstep) + c2*Thrcof(nfinp3-lstep)
        !
        IF ( lstep==nstep2 ) THEN
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
            DO n = nlim, nfin
              Thrcof(n) = ratio*Thrcof(n)
            ENDDO
            sumuni = sumfor + ratio*ratio*sumbac
          ELSE
            !
            DO n = 1, nlim
              Thrcof(n) = ratio*Thrcof(n)
            ENDDO
            sumuni = ratio*ratio*sumfor + sumbac
          ENDIF
          EXIT
        ELSE
          !
          Thrcof(nfinp1-lstep) = y
          sumbac = sum2
          sum2 = sum2 + (l1+l1-three)*y*y
          !
          !  See if last unnormalized 3j coefficient exceeds SRHUGE
          !
          IF ( ABS(y)>=srhuge ) THEN
            !
            !  This is reached if last 3j coefficient larger than SRHUGE,
            !  so that the recursion series THRCOF(NFIN), ... ,THRCOF(NFIN-LSTEP+1)
            !  has to be rescaled to prevent overflow
            !
            !     LSCALE = LSCALE + 1
            DO i = 1, lstep
              index = nfin - i + 1
              IF ( ABS(Thrcof(index))<srtiny ) Thrcof(index) = zero
              Thrcof(index) = Thrcof(index)/srhuge
            ENDDO
            sum2 = sum2/huge
            !
            !
            sumbac = sumbac/huge
          ENDIF
        ENDIF
      ELSE
        !
        !  If L1 = L1MAX + 1, the third term in the recursion formula vanishes
        !
        y = srtiny*c1
        Thrcof(nfin-1) = y
        sumbac = sum2
        !
        sum2 = sum2 + tiny*(l1+l1-three)*c1*c1
      ENDIF
    ENDDO
  ELSE
    !
    !
    !  If L1 = L1MIN + 1, the third term in the recursion equation vanishes,
    !  hence
    x = srtiny*c1
    Thrcof(2) = x
    sum1 = sum1 + tiny*(l1+l1+one)*c1*c1
    IF ( lstep/=nfin ) GOTO 100
    !
    sumuni = sum1
  ENDIF
  !
  !
  !  Normalize 3j coefficients
  !
  cnorm = one/SQRT(sumuni)
  !
  !  Sign convention for last 3j coefficient determines overall phase
  !
  sign1 = SIGN(one,Thrcof(nfin))
  sign2 = (-one)**INT(ABS(L2+M2-L3+M3)+eps)
  IF ( sign1*sign2<=0 ) cnorm = -cnorm
  !
  IF ( ABS(cnorm)<one ) THEN
    !
    thresh = tiny/ABS(cnorm)
    DO n = 1, nfin
      IF ( ABS(Thrcof(n))<thresh ) Thrcof(n) = zero
      Thrcof(n) = cnorm*Thrcof(n)
    ENDDO
    GOTO 99999
  ENDIF
  !
  DO n = 1, nfin
    Thrcof(n) = cnorm*Thrcof(n)
  ENDDO
  RETURN
  !
  99999 CONTINUE
  END SUBROUTINE RC3JJ
