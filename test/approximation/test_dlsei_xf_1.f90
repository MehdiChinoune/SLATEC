MODULE TEST_DLSEI_XF_1
  USE service, ONLY : DP
  IMPLICIT NONE

CONTAINS
  SUBROUTINE DLSEIT()
    USE approximation, ONLY : DLSEI
    !     .. Local Scalars ..
    REAL(DP) :: rnorml(1), rnorme
    INTEGER :: idigit, jdigit, ma, mdd, me, meap1, mep1, mg, mode, n, np1
    !     .. Local Arrays ..
    REAL(DP) :: d(11,6), prgopt(4), work(105), x(5)
    INTEGER :: ip(17)
    !     .. Data statements ..
    !
    !     Define the data arrays for the example.  The array A contains
    !     the least squares equations.  (There are no equality constraints
    !     in this example).
    !
    REAL(DP), PARAMETER :: a(6,5) = RESHAPE( [ &
      -74._DP, 80._DP, 18._DP, -11._DP, -4._DP,   14._DP, -69._DP, 21._DP, 28._DP, 0._DP, &
      66._DP, -72._DP, -5._DP, 7._DP, 1._DP,     -12._DP, 66._DP, -30._DP, -23._DP, 3._DP, &
      3._DP, 8._DP, -7._DP, -4._DP, 1._DP,        4._DP, -12._DP, 4._DP, 4._DP, 0._DP ], &
      [6,5], ORDER = [2,1] )
    !
    !     The array G contains the inequality constraint equations,
    !     written in the sense
    !     (row vector)*(solution vector) >= (given value).
    !
    REAL(DP), PARAMETER :: g(5,5) = RESHAPE( [ -1._DP, -1._DP, -1._DP, -1._DP, -1._DP, &
      10._DP, 10._DP, -3._DP, 5._DP, 4._DP,    -8._DP, 1._DP, -2._DP, -5._DP, 3._DP, &
      8._DP, -1._DP, 2._DP, 5._DP, -3._DP,    -4._DP, -2._DP, 3._DP, -5._DP, 1._DP ], &
      [5,5], ORDER = [2,1] )
    !
    !     Define the least squares right-side vector.
    !
    REAL(DP), PARAMETER :: f(6) = [ -5._DP, -9._DP, 708._DP, 4165._DP, &
      -13266._DP, 8409._DP ]
    !
    !     Define the inequality constraint right-side vector.
    !
    REAL(DP), PARAMETER :: h(5) = [ -5._DP, 20._DP, -40._DP, 11._DP, -30._DP ]
    !* FIRST EXECUTABLE STATEMENT  DDLSEIT
    !
    !  Define the matrix dimensions, number of least squares equations,
    !  number of equality constraints, total number of equations, and number
    !  of variables.  Set ME=0 to indicate there are no equality constraints.
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
    d(ma+1:mdd,1:n) = g
    d(1:ma,1:n) = a
    !
    !     Copy the right-side vectors into the work array in compatible order.
    !
    d(ma+1:mdd,np1) = h
    d(1:ma,np1) = f
    !
    !     Use default program options in DLSEI, and set matrix-vector
    !     printing accuracy parameters.
    !
    prgopt(1) = 1
    idigit = -4
    jdigit = -11
    !
    CALL DLSEI(d,0,me,ma,mg,n,prgopt,x,rnorme,rnorml(1),mode,work,ip)
    !
    RETURN
  END SUBROUTINE DLSEIT
  !
END MODULE TEST_DLSEI_XF_1
!
PROGRAM MAIN
  USE TEST_DLSEI_XF_1
  IMPLICIT NONE
  !
  CALL DLSEIT()
  !
END PROGRAM MAIN