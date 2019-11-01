MODULE TEST_LSEI_XF_2
  USE service, ONLY : SP
  IMPLICIT NONE

CONTAINS
  SUBROUTINE LSEIQX()
    !> Quick check for LSEI.
    !***
    ! **Description:**
    !
    !   The sample problem solved is from a paper by J. Stoer, in
    !   SIAM Journal of Numerical Analysis, June 1971.

    USE approximation, ONLY : LSEI
    !     .. Local Scalars ..
    REAL(SP) :: rnorme, rnorml(1)
    INTEGER :: idigit, jdigit, ma, mdd, me, meap1, mep1, mg, mode, n, np1
    !     .. Local Arrays ..
    REAL(SP) :: d(11,6), prgopt(4), work(105), x(5)
    INTEGER :: ip(17)
    !     .. Data statements ..
    !
    !  Define the data arrays for the example.  The array A contains the least
    !  squares equations.  (There are no equality constraints in this example).
    !
    REAL(SP), PARAMETER :: a(6,5) = RESHAPE( [ &
      -74._SP, 80._SP, 18._SP, -11._SP, -4._SP,    14._SP, -69._SP, 21._SP, 28._SP, 0._SP, &
      66._SP, -72._SP, -5._SP, 7._SP, 1._SP,      -12._SP, 66._SP, -30._SP, -23._SP, 3._SP, &
      3._SP, 8._SP, -7._SP, -4._SP, 1._SP,         4._SP, -12._SP, 4._SP, 4._SP, 0._SP ], &
      [6,5], ORDER = [2,1] )
    !
    !  The array G contains the inequality constraint equations, written in the sense
    !     (row vector)*(solution vector) >= (given value).
    !
    REAL(SP), PARAMETER :: g(5,5) = RESHAPE( [ -1._SP, -1._SP, -1._SP, -1._SP, -1._SP, &
      10._SP, 10._SP, -3._SP, 5._SP, 4._SP,    -8._SP, 1._SP, -2._SP, -5._SP, 3._SP, &
      8._SP, -1._SP, 2._SP, 5._SP, -3._SP,     -4._SP, -2._SP, 3._SP, -5._SP, 1._SP ], &
      [5,5], ORDER = [2,1] )
    !
    !     Define the least squares right-side vector.
    !
    REAL(SP), PARAMETER :: f(6) = [ -5._SP, -9._SP, 708._SP, 4165._SP, &
      -13266._SP, 8409._SP ]
    !
    !     Define the inequality constraint right-side vector.
    !
    REAL(SP), PARAMETER :: h(5) = [ -5._SP, 20._SP, -40._SP, 11._SP, -30._SP ]
    !* FIRST EXECUTABLE STATEMENT  LSEIQX
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
    !     Copy the right-side vectors into the work array in compatible
    !     order.
    !
    d(ma+1:mdd,np1) = h
    d(1:ma,np1) = f
    !
    !     Use default program options in LSEI, and set matrix-vector
    !     printing accuracy parameters.
    !
    idigit = -4
    jdigit = -11
    !
    !     Check calls to error processor.
    !
    prgopt(1) = -1
    CALL LSEI(d,mdd,me,ma,mg,n,prgopt,x,rnorme,rnorml(1),mode,work,ip)
    !
    RETURN
  END SUBROUTINE LSEIQX
  !
END MODULE TEST_LSEI_XF_2
!
PROGRAM MAIN
  USE TEST_LSEI_XF_2
  IMPLICIT NONE
  !
  CALL LSEIQX()
  !
END PROGRAM MAIN