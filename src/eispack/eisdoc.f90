!** EISDOC
SUBROUTINE EISDOC
  IMPLICIT NONE
  !>
  !***
  !  Documentation for EISPACK, a collection of subprograms for
  !            solving matrix eigen-problems.
  !***
  ! **Library:**   SLATEC (EISPACK)
  !***
  ! **Category:**  D4, Z
  !***
  ! **Type:**      ALL (EISDOC-A)
  !***
  ! **Keywords:**  EIGENVALUES, EIGENVECTORS, EISPACK
  !***
  ! **Author:**  Vandevender, W. H., (SNLA)
  !***
  ! **Description:**
  !
  !                 **********EISPACK Routines**********
  !
  ! single double complx
  ! ------ ------ ------
  !
  ! RS       -    CH     Computes eigenvalues and, optionally,
  !                      eigenvectors of real symmetric
  !                      (complex Hermitian) matrix.
  !
  ! RSP      -      -    Compute eigenvalues and, optionally,
  !                      eigenvectors of real symmetric matrix
  !                      packed into a one dimensional array.
  !
  ! RG       -    CG     Computes eigenvalues and, optionally,
  !                      eigenvectors of a real (complex) general
  !                      matrix.
  !
  ! BISECT   -      -    Compute eigenvalues of symmetric tridiagonal
  !                      matrix given interval using Sturm sequencing.
  !
  ! IMTQL1   -      -    Computes eigenvalues of symmetric tridiagonal
  !                      matrix implicit QL method.
  !
  ! IMTQL2   -      -    Computes eigenvalues and eigenvectors of
  !                      symmetric tridiagonal matrix using
  !                      implicit QL method.
  !
  ! IMTQLV   -      -    Computes eigenvalues of symmetric tridiagonal
  !                      matrix by the implicit QL method.
  !                      Eigenvectors may be computed later.
  !
  ! RATQR    -      -    Computes largest or smallest eigenvalues
  !                      of symmetric tridiagonal matrix using
  !                      rational QR method with Newton correction.
  !
  ! RST      -      -    Compute eigenvalues and, optionally,
  !                      eigenvectors of real symmetric tridiagonal
  !                      matrix.
  !
  ! RT       -      -    Compute eigenvalues and eigenvectors of
  !                      a special real tridiagonal matrix.
  !
  ! TQL1     -      -    Compute eigenvalues of symmetric tridiagonal
  !                      matrix by QL method.
  !
  ! TQL2     -      -    Compute eigenvalues and eigenvectors
  !                      of symmetric tridiagonal matrix.
  !
  ! TQLRAT   -      -    Computes eigenvalues of symmetric
  !                      tridiagonal matrix a rational variant
  !                      of the QL method.
  !
  ! TRIDIB   -      -    Computes eigenvalues of symmetric
  !                      tridiagonal matrix given interval using
  !                      Sturm sequencing.
  !
  ! TSTURM   -      -    Computes eigenvalues of symmetric tridiagonal
  !                      matrix given interval and eigenvectors
  !                      by Sturm sequencing.  This subroutine
  !                      is a translation of the ALGOL procedure
  !                      TRISTURM by Peters and Wilkinson. HANDBOOK
  !                      FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA,
  !                      418-439(1971).
  !
  ! BQR      -      -    Computes some of the eigenvalues of a real
  !                      symmetric matrix using the QR method with
  !                      shifts of origin.
  !
  ! RSB      -      -    Computes eigenvalues and, optionally,
  !                      eigenvectors of symmetric band matrix.
  !
  ! RSG      -      -    Computes eigenvalues and, optionally,
  !                      eigenvectors of symmetric generalized
  !                      eigenproblem: A*X=(LAMBDA)*B*X
  !
  ! RSGAB    -      -    Computes eigenvalues and, optionally,
  !                      eigenvectors of symmetric generalized
  !                      eigenproblem: A*B*X=(LAMBDA)*X
  !
  ! RSGBA    -      -    Computes eigenvalues and, optionally,
  !                      eigenvectors of symmetric generalized
  !                      eigenproblem: B*A*X=(LAMBDA)*X
  !
  ! RGG      -      -    Computes eigenvalues and eigenvectors
  !                      for real generalized eigenproblem:
  !                      A*X=(LAMBDA)*B*X.
  !
  ! BALANC   -    CBAL   Balances a general real (complex)
  !                      matrix and isolates eigenvalues whenever
  !                      possible.
  !
  ! BANDR    -      -    Reduces real symmetric band matrix
  !                      to symmetric tridiagonal matrix and,
  !                      optionally, accumulates orthogonal similarity
  !                      transformations.
  !
  ! HTRID3   -      -    Reduces complex Hermitian (packed) matrix
  !                      to real symmetric tridiagonal matrix by unitary
  !                      similarity transformations.
  !
  ! HTRIDI   -      -    Reduces complex Hermitian matrix to real
  !                      symmetric tridiagonal matrix using unitary
  !                      similarity transformations.
  !
  ! TRED1    -      -    Reduce real symmetric matrix to symmetric
  !                      tridiagonal matrix using orthogonal
  !                      similarity transformations.
  !
  ! TRED2    -      -    Reduce real symmetric matrix to symmetric
  !                      tridiagonal matrix using and accumulating
  !                      orthogonal transformations.
  !
  ! TRED3    -      -    Reduce  symmetric matrix stored in packed
  !                      form to symmetric tridiagonal matrix using
  !                      orthogonal transformations.
  !
  ! ELMHES   -    COMHES Reduces real (complex) general matrix to
  !                      upper Hessenberg form using stabilized
  !                      elementary similarity transformations.
  !
  ! ORTHES   -    CORTH  Reduces real (complex) general matrix to upper
  !                      Hessenberg form orthogonal (unitary)
  !                      similarity transformations.
  !
  ! QZHES    -      -    The first step of the QZ algorithm for solving
  !                      generalized matrix eigenproblems.  Accepts
  !                      a pair of real general matrices and reduces
  !                      one of them to upper Hessenberg and the other
  !                      to upper triangular form using orthogonal
  !                      transformations. Usually followed by QZIT,
  !                      QZVAL, QZ
  !
  ! QZIT     -      -    The second step of the QZ algorithm for
  !                      generalized eigenproblems.  Accepts an upper
  !                      Hessenberg and an upper triangular matrix
  !                      and reduces the former to quasi-triangular
  !                      form while preserving the form of the latter.
  !                      Usually preceded by QZHES and followed by QZVAL
  !                      and QZVEC.
  !
  ! FIGI     -      -    Transforms certain real non-symmetric
  !                      tridiagonal matrix to symmetric tridiagonal
  !                      matrix.
  !
  ! FIGI2    -      -    Transforms certain real non-symmetric
  !                      tridiagonal matrix to symmetric tridiagonal
  !                      matrix.
  !
  ! REDUC    -      -    Reduces generalized symmetric eigenproblem
  !                      A*X=(LAMBDA)*B*X, to standard symmetric
  !                      eigenproblem using Cholesky factorization.
  !
  ! REDUC2   -      -    Reduces certain generalized symmetric
  !                      eigenproblems standard symmetric eigenproblem,
  !                      using Cholesky factorization.
  !
  !   -      -    COMLR  Computes eigenvalues of a complex upper
  !                      Hessenberg matrix using the modified LR method.
  !
  !   -      -    COMLR2 Computes eigenvalues and eigenvectors of
  !                      complex upper Hessenberg matrix using
  !                      modified LR method.
  !
  ! HQR      -    COMQR  Computes eigenvalues of a real (complex)
  !                      upper Hessenberg matrix using the QR method.
  !
  ! HQR2     -    COMQR2 Computes eigenvalues and eigenvectors of
  !                      real (complex) upper Hessenberg matrix
  !                      using QR method.
  !
  ! INVIT    -    CINVIT Computes eigenvectors of real (complex)
  !                      Hessenberg matrix associated with specified
  !                      eigenvalues by inverse iteration.
  !
  ! QZVAL    -      -    The third step of the QZ algorithm for
  !                      generalized eigenproblems.  Accepts a pair
  !                      of real matrices, one quasi-triangular form
  !                      and the other in upper triangular form and
  !                      computes the eigenvalues of the associated
  !                      eigenproblem.  Usually preceded by QZHES,
  !                      QZIT, and followed by QZVEC.
  !
  ! BANDV    -      -    Forms eigenvectors of real symmetric band
  !                      matrix associated with a set of ordered
  !                      approximate eigenvalue by inverse iteration.
  !
  ! QZVEC    -      -    The optional fourth step of the QZ algorithm
  !                      for generalized eigenproblems.  Accepts
  !                      a matrix in quasi-triangular form and another
  !                      in upper triangular and computes the
  !                      eigenvectors of the triangular problem
  !                      and transforms them back to the original
  !                      coordinates Usually preceded by QZHES, QZIT,
  !                      QZVAL.
  !
  ! TINVIT   -      -    Eigenvectors of symmetric tridiagonal
  !                      matrix corresponding to some specified
  !                      eigenvalues, using inverse iteration.
  !
  ! BAKVEC   -      -    Forms eigenvectors of certain real
  !                      non-symmetric tridiagonal matrix from
  !                      symmetric tridiagonal matrix output from FIGI.
  !
  ! BALBAK   -    CBABK2 Forms eigenvectors of real (complex) general
  !                      matrix from eigenvectors of matrix output
  !                      from BALANC (CBAL).
  !
  ! ELMBAK   -    COMBAK Forms eigenvectors of real (complex) general
  !                      matrix from eigenvectors of upper Hessenberg
  !                      matrix output from ELMHES (COMHES).
  !
  ! ELTRAN   -      -    Accumulates the stabilized elementary
  !                      similarity transformations used in the
  !                      reduction of a real general matrix to upper
  !                      Hessenberg form by ELMHES.
  !
  ! HTRIB3   -      -    Computes eigenvectors of complex Hermitian
  !                      matrix from eigenvectors of real symmetric
  !                      tridiagonal matrix output from HTRID3.
  !
  ! HTRIBK   -      -    Forms eigenvectors of complex Hermitian
  !                      matrix from eigenvectors of real symmetric
  !                      tridiagonal matrix output from HTRIDI.
  !
  ! ORTBAK   -    CORTB  Forms eigenvectors of general real (complex)
  !                      matrix from eigenvectors of upper Hessenberg
  !                      matrix output from ORTHES (CORTH).
  !
  ! ORTRAN   -      -    Accumulates orthogonal similarity
  !                      transformations in reduction of real general
  !                      matrix by ORTHES.
  !
  ! REBAK    -      -    Forms eigenvectors of generalized symmetric
  !                      eigensystem from eigenvectors of derived
  !                      matrix output from REDUC or REDUC2.
  !
  ! REBAKB   -      -    Forms eigenvectors of generalized symmetric
  !                      eigensystem from eigenvectors of derived
  !                      matrix output from REDUC2
  !
  ! TRBAK1   -      -    Forms the eigenvectors of real symmetric
  !                      matrix from eigenvectors of symmetric
  !                      tridiagonal matrix formed by TRED1.
  !
  ! TRBAK3   -      -    Forms eigenvectors of real symmetric matrix
  !                      from the eigenvectors of symmetric tridiagonal
  !                      matrix formed by TRED3.
  !
  ! MINFIT   -      -    Compute Singular Value Decomposition
  !                      of rectangular matrix and solve related
  !                      Linear Least Squares problem.
  !
  !***
  ! **References:**  B. T. Smith, J. M. Boyle, J. J. Dongarra, B. S. Garbow,
  !                 Y. Ikebe, V. C. Klema and C. B. Moler, Matrix Eigen-
  !                 system Routines - EISPACK Guide, Springer-Verlag,
  !                 1976.
  !***
  ! **Routines called:**  (NONE)

  !* REVISION HISTORY  (YYMMDD)
  !   811101  DATE WRITTEN
  !   861211  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900723  PURPOSE section revised.  (WRB)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  
  !* FIRST EXECUTABLE STATEMENT  EISDOC
END SUBROUTINE EISDOC
