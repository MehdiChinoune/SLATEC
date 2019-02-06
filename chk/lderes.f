*DECK LDERES
      LOGICAL FUNCTION LDERES (TYPE, UPLO, M, N, AA, AS, LDA)
C***BEGIN PROLOGUE  LDERES
C***SUBSIDIARY
C***PURPOSE  Test if selected elements in two arrays are equal.
C***LIBRARY   SLATEC (BLAS)
C***AUTHOR  Du Croz, J. J., (NAG)
C           Hanson, R. J., (SNLA)
C***DESCRIPTION
C
C  Tests if selected elements in two arrays are equal.
C
C  TYPE is 'GE', 'SY' or 'SP'.
C
C  Auxiliary routine for test program for Level 2 Blas.
C***REFERENCES  (NONE)
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   870810  DATE WRITTEN
C   910620  Modified to meet SLATEC code and prologue standards.  (BKS)
C***END PROLOGUE  LDERES
C     .. Scalar Arguments ..
      INTEGER            LDA, M, N
      CHARACTER*1        UPLO
      CHARACTER*2        TYPE
C     .. Array Arguments ..
      DOUBLE PRECISION   AA( LDA, * ), AS( LDA, * )
C     .. Local Scalars ..
      INTEGER            I, IBEG, IEND, J
      LOGICAL            UPPER
C***FIRST EXECUTABLE STATEMENT  LDERES
      UPPER = UPLO.EQ.'U'
      IF( TYPE.EQ.'GE' )THEN
         DO 20 J = 1, N
            DO 10 I = M + 1, LDA
               IF( AA( I, J ).NE.AS( I, J ) )
     $            GO TO 70
   10       CONTINUE
   20    CONTINUE
      ELSE IF( TYPE.EQ.'SY' )THEN
         DO 50 J = 1, N
            IF( UPPER )THEN
               IBEG = 1
               IEND = J
            ELSE
               IBEG = J
               IEND = N
            END IF
            DO 30 I = 1, IBEG - 1
               IF( AA( I, J ).NE.AS( I, J ) )
     $            GO TO 70
   30       CONTINUE
            DO 40 I = IEND + 1, LDA
               IF( AA( I, J ).NE.AS( I, J ) )
     $            GO TO 70
   40       CONTINUE
   50    CONTINUE
      END IF
C
      LDERES = .TRUE.
      GO TO 80
   70 CONTINUE
      LDERES = .FALSE.
   80 RETURN
C
C     End of LDERES.
C
      END
