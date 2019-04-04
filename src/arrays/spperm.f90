!** SPPERM
SUBROUTINE SPPERM(X,N,Iperm,Ier)
  IMPLICIT NONE
  !>
  !***
  !  Rearrange a given array according to a prescribed
  !            permutation vector.
  !***
  ! **Library:**   SLATEC
  !***
  ! **Category:**  N8
  !***
  ! **Type:**      SINGLE PRECISION (SPPERM-S, DPPERM-D, IPPERM-I, HPPERM-H)
  !***
  ! **Keywords:**  APPLICATION OF PERMUTATION TO DATA VECTOR
  !***
  ! **Author:**  McClain, M. A., (NIST)
  !           Rhoads, G. S., (NBS)
  !***
  ! **Description:**
  !
  !         SPPERM rearranges the data vector X according to the
  !         permutation IPERM: X(I) <--- X(IPERM(I)).  IPERM could come
  !         from one of the sorting routines IPSORT, SPSORT, DPSORT or
  !         HPSORT.
  !
  !     Description of Parameters
  !         X - input/output -- real array of values to be rearranged.
  !         N - input -- number of values in real array X.
  !         IPERM - input -- permutation vector.
  !         IER - output -- error indicator:
  !             =  0  if no error,
  !             =  1  if N is zero or negative,
  !             =  2  if IPERM is not a valid permutation.
  !
  !***
  ! **References:**  (NONE)
  !***
  ! **Routines called:**  XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   901004  DATE WRITTEN
  !   920507  Modified by M. McClain to revise prologue text.

  INTEGER N, Iperm(*), i, Ier, indx, indx0, istrt
  REAL X(*), temp
  !* FIRST EXECUTABLE STATEMENT  SPPERM
  Ier = 0
  IF ( N<1 ) THEN
    Ier = 1
    CALL XERMSG('SLATEC','SPPERM',&
      'The number of values to be rearranged, N, is not positive.',Ier,1)
    RETURN
  END IF
  !
  !     CHECK WHETHER IPERM IS A VALID PERMUTATION
  !
  DO i = 1, N
    indx = ABS(Iperm(i))
    IF ( (indx>=1).AND.(indx<=N) ) THEN
      IF ( Iperm(indx)>0 ) THEN
        Iperm(indx) = -Iperm(indx)
        CYCLE
      END IF
    END IF
    Ier = 2
    CALL XERMSG('SLATEC','SPPERM',&
      'The permutation vector, IPERM, is not valid.',Ier,1)
    RETURN
  END DO
  !
  !     REARRANGE THE VALUES OF X
  !
  !     USE THE IPERM VECTOR AS A FLAG.
  !     IF IPERM(I) > 0, THEN THE I-TH VALUE IS IN CORRECT LOCATION
  !
  DO istrt = 1, N
    IF ( Iperm(istrt)<=0 ) THEN
      indx = istrt
      indx0 = indx
      temp = X(istrt)
      DO WHILE ( Iperm(indx)<0 )
        X(indx) = X(-Iperm(indx))
        indx0 = indx
        Iperm(indx) = -Iperm(indx)
        indx = Iperm(indx)
      END DO
      X(indx0) = temp
    END IF
  END DO
  !
END SUBROUTINE SPPERM
