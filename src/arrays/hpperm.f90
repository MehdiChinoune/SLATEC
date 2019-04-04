!** HPPERM
SUBROUTINE HPPERM(Hx,N,Iperm,Work,Ier)
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
  ! **Type:**      CHARACTER (SPPERM-S, DPPERM-D, IPPERM-I, HPPERM-H)
  !***
  ! **Keywords:**  APPLICATION OF PERMUTATION TO DATA VECTOR
  !***
  ! **Author:**  McClain, M. A., (NIST)
  !           Rhoads, G. S., (NBS)
  !***
  ! **Description:**
  !
  !         HPPERM rearranges the data vector HX according to the
  !         permutation IPERM: HX(I) <--- HX(IPERM(I)).  IPERM could come
  !         from one of the sorting routines IPSORT, SPSORT, DPSORT or
  !         HPSORT.
  !
  !     Description of Parameters
  !         HX - input/output -- character array of values to be
  !                 rearranged.
  !         N - input -- number of values in character array HX.
  !         IPERM - input -- permutation vector.
  !         WORK - character variable which must have a length
  !                   specification at least as great as that of HX.
  !         IER - output -- error indicator:
  !             =  0  if no error,
  !             =  1  if N is zero or negative,
  !             =  2  if work array is not long enough,
  !             =  3  if IPERM is not a valid permutation.
  !
  !***
  ! **References:**  (NONE)
  !***
  ! **Routines called:**  XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   901004  DATE WRITTEN
  !   920507  Modified by M. McClain to revise prologue text and to add
  !           check for length of work array.

  INTEGER N, Iperm(*), i, Ier, indx, indx0, istrt
  CHARACTER*(*) Hx(*), Work
  !* FIRST EXECUTABLE STATEMENT  HPPERM
  Ier = 0
  IF ( N<1 ) THEN
    Ier = 1
    CALL XERMSG('SLATEC','HPPERM',&
      'The number of values to be rearranged, N, is not positive.',Ier,1)
    RETURN
  END IF
  IF ( LEN(Work)<LEN(Hx(1)) ) THEN
    Ier = 2
    CALL XERMSG('SLATEC','HPPERM',&
      'The length of the work variable, WORK, is too short.',Ier,1)
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
    Ier = 3
    CALL XERMSG('SLATEC','HPPERM',&
      'The permutation vector, IPERM, is not valid.',Ier,1)
    RETURN
  END DO
  !
  !     REARRANGE THE VALUES OF HX
  !
  !     USE THE IPERM VECTOR AS A FLAG.
  !     IF IPERM(I) > 0, THEN THE I-TH VALUE IS IN CORRECT LOCATION
  !
  DO istrt = 1, N
    IF ( Iperm(istrt)<=0 ) THEN
      indx = istrt
      indx0 = indx
      Work = Hx(istrt)
      DO WHILE ( Iperm(indx)<0 )
        Hx(indx) = Hx(-Iperm(indx))
        indx0 = indx
        Iperm(indx) = -Iperm(indx)
        indx = Iperm(indx)
      END DO
      Hx(indx0) = Work
    END IF
  END DO
  !
END SUBROUTINE HPPERM
