!** DPPERM
SUBROUTINE DPPERM(Dx,N,Iperm,Ier)
  !>
  !  Rearrange a given array according to a prescribed
  !            permutation vector.
  !***
  ! **Library:**   SLATEC
  !***
  ! **Category:**  N8
  !***
  ! **Type:**      DOUBLE PRECISION (SPPERM-S, DPPERM-D, IPPERM-I, HPPERM-H)
  !***
  ! **Keywords:**  PERMUTATION, REARRANGEMENT
  !***
  ! **Author:**  McClain, M. A., (NIST)
  !           Rhoads, G. S., (NBS)
  !***
  ! **Description:**
  !
  !         DPPERM rearranges the data vector DX according to the
  !         permutation IPERM: DX(I) <--- DX(IPERM(I)).  IPERM could come
  !         from one of the sorting routines IPSORT, SPSORT, DPSORT or
  !         HPSORT.
  !
  !     Description of Parameters
  !         DX - input/output -- double precision array of values to be
  !                   rearranged.
  !         N - input -- number of values in double precision array DX.
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
  USE service, ONLY : XERMSG
  INTEGER N, Iperm(*), i, Ier, indx, indx0, istrt
  REAL(8) :: Dx(*), dtemp
  !* FIRST EXECUTABLE STATEMENT  DPPERM
  Ier = 0
  IF ( N<1 ) THEN
    Ier = 1
    CALL XERMSG('SLATEC','DPPERM',&
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
    CALL XERMSG('SLATEC','DPPERM',&
      'The permutation vector, IPERM, is not valid.',Ier,1)
    RETURN
  END DO
  !
  !     REARRANGE THE VALUES OF DX
  !
  !     USE THE IPERM VECTOR AS A FLAG.
  !     IF IPERM(I) > 0, THEN THE I-TH VALUE IS IN CORRECT LOCATION
  !
  DO istrt = 1, N
    IF ( Iperm(istrt)<=0 ) THEN
      indx = istrt
      indx0 = indx
      dtemp = Dx(istrt)
      DO WHILE ( Iperm(indx)<0 )
        Dx(indx) = Dx(-Iperm(indx))
        indx0 = indx
        Iperm(indx) = -Iperm(indx)
        indx = Iperm(indx)
      END DO
      Dx(indx0) = dtemp
    END IF
  END DO
  !
END SUBROUTINE DPPERM
