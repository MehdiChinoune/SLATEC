!*==HPPERM.f90  processed by SPAG 6.72Dc at 11:01 on  6 Feb 2019
!DECK HPPERM
SUBROUTINE HPPERM(Hx,N,Iperm,Work,Ier)
  IMPLICIT NONE
  !*--HPPERM5
  !***BEGIN PROLOGUE  HPPERM
  !***PURPOSE  Rearrange a given array according to a prescribed
  !            permutation vector.
  !***LIBRARY   SLATEC
  !***CATEGORY  N8
  !***TYPE      CHARACTER (SPPERM-S, DPPERM-D, IPPERM-I, HPPERM-H)
  !***KEYWORDS  APPLICATION OF PERMUTATION TO DATA VECTOR
  !***AUTHOR  McClain, M. A., (NIST)
  !           Rhoads, G. S., (NBS)
  !***DESCRIPTION
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
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  XERMSG
  !***REVISION HISTORY  (YYMMDD)
  !   901004  DATE WRITTEN
  !   920507  Modified by M. McClain to revise prologue text and to add
  !           check for length of work array.
  !***END PROLOGUE  HPPERM
  INTEGER N, Iperm(*), i, Ier, indx, indx0, istrt
  CHARACTER*(*) Hx(*), Work
  !***FIRST EXECUTABLE STATEMENT  HPPERM
  Ier = 0
  IF ( N<1 ) THEN
    Ier = 1
    CALL XERMSG('SLATEC','HPPERM',&
      'The number of values to be rearranged, N, is not positive.'&
      ,Ier,1)
    RETURN
  ENDIF
  IF ( LEN(Work)<LEN(Hx(1)) ) THEN
    Ier = 2
    CALL XERMSG('SLATEC','HPPERM',&
      'The length of the work variable, WORK, is too short.',Ier,&
      1)
    RETURN
  ENDIF
  !
  !     CHECK WHETHER IPERM IS A VALID PERMUTATION
  !
  DO i = 1, N
    indx = ABS(Iperm(i))
    IF ( (indx>=1).AND.(indx<=N) ) THEN
      IF ( Iperm(indx)>0 ) THEN
        Iperm(indx) = -Iperm(indx)
        CYCLE
      ENDIF
    ENDIF
    Ier = 3
    CALL XERMSG('SLATEC','HPPERM',&
      'The permutation vector, IPERM, is not valid.',Ier,1)
    RETURN
  ENDDO
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
      ENDDO
      Hx(indx0) = Work
    ENDIF
  ENDDO
  !
END SUBROUTINE HPPERM
