!** DAVINT
SUBROUTINE DAVINT(X,Y,N,Xlo,Xup,Ans,Ierr)
  !> Integrate a function tabulated at arbitrarily spaced
  !            abscissas using overlapping parabolas.
  !***
  ! **Library:**   SLATEC
  !***
  ! **Category:**  H2A1B2
  !***
  ! **Type:**      DOUBLE PRECISION (AVINT-S, DAVINT-D)
  !***
  ! **Keywords:**  INTEGRATION, QUADRATURE, TABULATED DATA
  !***
  ! **Author:**  Jones, R. E., (SNLA)
  !***
  ! **Description:**
  !
  !     Abstract
  !         DAVINT integrates a function tabulated at arbitrarily spaced
  !         abscissas.  The limits of integration need not coincide
  !         with the tabulated abscissas.
  !
  !         A method of overlapping parabolas fitted to the data is used
  !         provided that there are at least 3 abscissas between the
  !         limits of integration.  DAVINT also handles two special cases.
  !         If the limits of integration are equal, DAVINT returns a
  !         result of zero regardless of the number of tabulated values.
  !         If there are only two function values, DAVINT uses the
  !         trapezoid rule.
  !
  !     Description of Parameters
  !         The user must dimension all arrays appearing in the call list
  !              X(N), Y(N)
  !
  !         Input--
  !      X    - DOUBLE PRECISION array of abscissas, which must be in
  !             increasing order.
  !      Y    - DOUBLE PRECISION array of function values. i.e.,
  !                Y(I)=FUNC(X(I))
  !      N    - The integer number of function values supplied.
  !                N >= 2 unless XLO = XUP.
  !      XLO  - DOUBLE PRECISION lower limit of integration
  !      XUP  - DOUBLE PRECISION upper limit of integration.  Must have
  !              XLO<=XUP
  !
  !         Output--
  !      ANS  - Double Precision computed approximate value of integral
  !      IERR - A status code
  !           --Normal Code
  !                =1 Means the requested integration was performed.
  !           --Abnormal Codes
  !                =2 Means XUP was less than XLO.
  !                =3 Means the number of X(I) between XLO and XUP
  !                   (inclusive) was less than 3 and neither of the two
  !                   special cases described in the abstract occurred.
  !                   No integration was performed.
  !                =4 Means the restriction X(I+1)>X(I) was violated.
  !                =5 Means the number N of function values was < 2.
  !                   ANS is set to zero if IERR=2,3,4,or 5.
  !
  !    DAVINT is documented completely in SC-M-69-335
  !    Original program from *Numerical Integration* by Davis & Rabinowitz
  !    Adaptation and modifications by Rondall E Jones.
  !
  !***
  ! **References:**  R. E. Jones, Approximate integrator of functions
  !                 tabulated at arbitrarily spaced abscissas,
  !                 Report SC-M-69-335, Sandia Laboratories, 1969.
  !***
  ! **Routines called:**  XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   690901  DATE WRITTEN
  !   890831  Modified array declarations.  (WRB)
  !   890831  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  USE service, ONLY : XERMSG
  !
  INTEGER :: Ierr, N
  REAL(DP) :: Ans, Xlo, Xup, X(N), Y(N)
  INTEGER :: i, inlft, inrt, istart, istop
  REAL(DP) :: a, b, c, ca, cb, cc, fl, fr, r3, rp5, slope, summ, syl, syl2, syl3, &
    syu, syu2, syu3, term1, term2, term3, x1, x12, x13, x2, x23, x3
  !     BEGIN BLOCK PERMITTING ...EXITS TO 190
  !        BEGIN BLOCK PERMITTING ...EXITS TO 180
  !* FIRST EXECUTABLE STATEMENT  DAVINT
  Ierr = 1
  Ans = 0.0D0
  IF( Xlo>Xup ) THEN
    Ierr = 2
    !     ......EXIT
    CALL XERMSG('DAVINT',&
      'THE UPPER LIMIT OF INTEGRATION WAS NOT GREATER THAN THE LOWER LIMIT.',4,1)
  ELSEIF( Xlo/=Xup ) THEN
    IF( N>=2 ) THEN
      DO i = 2, N
        !        ............EXIT
        IF( X(i)<=X(i-1) ) GOTO 50
        !                 ...EXIT
        IF( X(i)>Xup ) EXIT
      END DO
      IF( N<3 ) THEN
        !
        !                    SPECIAL N=2 CASE
        slope = (Y(2)-Y(1))/(X(2)-X(1))
        fl = Y(1) + slope*(Xlo-X(1))
        fr = Y(2) + slope*(Xup-X(2))
        Ans = 0.5D0*(fl+fr)*(Xup-Xlo)
        !     ...............EXIT
        RETURN
      ELSEIF( X(N-2)<Xlo ) THEN
        Ierr = 3
        CALL XERMSG('DAVINT',&
          'THERE WERE LESS THAN THREE FUNCTION VALUES BETWEEN THE LIMITS OF INTEGRATION.',4,1)
        !     ...............EXIT
        RETURN
      ELSEIF( X(3)<=Xup ) THEN
        i = 1
        DO WHILE( X(i)<Xlo )
          i = i + 1
        END DO
        inlft = i
        i = N
        DO WHILE( X(i)>Xup )
          i = i - 1
        END DO
        inrt = i
        IF( (inrt-inlft)>=2 ) THEN
          istart = inlft
          IF( inlft==1 ) istart = 2
          istop = inrt
          IF( inrt==N ) istop = N - 1
          !
          r3 = 3.0D0
          rp5 = 0.5D0
          summ = 0.0D0
          syl = Xlo
          syl2 = syl*syl
          syl3 = syl2*syl
          !
          DO i = istart, istop
            x1 = X(i-1)
            x2 = X(i)
            x3 = X(i+1)
            x12 = x1 - x2
            x13 = x1 - x3
            x23 = x2 - x3
            term1 = Y(i-1)/(x12*x13)
            term2 = -Y(i)/(x12*x23)
            term3 = Y(i+1)/(x13*x23)
            a = term1 + term2 + term3
            b = -(x2+x3)*term1 - (x1+x3)*term2 - (x1+x2)*term3
            c = x2*x3*term1 + x1*x3*term2 + x1*x2*term3
            IF( i>istart ) THEN
              ca = 0.5D0*(a+ca)
              cb = 0.5D0*(b+cb)
              cc = 0.5D0*(c+cc)
            ELSE
              ca = a
              cb = b
              cc = c
            END IF
            syu = x2
            syu2 = syu*syu
            syu3 = syu2*syu
            summ = summ + ca*(syu3-syl3)/r3 + cb*rp5*(syu2-syl2) + cc*(syu-syl)
            ca = a
            cb = b
            cc = c
            syl = syu
            syl2 = syu2
            syl3 = syu3
          END DO
          syu = Xup
          Ans = summ + ca*(syu**3-syl3)/r3 + cb*rp5*(syu**2-syl2) + cc*(syu-syl)
        ELSE
          Ierr = 3
          !     ...............EXIT
          CALL XERMSG('DAVINT',&
            'THERE WERE LESS THAN THREE FUNCTION VALUES BETWEEN THE LIMITS OF INTEGRATION.',4,1)
        END IF
        RETURN
      ELSE
        Ierr = 3
        CALL XERMSG('DAVINT',&
          'THERE WERE LESS THAN THREE FUNCTION VALUES BETWEEN THE LIMITS OF INTEGRATION.',4,1)
        !     ...............EXIT
        RETURN
      END IF
    ELSE
      Ierr = 5
      CALL XERMSG('DAVINT',&
        'LESS THAN TWO FUNCTION VALUES WERE SUPPLIED.',4,1)
      !     ...............EXIT
      RETURN
    END IF
    50  Ierr = 4
    CALL XERMSG('DAVINT',&
      'THE ABSCISSAS WERE NOT STRICTLY INCREASING.  MUST HAVE X(I-1) < X(I) FOR ALL I.',4,1)
  END IF
  RETURN
END SUBROUTINE DAVINT
