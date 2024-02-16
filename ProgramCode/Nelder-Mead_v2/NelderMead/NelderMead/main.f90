!***********************************************************
!*   Multidimensional minimization of a function FUNC(X)   *
!*  where X is an NDIM-dimensional vector, by the downhill *
!*  simplex method of Nelder and Mead.                     *
!* ------------------------------------------------------- *
!* SAMPLE RUN: Find a minimum of function F(x,y):          *
!*             F=Sin(R)/R, where R = Sqrt(x*x+y*y).        *
!*                                                         *
!*  Number of iterations:          22                      *
!*                                                         *
!*  Best NDIM+1 points:                                    *
!*  0.41226859E+01  0.17869153E+01                         *
!*  0.41664774E+01  0.16826982E+01                         *
!*  0.41424544E+01  0.17411759E+01                         *
!*                                                         *
!*  Best NDIM+1 mimimum values:                            *
!* -0.21723363E+00                                         *
!* -0.21723363E+00                                         *
!* -0.21723363E+00                                         *
!*                                                         *
!* ------------------------------------------------------- *
!* REFERENCE: "Numerical Recipes, The Art of Scientific    *
!*             Computing by W.H. Press, B.P. Flannery,     *
!*             S.A. Teukolsky and W.T. Vetterling,         *
!*             Cambridge University Press, 1986"           *
!*             [BIBLI 08].                                 *
!*                                                         *
!*               Fortran 90 Release By J-P Moreau, Paris.  *
!*                           (www.jpmoreau.fr)             *
!*********************************************************** 
PROGRAM TEST_AMOEBA
IMPLICIT REAL*8 A-H,O-Z
PARAMETER(MP=21,NP=20)      

DIMENSION P(MP,NP), Y(MP), PT(MP)

  NDIM=2       ! 2 variables
  FTOL=1.D-8   ! Required tolerance

  !define NDIM+1 initial vertices (one by row)
  P(1,1)=1.d0;  P(1,2)=2.d0
  P(2,1)=-2.d0; P(2,2)=-3.d0
  P(3,1)=4.d0;  P(3,2)=2.d0

  !Initialize Y to the values of FUNC evaluated 
  !at the NDIM+1 vertices (rows) of P
  DO I=1, NDIM+1
    PT(1)=P(I,1); PT(2)=P(I,2); 
    Y(I)=FUNC(PT)
  END DO

  !call main subroutine
  CALL AMOEBA(P,Y,MP,NP,NDIM,FTOL,ITER)

  !print results
  print *,' '
  print *,' Number of iterations:', ITER
  print *,' '
  print *,' Best NDIM+1 points:'
  write(*,10) ((P(I,J),J=1,NDIM),I=1,NDIM+1)
  print *,' '
  print *,' Best NDIM+1 mimimum values:'
  write(*,20) (Y(I),I=1,NDIM+1) 
  print *,' '

10 format(2E16.8)
20 format(E16.8)

END

!user defined function to minimize
REAL*8 FUNCTION FUNC(P)
REAL*8 P(2),R
R=DSQRT(P(1)*P(1)+P(2)*P(2))
IF (DABS(R).LT.1.D-12) THEN
  FUNC=1.D0
ELSE
  FUNC=DSIN(R)/R
END IF
RETURN
END


SUBROUTINE AMOEBA(P,Y,MP,NP,NDIM,FTOL,ITER)
!-------------------------------------------------------------------
! Multidimensional minimization of the function FUNC(X) where X is
! an NDIM-dimensional vector, by the downhill simplex method of
! Nelder and Mead. Input is a matrix P whose NDIM+1 rows are NDIM-
! dimensional vectors which are the vertices of the starting simplex
! (Logical dimensions of P are P(NDIM+1,NDIM); physical dimensions
! are input as P(NP,NP)). Also input is the vector Y of length NDIM
! +1, whose components must be pre-initialized to the values of FUNC
! evaluated at the NDIM+1 vertices (rows) of P; and FTOL the fractio-
! nal convergence tolerance to be achieved in the function value. On
! output, P and Y will have been reset to NDIM+1 new points all within
! FTOL of a minimum function value, and ITER gives the number of ite-
! rations taken.
!--------------------------------------------------------------------
IMPLICIT REAL*8 A-H,O-Z
PARAMETER(NMAX=20,ALPHA=1.d0,BETA=0.5d0,GAMMA=2.d0,ITMAX=500)
! Expected maximum number of dimensions, three parameters which define
! the expansions and contractions, and maximum allowed number of
! iterations.
DIMENSION P(MP,NP), Y(MP), PR(NMAX), PRR(NMAX), PBAR(NMAX)
  MPTS=NDIM+1
  ITER=0
1 ILO=1
  IF(Y(1).GT.Y(2)) THEN
    IHI=1
    INHI=2
  ELSE
    IHI=2
    INHI=1
  ENDIF
  DO I=1, MPTS
    IF(Y(I).LT.Y(ILO)) ILO=I
    IF(Y(I).GT.Y(IHI)) THEN
      INHI=IHI
      IHI=I
    ELSE IF (Y(I).GT.Y(INHI)) THEN
      IF(I.NE.IHI) INHI=I
    END IF
  END DO
! Compute the fractional range from highest to lowest and return if
! satisfactory.
  RTOL=2.d0*DABS(Y(IHI)-Y(ILO))/(DABS(Y(IHI))+DABS(Y(ILO)))
  IF(RTOL.LT.FTOL) RETURN
  IF(ITER.EQ.ITMAX) Pause ' Amoeba exceeding maximum iterations.'
  ITER=ITER+1
  DO J=1, NDIM
    PBAR(J)=0.d0
  END DO
  DO I=1, MPTS
    IF(I.NE.IHI) THEN
      DO J=1,NDIM
        PBAR(J)=PBAR(J) + P(I,J)
      END DO
    END IF   
  END DO
  DO J=1, NDIM
    PBAR(J)=PBAR(J)/NDIM
    PR(J)=(1.d0+ALPHA)*PBAR(J) - ALPHA*P(IHI,J)
  END DO
  YPR=FUNC(PR)
  IF(YPR.LE.Y(ILO)) THEN
    DO J=1,NDIM
      PRR(J)=GAMMA*PR(J) + (1.d0-GAMMA)*PBAR(J)
    END DO
    YPRR=FUNC(PRR)
    IF(YPRR.LT.Y(ILO)) THEN
      DO J=1, NDIM
        P(IHI,J)=PRR(J)
      END DO
      Y(IHI)=YPRR
    ELSE
      DO J=1, NDIM
        P(IHI,J)=PR(J)
      END DO
      Y(IHI)=YPR	  
    END IF
  ELSE IF(YPR.GE.Y(INHI)) THEN
    IF(YPR.LT.Y(IHI)) THEN
      DO J=1, NDIM
        P(IHI,J)=PR(J)
      END DO
      Y(IHI)=YPR
    END IF
    DO J=1, NDIM
      PRR(J)=BETA*P(IHI,J) + (1.d0-BETA)*PBAR(J)
    END DO
    YPRR=FUNC(PRR)
    IF(YPRR.LT.Y(IHI)) THEN
      DO J=1, NDIM
        P(IHI,J)=PRR(J)
      END DO
      Y(IHI)=YPRR
    ELSE
      DO I=1, MPTS
        IF(I.NE.ILO) THEN
          DO J=1,NDIM
            PR(J)=0.5d0*(P(I,J) + P(ILO,J))
	    P(I,J)=PR(J)
	  END DO
          Y(I)=FUNC(PR)
	END IF
      END DO
    END IF
  ELSE
    DO J=1, NDIM
      P(IHI,J)=PR(J)
    END DO
    Y(IHI)=YPR
  END IF
  GO TO 1
  END
  
!end of file tamoeba.f90