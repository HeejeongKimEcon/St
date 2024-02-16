        !COMPILER-GENERATED INTERFACE MODULE: Wed Sep  8 12:09:59 2021
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE AMOEBA__genmod
          INTERFACE 
            SUBROUTINE AMOEBA(P,Y,MP,NP,NDIM,FTOL,FUNK,ITER)
              INTEGER(KIND=4) :: NP
              INTEGER(KIND=4) :: MP
              REAL(KIND=8) :: P(MP,NP)
              REAL(KIND=8) :: Y(MP)
              INTEGER(KIND=4) :: NDIM
              REAL(KIND=8) :: FTOL
              REAL(KIND=8) :: FUNK
              INTEGER(KIND=4) :: ITER
            END SUBROUTINE AMOEBA
          END INTERFACE 
        END MODULE AMOEBA__genmod
