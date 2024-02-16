        !COMPILER-GENERATED INTERFACE MODULE: Wed Sep  8 12:12:40 2021
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE AMOTRY__genmod
          INTERFACE 
            FUNCTION AMOTRY(P,Y,PSUM,MP,NP,NDIM,FUNK,IHI,FAC)
              INTEGER(KIND=4) :: NP
              INTEGER(KIND=4) :: MP
              REAL(KIND=8) :: P(MP,NP)
              REAL(KIND=8) :: Y(MP)
              REAL(KIND=8) :: PSUM(NP)
              INTEGER(KIND=4) :: NDIM
              REAL(KIND=8) :: FUNK
              EXTERNAL FUNK
              INTEGER(KIND=4) :: IHI
              REAL(KIND=8) :: FAC
              REAL(KIND=8) :: AMOTRY
            END FUNCTION AMOTRY
          END INTERFACE 
        END MODULE AMOTRY__genmod
