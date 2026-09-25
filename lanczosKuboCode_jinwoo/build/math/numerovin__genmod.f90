        !COMPILER-GENERATED INTERFACE MODULE: Thu Jun 13 16:27:35 2019
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE NUMEROVIN__genmod
          INTERFACE 
            FUNCTION NUMEROVIN(F,N,I,H) RESULT(Y)
              INTEGER(KIND=4), INTENT(IN) :: N
              REAL(KIND=8), INTENT(IN) :: F(N)
              INTEGER(KIND=4), INTENT(IN) :: I
              REAL(KIND=8), INTENT(IN) :: H
              REAL(KIND=8) :: Y
            END FUNCTION NUMEROVIN
          END INTERFACE 
        END MODULE NUMEROVIN__genmod
