        !COMPILER-GENERATED INTERFACE MODULE: Thu Jun 13 16:27:35 2019
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE TRAPEZOIDALINT__genmod
          INTERFACE 
            FUNCTION TRAPEZOIDALINT(F,N,H,K) RESULT(S)
              INTEGER(KIND=4), INTENT(IN) :: N
              REAL(KIND=8), INTENT(IN) :: F(N)
              REAL(KIND=8), INTENT(IN) :: H
              INTEGER(KIND=4), INTENT(IN) :: K
              REAL(KIND=8) :: S
            END FUNCTION TRAPEZOIDALINT
          END INTERFACE 
        END MODULE TRAPEZOIDALINT__genmod
