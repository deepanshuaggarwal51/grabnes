        !COMPILER-GENERATED INTERFACE MODULE: Thu Jun 13 16:27:35 2019
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE CROSSPROD__genmod
          INTERFACE 
            FUNCTION CROSSPROD(A,B) RESULT(AXB)
              REAL(KIND=8), INTENT(IN) :: A(3)
              REAL(KIND=8), INTENT(IN) :: B(3)
              REAL(KIND=8) :: AXB(3)
            END FUNCTION CROSSPROD
          END INTERFACE 
        END MODULE CROSSPROD__genmod
