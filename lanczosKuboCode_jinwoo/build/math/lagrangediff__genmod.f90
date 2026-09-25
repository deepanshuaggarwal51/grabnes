        !COMPILER-GENERATED INTERFACE MODULE: Thu Jun 13 16:27:34 2019
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE LAGRANGEDIFF__genmod
          INTERFACE 
            SUBROUTINE LAGRANGEDIFF(F,DF,K)
              INTEGER(KIND=4), INTENT(IN) :: K
              REAL(KIND=8), INTENT(IN) :: F(K+1)
              REAL(KIND=8), INTENT(OUT) :: DF(K+1)
            END SUBROUTINE LAGRANGEDIFF
          END INTERFACE 
        END MODULE LAGRANGEDIFF__genmod
