        !COMPILER-GENERATED INTERFACE MODULE: Thu Jun 13 16:27:33 2019
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE ADAMSIE__genmod
          INTERFACE 
            FUNCTION ADAMSIE(DF,N,I,K) RESULT(F)
              INTEGER(KIND=4), INTENT(IN) :: N
              REAL(KIND=8), INTENT(IN) :: DF(N)
              INTEGER(KIND=4), INTENT(IN) :: I
              INTEGER(KIND=4), INTENT(IN) :: K
              REAL(KIND=8) :: F
            END FUNCTION ADAMSIE
          END INTERFACE 
        END MODULE ADAMSIE__genmod
