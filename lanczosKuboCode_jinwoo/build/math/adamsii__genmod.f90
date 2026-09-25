        !COMPILER-GENERATED INTERFACE MODULE: Thu Jun 13 16:27:33 2019
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE ADAMSII__genmod
          INTERFACE 
            FUNCTION ADAMSII(DF,N,I,K) RESULT(F)
              INTEGER(KIND=4), INTENT(IN) :: N
              REAL(KIND=8), INTENT(IN) :: DF(N)
              INTEGER(KIND=4), INTENT(IN) :: I
              INTEGER(KIND=4), INTENT(IN) :: K
              REAL(KIND=8) :: F
            END FUNCTION ADAMSII
          END INTERFACE 
        END MODULE ADAMSII__genmod
