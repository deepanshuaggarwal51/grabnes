        !COMPILER-GENERATED INTERFACE MODULE: Thu Jun 13 16:27:33 2019
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE ADAMSIN__genmod
          INTERFACE 
            SUBROUTINE ADAMSIN(F,DF,D2F,C1,N,MN,H,K)
              INTEGER(KIND=4), INTENT(IN) :: N
              REAL(KIND=8), INTENT(INOUT) :: F(N)
              REAL(KIND=8), INTENT(INOUT) :: DF(N)
              REAL(KIND=8), INTENT(INOUT) :: D2F(N)
              REAL(KIND=8), INTENT(IN) :: C1(N)
              INTEGER(KIND=4), INTENT(IN) :: MN
              REAL(KIND=8), INTENT(IN) :: H
              INTEGER(KIND=4), INTENT(IN) :: K
            END SUBROUTINE ADAMSIN
          END INTERFACE 
        END MODULE ADAMSIN__genmod
