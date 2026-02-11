        !COMPILER-GENERATED INTERFACE MODULE: Fri Sep 26 16:45:25 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE PHOSORT__genmod
          INTERFACE 
            SUBROUTINE PHOSORT(NUMIC,IC,GAMMA,RANGE,EYY,PCO2Y)
              INTEGER(KIND=4), INTENT(IN) :: NUMIC
              INTEGER(KIND=4), INTENT(IN) :: IC
              REAL(KIND=8), INTENT(IN) :: GAMMA
              REAL(KIND=8), INTENT(IN) :: RANGE
              REAL(KIND=8), INTENT(INOUT) :: EYY(NUMIC)
              REAL(KIND=8), INTENT(INOUT) :: PCO2Y(NUMIC)
            END SUBROUTINE PHOSORT
          END INTERFACE 
        END MODULE PHOSORT__genmod
