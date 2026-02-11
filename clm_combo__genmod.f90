        !COMPILER-GENERATED INTERFACE MODULE: Fri Sep 26 16:45:32 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE CLM_COMBO__genmod
          INTERFACE 
            SUBROUTINE CLM_COMBO(DZ1,LIQ1,ICE1,TEMP1,DZ2,LIQ2,ICE2,TEMP2&
     &)
              REAL(KIND=8), INTENT(INOUT) :: DZ1
              REAL(KIND=8), INTENT(INOUT) :: LIQ1
              REAL(KIND=8), INTENT(INOUT) :: ICE1
              REAL(KIND=8), INTENT(INOUT) :: TEMP1
              REAL(KIND=8), INTENT(IN) :: DZ2
              REAL(KIND=8), INTENT(IN) :: LIQ2
              REAL(KIND=8), INTENT(IN) :: ICE2
              REAL(KIND=8), INTENT(IN) :: TEMP2
            END SUBROUTINE CLM_COMBO
          END INTERFACE 
        END MODULE CLM_COMBO__genmod
