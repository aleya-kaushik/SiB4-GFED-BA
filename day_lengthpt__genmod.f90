        !COMPILER-GENERATED INTERFACE MODULE: Fri Sep 26 16:45:33 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE DAY_LENGTHPT__genmod
          INTERFACE 
            SUBROUTINE DAY_LENGTHPT(PI180,TAN_DEC,LAT,DLENGTH,DLENGTHDT)
              REAL(KIND=8), INTENT(IN) :: PI180
              REAL(KIND=8), INTENT(IN) :: TAN_DEC
              REAL(KIND=4), INTENT(IN) :: LAT
              REAL(KIND=8), INTENT(INOUT) :: DLENGTH
              REAL(KIND=8), INTENT(INOUT) :: DLENGTHDT
            END SUBROUTINE DAY_LENGTHPT
          END INTERFACE 
        END MODULE DAY_LENGTHPT__genmod
