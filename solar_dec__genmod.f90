        !COMPILER-GENERATED INTERFACE MODULE: Fri Sep 26 16:45:33 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE SOLAR_DEC__genmod
          INTERFACE 
            SUBROUTINE SOLAR_DEC(DOY,CUREQNX,LONEARTH,SIN_DEC,COS_DEC,  &
     &TAN_DEC)
              INTEGER(KIND=4), INTENT(IN) :: DOY
              INTEGER(KIND=4), INTENT(IN) :: CUREQNX
              REAL(KIND=8), INTENT(INOUT) :: LONEARTH
              REAL(KIND=8), INTENT(INOUT) :: SIN_DEC
              REAL(KIND=8), INTENT(INOUT) :: COS_DEC
              REAL(KIND=8), INTENT(INOUT) :: TAN_DEC
            END SUBROUTINE SOLAR_DEC
          END INTERFACE 
        END MODULE SOLAR_DEC__genmod
