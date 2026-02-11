        !COMPILER-GENERATED INTERFACE MODULE: Fri Sep 26 16:45:33 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE ZENITH_ANGLE__genmod
          INTERFACE 
            SUBROUTINE ZENITH_ANGLE(LATSIB,COS_DEC,SIN_DEC,LOCTIME,COSZ)
              REAL(KIND=4), INTENT(IN) :: LATSIB
              REAL(KIND=8), INTENT(IN) :: COS_DEC
              REAL(KIND=8), INTENT(IN) :: SIN_DEC
              REAL(KIND=8), INTENT(IN) :: LOCTIME
              REAL(KIND=8), INTENT(INOUT) :: COSZ
            END SUBROUTINE ZENITH_ANGLE
          END INTERFACE 
        END MODULE ZENITH_ANGLE__genmod
