        !COMPILER-GENERATED INTERFACE MODULE: Fri Sep 26 16:46:04 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE LOCAL_SET__genmod
          INTERFACE 
            SUBROUTINE LOCAL_SET(PS,TC,CAPACC_LIQ,CAPACC_SNOW,CAPACG,   &
     &SNOW_GMASS,NSL,WWW_LIQ,WWW_ICE,TD)
              REAL(KIND=8), INTENT(IN) :: PS
              REAL(KIND=8), INTENT(IN) :: TC
              REAL(KIND=8), INTENT(IN) :: CAPACC_LIQ
              REAL(KIND=8), INTENT(IN) :: CAPACC_SNOW
              REAL(KIND=8), INTENT(IN) :: CAPACG
              REAL(KIND=8), INTENT(IN) :: SNOW_GMASS
              INTEGER(KIND=1), INTENT(IN) :: NSL
              REAL(KIND=8), INTENT(IN) :: WWW_LIQ(-4:10)
              REAL(KIND=8), INTENT(IN) :: WWW_ICE(-4:10)
              REAL(KIND=8), INTENT(IN) :: TD(-4:10)
            END SUBROUTINE LOCAL_SET
          END INTERFACE 
        END MODULE LOCAL_SET__genmod
