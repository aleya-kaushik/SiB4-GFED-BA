        !COMPILER-GENERATED INTERFACE MODULE: Fri Sep 26 16:45:26 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE POOL_AUTO_RESP__genmod
          INTERFACE 
            SUBROUTINE POOL_AUTO_RESP(POOLCONT,CLIM_ASSIM,CLIM_LAI,     &
     &ASSIMD,C13ASSIMD,LAI,TC,TD_LAY,ROOTF_LAY,POOLLT,RESP_SOIL,        &
     &RESP_SOIL_LAY,RESP_SOILC13,RESP_SOILC13_LAY,FRACT)
              USE MODULE_SIB, ONLY :                                    &
     &          SOIL_TYPE,                                              &
     &          POOLL_TYPE,                                             &
     &          FRACT_TYPE
              USE MODULE_PARAM, ONLY :                                  &
     &          POOL_PARAM
              TYPE (POOL_PARAM), INTENT(INOUT) :: POOLCONT
              REAL(KIND=8), INTENT(IN) :: CLIM_ASSIM
              REAL(KIND=8), INTENT(IN) :: CLIM_LAI
              REAL(KIND=8), INTENT(IN) :: ASSIMD
              REAL(KIND=8), INTENT(IN) :: C13ASSIMD
              REAL(KIND=8), INTENT(IN) :: LAI
              REAL(KIND=8), INTENT(IN) :: TC
              REAL(KIND=8), INTENT(IN) :: TD_LAY(10)
              REAL(KIND=8), INTENT(IN) :: ROOTF_LAY(10)
              TYPE (POOLL_TYPE), INTENT(INOUT) :: POOLLT
              REAL(KIND=8), INTENT(INOUT) :: RESP_SOIL
              REAL(KIND=8), INTENT(INOUT) :: RESP_SOIL_LAY(10)
              REAL(KIND=8), INTENT(INOUT) :: RESP_SOILC13
              REAL(KIND=8), INTENT(INOUT) :: RESP_SOILC13_LAY(10)
              TYPE (FRACT_TYPE), INTENT(IN) :: FRACT
            END SUBROUTINE POOL_AUTO_RESP
          END INTERFACE 
        END MODULE POOL_AUTO_RESP__genmod
