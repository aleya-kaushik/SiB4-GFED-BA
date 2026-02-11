        !COMPILER-GENERATED INTERFACE MODULE: Fri Sep 26 16:45:17 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE CFRAX_CALC__genmod
          INTERFACE 
            SUBROUTINE CFRAX_CALC(PHYSCONT,GPROGT,CO2T,FRACT,FLUXT,CO2M,&
     &CO2GAMMA,CO2_ASSIM)
              USE MODULE_SIB, ONLY :                                    &
     &          GPROG_TYPE,                                             &
     &          FRACT_TYPE,                                             &
     &          CO2_TYPE,                                               &
     &          FLUX_TYPE
              USE MODULE_PARAM, ONLY :                                  &
     &          PHYS_PARAM
              TYPE (PHYS_PARAM), INTENT(IN) :: PHYSCONT
              TYPE (GPROG_TYPE), INTENT(IN) :: GPROGT
              TYPE (CO2_TYPE), INTENT(IN) :: CO2T
              TYPE (FRACT_TYPE), INTENT(INOUT) :: FRACT
              TYPE (FLUX_TYPE), INTENT(IN) :: FLUXT
              REAL(KIND=8), INTENT(IN) :: CO2M
              REAL(KIND=8), INTENT(IN) :: CO2GAMMA
              REAL(KIND=8), INTENT(IN) :: CO2_ASSIM
            END SUBROUTINE CFRAX_CALC
          END INTERFACE 
        END MODULE CFRAX_CALC__genmod
