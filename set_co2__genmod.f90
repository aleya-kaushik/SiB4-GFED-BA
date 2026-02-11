        !COMPILER-GENERATED INTERFACE MODULE: Fri Sep 26 16:45:33 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE SET_CO2__genmod
          INTERFACE 
            SUBROUTINE SET_CO2(YEAR,MON,GPROGT)
              USE MODULE_SIBCONST, ONLY :                               &
     &          SUBCOUNT,                                               &
     &          VARCO2_SWITCH
              USE MODULE_SIB, ONLY :                                    &
     &          SIB,                                                    &
     &          GPROG_TYPE
              INTEGER(KIND=4) :: YEAR
              INTEGER(KIND=4) :: MON
              TYPE (GPROG_TYPE), INTENT(INOUT) :: GPROGT
            END SUBROUTINE SET_CO2
          END INTERFACE 
        END MODULE SET_CO2__genmod
