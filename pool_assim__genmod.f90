        !COMPILER-GENERATED INTERFACE MODULE: Fri Sep 26 16:45:26 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE POOL_ASSIM__genmod
          INTERFACE 
            SUBROUTINE POOL_ASSIM(SIBPT,LONSIB,LATSIB,PREF,HHTI,HLTI,   &
     &SHTI,SLTI,GR_FRAC,ASSIM,C13ASSIM,RSTFAC2,TM,POOLLT,CO2T,FRACT)
              USE MODULE_SIBCONST, ONLY :                               &
     &          NPOOLPFT,                                               &
     &          NPOOLLU
              USE MODULE_SIB, ONLY :                                    &
     &          POOLL_TYPE,                                             &
     &          CO2_TYPE,                                               &
     &          FRACT_TYPE
              INTEGER(KIND=4), INTENT(IN) :: SIBPT
              REAL(KIND=4), INTENT(IN) :: LONSIB
              REAL(KIND=4), INTENT(IN) :: LATSIB
              INTEGER(KIND=4), INTENT(IN) :: PREF
              REAL(KIND=4), INTENT(IN) :: HHTI
              REAL(KIND=4), INTENT(IN) :: HLTI
              REAL(KIND=4), INTENT(IN) :: SHTI
              REAL(KIND=4), INTENT(IN) :: SLTI
              REAL(KIND=4), INTENT(IN) :: GR_FRAC(NPOOLPFT)
              REAL(KIND=8), INTENT(IN) :: ASSIM
              REAL(KIND=8), INTENT(IN) :: C13ASSIM
              REAL(KIND=8), INTENT(IN) :: RSTFAC2
              REAL(KIND=8), INTENT(IN) :: TM
              TYPE (POOLL_TYPE), INTENT(INOUT) :: POOLLT
              TYPE (CO2_TYPE), INTENT(IN) :: CO2T
              TYPE (FRACT_TYPE), INTENT(IN) :: FRACT
            END SUBROUTINE POOL_ASSIM
          END INTERFACE 
        END MODULE POOL_ASSIM__genmod
