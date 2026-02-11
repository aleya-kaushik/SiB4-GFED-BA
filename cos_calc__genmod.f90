        !COMPILER-GENERATED INTERFACE MODULE: Fri Sep 26 16:45:17 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE COS_CALC__genmod
          INTERFACE 
            SUBROUTINE COS_CALC(GREF,PFTREF,LONSIB,LATSIB,I3,PCOSM,TCAS,&
     &APARKK,CO2_ASSIM,RSTFAC2,POROS,ZM,WOPTZM,WSAT,ROOTF3,DZ3,WWW_ICE3,&
     &WWW_LIQ3,SOILRESP_LAY,HRT_SOIL_TEMP,TCAN,CO2T,COST,PREF,SSCOLT)
              USE MODULE_SIB, ONLY :                                    &
     &          COS_TYPE,                                               &
     &          CO2_TYPE,                                               &
     &          SSCOL_TYPE
              INTEGER(KIND=4), INTENT(IN) :: I3
              INTEGER(KIND=4), INTENT(IN) :: GREF
              INTEGER(KIND=4), INTENT(IN) :: PFTREF
              REAL(KIND=4), INTENT(IN) :: LONSIB
              REAL(KIND=4), INTENT(IN) :: LATSIB
              REAL(KIND=8), INTENT(IN) :: PCOSM
              REAL(KIND=8), INTENT(IN) :: TCAS
              REAL(KIND=8), INTENT(IN) :: APARKK
              REAL(KIND=8), INTENT(IN) :: CO2_ASSIM
              REAL(KIND=8), INTENT(IN) :: RSTFAC2
              REAL(KIND=8), INTENT(IN) :: POROS
              REAL(KIND=8), INTENT(IN) :: ZM
              REAL(KIND=8), INTENT(IN) :: WOPTZM
              REAL(KIND=8), INTENT(IN) :: WSAT
              REAL(KIND=8), INTENT(IN) :: ROOTF3(I3)
              REAL(KIND=8), INTENT(IN) :: DZ3(I3)
              REAL(KIND=8), INTENT(IN) :: WWW_ICE3(I3)
              REAL(KIND=8), INTENT(IN) :: WWW_LIQ3(I3)
              REAL(KIND=8) :: SOILRESP_LAY(10)
              REAL(KIND=8), INTENT(IN) :: HRT_SOIL_TEMP
              REAL(KIND=8), INTENT(IN) :: TCAN
              TYPE (CO2_TYPE), INTENT(IN) :: CO2T
              TYPE (COS_TYPE), INTENT(INOUT) :: COST
              INTEGER(KIND=4), INTENT(IN) :: PREF
              TYPE (SSCOL_TYPE), INTENT(IN) :: SSCOLT
            END SUBROUTINE COS_CALC
          END INTERFACE 
        END MODULE COS_CALC__genmod
