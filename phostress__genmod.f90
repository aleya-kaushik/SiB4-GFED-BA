        !COMPILER-GENERATED INTERFACE MODULE: Fri Sep 26 16:45:25 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE PHOSTRESS__genmod
          INTERFACE 
            SUBROUTINE PHOSTRESS(PNUM,PHYSCONT,PS,ETC,TC,EACAS,RB,ECMASS&
     &,TD1,TD2,PAWFZW,TAWFRW,TCMIN,CO2T,VPD,ECT)
              USE MODULE_SIB, ONLY :                                    &
     &          CO2_TYPE
              USE MODULE_PARAM, ONLY :                                  &
     &          PHYS_PARAM
              INTEGER(KIND=4), INTENT(IN) :: PNUM
              TYPE (PHYS_PARAM), INTENT(IN) :: PHYSCONT
              REAL(KIND=8), INTENT(IN) :: PS
              REAL(KIND=8), INTENT(IN) :: ETC
              REAL(KIND=8), INTENT(IN) :: TC
              REAL(KIND=8), INTENT(IN) :: EACAS
              REAL(KIND=8), INTENT(IN) :: RB
              REAL(KIND=8), INTENT(IN) :: ECMASS
              REAL(KIND=8), INTENT(IN) :: TD1
              REAL(KIND=8), INTENT(IN) :: TD2
              REAL(KIND=8), INTENT(IN) :: PAWFZW
              REAL(KIND=8), INTENT(IN) :: TAWFRW
              REAL(KIND=8), INTENT(INOUT) :: TCMIN
              TYPE (CO2_TYPE), INTENT(INOUT) :: CO2T
              REAL(KIND=8), INTENT(IN) :: VPD
              REAL(KIND=8), INTENT(IN) :: ECT
            END SUBROUTINE PHOSTRESS
          END INTERFACE 
        END MODULE PHOSTRESS__genmod
