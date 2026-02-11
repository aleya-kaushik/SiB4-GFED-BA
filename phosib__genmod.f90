        !COMPILER-GENERATED INTERFACE MODULE: Fri Sep 26 16:45:25 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE PHOSIB__genmod
          INTERFACE 
            SUBROUTINE PHOSIB(PHYSCONT,GPROGT,GDIAGT,TC,TCAS,RA,RB,     &
     &RESP_LEAF,RESP_AUTO,RESP_HET,VEGT,CO2T)
              USE MODULE_SIB, ONLY :                                    &
     &          GPROG_TYPE,                                             &
     &          GDIAG_TYPE,                                             &
     &          VEG_TYPE,                                               &
     &          CO2_TYPE
              USE MODULE_PARAM, ONLY :                                  &
     &          PHYS_PARAM
              TYPE (PHYS_PARAM), INTENT(IN) :: PHYSCONT
              TYPE (GPROG_TYPE), INTENT(INOUT) :: GPROGT
              TYPE (GDIAG_TYPE), INTENT(IN) :: GDIAGT
              REAL(KIND=8), INTENT(IN) :: TC
              REAL(KIND=8), INTENT(IN) :: TCAS
              REAL(KIND=8), INTENT(IN) :: RA
              REAL(KIND=8), INTENT(IN) :: RB
              REAL(KIND=8), INTENT(IN) :: RESP_LEAF
              REAL(KIND=8), INTENT(IN) :: RESP_AUTO
              REAL(KIND=8), INTENT(IN) :: RESP_HET
              TYPE (VEG_TYPE), INTENT(INOUT) :: VEGT
              TYPE (CO2_TYPE), INTENT(INOUT) :: CO2T
            END SUBROUTINE PHOSIB
          END INTERFACE 
        END MODULE PHOSIB__genmod
