        !COMPILER-GENERATED INTERFACE MODULE: Fri Sep 26 16:45:22 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE SOIL_WATER__genmod
          INTERFACE 
            SUBROUTINE SOIL_WATER(INFIL,ECT,EGS,POROS,SATCO,BEE,PHSAT,  &
     &DZMM,ZMM,ROOTR,TD,VOL_LIQ,EFF_POROS,ROFF,DW_LIQ)
              REAL(KIND=8), INTENT(IN) :: INFIL
              REAL(KIND=8), INTENT(IN) :: ECT
              REAL(KIND=8), INTENT(IN) :: EGS
              REAL(KIND=8), INTENT(IN) :: POROS
              REAL(KIND=8), INTENT(IN) :: SATCO
              REAL(KIND=8), INTENT(IN) :: BEE
              REAL(KIND=8), INTENT(IN) :: PHSAT
              REAL(KIND=8), INTENT(IN) :: DZMM(10)
              REAL(KIND=8), INTENT(IN) :: ZMM(10)
              REAL(KIND=8), INTENT(IN) :: ROOTR(10)
              REAL(KIND=8), INTENT(IN) :: TD(10)
              REAL(KIND=8), INTENT(IN) :: VOL_LIQ(10)
              REAL(KIND=8), INTENT(IN) :: EFF_POROS(10)
              REAL(KIND=8), INTENT(INOUT) :: ROFF
              REAL(KIND=8), INTENT(OUT) :: DW_LIQ(1:10)
            END SUBROUTINE SOIL_WATER
          END INTERFACE 
        END MODULE SOIL_WATER__genmod
