        !COMPILER-GENERATED INTERFACE MODULE: Fri Sep 26 16:45:24 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE PHOCYCALC__genmod
          INTERFACE 
            SUBROUTINE PHOCYCALC(TOA_PAR,PAR,APARKK,VMAX,VMAXSS,ATHETA, &
     &BTHETA,GAMMA,RRKK,OMSS,C3,C4,PCO2I,OMC,OME,OMS,ASSIM)
              REAL(KIND=8), INTENT(IN) :: TOA_PAR
              REAL(KIND=8), INTENT(IN) :: PAR
              REAL(KIND=8), INTENT(IN) :: APARKK
              REAL(KIND=8), INTENT(IN) :: VMAX
              REAL(KIND=8), INTENT(IN) :: VMAXSS
              REAL(KIND=8), INTENT(IN) :: ATHETA
              REAL(KIND=8), INTENT(IN) :: BTHETA
              REAL(KIND=8), INTENT(IN) :: GAMMA
              REAL(KIND=8), INTENT(IN) :: RRKK
              REAL(KIND=8), INTENT(IN) :: OMSS
              REAL(KIND=8), INTENT(IN) :: C3
              REAL(KIND=8), INTENT(IN) :: C4
              REAL(KIND=8), INTENT(IN) :: PCO2I
              REAL(KIND=8), INTENT(INOUT) :: OMC
              REAL(KIND=8), INTENT(INOUT) :: OME
              REAL(KIND=8), INTENT(INOUT) :: OMS
              REAL(KIND=8), INTENT(INOUT) :: ASSIM
            END SUBROUTINE PHOCYCALC
          END INTERFACE 
        END MODULE PHOCYCALC__genmod
