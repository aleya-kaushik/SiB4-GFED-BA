        !COMPILER-GENERATED INTERFACE MODULE: Fri Sep 26 16:45:20 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE FLUX_RBRD__genmod
          INTERFACE 
            SUBROUTINE FLUX_RBRD(Z2,U2,TC,TA,TD,SNOW_CVFC,LAI,CCC1,CCC2,&
     &RBC,RDC,RB,RD)
              REAL(KIND=4), INTENT(IN) :: Z2
              REAL(KIND=8), INTENT(IN) :: U2
              REAL(KIND=8), INTENT(IN) :: TC
              REAL(KIND=8), INTENT(IN) :: TA
              REAL(KIND=8), INTENT(IN) :: TD
              REAL(KIND=8), INTENT(IN) :: SNOW_CVFC
              REAL(KIND=8), INTENT(IN) :: LAI
              REAL(KIND=8), INTENT(IN) :: CCC1
              REAL(KIND=8), INTENT(IN) :: CCC2
              REAL(KIND=8), INTENT(INOUT) :: RBC
              REAL(KIND=8), INTENT(INOUT) :: RDC
              REAL(KIND=8), INTENT(INOUT) :: RB
              REAL(KIND=8), INTENT(INOUT) :: RD
            END SUBROUTINE FLUX_RBRD
          END INTERFACE 
        END MODULE FLUX_RBRD__genmod
