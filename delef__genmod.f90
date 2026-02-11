        !COMPILER-GENERATED INTERFACE MODULE: Fri Sep 26 16:45:17 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE DELEF__genmod
          INTERFACE 
            SUBROUTINE DELEF(PSY,ROS,EM,EA,SOILRH,RD,RA,EC,EG,ES,FWS)
              REAL(KIND=8), INTENT(IN) :: PSY
              REAL(KIND=8), INTENT(IN) :: ROS
              REAL(KIND=8), INTENT(IN) :: EM
              REAL(KIND=8), INTENT(IN) :: EA
              REAL(KIND=8), INTENT(IN) :: SOILRH
              REAL(KIND=8), INTENT(IN) :: RD
              REAL(KIND=8), INTENT(IN) :: RA
              REAL(KIND=8), INTENT(INOUT) :: EC
              REAL(KIND=8), INTENT(INOUT) :: EG
              REAL(KIND=8), INTENT(INOUT) :: ES
              REAL(KIND=8), INTENT(INOUT) :: FWS
            END SUBROUTINE DELEF
          END INTERFACE 
        END MODULE DELEF__genmod
