        !COMPILER-GENERATED INTERFACE MODULE: Fri Sep 26 16:45:17 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE DELHF__genmod
          INTERFACE 
            SUBROUTINE DELHF(TM,BPS,ROS,TA,NSL,TSFC,FSS,HC,HG,HS,TC,RA, &
     &RB,RD)
              REAL(KIND=8), INTENT(IN) :: TM
              REAL(KIND=8), INTENT(IN) :: BPS(2)
              REAL(KIND=8), INTENT(IN) :: ROS
              REAL(KIND=8), INTENT(IN) :: TA
              INTEGER(KIND=1), INTENT(IN) :: NSL
              REAL(KIND=8), INTENT(IN) :: TSFC
              REAL(KIND=8), INTENT(INOUT) :: FSS
              REAL(KIND=8), INTENT(INOUT) :: HC
              REAL(KIND=8), INTENT(INOUT) :: HG
              REAL(KIND=8), INTENT(INOUT) :: HS
              REAL(KIND=8), INTENT(INOUT) :: TC
              REAL(KIND=8), INTENT(INOUT) :: RA
              REAL(KIND=8), INTENT(INOUT) :: RB
              REAL(KIND=8), INTENT(INOUT) :: RD
            END SUBROUTINE DELHF
          END INTERFACE 
        END MODULE DELHF__genmod
