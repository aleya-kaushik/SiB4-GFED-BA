        !COMPILER-GENERATED INTERFACE MODULE: Fri Sep 26 16:45:21 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE FLUX_VMF__genmod
          INTERFACE 
            SUBROUTINE FLUX_VMF(ZZWIND,ZZTEMP,Z0,ROS,SPDM,SH,THM,SHA,THA&
     &,USTAR,CUNI,CU,CT,VENTMF)
              REAL(KIND=8), INTENT(IN) :: ZZWIND
              REAL(KIND=8), INTENT(IN) :: ZZTEMP
              REAL(KIND=8), INTENT(IN) :: Z0
              REAL(KIND=8), INTENT(IN) :: ROS
              REAL(KIND=8), INTENT(IN) :: SPDM
              REAL(KIND=8), INTENT(IN) :: SH
              REAL(KIND=8), INTENT(IN) :: THM
              REAL(KIND=8), INTENT(IN) :: SHA
              REAL(KIND=8), INTENT(IN) :: THA
              REAL(KIND=8), INTENT(INOUT) :: USTAR
              REAL(KIND=8), INTENT(INOUT) :: CUNI
              REAL(KIND=8), INTENT(INOUT) :: CU
              REAL(KIND=8), INTENT(INOUT) :: CT
              REAL(KIND=8), INTENT(INOUT) :: VENTMF
            END SUBROUTINE FLUX_VMF
          END INTERFACE 
        END MODULE FLUX_VMF__genmod
