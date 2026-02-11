        !COMPILER-GENERATED INTERFACE MODULE: Fri Sep 26 16:45:27 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE POOL_AUTO_TRAN__genmod
          INTERFACE 
            SUBROUTINE POOL_AUTO_TRAN(POOLCONT,DAYLEN,DAYLENDT,DAYLENMAX&
     &,TC,PAWFRW,ROOTF_LAY,POOLLT,GAIN_TRANSL_LAY,POOLLU_DGAIN,FRACT)
              USE MODULE_SIBCONST, ONLY :                               &
     &          NPOOLPFT,                                               &
     &          NPOOLLU,                                                &
     &          NTPOOL,                                                 &
     &          NSOIL
              USE MODULE_SIB, ONLY :                                    &
     &          POOLL_TYPE,                                             &
     &          FRACT_TYPE
              USE MODULE_PARAM, ONLY :                                  &
     &          POOL_PARAM
              TYPE (POOL_PARAM), INTENT(IN) :: POOLCONT
              REAL(KIND=8), INTENT(IN) :: DAYLEN
              REAL(KIND=8), INTENT(IN) :: DAYLENDT
              REAL(KIND=4), INTENT(IN) :: DAYLENMAX
              REAL(KIND=8), INTENT(IN) :: TC
              REAL(KIND=8), INTENT(IN) :: PAWFRW
              REAL(KIND=8), INTENT(IN) :: ROOTF_LAY(10)
              TYPE (POOLL_TYPE), INTENT(INOUT) :: POOLLT
              REAL(KIND=8), INTENT(INOUT) :: GAIN_TRANSL_LAY(NPOOLLU,10)
              REAL(KIND=8), INTENT(INOUT) :: POOLLU_DGAIN(NPOOLLU,10)
              TYPE (FRACT_TYPE), INTENT(IN) :: FRACT
            END SUBROUTINE POOL_AUTO_TRAN
          END INTERFACE 
        END MODULE POOL_AUTO_TRAN__genmod
