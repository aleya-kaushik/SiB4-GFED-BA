        !COMPILER-GENERATED INTERFACE MODULE: Fri Sep 26 16:45:18 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE EQUIPOOLS_RESTART__genmod
          INTERFACE 
            SUBROUTINE EQUIPOOLS_RESTART(PNUM,POOLPFT_END,POOLPFT_EQUIB,&
     &POOLPFT_FLAY,POOLPFT_OUT)
              USE MODULE_SIBCONST, ONLY :                               &
     &          NPOOLPFT,                                               &
     &          NSOIL,                                                  &
     &          SPINUP
              INTEGER(KIND=4), INTENT(IN) :: PNUM
              REAL(KIND=8), INTENT(IN) :: POOLPFT_END(NPOOLPFT)
              REAL(KIND=8), INTENT(IN) :: POOLPFT_EQUIB(NPOOLPFT)
              REAL(KIND=8), INTENT(IN) :: POOLPFT_FLAY(NPOOLPFT,10)
              REAL(KIND=8), INTENT(INOUT) :: POOLPFT_OUT(NPOOLPFT,10)
            END SUBROUTINE EQUIPOOLS_RESTART
          END INTERFACE 
        END MODULE EQUIPOOLS_RESTART__genmod
