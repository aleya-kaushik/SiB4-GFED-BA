        !COMPILER-GENERATED INTERFACE MODULE: Fri Sep 26 16:45:26 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE POOL_ALLOC__genmod
          INTERFACE 
            SUBROUTINE POOL_ALLOC(SIBPT,LONSIB,LATSIB,PREF,HHTI,HLTI,   &
     &SHTI,SLTI,RSTFAC2,TM,ADJ_MOIST,ADJ_TEMP,ALLOC_PHEN,ALLOC_MOIST,   &
     &ALLOC_TEMP,ALLOC)
              USE MODULE_SIBCONST, ONLY :                               &
     &          NPOOLPFT
              INTEGER(KIND=4), INTENT(IN) :: SIBPT
              REAL(KIND=4), INTENT(IN) :: LONSIB
              REAL(KIND=4), INTENT(IN) :: LATSIB
              INTEGER(KIND=4), INTENT(IN) :: PREF
              REAL(KIND=4), INTENT(IN) :: HHTI
              REAL(KIND=4), INTENT(IN) :: HLTI
              REAL(KIND=4), INTENT(IN) :: SHTI
              REAL(KIND=4), INTENT(IN) :: SLTI
              REAL(KIND=8), INTENT(IN) :: RSTFAC2
              REAL(KIND=8), INTENT(IN) :: TM
              LOGICAL(KIND=4), INTENT(IN) :: ADJ_MOIST
              LOGICAL(KIND=4), INTENT(IN) :: ADJ_TEMP
              REAL(KIND=8), INTENT(IN) :: ALLOC_PHEN(NPOOLPFT)
              REAL(KIND=8), INTENT(INOUT) :: ALLOC_MOIST(NPOOLPFT)
              REAL(KIND=8), INTENT(INOUT) :: ALLOC_TEMP(NPOOLPFT)
              REAL(KIND=8), INTENT(INOUT) :: ALLOC(NPOOLPFT)
            END SUBROUTINE POOL_ALLOC
          END INTERFACE 
        END MODULE POOL_ALLOC__genmod
