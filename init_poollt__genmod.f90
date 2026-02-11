        !COMPILER-GENERATED INTERFACE MODULE: Fri Sep 26 16:46:52 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE INIT_POOLLT__genmod
          INTERFACE 
            SUBROUTINE INIT_POOLLT(NSOIL,NPOOLCAN,NPOOLCANC13,NPOOLPFT, &
     &POOLLT)
              USE MODULE_SIB, ONLY :                                    &
     &          POOLL_TYPE
              INTEGER(KIND=4), INTENT(IN) :: NSOIL
              INTEGER(KIND=4), INTENT(IN) :: NPOOLCAN
              INTEGER(KIND=4), INTENT(IN) :: NPOOLCANC13
              INTEGER(KIND=4), INTENT(IN) :: NPOOLPFT
              TYPE (POOLL_TYPE), INTENT(INOUT) :: POOLLT
            END SUBROUTINE INIT_POOLLT
          END INTERFACE 
        END MODULE INIT_POOLLT__genmod
