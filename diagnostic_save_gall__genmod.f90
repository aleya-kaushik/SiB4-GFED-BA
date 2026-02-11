        !COMPILER-GENERATED INTERFACE MODULE: Tue Jan 13 14:42:01 2026
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE DIAGNOSTIC_SAVE_GALL__genmod
          INTERFACE 
            SUBROUTINE DIAGNOSTIC_SAVE_GALL(GNUM,SIBG,GPROGT)
              USE MODULE_SIB, ONLY :                                    &
     &          GRIDCELL_TYPE,                                          &
     &          GPROG_TYPE
              INTEGER(KIND=4), INTENT(IN) :: GNUM
              TYPE (GRIDCELL_TYPE), INTENT(IN) :: SIBG
              TYPE (GPROG_TYPE), INTENT(IN) :: GPROGT
            END SUBROUTINE DIAGNOSTIC_SAVE_GALL
          END INTERFACE 
        END MODULE DIAGNOSTIC_SAVE_GALL__genmod
