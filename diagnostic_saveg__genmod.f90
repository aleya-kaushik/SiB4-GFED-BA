        !COMPILER-GENERATED INTERFACE MODULE: Tue Jan 13 14:42:01 2026
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE DIAGNOSTIC_SAVEG__genmod
          INTERFACE 
            SUBROUTINE DIAGNOSTIC_SAVEG(NVARS,NSAVE,REFVARS,REFSAVE,SIBG&
     &,OUTSAVE,GPROGT)
              USE MODULE_SIB, ONLY :                                    &
     &          GRIDCELL_TYPE,                                          &
     &          GPROG_TYPE
              INTEGER(KIND=4), INTENT(IN) :: NSAVE
              INTEGER(KIND=4), INTENT(IN) :: NVARS
              INTEGER(KIND=4), INTENT(IN) :: REFVARS(NVARS)
              INTEGER(KIND=4), INTENT(IN) :: REFSAVE
              TYPE (GRIDCELL_TYPE), INTENT(IN) :: SIBG
              REAL(KIND=8), INTENT(INOUT) :: OUTSAVE(NVARS,NSAVE)
              TYPE (GPROG_TYPE), INTENT(IN) :: GPROGT
            END SUBROUTINE DIAGNOSTIC_SAVEG
          END INTERFACE 
        END MODULE DIAGNOSTIC_SAVEG__genmod
