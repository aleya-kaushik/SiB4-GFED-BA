        !COMPILER-GENERATED INTERFACE MODULE: Fri Sep 26 16:45:33 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE UPDATE_LST__genmod
          INTERFACE 
            SUBROUTINE UPDATE_LST(SUBCOUNT,DTDAYFRAC,DAYFRAC_LST,       &
     &NEW_DAY_LST)
              INTEGER(KIND=4), INTENT(IN) :: SUBCOUNT
              REAL(KIND=8), INTENT(IN) :: DTDAYFRAC
              REAL(KIND=8), INTENT(INOUT) :: DAYFRAC_LST(SUBCOUNT)
              LOGICAL(KIND=4), INTENT(INOUT) :: NEW_DAY_LST(SUBCOUNT)
            END SUBROUTINE UPDATE_LST
          END INTERFACE 
        END MODULE UPDATE_LST__genmod
