        !COMPILER-GENERATED INTERFACE MODULE: Fri Sep 26 16:45:33 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE NC_ENSURE_VAR__genmod
          INTERFACE 
            SUBROUTINE NC_ENSURE_VAR(NCID,VARNAME,VARID,FILE,LINE)
              INTEGER(KIND=4), INTENT(IN) :: NCID
              CHARACTER(*), INTENT(IN) :: VARNAME
              INTEGER(KIND=4), INTENT(OUT) :: VARID
              CHARACTER(*), INTENT(IN) :: FILE
              INTEGER(KIND=4), INTENT(IN) :: LINE
            END SUBROUTINE NC_ENSURE_VAR
          END INTERFACE 
        END MODULE NC_ENSURE_VAR__genmod
