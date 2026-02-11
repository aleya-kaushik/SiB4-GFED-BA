        !COMPILER-GENERATED INTERFACE MODULE: Fri Sep 26 16:45:33 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE DAY_LENGTH__genmod
          INTERFACE 
            SUBROUTINE DAY_LENGTH(SUBCOUNT,PI180,TAN_DEC,LAT,DLENGTH,   &
     &DLENGTHDT)
              INTEGER(KIND=4), INTENT(IN) :: SUBCOUNT
              REAL(KIND=8), INTENT(IN) :: PI180
              REAL(KIND=8), INTENT(IN) :: TAN_DEC
              REAL(KIND=4), INTENT(IN) :: LAT(SUBCOUNT)
              REAL(KIND=8), INTENT(INOUT) :: DLENGTH(SUBCOUNT)
              REAL(KIND=8), INTENT(INOUT) :: DLENGTHDT(SUBCOUNT)
            END SUBROUTINE DAY_LENGTH
          END INTERFACE 
        END MODULE DAY_LENGTH__genmod
