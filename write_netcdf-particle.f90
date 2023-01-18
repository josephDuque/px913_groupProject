! Structure of the code heavily based on the example from:
! https://people.sc.fsu.edu/~jburkardt/f_src/netcdf/netcdf.html
! The page is available via the Wayback machine at
! https://web.archive.org/web/20190623025346/http://people.sc.fsu.edu/~jburkardt/f_src/netcdf/netcdf.html (accessedNov 2021)

MODULE write_netcdf
  USE ISO_FORTRAN_ENV
  USE netcdf

  IMPLICIT NONE
  
  !Creating a derived type
  TYPE :: run_data
    CHARACTER(LEN=50) :: code_filename, initialState
    INTEGER :: Nx, Ny, iter
  END TYPE 
 
  CONTAINS
  SUBROUTINE check(stat, operation)
    USE netcdf
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: stat
    CHARACTER(LEN=*), INTENT(IN) :: operation
    IF (stat == NF90_NOERR) RETURN
    PRINT *, "Error encountered during ", operation
    PRINT *, nf90_strerror(stat)
    STOP 1
  END SUBROUTINE check

  IMPURE SUBROUTINE write_sub(positions, velocities, accels, rho_grid, phi_grid, &
      & Ex, Ey, Nx, Ny, iter, filename, neterr, code_name, init)
    INTEGER, INTENT(IN) :: Nx, Ny, iter
    REAL(REAL64), DIMENSION(2, 0:iter) :: positions, velocities, accels
    REAL(REAL64), DIMENSION(Ny, Nx) :: rho_grid, phi_grid, Ex, Ey
    INTEGER, PARAMETER :: numDims = 4
    INTEGER, DIMENSION(numDims) :: sizes, dim_ids
    CHARACTER(LEN=20), DIMENSION(numDims) :: dims
    CHARACTER(LEN=*), INTENT(IN) :: filename, code_name, init
    INTEGER, DIMENSION(2) :: posShape, gridShape
    INTEGER, DIMENSION(7) :: var_ids
    INTEGER :: neterr, file_id, i
    TYPE(run_data) :: runDat
    
    dims(1) = "y"
    dims(2) = "x"
    dims(3) = "x-y components"
    dims(4) = "iteration"

    runDat%code_filename = code_name
    runDat%initialState = init
    runDat%Nx = Nx
    runDat%Ny = Ny
    runDat%iter = iter
    
    print *, 'subroutine starts...'

    posShape = SHAPE(positions)
    gridShape = SHAPE(rho_grid)
    
    ! Concatenating dimension lengths of the two arrays to one array:
    sizes(1) = gridShape(1)
    sizes(2) = gridShape(2)
    sizes(3) = posShape(1)
    sizes(4) = posShape(2)

    print *, 'creating dataset...'
    neterr = nf90_create(filename, NF90_CLOBBER, file_id)
    call check(neterr, 'open')
    ! Checking that creating the dataset produced no error:
    IF (neterr /= nf90_noerr) THEN
      PRINT *, TRIM(nf90_strerror(neterr))
      RETURN
    END IF

    print *, 'assigning attributes...'
    neterr = nf90_put_att(file_id, NF90_GLOBAL, 'Code filename', runDat%code_filename)
    IF (neterr /= nf90_noerr) THEN
      PRINT *, TRIM(nf90_strerror(neterr))
      RETURN
    END IF

    neterr = nf90_put_att(file_id, NF90_GLOBAL, 'Initial charge density state', runDat%initialState)
    IF (neterr /= nf90_noerr) THEN
      PRINT *, TRIM(nf90_strerror(neterr))
      RETURN
    END IF
    
    neterr = nf90_put_att(file_id, NF90_GLOBAL, 'Array column length', runDat%Ny)
    IF (neterr /= nf90_noerr) THEN
      PRINT *, TRIM(nf90_strerror(neterr))
      RETURN
    END IF

    neterr = nf90_put_att(file_id, NF90_GLOBAL, 'Array row length', runDat%Nx)
    IF (neterr /= nf90_noerr) THEN
      PRINT *, TRIM(nf90_strerror(neterr))
      RETURN
    END IF

    neterr = nf90_put_att(file_id, NF90_GLOBAL, 'Number of iterations', runDat%iter)
    IF (neterr /= nf90_noerr) THEN
      PRINT *, TRIM(nf90_strerror(neterr))
      RETURN
    END IF

    print *, 'dimension ids...'
    ! Adding dimensions to the opened netCDF dataset, and assigning the dimension IDs:
    DO i = 1, numDims
      neterr = nf90_def_dim(file_id, dims(i), sizes(i), dim_ids(i))
      IF (neterr /= nf90_noerr) THEN
        PRINT *, TRIM(nf90_strerror(neterr))
        RETURN
      END IF
    END DO

    print *, 'desnity rho array id...'
    ! Define the positions array variable type, and assigning variable ID:
    neterr = nf90_def_var(file_id, "Density rho", NF90_REAL, dim_ids(1:2), var_ids(1))
    IF (neterr /= nf90_noerr) THEN
      PRINT *, TRIM(nf90_strerror(neterr))
      RETURN
    END IF

    print *, 'potential phi array id...'
    ! Define the positions array variable type, and assigning variable ID:
    neterr = nf90_def_var(file_id, "Potential phi", NF90_REAL, dim_ids(1:2), var_ids(2))
    IF (neterr /= nf90_noerr) THEN
      PRINT *, TRIM(nf90_strerror(neterr))
      RETURN
    END IF

    print *, 'Ex array id...'
    ! Define the positions array variable type, and assigning variable ID:
    neterr = nf90_def_var(file_id, "Electric field x-component", NF90_REAL, dim_ids(1:2), var_ids(3))
    IF (neterr /= nf90_noerr) THEN
      PRINT *, TRIM(nf90_strerror(neterr))
      RETURN
    END IF

    print *, 'Ey array id...'
    ! Define the positions array variable type, and assigning variable ID:
    neterr = nf90_def_var(file_id, "Electric field y-component", NF90_REAL, dim_ids(1:2), var_ids(4))
    IF (neterr /= nf90_noerr) THEN
      PRINT *, TRIM(nf90_strerror(neterr))
      RETURN
    END IF

    print *, 'positions array id...'
    ! Define the positions array variable type, and assigning variable ID:
    neterr = nf90_def_var(file_id, "Positions", NF90_REAL, dim_ids(3:4), var_ids(5))
    IF (neterr /= nf90_noerr) THEN
      PRINT *, TRIM(nf90_strerror(neterr))
      RETURN
    END IF

    print *, 'velocities array id...'
    ! Define the positions array variable type, and assigning variable ID:
    neterr = nf90_def_var(file_id, "Velocities", NF90_REAL, dim_ids(3:4), var_ids(6))
    IF (neterr /= nf90_noerr) THEN
      PRINT *, TRIM(nf90_strerror(neterr))
      RETURN
    END IF

    print *, 'accelerations array id...'
    ! Define the positions array variable type, and assigning variable ID:
    neterr = nf90_def_var(file_id, "Accelerations", NF90_REAL, dim_ids(3:4), var_ids(7))
    IF (neterr /= nf90_noerr) THEN
      PRINT *, TRIM(nf90_strerror(neterr))
      RETURN
    END IF

    print *, 'switching to write mode...'
    ! Finish defining metadata, switch from define mode to write mode:
    neterr = nf90_enddef(file_id)
    IF (neterr /= nf90_noerr) THEN
      PRINT *, TRIM(nf90_strerror(neterr))
      RETURN
    END IF

    print *, 'writing density rho array...'
    ! Write the variables:
    neterr = nf90_put_var(file_id, var_ids(1), rho_grid)
    IF (neterr /= nf90_noerr) THEN
      PRINT *, TRIM(nf90_strerror(neterr))
      RETURN
    END IF

    print *, 'writing potential phi array...'
    ! Write the variables:
    neterr = nf90_put_var(file_id, var_ids(2), phi_grid)
    IF (neterr /= nf90_noerr) THEN
      PRINT *, TRIM(nf90_strerror(neterr))
      RETURN
    END IF

    print *, 'writing Ex array...'
    ! Write the variables:
    neterr = nf90_put_var(file_id, var_ids(3), Ex)
    IF (neterr /= nf90_noerr) THEN
      PRINT *, TRIM(nf90_strerror(neterr))
      RETURN
    END IF

    print *, 'writing Ey array...'
    ! Write the variables:
    neterr = nf90_put_var(file_id, var_ids(4), Ey)
    IF (neterr /= nf90_noerr) THEN
      PRINT *, TRIM(nf90_strerror(neterr))
      RETURN
    END IF

    print *, 'writing positions array...'
    ! Write the variables:
    neterr = nf90_put_var(file_id, var_ids(5), positions)
    IF (neterr /= nf90_noerr) THEN
      PRINT *, TRIM(nf90_strerror(neterr))
      RETURN
    END IF

    print *, 'writing velocities array...'
    ! Write the variables:
    neterr = nf90_put_var(file_id, var_ids(6), velocities)
    IF (neterr /= nf90_noerr) THEN
      PRINT *, TRIM(nf90_strerror(neterr))
      RETURN
    END IF

    print *, 'writing accelerations array...'
    ! Write the variables:
    neterr = nf90_put_var(file_id, var_ids(7), accels)
    IF (neterr /= nf90_noerr) THEN
      PRINT *, TRIM(nf90_strerror(neterr))
      RETURN
    END IF
    
    print *, 'closing file...'
    ! Close file:
    neterr = nf90_close(file_id)
    IF (neterr /= nf90_noerr) THEN
      PRINT *, TRIM(nf90_strerror(neterr))
      RETURN
    END IF

    print *, 'subroutine ended...'
  END SUBROUTINE write_sub 
  
END MODULE write_netcdf
