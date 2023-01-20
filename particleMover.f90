MODULE subroutines
    USE iso_fortran_env
    USE command_line
    IMPLICIT NONE
    CONTAINS
    
    SUBROUTINE params(Ny, Nx, rho)
        LOGICAL :: success, exists
        INTEGER, INTENT(OUT) :: Nx,Ny
        CHARACTER(len=25) :: rho

        CALL parse_args !CALL parse_args function to 'grab and store' arguments specified in command line.
        success = get_arg("Nx", Nx,exists=exists)
        IF(success) THEN
            IF(Nx <= 0_INT32) THEN    !Number of grid cells in x-direction must be a positive integer. If Nx <= 0 is supplied, code stops and error prints.
                PRINT*,'Error: Nx must be set to a positive value'
                STOP
            END IF
        ELSE
            IF(exists) THEN !Notify user of incorrect type for Nx. Stops code.
                PRINT*,"Value for 'Nx' of incorrect type/kind passed, ensure command line input for 'Nx' is suitable for type and&
                 & kind specified by INT32."
                STOP
            ELSE        !Notify user no value for Nx was passed. Stops code.
                PRINT*,"Value for 'Nx' not parsed, check command line input." 
                STOP
            END IF
        END IF

        success = get_arg("Ny", Ny,exists=exists)
        IF(success) THEN
            IF(Ny <= 0_INT32) THEN    !Number of grid cells in y-direction must be a positive integer. If Ny <= 0 is supplied, code stops and error prints.
                PRINT*,'Error: Ny must be set to a positive value'
                STOP
            END IF
        ELSE
            IF(exists) THEN !Notify user of incorrect type for Ny. Stops code.
                PRINT*,"Value for 'Ny' of incorrect type/kind passed, ensure command line input for 'Ny' is suitable for type and&
                 & kind specified by INT32."
                STOP
            ELSE        !Notify user no value for Nx was passed. Stops code.
                PRINT*,"Value for 'Ny' not parsed, check command line input."
                STOP
            END IF
        END IF

        success = get_arg("rho", rho, exists=exists)
        IF(success) THEN
        ELSE
            IF(exists) THEN !Notify user of incorrect type for rho. Stops code.
                PRINT*,"Value for 'rho' of incorrect type passed, ensure command line input for 'init' is one of &
                &'null','single' or 'double'"
                STOP
            ELSE        !Notify user no value for rho was passed. Stops code.
                PRINT*,"Value for 'rho' not parsed, check command line input"
                STOP
            END IF
        END IF

    END SUBROUTINE params


    SUBROUTINE null(lower, rho_grid, initPos, initVel)
        INTEGER :: lower
        REAL(REAL64), DIMENSION(2), INTENT(OUT) :: initPos, initVel
        REAL(REAL64), DIMENSION(lower:,lower:), INTENT(INOUT) :: rho_grid
                
        rho_grid = 0.0_REAL64       !Generate density field for rho=null condition.

        ! Setting initial position and velocity:
        initPos = 0.0_REAL64
        initVel = 0.1_REAl64

    END SUBROUTINE null


    SUBROUTINE single(lower, rho_grid,x_axis,y_axis,Nx, Ny, initPos, initVel)
        INTEGER, INTENT(IN) :: lower, Nx, Ny
        REAL(REAL64), DIMENSION(2), INTENT(OUT) :: initPos, initVel
        REAL(REAL64), DIMENSION(lower:,lower:), INTENT(INOUT) :: rho_grid
        REAL(REAL64), DIMENSION(lower:), INTENT(IN) :: x_axis, y_axis
        REAL(REAL64) :: x_div, y_div, z_1
        INTEGER :: i,j

        
        z_1 = 0.1_REAL64    !Variable assigned to minimise risk of inputting incorrect float into equation.

        !Nested loop used to generate density field for rho=single condition.      
        !Code optimized with column-major favouring loop.
        DO i=1,Nx
            x_div = (x_axis(i)/z_1)**2       
            DO j=1,Ny
                y_div = (y_axis(j)/z_1)**2
                rho_grid(j,i) = EXP(-1.0_REAL64*(x_div + y_div))
            
            END DO
        
        END DO
        

        ! Setting initial position and velocity:
        initPos(1) = 0.1_REAL64; initPos(2) = 0.0_REAL64
        initVel = 0.0_REAl64
          

    END SUBROUTINE single


    SUBROUTINE double(lower, rho_grid,x_axis,y_axis,Nx, Ny, initPos, initVel)
        INTEGER, INTENT(IN) :: lower, Nx, Ny
        REAL(REAL64), DIMENSION(2), INTENT(OUT) :: initPos, initVel
        REAL(REAL64), DIMENSION(lower:,lower:), INTENT(INOUT) :: rho_grid
        REAL(REAL64), DIMENSION(lower:), INTENT(IN) :: x_axis, y_axis
        REAL(REAL64) :: x_div, y_div, extra_divx, extra_divy, z_1, z_2
        INTEGER :: i, j
        !Variables assigned to minimise risk of inputting incorrect float into equation.
        z_1 = 0.1_REAL64    
        z_2 = 0.2_REAL64    

        !Nested loop used to generate density field for rho=double condition.      
        !Code optimized with column-major favouring loop.
        DO i=1,Nx
            x_div = ((x_axis(i)+0.25_REAL64)/z_1)

            DO j=1,Ny
                y_div = ((y_axis(j)+0.25_REAL64)/z_1)
                extra_divy = ((y_axis(j)-0.75_REAL64)/z_2)
                extra_divx = ((x_axis(i)-0.75_REAL64)/z_2)
                rho_grid(j,i) = EXP(-(x_div**2)-(y_div**2)) + EXP(-(extra_divx**2)-(extra_divy**2))
            
            END DO

        END DO
        
        ! Setting initial position and velocity:
        initPos(1) = 0.0_REAL64; initPos(2) = 0.5_REAL64
        initVel = 0.0_REAl64

    END SUBROUTINE double

    !Subroutine generates phi_grid that satisfies convergence criteria.
    SUBROUTINE phi(lower, phi_grid, Nx, Ny, rho_grid, dx, dy)
        INTEGER, INTENT(IN) :: lower, Nx, Ny
        REAL(REAL64), INTENT(IN) :: dx, dy
        REAL(REAL64), DIMENSION(lower:,lower:), INTENT(INOUT) :: rho_grid, phi_grid
        REAL(REAL64) :: dens, div1, div2, two, denom, N
        REAL(REAL64) :: etot, conv1, conv2, drms_sum, drms
        INTEGER :: i, j

        two = 2.0_REAL64 !Variable assigned to minimise risk of inputting incorrect float into equation.

        !Initialise values of convergence criteria to allow DO WHILE loop to start,
        etot = 10.0_REAL64 
        drms = 1.0_REAL64

        
        DO i=1,Nx                                                           !Calculate one initial set of elements for phi grid.
            DO j=1,Ny                                                       !ABS() in final line of DO loop initializes grid to give
                dens=rho_grid(j,i)                                          !positive values of 'dmrs', allowing convergence criteria to work.
                div1=((phi_grid(j,i+1) + phi_grid(j,i-1))/(dx**2))
                div2=((phi_grid(j+1,i) + phi_grid(j-1,i))/(dy**2))
                denom=((two/(dx**2))+(two/(dy**2)))
                phi_grid(j,i)=(-(ABS(dens-div1-div2)/denom))            
                                                                                   
            END DO                                                      

        END DO

        !Loop until convergence criteria are met.
        DO WHILE(etot/drms > 0.00001_REAL64)
            
            DO i=1,Nx
                DO j=1,Ny                                                       !Do one cycle of gauss-seidel method to generate new values for phi grid.
                    dens=rho_grid(j,i)                                          
                    div1=((phi_grid(j,i+1) + phi_grid(j,i-1))/(dx**2))
                    div2=((phi_grid(j+1,i) + phi_grid(j-1,i))/(dy**2))
                    denom=((two/(dx**2))+(two/(dy**2)))
                    phi_grid(j,i)=-((dens-div1-div2)/denom)

                END DO

            END DO

            N = REAL(Nx, kind=REAL64)*REAL(Ny, kind=REAL64)    !Make REAL to prevent REAL-INTEGER arithmetic below..

            !(Re-)initialise values for convergence criteria.
            etot = 0.0_REAL64
            drms_sum = 0.0_REAL64       
            drms = 0.0_REAL64

            !Loop through phi_grid to generate values for convergence criteria.
            !Generation of new values for convergence criteria after each Gauss-Seidel cycle allows exit from loop as soon as criteria is met.
            DO i=1,Nx
                DO j=1,Ny
                    
                    dens=rho_grid(j,i)
                    conv1 = ((phi_grid(j,i-1) - (two*(phi_grid(j,i))) + phi_grid(j,i+1))/(dx**2))       
                    conv2 = ((phi_grid(j-1,i) - (two*(phi_grid(j,i))) + phi_grid(j+1,i))/(dy**2))
                    
                    etot = etot + ABS(conv1 + conv2 - dens) !Use loop to calculate sum for etot.

                    drms_sum = drms_sum + conv1 + conv2 !Use loop to calculate the sum part of drms.
                
                END DO
            
            END DO

            drms = SQRT((1/N)*((drms_sum)))    !After sum part of drms calculated in loop, calculate full drms outside of loop.
            
        
        END DO


    END SUBROUTINE phi


    SUBROUTINE E_fields(E_field_x, E_field_y, phi_grid, dx, dy, Nx, Ny,lower)
        INTEGER, INTENT(IN) :: Nx, Ny, lower
        REAL(REAL64), DIMENSION(:,:), INTENT(INOUT) :: E_field_x, E_field_y
        REAL(REAL64), DIMENSION(lower:,lower:), INTENT(IN) :: phi_grid
        REAL(REAL64), INTENT(IN) :: dx, dy
        REAL(REAL64) :: two
        INTEGER :: i, j

        two = 2.0_REAL64    !Variable assigned to minimise risk of inputting incorrect float into equation.
        !Use central difference to find gradients of phi in x and y directions (x and y electric fields).
        DO i=1,Nx
            Do j=1,Ny
                E_field_x(j,i) = ((phi_grid(j,i+1) - phi_grid(j,i-1))/(two*dx))     
                E_field_y(j,i) = ((phi_grid(j+1,i) - phi_grid(j-1,i))/(two*dy))

            END DO
        
        END DO
    END SUBROUTINE

!-------------------------------------------------------------
! Subroutine to convert a particles's real position into
! an integer value corresponding to the cell the particle is
! in.
    SUBROUTINE convertPos(particlePos, cellPos, dx, dy)
    INTEGER, DIMENSION(2), INTENT(OUT) :: cellPos
    REAL(REAL64), INTENT(IN) :: dx, dy
    REAL(REAL64), DIMENSION(2), INTENT(IN) :: particlePos
    
    cellPos(1) = FLOOR((particlePos(1) - 1.0_REAL64)/dx) + 1
    cellPos(2) = FLOOR((particlePos(2) - 1.0_REAL64)/dy) + 1
    
    END SUBROUTINE convertPos
!-------------------------------------------------------------
    
    
!-------------------------------------------------------------
! Verlet algorithm subroutine:
    SUBROUTINE verletAlg(particlePos, particleVel, particleAcc, Ex, Ey, q, dx, dy, dt)
    INTEGER, DIMENSION(2) :: indexPos
    REAL(REAL64), INTENT(IN) :: q, dx, dy, dt
    REAL(REAL64), DIMENSION(2) :: tempPos, tempVel, tempAcc
    REAL(REAL64), DIMENSION(2), INTENT(INOUT) :: particlePos, particleVel, particleAcc
    REAL(REAL64), DIMENSION(:, :), INTENT(IN), ALLOCATABLE :: Ex, Ey
    
    !-------------------------------
    ! Update positions into a temporary array:
    tempPos(1) = particlePos(1) + particleVel(1)*dt + 0.5_REAL64*particleAcc(1)*(dt**2)
    tempPos(2) = particlePos(2) + particleVel(2)*dt + 0.5_REAL64*particleAcc(2)*(dt**2)
    !-------------------------------
    
    !-------------------------------
    ! Convert positions into indexes:
    CALL convertPos(particlePos, indexPos, dx, dy)
    !-------------------------------
    
    !-------------------------------
    ! Update acceleration into a temporary array:
    tempAcc(1) = q * Ex(indexPos(2), indexPos(1))
    tempAcc(2) = q * Ey(indexPos(2), indexPos(1))
    !-------------------------------
    
    !-------------------------------
    ! Update velocity into a temporary array:
    tempVel(1) = particleVel(1) + 0.5_REAL64*dt*(tempAcc(1) + particleAcc(1))
    tempVel(2) = particleVel(2) + 0.5_REAL64*dt*(tempAcc(2) + particleAcc(2))
    !-------------------------------
    
    !-------------------------------
    ! Save the updated values:
    particlePos(1) = tempPos(1)
    particlePos(2) = tempPos(2)
    
    particleAcc(1) = tempAcc(1)
    particleAcc(2) = tempAcc(2)
    
    particleVel(1) = tempVel(1)
    particleVel(2) = tempVel(2)
    !-------------------------------
    
    END SUBROUTINE verletAlg
!-------------------------------------------------------------


END MODULE subroutines


PROGRAM main                !REMEMBER KING, I HAVE INDEXED IN STRANGE WAYS (J,I)
    USE iso_fortran_env
    USE subroutines
    USE domain_tools
    USE write_netcdf
    IMPLICIT NONE
    INTEGER :: Ny, Nx, i, lower
    INTEGER, PARAMETER :: iterations = 1000
    INTEGER, DIMENSION(2) :: pPosIndex
    REAL(REAL64) :: dx, dy
    REAL(REAL64), PARAMETER :: dt = (1.0_REAL64)*1.0e-2, eCharge = -1.0_REAL64 !eCharge = (1.602176634_REAL64)*1.0e-19
    REAL(REAL64), DIMENSION(2) :: pPos, pVel, pAcc, axis_range, axis_range2
    REAL(REAL64), DIMENSION(2, 0:iterations) :: positions, velocities, accels
    REAL(REAL64), DIMENSION(:,:), ALLOCATABLE :: grid, rho_grid, phi_grid, E_field_x, E_field_y
    REAL(REAL64), DIMENSION(:), ALLOCATABLE :: x_axis, y_axis
    CHARACTER(len=25) :: rho
    CHARACTER(LEN=*), PARAMETER :: netcdf_filename = 'particle-move.nc'
    CHARACTER(LEN=*), PARAMETER :: code_filename = 'particleMover.f90'
    INTEGER :: neterr
    
    CALL params(Ny, Nx, rho) !This subroutine retrieves values inputted into command line for Nx and Ny parameters.

    ALLOCATE(grid(0:Ny+1,0:Nx+1))
!-------------------------------------------------Initialise All Cells, Then Set Ghost Cells-----------------------------------------------------
    grid=1.0_REAL64
    
    DO i=0,Nx+1
        grid(0,i) = 0.0_REAL64
        grid((Ny+1),i) = 0.0_REAL64
    END DO

    DO i=0,Ny+1
        grid(i,0) = 0.0_REAL64
        grid(i,(Nx+1)) = 0.0_REAL64
    END DO

!-------------------------------------------------Set axes to be indexed later-------------------------------------------

    axis_range(1)=-1.0_REAL64
    axis_range(2)=1.0_REAL64        !Use of two axis ranges allows correct orientation of axes in plots.
    axis_range2(1)=1.0_REAL64
    axis_range2(2)=-1.0_REAL64
    CALL create_axis(x_axis, Nx, axis_range, nghosts=1, delta=dx)       !Use create_axis subroutine to generate axes with points aligning with centres of elements in grids.
    CALL create_axis(y_axis, Ny, axis_range2, nghosts=1, delta=dy)      !Used added functionality of calculating dx and dy, as well as including positions of ghost nodes.

!--------------------------------------------------Generate rho_grid-----------------------------------------------------
    
    lower=0_INT32
    ALLOCATE(rho_grid(0:Ny+1,0:Nx+1))
    rho_grid = grid     !Initialise rho_grid with values. Does not matter element values subsequently overwritten.
    IF(rho == 'null' .OR. rho == 'NULL' .OR. rho == 'Null') THEN !It is plausible that users will use upper and lower case versions of state names interchangeably. 
      CALL null(lower,rho_grid, pPos, pVel)                      
        
        
    ELSE IF(rho == 'single' .OR. rho == 'SINGLE' .OR. rho == 'Single' ) THEN
      CALL single(lower, rho_grid,x_axis,y_axis, Nx, Ny, pPos, pVel)
            
    
    ELSE IF(rho == 'double' .OR. rho == 'DOUBLE' .OR. rho=='Double') THEN
      CALL double(lower, rho_grid,x_axis,y_axis,Nx, Ny, pPos, pVel)
    
            
    ELSE
      PRINT*,"Error: Invalid state initialisation requested, choose one of 'null', 'single' or 'double'. &
        &Write in lowercase only."           !If no valid state of initialisation is requested, this program will not work as 
      STOP                                   !intended. In response to an invalid input, an error prints and the code stops running. 
                                             !An explanation of allowed states is provided.
    END IF
    
    !---------------------------------------Generate phi_grid------------------------------------------------------------
    ALLOCATE(phi_grid(0:Ny+1,0:Nx+1))
    phi_grid = grid     !Initialise phi_grid with boundary of 0 values, and +1 inside. 
                        !It was suspected this would reduce Gauss-Seidel cycles required for convergence in phi subroutine, in comparison to array of zeros.
    ALLOCATE(E_field_x(Ny,Nx))
    ALLOCATE(E_field_y(Ny,Nx))
    
    
    CALL phi(lower, phi_grid, Nx, Ny, rho_grid, dx, dy)     !Call subroutine to generate phi_grid that satisfies convergence criteria.
    CALL E_fields(E_field_x, E_field_y, phi_grid, dx, dy, Nx, Ny,lower)     !Call subroutine to generate arrays containing values for electric fields.


    !-----------------------------------------------------------------
    ! Using the electric fields to set the inital acceleration
    CALL convertPos(pPos, pPosIndex, dx, dy)

    !---------------------------------------------------
    ! Initial acceleration:
    pAcc(1) = eCharge * E_field_x(pPosIndex(2), pPosIndex(1))
    pAcc(2) = eCharge * E_field_y(pPosIndex(2), pPosIndex(1))
    !---------------------------------------------------
    !-----------------------------------------------------------------
    
    
    !-----------------------------------------------------------------
    ! Saving the inital position, velocity and acceleration into
    ! the arrays that are outputted to netcdf.
    positions(:, 0) = pPos; velocities(:, 0) = pVel; accels(:, 0) = pAcc
    !-----------------------------------------------------------------


    !----------------------------------------------------------------- 
    ! Perform the verlet algorithm for 1000 iterations, and at each step
    ! save the positions, velocities and accelerations into the netcdf
    ! arrays:
    DO i = 1, iterations
      CALL verletAlg(pPos, pVel, pAcc, E_field_x, E_field_y, eCharge, dx, dy, dt)
      positions(:, i) = pPos; velocities(:, i) = pVel; accels(:, i) = pAcc
    END DO
    !-----------------------------------------------------------------


    !----------------------------------------------------------------- 
    ! Call the netcdf writer subroutine:
    CALL write_sub(positions, velocities, accels, rho_grid(1:Ny, 1:Nx), phi_grid(1:Ny, 1:Nx), E_field_x, E_field_y, Nx, Ny, &
     & iterations, netcdf_filename, neterr, code_filename, rho)
    !----------------------------------------------------------------- 
    

    DEALLOCATE(rho_grid)
    DEALLOCATE(grid)
    DEALLOCATE(phi_grid)
    DEALLOCATE(E_field_x)
    DEALLOCATE(E_field_y)




END PROGRAM main




