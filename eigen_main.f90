PROGRAM eigen_main

   ! This is a template for calling the pdsyevd scalapack routine.
   ! The matrix is random, but a module TOMATO could be added
   ! Based on the c program from mhadi@cray.com
   ! written in fortran by pcarrier@cray.com
   use mpi
   use parameters
   use matrix_module
   use pdsyevd_module
   use elpa_module

   implicit none

   real(kind=8), allocatable :: A(:,:), Z(:,:), W(:)
   real(kind=8), allocatable :: globA(:,:), globZ(:,:)

   type (input)     :: input_var
   type (scalapack) :: scalapack_var
   type (elpa)      :: elpa_var

   call read_input(input_var)

   !allocate space and generate random matrix
   allocate( globA(input_var%N, input_var%N) )
   allocate( globZ(input_var%N, input_var%N) )
   allocate( W(input_var%N) )
   call build_random_matrix(input_var%N, globA, globZ)

   select case (input_var%solver_type)
      case ("PDSYEVD")
        !initialize
         call initialize_scalapack(input_var, scalapack_var)
         allocate( A(scalapack_var%rowA, scalapack_var%colA) )
         allocate( Z(scalapack_var%rowA, scalapack_var%colA) )
         !distribute matrix
         call distribute_matrix_scalapack(input_var, scalapack_var, globA, globZ, A, Z)
         !solve
         call pdsyevd_solver(input_var%N, scalapack_var, A, Z, W)
         !gather
         !call gather_eigenvectors_scalapack(input_var, scalapack_var, globZ, Z)
         !finalize
         call finalize_scalapack(scalapack_var)
      case ("ELPA")
        !initialize
         call initialize_elpa(input_var, elpa_var)
      case default
         write(6,*) "<solver_type> = PDSYEVD or ELPA only (for now)"
   end select

   !deallocate
   deallocate( A, Z, W )
   deallocate( globA, globZ )
END
