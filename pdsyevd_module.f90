
module pdsyevd_module
   use parameters
   use mpi
   implicit none


contains

subroutine initialize_scalapack(input_var, scalapack_var)
   type (input),     intent(in)  :: input_var
   type (scalapack), intent(out) :: scalapack_var

!  .. External Functions ..
   integer, external ::   numroc

   integer :: info = 0
   character(len=9) :: procOrder = "Row-major"
   integer :: moneI = -1, zeroI = 0, oneI = 1

   call blacs_get(moneI, zeroI, scalapack_var%ictxt)
   call blacs_gridinit(scalapack_var%ictxt, procOrder, input_var%nprow, input_var%npcol);
   call blacs_gridinfo(scalapack_var%ictxt, scalapack_var%dnprow, scalapack_var%dnpcol, scalapack_var%myrow, scalapack_var%mycol);

   scalapack_var%rowA = numroc(input_var%N, input_var%NBrow, scalapack_var%myrow, zeroI, input_var%nprow)
   scalapack_var%colA = numroc(input_var%N, input_var%NBcol, scalapack_var%mycol, zeroI, input_var%npcol)

   call descinit(scalapack_var%descA, input_var%N, input_var%N, input_var%NBrow, input_var%NBcol, zeroI, zeroI, scalapack_var%ictxt, scalapack_var%rowA, info)
   call descinit(scalapack_var%descZ, input_var%N, input_var%N, input_var%NBrow, input_var%NBcol, zeroI, zeroI, scalapack_var%ictxt, scalapack_var%rowA, info)

   return
end subroutine initialize_scalapack

!scalapack assumes block cyclic layout:
subroutine distribute_matrix_scalapack(input_var, scalapack_var, GA, GZ, A, Z)
   type (input),                 intent(in)  :: input_var
   type (scalapack),             intent(in)  :: scalapack_var
   real(kind=8),                 intent(in)  :: GA(input_var%N, input_var%N)
   real(kind=8),                 intent(in)  :: GZ(input_var%N, input_var%N)
   real(kind=8),                 intent(out) :: A(scalapack_var%rowA, scalapack_var%colA)
   real(kind=8),                 intent(out) :: Z(scalapack_var%rowA, scalapack_var%colA)

   !local variables
   !this assumes that all nodes get the same share, but that is right because earlier we error out if that is not the case
   integer :: my_rank, ierror, ndim
   integer, dimension(2) :: sizes, subsizes, starts
   integer :: global_block_type, local_block_type

   !get rank info
   call MPI_Comm_rank(MPI_COMM_WORLD, my_rank, ierror);

   !scatter the matrix. Create a vector datatype for that:
   ndim=2
   sizes(1)=input_var%N
   sizes(2)=input_var%N
   subsizes(1)=scalapack_var%rowA
   subsizes(2)=scalapack_var%colA
   starts=0
   call MPI_Type_create_subarray(ndim,sizes,subsizes,starts,MPI_ORDER_FORTRAN,MPI_DOUBLE,global_block_type,ierror)
   call MPI_Type_create_subarray(ndim,subsizes,subsizes,starts,MPI_ORDER_FORTRAN,MPI_DOUBLE,local_block_type,ierror)
   call MPI_Type_commit(global_block_type,ierror)
   call MPI_Type_commit(local_block_type,ierror)

   !scatter the matrix
   call MPI_Scatter(GA,1,global_block_type,A,1,local_block_type,0,MPI_COMM_WORLD,ierror)
   call MPI_Scatter(GZ,1,global_block_type,Z,1,local_block_type,0,MPI_COMM_WORLD,ierror)

   !clean up
   call MPI_Type_free(global_block_type,ierror)
   call MPI_Type_free(local_block_type,ierror)

   return
end subroutine distribute_matrix_scalapack

!call the actual solver
subroutine pdsyevd_solver(N, scalapack_var, A, Z, W)
   integer,          intent(in)  :: N
   type (scalapack), intent(in)  :: scalapack_var
   real(kind=8),     intent(in)  :: A(scalapack_var%rowA, scalapack_var%colA)
   real(kind=8),     intent(in)  :: Z(scalapack_var%rowA, scalapack_var%colA)
   real(kind=8),     intent(out) :: W(N)

   character :: jobz = 'V'
   character :: uplo = 'L'
   integer :: ia=1, ja=1, iz=1, jz=1
   real(kind=8), allocatable ::  work(:)
   integer,      allocatable :: iwork(:)
   real(kind=8) :: temp(2)
   integer :: liwork, lwork
   integer :: info = 0
   integer :: moneI = -1
   integer :: my_rank, j, ierror



   write(6,*) "PDSYEVD -- work query start"

   call pdsyevd(jobz, uplo, N, A, ia, ja, scalapack_var%descA, W, Z, iz, jz, scalapack_var%descZ, temp, moneI, liwork, moneI, info);

   write(6,*) "PDSYEVD -- work query end"

   lwork    = temp(1)
   allocate(  work( lwork) )
   allocate( iwork(liwork) )

   write(6,*) "PDSYEVD -- call start"

   call pdsyevd(jobz, uplo, N, A, ia, ja, scalapack_var%descA, W, Z, iz, jz, scalapack_var%descZ, work, lwork, iwork, liwork, info)

   write(6,*) "PDSYEVD -- call end"

   ! call finalize_lapack
   call MPI_Comm_rank(MPI_COMM_WORLD, my_rank, ierror);

   if (my_rank == 0 ) then
      do j=1, 5
         write(6,*) "W(", j, ")=", W(j), "Z(", j, ",1)=", Z(j,1), " LWORK=", lwork, " LIWORK=", liwork
      end do
   end if

   deallocate( work, iwork )

   return
end subroutine pdsyevd_solver

subroutine finalize_scalapack(scalapack_var)
   type (scalapack), intent(in)  :: scalapack_var

   integer :: zeroI = 0

   call blacs_gridexit(scalapack_var%ictxt)
   call blacs_exit(zeroI)

   return
end subroutine finalize_scalapack

end module pdsyevd_module
