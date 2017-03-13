
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
   integer :: my_rank, num_ranks, ierror, prow, pcol, colextent, doubleextent
   integer :: coltype, global_block_type, local_block_type
   integer, dimension(:), allocatable :: rec_requests, snd_requests

   !get rank info
   call MPI_Comm_size(MPI_COMM_WORLD, num_ranks, ierror);
   call MPI_Comm_rank(MPI_COMM_WORLD, my_rank, ierror);

   !allocate buffers
   allocate(rec_requests(2),snd_requests(2*num_ranks))

   !create vector type for column
   call MPI_Type_extent(MPI_DOUBLE,doubleextent,ierror)
   call MPI_Type_vector(scalapack_var%rowA, 1, input_var%nprow, MPI_DOUBLE, coltype, ierror)
   !we need to know the exact size of the datatype so that we can calculate the stride in the next task
   call MPI_Type_extent(coltype,colextent,ierror)
   colextent=colextent/doubleextent
   !with the known vector type extent, we can now precisely compute the offsets: we do that in bytes (using hvector) to avoid fractional strides
   call MPI_Type_hvector(scalapack_var%colA, 1, input_var%N*input_var%npcol*doubleextent, coltype, global_block_type, ierror)
   call MPI_Type_contiguous(scalapack_var%rowA*scalapack_var%colA, MPI_DOUBLE, local_block_type, ierror)

   !commit types
   call MPI_Type_commit(global_block_type,ierror)
   call MPI_Type_commit(local_block_type,ierror)

   !iterate over the processes and do sendrecieves:
   !everybody receives from rank 0
   call MPI_Irecv(A,1,local_block_type,0,0,MPI_COMM_WORLD,rec_requests(1),ierror)
   call MPI_Irecv(Z,1,local_block_type,0,0,MPI_COMM_WORLD,rec_requests(2),ierror)
   if(my_rank==0) then
    do prow=1,input_var%nprow
      do pcol=1,input_var%npcol
        call MPI_Isend(GA(prow,pcol),1,global_block_type,(pcol-1)+input_var%npcol*(prow-1),0,MPI_COMM_WORLD,snd_requests(2*(pcol-1+input_var%npcol*(prow-1))+1),ierror)
        call MPI_Isend(GZ(prow,pcol),1,global_block_type,(pcol-1)+input_var%npcol*(prow-1),0,MPI_COMM_WORLD,snd_requests(2*(pcol-1+input_var%npcol*(prow-1))+2),ierror)
      end do
    end do
   end if

   !wait for the sends to finish
   if(my_rank==0) call MPI_Waitall(2*num_ranks,snd_requests,MPI_STATUSES_IGNORE,ierror)
   call MPI_Waitall(2,rec_requests,MPI_STATUSES_IGNORE,ierror)
   call MPI_Barrier(MPI_COMM_WORLD,ierror)

   !clean up
   deallocate(rec_requests,snd_requests)
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
