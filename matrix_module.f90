module matrix_module
   use parameters
   use mpi
   implicit none


contains

subroutine build_random_matrix(N, A, Z)
   integer,          intent(in)  :: N
   real(kind=8),     intent(out) :: A(N,N)
   real(kind=8),     intent(out) :: Z(N,N)
   integer :: oneI = 1, i, j

   integer :: ISEED(4)
   integer :: numEle

   ISEED(1) = 1
   ISEED(2) = 2
   ISEED(3) = 3
   ISEED(4) = 4

   numEle = N*N;
   call dlarnv(oneI, ISEED, numEle, A);
   call dlarnv(oneI, ISEED, numEle, Z);

   !symmetrize
   do i=1,N
     do j=i+1,N
       A(i,j)=A(j,i)
     end do
   end do

   return
end subroutine build_random_matrix

end module matrix_module
