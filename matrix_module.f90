module matrix_module
   use parameters
   use mpi
   implicit none


contains

subroutine build_random_matrix(N, A, Z)
   integer,          intent(in)  :: N
   real(kind=8),     intent(out) :: A(N,N)
   real(kind=8),     intent(out) :: Z(N,N)
   integer :: oneI = 1

   integer :: ISEED(4)
   integer :: numEle

   ISEED(1) = 1
   ISEED(2) = 2
   ISEED(3) = 3
   ISEED(4) = 4

   numEle = N*N;
   call dlarnv(oneI, ISEED, numEle, A);
   call dlarnv(oneI, ISEED, numEle, Z);
   return
end subroutine build_random_matrix

end module matrix_module
