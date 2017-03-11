module parameters

   implicit none

   type input
      integer :: nprow, npcol
      integer :: N
      integer :: NBrow, NBcol
      character(len=8) :: solver_type
   end type input

   type scalapack
      integer :: ictxt
      integer :: dnprow, dnpcol
      integer :: myrow, mycol
      integer :: descA(15), descZ(15)
      integer :: rowA, colA
   end type scalapack

   type elpa
      integer :: something
   end type elpa

contains

subroutine read_input(input_var)
   type (input), intent(out) :: input_var

   integer :: i, b, c
   character(len=32) :: arg


   call getarg(1, arg)
   read(arg,*) b
   input_var%nprow = b
   call getarg(2, arg)
   read(arg,*) b
   input_var%npcol = b
   call getarg(3, arg)
   read(arg,*) b
   input_var%N = b
   call getarg(4, arg)
   input_var%solver_type = arg

   ! Define the block size, based on the number of MPI ranks
   c = input_var%N/input_var%nprow
   if (c*input_var%nprow /= input_var%N) then
      write(6,*) "Stopping: N (",input_var%N,") is not divisible by nprow (",input_var%nprow,")"
      STOP
   end if
   input_var%NBrow =c
   c = input_var%N/input_var%npcol
   if (c*input_var%npcol /= input_var%N) then
      write(6,*) "Stopping: N (",input_var%N,") is not divisible by npcol (",input_var%npcol,")"
      STOP
   end if
   input_var%NBcol =c

   return
end subroutine read_input

end module parameters
