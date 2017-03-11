module elpa_module
   use parameters
   use mpi
   implicit none

contains

subroutine initialize_elpa(input_var, elpa_var)
   type (input), intent(in)  :: input_var
   type (elpa),  intent(out) :: elpa_var

   write(6,*) "Transfer ELSI work here"
   write(6,*) "Matrix size is ", input_var%N
   return
end subroutine initialize_elpa

end module elpa_module
