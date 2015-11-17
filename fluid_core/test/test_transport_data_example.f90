program test_transport_data_example
   use m_transport_data
   use m_lookup_table

   implicit none

   integer, parameter          :: dp         = kind(0.0d0)
   integer, parameter          :: num_rows_lookup = 1000 ! Size of lookup table
   integer, parameter          :: test_size  = 10   ! Size of test array
   character(len=*), parameter :: td_file    = "test_td_input.txt"
   character(len=*), parameter :: gas_name   = "N2"

   integer                     :: ix, col_ix
   real(dp)                    :: rand_vals(test_size)
   real(dp)                    :: test_results(test_size)
   real(dp), parameter         :: x_min      = 0, x_max = 1e7_dp
   real(dp), allocatable       :: x_data(:), y_data(:)

   type(LT_col_t)              :: lkp_tbl
   type(LT_mcol_t)             :: lkp_tbl_m

   ! First example: use only one column
   ! ----------------------------------

   ! Create a lookup table
   lkp_tbl = LT_create_col(x_min, x_max, num_rows_lookup)

   ! Get the data from the input file
   call TD_get_td_from_file(td_file, gas_name, &
        "efield[V/m]_vs_energy[eV]", x_data, y_data)

   ! Store it in the lookup table
   call LT_set_col(lkp_tbl, x_data, y_data)

   print *, "Here are some random values from the lookup table"
   ! Get random locations in range x_min : x_max
   call random_number(rand_vals)
   rand_vals = x_min + rand_vals * (x_max - x_min)

   do ix = 1, test_size
      write(*,"(E8.2,A,E8.2)") rand_vals(ix), " | ", LT_get_col(lkp_tbl, rand_vals(ix))
   end do

   print *, "You can also get multiple values at the same time"
   test_results = LT_get_col(lkp_tbl, rand_vals)

   do ix = 1, test_size
      write(*,"(E8.2,A,E8.2)") rand_vals(ix), " | ", test_results(ix)
   end do

   ! Second example: multiple columns
   ! ----------------------------------

   ! Create lookup table with 3 columns
   lkp_tbl_m = LT_create_mcol(x_min, x_max, num_rows_lookup, 3)

   ! Set first column
   call TD_get_td_from_file(td_file, gas_name, &
        "efield[V/m]_vs_energy[eV]", x_data, y_data)
   call LT_set_mcol(lkp_tbl_m, 1, x_data, y_data)

   ! Set second column
   call TD_get_td_from_file(td_file, gas_name, &
        "efield[V/m]_vs_mu[m2/Vs]", x_data, y_data)
   call LT_set_mcol(lkp_tbl_m, 2, x_data, y_data)

   ! Set third column
   call TD_get_td_from_file(td_file, gas_name, &
        "efield[V/m]_vs_dif[m2/s]", x_data, y_data)
   call LT_set_mcol(lkp_tbl_m, 3, x_data, y_data)

   print *, "Here are some random values from the lookup table"
   call random_number(rand_vals)
   rand_vals = x_min + rand_vals * (x_max - x_min)

   do ix = 1, test_size
      write(*,"(E8.2,A,3E9.2)") rand_vals(ix), " | ", LT_get_mcol(lkp_tbl_m, rand_vals(ix))
   end do

   print *, "You can also get them at the same time like this"
   do col_ix = 1, 3
      print *, "Column", col_ix
      test_results = LT_get_1col(lkp_tbl_m, col_ix, rand_vals)
      do ix = 1, test_size
         write(*,"(E8.2,A,E8.2)") rand_vals(ix), " | ", test_results(ix)
      end do
   end do

   print *, "Note that LT_get_col and LT_get_1col are elemental!"

end program test_transport_data_example
