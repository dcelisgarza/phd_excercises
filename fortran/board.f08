program main
  implicit none
  integer, parameter :: dp = kind(1.0D0)
  integer, parameter :: i4b =  selected_int_kind(9)
  integer(i4b) :: i, max_it
  real(dp), allocatable :: denum(:)
  real(dp) :: tmp, x
  real(dp) :: result

  INTEGER :: icount, fcount


  write(*,*) " enter maximum number of iterations ", new_line(' ')
  read(*,*) max_it

  allocate(denum(max_it))

  CALL SYSTEM_CLOCK(icount)
  do i = 1, max_it
    denum(i) = (i - 0.5)**2/max_it**2
  end do
  denum = 1 + denum
  result = 4._dp*sum(1._dp/denum)/max_it
  CALL SYSTEM_CLOCK(fcount)
  write(*,*) fcount-icount, result

  CALL SYSTEM_CLOCK(icount)
  do i = max_it, 1, -1
    denum(i) = (i - 0.5)**2/max_it**2
  end do
  denum = 1 + denum
  result = 4._dp*sum(1._dp/denum)/max_it
  CALL SYSTEM_CLOCK(fcount)
  write(*,*) fcount-icount, result




end program main
