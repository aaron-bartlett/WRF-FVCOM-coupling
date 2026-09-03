
program testTimeUtil

  use TimeUtil_mod

  implicit none

  integer :: year,dayOfYear
  integer :: month,dayOfMonth

  year = 2003
  do dayOfYear = 1,366
     call calcMonthDayOfMonth(year,dayOfYear,month,dayOfMonth)
  end do

end program testTimeUtil
