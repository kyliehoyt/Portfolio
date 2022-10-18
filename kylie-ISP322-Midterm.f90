module iteration
type trial
        complex :: startx
        integer :: iters
end type trial
end module

program NewtonMethod
use iteration
implicit none
real, parameter :: epsilon = 0.0001
complex :: xstart
real :: s, t, realpart, imagpart
integer, parameter :: Max_Itns = 100
integer :: i
type(trial), dimension(160000) :: trials

i = 1
do s=-2,2,0.01
        do t=-2,2,0.01
                xstart = cmplx(s, t)
                trials(i)%startx = xstart
                trials(i)%iters = converge(xstart, Max_Itns) ! Function Call
                realpart = real(xstart)
                imagpart = aimag(xstart)
                print *, realpart, ", ", imagpart, ", ", trials(i)%iters ! csv
                i = i+1
        end do
end do


contains
        complex function f(x)
                implicit none
                complex :: x
                f = x**3 - 1
        end function
        complex function p(x)
                implicit none
                complex :: x
                p = 3*x**2
        end function
        complex function x_new(x)
                implicit none
                complex :: x
                x_new = x - f(x)/p(x)
        end function
        integer function converge(x0, MaxIter)
                implicit none
                complex :: x0, x_n, x
                integer :: iter, MaxIter
                x_n = x_new(x0)
                x = x0
                iter = 1
                do while ((abs(x_n - x) > epsilon) .and. iter<MaxIter)
                        x = x_n
                        if(p(x_n) /= 0) then ! if denominator is not zero
                                x_n = x_new(x_n) ! get new x_n
                        else
                                iter = 100
                                exit
                        end if 
                        iter = iter+1
                end do 
                converge = iter
        end function
end program
