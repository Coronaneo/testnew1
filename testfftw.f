	program test
	implicit none
	complex*16 in1(8,2),out1(8,2)
	integer*16 plan,i
	real*16 pi
	parameter (pi=3.141592653589793238462643383279502884197d0)

	call dfftw_plan_dft_2d(plan,8,in1,out1,-1,0)

	do i = 1,8
	   in1(i,1)=pi*i/8
           in1(i,2)=-pi*i/8
	enddo

	call dfftw_execute_dft(plan, in1, out1)

	print *,out1

	end program
