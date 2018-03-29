	program main

	    implicit none
	    complex*16 :: a(2,2),b(2,1),c(2,1)
	    a=reshape((/dcmplx(1,1),dcmplx(3,0),dcmplx(1,2),dcmplx(3,0)/),
     &      (/2,2/))
	    b=reshape((/dcmplx(3,0),dcmplx(5,1)/),(/2,1/))

	    call zgemm('N','N',2,1,2,dcmplx(1,0),a,2,b,2,dcmplx(0,0),c,2)
	    print *,c

	end program
