	subroutine nufft1dI(nj,x,iflag,kflag,ns,rt,tol,U,V,xsub,r)
	implicit none

	integer :: ns,rt,iflag,nj,k,j,r,xsub(nj),kflag
	real  :: tol
	real*8 pi,x(nj)
	parameter (pi=3.141592653589793238462643383279502884197d0)
	complex*16 fftconst,U(ns,r),V(nj,r)
	complex*16 M(ns,nj)
        double complex in1, out1
	dimension in1(nj), out1(nj)
	integer*8 :: plan

        call dfftw_plan_dft_1d(plan,nj,in1,out1,iflag,0)
	fftconst = iflag*dcmplx(0.0,1.0)/ns*2*pi

        if (kflag .ge. 0) then
	   do k = 1,ns
	      do j = 1,nj
	         M(k,j) = exp(fftconst*(k-1)*(x(j)-floor(x(j)+0.5)))
	      enddo
	   enddo
        else 
           do k = -ns/2,ns/2-1
              do j = 1,nj
	         M(k,j) = exp(fftconst*k*(x(j)-floor(x(j)+0.5)))
	      enddo
	   enddo
        endif
        
	!call lowrankfac(M,tol,rt,rt,U,V)

	r = size(V,2)

	xsub = mod(floor(x+0.5),ns)+1

	return

	end subroutine



        subroutine nufft1dIapp(nj,plan,c,U,V,xsub,ns,kflag,r,S)
	implicit none
	integer  r,i,j,k,nj,ns,iflag,num,kflag
        integer mm
	integer xsub(nj)
	complex*16 M(nj,r),N(ns,r),S(ns),c(nj),U(ns,r),V(nj,r)
        complex*16,allocatable :: NN(:,:)
	double complex in1, out1
        real*16  time_begin,time_end,countrage,countmax
	dimension in1(nj), out1(ns)
	integer*8  plan

	M=0
	do i = 1,nj
	   do k = 1,r
	      M(xsub(i),k) = M(xsub(i),k)+conjg(V(i,k))*c(i)
	   enddo
	enddo
        !print *,'M(1,:)=',M(1,:)
	do i = 1,r
	   in1 = M(:,i)
           call dfftw_execute_dft(plan, in1, out1)
	   N(:,i) = out1
	enddo

       
        if (kflag .lt. 0) then
           mm=floor(ns/2.0+0.6)
           allocate(NN(mm,r))
           NN=N(1:mm,:)
           N(1:ns-mm,:)=N(mm+1:ns,:)
           N(ns-mm+1:ns,:)=NN
        endif

        N=U*N
	S = sum(N,2)

        
	end
