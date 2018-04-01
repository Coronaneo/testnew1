	subroutine nufft1dII(nj,x,iflag,ns,rt,tol,U,V,xsub,r)
	implicit none

	integer :: ns,rt,iflag,nj,k,j,r,xsub(nj)
	real  :: tol
	real*8 pi,x(nj)
	parameter (pi=3.141592653589793238462643383279502884197d0)
	complex*16 fftconst,U(nj,r),V(ns,r)
	complex*16 M(nj,ns)

	fftconst = iflag*dcmplx(0,1)/ns*2*pi

	do k = 1,nj
	   do j = 1,ns
	      M(k,j) = exp(fftconst*(j-1)*(x(k)-floor(x(k)+0.5)))
	   enddo
	enddo

	!call lowrankfac(M,tol,rt,rt,U,V)
	xsub = mod(floor(x+0.5),ns)+1
	r = size(V,2)

	end subroutine

	subroutine nufft1dIIapp(nj,plan,c,U,V,xsub,ns,iflag,r,S)
	integer  r,i,j,k,nj,ns,iflag,num
	integer mm,xsub(nj)
	complex*16 M(r,ns),N(r,ns),S(nj),c(ns),U(r,nj),V(r,ns)
	complex*16 NN(nj,r),CC(nj,r),w1,w2,Idx(nj,ns)
        complex*16,allocatable :: MMM(:,:)
	double complex in1, out1
	dimension in1(ns), out1(nj)
	integer*8 :: plan
        !character*8 date
        !character*10 time
        !character*5 zone 
        !integer*4 values1(8),values2(8)
        !real*16 arr(4),time1
        !num=1000
        !arr(1)=3600
        !arr(2)=60
        !arr(3)=1
        !arr(4)=0.001

	do k = 1,ns
	   M(:,k) = V(:,k)*c(k)
	enddo

        
        mm=floor(ns/2.0+0.6)
        allocate(MMM(r,mm))
        MMM=M(:,1:mm)
        M(:,1:ns-mm)=M(:,mm+1:ns)
        M(:,ns-mm+1:ns)=MMM

	do i = 1,r
	   in1 = M(i,:)
	   call dfftw_execute_dft(plan, in1, out1)
	   N(i,:) = out1
	enddo

        !do i = 1,nj
        !   CC(i,:)=N(xsub(i),:)
        !enddo
	S = sum(U*N(:,xsub),1)

	end subroutine
