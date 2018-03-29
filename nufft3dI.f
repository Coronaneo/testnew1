subroutine nufft3dI(nj,x,iflag,ns,rt,tol,U,V,xxsub,r)
implicit none

integer :: ns,rt,iflag,nj,k(ns*ns*ns,3),j,r,xsub(nj,3),i,xxsub(nj),w
real  :: tol
real*8 pi,x(nj,3)
parameter (pi=3.141592653589793238462643383279502884197d0)
complex*16 fftconst,U(ns*ns*ns,r),V(nj,r)
complex*16 M(ns*ns*ns,nj)

do j = 1,ns
   do i = 1,ns
      do w = 1,ns
	 k((j-1)*ns*ns+(i-1)*ns+w,3) = j
	 k((j-1)*ns*ns+(i-1)*ns+w,2) = i
         k((j-1)*ns*ns+(i-1)*ns+w,1) = w
      enddo
   enddo
enddo

fftconst = iflag*dcmplx(0,1)/ns*2*pi

do i = 1,ns
   do j = 1,nj
      M(i,j) = exp(fftconst*(k(i,1)*(x(j,1)-floor(x(j,1)+0.5))+
&     k(i,2)*(x(j,2)-floor(x(j,2)+0.5))+k(i,3)*(x(j,3)-floor(x(j,3)+0.5))))
   enddo
enddo

call lowrankfac(M,tol,rt,rt,U,V)

xsub = mod(floor(x+0.5),ns)+1
do i = 1,nj
   xxsub(i) = xsub(i,3)*ns*ns-ns*ns+xsub(i,2)*ns-ns+xsub(i,1)
enddo

r = size(V,2)

end subroutine

subroutine nufft3dIapp(nj,c,U,V,xxsub,ns,iflag,r,S,plan)
implicit none
integer  r,i,j,k,nj,ns,iflag,num
integer mm
integer xxsub(nj)
complex*16 M(nj,r),N(ns,ns,r),S(ns*ns*ns),c(nj),U(ns*ns*ns,r),V(nj,r)
complex*16 NN(ns*ns*ns,r)
double complex in1, out1
real*16  time_begin,time_end,countrage,countmax
dimension in1(ns,ns,ns), out1(ns,ns,ns)
integer*8 :: plan
integer FFTW_FORWARD,FFTW_MEASURE
parameter (FFTW_FORWARD=-1)
parameter (FFTW_MEASURE=0)


M = 0

do i = 1,nj
   do k = 1,r
      M(xsub(i),k) = M(xsub(i),k)+V(i,k)*c(i)
   enddo
enddo

do i = 1,r
   in1 = reshape(conjg(M(:,i)),ns,ns,ns)

   call dfftw_execute_dft(plan, in1, out1)
   N(:,i) = reshape(out1,ns*ns*ns,1)
enddo

S = (U*N,2)
end subroutine
