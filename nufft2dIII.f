subroutine nufft2dIII(nj,k,x,nk,iflag,ns,rt,tol,U,V,kksub,r)
implicit none

integer :: ns,rt,iflag,nj,j,r,xsub(nj),i,nk
real  :: tol
real*8 pi,x(nj),k(nk)
parameter (pi=3.141592653589793238462643383279502884197d0)
complex*16 fftconst,U(nk,r),V(nj,r)
complex*16 M(nk,nj)

fftconst = iflag*dcmplx(0,1)/ns*2*pi

do i = 1,nk
   do j = 1,nj
      M(i,j) = exp(fftconst*(x(j,1)*(k(i,1)-floor(k(i,1)+0.5))+
&     x(j,2)*(k(i,2)-floor(k(i,2)+0.5))))
   enddo
enddo

call lowrankfac(M,tol,rt,rt,U,V)

ksub = mod(floor(k+0.5),ns)+1
do i = 1,nk
   kksub(i) = ksub(i,2)*ns-ns+ksub(i,1)
enddo

r = size(V,2)

end subroutine

subroutine nufft1dIIIapp(nj,nk,c,U,V,kksub,ns,iflag,r,S,plan)
integer  r,i,j,k,nj,ns,iflag,num,nk
integer mm
integer xsub(nj),Ksub(nk,nj)
complex*16 M(nj,r),N(nj,r),S(nk),c(nj),U(nk,r),V(nj,r)
complex*16 NN(nj,r)
double complex in1, out1
real*16  time_begin,time_end,countrage,countmax
dimension in1(ns), out1(nj)
integer*8 :: plan
integer FFTW_FORWARD,FFTW_MEASURE
parameter (FFTW_FORWARD=-1)
parameter (FFTW_MEASURE=0)

M=0

do i = 1,nj
   do k = 1,r
      M(i,k) = V(i,k)*c(i)
   enddo
enddo

end subroutine

