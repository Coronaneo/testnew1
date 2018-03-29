	program test
	implicit none
        integer i,iflag,xsub(128),ier,num,ksub(128)
        integer nj,ns,r,nk
        real*16 begin1,end1
        integer*8  time_begin,time_end,countrage,countmax
        real*16 U1(50,128),V1(50,128),U2(50,128),V2(50,128)
        real*16 re1(128),re2(128),x(128),xsubb(128),pi,time1,time2
        real*16 arr(4)
        parameter (pi=3.141592653589793d0)
        complex*16 U(128,50),V(128,50),c(128),S(128),re(128)
        complex*16 fk(128),Idx(128,128)
        real*8 k(128),eps
        real*8 error,x1(128)
        double complex in1, out1
        dimension in1(128), out1(128)
	integer*8 :: plan
        integer FFTW_FORWARD,FFTW_MEASURE
        parameter (FFTW_FORWARD=-1)
        parameter (FFTW_MEASURE=0)
    
        character*8 date
        character*10 time
        character*5 zone 
        integer*4 values1(8),values2(8)
        
        arr(1)=3600
        arr(2)=60
        arr(3)=1
        arr(4)=0.001
        nj=128
        ns=128
        nk=128
        r=50
        iflag=-1
        eps=1d-12
        num=1

        open(unit = 10,file = 'Rer.txt')
        read(10,*) re1
        open(unit = 10,file = 'Rei.txt')
        read(10,*) re2

        print *,'start 1D type 3 testing:','nj  =',nj,'ns  =',ns
        print *,'eps            =',eps
       
        U=dcmplx(transpose(U1),transpose(U2))
        V=dcmplx(transpose(V1),transpose(V2))
        re=dcmplx(re1,re2)


        do i = 1,128
           k(i) = (i-65)*3.0/7.0
        enddo

        print *,'k(1:5)=',k(66:75)
        do i = 1,128
           c(i) = exp(-dcmplx(0,1)*i/ns)
        enddo
        !print *,'c =',c(-64:-59)
        do i = 1,128
           x1(i) = i*pi*2*pi/(8*nj)
        enddo




        call date_and_time(date,time,zone,values1)
        do i=1,num
        call nufft1d3f90(nj,x1,c,iflag,eps,nk,k,fk,ier)
        enddo
        call date_and_time(date,time,zone,values2)
        time2=sum((values2(5:8)-values1(5:8))*arr)

        !print *,'ier=',ier
        !print *,'fk(1:5)=',fk(65:75)
        print *,' T_nyu         = ',time2/num

        do i = 1,128
           print *,i,re(i)-fk(i)
        enddo

        error=sqrt(real(sum((fk(1:128)-re)*conjg(fk(1:128)-re))/
     &  sum(re*conjg(re))))
        print *,' relative error2= ',error

        

	end program
