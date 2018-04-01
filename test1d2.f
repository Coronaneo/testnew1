	program test
	implicit none
        integer i,iflag,xsub(128),ier,num
        integer nj,ns,r
        real*16 begin1,end1
        integer*8  time_begin,time_end,countrage,countmax
        real*16 U1(12,128),V1(12,128),U2(12,128),V2(12,128)
        real*16 re1(128),re2(128),x(128),xsubb(128),pi,time1,time2
        real*16 arr(4)
        parameter (pi=3.141592653589793238462643383279502884197d0)
        complex*16 U(12,128),V(12,128),c(-64:63),S(128),re(128)
        complex*16 fk(128),Idx(128,128)
        real*8 x1(128),eps,error
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
        r=12
        iflag=1
        eps=1E-12
        num=10000
        open(unit = 10,file = 'Ur2.txt')
        read(10,*) U1
        open(unit = 20,file = 'Vr2.txt')
        read(20,*) V1
        open(unit = 10,file = 'Ui2.txt')
        read(10,*) U2
        open(unit = 10,file = 'Vi2.txt')
        read(10,*) V2
        call dfftw_plan_dft_1d(plan,nj,in1,out1,FFTW_FORWARD,0)
        print *,'start 1D type 2 testing:','nj  =',nj,'ns  =',ns
        print *,'eps             =',eps
        re=dcmplx(re1,re2)        
        U=dcmplx(U1,U2)
        V=dcmplx(V1,V2)
        V=conjg(V)
        !print *,V(1,:)
        !print *,U(1,:)
        do i = 1,128
           x(i) = i*pi/8
        enddo
        xsub=mod(floor(x+0.5),ns)+1
        do i = 1,128
           c(i-65) = exp(-dcmplx(0,1)*i/ns)
        enddo
        !print *,c(-64:-59)
        do i = 1,128
           x1(i) = i*pi*2*pi/(8*nj)
        enddo
        !print *,'ok'
        Idx=0
	do i = 1,nj
	   Idx(i,xsub(i))=1
	enddo

        call date_and_time(date,time,zone,values1)
        do i=1,num
        call nufft1dIIapp(nj,plan,c,U,V,xsub,ns,-1,r,S)
        enddo
        call date_and_time(date,time,zone,values2)
        !print *,'S(1:5)=',S(1:5)
        time1=sum((values2(5:8)-values1(5:8))*arr)
        print *,' T_our         = ',time1/num

        call date_and_time(date,time,zone,values1)
        do i=1,num
        call nufft1d2f90(nj,x1,fk,-1,eps,ns,c,ier)
        enddo
        call date_and_time(date,time,zone,values2)
        time2=sum((values2(5:8)-values1(5:8))*arr)
        !print *,'ier=',ier
        !print *,'fk(1:5)=',fk(1:5)
        print *,' T_nyu         = ',time2/num
        print *,' T_our/T_nyu   = ',time1/time2
        error=sqrt(real(sum((S-fk)*conjg(S-fk))/
     &  sum(S*conjg(S))))
        print *,' relative error= ',error
        call dfftw_destroy_plan(plan)
        

	end program

