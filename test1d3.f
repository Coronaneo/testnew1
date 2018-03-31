	program test
	implicit none
        integer i,iflag,xsub(128),ier,num,ksub(128)
        integer nj,ns,r,nk
        real*16 begin1,end1
        integer*8  time_begin,time_end,countrage,countmax
        real*16 U1(57,128),V1(57,128),U2(57,128),V2(57,128)
        real*16 re1(128),re2(128),x(128),xsubb(128),time1,time2
        real*16 arr(4)
        real*8 pi
        parameter (pi=3.141592653589793238462643383279502884197d0)
        complex*16 U(128,57),V(128,57),c(128),S(128),re(128)
        complex*16 fk(128),Idx(128,128),Idk(128,128)
        real*8 x1(128),eps,error,k(128)
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
        r=57
        iflag=-1
        eps=1d-12
        num=10000
        open(unit = 10,file = 'Ur3.txt')
        read(10,*) U1
        open(unit = 20,file = 'Vr3.txt')
        read(20,*) V1
        open(unit = 10,file = 'Ui3.txt')
        read(10,*) U2
        open(unit = 10,file = 'Vi3.txt')
        read(10,*) V2
        open(unit = 10,file = 'Rer.txt')
        read(10,*) re1
        open(unit = 10,file = 'Rei.txt')
        read(10,*) re2
        call dfftw_plan_dft_1d(plan,nj,in1,out1,FFTW_FORWARD,0)
        print *,'start 1D type 3 testing:','nj  =',nj,'ns  =',ns
        print *,'eps            =',eps
      
        U=dcmplx(transpose(U1),transpose(U2))
        U=conjg(U)
        V=dcmplx(transpose(V1),transpose(V2))
        re=dcmplx(re1,re2)
        !print *,'V(1,1:5)=',V(1,1:5)
        !print *,'U(1,1:5)=',U(1,1:5)
        do i = 1,128
           x(i) = dsin(-pi*i/nj)/2*nj
        enddo
        xsub=mod(floor(x+0.5)+ns,ns)+1
        !print *,'xsub(1:5)=',xsub(65:74)
        do i = 1,128
           k(i) = 48*dcos(i*pi/ns)
        enddo
        ksub=mod(floor(k+0.5)+ns,ns)+1
        !print *,'ksub(1:5)=',ksub(65:74)
        do i = 1,128
           c(i) = dcmplx( dcos(pi*i/nj), -dsin(pi*i/nj))
        enddo
        !print *,'c =',c(-64:-59)
        do i = 1,128
           x1(i) = pi * dsin(-pi*i/nj)
        enddo
        !print *,'ok'
        Idx=0
	do i = 1,nj
	   Idx(xsub(i),i)=1
	enddo
        do i = 1,nk
           Idk(i,ksub(i))=1
        enddo

        call date_and_time(date,time,zone,values1)
        do i=1,num
        call nufft1dIIIapp(nj,nk,plan,c,U,V,xsub,ksub,ns,-1,r,S)
        enddo
        call date_and_time(date,time,zone,values2)
        !print *,'S(1:5)=',S(65:75)
        time1=sum((values2(5:8)-values1(5:8))*arr)
        print *,' T_our         = ',time1/num

        call date_and_time(date,time,zone,values1)
        do i=1,num
           call nufft1d3f90(nj,x1,c,iflag,1d-12,nk,k,fk,ier)
        enddo
        call date_and_time(date,time,zone,values2)
        time2=sum((values2(5:8)-values1(5:8))*arr)
     
        !print *,'ier=',ier
        !print *,'fk(1:5)=',fk(1:5)
        print *,' T_nyu         = ',time2/num
        print *,' T_our/T_nyu   = ',time1/time2
c        do i = 1,128
c           print *,i,re(i)-fk(i)
c        enddo
c        error=sqrt(real(sum((S-re)*conjg(S-re))/
c     &  sum(re*conjg(re))))
c        print *,' relative error1= ',error
c        error=sqrt(real(sum((fk-re)*conjg(fk-re))/
c     &  sum(re*conjg(re))))
c        print *,' relative error2= ',error
        error=sqrt(real(sum((S-fk)*conjg(S-fk))/
     &  sum(S*conjg(S))))
        print *,' relative error= ',error
        call dfftw_destroy_plan(plan)
        

	end program
