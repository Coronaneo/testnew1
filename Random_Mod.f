	Module Random_Mod

	  Implicit None

	  !// 使用内部函数获得int64的Kind值

	  Integer , parameter :: int64 = Selected_Int_Kind( 10 )

	contains
	  Subroutine Init_Random_Seed()   

	    integer :: ised , i , pid

	    integer(int64) :: t

	    integer , allocatable :: sed(:)

	    call random_seed( size = ised ) !// 获得种子大小

	    allocate( sed(ised) ) !// 分配种子

	    call system_clock(t) !// 获得时间

	    pid = getpid() !// 获得处理器ID

	    t = ieor(t, int(pid, kind(t))) !// 用 pid 和日期做XOR运算

	    do i = 1, ised

		sed(i) = lcg(t) !// 用线性同余计算种子

	    end do

	    call random_seed( put=sed ) !// 给定种子 

	  End Subroutine Init_Random_Seed 
	  Function lcg(s) !// 线性同余算法 Linear congruential generator

	    integer :: lcg

	    integer(int64) :: s

	    if (s == 0) then

	       s = 104729

	    else

	       s = mod(s, 4294967296_int64)

	    end if

	    s = mod(s * 279470273_int64, 4294967291_int64)

	    lcg = int(mod(s, int(huge(0), int64)), kind(0))

	  End Function lcg
	End Module Random_Mod

	subroutine randsample(N,t,r)
	use Random_Mod

	integer :: N,t
	integer r(t)
        real x(t)


	call Init_Random_Seed()
	call random_number(x)
	r = floor(x*N)+1
	return
	end

       
c	Program www_fcode_cn

c	  use Random_Mod

c	  Real :: x(8)

c	  call Init_Random_Seed()

c	  call random_number(x)

c	  write(*,*) x

c	End Program www_fcode_cn
