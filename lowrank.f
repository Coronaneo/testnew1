subroutine lowrankfac(M,tol,tR,mR,U,V)
implicit none
integer :: tR,mR,Nx,Ny,rs(tR),INFO
real :: tol
complex*16 TAU1(min(Nx,Ny)),
complex*16,allocatable :: M(:,:)


Nx = size(M,1)
Ny = size(M,2)

if (Nx==0 .or. Np==0)
   U = 0
   V = 0
   return
endif

if (tR<Np .and. tR<Nx)
c  get columns
   call randsample(Nx,tR,rs)
   M2 = M(rs,:)
   call zgeqrf(tR,Ny,M2,tR,TAU1,WORK1,-1,INFO) ![~,R2,E2]=qr(M2,0)
   where(abs(diag(R2))>tol*abs(R2(1))) Cidx = E2(1:tR)
   
c  get rows
   call randsample(Ny,3,cs)
   !cs = unique([cs' Cidx])
   M1 = M(:,cs)
   call qr(M1,R1,E1) ![~,R2,E2]=qr(M2,0)
   where(abs(diag(R1))>tol*abs(R1(1))) Ridx = E1(1:tR)

c  get columns again
   call randsample(Nx,3,rs)
   !rs = unique([rs' Ridx])
   M2 = M(rs,:)
   call qr(M2,R2,E2) ![~,R2,E2]=qr(M2,0)
   where(abs(diag(R2))>tol*abs(R2(1))) Cidx = E2(1:tR)

else
   Ridx = 1:Nx
   Cidx = 1:Ny
endif

c get rows
MR = M(Ridx,y)

c get middle matrix
MC = M(x,Cidx)

c get middle matrix
call qr(MC,QC)
call qr(transpose(conjg(MR)),QR)

if (tR<Ny .and. tR<Nx)
   call randsample(Ny,tR,cs)
   !cs = unique([cs' Cidx])
   call randsample(Nx,tR,rs)
   !rs = unique([rs' Ridx])
else
   cs = 1:Ny
   rs = 1:Nx
end

M1 = QC(rs,:)
M2 = QR(cs,:)
M3 = fun(x(rs,:),y(cs,:))
MD = pinv(M1)*(M3*oinv(transpose(conjg(M2))))
call svdtrunc(MD,mR,tol,U,S,V)
U = QC*U*sqrt(S)
V = QR*V*sqrt(S)
end subroutine



function diag(R)
   implicit none
   complex*16,allocatable :: R(:,:),diag(:)
   integer :: N,i
   N=size(R,1)
   do i = 1,N
      diag(i) = R(i,i)
   enddo
   return
end

subroutine svdtrunc(A,r,tol)

ms = min(size(A))
call gesvd(A,U,S,V)
if (ms>r)
   !idx = find(find(diag(S)>tol*S(1,1))<=r)
   U = U(:,idx)
   S = S(idx,idx)
   V = V(:,idx)
end
end subroutine




