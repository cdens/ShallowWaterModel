function [X,curresidual] = iterate_cgm(C,Xo,b,nx,ny,xbound,minresidual,maxiterations)

%pull A and b from C matrix
A = C(:,:,1:5);
b = b - squeeze(C(:,:,6));

%initialize X 
X = Xo; 

%get initial residual (r) 
ri = solve_residual(A,X,b,nx,ny,xbound);

%initial d = r
di = ri;

%begin looping towards solution
nn = 0;
curresidual = minresidual + 1;
while nn <= maxiterations && curresidual > minresidual
    nn = nn + 1;
    
    %1: solve for step size alpha_i = ri^T*ri/di^T*A*di
    riTri = sum(ri(:).^2);
    di_bound = getboundaries(di,nx,ny,xbound);
    Adi = NaN.*ones(nx,ny); %preallocating
    for i = 1:nx 
        ic = i + 1;
        for j = 1:ny
            jc = j + 1;
            Adi(i,j) = A(i,j,1)*di_bound(ic,jc) + A(i,j,2)*di_bound(ic+1,jc) + A(i,j,3)*di_bound(ic,jc+1) + ...
                A(i,j,4).*di_bound(ic-1,jc) + A(i,j,5).*di_bound(ic,jc-1);
        end
    end
    diTAdi = sum(di(:).*Adi(:)); 
    if riTri ~= 0 || diTAdi ~= 0
        alpha = riTri/diTAdi;
    else
        alpha = 0;
    end
    
    %2: iterate profile forward Xi+1 = Xi + alpha_i*di
    X = X + alpha.*di;
    
    %3: solve new residual ri+1 = b - A*Xi+1
    rip1 = solve_residual(A,X,b,nx,ny,xbound);
    
    %4: solve new step direction di+1 = ri+1 + di*beta = ri+1 + di*(ri+1^T*ri+1/ri^T*ri)
    dip1 = rip1 + di*sum(rip1(:).^2)/sum(ri(:).^2);
    
    %5: move new index (_i+1) to current index (_i) for r and d
    ri = rip1;
    di = dip1;
    
    %determine current residual magnitude
    curresidual = sum(abs(ri(:)));
    
end
