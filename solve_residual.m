function ri = solve_residual(A,X,b,nx,ny,xbound)

% X_bounded = getiterationprof(X,nx,ny);
X_bounded = getboundaries(X,nx,ny,xbound);
ri = zeros(nx,ny);

for i = 1:nx 
    ic = i + 1;
    for j = 1:ny
        jc = j + 1;
        ri(i,j) = -1*(A(i,j,1)*X_bounded(ic,jc) + A(i,j,2)*X_bounded(ic+1,jc) + A(i,j,3)*X_bounded(ic,jc+1) + ...
            A(i,j,4).*X_bounded(ic-1,jc) + A(i,j,5).*X_bounded(ic,jc-1) - b(i,j));
    end
end