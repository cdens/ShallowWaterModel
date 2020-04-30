function out = iterate2Dexplicit(in,nx,ny,C,xbound)

iter = getboundaries(in,nx,ny,xbound);

out = zeros(nx,ny);
for i = 1:nx
    ic = i + 1;
    
    for j = 1:ny %not solving for boundary points along y axis
        jc = j+1; %j correction for calling indices of P
        
        out(i,j) = C(i,j,1)*iter(ic,jc) + C(i,j,2).*iter(ic+1,jc) + C(i,j,3).*iter(ic,jc+1) + ...
            C(i,j,4).*iter(ic-1,jc) + C(i,j,5).*iter(ic,jc-1) + C(i,j,6);
    end
end





