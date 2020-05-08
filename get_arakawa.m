function Jpz = get_arakawa(psi,zeta,nx,ny,dx,dy)

Jpz = zeros(nx,ny);
coeff = 1/(3*dx*dy);

% for i = 1:nx
%     ic = i + 1;
%     for j = 1:ny
%         jc = j + 1;
     
%         Jpz(i,j) = + zeta(ic-1,jc+1)*(2*psi(ic,jc+1) - psi(ic-1,jc) - psi(ic,jc-1)) + ... %z0
%                    + zeta(ic-1,jc)*(psi(ic-1,jc+1) - psi(ic-1,jc-1)) + ... %z1
%                    + zeta(ic-1,jc-1)*(psi(ic-1,jc) - psi(ic,jc-1)) + ... %z2
%                    + zeta(ic,jc-1)*(psi(ic-1,jc-1) - psi(ic+1,jc-1)) + ... %z3
%                    + zeta(ic+1,jc-1)*(psi(ic,jc-1) + psi(ic-1,jc) - 2*psi(ic+1,jc)) + ... %z4
%                    + zeta(ic,jc+1)*(psi(ic+1,jc+1) - psi(ic-1,jc+1)) + ... %z6
%                    + zeta(ic+1,jc+1)*(2*psi(ic+1,j) - 2*psi(ic,jc+1) - psi(ic-1,jc)) + ... %z7
%                    + zeta(ic+1,jc)*(psi(ic+1,jc-1) - psi(ic+1,jc+1)); %z8


for ic = 1:nx
    i = ic + 1;
    for jc = 1:ny
        j = jc + 1;
        Jpz(ic,jc)= coeff*( ...
            (psi(i,j-1)+psi(i+1,j-1)-psi(i,j+1)-psi(i+1,j+1))*zeta(i+1,j) ...
            - (psi(i-1,j-1)+psi(i,j-1)-psi(i-1,j+1)-psi(i,j+1))*zeta(i-1,j) ...
            + (psi(i+1,j)+psi(i+1,j+1)-psi(i-1,j)-psi(i-1,j+1))*zeta(i,j+1) ...
            - (psi(i+1,j-1)+psi(i+1,j)-psi(i-1,j-1)-psi(i-1,j))*zeta(i,j-1) ...
            + (psi(i+1,j)-psi(i,j+1))*zeta(i+1,j+1) ...
            - (psi(i,j-1)-psi(i-1,j))*zeta(i-1,j-1) ...
            + (psi(i,j+1)-psi(i-1,j))*zeta(i-1,j+1) ...
            - (psi(i+1,j)-psi(i,j-1))*zeta(i+1,j-1));
        
    end
end