function [u,v] = getvelocityfrompsi(psi,dx,dy,gridtype)

% Solves velocities from psi on C-grid, handle boundary conditions

%U = -dPsi/dy
%V = dPsi/dx

%C grid- U on E/W faces, V on N/S faces
%D grid- U on N/S faces, V on E/W faces

%dimensions
nx = size(psi,1);
ny = size(psi,2);

switch(lower(gridtype))

    case 'a'

        u = -0.5.*(psi(2:end-1,3:end) - psi(2:end-1,1:end-2))./dy;
        v = 0.5.*(psi(3:end,2:end-1) - psi(1:end-2,2:end-1))./dx;
        
        
    case 'c'

        uf = -diff(psi,1,2)./dy;
        vf = diff(psi,1,1)./dx;
        
        %interpolating to grid centers
        uc = 0.5.*(uf(1:end-1,:) + uf(2:end,:));
        vc = 0.5.*(vf(:,1:end-1) + vf(:,2:end));
        
        %interpolating to grid faces
        u = 0.5.*(uc(:,1:end-1) + uc(:,2:end));
        v = 0.5.*(vc(1:end-1,:) + vc(2:end,:));
        

    case 'd'

        u = -1*diff(psi(2:end-1,:),1,2)./dy;
        v = diff(psi(:,2:end-1),1,1)./dx;
        
end