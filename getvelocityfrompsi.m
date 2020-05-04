function [u,v] = getvelocityfrompsi(psi,dx,dy,xbound,gridtype)

% Solves velocities from psi on C-grid, handle boundary conditions

%U = -dPsi/dy
%V = dPsi/dx

%C grid- U on E/W faces, V on N/S faces
%D grid- U on N/S faces, V on E/W faces

%dimensions
nx = size(psi,1);
ny = size(psi,2);

switch(lower(gridtype))
    case 'c'
%         uf = zeros(nx+2,ny+1);
%         uf(2:end-1,2:end-1) = -diff(psi,1,2)./dy;
%         if strcmpi(xbound,'periodic')
%             uf(1,:) = uf(end-1,:);
%             uf(end,:) = uf(2,:);
%         end
%         
%         vf = zeros(nx+1,ny+2);
%         vf(2:end-1,2:end-1) = diff(psi,1,1)./dx;
%         
%         %interpolating to grid centers
%         uc = 0.5.*(uf(1:end-1,:) + uf(2:end,:));
%         vc = 0.5.*(vf(:,1:end-1) + vf(:,2:end));
%         
%         %interpolating to grid faces
%         u = 0.5.*(uc(:,1:end-1) + uc(:,2:end));
%         v = 0.5.*(vc(1:end-1,:) + vc(2:end,:));

        uf = -diff(psi,1,2)./dy;
        vf = diff(psi,1,1)./dx;
        
        %interpolating to grid centers
        uc = 0.5.*(uf(1:end-1,:) + uf(2:end,:));
        vc = 0.5.*(vf(:,1:end-1) + vf(:,2:end));
        
        %interpolating to grid faces
        u = 0.5.*(uc(:,1:end-1) + uc(:,2:end));
        v = 0.5.*(vc(1:end-1,:) + vc(2:end,:));

    case 'd'

%         u = zeros(nx,ny+1);
%         u(:,2:end-1) = -1*diff(psi,1,2)./dy;
%         if strcmpi(xbound,'periodic')
%             u(1,:) = u(end-1,:);
%             u(end,:) = u(2,:);
%         end
% 
%         v = zeros(nx+1,ny);
%         v(2:end-1,:) = diff(psi,1,1)./dx;

        u = -1*diff(psi(2:end-1,:),1,2)./dy;
        v = diff(psi(:,2:end-1),1,1)./dx;
        
end