function [u,v] = getvelocityfrompsi(psi,dx,dy,xbound)

% Solves velocities from psi on C-grid setup, includes boundary conditions
% for velocities

%dimensions
nx = size(psi,1);
ny = size(psi,2);

%U = -dPsi/dy
u = zeros(nx+1,ny);
u(2:end-1,:) = diff(psi,1,1)./dy;
if strcmpi(xbound,'periodic')
    u(1,:) = u(end-1,:);
    u(end,:) = u(2,:);
end

%V = dPsi/dx
v = zeros(nx,ny+1);
v(:,2:end-1) = diff(psi,1,2)./dx;
