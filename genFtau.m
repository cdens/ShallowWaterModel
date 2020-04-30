function Ftau = genFtau(taux,tauy,dx,dy,rho_o,H)

if ismatrix(taux)
    if size(taux,1) == 1 && size(taux,2) == 1
        taux(1:nx,1:ny+1) = taux;
    elseif size(taux,1) ~= nx+1 || size(taux,2) ~= ny
        error('Invalid initial taux matrix')
    end
    
    taux2D = taux; clear taux;
    for i = 1:nt
        taux(:,:,i) = taux2D;
    end
end
if ismatrix(tauy)
    if size(tauy,1) == 1 && size(tauy,2) == 1
        tauy(1:nx+1,1:ny) = tauy;
    elseif size(tauy,1) ~= nx || size(tauy,2) ~= ny+1
        error('Invalid initial tauy matrix')
    end
    
    tauy2D = tauy; clear tauy;
    for i = 1:nt
        tauy(:,:,i) = tauy2D;
    end
end

%calculating stress forcing for all times
Ftau = (1/(rho_o*H)).*(diff(tauy,1,1)./dx - diff(taux,1,2)./dy);