function out = getboundaries(in,nx,ny,xbound)

%adds boundaries to values
out = zeros(nx+2,ny+2);
out(2:end-1,2:end-1) = in;

%if boundtype = 0, do nothing as its correct
if strcmpi(xbound,'periodic') %configuring periodic boundary for x
    out(1,:) = out(end-1,:);
    out(end,:) = out(2,:);
end
