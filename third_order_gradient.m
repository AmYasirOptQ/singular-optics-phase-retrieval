function [dudx, dudy] = third_order_gradient(u, dx)
% THIRD_ORDER_GRADIENT  Sixth-order central finite-difference gradient.
%
%   [dudx, dudy] = third_order_gradient(u, dx)
%
%   Stencil (O(dx^6)):
%     du/dx|_j = [ (3/4)(u_{j+1}-u_{j-1}) - (3/20)(u_{j+2}-u_{j-2})
%                  + (1/60)(u_{j+3}-u_{j-3}) ] / dx
%
%   Boundary columns/rows (within 3 pixels of edge) are left as zero.

N    = size(u, 1);
dudx = zeros(N);
dudy = zeros(N);

for j = 4 : N-3
    dudx(:,j) = ( (3/4) *(u(:,j+1) - u(:,j-1)) ...
                - (3/20)*(u(:,j+2) - u(:,j-2)) ...
                + (1/60)*(u(:,j+3) - u(:,j-3)) ) / dx;

    dudy(j,:) = ( (3/4) *(u(j+1,:) - u(j-1,:)) ...
                - (3/20)*(u(j+2,:) - u(j-2,:)) ...
                + (1/60)*(u(j+3,:) - u(j-3,:)) ) / dx;
end
end
