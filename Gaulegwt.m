function [x, w]=Gaulegwt(x1, x2, n)
% Gaulegwt.m
%
% This script is for computing definite integrals using Legendre-Gauss 
% Quadrature. Computes the Legendre-Gauss nodes and weights  on an interval
% [x1, x2] with truncation order N
%
% Suppose you have a continuous function f(x) which is defined on [x1, x2]
% which you can evaluate at any x in [x1, x2]. Simply evaluate it at all of
% the values contained in the x vector to obtain a vector f. Then compute
% the definite integral using sum(f.*w);
%
% Written by Lateef Kareem - 12/02/2018
x = zeros(1,n); w = x;
EPS = 1.0e-14; % EPS is the relative precision.
m  = floor((n + 1)/2); % The roots are symmetric in the interval, so
xm = 0.5*(x2 + x1);
xl = 0.5*(x2 - x1);
for i = 0:m-1 % Loop over the desired roots.
    z = cos(3.141592654*(i+0.75)/(n+0.5));
    % Starting with this approximation to theith root, we enter the main loop of refinement by Newton? method.
    error = 1;
    while(error >EPS)
        p1=1.0; p2=0.0;
        for j = 0:n-1; % Loop up the recurrence relation to get the Legendre polynomial evaluated atz. 
            p3 = p2;
            p2 = p1;
            p1 = ((2.0*j+1.0)*z*p2-j*p3)/(j+1);
        end
        %p1 is now the desired Legendre polynomial. We next compute pp, its derivative,
        % by a standard relation involving alsop2, the polynomial of one lower order.
        pp=n*(z*p1-p2)/(z*z-1.0);
        z1=z;
        z=z1-p1/pp; %Newton's method.
        error = abs(z-z1);
    end
    x(i + 1) = xm - xl*z; % Scale the root to the desired interval,
    x(n-i)   = xm + xl*z; % and put in its symmetric counterpart.
    w(i + 1) =2.0*xl/((1.0-z*z)*pp*pp); %Compute the weight
    w(n-i)   =w(i + 1);% and its symmetric counterpart.
end