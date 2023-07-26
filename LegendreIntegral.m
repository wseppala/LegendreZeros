function res1 = LegendreIntegral(nu,m,x)
%------------------------------------------------------%
%   Function to compute the value of the unscaled associated
% Legendre function with integer order 'm', complex degree 'nu' 
% and Re(x) > 0 to precision given by default tolerances of
% the 'integral'-method.

%   The method makes use of an integral form of the function,
% see e.g. NIST DLMF §14.12.7
%------------------------------------------------------%
res1 = integral(@(t) (x + (x^2-1)^(1/2) ...
    .*cos(t)).^(nu).*cos(m.*t), 0, pi);
end