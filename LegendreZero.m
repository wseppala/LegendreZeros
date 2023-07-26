%------------------------------------------------------%
%   Function to compute the n:th k-zero of the associated Legendre
% function with complex degree -1/2 + ik, integer order m
% and Re(x) > 0.

% The function is just a convenient wrapper for LegendreZeros.
%------------------------------------------------------%
function res = LegendreZero(m, n, x, int_len)
    vec = LegendreZeros(m, n, x, int_len);
    res = vec(end);
end 