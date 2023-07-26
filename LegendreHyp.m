function a = LegendreHyp(l,m,z)
%-----------------------------------------------%
%   Function to compute the value of the hypergeometric 
% function present in the definition, see e.g. NIST DLMF 
% §14.3.6. with no parameter restrictions except that 'l'
% cannot be negative integer, l != -1, -2, -3, ... .

%   For integer order 'l', this method is much slower than 
% the methods in this package based on quadrature methods of 
% integral representations. 
    a = hypergeom([-l, l+1], 1+m , 1/2*(1-z));
end