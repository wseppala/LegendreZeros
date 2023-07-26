%------------------------------------------------------%
%   Function to compute the first 'n' k-zeros of the associated Legendre
% function with complex degree -1/2 + ik, integer order mu
% and Re(x) > 0 for all orders contained in input vector 'mus'

%   Search parameter 'int_len' specifies the length of the search interval.
% The free parameter is advised to be chosen not too large, e.g. 10.

% Constraints on input: 
%   elements in 'mus' must be non-negative integers
%   Re(x) > 0 and scalar
%   n must be integer scalar and at least 1
%   int_len must be non-negative scalar floating point number.

%   Return value of function is a 'n' by 'length(mus)' matrix, with the k:th
% column corresponds to zeros of the k:th element of 'mus' and the k:th row
% corresponding to the k:th roots for each of the input orders in 'mus'.

%   The method used for function evaluation in interpolation can be varied;
% package 'LegendreCC' is default, but also 'LegendreIntegral' (fast) and 
% 'LegendreHyp' (slow) are available.

% This method is a variant of 'LegendreZeros' which does not take advantage
% of the theoretical lower limit, but which computes the beginning points
% for the search by using the 'FirstRoots'-method.
%------------------------------------------------------%
function zero = LegendreZerosF(mus, n, x, int_len)
arguments
    mus {mustBeVector, mustBeInteger}
    n (1,1) {mustBeInteger, mustBePositive}
    x (1,1) 
    int_len (1,1) {mustBePositive}
end
q = int_len;
res = zeros(n,length(mus)); 

climbedroots = FirstRoots(mus, x, int_len);
begpoints = climbedroots( end-length(mus)+1 : end) - 0.1;

for i = 1:length(mus)
    mu = mus(i);
    foundcount = 0;
    inter = [begpoints(i) begpoints(i)+q];
    while foundcount < n
        f = chebfun(@(k) LegendreCC(-1/2 + 1i*k, mu, x, 10^-9), inter);
        r = roots(f,'complex'); 
        r = sort(real(r));
        ind1 = inter(1) < r; ind2 = r <= inter(2);
        r = r(ind1 & ind2);     

        newcount = length(r);        
        needed  = n - foundcount;   
        addcount = min(needed, newcount);   
        if addcount > 0
            res((1+foundcount):(foundcount+addcount),i) = r(1:addcount);
            foundcount = foundcount + addcount;
        end
        inter = inter + q;
    end 
end
zero = res;
end 