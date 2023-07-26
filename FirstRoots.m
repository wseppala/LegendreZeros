%------------------------------------------------------%
%   Function to compute the first k-zeros of the associated
% Legendre function with complex degree -1/2 + i*k,
% integer order 'mus and main argument 's' for all orders
% 0...mus(end), ie. from order 0 until the last element in
% the input vector 'mus'. 

% NOTE: it is advised to only input orders in ascending order.

% The method in based on the first roots increasing with increasing
% order and does not take advantage of the theoretical lower limit for the
% first zeros, but moves iteratively from root to the next.
%------------------------------------------------------%
function res1 = FirstRoots(mus, x, int_len)
res1 = zeros(1,mus(end)+1);             % reserve space for result
mu = 0; start = 0;
tol = 10^-7;                            % default tolerance
inter = [start, start+int_len];
    while mu <= mus(end)
    f = chebfun(@(k) LegendreCC(-1/2 + 1i*k, mu, x, tol), inter);
    r = roots(f,'complex'); 
    r = sort(real(r));                  % project to reals and resort
        if isempty(r) && mu == 0        
            inter = inter + int_len;
            continue
        end
        if mu == 0                      % search interval after having 
            int_len = r(1);             % ... first root at mu = 0.
        end
        if isempty(r)                   % fail-safe for failed search
             f = chebfun(@(k) LegendreHyp(-1/2 + 1i*k, -mu, x), inter);
             f = projection(f,mu); r = roots(f);
             fprintf("first search failed at mu= %d\n",mu)
        end                             % dynamically check that roots fulfill sought-after form
        if mu > 2 && (r(1)-res1(mu)) > (res1(mu-1) - res1(mu-2))
            fprintf("Climb failed at mu= %d, no root matching expected pattern could be found\n",mu)
            break
        end
    mu = mu + 1;
    res1(mu) = real(r(1));
    inter = [r(1), r(1)+int_len];
    end
end