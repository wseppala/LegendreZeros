function res1 = LegendreCC(nu,m,x,tol)
%------------------------------------------------------%
%   Function to compute the value of the unscaled associated
% Legendre function with integer order 'm', complex degree 'nu' 
% and Re(x) > 0 to an absolute precision of at least 'tol'.

%   The method makes use of an integral form of the function,
% see e.g. NIST DLMF §14.12.7

%   The quadrature method used is a nested form of Clenshaw-Curtis
% -quadrature with the lower grid size chosen according to 'm' 
% until a (default) maximum of 400.
%------------------------------------------------------%
f = @(t) (x + (x^2-1)^(1/2).*cos(t)).^(nu).*cos(m*t);
n = min(m+50,400); k = 2*n-1;
[xn0, wn0] = chebpts(n);                    % compute initial weights and nodes to be rescaled
[xk0, wk0] = chebpts(k);                    % ... from the interval [-1,1].
    function res2 = ccadpt(a, b, tol)
        xn = scaleNodes(xn0, [a b]);        % scale nodes and weights
        wn = scaleWeights(wn0, [a b]);
        Qn = wn*f(xn);                      % apply quadrature rule
        xk = scaleNodes(xk0, [a b]);
        wk = scaleWeights(wk0, [a b]);
        Qk = wk*f(xk);

        if abs(Qn - Qk) < tol               % approximation is good enough
            res2 = Qk;  
            return
        else                                % else recurse
            q1 = ccadpt(a, (a+b)/2, tol);
            q2 = ccadpt((a+b)/2, b, tol);
            res2 = q1+q2;
            return
        end

    end 
res1 = ccadpt(0, pi, tol);                  %Start recursion.
end

%   These functions are taken directly from the source code of
% chebfun/chebpts.
function y = scaleNodes(x, dom)
%SCALENODES   Scale the Chebyshev nodes X from [-1,1] to DOM.
% TODO: Deal with unbounded domains
if ( dom(1) == -1 && dom(2) == 1 )
    % Nodes are already on [-1, 1];
    y = x;
    return
end
% Scale the nodes:
y = dom(2)*(x + 1)/2 + dom(1)*(1 - x)/2;
end
function w = scaleWeights(w, dom)
%SCALEWEIGHTS   Scale the Chebyshev weights W from [-1,1] to DOM.
% TODO: Deal with unbounded domains
if ( dom(1) == -1 && dom(2) == 1 )
    % Nodes are already on [-1, 1];
    return
end
% Scale the weights:
w = (diff(dom)/2)*w;
end