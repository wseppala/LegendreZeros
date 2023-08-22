%------------------------------------------------------%
%   Function to compute the n:th k-zeros for a fixed main argument
% z > 1 among all orders m = 0, 1, 2, ... . So whereas 'LegendreZeros'
% computes the n:th zero for fixed z > 1 AND fixed order m, this function
% finds the n:th roots with only the main argument z fixed. 

%   The return value is a struct with fields:
% - orders: the orders of the Legendre function at given zeros
% - placement_within_order: this positive integer tells that if you
%   fix the order, then the returned zero is the
%   placement_within_order:th zero. 
% - value: the k-value of the zero.
%   **  the k:th element in each of these fields corresponds to the k:th
%       element in the input 'n'.

%   On completion, the function also prints the figure showing where the
% sought-after root is situated (marked with a black circle) in the 
% grid when the zeros are plotted as functions of order m. This is only
% done if the input is scalar.

%   Constraints on input: 
% - input 'n' must be of the form 'a:b' (linearly spaced integer vector)
%   with a >=1 and b >=a.
% - input 'z' must be a positive real scalar

%   The method used for function evaluation in interpolation can be varied;
% package 'LegendreCC' is default, but also 'LegendreIntegral' (fast) and 
% 'LegendreHyp' (slow) are available. 
% 
%   The default search interval length is 8 units, but this is naturally 
% variable, though smallish values are
% strongly recommended.
%------------------------------------------------------%
function returnedstruct = NthAbsoluteZero(n,z)
max_det_n = 0; max_order = 0;
n_inc = 3; % we are free to choose how much we wish to increment n by if 
               % we haven't found enough roots yet
while(max_det_n < max(n))
    n_inc = n_inc + 2;
    abs_max_root = LegendreZero(0,n_inc, z ,7);
    while(LegendreZero(max_order,1, z, 8) < abs_max_root )
            max_order = max_order + 1;
    end
    range_restrict = LegendreZeros(0:max_order, n_inc,z,8);
    max_det_n = sum(range_restrict <= abs_max_root,"all");
end
finaldata = LegendreZeros(0:max_order,n_inc,z,8);
sorted = sort(finaldata(finaldata <= abs_max_root));

returnval = sorted(n);
placements = []; orders = [];
for i = 1:length(n)
[t,indexx] = ismember(returnval(i),finaldata);
if t ==1
    [placement, order] = ind2sub([height(finaldata),length(finaldata)],indexx);
    
end
placements = [placements, placement]; orders = [orders, order-1];
end
returnedstruct = struct('orders',orders,'placement_within_orders',placements,'values',returnval);
% Comment out this last if-statement you don't wish to see the plots after function calls. 
if(length(n)==1)
    figure, hold on
    scatter(0:max_order,finaldata,15,'filled')
    scatter(order-1, returnval,50, 'black')
    plot(0:max_order,ones(max_order+1).*abs_max_root,'--')
    grid on, axis tight
    title("Circled n:th root, with n="+n)
    xlabel("Order, m"), ylabel("k-zero, k_{nm}")
    hold off
end
end
