%------------------------------------------------------%
%   Function to compute the n:th k-zero for a fixed main argument
% z > 1 among all orders m = 0, 1, 2, ... . So whereas 'LegendreZeros'
% computes the n:th zero for fixed z > 1 AND fixed order m, this function
% finds the n:th root with only the main argument z fixed. 

%   The return value is a struct with fields:
% - order: the order of the Legendre function at that zero
% - placement_within_order: this positive integer tells that if you
%   fix the order, then the returned zero is the
%   placement_within_order:th zero. 
% - value: the k-value of the zero.

%   On completion, the function also prints the figure showing where the
% sought-after root is situated (marked with a black circle) in the 
% grid when the zeros are plotted as functions of order m.

%   Constraints on input: 
% - input 'n' must be a positive integer scalar
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
while(max_det_n < n)
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
[t,indexx] = ismember(returnval,finaldata);
if t ==1
    [placement, order] = ind2sub([height(finaldata),length(finaldata)],indexx);
end
returnedstruct = struct('order',order,'placement_within_order',placement,'value',returnval);
figure, hold on
scatter(0:max_order,finaldata,15,'filled')
scatter(order-1, returnval,50, 'black')
plot(0:max_order,ones(max_order+1).*abs_max_root,'--')
grid on, axis tight
title("Circled n:th root, with n="+n)
xlabel("Order, m"), ylabel("k-zero, k_{nm}")
hold off
end