function [P_span, spline] = b_spline_basis_functions(M,d,max_P)
% %basis function set up
% d = 2; 
knot_width = max_P/(M-d);
knots = -(d)*knot_width:knot_width:(M+d-1)*knot_width;

P_span = [knots(1):0.01:knots(end)]';
phi = zeros(length(P_span),M,d+1);
%zeroth order
for ii = 1:length(knots) -1; %M+d+1;
    aa = find(P_span>=knots(ii) & P_span<knots(ii+1));
    if ~isempty(aa)
    phi(aa,ii,1) = ones(length(aa),1);
    end
end
spline = phi(:,:,1);

% plot(P_span,spline)
%other order
for dd = 1:d
    for ii = 1:M+(d-dd)
        phi(:,ii,dd+1) = ((P_span - knots(ii))/(knots(ii+dd)-knots(ii))).*phi(:,ii,dd)+ ((knots(ii+dd+1)-P_span)/(knots(ii+dd+1)-knots(ii+1))).*phi(:,ii+1,dd);
    end
% spline = phi(:,:,(dd+1));
% plot(P_span,spline)
end

spline = phi(:,:,end);
aa = find(P_span>=0 & P_span<=max_P);
P_span = P_span(aa);
spline = spline(aa,1:M);
% plot(P_span,spline)

        
        