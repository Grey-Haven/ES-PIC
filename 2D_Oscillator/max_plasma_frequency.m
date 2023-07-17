function w_pe = max_plasma_frequency(electrons,x,y,w_e,eps_0)
    dx = x(2)-x(1);
    dy = y(2)-y(1);
    q = electrons(1,5);
    m = electrons(1,6);
%     max_n_bar = -inf;
%     for i = 1:length(x)-1
%         for j = 1:length(y)-1
%             n_bar = size(find(electrons(:,1) >= x(i) & electrons(:,1) <= x(i+1) & electrons(:,2) >= y(j) & electrons(:,2) <= y(j+1)),1)/(dx*dy);
% %             disp(n_bar)
%             max_n_bar = max(max_n_bar,n_bar);
%         end
%     end
%     bins0 = histogram2(electrons(:,1),electrons(:,2),x,y);
    bins = histcounts2(electrons(:,1),electrons(:,2),x,y);
%     [C,I] = max(bins.Values)
%     max_n_bar0 = max(max(bins0.Values))/(dx*dy);
    max_n_bar = max(max(bins))/(dx*dy);
%     assert(max_n_bar0 == max_n_bar);
    w_pe = sqrt(max_n_bar*w_e*q^2/(m*eps_0));
end