function y = test_holm_bonferroni(PVAL)
%% used to counteract the problem of multiple comparisons
% PVAL = a Nx1 array, a series of p values.

% Order p values from smallest to greatest
[x,I] = sort(PVAL);

n = numel(PVAL);
alpha = 0.05;

c = zeros(size(I));
for rank = 1:n
    
    %% The significance level for this rank
    HB = alpha/(n-rank+1);
    
    %% Accept if pvalue <= HB
    if x(rank)<=HB
        c(rank) = 1;
    end
    
    y = zeros(size(PVAL));
    y(I(find(c==1)))=1;
end

%% y = a N*1 array, h 
