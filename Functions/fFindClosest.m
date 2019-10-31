function [Y,idx] = fFindClosest(X,V)
% X = [1,2,3,4,5];
% V = 0.1;

[~,idx] = min(abs(X-V));
Y = X(idx);
end