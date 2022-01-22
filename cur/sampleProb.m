function T = sampleProb(U, pr, k)
%   U = universal set = [n*d] of n items, each of size d
%   pr = [n*1] probability of choosing each of the n items
%   k = [int] sample size
%   T =  [k*d] sample
    cumpr = [0;reshape(cumsum(pr),[],1)];
    idx = sum(cumpr < rand(1,k), 1);
    T = U(idx,:);
end