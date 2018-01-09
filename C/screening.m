%% screening for fused multiple graphical lasso
% 
%  Input: 
%       S      -  emprical covariance matrices (N x N x K)
%       lambda -  l1 regularization parameter
%       rho    -  fused regularization parameter
%  Output:
%       adj    -  the constructed adjacency matrix
%  
%  How to use?
%  
%  adj = screening(S,lambda,rho);
% [SC, ID] = graphconncomp(sparse(adj));
% %[SC,ID] includes the indices of connected components