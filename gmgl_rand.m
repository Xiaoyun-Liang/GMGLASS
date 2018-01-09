% Reference: The Matlab code is based on the following paper:
%            Xiaoyun Liang, David N. Vaughan, Alan Connelly, Fernando Calamante. A novel
%            group-fused sparse partial correlation method for simultaneous
%            estimation of functional networks.in group comparison studies. Brain Topography, 12/2017; DOI:10.1007/s10548-017-0615-6.
%
% Copyright 2017 Florey Institute of Neuroscience and Mental Health, Melbourne, Australia
% Written by Xiaoyun Liang (Email: x.liang@brain.org.au)
% This software is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied  
% warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 


%This implementation is adapted from the follwing paper:
%Yang S, Lu Z, Shen X, Wonka P, Ye J (2015) Fused multiple graphical lasso SIAM J Optim 25:916–943.

function [P,funVal] = gmgl_rand(S,para_struc)
%% second order method for fused multiple graphical lasso
% This function solves the following problem:
% min \sum{-log det(P(:,:,k) + tr(P(:,:,k)*S(:,:,k))} + alpha*||P||_1
%      + beta*\sum||P(:,:,k)-P(:,:,k+1)||_1
% Input:
%     S    - emprical covariance matrices
%     para_struc:
%          - alpha:  hyper-parameter - l1 regularization parameter
%          - beta:   hyper-parameter - fused regularization parameter
%          - maxlineiter: the maximum number of line search
%          - Newtontol:  tolerence 
%          - maxiter:   the maximum iteration number
%          - SPGreltol: the tolerence for SPG
%          - SPGmaxiter: the maxinum iteration number for SPG
%          - sigma:     the line search parameter
global marks

if nargin < 2
  error('para_struc.alpha and para_struc.beta should be specified!');
end

if isfield(para_struc,'alpha')
    alpha = para_struc.alpha;
else
   error('para_struc.beta should be specified!');
end

if isfield(para_struc,'beta')
    beta = para_struc.beta;
else
   error('para_struc.beta should be specified!');
end

if ~isfield(para_struc,'maxlineiter')
    para_struc.maxlineiter = 50;
end
if ~isfield(para_struc,'SPGreltol')
    para_struc.SPGreltol = 1e-6;
end
if ~isfield(para_struc, 'SPGmaxiter')
    para_struc.SPGmaxiter = 20;
end
if ~isfield(para_struc,'Newtontol')
   para_struc.Newtontol = 1e-6;
end

if ~isfield(para_struc,'maxiter')
   para_struc.maxiter = 1500;
end

if ~isfield(para_struc,'sigma')
    para_struc.sigma = 0.001;
end



[n,n,K] = size(S);
%randam permutation
order=randperm(K/2);
[a, rank]=sort(order);
S=S(:,:,[order K/2+1:K]);

if isfield(para_struc, 'P0')
    P = para_struc.P0;
else
    P = zeros(n,n,K);
    for k = 1:K
        P(:,:,k) = diag(1./diag(S(:,:,k)));
    end
end

marks = zeros(n,n,K);
for k = 1:K
    marks(:,:,k) = 1 - eye(n,n);
end

funVal = [];
R = zeros(n,n,K);
B = R;
for k = 1:K
    [R(:,:,k),info] = chol(P(:,:,k));
    if info > 0
        error('%d not positive definite', k);
    end
    invR = inv(R(:,:,k));
    B(:,:,k) = invR*invR';  %inv(P)
end

[fobjp, l1normX] = computLogDet(P, S, R, alpha, beta);
funVal(1) = fobjp;
fx = zeros((n-1)*n/2,1);
fy = zeros((n-1)*n/2,1);
tev = zeros(K,1);
iter = 1;
while (1)
    iter = iter+1;
    G = S - B;
    numF = findFreeSet(G,P,fx,fy,alpha,beta,tev);
    para_struc.NumF = numF;
    D = fmgl_spg(S,P*1.0,para_struc,B,fx,fy,max(para_struc.SPGreltol,0));

    if (norm(D(:),'inf') < 5e-5)
        break;
    end
    
   %linesearch
    LS_alpha = 1;
    PW = (P+D).*marks;
    l1normXD = sum(abs(PW(:)));
    

     PW = PW(:,:,1:K/2) - PW(:,:,K/2+1:K);   



    l1normXD = alpha*l1normXD + beta*sum(abs(PW(:)));   
    trdg = 0;
    
    for k = 1:K
        t1 = G(:,:,k);
        t2 = D(:,:,k);
        trdg = trdg + t1(:)'*t2(:);
    end
    
    fobjp1 = 1e15;
    for liter = 1 : para_struc.maxlineiter
        W = P+LS_alpha*D;
        flag = 0;
        for k = 1:K
            [Tz,info] = chol(W(:,:,k));
            if info > 0
                flag = 1;
                break;
            end
            R(:,:,k) = Tz;
        end
        
        if flag == 1
            LS_alpha = LS_alpha/2;
            continue;
        end
        
        [fobj, l1normX1] = computLogDet(W, S, R, alpha, beta);
        
        if fobj <= fobjp + LS_alpha*para_struc.sigma*(trdg + l1normXD - l1normX);
            l1normX = l1normX1;
            break;
        end
        
        if fobjp1 < fobj
            l1normX = l1normX1;
            break
        end
        fobjp1 = fobj;
        LS_alpha = LS_alpha/2;
    end
    P = W;
    
    funVal(iter) = fobj;
    
    fobjp = fobj;
    if iter == 1
        for k = 1:K
            invR = inv(R(:,:,k));
            B(:,:,k) = invR*invR';  %inv(P)
        end
        continue
    end
    
    if(abs(funVal(iter) - funVal(iter-1)) < para_struc.Newtontol*abs(funVal(iter-1))) || norm(D(:),'inf') < 5e-5
        break;
    end
    
    if iter>para_struc.maxiter
        break;
    end
    
    for k = 1:K
        invR = inv(R(:,:,k));
        B(:,:,k) = invR*invR';  %inv(P)
    end
end
    P(:,:,1:K/2)= P(:,:,rank);
end

function [fobj, l1normX] = computLogDet(W, S, R, alpha, beta)
global marks
[n,n,K] = size(S);
logdet = 0;
trdg = 0;
for k = 1:K
    logdet = logdet + 2*sum(log(diag(R(:,:,k))));
    t1 = S(:,:,k);
    t2 = W(:,:,k);
    trdg = trdg + t1(:)'*t2(:);
end
W = W.*marks;
l1normX = sum(abs(W(:)));

PW = W(:,:,1:K/2) - W(:,:,K/2+1:K);   % 26/02/15

l1normX = alpha*l1normX + beta*sum(abs(PW(:)));  

fobj = trdg - logdet + l1normX;
end


function opterr = objcompute(D,hatP, B, S, K, alpha, beta)
global marks
hatP = hatP.*marks;
opterr = 0;
for k = 1:K
    Wz = B(:,:,k)*D(:,:,k);
    opterr = opterr + Wz(:)'*Wz(:)/2 + ...
        trace((S(:,:,k) - B(:,:,k))*D(:,:,k));
end
opterr = opterr + alpha*sum(abs(hatP(:)));
% Temp = hatP(:,:,1:K-1) - hatP(:,:,2:K);   
Temp = hatP(:,:,1:K/2) - hatP(:,:,K/2+1:K); % 26/02/15

opterr = opterr + beta*sum(abs(Temp(:)));
end

function [W,funVal,lambdak] = fmgl_spg(S,P,para_struc,B,fx,fy,termtol)
lambdamin = 1e-30;
lambdamax = 1e5;
gamma = 1e-4;
eta = 2;
M = 10;
f_1 = -inf*ones(M,1);

K = size(S,3);
l1reg = para_struc.alpha;
freg = para_struc.beta;

SB = S - 2*B;
grad = SB;
hatP = 1.0*P; P_start = 1.0*P;
iter = 1;
for k = 1:K
    grad(:,:,k) = B(:,:,k)*hatP(:,:,k)*B(:,:,k);
end
grad  = grad + SB;
z = hatP - grad;
fmgl_subfusedLasso(P,z,fx,fy,l1reg, freg, para_struc.NumF);
f_1(iter) = objcompute(hatP-P_start,hatP, B, S, K, l1reg, freg);
funVal(iter) = f_1(iter);
pomga = P - hatP;
err = norm(pomga(:),'inf');
if err <= termtol
    W = P - P_start;
    return;
end

if isfield(para_struc,'SPGlambdak')
    lambdak =  para_struc.SPGlambdak;
else
    lambdak = 1/(err + 1e-15);
end
P = 1.0*hatP;
iter = iter+1;

while (iter < para_struc.SPGmaxiter)
    
    cont_inner = 1;
    while cont_inner
        z = hatP - 1/lambdak*grad;
        fmgl_subfusedLasso(P,z,fx,fy,l1reg/lambdak, freg/lambdak, para_struc.NumF);
        D = P - hatP;
        fnew = objcompute(P-P_start,P, B, S, K, l1reg, freg);
        ddr = D(:)'*D(:);
        f_threshold = max(f_1) - 0.5*gamma*lambdak*ddr;
        
        
        if fnew <= f_threshold
            cont_inner=0;
        else
            lambdak = eta*lambdak;
            
        end
    end
    f_1(mod(iter,M)+1) = fnew;
    funVal(iter) = fnew;
	iter = iter+1;
    hatgrad = grad;
    for k = 1:K
        grad(:,:,k) = B(:,:,k)*P(:,:,k)*B(:,:,k);
    end
    grad  = grad + SB;
    gk = grad - hatgrad;
    td = D(:)'*gk(:);
    if td <=0
        lambdak = lambdamax;
        continue
    end    
    td = td/ddr;
    lambdak = max(lambdamin,min(td,lambdamax));
    
    pomga = P - hatP;
    err = norm(pomga(:),'inf')/norm(P(:),'inf');
    
    if err <= termtol
        break;
    end
    
    hatP = 1.0*P;   
end
W = P - P_start;
end    

