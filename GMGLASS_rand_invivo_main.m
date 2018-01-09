%%%%%%%%  Main program
%        Simultaneous estimation of functional networks from 2 groups
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input data:(1) Mean regional tiem series from two groups of subejcts;
%            (2) Time series should be normalized;
%            (3) For each subject, the dataset is equally divided into n
%            blocks (n=20), with each block having 10 (200/20) time points; 

% Output data:(1) S: Stable matrix across 200 subsamples (p*p*lambda*rho)
%            (2) Y: Estimated networks at group-level (size: p*p*w*lambda*rho)
%            (3) I: Estimated networks at individual-level (size:p*p*k*w*lambda*rho):I1PT, I2NC



%Note:(1) Alpha and beta can be adjusted;
%         
%     (2) The ranges of alpha and beta can generally be chosen following 2 empirical observations: 
%           (a) Parameters chosen can achieve as dense networks as possible to potentially include all true edges
%           (b) Parameters chosen should also achieve as sparse networks as possible, but avoid empty networks.
%     (3) Stability selection is employed by subsampling data 100 times.
%     (4) Average number of selected connections across the set /\ needs to be estimated at both group- and individual-level, which is then
%         employed to calculate Pthr.
%      
% Reference: The Matlab code is based on the following paper:
%            Xiaoyun Liang, David N. Vaughan, Alan Connelly, Fernando Calamante. A novel
%            group-fused sparse partial correlation method for simultaneous
%            estimation of functional networks.in group comparison studies. Brain Topography, 12/2017; DOI:10.1007/s10548-017-0615-6.
%            The fused multiple graphical lasso problem is solved by employing a second-order method proposed by Yang et. al.: 
%            Yang S, Lu Z, Shen X, Wonka P, Ye J (2015) Fused multiple.
%            
%            graphical lasso SIAM J Optim 25:916–943
% Copyright 2017 Florey Institute of Neuroscience and Mental Health, Melbourne, Australia
% Written by Xiaoyun Liang (Email: x.liang@brain.org.au)
% This software is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied  
% warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 

%%
% Note: In those places with .../, please insert your specified directory
% of either input data or output results.
clc
clear all
addpath('C');
addpath('Matlab');
K = 26;
n_array = [90];

para_struc.alpha=0;
para_struc.beta=0;
para_struc.sigma = 1e-3;
para_struc.maxlineiter = 100;
para_struc.maxiter = 100;
para_struc.SPGmaxiter = 100;


n = n_array(1);
A = randn(n,n);
S = zeros(n,n,K);
%%%%%%%%%%%%
P =  zeros(n,n,K);
trueSparsity  = 0;

%% Patients
I1=load('.../mean_ts1_PT.mat');
I2=load('.../mean_ts2_PT.mat');
I3=load('.../mean_ts3_PT.mat');
I4=load('.../mean_ts4_PT.mat');
I5=load('.../mean_ts5_PT.mat');
I6=load('.../mean_ts6_PT.mat');
I7=load('.../mean_ts7_PT.mat');
I8=load('.../mean_ts8_PT.mat');
I9=load('.../mean_ts9_PT.mat');
I10=load('.../mean_ts10_PT.mat');
I11=load('.../mean_ts11_PT.mat');
I12=load('.../mean_ts12_PT.mat');
I13=load('.../mean_ts13_PT.mat');


T=200; %time points
p=90; %nodes
% k=10; %number of subjects
%abnormal
D=zeros(T,n,K);
D(:,:,1)=normalize(I1.ts1(1:200,:));
D(:,:,2)=normalize(I2.ts1(1:200,:));
D(:,:,3)=normalize(I3.ts1(1:200,:));
D(:,:,4)=normalize(I4.ts1(1:200,:));
D(:,:,5)=normalize(I5.ts1(1:200,:));
D(:,:,6)=normalize(I6.ts1(1:200,:));
D(:,:,7)=normalize(I7.ts1(1:200,:));
D(:,:,8)=normalize(I8.ts1(1:200,:));
D(:,:,9)=normalize(I9.ts1(1:200,:));
D(:,:,10)=normalize(I10.ts1(1:200,:));
D(:,:,11)=normalize(I11.ts1(1:200,:));
D(:,:,12)=normalize(I12.ts1(1:200,:));
D(:,:,13)=normalize(I13.ts1(1:200,:));

%% Normal controls
I1=load('.../NC_mean_ts1.mat');
I2=load('.../NC_mean_ts2.mat');
I3=load('.../NC_mean_ts3.mat');
I4=load('.../NC_mean_ts4.mat');
I5=load('.../NC_mean_ts5.mat');
I6=load('.../NC_mean_ts6.mat');
I7=load('.../NC_mean_ts7.mat');
I8=load('.../NC_mean_ts8.mat');
I9=load('.../NC_mean_ts9.mat');
I10=load('.../NC_mean_ts10.mat');
I11=load('.../NC_mean_ts11.mat');
I12=load('.../NC_mean_ts12.mat');
I13=load('.../NC_mean_ts13.mat');

%% Normalized time series
D(:,:,14)=normalize(I1.ts1(1:200,:));
D(:,:,15)=normalize(I2.ts1(1:200,:));
D(:,:,16)=normalize(I3.ts1(1:200,:));
D(:,:,17)=normalize(I4.ts1(1:200,:));
D(:,:,18)=normalize(I5.ts1(1:200,:));
D(:,:,19)=normalize(I6.ts1(1:200,:));
D(:,:,20)=normalize(I7.ts1(1:200,:));
D(:,:,21)=normalize(I8.ts1(1:200,:));
D(:,:,22)=normalize(I9.ts1(1:200,:));
D(:,:,23)=normalize(I10.ts1(1:200,:));
D(:,:,24)=normalize(I11.ts1(1:200,:));
D(:,:,25)=normalize(I12.ts1(1:200,:));
D(:,:,26)=normalize(I13.ts1(1:200,:));

%% Each dataset is divided into 20 blocks
B1(1:10,:,:)=D(1:10,:,:);
B2(1:10,:,:)=D(11:20,:,:);
B3(1:10,:,:)=D(21:30,:,:);
B4(1:10,:,:)=D(31:40,:,:);
B5(1:10,:,:)=D(41:50,:,:);
B6(1:10,:,:)=D(51:60,:,:);
B7(1:10,:,:)=D(61:70,:,:);
B8(1:10,:,:)=D(71:80,:,:);
B9(1:10,:,:)=D(81:90,:,:);
B10(1:10,:,:)=D(91:100,:,:);
B11(1:10,:,:)=D(101:110,:,:);
B12(1:10,:,:)=D(111:120,:,:);
B13(1:10,:,:)=D(121:130,:,:);
B14(1:10,:,:)=D(131:140,:,:);
B15(1:10,:,:)=D(141:150,:,:);
B16(1:10,:,:)=D(151:160,:,:);
B17(1:10,:,:)=D(161:170,:,:);
B18(1:10,:,:)=D(171:180,:,:);
B19(1:10,:,:)=D(181:190,:,:);
B20(1:10,:,:)=D(191:200,:,:);

%% Stability selection
for i=10:10 %10
    for j=5:5  %5
          para_struc.alpha=0.065*i-0.03;   %0.01-0.4
          para_struc.beta=0.035*j+0.15;       %0.2-0.6

        %random permutation and subsample floor(N/2) observations
        for w=1:1
          vec(w,:)=randperm(20,10);
          ind1=sprintf('B%d',vec(w,1));
          ind2=sprintf('B%d',vec(w,2));
          ind3=sprintf('B%d',vec(w,3));
          ind4=sprintf('B%d',vec(w,4));
          ind5=sprintf('B%d',vec(w,5));
          ind6=sprintf('B%d',vec(w,6));
          ind7=sprintf('B%d',vec(w,7));
          ind8=sprintf('B%d',vec(w,8));
          ind9=sprintf('B%d',vec(w,9));
          ind10=sprintf('B%d',vec(w,10));
        %   temp1=eval(ind1);
           data=cat(1,eval(ind1),eval(ind2),eval(ind3),eval(ind4),eval(ind5),eval(ind6),eval(ind7),eval(ind8),eval(ind9),eval(ind10));  




          S=zeros(n,n,K);
          for k=1:K
                A = data(:,:,k);
                S(:,:,k) = cov(zscore(A));
                P(:,:,k) = diag(1./diag(S(:,:,k)));
          end



            para_struc.P0 = P;
            tic
            [Ps,funVals] = gmgl_rand(S,para_struc);
%             


            
            for a=1:n
                for b=1:n
                   if Ps(a,b,1)&&Ps(a,b,2)&&Ps(a,b,3)&&Ps(a,b,4)&&Ps(a,b,5)&&Ps(a,b,6)&&Ps(a,b,7)&&Ps(a,b,8)&&Ps(a,b,9)&&Ps(a,b,10)&&Ps(a,b,11)&&Ps(a,b,12)&&Ps(a,b,13)
                        temp1(a,b)=1;
                    else
                        temp1(a,b)=0;
                    end
                end
                temp1(a,a)=0;
            end
%             figure(9)
%             imagesc(temp1)
            Y1(:,:,w,i,j)=temp1(:,:);
            
            I1PT(:,:,:,w,i,j)=Ps(:,:,1:13);
            
             for a=1:n
                for b=1:n
                   if Ps(a,b,14)&&Ps(a,b,15)&&Ps(a,b,16)&&Ps(a,b,17)&&Ps(a,b,18)&&Ps(a,b,19)&&Ps(a,b,20)&&Ps(a,b,21)&&Ps(a,b,22)&&Ps(a,b,23)&&Ps(a,b,24)&&Ps(a,b,25)&&Ps(a,b,26)
                        temp2(a,b)=1;
                    else
                        temp2(a,b)=0;
                   end
                    
                   
                end
                temp2(a,a)=0;
            end
%             figure(10)
%             imagesc(temp2)
            Y2(:,:,w,i,j)=temp2(:,:);
            I2NC(:,:,:,w,i,j)=Ps(:,:,14:26);
            

            
            
            for a=1:n
                for b=1:n
                   if Ps(a,b,1)<0&&Ps(a,b,2)<0&&Ps(a,b,3)<0&&Ps(a,b,4)<0&&Ps(a,b,5)<0&&Ps(a,b,6)<0&&Ps(a,b,7)<0&&Ps(a,b,8)<0&&Ps(a,b,9)<0&&Ps(a,b,10)<0&&Ps(a,b,11)<0&&Ps(a,b,12)<0&&Ps(a,b,13)<0
                        temp1(a,b)=1;
                    else
                        temp1(a,b)=0;
                    end
                end
                temp1(a,a)=0;
            end
%             figure(11)
%             imagesc(temp1)
            Y3(:,:,w,i,j)=temp1(:,:);
            
%             I1PT(:,:,:,w,i,j)=Ps(:,:,1:13);
            
             for a=1:n
                for b=1:n
                   if Ps(a,b,14)<0&&Ps(a,b,15)<0&&Ps(a,b,16)<0&&Ps(a,b,17)<0&&Ps(a,b,18)<0&&Ps(a,b,19)<0&&Ps(a,b,20)<0&&Ps(a,b,21)<0&&Ps(a,b,22)<0&&Ps(a,b,23)<0&&Ps(a,b,24)<0&&Ps(a,b,25)<0&&Ps(a,b,26)<0
                        temp2(a,b)=1;
                    else
                        temp2(a,b)=0;
                    end
                end
                temp2(a,a)=0;
            end
%             figure(12)
%             imagesc(temp2)
            Y4(:,:,w,i,j)=temp2(:,:);
%             I2PT(:,:,:,w,i,j)=Ps(:,:,14:26);
            
            
            
            
        end

%       
      for a=1:n
          for b=1:n
              W1(a,b,i,j)=sum(Y1(a,b,:,i,j))/200;
          end
      end
      for a=1:n
          for b=1:n
              W2(a,b,i,j)=sum(Y2(a,b,:,i,j))/200;
          end
      end
      
      
      for a=1:n
          for b=1:n
              W3(a,b,i,j)=sum(Y3(a,b,:,i,j))/200;
          end
      end
      for a=1:n
          for b=1:n
              W4(a,b,i,j)=sum(Y4(a,b,:,i,j))/200;
          end
      end
      
      
    end
    
    
end


save('.../W1_FMGL_GGL_SS200_50_invivo_indiviudal_rand_200.mat','W1');
save('.../W2_FMGL_GGL_SS200_50_invivo_indiviudal_rand_200.mat','W2');
save('.../Y1_FMGL_GGL_SS200_50_invivo_indiviudal_rand_200.mat','Y1');
save('.../Y2_FMGL_GGL_SS200_50_invivo_indiviudal_rand_200.mat','Y2');
save('.../W3_FMGL_GGL_SS200_50_invivo_indiviudal_rand_200.mat','W3');
save('.../W4_FMGL_GGL_SS200_50_invivo_indiviudal_rand_200.mat','W4');
save('.../Y3_FMGL_GGL_SS200_50_invivo_indiviudal_rand_200.mat','Y3');
save('.../Y4_FMGL_GGL_SS200_50_invivo_indiviudal_rand_200.mat','Y4');


save('.../I1_FMGL_GGL_SS200_50_invivo_individual_rand_200_PT2016.mat','I1PT','-v7.3');
save('.../I1_FMGL_GGL_SS200_50_invivo_individual_rand_200_NC2016.mat','I2NC','-v7.3');
% 
