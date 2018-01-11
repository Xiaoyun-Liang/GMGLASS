% Reference: The Matlab code is based on the following paper:
%            Xiaoyun Liang, David N. Vaughan, Alan Connelly, Fernando Calamante. A novel
%            group-fused sparse partial correlation method for simultaneous
%            estimation of functional networks.in group comparison studies. Brain Topography, 12/2017; DOI:10.1007/s10548-017-0615-6.
%
% Copyright 2017 Florey Institute of Neuroscience and Mental Health, Melbourne, Australia
% Written by Xiaoyun Liang (Email: x.liang@brain.org.au)
% This software is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied  
% warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
%%
%%% Calculating Calculating individual-level networks for all subjects of NC group

%%
clear all
load('.../I1_FMGL_GGL_SS200_50_invivo_individual_rand_200_NC2016.mat')
Y=I2NC;
for lambda1=1:10
    for lambda2=1:5
       for k=1:13 
          for i=1:90
             for j=1:90
                 for m=1:100
                     if Y(i,j,k,m,lambda1,lambda2)<0
                         I(i,j,k,m,lambda1,lambda2)=1;
                     else
                         I(i,j,k,m,lambda1,lambda2)=0;
                     end
                 
                 end
             S(i,j,k,lambda1,lambda2)=sum(I(i,j,k,:,lambda1,lambda2))/100;
          end
       end
       end
    end
end

save('.../S_indiviudal_negative_NC.mat','S');


sum_S=zeros(90,90,10,5);
for i=1:13
    sum_S=sum_S+squeeze(S(:,:,i,:,:));
end
M_reshape=reshape(sum_S,90,90,50);
temp=M_reshape(:,:,1:50);
% 
% M_reshape=reshape(S(:,:,13,:,:),90,90,50);
% temp=M_reshape(:,:,1:50);
[C,mask]=max(temp,[],3);
M=S(:,:,1,:,:);
M_reshape=reshape(M,90,90,50);
for i=1:90
    for j=1:90
        M1(i,j)=M_reshape(i,j,mask(i,j));
    end
end

M=S(:,:,2,:,:);
M_reshape=reshape(M,90,90,50);
for i=1:90
    for j=1:90
        M2(i,j)=M_reshape(i,j,mask(i,j));
    end
end

M=S(:,:,3,:,:);
M_reshape=reshape(M,90,90,50);
for i=1:90
    for j=1:90
        M3(i,j)=M_reshape(i,j,mask(i,j));
    end
end

M=S(:,:,4,:,:);
M_reshape=reshape(M,90,90,50);
for i=1:90
    for j=1:90
        M4(i,j)=M_reshape(i,j,mask(i,j));
    end
end

M=S(:,:,5,:,:);
M_reshape=reshape(M,90,90,50);
for i=1:90
    for j=1:90
        M5(i,j)=M_reshape(i,j,mask(i,j));
    end
end

M=S(:,:,6,:,:);
M_reshape=reshape(M,90,90,50);
for i=1:90
    for j=1:90
        M6(i,j)=M_reshape(i,j,mask(i,j));
    end
end

M=S(:,:,7,:,:);
M_reshape=reshape(M,90,90,50);
for i=1:90
    for j=1:90
        M7(i,j)=M_reshape(i,j,mask(i,j));
    end
end

M=S(:,:,8,:,:);
M_reshape=reshape(M,90,90,50);
for i=1:90
    for j=1:90
        M8(i,j)=M_reshape(i,j,mask(i,j));
    end
end

M=S(:,:,9,:,:);
M_reshape=reshape(M,90,90,50);
for i=1:90
    for j=1:90
        M9(i,j)=M_reshape(i,j,mask(i,j));
    end
end

M=S(:,:,10,:,:);
M_reshape=reshape(M,90,90,50);
for i=1:90
    for j=1:90
        M10(i,j)=M_reshape(i,j,mask(i,j));
    end
end

M=S(:,:,11,:,:);
M_reshape=reshape(M,90,90,50);
for i=1:90
    for j=1:90
        M11(i,j)=M_reshape(i,j,mask(i,j));
    end
end

M=S(:,:,12,:,:);
M_reshape=reshape(M,90,90,50);
for i=1:90
    for j=1:90
        M12(i,j)=M_reshape(i,j,mask(i,j));
    end
end

M=S(:,:,13,:,:);
M_reshape=reshape(M,90,90,50);
for i=1:90
    for j=1:90
        M13(i,j)=M_reshape(i,j,mask(i,j));
    end
end

%% Thr (P_thr) is calculated for individual subject by using equation 9. For different data, one should recalculate these values accordingly.
SV=zeros(90,90); % stable variable
Thr=0.71;
for i=1:90
    for j=1:90
        if M1(i,j)>=Thr
            SV(i,j)=1;
        end
        if i==j
            SV(i,j)=0;
        end
    end
end
figure(1)
imagesc(SV)
save('C:/Work/projects/FMGL/individual/NC_new/SV1.mat','SV')

SV=zeros(90,90); % stable variable
Thr=0.72;
for i=1:90
    for j=1:90
        if M2(i,j)>=Thr
            SV(i,j)=1;
        end
        if i==j
            SV(i,j)=0;
        end
    end
end
figure(2)
imagesc(SV)
save('C:/Work/projects/FMGL/individual/NC_new/SV2.mat','SV')

SV=zeros(90,90); % stable variable
Thr=0.74;
for i=1:90
    for j=1:90
        if M3(i,j)>=Thr
            SV(i,j)=1;
        end
        if i==j
            SV(i,j)=0;
        end
    end
end
figure(3)
imagesc(SV)
save('C:/Work/projects/FMGL/individual/NC_new/SV3.mat','SV')

SV=zeros(90,90); % stable variable
Thr=0.72;
for i=1:90
    for j=1:90
        if M4(i,j)>=Thr
            SV(i,j)=1;
        end
        if i==j
            SV(i,j)=0;
        end
    end
end
figure(4)
imagesc(SV)
save('C:/Work/projects/FMGL/individual/NC_new/SV4.mat','SV')

SV=zeros(90,90); % stable variable
Thr=0.67;
for i=1:90
    for j=1:90
        if M5(i,j)>=Thr
            SV(i,j)=1;
        end
        if i==j
            SV(i,j)=0;
        end
    end
end
figure(5)
imagesc(SV)
save('C:/Work/projects/FMGL/individual/NC_new/SV5.mat','SV')

SV=zeros(90,90); % stable variable
Thr=0.64;
for i=1:90
    for j=1:90
        if M6(i,j)>=Thr
            SV(i,j)=1;
        end
        if i==j
            SV(i,j)=0;
        end
    end
end
figure(6)
imagesc(SV)
save('C:/Work/projects/FMGL/individual/NC_new/SV6.mat','SV')

SV=zeros(90,90); % stable variable
Thr=0.62;
for i=1:90
    for j=1:90
        if M7(i,j)>=Thr
            SV(i,j)=1;
        end
        if i==j
            SV(i,j)=0;
        end
    end
end
figure(7)
imagesc(SV)
save('C:/Work/projects/FMGL/individual/NC_new/SV7.mat','SV')

SV=zeros(90,90); % stable variable
Thr=0.62;
for i=1:90
    for j=1:90
        if M8(i,j)>=Thr
            SV(i,j)=1;
        end
        if i==j
            SV(i,j)=0;
        end
    end
end
figure(8)
imagesc(SV)
save('C:/Work/projects/FMGL/individual/NC_new/SV8.mat','SV')

SV=zeros(90,90); % stable variable
Thr=0.61;
for i=1:90
    for j=1:90
        if M9(i,j)>=Thr
            SV(i,j)=1;
        end
        if i==j
            SV(i,j)=0;
        end
    end
end
figure(9)
imagesc(SV)
save('C:/Work/projects/FMGL/individual/NC_new/SV9.mat','SV')

SV=zeros(90,90); % stable variable
Thr=0.61;
for i=1:90
    for j=1:90
        if M10(i,j)>=Thr
            SV(i,j)=1;
        end
        if i==j
            SV(i,j)=0;
        end
    end
end
figure(10)
imagesc(SV)
save('C:/Work/projects/FMGL/individual/NC_new/SV10.mat','SV')

SV=zeros(90,90); % stable variable
Thr=0.62;
for i=1:90
    for j=1:90
        if M11(i,j)>=Thr
            SV(i,j)=1;
        end
        if i==j
            SV(i,j)=0;
        end
    end
end
figure(11)
imagesc(SV)
save('C:/Work/projects/FMGL/individual/NC_new/SV11.mat','SV')

SV=zeros(90,90); % stable variable
Thr=0.65;
for i=1:90
    for j=1:90
        if M12(i,j)>=Thr
            SV(i,j)=1;
        end
        if i==j
            SV(i,j)=0;
        end
    end
end
figure(12)
imagesc(SV)
save('C:/Work/projects/FMGL/individual/NC_new/SV12.mat','SV')

SV=zeros(90,90); % stable variable
Thr=0.65;
for i=1:90
    for j=1:90
        if M13(i,j)>=Thr
            SV(i,j)=1;
        end
        if i==j
            SV(i,j)=0;
        end
    end
end
figure(13)
imagesc(SV)
save('C:/Work/projects/FMGL/individual/NC_new/SV13.mat','SV')


