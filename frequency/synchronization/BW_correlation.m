function [r_bins,Cf] = BW_correlation(movie,NR )
% Nframes=metarray(1);dt=metarray(2);Ndt=metarray(3)

if nargin < 1 || isempty(movie)
    movie='P0';
end

d=dir(strcat('*',movie,'*movie'));
mo=moviereader(d(1).name);
Res.mo=mo;

disp(mo);

fs=mo.read();
N_frames=size(fs,3);
FR=mo.FrameRate;

Ntimes=200;
fsbw= double(fs(:,:,1:Ntimes))- mean(fs(:,:,1:Ntimes),3);
for jj=1:Ntimes; fsbw(:,:,jj)=medfilt2(fsbw(:,:,jj),[3,3]);end
k=0;scroll_stack(fsbw)
nc=1;
while k==0;
    k = waitforbuttonpress;    
    BW_temp{nc}=roipoly();
    nc=nc+1;
end
cc=1;
for nc=1:numel(BW_temp)
if ~isempty(BW_temp{nc}); BW{cc}=BW_temp{nc}; cc=cc+1;end
end;close all;
    
BWT=BW_temp{1};
[height, width, times] = size(fs);
[X,Y]=meshgrid(1:width,1:height);


X=X(BWT);
Y=Y(BWT);

R = sqrt(((repmat(X,[1,numel(X)])) - (repmat(X,[1,numel(X)]))').^2 ...
       + ((repmat(Y,[1,numel(Y)])) - (repmat(Y,[1,numel(Y)]))').^2 );
TH= (atan2d( (repmat(Y,[1,numel(Y)])) - (repmat(Y,[1,numel(Y)]))',(repmat(X,[1,numel(X)])) - (repmat(X,[1,numel(X)]))'));    
TH(TH<0)=TH(TH<0)+180;   

Rl=(min(R(:))+1);
Ru=(max(R(:)))/2;
dr= floor((Ru-Rl)/NR);
if dr<1; dr=1; end
dth=floor(180/Ntheta);
%th_bin= (TH_bin(1:end-1)+TH_bin(2:end))/2;
%r_bin= (R_bin(1:end-1)+R_bin(2:end))/2;
R =1+tril(1+(floor((R+Rl)/dr)));
TH= tril(1+floor(TH/dth));
THR=TH+R;
r_bins= unique(R);r_bins=r_bins(1:end-1)*dr;
th_bins= unique(TH);th_bins=th_bins(1:end-1)*dth;
%ff=reshape(fsbw(repmat(logical(BWT),[1,Ntimes])),sum(BWT(:)),Ntimes);

Cf=zeros([floor(Ntimes/2),numel(r_bins)]);
for tau=1:floor(Ntimes/2)
    disp(tau)
    counts_sum=zeros([numel(r_bins)+1,1]);
    C_sum=zeros([numel(r_bins)+1,1]);
    for tt=1:Ntimes-tau
         f1=fsbw(:,:,tt);f2=fsbw(:,:,tt+tau);
         Mat1=f1(BWT);
         Mat2=f2(BWT);
         C= (Mat1).*Mat2';
         counts_sum = counts_sum + accumarray(R(:),ones(size(R,1)^2,1));
         C_sum = C_sum+ accumarray(R(:),C(:));
    end
    Cf(tau,:)=C_sum(2:end)./counts_sum(2:end);
end    
    
end
