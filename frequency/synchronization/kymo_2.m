function [r_bins,Cf] = kymo_2(movie,NR )
% NOT finished!!!!
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
    [cx,cy,c,xi,yi] = improfile()
%    nc=nc+1;
end
kymo=zeros([numel(c),Ntimes]);


for jj=1:Ntimes
    kymo(:,jj)=improfile(fsbw(:,:,jj),xi,yi);    
end


    
end
