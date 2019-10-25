function [Centrin] = Centrin_fun()
%%%%% measure the number of cilia by measuring the area of centriols in the
%%%%% cells

d=dir('*centrin*movie');
Centrin.movie='Centrin';
Objective=60;
px2mu= 0.134*40/60;


if isempty(d)==0;

mo=moviereader(d(1).name);
fs=mo.read();
fm=mean(fs,3);
figure();imagesc(fm);BW=roipoly();
close all

gfm=imgaussfilt(fm,3);  %%% removing noise camera
gfm=mat2gray(gfm);
%gfm=gfm/max(gfm(:));  %%% normalising to 1;
[level,EM]=graythresh(gfm(BW));


thresh=multithresh(fm,2);
seg_BW = imquantize(fm,thresh);
th_mask = seg_BW>2

%%% using Otsu's method to find a threeshold
%th_mask= zeros(size(gfm));
%th_mask(gfm>=level)=1;

mask= BW & th_mask;

[G,id]=findgroups(mask(:));
CC = bwconncomp(mask)
S = regionprops(CC, 'Area');
%Remove small objects:
L = labelmatrix(CC);
Area=[];for fff=1:numel(S); Area(fff)=S(fff).Area;end
max_Area= max(Area);    
mask = ismember(L, find([S.Area] >= max_Area-1))



area= sum(mask(:));
area_mu=area*(px2mu)^2;

mean_int= mean(fm(mask(:)))/(2^16-1); 

Centrin.size=area; Centrin.size_mu=area_mu;
Centrin.fm=fm;Centrin.gfm=gfm;Centrin.level=level;Centrin.px2mu=px2mu;
Centrin.BW=BW; Centrin.mask=mask; Centrin.th_mask=th_mask;
Centrin.mean_int=mean_int; 


%%% make plots
subplot(1,2,1);
imagesc(fm);
title('centrin average')
subplot(1,2,2);
imagesc(mask);
title('centrin postpro');



else Centrin.movie=[];

end
end