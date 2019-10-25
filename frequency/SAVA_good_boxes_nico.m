function [BW]= SAVA_good_boxes_nico(cilia,bsz,N_max_px);
cc=1;
smooth_s= mat2gray(medfilt2(cilia.lstd_sfd,[5,5]));

[level,EM] = graythresh(smooth_s);
bsz=4;
cc=1;
pp=1;
while cc==1 
    level_pp=level*pp;
BW = imresize(imbinarize(smooth_s,level_pp),1/bsz);
    if sum(BW(:))<N_max_px; cc=0;
    else cc=1; pp=pp+0.1;
    end
end


end