function [] = movie2avi_nicola_timestamp(minFrame,maxFrame,times_slower,skip_frames,flag_crop)
%movie2mp4 Interactive program to save a .movie file to an .mp4 video

[filename, pathname] = uigetfile('*.movie');

mo = moviereader(fullfile(pathname,filename));
ind=minFrame:skip_frames:maxFrame;
fs=zeros([mo.height,mo.width,numel(ind)]);

kk=1;
for i= ind;
    fs(:,:,kk)=mo.read(i);kk=kk+1;
end
fs=fs-mean(fs,3);
fs=(fs-min(fs(:)))/(max(fs(:))-min(fs(:)));
%fs=imadjustn(fs,stretchlim(fs(:)));

if flag_crop==1
    k=0;nc=1;
    scroll_stack(fs)
    while k==0;
    k = waitforbuttonpress;    
    Rect{nc}=getrect();
    end
    rect=Rect{1};
    [rows, cols]= rect2sub(rect)
else
    rows=1:mo.height;
    cols=1:mo.width;
end

fs=fs(rows,cols,:);  %%%% select only the cropped 
fs=(fs-min(fs(:)))/(max(fs(:))-min(fs(:)));  %%%% normalise intensity to 1


outputName = fullfile(pathname,filename);
outputName(end-5:end) = [];
outputName = [outputName,'.avi'];

v = VideoWriter(outputName);


v.FrameRate= mo.FrameRate/(skip_frames*times_slower);
v.Quality = 95;
open(v);
ind_frames= minFrame:times_slower:maxFrame;
map = colormap(gray(256));

for i=1:size(fs,3);
    time_stamp= ind(i)/mo.FrameRate
   fig = figure;
   t = text(.05,.1,strcat(num2str(time_stamp,'%1.3f'),' s'),'FontSize',20,'FontWeight','bold');
   F = getframe(gca,[10 10 200 200]);
   close(fig)  
   c = F.cdata(:,:,1);
   [ii,j] = find(c==0);
   frame= medfilt2(fs(:,:,i),[2,2]);
   ind_stamp = sub2ind(size(frame),ii,j);
   frame(ind_stamp) = 1;
   writeVideo(v,frame);
end

close(v)
close all
end