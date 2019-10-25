close all;
clear all

%%% looking for a movie without external flow and select a mask from it %%%%
d=dir('*P0*'); %%%% P0 stands for period 0

cd(d(end).name); %%%%% move in the directotry of the file

%%%% load all the images in a matrix fs
A= dir('*.tiff');
for kk=1:numel(A) 
    fs(:,:,kk) = imread(A(kk).name);
end
cd('..');

N_frames= size(fs,3);
s=std(double(fs(:,:,2:300)),[],3); %%% calculate standard deviation map
s=mat2gray(s);

%%%%% remove background with a moving mean
fsbw= double(fs(:,:,2:300))- mean(fs(:,:,2:300),3); 

%%%%% script to select the cell, you need to close the polygon and then
%%%%% press k and then close the window
k=0;
scroll_stack(fsbw)
%nc=1;
nc=11;
while k==0;
k = waitforbuttonpress;    
BW_temp{nc}=roipoly();
nc=nc+1;
end


%%%%%% create a mask that will be BW and do an image of the mask over the
%%%%%% standard deviation
cc=1;
for nc=1:numel(BW_temp)
if ~isempty(BW_temp{nc}); BW{cc}=BW_temp{nc}; cc=cc+1;end
end
%close all;
mkdir('plots')
%figure();
Reds= cat(3,s(:,:,1),0*s(:,:,1),0*s(:,:,1));  
    imshow(imadjust(s));
    BWT=zeros(size(s));
    for nc=1:numel(BW);
        BWT= BWT+BW{nc};
        [maskx,masky]=find(BW{nc});
        text(masky(1)-30,maskx(1)-30,num2str(nc),'Color','red','FontSize',24)
    end
    hold on
    h = imshow(Reds); % Save the handle; we'll need it later
    hold off        
    set(h, 'AlphaData', BWT); 
    saveas(gca,'cell_analysed.jpg')


%%
Ncell=numel(BW); %%%% in our case should be 1
cd  ()
%%%%% now load all the other videos at different applied flow the we
%%%%% calculate the periodogram over all the selected pixel in the mask
%%%%% all the information are stored in the variable R

d=dir('*_V*movie');
for jj=1:numel(d);
    disp(jj)
    
%%%%% load data and the roi    
    filename=d(jj).name;
    cd(filename)
    A= dir('*.tiff');
    for kk=1:numel(A) 
        fs(:,:,kk) = imread(A(kk).name);
    end
    cd('..')
    N_frames= size(fs,3);
    FR =mo.FrameRate;
    
%%% store in a variable RC all the useful information such as the period of the
%%% external flow applied, the voltage (that correspond to a claibrated flow velocity)
    for nc=1:Ncell
    
     R{nc,jj}.filename=filename;
    
    %%%%% save voltage %%%%
    Vind=strfind(R{nc,jj}.filename,'V');
    Vend=strfind(R{nc,jj}.filename,'_O');
    R{nc,jj}.Volt= str2num(R{nc,jj}.filename(Vind+1:Vend-1));  
    %%%%% save offset %%%%%
    Offind=strfind(filename,'O');
    R{nc,jj}.Offset= str2num(filename(Offind+1));  
    %%%%%% save Period %%%%%%%
    Pind=strfind(filename,'P');
    Endind=strfind(R{nc,jj}.filename,'.27Nov');
    R{nc,jj}.Period= str2num(filename(Pind+1:Endind-1));
    %%%%%% save Framerate %%%%%%    
    R{nc,jj}.FR = FR;
    
    %%%%% place in the variale roi all the intensity pixel over time that
    %%%%% where selected by the mask
    temp_BW=BW{nc};
    fs_roi= fs(repmat(temp_BW,[1,1,N_frames]));
    fs_roi=reshape(fs_roi,[sum(temp_BW(:)),N_frames]);
    roi=double(fs_roi)- movmean(fs_roi,30,2);

    %%%% Calculate the periodogram %%%%%%%%%%
    %%%% periodogram needs time on the row and pixel on the coloun; 
    window = hann(floor(N_frames));
    window= repmat(window,[1,size(roi,1)])';
    n= floor(N_frames/2);
    if mod(n,2)==0; n= n-1;end
    pxx= abs(fft(double(roi).*window,n,2)).^2;
    R{nc,jj}.m_pxx= mean(pxx(:,1:floor(n/2)),1);
    fq= (0:(FR./n):(FR./2-FR./n));    
    R{nc,jj}.fq= fq;     
    end
end
    
    
    save('fft_results_fft.mat','R','BW','s');
    