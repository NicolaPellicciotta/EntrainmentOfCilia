function [Meta] = spatial_noise_fft(f_guess,box_size,movie,rect_out)
% 
if nargin < 3 || isempty(movie)
    movie='P0';
    d=dir(strcat('*',movie,'*movie'));
else
    d=dir(movie);
  end

mo=moviereader(d(1).name);
Meta.mo=mo;

disp(mo)

fq_max=f_guess+7;fq_min=f_guess-7;


fs=mo.read();
N_frames=size(fs,3);
FR=mo.FrameRate;
s= std(double(fs(:,:,1:300)),[],3);
s=medfilt2(s,[4,4]);
s=mat2gray(s);
%[level,EM] = graythresh(s);
thresh=multithresh(s,2);
seg_BW = imquantize(s,thresh);
%BW=imbinarize(s,level);
if isempty(rect_out);
    imagesc(s);colorbar();axis equal
    rect_out= getrect();    
end
BW= seg_BW; BW(BW==1)=0; BW(BW>1)=1; BW=logical(BW);

Meta.BW=BW;
Meta.seg_BW= seg_BW;

[cols,rows]= rect2sub(rect_out);  %%% get the doordinate of the rect
fs_roi= fs(cols(1):cols(end),rows(1):rows(end),:); %%% crop the frames_stack
fsm_roi=double(fs_roi)-movmean(fs_roi,100,3); %%%% average the frame_stack
row_dim= size(fsm_roi,2);
col_dim= size(fsm_roi,1);

s_roi=s(cols(1):cols(end),rows(1):rows(end)); %%% crop the std with rec
s_bin=BW(cols(1):cols(end),rows(1):rows(end)); %%% crop the binarised std
mask_tot= zeros(size(s_roi));  


%box_size=8;   %%%% setting the box size 
Nbox=floor((row_dim*col_dim)/box_size^2);  %%%% number of boxes based on the total area
nboxes_row=floor(row_dim/box_size);
nboxes_col=floor(col_dim/box_size);
%%% coordinate of all the boxes in the rectangle
[X,Y]=meshgrid(1:box_size:(nboxes_row)*box_size,1: box_size: (nboxes_col)*box_size);



F= zeros([size(X)]);
P= zeros([size(X)]);
P_meanf= zeros([size(X)]);

for xx=1:(size(X,2)-1);
    for yy=1:(size(X,1)-1);    

        std_value= s_bin(Y(yy,xx):(Y(yy+1,xx+1)-1),X(yy,xx):(X(yy+1,xx+1)-1));
 
        %%%% if the cropped standard deviation has enough one, calculate
        %%%% the fft
        if mean(std_value(:))> 0.8; good=1;
            mask_tot(Y(yy,xx):(Y(yy+1,xx+1)-1),X(yy,xx):(X(yy+1,xx+1)-1))=1;
        else good=0;
        end
        Meta.Box(yy,xx).good=good;
        
        roi=fsm_roi(Y(yy,xx):(Y(yy+1,xx+1)-1),X(yy,xx):(X(yy+1,xx+1)-1),:);
        roi2=reshape(roi,[size(roi,1)*size(roi,2),N_frames]);
        roi=roi2;
        Meta.Box(yy,xx).vec_pos=[Y(yy,xx):(Y(yy+1,xx+1)-1),X(yy,xx):(X(yy+1,xx+1)-1)];
        


        %%%%%% calculate fft on the roi
        window = hann(floor(N_frames));
        window= repmat(window,[1,size(roi,1)])';
        %n= floor(N_frames/2);
        n=N_frames;
        if mod(n,2)==0; n= n-1;end
        
        f_spectrum= fft(double(roi).*window,n,2);
        phase=angle(f_spectrum);
        pxx= abs(f_spectrum).^2;
        m_pxx= mean(pxx(:,1:floor(n/2)),1);
        fq= (0:(FR./n):(FR./2-FR./n));    
        Meta.fq=fq;
        Meta.Box(yy,xx).m_pxx=m_pxx; 
        Meta.Box(yy,xx).f_spectrum=f_spectrum;
        
        %%%%% calculate the peak in frequency
        
        f_range= fq> (fq_min) & fq<(fq_max);
        baseline= min((m_pxx(f_range)));
        [pks,locs,w,p] = findpeaks((m_pxx(f_range)),fq(f_range));%%%% find the peaks frequency in the selected freq range
        [~,ind_sort]= sort(pks);                                           %%%% sort peaks and get an index 
        pks=pks(ind_sort); locs=locs(ind_sort);w=w(ind_sort);p=p(ind_sort);  %%% order all the variables with the same index
        Meta.Box(yy,xx).pks=pks; Meta.Box(yy,xx).locs=locs; Meta.Box(yy,xx).w=w; Meta.Box(yy,xx).p=p;     %%%%% load results in Rc
        Meta.Box(yy,xx).baseline=baseline;
        [~,where] = max(pks);       
        if numel(locs)==0 | good==0; 
            Meta.Box(yy,xx).fp1 =nan;   %%% fp1c is the frequency with higher peak. 
            F(yy,xx)=nan;
            P(yy,xx)=nan;
            P_meanf(yy,xx)=nan;
        else
            F(yy,xx)=locs(end);
            P(yy,xx)=phase(fq==locs(end));
            [~,ind_frest]=min(abs(fq-f_guess));
            P_meanf(yy,xx)=phase(ind_frest);
        end
        
        Meta.Box(yy,xx).fp1 = F(yy,xx);
        Meta.Box(yy,xx).phase = P(yy,xx);
        Meta.Box(yy,xx).phase_meanf = P_meanf(yy,xx);
end
end

Meta.F=F;Meta.P=P;Meta.P_meanf=P_meanf; Meta.s=s;Meta.s_bin=s_bin;Meta.s_roi=s_roi;Meta.N_frames=N_frames;
Meta.box_size=box_size;
Meta.rect_out=rect_out;
Meta.median_freq= nanmedian(F(:));
Meta.std_freq= nanstd(F(:));
Meta.X=X;
Meta.Y=Y;
Meta.mask_tot=mask_tot;
end