function [Meta] = spatial_noise_fft_meta(f_guess,box_size,movie,rect_out,metarray )
% Nframes=metarray(1);dt=metarray(2);Ndt=metarray(3)

if nargin < 3 || isempty(movie)
    movie='P0';
end

d=dir(strcat('*',movie,'*movie'));
mo=moviereader(d(1).name);
Meta.mo=mo;

if nargin < 5 || isempty(metarray)
    Nframes=mo.NumberOfFrames;
    dt=1;
    Ndt=1;
end

disp(mo)
Nframes=metarray(1);dt=metarray(2);Ndt=metarray(3);

fq_max=f_guess+7;fq_min=f_guess-7;


fs=mo.read();
N_frames=size(fs,3);
FR=mo.FrameRate;
s= std(double(fs(:,:,1:300)),[],3);
s=medfilt2(s,[3,3]);
[level,EM] = graythresh(s);
s_bin = imbinarize(s,level);   %%%%% standard deviation of the movie binarised between 1 and 0


[level,EM]= graythresh(s);
BW=imbinarize(s,level);
if isempty(rect_out);
    imagesc(s);
    rect_out= getrect();    
end

[cols,rows]= rect2sub(rect_out);
fs_roi= fs(cols(1):cols(end),rows(1):rows(end),:);
fsm_roi=double(fs_roi)-movmean(fs_roi,100,3);
row_dim= size(fsm_roi,2);
col_dim= size(fsm_roi,1);

s_roi=s(cols(1):cols(end),rows(1):rows(end));
s_bin=s_bin(cols(1):cols(end),rows(1):rows(end));
mask_tot= zeros(size(s_roi));
[level,EM] = graythresh(s_roi);


%box_size=8;   %%%% setting the box size 
Nbox=floor((row_dim*col_dim)/box_size^2);  %%%% number of boxes based on the total area
nboxes_row=floor(row_dim/box_size);
nboxes_col=floor(col_dim/box_size);
[X,Y]=meshgrid(1:box_size:(nboxes_row)*box_size,1: box_size: (nboxes_col)*box_size);

F= zeros([size(X),Ndt]);
P= zeros([size(X),Ndt]);
P_meanf= zeros([size(X),Ndt]);


times=1:dt:Ndt*dt;
for tt=1:Ndt;

for xx=1:(size(X,2)-1);
    for yy=1:(size(X,1)-1);    

        std_value= s_bin(Y(yy,xx):(Y(yy+1,xx+1)-1),X(yy,xx):(X(yy+1,xx+1)-1));
        if mean(std_value(:))> 0.7; good=1;
            mask_tot(Y(yy,xx):(Y(yy+1,xx+1)-1),X(yy,xx):(X(yy+1,xx+1)-1))=1;
        else good=0;
        end
        Meta.Box(yy,xx).good=good;
        
        roi=fsm_roi(Y(yy,xx):(Y(yy+1,xx+1)-1),X(yy,xx):(X(yy+1,xx+1)-1), times(tt):times(tt)+Nframes);
        
        
        roi2=reshape(roi,[size(roi,1)*size(roi,2),Nframes+1]);
        roi=roi2;
        Meta.Box(yy,xx).vec_pos=[Y(yy,xx):(Y(yy+1,xx+1)-1),X(yy,xx):(X(yy+1,xx+1)-1)];
        


        %%%%%% calculate fft on the roi
        window = hann(floor(Nframes+1));
        window= repmat(window,[1,size(roi,1)])';
        %n= floor(N_frames/2);
        n=Nframes+1;
        if mod(n,2)==0; n= n-1;end
        
        f_spectrum= fft(double(roi).*window,n,2);
        phase=  angle(f_spectrum);
        pxx= abs(f_spectrum).^2;
        m_pxx= mean(pxx(:,1:floor(n/2)),1);
        m_phase= mean(exp(i*phase),1);
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
            F(yy,xx,tt)=nan;
            P(yy,xx,tt)=nan;
            P_meanf(yy,xx,tt)=nan;
        else
            F(yy,xx,tt)=locs(end);
            P(yy,xx,tt)=m_phase(fq==locs(end));
            [~,ind_frest]=min(abs(fq-f_guess));
            P_meanf(yy,xx,tt)=phase(ind_frest);
        end
        
        Meta.Box(yy,xx).fp1 = F(yy,xx);
        Meta.Box(yy,xx).phase = P(yy,xx);
        Meta.Box(yy,xx).phase_meanf = P_meanf(yy,xx);
end
end

end


Meta.F=F;Meta.P=P;Meta.P_meanf=P_meanf; Meta.s=s;Meta.s_bin=s_bin;Meta.s_roi=s_roi;Meta.N_frames=N_frames;Meta.dt=dt;Meta.Ndt=Ndt;Meta.Nframes=Nframes;
Meta.box_size=box_size;
Meta.rect_out=rect_out;
Meta.median_freq= nanmedian(F(:));
Meta.std_freq= nanstd(F(:));
Meta.X=X;
Meta.Y=Y;
Meta.mask_tot=mask_tot;
end