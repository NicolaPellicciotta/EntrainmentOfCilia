%%%% frequency, area and number of cilia script.
%data_dir='/media/np451/Seagate Backup Plus Drive/DATA/1.11.18/'  %%% experiment in DMEM P/S

cd('/media/np451/Seagate Backup Plus Drive/DATA/synchronisation_paper/26.11.18/P1/')

%cd(data_dir);
d=dir('*P0*movie');
for dd=1:1%numel(d)
    mo=moviereader(d(dd).name);
    fs=mo.read();
    fsbw= double(fs(:,:,2:300))- movmean(fs(:,:,2:300),100,3);
    for kk=1:size(fsbw,3)
    fsbw(:,:,kk)= wiener2(fsbw(:,:,kk),[5,5]);
    end
    %%%%%%% count the number of cilia
    scroll_stack(fsbw)
    prompt= 'number of cilia';
    Ncilia = (input(prompt));
    FC.Res(dd).Ncilia=Ncilia;
    close all
    %%%%% define a region for the spatial fft calculation
    s= std(double(fs(:,:,1:300)),[],3)
    imagesc(s);
    rect_out= getrect();
    close all
    mask= zeros(size(fs(:,:,1)));
    [rows,cols]= rect2sub(rect_out);
    mask(rows(1):rows(end),cols(1):cols(end),:)=1;
    FC.Res(dd).rect_out= rect_out;
    FC.Res(dd).mask= mask;
   %%%%%%% calculate spectra to find fguess
    [fq,m_pxx]=average_fft(fs,logical(mask),mo.FrameRate);
    plot(fq,m_pxx);
    prompt='frequency guess'
    f_guess = (input(prompt));
    
    FC.Res(dd).f_guess=f_guess;
    FC.Res(dd).fq=fq;
    FC.Res(dd).m_pxx = m_pxx;
    FC.Res(dd).filename = d(dd).name;
    
end
%% find fguess and frequency peaks

for dd=1:numel(FC.Res)
    box_size=6;
    movie= FC.Res(dd).filename;
    movie=movie(end-8-6:end-6);
    f_guess=FC.Res(dd).f_guess;
    rect_out=FC.Res(dd).rect_out;
    [Meta] = spatial_noise_fft(f_guess,box_size,movie,rect_out);
    FC.Res(dd).F=Meta.F;
    FC.Res(dd).s_roi=Meta.s_roi;
    FC.Res(dd).area= sum(~isnan(Meta.F(:)));
    FC.Res(dd).medf=nanmedian(Meta.F(:));
    FC.Res(dd).medferr= nanmedian(abs(Meta.F(:)-nanmedian(Meta.F(:))));
    FC.Res(dd).stdf=nanstd(Meta.F(:));
    
end


%% find spatial noise






