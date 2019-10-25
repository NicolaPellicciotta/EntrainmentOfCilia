function [r_bin,th_bin,fph ] = phase_fft( ph,Ntheta,NR )

[height, width, times] = size(ph);
[X,Y]=meshgrid(1:width,1:height);
max_qy= height/2;
max_qx= width/2;

if NR> 2/3*min(height/2,width/2);
disp('requested resolution is too high, it will be cutted');
NR= floor(2/3*min(height/2,width/2))
end
%% alternative map from Nicola
R = (round(sqrt((X-max_qx).*(X-max_qx)+(Y-max_qy).*(Y-max_qy)))+1);

TH= (atan2d( -(Y-max_qy),(X-max_qx)));    
TH(TH<0)=TH(TH<0)+180;
Rl=(min(R(:)));Ru=(max(R(:)))*(2/3);
dr= ((Ru-Rl)./NR); R_bin=floor(Rl:dr:(NR*dr));
dth=floor(180/Ntheta);TH_bin=0:dth:180;
th_bin= (TH_bin(1:end-1)+TH_bin(2:end))/2;
r_bin= (R_bin(1:end-1)+R_bin(2:end))/2;


%% actual fft of the phase
%(difference of actual images, then |FFT|^2, average and azimuthal average)

    for tt=1:times
         pht=ph(:,:,tt);
         pht(isnan(pht))=0;
         ham2=window2(height,width,@hamming);
 
         fp= medfilt2(abs(fft2(pht.*ham2)),[3,3]);
         fp=fftshift(fp);
         
         for ii=1:(numel(R_bin)-1)
             for jj=1:(numel(TH_bin)-1)
                    slot=((R>=R_bin(ii) & R < (R_bin(ii+1)) ) & (TH>=TH_bin(jj) & TH < (TH_bin(jj+1)) )) ;
                    fph(ii,jj,tt)=nanmean(fp(slot));
             end
         end
    end
end