cd '/media/np451/Seagate Backup Plus Drive/DATA/synchronisation_paper/14.12.18/frequency/'
load('FC.mat')
%%
N_cells=numel(FC.Res)
for nc=9:N_cells
f_guess= FC.Res(nc).f_guess;
movie=FC.Res(nc).filename(1:(end-6));
disp(FC.Res(nc).Ncilia);
box_size=4;
[Meta] = spatial_noise_fft_mask(f_guess,box_size,movie,[])
FC.Res(nc).Sp.F=Meta.F;
FC.Res(nc).Sp.P=Meta.P;
FC.Res(nc).Sp.P_meanf=Meta.P_meanf;
FC.Res(nc).Sp.box_size=box_size;
FC.Res(nc).Sp.area= sum(~isnan(Meta.F(:)))*box_size^2;
FC.Res(nc).Sp.medf=nanmedian(Meta.F(:));
FC.Res(nc).Sp.medferr= nanmedian(abs(Meta.F(:)-nanmedian(Meta.F(:))));
FC.Res(nc).Sp.stdf=nanstd(Meta.F(:));
FC.Res(nc).Sp.mask=Meta.mask;

end

%% add this info to the other data
cd '/media/np451/Seagate Backup Plus Drive/DATA/synchronisation_paper/14.12.18/frequency/'
load('FC_updated.mat')
cc=58;
for nc=1:N_cells 
    FREQ(cc)=FC.Res(nc).f_guess;
    MED(cc)=FC.Res(nc).Sp.medf;
    NCILIA(cc)= FC.Res(nc).Ncilia;
    EMED(cc)=FC.Res(nc).Sp.medferr;
    ESTD(cc)=FC.Res(nc).Sp.stdf;
    
%    leg{cc}=num2str(nc); ;
cc=cc+1;
end
%%
cd '/media/np451/Seagate Backup Plus Drive/DATA/synchronisation_paper/14.12.18/frequency/'
load('FC_POST.mat')
cc=58;
for nc=1:N_cells 
    FREQ(cc)=FC.Res(nc).f_guess;
    MED(cc)=FC.Res(nc).medf;
    NCILIA(cc)= FC.Res(nc).Ncilia;
    EMED(cc)=FC.Res(nc).medferr;
    ESTD(cc)=FC.Res(nc).stdf;
    ESTD_M(cc)=FC.Res(nc).stdf/sqrt(sum(~isnan(FC.Res(nc).F(:))));
%    leg{cc}=num2str(nc); ;
cc=cc+1;
end