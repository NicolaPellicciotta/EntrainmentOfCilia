%%%%%%% script to find out if frequency increase is similar in cell treated
%%%%%%% with cyto-D
%%%%%% I need to find number of cilia, frequency map (and std) and area

folder= '/media/np451/Seagate Backup Plus Drive1/DATA/Cyto_2.4.19/40hr/'
subfolders={'control','Cyto1','Cyto2'};

datadir=strcat(folder,subfolders{3})
cd(datadir);
d=dir('60X_P1*.movie');
for nc=1:numel(d)
movie=d(nc).name(1:(end-6));
box_size=4;
[Meta] = spatial_noise_fft_mask_ncilia_fguess(box_size,movie,[])
Cyto(nc).F=Meta.F;
Cyto(nc).P=Meta.P;
Cyto(nc).P_meanf=Meta.P_meanf;
Cyto(nc).box_size=box_size;
Cyto(nc).area= sum(~isnan(Meta.F(:)))*box_size^2;
Cyto(nc).medf=nanmedian(Meta.F(:));
Cyto(nc).medferr= nanmedian(abs(Meta.F(:)-nanmedian(Meta.F(:))));
Cyto(nc).stdf=nanstd(Meta.F(:));
Cyto(nc).mask=Meta.mask;
Cyto(nc).s=Meta.s;
Cyto(nc).s_bin=Meta.s_bin;
Cyto(nc).filename= d(nc).name;
Cyto(nc).Ncilia=Meta.Ncilia;

end

save('Cyto2.mat','Cyto');

%% for large field of view only
folder= '/media/np451/Seagate Backup Plus Drive1/DATA/Cyto_2.4.19/40hr/'
subfolders={'control','Cyto1','Cyto2'};

datadir=strcat(folder,subfolders{2})
cd(datadir);
d=dir('*Large*.movie');
cc=1;
for nc=1:numel(d)
movie=d(nc).name(1:(end-6));
for kk=1:5;
box_size=4;
[Meta] = spatial_noise_fft_mask_ncilia_fguess(box_size,movie,[])
Cyto(cc).F=Meta.F;
Cyto(cc).P=Meta.P;
Cyto(cc).P_meanf=Meta.P_meanf;
Cyto(cc).box_size=box_size;
Cyto(cc).area= sum(~isnan(Meta.F(:)))*box_size^2;
Cyto(cc).medf=nanmedian(Meta.F(:));
Cyto(cc).medferr= nanmedian(abs(Meta.F(:)-nanmedian(Meta.F(:))));
Cyto(cc).stdf=nanstd(Meta.F(:));
Cyto(cc).mask=Meta.mask;
Cyto(cc).s=Meta.s;
Cyto(cc).s_bin=Meta.s_bin;
Cyto(cc).filename= d(nc).name;
Cyto(cc).Ncilia=Meta.Ncilia;
cc=cc+1;
end
end

save('CytoLarge.mat','Cyto');


%% plot all data for frequecny graph

%% Number of cilia and Frequency
close all
set(0,'defaulttextInterpreter','latex') ;

cc=1;FREQ=[];NCILIA=[];
folder= '/media/np451/Seagate Backup Plus Drive1/DATA/Cyto_2.4.19/40hr/'
subfolders={'control','Cyto1','Cyto2'};
lw=1;
mks=14
   

for subf=[2,3]
    datadir=strcat(folder,subfolders{subf})
    cd(datadir);load('Cyto2.mat');
for nc=1:numel(Cyto);
    %plot(Cell(nc).Ncilia, Cell(nc).F_rest,'bo','LineWidth',lw,'MarkerSize',mks);hold on;
    data= Cyto(nc).F(~isnan(Cyto(nc).F));
    FREQ(cc)=mean(data(:));
    MED(cc)=Cyto(nc).medf;
    NCILIA(cc)= Cyto(nc).Ncilia;
    EMED(cc)=Cyto(nc).medferr;
    ESTD(cc)=Cyto(nc).stdf;
    ESTD_M(cc)=Cyto(nc).stdf/sqrt(numel(data));
    cc=cc+1;
end
if ~isempty(dir('CytoLarge.mat'))
    load('CytoLarge.mat');end
for nc=1:numel(Cyto);
    %plot(Cell(nc).Ncilia, Cell(nc).F_rest,'bo','LineWidth',lw,'MarkerSize',mks);hold on;
    data= Cyto(nc).F(~isnan(Cyto(nc).F));
    FREQ(cc)=mean(data(:));
    MED(cc)=Cyto(nc).medf;
    NCILIA(cc)= Cyto(nc).Ncilia;
    EMED(cc)=Cyto(nc).medferr;
    ESTD(cc)=Cyto(nc).stdf;
    ESTD_M(cc)=Cyto(nc).stdf/sqrt(numel(data));
    cc=cc+1;
end


end
NCILIA(NCILIA==14 & MED<10) = 0;

Eper= (ESTD./FREQ);
%onlyfirst= ones([1,numel(NCILIA)]);
%onlyfirst(57:end)=0;
in=~isinf(Eper) & ~isnan(Eper) & (NCILIA<25) & (NCILIA~=0) & FREQ>3;% & onlyfirst %./sqrt(NCILIA);
%Eper= (100*EMED./FREQ);
ind1= in & Eper<=nanmean(Eper(in))+0.01 ;
ind2= in &  Eper>nanmean(Eper(in))+0.01;


nbins=6;
figure(99);
scatter(NCILIA(in),MED(in),20,Eper(in),'filled');
colormap([0 0 1;0 1 0; 1 0 0])
c = colorbar();
c.Label.String = '\eta_{s}';


nbins=5;
%[histw,ehistw,bins]=hist_nico(NCILIA(ind1),MED(ind1),nbins,[1:5:26]);hold on;
[histw,ehistw,bins]=hist_nico_std(NCILIA(ind1),MED(ind1),nbins,[]);hold on;

errorbar(bins,histw,ehistw,'--ko','LineWidth',lw,'MarkerSize',mks,'MarkerFaceColor','b');
xlabel('Number of cilia per cell','FontSize',20);
ylabel('CBF [Hz]','FontSize',20);
ylim([5,25]);


%[histw,ehistw,bins]=hist_nico(NCILIA(ind2),MED(ind2),nbins,[1:5:26]);
[histw,ehistw,bins]=hist_nico_std(NCILIA(ind2),MED(ind2),nbins,[]);
errorbar(bins(1:end-1),histw(1:end-1),ehistw(1:end-1),'--ro','LineWidth',lw,'MarkerSize',mks,'MarkerFaceColor','r');
xlabel('Number of cilia per cell','FontSize',20);
ylabel('CBF [Hz]','FontSize',20);
ylim([5,25]);

title('cell treated with cytochalasin D')


ylim([8,32])
xlim([0,20])
legend({'data','$\eta_{s}<\overline{\eta_{s}}$','$\eta_{s}> \overline{\eta_{s}}$'},'Location', 'nw','Interpreter','latex')

x0=0;
y0=0;
width=500;
height=400;
set(gcf,'position',[x0,y0,width,height])
set(gca,'FontSize',15);
%set(0,'defaulttextInterpreter','latex') ;

title(strcat('mean spatial noise used $\overline{\eta_{s}}=',num2str(mean(Eper(in)),3),'$'...,'\pm',num2str(eint_d_all,1)
    ),'Interpreter','latex','FontSize',9)

fig=figure(99);
saveas(fig,'/u/homes/np451/Desktop/sync_paper/figures/cytoCBF/cyto_cbf2.pdf')



