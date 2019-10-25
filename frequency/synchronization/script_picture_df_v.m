cd '/media/np451/Seagate Backup Plus Drive/DATA/synchronisation_paper';
load('allvariables_area.mat');
%%
N_cells=numel(Cell)
%for nc=1:N_cells
%%
nc= 55;
f_guess= Cell(nc).F_rest;
keystr='Drive/';
keystr2='DATA/';

ind_directory= strfind(Cell(nc).Pos,keystr2)+numel(keystr2);
if isempty(ind_directory)
ind_directory= strfind(Cell(nc).Pos,keystr)+numel(keystr)
end
moviedirectory= Cell(nc).Pos(ind_directory:end); 
cd(strcat('/media/np451/Seagate Backup Plus Drive/DATA/synchronisation_paper/',moviedirectory));
box_size=4;
[Meta] = spatial_noise_fft_mask(f_guess,box_size,[],[],[]);
mask_input=Meta.mask_tot;
Res(1).F=Meta.F;
data=Meta.F(~isnan(Meta.F));
M(1)=nanmedian(data);
%medferr= nanmedian(abs(Meta.F(:)-nanmedian(Meta.F(:))));
STD(1)=nanstd(data);
MAD(1)=mad(data);
STD_M(1)=nanstd(data)/ sqrt(numel(data));
Period_array(1)=0;
mask=Meta.BW;
%%

for Volt=2:5;
%%% look for a synchronised cell and caluclate meta for it
ind=[];cc=1;
for rc=1:numel(Cell(nc).Rc)
if abs(Cell(nc).Rc{rc}.Period - 1000/M(1))<0.5 & Cell(nc).Rc{rc}.Volt==Volt;
    ind(cc)= rc;cc=cc+1;
end
end
movie=strcat('V',num2str(Volt),'*P',num2str(Cell(nc).Rc{ind(1)}.Period));
Period_array(Volt)= Cell(nc).Rc{ind(1)}.Period;
[Meta] = spatial_noise_fft_mask(1000/Cell(nc).Rc{ind(1)}.Period,box_size,movie,mask,mask_input);
Res(Volt).F=Meta.F;
M(Volt)=nanmedian(Meta.F(:));
%medferr= nanmedian(abs(Meta.F(:)-nanmedian(Meta.F(:))));
data=Meta.F(~isnan(Meta.F));
MAD(Volt)=mad(data) ;
STD(Volt)=nanstd(data);
STD_M(Volt)=nanstd(data)/ sqrt(numel(data));
mask=Meta.mask;
end


%end
save('nc_55_df_v.mat','Res','STD','mask','STD_M','M');


%% try to do a nice plot with the distribution
boxF=[];g=[];
for Volt=1:5
   boxF=cat(1,boxF,Res(Volt).F(:));
   g=cat(1,g, Volt*(~isnan(Res(Volt).F(:))));
end
boxplot(boxF,g)


%%%% plottting multiple histograms
plot_y= 0.05 +(0.9/6)*(1:5);
dy=diff(plot_y);
cleg=1;
figure('Renderer', 'painters', 'Position', [10 10 500 500])

 for kk=1:5;    
 %subplot(numel(Period_toplot)+1,1,kk);
 subplot('Position',[0.1 plot_y(kk) 0.8 dy(1) ]); 
histogram(Res(kk).F(~isnan(Res(kk).F)),[22:0.25:25]); 
 end

%% plot with multiple map of frequency
cc=1;
figure()
plot_y= 0.05 +(0.9/5)*(1:5);

plot_y= (1./5)*(1:5);
dy=diff(plot_y);
for kk=[1,3,4,5];    
 %subplot(numel(Period_toplot)+1,1,kk);
 %subplot(2,2,cc);cc=cc+1; 
 subplot('Position',[0.3 plot_y(cc) 0.6 dy(1) ]); cc=cc+1;
imagesc(kk,kk,Res(kk).F,[18,25]); 
text(-50,10,strcat('std= ',num2str(STD(kk)/M(kk),2)),'FontSize',10,'FontWeight','bold')
text(-50,20,strcat('Volts= ',num2str(Volts(kk))),'FontSize',10,'FontWeight','bold')
text(-50,30,strcat('f_ex= ',num2str(1000./Period_array(kk),3)),'FontSize',10,'FontWeight','bold')
;colorbar()
axis equal
axis off;
 end
%colorbar()
%% final plot (I hope)


mks=10;
lw=2;
Cal= load('/media/np451/Seagate Backup Plus Drive/DATA/calibration/calibration_matrix.mat');
load('nc_55_df_v.mat')
cd '/media/np451/Seagate Backup Plus Drive/DATA/synchronisation_paper';
Flow_mean= mean(Cal.Flow_int,2);
eps_plot=[1,3,5,7];
figure(76);
flow= cat(1,[0],Flow_mean(eps_plot));
medF=[];
for kk=[1,2,3,4,5]; medF(kk)= nanmedian(Res(kk).F(:));end

plot(flow,STD(:)./medF,'ko','MarkerSize',mks,'LineWidth',lw,'MarkerFaceColor','w');
xlabel('Vex mm[s]');
ylabel('stdf/f');
xlim([0,1.8])
x0=10;
y0=10;
width=400;
height=200
set(gcf,'position',[x0,y0,width,height])
set(gca,'FontSize',14);
%You can specify other units (inches, centimeters, normalized, points, or characters). For example:
%set(gcf,'normalized','points','position',[x0,y0,0.5,0.2])
fig=figure(76)
saveas(fig,'/u/homes/np451/Desktop/sync_paper/figures/figure3_new/df_v_nc55.pdf')





figure()
cc=1;
plot_x= (1./6)*(0:6);
dx=diff(plot_x);
for kk=[1,2,3,4,5];    
 %subplot(numel(Period_toplot)+1,1,kk);
 %subplot(2,2,cc);cc=cc+1; 
 subplot('Position',[plot_x(cc) 0.0 dx(1) 1 ]); cc=cc+1;
imagesc(kk,kk,Res(kk).F,[18,25]); colorbar('southoutside');
%text(-50,10,strcat('std= ',num2str(STD(kk)/M(kk),2)),'FontSize',10,'FontWeight','bold')
%text(-50,20,strcat('Volts= ',num2str(Volts(kk))),'FontSize',10,'FontWeight','bold')
%text(-50,30,strcat('f_ex= ',num2str(1000./Period_array(kk),3)),'FontSize',10,'FontWeight','bold')
axis equal
axis off;
 end
width=350;
height=100
set(gcf,'position',[x0,y0,width,height])


