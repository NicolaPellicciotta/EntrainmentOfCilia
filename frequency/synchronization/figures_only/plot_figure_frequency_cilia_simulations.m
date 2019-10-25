%%%% to make the figure here you ned first to load the data from
%%%% experimental data. For example you can use
%%%% figure_CBF_Ncilia_spnoise_PNAS.mat

load('/u/homes/np451/Desktop/sync_paper/evelyn/simulations/no external flow/FrequencyData_extV0_d1.200000e+00.mat')
load('/u/homes/np451/Desktop/sync_paper/evelyn/simulations/no external flow/FrequencyData_extV0_d4.000000e-01.mat')
load('/u/homes/np451/Desktop/sync_paper/evelyn/simulations/no external flow/FrequencyData_extV0_d6.000000e-01.mat')

B=[1,2,3,4,6,8,10,20,30];

%load('/media/np451/Seagate Backup Plus Drive/DATA/26.11.18/all_loaded_variable'); %

lw=0.3;
mks=14

%cd '/media/np451/Seagate Backup Plus Drive/DATA/synchronisation_paper';
%load('allvariables_area.mat');
%%% Number of cilia and Frequency with cells used with external flow
cc=1;FREQ=[];NCILIA=[];

lw=1;
mks=14
% 
% for nc=1:Ncells;
%     %plot(Cell(nc).Ncilia, Cell(nc).F_rest,'bo','LineWidth',lw,'MarkerSize',mks);hold on;
%     leg{cc}=Cell(nc).Pos(end-2:end); 
%     FREQ(cc)=Cell(nc).F_rest;
%     MED(cc)=Cell(nc).Sp.medf;
%     NCILIA(cc)= Cell(nc).Ncilia;
%     EMED(cc)=Cell(nc).Sp.medferr;
%     ESTD(cc)=Cell(nc).Sp.stdf;
%     ESTD_M(cc)=Cell(nc).Sp.stdf/sqrt(sum(~isnan(Cell(nc).Sp.F(:))));
%     leg{cc}=num2str(nc); ;cc=cc+1;
% end
% 
% % add this info to the other experiments with cells not use for flow experiments
% cd '/media/np451/Seagate Backup Plus Drive/DATA/synchronisation_paper/14.12.18/frequency/'
% load('FC_POST.mat')
% cc=58;
% for nc=1:N_cells 
%     FREQ(cc)=FC.Res(nc).f_guess;
%     MED(cc)=FC.Res(nc).medf;
%     NCILIA(cc)= FC.Res(nc).Ncilia;
%     EMED(cc)=FC.Res(nc).medferr;
%     ESTD(cc)=FC.Res(nc).stdf;
%     ESTD_M(cc)=FC.Res(nc).stdf/sqrt(sum(~isnan(FC.Res(nc).F(:))));
% %    leg{cc}=num2str(nc); ;
% cc=cc+1;
% end
% 
% 
% % start plot
% 
% Eper= (ESTD./FREQ);
% %onlyfirst= ones([1,numel(NCILIA)]);
% %onlyfirst(57:end)=0;
% in=~isinf(Eper) & ~isnan(Eper) & (NCILIA<25) & (NCILIA~=0)% & onlyfirst %./sqrt(NCILIA);
% %Eper= (100*EMED./FREQ);
% ind1= in & Eper<mean(Eper(in)) ;
% ind2= in &  Eper>mean(Eper(in));
% nbins=6;
% figure(30);
% % scatter(NCILIA(in),MED(in),20,Eper(in),'filled');
% colormap( [0 0 1; 0.5 0  0.5 ; 1 0 0]);
% c = colorbar();
% c.Label.String = '\eta_{s}';

%[histw,ehistw,vinterval]=hist_nico(NCILIA(ind1),MED(ind1),nbins);hold on;
%errorbar(vinterval,histw,ehistw,'--ko','LineWidth',lw,'MarkerSize',mks,'MarkerFaceColor','k');hold on;

% [histw,ehistw,vinterval2]=hist_nico(NCILIA(ind2),MED(ind2),nbins,vinterval);
% errorbar(vinterval(1:end-1),histw(1:end-1),ehistw(1:end-1),'--ro','LineWidth',lw,'MarkerSize',mks,'MarkerFaceColor','r');
% xlabel('Number of cilia per cell','FontSize',20);
% ylabel('CBF [Hz]','FontSize',20);
% ylim([8,32])
% xlim([0,23])
% legend({'data','$\eta_{s}<\overline{\eta_{s}}$','$\eta_{s}> \overline{\eta_{s}}$'},'Location', 'nw','Interpreter','latex')
% 
% 

%figure();
%[histw,ehistw,vinterval]=hist_nico(NCILIA(NCILIA<23),FREQ(NCILIA<23),nbins);
close all
figure(27)
errorbar(vinterval,histw,ehistw,'--kd','LineWidth',lw,'MarkerSize',mks*1.3,'MarkerFaceColor','k');hold on;
hold on;


plot(B,15*freq1d2,'--m^','LineWidth',lw,'MarkerSize',mks,'MarkerFaceColor','m');
plot(B,15*freq0d6,'--g>','LineWidth',lw,'MarkerSize',mks,'MarkerFaceColor','g');
plot(B,15*freq0d4,'--rs','LineWidth',lw,'MarkerSize',mks,'MarkerFaceColor','r'); 
%legend({'EXP $\eta_{s}<\overline{\eta_{s}}$','d=1.2','d=0.6','d=0.4'},'Location', 'nw','Interpreter','latex');
legend({'exp','d=1.2','d=0.6','d=0.4'},'Location', 'nw','Interpreter','latex');

set(gca,'FontSize',14);
%xlabel('Number of Cilia/ N','FontSize',20);
%ylabel('$CBF $ [Hz]','FontSize',20);
ylim([12,30]);xlim([0,31]); box on;
x0=0;
y0=0;
width=500;
height=400;
set(gcf,'position',[x0,y0,width,height])
set(gca,'FontSize',15);

fig=figure(27);
saveas(fig,'/u/homes/np451/Dropbox/Synchonisation of mammalian cilia with hydrodynamic forces/PNAS/figures/v3/figure_sim/CBF_simulations.pdf')
%saveas(gcf,'/u/homes/np451/Desktop/sync_paper/figures/figure4/frequency_increase_simulations.pdf')
%saveas(gcf,'frequency_increase_simulations.pdf');
