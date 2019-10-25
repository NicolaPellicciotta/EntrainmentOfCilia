clear all
%%%% set data for plot simulation Evelyn
B=[1,2,4,6,8,10,20,30];
%%%% load simulation with chains 0.6 and 1.2 um distance
load('/u/homes/np451/Desktop/sync_paper/evelyn/simulations/external_flow/PlateauWidths_d0d6_1d2_date181218.mat');
pw6= mean(pwidth0d6,2);
pw6E= pwidth0d6ERR;
pw12= mean(pwidth1d2,2);
pw12E = pwidth1d2ERR;

%%%%% load data for 0.4 um distance rowers
%here is the data for the high density case. The standard definitions apply for STRONG variable B = [2 3 4 6]
%For the WEAK data B=[8 10 20]
%CASE one is (:,:.1) is for very weak definition, the number of states with more than 50% synchronised
%CASE two is (:,:,2) is for the moderate definition, the number of states with more than 80% synchronised.
load('/u/homes/np451/Desktop/sync_paper/evelyn/simulations/external_flow/PlateauWidths_d0d4_date191218.mat')

B4=[1,2,4,6,8,10,20];
pw4weak=cat(1,mean(pwidthsSTRONG,2),mean(pwidthsWEAK(:,:,1),2));
pw4mod=cat(1,mean(pwidthsSTRONG,2),mean(pwidthsWEAK(:,:,2),2));
pw4weakE=cat(1,pwidthErr,mean(pwidthsWEAKerr(:,:,1),2));
pw4modE=cat(1,pwidthErr,mean(pwidthsWEAKerr(:,:,2),2));


load('/u/homes/np451/Desktop/sync_paper/evelyn/simulations/external_flow/PlateauWidths_d0d5_date020119.mat');


pw5weak=cat(1,mean(pwidthsSTRONG0d5,2),mean(pwidthsWEAK0d5(:,:,1),2));
pw5mod=cat(1,mean(pwidthsSTRONG0d5,2),mean(pwidthsWEAK0d5(:,:,2),2));
pw5weakE=cat(1,pwidthErr0d5,mean(pwidthsWEAKerr0d5(:,:,1),2));
pw5modE=cat(1,pwidthErr0d5,mean(pwidthsWEAKerr0d5(:,:,2),2));


%%%%%load data from experiments

%load('/media/np451/Seagate Backup Plus Drive/DATA/26.11.18/all_loaded_variable'); 
load('/u/homes/np451/Desktop/sync_paper/evelyn/simulations/external_flow/EPS_CILIA_VOLT3.mat');
cd '/media/np451/Seagate Backup Plus Drive/DATA/synchronisation_paper';
load('allvariables_area.mat');


%%%%%%% load experimental data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Volts_plot=[3];
lw=0.3;
mks=12

color_scatter_array={'ko','bo','ro','mo'};
%color_scatter_array={'ko','ko','ko','ko'};
color_median_array={'k--o','r--s','b--<','m--^'};
color_face_array={'k','r','b','m'};
clear l;cleg=1;
  
Flow_mean= mean(Cal.Flow_int,2);


figure();
whichVolt=Volts_plot;
v_eps= (Volts == whichVolt);
l{cleg}= strcat(num2str(Flow_mean(v_eps),'%-5.2f'),' mm/s'); cleg=cleg+1

cc=1;
EPS=[]; NCILIA=[];
ms=11;
%lw=0.5;
lw=1;jj=1;
color_scatter=color_scatter_array{jj};
color_median=color_median_array{jj};
color_face=color_face_array{jj};

figure(jj)
for nc=1:Ncells;
    if Cell(nc).good
        if whichVolt<4 & nc>28
%        plot(Cell(nc).Ncilia,Cell(nc).eps(v_eps),color_scatter,'MarkerSize',ms,'LineWidth',lw,'MarkerFaceColor',[1,1,1]);hold on;
%        leg{cc}=num2str(nc); 
        EPS(cc)=Cell(nc).eps(v_eps);
        NCILIA(cc)=Cell(nc).Ncilia;
        cc=cc+1;
        end

        if whichVolt>=4 
%        plot(Cell(nc).Ncilia,Cell(nc).eps(v_eps),color_scatter,'MarkerSize',ms,'LineWidth',lw,'MarkerFaceColor',[1,1,1]);hold on;
%        leg{cc}=num2str(nc); 
        EPS(cc)=Cell(nc).eps(v_eps);
        NCILIA(cc)=Cell(nc).Ncilia;
        cc=cc+1;
        end

    end
end
%legend(leg)
%xlabel('$N_{c}$ number of cilia','Interpreter','latex', 'FontSize',20);
%ylabel('$\epsilon [Hz]$ ','Interpreter','latex', 'FontSize',20);
%title(strcat('average flow=  ',num2str(Flow_mean(v_eps),'%.1f'),'[mm/s]'),'Interpreter','latex', 'FontSize',15);
%saveas(gca,strcat('average flow  ',num2str(Flow_mean(v_eps)),'.jpg'));


nbins=5;
figure(5);

[histw,ehistw,vinterval]=hist_nico_std(NCILIA,EPS,nbins,[]);
histw(histw>9)=nan;  %%%% to remove from statistic the one the went out of range 
%errorbar(vinterval,histw,ehistw,color_median,'MarkerSize',mks,'LineWidth',lw,'MarkerFaceColor',color_face);hold on;
%errorbar(vinterval,histw./(Flow_mean(v_eps)),ehistw./(Flow_mean(v_eps)),k,'MarkerSize',mks,'LineWidth',lw,'MarkerFaceColor',k);hold on;

nbins=5;
[histw,ehistw,vinterval]=hist_nico(NCILIA,EPS,nbins,[1:5:30]);
histw(histw>9)=nan;  %%%% to remove from statistic the one the went out of range 
ehistw(ehistw<0.5)=0.5

errorbar(vinterval,histw,ehistw,'--kd','LineWidth',lw,'MarkerSize',mks*1.3,'MarkerFaceColor','k');hold on;
errorbar(B,pw12,pw12E,'--m^','LineWidth',lw,'MarkerSize',mks,'MarkerFaceColor','m');
errorbar(B,pw6,pw6E,'--g>','LineWidth',lw,'MarkerSize',mks,'MarkerFaceColor','g');
errorbar(B,pw5mod,pw5modE,'--b<','LineWidth',lw,'MarkerSize',mks,'MarkerFaceColor','b');
errorbar(B4,pw4mod,pw4modE,'--rs','LineWidth',lw,'MarkerSize',mks,'MarkerFaceColor','r'); hold on;


legend({'EXP','d=1.2','d=0.6','d=0.5','d=0.4'},'Interpreter','latex');
%title('80% sync');
x0=0;
y0=0;
width=500;
height=400;
set(gcf,'position',[x0,y0,width,height])
set(gca,'FontSize',15);

fig=figure(5);
%saveas(fig,'/u/homes/np451/Dropbox/Synchonisation of mammalian cilia with hydrodynamic forces/PNAS/figures/v3/figure_sim/epsilon_simulations.pdf')

%xlabel('Number of Cilia/ N','FontSize',20);
%ylabel('$\epsilon $ [Hz]','FontSize',20);

%saveas(gcf,'/u/homes/np451/Desktop/sync_paper/figures/figure4/80pc.pdf')
%saveas(gcf,'80pc.pdf');



