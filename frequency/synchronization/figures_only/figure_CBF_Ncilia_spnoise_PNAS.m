%%%%% this script make figure for frequency vs N
%%%%% taking into account spatial noise and other staff

%%frequency vs noise and N

clear all
cd '/media/np451/Seagate Backup Plus Drive/DATA/synchronisation_paper';
%load('allvariables_area.mat');
load('all_variables_area_PNAS.mat');

%%% Number of cilia and Frequency with cells used with external flow
cc=1;FREQ=[];NCILIA=[];

lw=1;
mks=14
px2mu = 0.13*40/60;
nogood=[2];
areamin=15;
Ncells=numel(Cell)
for nc=1:Ncells;
    if any(nc==nogood); NCILIA(cc)= 0;
    else NCILIA(nc)= Cell(nc).Ncilia;
    end
    %plot(Cell(nc).Ncilia, Cell(nc).F_rest,'bo','LineWidth',lw,'MarkerSize',mks);hold on;
    leg{nc}=Cell(nc).Pos(end-2:end); 
    FREQ(nc)=Cell(nc).F_rest;
    AREA(nc)= Cell(nc).Sp.area_debris2*(px2mu^2); 
    
    
    
    F_temp = Cell(nc).Sp.F_debris2 ;
    [F,BW] = remove_debris(F_temp,areamin,[nanmedian(F_temp(:))*0.6,nanmedian(F_temp(:))*1.4]);
    
    MED(nc)=nanmedian(F(:)); %Cell(nc).Sp.medf;
    ESTD(nc)= nanstd(F(:));
    %    EMED(cc)=Cell(nc).Sp.medferr;

    %MED(cc)=Cell(nc).Sp.medf_debris2 %Cell(nc).Sp.medf;
%    ESTD(cc)=nanstd(F(:))
%    ESTD(cc)=Cell(nc).Sp.stdf_debris2%Cell(nc).Sp.stdf;
    %    ESTD_M(cc)=Cell(nc).Sp.stdf/sqrt(sum(~isnan(Cell(nc).Sp.F(:))));
    leg{nc}=num2str(nc); cc=cc+1;
end

% add this info to the other experiments with cells not use for flow experiments
cd '/media/np451/Seagate Backup Plus Drive/DATA/synchronisation_paper/14.12.18/frequency/'
load('FC_PNAS.mat');
nogood=[3,7,8,12,15,20,22,26,30,29,38,47];
%nogood=0;

cc=58;
for nc=1:numel(FC.Res)
    if any(nc==nogood); NCILIA(cc)= 0;
    else NCILIA(cc)= FC.Res(nc).Ncilia;
    end
    F_temp = FC.Res(nc).F_debris2; 
    [F,BW] = remove_debris(F_temp,areamin,[nanmedian(F_temp(:))*0.6,nanmedian(F_temp(:))*1.4]);
    MED(cc)=nanmedian(F(:)); %Cell(nc).Sp.medf;
    ESTD(cc)= nanstd(F(:));
    AREA(cc)= FC.Res(nc).area_debris2*(px2mu^2); 
    
%     FREQ(cc)=FC.Res(nc).f_guess
%     MED(cc)=FC.Res(nc).medf_debris2;
%     EMED(cc)=FC.Res(nc).medferr;
%     ESTD(cc)=FC.Res(nc).stdf_debris2;
%     ESTD_M(cc)=FC.Res(nc).stdf/sqrt(sum(~isnan(FC.Res(nc).F(:))));
%    leg{cc}=num2str(nc); ;
cc=cc+1;
end


% start plot

Eper= (ESTD./MED);
%Eper = AREA;
in=~isinf(Eper) & ~isnan(Eper) & (NCILIA<25) & (NCILIA~=0)% & onlyfirst %./sqrt(NCILIA);

half=median(Eper(in));
indice=1:numel(Eper)
%half= (max(Eper)+min(Eper))/2;
ind1= in & Eper< half ;
%ind1= in & indice<58
ind2= in &  Eper>half;
%ind2= in & indice>58
nbins=6;
figure(30);
scatter(NCILIA(in),MED(in),60,Eper(in),'filled');
colorbar()
%colormap( [0 0 1; 0.5 0  0.5 ; 1 0 0]);
c = colorbar();
c.Label.String = '\eta_{s}';



nbins=5;
[histw,ehistw,vinterval]=hist_nico_std(NCILIA(ind1),MED(ind1),nbins,[]);hold on;
errorbar(vinterval,histw,ehistw,'--kd','LineWidth',lw,'MarkerSize',mks,'MarkerFaceColor','k');
xlabel('Number of cilia per cell','FontSize',20);
ylabel('CBF [Hz]','FontSize',20);
ylim([10,32])

[histw,ehistw,vinterval2]=hist_nico_std(NCILIA(ind2),MED(ind2),nbins,[])%vinterval);
errorbar(vinterval(1:end-1),histw(1:end-1),ehistw(1:end-1),'--k^','LineWidth',lw,'MarkerSize',mks,'MarkerFaceColor','w');
xlabel('Number of cilia per cell','FontSize',20);
ylabel('CBF [Hz]','FontSize',20);
ylim([8,32])
xlim([0,23])
legend({'data','$\eta_{s}<\overline{\eta_{s}}$','$\eta_{s}> \overline{\eta_{s}}$'},'Location', 'nw','Interpreter','latex')

x0=0;
y0=0;
width=500;
height=400;
set(gcf,'position',[x0,y0,width,height])
set(gca,'FontSize',15);
%set(0,'defaulttextInterpreter','latex') ;

title(strcat('mean spatial noise used $\overline{\eta_{s}}=',num2str(half,3),'$'...,'\pm',num2str(eint_d_all,1)
    ),'Interpreter','latex','FontSize',9)

fig=figure(30);
%saveas(fig,'/u/homes/np451/Dropbox/Synchonisation of mammalian cilia with hydrodynamic forces/PNAS/figures_review/figure_CBF_cilia/CBF_spatial_noise.pdf')

%saveas(fig,'/u/homes/np451/Desktop/sync_paper/figures/figure2_new/frequency_df2.pdf')
%saveas(fig,'/u/homes/np451/Desktop/sync_paper/figures/figure3_new/frequency_df2.pdf')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% make unique average %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% This is for le plot of the figure in the paper %%
%-------- CBF vs number of cilia with theoretical line from hydrodynamic
%screening ------------------------------------------------------------
in=~isinf(Eper) & ~isnan(Eper) & (NCILIA~=0) & (NCILIA<25)% & onlyfirst %./sqrt(NCILIA);
%Eper= (100*EMED./FREQ);
%half=median(Eper(in));
%half= (max(Eper)+min(Eper))/2;

nbins=5;
figure(55);
scatter(NCILIA(in),MED(in),60,Eper(in),'filled');
colorbar()
%colormap( [0 0 1; 0.5 0  0.5 ; 1 0 0]);
c = colorbar();
%c.Label.String = '\eta_{s}';
hold on;


nbins=6;
[histw,ehistw,vinterval]=hist_nico_std(NCILIA(in),MED(in),nbins,[]);hold on;
errorbar(vinterval,histw,ehistw,'kd','LineWidth',lw,'MarkerSize',mks,'MarkerFaceColor','k');

%%%%%%%%%%%%% linear fit %%%%%%%%%%%%%%%%%%%%%%%%
[p,S]=polyfit(NCILIA(in),MED(in),3);
%plot(1:25,polyval(p,1:25),'g-','LineWidth',lw)
y=MED(in);
r_squared= 1 - (S.normr/norm(y - mean(y)))^2


% x=NCILIA(in);
% y=MED(in);
% ft_CBF = fittype('a*x+b','independent','x');
% fo_CBF = fitoptions('Method','NonLinearLeastSquares',...
%                'Lower',[0,0],...
%                'Upper',[3,20],...
%                'StartPoint',[0.3,13],'Robust','Bisquare');
% [fit_out_CBF,gof,output] = fit(x',y',ft_CBF,fo_CBF);
% efit= confint(fit_out_CBF); ea= abs(efit(1)-efit(2));
% %ea= abs(efit(1,2)-efit(2,2));
% plot(fit_out_CBF)



%%% model with N
N=1:25;
C= 0.07*10^-6;   %%% um
Fdr= 6*10^-12;  %%%% pN

r_e=C.*sqrt(N);
L=11*10^-6;   %%%%% lenght cilium
a=(1*10^-6);
r0= (1*10^-7)
eta=0.8*10^-3;%%%%% viscosity water
fd=15;        %%%%% frequency single cilium
gamma_0= (4*pi*eta*L)/(log(L./r0)-(1/2));
Amp=10*10^-6;
gamma_N= (4*pi*eta*L)./(log(L./r_e)-(1/2));
%CBF= N*Fdr./(2*Amp*gamma_N);
CBF= (N*fd*gamma_0)./(gamma_N);

hold on;
plot(N,CBF,'r-','LineWidth',3);

%%%plot trend of N
%plot(5:25,polyval([1,10],5:25),'k-','LineWidth',2);
%plot(5:25,polyval([1,0],5:25),'k-','LineWidth',2);


%legend({'data','binning average','linear fit','IER model'},'Location', 'ne','Interpreter','latex');
legend({'data','binning average','IER model'},'Location', 'ne','Interpreter','latex');
%xlabel('Number of cilia per cell','FontSize',20);
%ylabel('CBF [Hz]','FontSize',20);
ylim([10,40])
xlim([0,26])
x0=0;
y0=0;
width=500;
height=400;
set(gcf,'position',[x0,y0,width,height])
set(gca,'FontSize',15);
%set(0,'defaulttextInterpreter','latex') ;
% set(gca,'XScale', 'log', 'YScale', 'log');
set(gca,'XScale', 'Lin', 'YScale', 'Lin');

%title(strcat('average CBF Ncilia'...,'\pm',num2str(eint_d_all,1)
%    ),'Interpreter','latex','FontSize',9)

fig=figure(55);
%saveas(fig,'/u/homes/np451/Dropbox/Synchonisation of mammalian cilia with hydrodynamic forces/PNAS/figures/v3/figure_CBF/CBF_modelcompare.pdf')
%saveas(fig,'/u/homes/np451/Desktop/sync_paper/figures/figure3_new/frequency_df2.pdf')
%%
%%%%%%%%%%%%%%%%%%%%%%% PLOT for REFEREE 1 %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% FREQUENCY vs SPATIAL NOIsE %%%%%%%%%%%%%%%%%%%%%

%onlyfirst(57:end)=0;
in=~isinf(Eper) & ~isnan(Eper) & (NCILIA~=0) & (NCILIA<25)% & onlyfirst %./sqrt(NCILIA);
%Eper= (100*EMED./FREQ);
%half=median(Eper(in));
%half= (max(Eper)+min(Eper))/2;

nbins=6;
figure(56);
scatter(Eper(in),MED(in),60,NCILIA(in),'filled');
colorbar()
%colormap( [0 0 1; 0.5 0  0.5 ; 1 0 0]);
c = colorbar();
c.Label.String = 'cilia number N';
hold on;

nbins=6;
[histw,ehistw,vinterval]=hist_nico_std(Eper(in),MED(in),nbins,[]);hold on;
errorbar(vinterval,histw,ehistw,'kd','LineWidth',lw,'MarkerSize',mks,'MarkerFaceColor','k');

%title(strcat('average CBF vs spatial noise'...,'\pm',num2str(eint_d_all,1)
%    ),'Interpreter','latex','FontSize',9)

xlabel('Spatial noise','FontSize',15);
ylabel('CBF [Hz]','FontSize',15);
legend('raw data','binning median');

fig=figure(56);
%saveas(fig,'/u/homes/np451/Dropbox/Synchonisation of mammalian cilia with hydrodynamic forces/PNAS/figures/v3/CBF_vs_spatial_noise.png')
%saveas(fig,'/u/homes/np451/Desktop/sync_paper/figures/figure3_new/frequency_df2.pdf')

indice=1:numel(Eper);
 m15=median(Eper(in & NCILIA>20 & indice <58));
 em15=std(Eper(in & NCILIA>20 & indice <58));
 disp(strcat('spatial noise for N>15 =',num2str(m15),' pom',num2str(em15)));
 
 m15=median(Eper(in & NCILIA<20 & indice <58));
 em15=std(Eper(in & NCILIA<20 & indice <58));
 disp(strcat('spatial noise for N<15 =',num2str(m15),' pom',num2str(em15)));
  
