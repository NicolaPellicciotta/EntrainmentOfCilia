%%%% scrip for figure 2 new

%% plot effective radius vs number of cilia
 load('/media/np451/Seagate Backup Plus Drive/DATA/synchronisation_paper/allvariables_area.mat');
set(0,'defaulttextInterpreter','latex') ;
lw=1;
mks=7;

px2mu=0.0973;
AREA=[]
for nc=1:numel(Cell)
    
    
    
   AREA(nc)= Cell(nc).Sp.area*(px2mu^2);
   
   
   %plot(Cell(nc).Ncilia,Cell(nc).Sp.area*px2mu^2/10,'o');hold on   
end

Norm= AREA(NCILIA==1);
AREA=AREA;
R_eff= sqrt(AREA./pi);
Norm= R_eff(NCILIA==1);
R_eff=0.1*R_eff/Norm

figure(10)
plot(NCILIA,AREA,'ko','MarkerSize',mks,'LineWidth',lw);
%%%%% fit with polyfit
p=polyfit(log(NCILIA),log(AREA),1);hold on;
plot(1:25,exp(polyval(p,log(1:25))),'r-','MarkerSize',mks,'LineWidth',lw)
title('AREA')
xlabel('Ncilia');
ylabel('Area [um2]');
legend('exp','y=Ax^{\alpha}');
set(gca,'XScale', 'log', 'YScale', 'log');

%%%%% plot effective radius

figure(11)
plot(NCILIA,R_eff,'ko','MarkerSize',mks,'LineWidth',lw);
%%%%% fit with polyfit
p_radius=polyfit(log(NCILIA),log(R_eff),1);hold on;
%plot(1:25,exp(polyval(p_radius,log(1:25))),'r-','MarkerSize',mks,'LineWidth',lw)

%%%%% fit with matlab fit
ft_radius = fittype('(1/2)*x+C','independent','x');
fo_radius = fitoptions('Method','NonLinearLeastSquares',...
               'Lower',[-10],...
               'Upper',[0],...
               'StartPoint',[-2.6]);
fit_out_radius = fit(log(NCILIA)',log(R_eff)',ft_radius,fo_radius);
plot(1:25,exp(fit_out_radius(log(1:25))),'r-','MarkerSize',mks,'LineWidth',lw)
int_c_exp= exp(fit_out_radius.C)
eC= confint(fit_out_radius); eC= abs(eC(1)-eC(2));
eint_c_exp= abs(int_c_exp*fit_out_radius.C*eC);


xlabel('Ncilia'); ylim([0.05,1]);
ylabel('Effective radius[um]');
legend({'data','$y=CN^{1/2}$'},'Interpreter','latex');
%title(strcat('A=',num2str(exp(p_radius(2))),'   and $\alpha$= ',num2str(p_radius(1))),'Interpreter','latex', 'FontSize',15);
%title(strcat('$A=',num2str(exp(fit_out_radius.C),2),'\mum$   and $\alpha$= ',num2str(p_radius(1))),'Interpreter','latex', 'FontSize',15);
title(strcat('$A=',num2str(exp(fit_out_radius.C),1),'\pm',num2str(eint_c_exp,1),'\mu m$'),'Interpreter','latex', 'FontSize',15);
set(gca,'XScale', 'log', 'YScale', 'log');
x0=0;
y0=0;
width=500;
height=400;
set(gcf,'position',[x0,y0,width,height])
set(gca,'FontSize',15);

fig=figure(11)
%saveas(fig,'/u/homes/np451/Desktop/sync_paper/figures/figure2_new/radius_eff_area.pdf')
%% plot epsilon/v vs N with fit
%  epsilon/v_ex fitted using all data instead of the average
%Cal=load('/home/np451/Documents/MATLAB/frequency/synchronization/calibration/calibration_matrix.mat')
load('/media/np451/Seagate Backup Plus Drive/DATA/synchronisation_paper/allvariables_area.mat');
close all
set(0,'defaulttextInterpreter','latex') ;
Volts_plot=[2,3,4,5]; lw=0.3; mks=14;

color_scatter_array={'ko','bo','ro','mo'};
color_median_array={'k--o','r--s','b--<','m--^'};
color_face_array={'k','r','b','m'};
clear l;cleg=1;

HISTW=nan(4,5);
EHISTW=nan(4,5);
cc=1;
EPS_V=[]; NCILIA_ALL=[]; VOLTAGE=[];
 for jj=1:4 
Flow_mean= mean(Cal.Flow_int,2);
whichVolt=Volts_plot(jj);
v_eps= (Volts == whichVolt);
l{cleg}= strcat(num2str(Flow_mean(v_eps),'%-5.2f'),' mm/s'); cleg=cleg+1;

ms=11;lw=1;
color_scatter=color_scatter_array{jj};color_median=color_median_array{jj};color_face=color_face_array{jj};

for nc=1:Ncells;
    if Cell(nc).good
        if whichVolt<4 & nc>28
%        plot(Cell(nc).Ncilia,Cell(nc).eps(v_eps),color_scatter,'MarkerSize',ms,'LineWidth',lw,'MarkerFaceColor',[1,1,1]);hold on;
%        leg{cc}=num2str(nc); 
        if Cell(nc).eps(v_eps)<12; EPS_V(cc)=Cell(nc).eps(v_eps)/(Flow_mean(v_eps));
        else EPS_V(cc)=nan;end
        NCILIA_ALL(cc)=Cell(nc).Ncilia;
        VOLTAGE(cc)=jj;
        cc=cc+1;
        end

        if whichVolt>=4 
%        plot(Cell(nc).Ncilia,Cell(nc).eps(v_eps),color_scatter,'MarkerSize',ms,'LineWidth',lw,'MarkerFaceColor',[1,1,1]);hold on;
%        leg{cc}=num2str(nc);
        if Cell(nc).eps(v_eps)<12; EPS_V(cc)=Cell(nc).eps(v_eps)/(Flow_mean(v_eps));
        else EPS_V(cc)=nan;end
        NCILIA_ALL(cc)=Cell(nc).Ncilia;
        VOLTAGE(cc)=jj;
        cc=cc+1;
        end
        
        

    end
end

end

Volts_plot=[2,3,4,5];
lw=0.3;
mks=14;

color_scatter_array={'ko','bo','ro','mo'};
%color_scatter_array={'ko','ko','ko','ko'};
color_median_array={'ko','rs','b<','m^'};
color_face_array={'k','r','b','m'};
clear l;cleg=1;

HISTW=nan(4,5);
EHISTW=nan(4,5);
 for jj=1:4 
 
Flow_mean= mean(Cal.Flow_int,2);

whichVolt=Volts_plot(jj);
v_eps= (Volts == whichVolt);
l{cleg}= strcat('$',num2str(Flow_mean(v_eps),'%-5.2f'),'$ mm/s'); cleg=cleg+1;
cc=1;
EPS=[]; NCILIA=[]; 
ms=11;
%lw=0.5;
lw=1;
color_scatter=color_scatter_array{jj};
color_median=color_median_array{jj};
color_face=color_face_array{jj};

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
nbins=5;
figure(12);
%[histw,ehistw,vinterval]=hist_nico(NCILIA,EPS,nbins,[1:5:30]);
[histw,ehistw,vinterval]=hist_nico_std(NCILIA,EPS,[],[1:5:30])

histw(histw>9)=nan;
ehistw(histw>9)=nan;
ehistw(ehistw<0.4)=0.4;
%%%% to remove from statistic the one the went out of range 
%errorbar(vinterval,histw,ehistw,color_median,'MarkerSize',mks,'LineWidth',lw,'MarkerFaceColor',color_face);hold on;
errorbar(vinterval,histw./(Flow_mean(v_eps)),ehistw./(Flow_mean(v_eps)),color_median,'MarkerSize',mks,'LineWidth',lw,'MarkerFaceColor',color_face);hold on;
disp(jj)
HISTW(jj,:)= histw'./(Flow_mean(v_eps));
EHISTW(jj,:)=ehistw'./(Flow_mean(v_eps))
end
 
L=11*10^-6;   %%%%% lenght cilium
a=(1*10^-6)
eta=0.8*10^-3;%%%%% viscosity water
fd=15%15;        %%%%% frequency single cilium
Fd=33*10^-12;  %%%%% force acting on cilium in piconewton
Amp=10*10^-6;

A0= 2*(fd/Fd)*4*pi*eta*L*10^-3;     %%%% in mm
A0=10;


%%%%%%%% setting the reff constant
%C0=0.07*10^(-6);
C0=0.67*10^(-6); %%%% with reviewer model 
%C0=0.52*10^(-6); %%%% with centrin model and intercept different from 0. 
C0=0.11*10^(-6); %%%% ideal case when all cilia are in contact


em= exp(1)^(1/2)
B0= (L/(em*C0))^-2; 
B0l= (L/(em*C0*0.8))^-2; 
B0u= (L/(em*C0*1.2))^-2; 

eEPS_V=ones(size(EPS_V(~isnan(EPS_V))));

ft_all = fittype('A/(-N*(0.5*log(B*N)))','independent','N');
fo_all = fitoptions('Method','NonLinearLeastSquares',...
               'Lower',[A0/10,B0l],...
               'Upper',[A0*100,B0u],...
               'StartPoint',[A0,B0],'Weights',eEPS_V,'Robust','Bisquare');
fit_out_all = fit(NCILIA_ALL(~isnan(EPS_V))',EPS_V((~isnan(EPS_V)))',ft_all,fo_all);

efit_out_all=confint(fit_out_all);
eA=abs(efit_out_all(1,1)-efit_out_all(2,1));
eB=abs(efit_out_all(1,2)-efit_out_all(2,2));
if isnan(eB); eB=abs(B0u-B0l); end

nbins=5;
%errorbar(vinterval,histw,ehistw*2,color_median,'MarkerSize',mks,'LineWidth',lw,'MarkerFaceColor',color_face);hold on;
%plot(NCILIA,EPS_V,'ko','MarkerSize',mks-3,'LineWidth',lw);hold on;
%errorbar(vinterval,y,ey,'ko','MarkerSize',mks,'LineWidth',lw); hold on;
plot(2:30,fit_out_all(2:30),'r-','MarkerSize',mks,'LineWidth',lw+2);
set(gca,'XScale', 'log', 'YScale', 'log');

%int_d_all= L/(sqrt(fit_out_all.B)*exp(1)^(1/2));
%eint_d_all= (10^6)*(L/(2*exp(1)^(1/2)))*fit_out_all.B^(-3/2)* eB ;
int_d_all= (L*sqrt(fit_out_all.B))/(em);
eint_d_all= (10^6)*(L/(2*exp(1)^(1/2)))*fit_out_all.B^(-3/2)* eB ;


Fd_all= 2*fd*4*pi*eta*L/(fit_out_all.A*10^3);
eFd_all = (10^12)*((2*fd*4*pi*eta*L/10^3)/fit_out_all.A^2)* eA;  %%% error in pN



A_exp=(log( L/a )-1/2)/fit_out_all.A;
Fd_theo= (Amp*2*fd)*4*pi*eta*L /(log(L/a)-1/2);


%l{cleg}='$\frac{A}{\frac{N}{2}\log{(BN)}}$';
%l{cleg}='$\frac{A}{ N[\log{L/r_{e}(N)}-1/2}$';
l{cleg}='$\frac{[2 (CBF_{0})/F_{dr}]4\pi\eta L}{ N[\log{(L/r_{e}(N))}-1/2]}$';

%plot the -1 power line
nmeno1x=2:10;
nmeno1y=5*nmeno1x.^(-1)
plot(nmeno1x,nmeno1y,'k--','MarkerSize',mks,'LineWidth',lw+2)

legend(l,'FontSize',13,'Interpreter','latex');
set(gca,'XScale', 'log', 'YScale', 'log');
%xlabel('$N_{c}$ number of cilia','Interpreter','latex', 'FontSize',20);
%ylabel('$\epsilon/v_{EX} [mm]^{-1}$ ','Interpreter','latex', 'FontSize',20);

x0=100;
y0=100;
width=500;
height=400;
set(gcf,'position',[x0,y0,width,height])
set(gca,'FontSize',15);
title(strcat('ciliary distance $C=',num2str(int_d_all*10^6,2)...,'\pm',num2str(eint_d_all,1)
    ,'\mu m$ and $F_{d}=',num2str(Fd_all*10^12,2),'\pm',num2str(eFd_all,2),'pN$'),'Interpreter','latex','FontSize',9)
fig=figure(12);
%saveas(fig,'/u/homes/np451/Desktop/sync_paper/figures/figure3_new/epsilon_v_n_std_lin.pdf')

%% figure for supplementaries reviwer1 pnas
figure(15)
mks=8;
lw=1;
color_median_array={'ko','rs','b<','m^'};
for tt=1:4
    COLOR=(VOLTAGE==tt);
    color_median=color_median_array{tt};
    color_face=color_face_array{tt};
    
%errorbar(NCILIA_ALL(~isnan(EPS_V) & COLOR),EPS_V((~isnan(EPS_V) & COLOR)),color_median,'MarkerSize',mks,'LineWidth',lw,'MarkerFaceColor',color_face)
plot(NCILIA_ALL(~isnan(EPS_V) & COLOR),EPS_V((~isnan(EPS_V) & COLOR)),color_median,'MarkerSize',mks,'LineWidth',lw,'MarkerFaceColor',color_face)

hold on;
end


x0=100;
y0=100;
width=700;
height=600;
set(gcf,'position',[x0,y0,width,height])
set(gca,'FontSize',15);
plot(1:30,fit_out_all(1:30),'r-','MarkerSize',mks,'LineWidth',lw+2);

%plot the -1 power line
nmeno1x=2:10;
nmeno1y=5*nmeno1x.^(-1)
plot(nmeno1x,nmeno1y,'k-','MarkerSize',mks,'LineWidth',lw+2)


%set(gca,'XScale', 'log', 'YScale', 'log');
legend(l,'FontSize',13,'Interpreter','latex');

x0=0;
y0=0;
width=500;
height=400;
set(gcf,'position',[x0,y0,width,height])
set(gca,'FontSize',15);


%title(strcat('ciliary distance $C=',num2str(int_d_all*10^6,2)...,'\pm',num2str(eint_d_all,1)
%    ,'\mu m$ and $F_{d}=',num2str(Fd_all*10^12,2),'\pm',num2str(eFd_all,2),'pN$'),'Interpreter','latex','FontSize',9)
fig=figure(15);
%saveas(fig,'/u/homes/np451/Desktop/sync_paper/figures/figure3_new/suppl_epsilon_v_n.pdf')
saveas(fig,'/u/homes/np451/Desktop/sync_paper/figures/figure3_new/suppl_epsilon_v_n_a_1.2um.pdf')



%% frequency vs noise and N

clear all
cd '/media/np451/Seagate Backup Plus Drive/DATA/synchronisation_paper';
%load('allvariables_area.mat');
load('all_variables_area_PNAS.mat');

%%% Number of cilia and Frequency with cells used with external flow
cc=1;FREQ=[];NCILIA=[];

lw=1;
mks=14

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
    
%     FREQ(cc)=FC.Res(nc).f_guess;
%     MED(cc)=FC.Res(nc).medf_debris2;
%     EMED(cc)=FC.Res(nc).medferr;
%     ESTD(cc)=FC.Res(nc).stdf_debris2;
%     ESTD_M(cc)=FC.Res(nc).stdf/sqrt(sum(~isnan(FC.Res(nc).F(:))));
%    leg{cc}=num2str(nc); ;
cc=cc+1;
end


% start plot

Eper= (ESTD./MED);
%onlyfirst= ones([1,numel(NCILIA)]);
%onlyfirst(57:end)=0;
in=~isinf(Eper) & ~isnan(Eper) & (NCILIA<25) & (NCILIA~=0)% & onlyfirst %./sqrt(NCILIA);
%Eper= (100*EMED./FREQ);
half=median(Eper(in));
%half= (max(Eper)+min(Eper))/2;
ind1= in & Eper<half ;

ind2= in &  Eper>half;
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
%saveas(fig,'/u/homes/np451/Desktop/sync_paper/figures/figure2_new/frequency_df2.pdf')
%saveas(fig,'/u/homes/np451/Desktop/sync_paper/figures/figure3_new/frequency_df2.pdf')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% make unique average %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%onlyfirst(57:end)=0;
in=~isinf(Eper) & ~isnan(Eper) & (NCILIA~=0) & (NCILIA<25)% & onlyfirst %./sqrt(NCILIA);
%Eper= (100*EMED./FREQ);
%half=median(Eper(in));
%half= (max(Eper)+min(Eper))/2;

nbins=6;
figure(55);
scatter(NCILIA(in),MED(in),60,Eper(in),'filled');
colorbar()
%colormap( [0 0 1; 0.5 0  0.5 ; 1 0 0]);
c = colorbar();
c.Label.String = '\eta_{s}';



nbins=5;
[histw,ehistw,vinterval]=hist_nico_std(NCILIA(in),MED(in),nbins,[]);hold on;
errorbar(vinterval,histw,ehistw,'kd','LineWidth',lw,'MarkerSize',mks,'MarkerFaceColor','k');
xlabel('Number of cilia per cell','FontSize',20);
ylabel('CBF [Hz]','FontSize',20);


%%% model with N
N=1:25;
C= 0.07*10^-6;   %%% um
Fdr= 6*10^-12;  %%%% pN

r_e=C.*sqrt(N);
L=11*10^-6;   %%%%% lenght cilium
a=(1*10^-6)
eta=0.8*10^-3;%%%%% viscosity water
fd=15;        %%%%% frequency single cilium
gamma_0= (4*pi*eta*L)/(log(L./r_e(1))-(1/2));
Amp=10*10^-6;
gamma_N= (4*pi*eta*L)./(log(L./r_e)-(1/2));
%CBF= N*Fdr./(2*Amp*gamma_N);
CBF= (N*fd*gamma_0)./(gamma_N);

hold on;
plot(N,CBF,'r-','LineWidth',3);

%%%plot trend of N
%plot(5:25,polyval([1,10],5:25),'k-','LineWidth',2);
%plot(5:25,polyval([1,0],5:25),'k-','LineWidth',2);


legend({'data','binning average','IR model'},'Location', 'nw','Interpreter','latex')
ylim([10,100])
xlim([0,26])
x0=0;
y0=0;
width=500;
height=400;
set(gcf,'position',[x0,y0,width,height])
set(gca,'FontSize',15);
%set(0,'defaulttextInterpreter','latex') ;
 set(gca,'XScale', 'log', 'YScale', 'log');

title(strcat('average CBF Ncilia'...,'\pm',num2str(eint_d_all,1)
    ),'Interpreter','latex','FontSize',9)

fig=figure(55);
%saveas(fig,'/u/homes/np451/Desktop/sync_paper/figures/figure2_new/frequency_df2.pdf')
%saveas(fig,'/u/homes/np451/Desktop/sync_paper/figures/figure3_new/frequency_df2.pdf')

