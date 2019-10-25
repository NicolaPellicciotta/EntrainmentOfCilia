%%%%%%%%%%%%%%%%%% this is the script to make figure for PNAS
%%%%%%%%%%%%%%%%%% epilon/v fitted with different effective radius

%% plot epsilon/v vs N with fit
%  epsilon/v_ex fitted using all data instead of the average
%Cal=load('/home/np451/Documents/MATLAB/frequency/synchronization/calibration/calibration_matrix.mat')

clear all;
C0_array=[0.67*10^(-6),...  %%%%% amplitude model C=a/sqrt(pi)
    0.11*10^(-6)];          %%%% touching cilia   
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

%%%%%%%%%%%%load data %%%%%%%%%%%%%%%%%%%%%%%%%

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
 

%%%%%%%%%%%%%%%%%%% set data for fit %%%%%%%%%%%%%%%%%%%%%%%%% 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
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

%%%%%%%%%%%%%%%% FIRST FIT MODEL %%%%%%%%%%%%%%%%%%%%%%%

C0=C0_array(1);
%C0=0.67*10^(-6); %%%% with reviewer model 
%C0=0.52*10^(-6); %%%% with centrin model and intercept different from 0. 
%C0=0.11*10^(-6); %%%% ideal case when all cilia are in contact

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
[fit_out_all,gof,output] = fit(NCILIA_ALL(~isnan(EPS_V))',EPS_V((~isnan(EPS_V)))',ft_all,fo_all);
GOF{1}=gof;



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


Fd_all(1)= 2*fd*4*pi*eta*L/(fit_out_all.A*10^3);
eFd_all(1) = (10^12)*((2*fd*4*pi*eta*L/10^3)/fit_out_all.A^2)* eA;  %%% error in pN



%l{cleg}='$\frac{[2 (CBF_{0})/F_{dr}]4\pi\eta L}{ N[\log{(L/r_{e}(N))}-1/2]}$';
l{cleg}='IER model $r_{e}$(N)';
cleg=cleg+1;
%%%%%%%%%%%%%%%%%%%%% second plot with other model %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

C0=C0_array(2);
%C0=0.67*10^(-6); %%%% with reviewer model 
%C0=0.52*10^(-6); %%%% with centrin model and intercept different from 0. 
%C0=0.11*10^(-6); %%%% ideal case when all cilia are in contact


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
[fit_out_all,gof,output] = fit(NCILIA_ALL(~isnan(EPS_V))',EPS_V((~isnan(EPS_V)))',ft_all,fo_all);
GOF{2}=gof;
efit_out_all=confint(fit_out_all);
eA=abs(efit_out_all(1,1)-efit_out_all(2,1));
eB=abs(efit_out_all(1,2)-efit_out_all(2,2));
if isnan(eB); eB=abs(B0u-B0l); end

nbins=5;
%errorbar(vinterval,histw,ehistw*2,color_median,'MarkerSize',mks,'LineWidth',lw,'MarkerFaceColor',color_face);hold on;
%plot(NCILIA,EPS_V,'ko','MarkerSize',mks-3,'LineWidth',lw);hold on;
%errorbar(vinterval,y,ey,'ko','MarkerSize',mks,'LineWidth',lw); hold on;
plot(2:30,fit_out_all(2:30),'g-','MarkerSize',mks,'LineWidth',lw+2);
set(gca,'XScale', 'log', 'YScale', 'log');

%int_d_all= L/(sqrt(fit_out_all.B)*exp(1)^(1/2));
%eint_d_all= (10^6)*(L/(2*exp(1)^(1/2)))*fit_out_all.B^(-3/2)* eB ;
int_d_all= (L*sqrt(fit_out_all.B))/(em);
eint_d_all= (10^6)*(L/(2*exp(1)^(1/2)))*fit_out_all.B^(-3/2)* eB ;


Fd_all(2)= 2*fd*4*pi*eta*L/(fit_out_all.A*10^3);
eFd_all(2) = (10^12)*((2*fd*4*pi*eta*L/10^3)/fit_out_all.A^2)* eA;  %%% error in pN


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A_exp=(log( L/a )-1/2)/fit_out_all.A;
Fd_theo= (Amp*2*fd)*4*pi*eta*L /(log(L/a)-1/2);


%l{cleg}='$\frac{A}{\frac{N}{2}\log{(BN)}}$';
%l{cleg}='$\frac{A}{ N[\log{L/r_{e}(N)}-1/2}$';
%l{cleg}='$\frac{[2 (CBF_{0})/F_{dr}]4\pi\eta L}{ N[\log{(L/r2_{e}(N))}-1/2]}$';
l{cleg}='IER model $r_{lim}$(N)';


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
%title(strcat('ciliary distance $C=',num2str(int_d_all*10^6,2)...,'\pm',num2str(eint_d_all,1)
%    ,'\mu m$ and $F_{d}=',num2str(Fd_all(1)*10^12,2),'\pm',num2str(eFd_all(1),2),'pN$',...
%    ' $F_{d}=',num2str(Fd_all(2)*10^12,2),'\pm',num2str(eFd_all(2),2),'pN$'),'Interpreter','latex','FontSize',9)
set(gca,'FontSize',15);
fig=figure(12);
 disp(strcat('r-squared for centrin =',num2str(GOF{1}.rsquare),'r-squared for very close data =',num2str(GOF{2}.rsquare)) )


%saveas(fig,'/u/homes/np451/Dropbox/Synchonisation of mammalian cilia with hydrodynamic forces/PNAS/figures/v3/figure3/epsilon_v_N.pdf')

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
%plot(nmeno1x,nmeno1y,'k-','MarkerSize',mks,'LineWidth',lw+2)


%set(gca,'XScale', 'log', 'YScale', 'log');
legend(l{1:end-1},'FontSize',13,'Interpreter','latex');

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

