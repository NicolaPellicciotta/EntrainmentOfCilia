%%%%%%%% this is a scrip to find the cell area and density
N_cells=numel(Cell)
for nc=1:N_cells
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
metarray=[1000,1,2]
[Meta] = spatial_noise_fft_meta(f_guess,box_size,[],[],metarray);

Cell(nc).Sp.F=Meta.F;
Cell(nc).Sp.P=Meta.P;
Cell(nc).Sp.P_meanf=Meta.P_meanf;
Cell(nc).Sp.box_size=box_size;
Cell(nc).Sp.s_roi=Meta.s_roi
Cell(nc).Sp.area= sum(~isnan(Meta.F(:)))*box_size^2;
Cell(nc).Sp.medf=nanmedian(Meta.F(:));
Cell(nc).Sp.medferr= nanmedian(abs(Meta.F(:)-nanmedian(Meta.F(:))));
Cell(nc).Sp.stdf=nanstd(Meta.F(:));
Cell(nc).Sp.rect_out=Meta.rect_out;



%%% look for a synchronised cell and caluclate meta for it
ind=[];cc=1;
for rc=1:numel(Cell(nc).Rc)
if Cell(nc).Rc{rc}.Period~=0 & Cell(nc).Rc{rc}.syn==1 & Cell(nc).Rc{rc}.Volt==4;
    ind(cc)= rc;cc=cc+1;
end
end

movie=strcat('V4*P',num2str(Cell(nc).Rc{ind(1)}.Period));
[Meta_syn] = spatial_noise_fft_meta(1000/Cell(nc).Rc{ind(1)}.Period,box_size,movie,Cell(nc).Sp.rect_out,metarray);
  
end

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
legend('data','y=A N^{1/2}');
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

fig=figure()
%saveas(fig,'/u/homes/np451/Desktop/sync_paper/figures/figure2_new/radius_eff_area.pdf')

%%
%%%%% plot epsilon/velocity as a function of the effective radius
%%%%% epsilon vs N_cilia

set(0,'defaulttextInterpreter','latex') ;

Volts_plot=[2,3,4,5];
lw=0.3;
mks=14;

color_scatter_array={'ko','bo','ro','mo'};
%color_scatter_array={'ko','ko','ko','ko'};
color_median_array={'k--o','r--s','b--<','m--^'};
color_face_array={'k','r','b','m'};
clear l;cleg=1;

HISTW=nan(4,5);
EHISTW=nan(4,5);
 for jj=1:4 
 
Flow_mean= mean(Cal.Flow_int,2);

whichVolt=Volts_plot(jj);
v_eps= (Volts == whichVolt);
l{cleg}= strcat(num2str(Flow_mean(v_eps),'%-5.2f'),' mm/s'); cleg=cleg+1;
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
figure(5);
[histw,ehistw,vinterval]=hist_nico(NCILIA,EPS,nbins,[1:5:30]);
histw(histw>9)=nan;
ehistw(histw>9)=nan;
%%%% to remove from statistic the one the went out of range 
%errorbar(vinterval,histw,ehistw,color_median,'MarkerSize',mks,'LineWidth',lw,'MarkerFaceColor',color_face);hold on;
ehistw(ehistw<0.4)=0.4;
errorbar(vinterval,histw./(Flow_mean(v_eps)),ehistw./(Flow_mean(v_eps)),color_median,'MarkerSize',mks,'LineWidth',lw,'MarkerFaceColor',color_face);hold on;
disp(jj)
HISTW(jj,:)= histw'./(Flow_mean(v_eps));
EHISTW(jj,:)=ehistw'./(Flow_mean(v_eps))
end
 
y= nanmean(HISTW,1); 
ey=nanmean(EHISTW,1)/sqrt(4);

L=11*10^-6;   %%%%% lenght cilium
a=(1*10^-6)
eta=0.8*10^-3;%%%%% viscosity water
fd=15;        %%%%% frequency single cilium
Fd=33*10^-12;  %%%%% force acting on cilium in piconewton
Amp=10*10^-6;

A0= 2*(fd/Fd)*4*pi*eta*L*10^-3;
A0=10;
B0= (L/(10^-7))^(1/0.5);

ft = fittype('A/(N*(0.5*log(B*N)))','independent','N');
fo = fitoptions('Method','NonLinearLeastSquares',...
               'Lower',[A0/10,(B0-B0*0.7)],...
               'Upper',[A0*100,(B0+B0*0.7)],...
               'StartPoint',[A0,B0],'Weights',1./(ey(1:end).^2),'Robust','Bisquare');
fit_out = fit(vinterval(1:end)',y(1:end)',ft,fo);
efit_out=confint(fit_out);
eA=abs(efit_out(1,1)-efit_out(2,1));
eB=abs(efit_out(1,2)-efit_out(2,2));
if isnan(eB); eB=B0*1.4; end
%errorbar(vinterval,y,ey,'ko','MarkerSize',mks,'LineWidth',lw); hold on;
plot(2:30,fit_out(2:30),'k-','MarkerSize',mks,'LineWidth',lw+2);
%set(gca,'XScale', 'log', 'YScale', 'log');

int_d_exp= L/(sqrt(fit_out.B)*exp(1)^(1/2));
eint_d_exp= (10^6)*(L/(2*exp(1)^(1/2)))*fit_out.B^(-3/2)* eB ;

Fd_exp= 2*fd*4*pi*eta*L/(fit_out.A*10^3);
eFd_exp = (10^12)*((2*fd*4*pi*eta*L/10^3)/fit_out.A^2)* eA;  %%% error in pN

A_exp=(log( L/a )-1/2)/fit_out.A;
Fd_theo= (Amp*2*fd)*4*pi*eta*L /(log(L/a)-1/2);

% errorbar(vinterval,y,ey,'ko','MarkerSize',mks,'LineWidth',lw); hold on;
% p_value= polyfit(log(vinterval'),log(y'),1)
% plot(2:30,exp(polyval(p_value,log(2:30))),'r-','MarkerSize',mks,'LineWidth',lw)
% set(gca,'XScale', 'log', 'YScale', 'log');
l{cleg}='$\frac{A}{\frac{N}{2}\log{(BN)}}$';
 legend(l,'FontSize',13,'Interpreter','latex');
 set(gca,'XScale', 'log', 'YScale', 'log');
xlabel('$N_{c}$ number of cilia','Interpreter','latex', 'FontSize',20);
ylabel('$\epsilon/v_{EX} [mm]^{-1}$ ','Interpreter','latex', 'FontSize',20);
title(strcat('ciliary distance $C=',num2str(int_d_exp*10^6,1),'\pm',num2str(eint_d_exp,1),'\mu m$ and $F_{d}=',num2str(Fd_exp*10^12,2),'\pm',num2str(eFd_exp,2),'pN$'),'Interpreter','latex')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  epsilon/v_ex fitted using all data instead of the average

set(0,'defaulttextInterpreter','latex') ;
Volts_plot=[2,3,4,5]; lw=0.3; mks=14;

color_scatter_array={'ko','bo','ro','mo'};
color_median_array={'k--o','r--s','b--<','m--^'};
color_face_array={'k','r','b','m'};
clear l;cleg=1;

HISTW=nan(4,5);
EHISTW=nan(4,5);
cc=1;
EPS_V=[]; NCILIA_ALL=[]; 
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
        cc=cc+1;
        end

        if whichVolt>=4 
%        plot(Cell(nc).Ncilia,Cell(nc).eps(v_eps),color_scatter,'MarkerSize',ms,'LineWidth',lw,'MarkerFaceColor',[1,1,1]);hold on;
%        leg{cc}=num2str(nc);
        if Cell(nc).eps(v_eps)<12; EPS_V(cc)=Cell(nc).eps(v_eps)/(Flow_mean(v_eps));
        else EPS_V(cc)=nan;end
        NCILIA_ALL(cc)=Cell(nc).Ncilia;
        cc=cc+1;
        end

    end
end

end
% legend(l,'FontSize',13);
% set(gca,'XScale', 'log', 'YScale', 'log');
%%
L=11*10^-6;   %%%%% lenght cilium
a=(1*10^-6)
eta=0.8*10^-3;%%%%% viscosity water
fd=15;        %%%%% frequency single cilium
Fd=33*10^-12;  %%%%% force acting on cilium in piconewton
Amp=10*10^-6;

A0= 2*(fd/Fd)*4*pi*eta*L*10^-3;
A0=10;
B0= (L/(10^-6))^(1/0.4);

ft_all = fittype('A/(N*(0.5*log(B*N)))','independent','N');
fo_all = fitoptions('Method','NonLinearLeastSquares',...
               'Lower',[A0/10,B0/10],...
               'Upper',[A0*100,B0*10],...
               'StartPoint',[A0,B0]);
fit_out_all = fit(NCILIA_ALL(~isnan(EPS_V))',EPS_V((~isnan(EPS_V)))',ft_all,fo_all);
nbins=5;
[histw,ehistw,vinterval]=hist_nico(NCILIA_ALL,EPS_V,nbins,[1:5:30]);
%histw(histw>9)=nan;  %%%% to remove from statistic the one the went out of range 
%errorbar(vinterval,histw,ehistw,color_median,'MarkerSize',mks,'LineWidth',lw,'MarkerFaceColor',color_face);hold on;
figure(5);
errorbar(vinterval,histw,ehistw*2,color_median,'MarkerSize',mks,'LineWidth',lw,'MarkerFaceColor',color_face);hold on;
%plot(NCILIA,EPS_V,'ko','MarkerSize',mks-3,'LineWidth',lw);hold on;
%errorbar(vinterval,y,ey,'ko','MarkerSize',mks,'LineWidth',lw); hold on;
plot(2:30,fit_out_all(2:30),'r-','MarkerSize',mks,'LineWidth',lw+2);
set(gca,'XScale', 'log', 'YScale', 'log');

int_d_exp= L/(sqrt(fit_out_all.B)*exp(1)^(1/2));
Fd_exp= 2*fd*4*pi*eta*L/(fit_out_all.A*10^3);

A_exp=(log( L/a )-1/2)/fit_out.A;

Fd_theo= (Amp*2*fd)*4*pi*eta*L /(log(L/a)-1/2);

title(strcat('FITALLDATA ciliary distance $C=',num2str(int_d_exp*10^6,1),'\mu m$ and $F_{d}=',num2str(Fd_exp*10^12,2),'pN$'),'Interpreter','latex')
% errorbar(vinterval,y,ey,'ko','MarkerSize',mks,'LineWidth',lw); hold on;
% p_value= polyfit(log(vinterval'),log(y'),1)
% plot(2:30,exp(polyval(p_value,log(2:30))),'r-','MarkerSize',mks,'LineWidth',lw)
% set(gca,'XScale', 'log', 'YScale', 'log');
l{1}='exp';
l{2}='$\frac{A}{\frac{N}{2}\log{(BN)}}$';
 legend(l,'FontSize',13,'Interpreter','latex');
 set(gca,'XScale', 'log', 'YScale', 'log');
xlabel('$N_{c}$ number of cilia','Interpreter','latex', 'FontSize',20);
ylabel('$\epsilon/v_{EX} [mm]^{-1}$ ','Interpreter','latex', 'FontSize',20);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot the 4 lines + fit obtained from the all data fit


set(0,'defaulttextInterpreter','latex') ;

Volts_plot=[2,3,4,5];
lw=0.3;
mks=14;

color_scatter_array={'ko','bo','ro','mo'};
%color_scatter_array={'ko','ko','ko','ko'};
color_median_array={'k--o','r--s','b--<','m--^'};
color_face_array={'k','r','b','m'};
clear l;cleg=1;

HISTW=nan(4,5);
EHISTW=nan(4,5);
 for jj=1:4 
 
Flow_mean= mean(Cal.Flow_int,2);

whichVolt=Volts_plot(jj);
v_eps= (Volts == whichVolt);
l{cleg}= strcat(num2str(Flow_mean(v_eps),'%-5.2f'),' mm/s'); cleg=cleg+1;
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
[histw,ehistw,vinterval]=hist_nico(NCILIA,EPS,nbins,[1:5:30]);
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
fd=15;        %%%%% frequency single cilium
Fd=33*10^-12;  %%%%% force acting on cilium in piconewton
Amp=10*10^-6;

A0= 2*(fd/Fd)*4*pi*eta*L*10^-3;
A0=10;

C0=0.07*10^(-6);
em= exp(1)^(1/2)
%B0= L/(( exp(1)*C0 )^2); 
%B0u= L/(( exp(1)*C0*0.8 )^2);
%B0l= L/(( exp(1)*C0*1.2 )^2);
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

title(strcat('ciliary distance $C=',num2str(int_d_all*10^6,2)...,'\pm',num2str(eint_d_all,1)
    ,'\mu m$ and $F_{d}=',num2str(Fd_all*10^12,2),'\pm',num2str(eFd_all,2),'pN$'),'Interpreter','latex')
l{cleg}='$\frac{A}{\frac{N}{2}\log{(BN)}}$';
legend(l,'FontSize',13,'Interpreter','latex');
set(gca,'XScale', 'log', 'YScale', 'log');
xlabel('$N_{c}$ number of cilia','Interpreter','latex', 'FontSize',20);
ylabel('$\epsilon/v_{EX} [mm]^{-1}$ ','Interpreter','latex', 'FontSize',20);

x0=0;
y0=0;
width=500;
height=400;
set(gcf,'position',[x0,y0,width,height])
set(gca,'FontSize',15);

%% epsilon vs area normalised with velocity, much more noisy than the other one.
 
 
 
Volts_plot=[2,3,4,5];
lw=0.3;
mks=14

color_scatter_array={'ko','bo','ro','mo'};
%color_scatter_array={'ko','ko','ko','ko'};
color_median_array={'k--o','r--s','b--<','m--^'};
color_face_array={'k','r','b','m'};
clear l;cleg=1;
 for jj=1:4 
 
Flow_mean= mean(Cal.Flow_int,2);

figure();
whichVolt=Volts_plot(jj);
v_eps= (Volts == whichVolt);
l{cleg}= strcat(num2str(Flow_mean(v_eps),'%-5.2f'),' mm/s'); cleg=cleg+1
cc=1;
EPS=[];REFF=[];
ms=11;
%lw=0.5;
lw=1;
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
        REFF(cc)=R_eff(nc);
        cc=cc+1;
        end

        if whichVolt>=4 
%        plot(Cell(nc).Ncilia,Cell(nc).eps(v_eps),color_scatter,'MarkerSize',ms,'LineWidth',lw,'MarkerFaceColor',[1,1,1]);hold on;
%        leg{cc}=num2str(nc); 
        EPS(cc)=Cell(nc).eps(v_eps);
        REFF(cc)=R_eff(nc);
        cc=cc+1;
        end

    end
end
%legend(leg)
xlabel('$N_{c}$ number of cilia','Interpreter','latex', 'FontSize',20);
ylabel('$\epsilon [Hz]$ ','Interpreter','latex', 'FontSize',20);
title(strcat('average flow=  ',num2str(Flow_mean(v_eps),'%.1f'),'[mm/s]'),'Interpreter','latex', 'FontSize',15);
saveas(gca,strcat('average flow  ',num2str(Flow_mean(v_eps)),'.jpg'));


nbins=5;
figure(5);
[histw,ehistw,vinterval]=hist_nico(REFF,EPS,nbins)%,[1:5:30]*max(R_eff));
histw(histw>9)=nan;  %%%% to remove from statistic the one the went out of range 
%errorbar(vinterval,histw,ehistw,color_median,'MarkerSize',mks,'LineWidth',lw,'MarkerFaceColor',color_face);hold on;
errorbar(vinterval,histw./(Flow_mean(v_eps)),ehistw./(Flow_mean(v_eps)),color_median,'MarkerSize',mks,'LineWidth',lw,'MarkerFaceColor',color_face);hold on;

xlabel('Number of cilia per cell','FontSize',20);
ylabel('$\epsilon $ [Hz]','FontSize',20);
%ylim([0,20]);
set(gca,'FontSize',16);
%title(strcat('average flow ',num2str(Flow_mean(v_eps)),'[mm/s]'));



end
 legend(l,'FontSize',13);
 set(gca,'XScale', 'log', 'YScale', 'log');

