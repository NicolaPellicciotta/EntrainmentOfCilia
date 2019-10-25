%%%%%%%% this plot the CBF increase with N for wild cells and treated with
%%%%%%%% Cyto_D%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%load Cyto-D %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear all
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

Eper= (ESTD./MED);
med_Eper= nanmedian(Eper)
in=  Eper<(med_Eper+ 1*nanstd(Eper)) & ~isinf(Eper) & ~isnan(Eper) & (NCILIA<20) & (NCILIA~=0) & FREQ>3;% & onlyfirst %./sqrt(NCILIA);


nbins=6;
figure(99);
plot(NCILIA(in),MED(in),'ro','MarkerFaceColor','r');
%colormap([0 0 1;0 1 0; 1 0 0])
%c = colorbar();
%c.Label.String = '\eta_{s}';
NCILIA_cyto =NCILIA(in)
MED_cyto=MED(in);
Eper_cyto=Eper(in);
hold on;

[fit_cyto,gof_cyto,~]= fit(NCILIA(in)', MED(in)','poly1')
%[p_cyto,S]=polyfit(NCILIA(in),MED(in),1);
plot(fit_cyto,'r-')%,'LineWidth',lw)
%y=MED(in);
%r_squared_cyto= 1 - (S.normr/norm(y - mean(y)))^2

%

%%%%%%%%%%%%%%%%% now plor results from wild type CCF%%%%%%%%%%%%%%%%%%

clearvars -except NCILIA_cyto MED_cyto Eper_cyto fit_cyto gof_cyto
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
    
%     FREQ(cc)=FC.Res(nc).f_guess;
%     MED(cc)=FC.Res(nc).medf_debris2;
%     EMED(cc)=FC.Res(nc).medferr;
%     ESTD(cc)=FC.Res(nc).stdf_debris2;
%     ESTD_M(cc)=FC.Res(nc).stdf/sqrt(sum(~isnan(FC.Res(nc).F(:))));
%    leg{cc}=num2str(nc); ;
cc=cc+1;
end


Eper= (ESTD./MED);
in=~isinf(Eper) & ~isnan(Eper) & (NCILIA<25) & (NCILIA~=0)% & onlyfirst %./sqrt(NCILIA);
nbins=6;
figure(99);
plot(NCILIA(in),MED(in),'bd','MarkerFaceColor','b');
%colorbar()
%colormap( [0 0 1; 0.5 0  0.5 ; 1 0 0]);
%c = colorbar();
%c.Label.String = '\eta_{s}';
hold on;


[fit,gof,~]= fit(NCILIA(in)', MED(in)','poly1')
plot(fit,'b-')%,'LineWidth',lw)

%[p,S]=polyfit(NCILIA(in),MED(in),1);
%plot(1:25,polyval(p,1:25),'b-','LineWidth',lw)
%y=MED(in);
%r_squared= 1 - (S.normr/norm(y - mean(y)))^2;


x0=0;
y0=0;
width=500;
height=400;
set(gcf,'position',[x0,y0,width,height])
set(gca,'FontSize',12);
%set(0,'defaulttextInterpreter','latex') ;
legend({ 'CytoChalasin-D','linear fit','control','linear fit'})

e_fit=confint(fit);
ep1= abs(e_fit(2,1)-e_fit(1,1))/2;

e_fit_cyto=confint(fit_cyto);
ep1_cyto= abs(e_fit_cyto(2,1)-e_fit_cyto(1,1))/2;
title({ strcat('m_control =',num2str(fit.p1,2),'$\pm$',num2str(ep1,2),' m_cyto =',num2str(fit_cyto.p1,2),'$\pm$',num2str(ep1_cyto,2))...
    },'Interpreter','latex','FontSize',9)
xlabel([]);
ylabel([]);

fig=figure(99);
saveas(fig,'/u/homes/np451/Dropbox/Synchonisation of mammalian cilia with hydrodynamic forces/PNAS/figures/v3/supp/CBF_cyto_compare/CBF_N_cyto_wild.pdf')

disp(strcat(num2str(fit.p2),' e ',num2str(abs(e_fit(2,2)-e_fit(1,2))/2)))

disp(strcat(num2str(fit_cyto.p2),' e ',num2str(abs(e_fit_cyto(2,2)-e_fit_cyto(1,2))/2))) 































%%
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