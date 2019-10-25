clear all
cd '/media/np451/Seagate Backup Plus Drive/DATA/synchronisation_paper';
load('allvariables_area.mat');
%%% Number of cilia and Frequency with cells used with external flow
cc=1;FREQ=[];NCILIA=[];

lw=1;
mks=14

for nc=1:Ncells;
    %plot(Cell(nc).Ncilia, Cell(nc).F_rest,'bo','LineWidth',lw,'MarkerSize',mks);hold on;
    leg{cc}=Cell(nc).Pos(end-2:end); 
    FREQ(cc)=Cell(nc).F_rest;
    MED(cc)=Cell(nc).Sp.medf;
    NCILIA(cc)= Cell(nc).Ncilia;
    EMED(cc)=Cell(nc).Sp.medferr;
    ESTD(cc)=Cell(nc).Sp.stdf;
    ESTD_M(cc)=Cell(nc).Sp.stdf/sqrt(sum(~isnan(Cell(nc).Sp.F(:))));
    leg{cc}=num2str(nc); ;cc=cc+1;
end

%% add this info to the other experiments with cells not use for flow experiments
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


%%

Eper= (100*ESTD./FREQ);
%onlyfirst= ones([1,numel(NCILIA)]);
%onlyfirst(57:end)=0;
in=~isinf(Eper) & ~isnan(Eper) & (NCILIA<25) & (NCILIA~=0)% & onlyfirst %./sqrt(NCILIA);
%Eper= (100*EMED./FREQ);
ind1= in & Eper<mean(Eper(in)) ;
ind2= in &  Eper>mean(Eper(in));
nbins=6;
figure();
scatter(NCILIA(in),MED(in),20,Eper(in),'filled');
colormap( [0 0 1; 0.5 0  0.5 ; 1 0 0]);
c = colorbar();
c.Label.String = '(\sigma f)/ f';




[histw,ehistw,vinterval]=hist_nico(NCILIA(ind1),MED(ind1),nbins);hold on;
errorbar(vinterval,histw,ehistw,'--ko','LineWidth',lw,'MarkerSize',mks,'MarkerFaceColor','b');
xlabel('Number of cilia per cell','FontSize',20);
ylabel('CBF [Hz]','FontSize',20);
ylim([10,32])

[histw,ehistw,vinterval2]=hist_nico(NCILIA(ind2),MED(ind2),nbins,vinterval);
errorbar(vinterval(1:end-1),histw(1:end-1),ehistw(1:end-1),'--ro','LineWidth',lw,'MarkerSize',mks,'MarkerFaceColor','r');
xlabel('Number of cilia per cell','FontSize',20);
ylabel('CBF [Hz]','FontSize',20);
ylim([10,32])

legend('alldata','low noise', 'high noise')



% [histw,ehistw,vinterval]=hist_nico(NCILIA(in),MED(in),nbins);
% errorbar(vinterval,histw,ehistw,'--ko','LineWidth',lw,'MarkerSize',mks,'MarkerFaceColor','k');
% xlabel('Number of cilia per cell','FontSize',20);
% ylabel('CBF [Hz]','FontSize',20);
% ylim([10,32])
%% box plot
per= (100*ESTD./FREQ);
%onlyfirst= ones([1,numel(NCILIA)]);
%onlyfirst(57:end)=0;
in=~isinf(Eper) & ~isnan(Eper) & (NCILIA<25) & (NCILIA~=0)% & onlyfirst %./sqrt(NCILIA);
%Eper= (100*EMED./FREQ);
ind1= in & Eper<=nanmean(Eper(in)) ;
ind2= in &  Eper>nanmean(Eper(in));
nbins=6;
figure();

Box1= hist_box(NCILIA(ind2),MED(ind2),nbins);
boxplot(Box1.values(1:end-1),Box1.g(1:end-1),'Labels',{'1-4','5-9','10-14','15-20'},'Colors','r');

Box1= hist_box(NCILIA(ind1),MED(ind1),nbins);hold on;
boxplot(Box1.values,Box1.g,'Labels',{'1-4','5-9','10-14','15-20','20-25'},'Colors','k'); hold on;



xlabel('Number of cilia per cell','FontSize',20);
ylabel('CBF [Hz]','FontSize',20);
ylim([10,32])
legend(['low noise', 'high noise'])




%set(gca,'XScale', 'log', 'YScale', 'log');
%%


p_cbf= polyfit(log(NCILIA(NCILIA<23)),log(FREQ(NCILIA<23)),1);
plot(2:30,exp(polyval(p_cbf,log(2:30))),'r-','MarkerSize',mks,'LineWidth',lw)
set(gca,'XScale', 'log', 'YScale', 'log');
legend({'single cell','binnig average','$y=Ax^{\alpha}$'},'Interpreter','latex')
title(strcat('CBF vs cilia $\alpha= ',num2str(p_cbf(1),2),'$'),'Interpreter','latex');



% %%%% fitting with model
% 
% L=11*10^-6;   %%%%% lenght cilium
% a=(1*10^-6)
% eta=0.8*10^-3;%%%%% viscosity water
% fd=15;%%%%% frequency single cilium
% Amp=10*10^-6;
% Fd=10^-11;
% 
% A0= Fd/(Amp*4*pi*eta*L);
% B0= (L/(2*10^-4))^2;
% B0u= (L/(2*10^-5))^2;
% B0l= (L/(2*10^-3))^2;
% 
% ft_cbf = fittype('A*(N*(0.5*log(B*N)))','independent','N');
% fo_cbf = fitoptions('Method','NonLinearLeastSquares',...
%                'Lower',[A0/2,B0],...
%                'Upper',[A0*2, B0],...
%                'StartPoint',[A0,B0],'Robust','Bisquare');
% fit_out_cbf = fit(NCILIA',FREQ',ft_cbf,fo_cbf);
% 
% efit_out_cbf=confint(fit_out_cbf);
% eA=abs(efit_out_cbf(1,1)-efit_out_cbf(2,1));
% eB=abs(efit_out_cbf(1,2)-efit_out_cbf(2,2));
% if isnan(eB); eB=abs(B0u-B0l); end
% 
% plot(2:30,fit_out_cbf(2:30),'r-','MarkerSize',mks,'LineWidth',lw+2);
% set(gca,'XScale', 'log', 'YScale', 'log');
% 
% 
% 
% int_d_cbf= L/(sqrt(fit_out_cbf.B)*exp(1)^(1/2));
% eint_d_cbf= (10^6)*(L/(2*exp(1)^(1/2)))*fit_out_cbf.B^(-3/2)* eB ;
% 
% Fd_cbf= fit_out_cbf.A*Amp*4*pi*L*eta;
% eFd_cbf = Amp*4*pi*L*eta*eA;  %%% error in pN
% 
% Fd_theo= (Amp*2*fd)*4*pi*eta*L /(log(L/a)-1/2);
% 
% title(strcat('ciliary distance $C=',num2str(int_d_cbf*10^6,2),'\pm',num2str(eint_d_cbf,1),'\mu m$ and $F_{d}=',num2str(Fd_cbf*10^12,2),'\pm',num2str(eFd_cbf,2),'pN$'),'Interpreter','latex')
% l{cleg}='$\frac{A}{\frac{N}{2}\log{(BN)}}$';
% legend(l,'FontSize',13,'Interpreter','latex');
% set(gca,'XScale', 'log', 'YScale', 'log');
% xlabel('$N_{c}$ number of cilia','Interpreter','latex', 'FontSize',20);
% ylabel('$CBF [mm]^{-1}$ ','Interpreter','latex', 'FontSize',20);

%% this is a part to find out about metachronicity in the cell
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
metarray=[200,20,20]
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
if Cell(nc).Rc{rc}.Period~=0 & Cell(nc).Rc{rc}.syn==1 & Cell(nc).Rc{rc}.Volt==5;
    ind(cc)= rc;cc=cc+1;
end
end

movie=strcat('V5*P',num2str(Cell(nc).Rc{ind(2)}.Period));
[Meta_syn] = spatial_noise_fft_meta(1000/Cell(nc).Rc{ind(2)}.Period,box_size,movie,Cell(nc).Sp.rect_out,metarray);
std_p_syn=abs(nanmean( exp(i*Meta_syn.P_meanf(:)) ));
std_p=abs(nanmean( exp(i*Meta.P_meanf(:)) ))
std_f_syn=nanstd(Meta_syn.F(:));
std_f=nanstd(Meta.F(:));


%% for Meta at rest
figure(1)
pht=Meta.P_meanf;
Ntheta=6;NR=10;
%[ r_bin,th_bin,fph ] = phase_fft( pht-pht(:,:,1),Ntheta,NR );
%%%% spatial phase correlation
[ r_bin,th_bin,cf ] = phase_correlation_function_2d(pht-pht(:,:,1),Ntheta,NR );
%%%% time correlation of the spatial phase difference
clear cf_t
for tau=1:size(cf,3);
 cf_t(:,:,tau)= angle(nanmean(exp(i*cf(:,:,1:end-tau+1)).* conj(exp(i*cf(:,:,tau:end))),3));
end

%subplot(1,2,1);
figure(2)
px2mu=0.14*40/60
for th=1:size(cf,2)
    plot(r_bin*4*px2mu,unwrap(cf(:,th,10)),'-o'); hold on;
    l{th}=num2str(th_bin(th));
end
legend(l)
xlabel('distance [um]');
ylabel('phase difference [rad]');
%subplot(1,2,2);
%imagesc((pht(:,:,1)))
%axis equal
%% for Meta_syn
pht=Meta_syn.P_meanf;
Ntheta=6;NR=10;
%[ r_bin,th_bin,fph ] = phase_fft( pht-pht(:,:,1),Ntheta,NR );
%%%% spatial phase correlation
[ r_bin,th_bin,cf_syn ] = phase_correlation_function_2d(pht-pht(:,:,1),Ntheta,NR );
%%%% time correlation of the spatial phase difference

%subplot(1,2,1);
figure()
px2mu=0.14*40/60
for th=1:size(cf_syn,2)
    plot(r_bin*4*px2mu,(cf_syn(:,th,10)),'-o'); hold on;
    l{th}=num2str(th_bin(th));
end
legend(l)
xlabel('distance [um]');
ylabel('phase difference [rad]');
%subplot(1,2,2);
%imagesc((pht(:,:,1)))
%axis equal
%%

end
%% nuovo approccio per trovare metacronal wave based on black and white.
N_cells=numel(Cell)
%for nc=1:N_cells
nc= 55; %%%% my favourite cell
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
metarray=[200,20,20]
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
if Cell(nc).Rc{rc}.Period~=0 & Cell(nc).Rc{rc}.syn==1 & Cell(nc).Rc{rc}.Volt==5;
    ind(cc)= rc;cc=cc+1;
end

end
%end

