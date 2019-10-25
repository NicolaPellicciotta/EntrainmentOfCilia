%%% in this scrip I
%%% do the plot of the Area vs Number of cilia + fit for reviewrs PNAS

%%% I calculate the Area and frequency spatial noise for the cells treated
%%% with flow, 

%%% And area for cells not treated with flow

%%% I calculate the centrin area
%%%% I plot the centrin area vs numer of cilia + fit

%% plot effective radius vs number of cilia
load('/media/np451/Seagate Backup Plus Drive/DATA/synchronisation_paper/all_variables_area_PNAS.mat');
set(0,'defaulttextInterpreter','latex') ;
lw=2;
mks=10;
px2mu=0.0973;

AREA=[]
for nc=1:numel(Cell)  
   AREA(nc)= Cell(nc).Sp.area_debris2*(px2mu^2);   
end

Norm= AREA(NCILIA==1);
AREA=AREA;
R_eff= sqrt(AREA./pi);
Norm= R_eff(NCILIA==1);
R_eff=0.1*R_eff/Norm

%%
%%% making plot with fit from model of referee PNAS model
%%% of nodal cilia 
%%% dynamic area pi*[r+R ]^2 = pi*[ a*(N/pi)^(1/2) + R]^2. 

figure(10)

plot(NCILIA,AREA,'ko','MarkerSize',mks,'LineWidth',lw);
%%%%% fit with polyfit

ft_area = fittype('pi*[a*(x/pi)^(1/2)+R]^2','independent','x');
fo_area = fitoptions('Method','NonLinearLeastSquares',...
               'Lower',[0,0],...
               'Upper',[10,2],...
               'StartPoint',[1,2]);
fit_out_area = fit(NCILIA',AREA',ft_area,fo_area);
efit= confint(fit_out_area); eR= abs(efit(1,1)-efit(2,1));
ea= abs(efit(1,2)-efit(2,2));



hold on;
plot(fit_out_area)%,'r-','MarkerSize',mks,'LineWidth',lw)
title(strcat('fit results R=',num2str(fit_out_area.R,2),'$\pm$',num2str(eR,2),'$\mu$m and a=',num2str(fit_out_area.a,2),'$\pm$',num2str(ea,2),'$\mu$m'),'Interpreter','latex')  
xlabel('N [number of cilia]','Interpreter','latex');
ylabel({'Area [$\mu$m$^{2}$]'},'Interpreter','latex');
legend({'data','$\pi*[a (\frac{N}{\pi})^{1/2}+R]^{2}$ '},'Interpreter','latex');
set(gca,'FontSize',15);
fig=figure(10)
saveas(fig,'/u/homes/np451/Dropbox/Synchonisation of mammalian cilia with hydrodynamic forces/PNAS/figures_review/reviewer2_nodalmodel.pdf')

%% figure with model with amplitude

figure(11)
%%% making plot with fit from model of referee PNAS model
%%% of nodal cilia 
%%% dynamic area pi*[r+R ]^2 = pi*[ a*(N/pi)^(1/2) + R]^2. 
plot(NCILIA,AREA,'ko','MarkerSize',mks,'LineWidth',lw);
%%%%% fit with polyfit

ft_area = fittype(' 2*A*a*(x/pi)^(1/2) + x*a^(2)','independent','x');
fo_area = fitoptions('Method','NonLinearLeastSquares',...
                'Robust','Bisquare',...
               'Lower',[10,0.2],...
               'Upper',[10,2],...
               'StartPoint',[1,2]);
fit_out_area = fit(NCILIA',AREA',ft_area,fo_area);hold on;
efit= confint(fit_out_area); eA= abs(efit(1,1)-efit(2,1));
ea= abs(efit(1,2)-efit(2,2));



plot(fit_out_area)%,'r-','MarkerSize',mks,'LineWidth',lw)
title(strcat('fit results A=',num2str(fit_out_area.A,2),'$\pm$',num2str(eA,2),'$\mu$m and a=',num2str(fit_out_area.a,2),'$\pm$',num2str(ea,2),'$\mu$m'),'Interpreter','latex')  
xlabel('N [number of cilia]','Interpreter','latex');
ylabel({'Area [$\mu$m$^{2}$]'},'Interpreter','latex');
legend({'data','$ 2Aa(\frac{N}{\pi})^{1/2} + Na^{2}$ '},'Interpreter','latex');
x0=0;
y0=0;
width=500;
height=400;
set(gcf,'position',[x0,y0,width,height])
set(gca,'FontSize',15);

set(gca,'FontSize',15);
fig=figure(11)
saveas(fig,'/u/homes/np451/Dropbox/Synchonisation of mammalian cilia with hydrodynamic forces/PNAS/figures_review/reviewer2_amplitude_model_fixedamplitude10um.pdf')


%%
p=polyfit(log(NCILIA),log(AREA),1);hold on;
plot(1:25,exp(polyval(p,log(1:25))),'r-','MarkerSize',mks,'LineWidth',lw)
title('AREA')
xlabel('Ncilia');
ylabel('Area [um2]');
legend('exp','y=Ax^{\alpha}');
%set(gca,'XScale', 'log', 'YScale', 'log');

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










%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%        Here I calculate Area and frequency %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%  
%%%%%%%%%%%%%%%%%%%%%%%%        distribution for cells treated with flow
%%
%nc= 41; %%%% test with a cell with 1 cilium

clear all
cd '/media/np451/Seagate Backup Plus Drive/DATA/synchronisation_paper';
load('allvariables_area.mat');

for nc=25:N_cells
f_guess= Cell(nc).F_rest;
keystr='Drive/';
keystr2='DATA/';

ind_directory= strfind(Cell(nc).Pos,keystr2)+numel(keystr2);
if isempty(ind_directory)
ind_directory= strfind(Cell(nc).Pos,keystr)+numel(keystr);
end


moviedirectory= Cell(nc).Pos(ind_directory:end); 
cd(strcat('/media/np451/Seagate Backup Plus Drive/DATA/synchronisation_paper/',moviedirectory));

%%%% here we calculate the area and frequency of the cell. 
box_size=2;

disp(nc); disp(Cell(nc).Ncilia)
figure();
[Meta] = spatial_noise_fft(f_guess,box_size,[],[]);

%%%% this is to remove the possible strange frequencies 
[F_debris1,BW_debris1] = remove_debris(Meta.F,3,[12,35]);  %%% remove debris frequencies

%%% this to keep only the largest connected figure;
[G,id]=findgroups(BW_debris1(:));
CC = bwconncomp(BW_debris1)
S = regionprops(CC, 'Area');
%Remove small objects:
L = labelmatrix(CC);
Area=[];for fff=1:numel(S); Area(fff)=S(fff).Area;end
max_Area= max(Area);    
    
BW_debris2 = ismember(L, find([S.Area] >= max_Area-1));
F_debris2=F_debris1; F_debris2(~BW_debris2)=nan; 


%%% let s do a plot of the result before the debris
close all;
figure(31);
s = imadjust(Meta.s);
imshow(s)
hold on;

[B,L] = bwboundaries(Meta.BW);
for k=1:length(B);
    k=2;
    boundary = B{k};
plot(boundary(:,2), boundary(:,1), 'r', 'LineWidth', 2)
end
truesize
fig=figure(31);saveas(fig,'dynamic_area.pdf')

%red_image= s; red_image(Meta.BW)=1;
%s_rgb = cat(3,red_image, s,s);
%paintcolor = [1 0 0];   %red

%imagesc(s_rgb);
%rectangle('Position',Meta.rect_out, 'EdgeColor', paintcolor)

figure(2);
F_plot= F_debris2; F_plot(~isnan(F_plot))=0;
imagesc(F_plot);


%%%% this is to save all in the Cell variable
Cell(nc).Sp.F=Meta.F;
Cell(nc).Sp.F_debris2=F_debris2;
Cell(nc).Sp.BW_debris2=BW_debris2;
Cell(nc).Sp.BW=Meta.BW;

Cell(nc).Sp.box_size=box_size;
Cell(nc).Sp.s_roi=Meta.s_roi
Cell(nc).Sp.area= sum(~isnan(Meta.F(:)))*box_size^2;
Cell(nc).Sp.area_debris2= sum(~isnan(F_debris2(:)))*box_size^2;

Cell(nc).Sp.medf=nanmedian(Meta.F(:));
Cell(nc).Sp.medf_debris2=nanmedian(F_debris2(:));

Cell(nc).Sp.medferr= nanmedian(abs(Meta.F(:)-nanmedian(Meta.F(:))));
Cell(nc).Sp.stdf=nanstd(Meta.F(:));
Cell(nc).Sp.stdf_debris2=nanstd(F_debris2(:));

Cell(nc).Sp.rect_out=Meta.rect_out;




end

%save('all_variables_area_PNAS.mat','Cell')

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%        Here I calculate Area and frequency %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%  
%%%%%%%%%%%%%%%%%%%%%%%%        distribution for cells not treated with flow
%%
clear all
close all
data_dir='/media/np451/Seagate Backup Plus Drive/DATA/synchronisation_paper/14.12.18/frequency';
cd(data_dir);
%d=dir('60X_2*movie');
load('FC.mat')

N_cells=numel(FC.Res)
%%
for nc=34:N_cells
   
%%%% here we calculate the area and frequency of the cell. 
box_size=2;

disp(nc); disp(FC.Res(nc).Ncilia)
figure();
[Meta] = spatial_noise_fft(FC.Res(nc).f_guess,box_size,FC.Res(nc).filename,[]);

%%%% this is to remove the possible strange frequencies 
[F_debris1,BW_debris1] = remove_debris(Meta.F,3,[8,35]);  %%% remove debris frequencies

%%% this to keep only the largest connected figure;
[G,id]=findgroups(BW_debris1(:));
CC = bwconncomp(BW_debris1)
S = regionprops(CC, 'Area');
%Remove small objects:
L = labelmatrix(CC);
Area=[];for fff=1:numel(S); Area(fff)=S(fff).Area;end
max_Area= max(Area);    
    
BW_debris2 = ismember(L, find([S.Area] >= max_Area-1));
F_debris2=F_debris1; F_debris2(~BW_debris2)=nan; 


%%% let s do a plot of the result before the debris
close all;
figure(1);
s = imadjust(Meta.s);
red_image= s; red_image(Meta.BW)=1;
s_rgb = cat(3,red_image, s,s);
paintcolor = [1 0 0];   %red

imagesc(s_rgb);
rectangle('Position',Meta.rect_out, 'EdgeColor', paintcolor)

figure(2);
F_plot= F_debris2; F_plot(~isnan(F_plot))=0;
imagesc(F_plot); axis equal


%%%% this is to save all in the Cell variable
FC.Res(nc).F=Meta.F;
FC.Res(nc).F_debris2=F_debris2;
FC.Res(nc).BW_debris2=BW_debris2;
FC.Res(nc).BW=Meta.BW;

FC.Res(nc).box_size=box_size;
FC.Res(nc).s_roi=Meta.s_roi
FC.Res(nc).area= sum(~isnan(Meta.F(:)))*box_size^2;
FC.Res(nc).area_debris2= sum(~isnan(F_debris2(:)))*box_size^2;

FC.Res(nc).medf=nanmedian(Meta.F(:));
FC.Res(nc).medf_debris2=nanmedian(F_debris2(:));

FC.Res(nc).medferr= nanmedian(abs(Meta.F(:)-nanmedian(Meta.F(:))));
FC.Res(nc).stdf=nanstd(Meta.F(:));
FC.Res(nc).stdf_debris2=nanstd(F_debris2(:));

FC.Res(nc).rect_out=Meta.rect_out;    
    
end

 save('FC_PNAS.mat','FC');
 
 %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%        Here I calculate Cenrin Area for the one %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%  
%%%%%%%%%%%%%%%%%%%%%%%%        that they have it

%%
%nc= 41; %%%% test with a cell with 1 cilium

clear all
cd '/media/np451/Seagate Backup Plus Drive/DATA/synchronisation_paper';
load('all_variables_area_PNAS.mat');
N_cells=numel(Cell)

%%
for nc=25:N_cells
f_guess= Cell(nc).F_rest;
keystr='Drive/';
keystr2='DATA/';

ind_directory= strfind(Cell(nc).Pos,keystr2)+numel(keystr2);
if isempty(ind_directory);
ind_directory= strfind(Cell(nc).Pos,keystr)+numel(keystr);
end


moviedirectory= Cell(nc).Pos(ind_directory:end); 
cd(strcat('/media/np451/Seagate Backup Plus Drive/DATA/synchronisation_paper/',moviedirectory));

cen_file= dir('*entrin*movie');
if ~isempty(cen_file)
    disp(strcat('nc=',num2str(nc)));
    disp(strcat('cilia =',num2str(Cell(nc).Ncilia)));
    
box_size=2;
figure();
[Meta] = spatial_noise_fft(f_guess,box_size,[],[]);    
Cell(nc).Centrin= Centrin_fun();
close all;

%%%%% plot centrin image
figure(21);
fm=(mat2gray(Cell(nc).Centrin.fm));
imshow(fm);hold on;
[B,L] = bwboundaries(Cell(nc).Centrin.mask,'noholes');
 boundary = B{1};
   plot(boundary(:,2), boundary(:,1), 'r', 'LineWidth', 2)
%red_image= fm; red_image(Cell(nc).Centrin.mask)=1;
%s_rgb_cen = cat(3,red_image, fm,fm);
%imagesc(s_rgb_cen);
truesize

fig=figure(21);saveas(fig,'centrin_threeshold.pdf')

%%%%% plot centrin on standard deviation
figure(22)
s = imadjust(Meta.s);
red_image= s; red_image(Cell(nc).Centrin.mask)=1;
s_rgb = cat(3,red_image, s,s);
imagesc(s_rgb);
truesize
title(strcat('Number of cilia= ',num2str(Cell(nc).Ncilia), ' nc=',num2str(nc)));
hold on;
fig=figure(22);saveas(fig,'centrin_std.pdf')

else
Cell(nc).Centrin= [];
end

end

%%save('all_variables_area_PNAS_centr.mat','Cell');

%% plot centrin area as a function of number of cilia

load('/media/np451/Seagate Backup Plus Drive/DATA/synchronisation_paper/all_variables_area_PNAS_centr.mat');
set(0,'defaulttextInterpreter','latex') ;
lw=2;
mks=10;
px2mu=0.0973;

AREA=[];
CEN=[]
for nc=1:numel(Cell) 
    
    if ~isempty(Cell(nc).Centrin) & nc~=36 & nc~=41 
        CEN(nc)= sum(Cell(nc).Centrin.BW(:))*(px2mu^2);
    else
        CEN(nc) = nan;
    end
    AREA(nc)= Cell(nc).Sp.area_debris2*(px2mu^2);
    NCILIA(nc)=Cell(nc).Ncilia;
end


figure(10)

plot(NCILIA,CEN,'ko','MarkerSize',mks,'LineWidth',lw);
hold on;
%%
%%%%% fit with polyfit
fit_n= NCILIA(~isnan(CEN))
fit_c= CEN(~isnan(CEN))

[p,s]=polyfit(fit_n,fit_c,1);
ste = sqrt(diag(inv(s.R)*inv(s.R')).*s.normr.^2./s.df);

plot(polyval(p,1:25),'r','LineWidth',lw);


title(strcat('fit results a$^2$=',num2str(p(1),2),'$\pm$',num2str(ste(1),2),'$\mu$m$^2$ and C=',num2str(p(2),2),'$\pm$',num2str(ste(2),2),'$\mu$m$^2$'),'Interpreter','latex')  
xlabel('N [number of cilia]','Interpreter','latex');
ylabel({'Area [$\mu$m$^{2}$]'},'Interpreter','latex');
legend({'data','y= Na$^2+C$ '},'Interpreter','latex');
set(gca,'FontSize',15);
x0=0;
y0=0;
width=500;
height=400;
set(gcf,'position',[x0,y0,width,height])
set(gca,'FontSize',15);
fig=figure(10)
saveas(fig,'/u/homes/np451/Dropbox/Synchonisation of mammalian cilia with hydrodynamic forces/PNAS/figures_review/reviewer2_centrin_area_ncilia.pdf')


%%
%%%%% fit with fit matlab forcing C=0
fit_n= NCILIA(~isnan(CEN))
fit_c= CEN(~isnan(CEN))

ft_cen = fittype('a2*x','independent','x');
fo_cen = fitoptions('Method','NonLinearLeastSquares',...
               'Lower',[0],...
               'Upper',[2],...
               'StartPoint',[0.8],'Robust','Bisquare');
[fit_out_cen,gof,output] = fit(fit_n',fit_c',ft_cen,fo_cen);
efit= confint(fit_out_cen); ea2= abs(efit(1)-efit(2));
%ea= abs(efit(1,2)-efit(2,2));

figure(111)
plot(fit_n,fit_c,'ko','MarkerSize',mks,'LineWidth',lw);
hold on;
plot(fit_out_cen);


%title(strcat('fit results a$^2$=',num2str(fit_out_cen.a2,2),'$\pm$',num2str(ea2,2),'$\mu$m$^2$'),'Interpreter','latex')  
%xlabel('N [number of cilia]','Interpreter','latex');
%ylabel({'Area [$\mu$m$^{2}$]'},'Interpreter','latex');
legend({'data','linear fit '},'Interpreter','latex');
x0=0;
y0=0;
width=500;
height=400;
set(gcf,'position',[x0,y0,width,height])
set(gca,'FontSize',15);
xlabel('');ylabel('')
fig=figure(111)
%saveas(fig,'/u/homes/np451/Dropbox/Synchonisation of mammalian cilia with hydrodynamic forces/PNAS/figures_review/reviewer2_centrin_area_ncilia_0intercept.pdf')



