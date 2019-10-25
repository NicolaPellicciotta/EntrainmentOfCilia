%%%% script DDM for frequency (on boxes) %%%
%%% to use this script you should load only Documents/MATLAB/DDM_cilia and Documents/MATLAB/DDM/MatLab_Common_Function 
%% 
path = '/home/np451/Desktop/Mouse/synchronisation/23.5.18/normal/';
cd(path);
 mkdir('analysis_multiDDM'); 
ana_dir= strcat(path,'analysis_multiDDM/');
%store_dir= '/run/user/10704/gvfs/smb-share:server=sf3.bss.phy.private.cam.ac.uk,share=space/np451/6.4.18/beads/';
store_dir= '/home/np451/Desktop/Mouse/synchronisation/23.5.18/normal/data'
cd(store_dir);
d=dir('*.movie');
box_size=[32];  %%% box size for ddm


for nf=1:size(d,1)
    filename= d(nf).name;
    cd(ana_dir);
    if exist(filename)==0;    
    cd(store_dir);    %%% moving to store to load frames

    cilia = DDM_Analysis_nico(filename);

    cd(ana_dir); mkdir(filename); cd(filename);  %%%% moving to the analysis dit

   
    cilia.VariableBoxSize_Analysis(box_size);
    save([cilia.Filename(1:end-5),'mat'],'cilia');  %%% save sruslts

    % plot good boxes with std %
    end
end
%% extract frequency and plots from DDM
%ana_dir = '/home/np451/Desktop/Mouse/MAR2018/analysis_multiDDM/'
%cd(ana_dir);

box_size=64;

analysis_folder = {'/run/user/10704/gvfs/smb-share:server=sf3.bss.phy.private.cam.ac.uk,share=space/np451/ependymalJune/DDM_sigmoids/Analysis_18.6_T/',...
    '/run/user/10704/gvfs/smb-share:server=sf3.bss.phy.private.cam.ac.uk,share=space/np451/ependymalJune/DDM_sigmoids/Analysis_18.6_NT',...
    '/run/user/10704/gvfs/smb-share:server=sf3.bss.phy.private.cam.ac.uk,share=space/np451/ependymalJune/DDM_sigmoids/Analysis_13.6_MT',...
    '/run/user/10704/gvfs/smb-share:server=sf3.bss.phy.private.cam.ac.uk,share=space/np451/ependymalJune/DDM_sigmoids/Analysis_18.6_MT',...
    '/run/user/10704/gvfs/smb-share:server=sf3.bss.phy.private.cam.ac.uk,share=space/np451/ependymalJune/DDM_sigmoids/Analysis_13.6_T',...
    '/run/user/10704/gvfs/smb-share:server=sf3.bss.phy.private.cam.ac.uk,share=space/np451/ependymalJune/DDM_sigmoids/Analysis_13.6_NT'};

for dd=1:numel(analysis_folder)
 cd(analysis_folder{dd})
 d=dir('20X*.mat');
 clear CF; clear NF
for nf=1:size(d,1)
%    filename= d(nf).name;
%    cd(ana_dir);cd(filename);
%    load(strcat(filename(1:end-6),'.mat'));
    load(d(nf).name)
    
    boxes=[];for bb= 1:numel(cilia.Results); boxes(bb)= cilia.Results(bb).BoxSize;
    end
    bsi = find(boxes == box_size);
        
    new_std = cilia.std_fs;
    f_y=cilia.width;
    f_x=cilia.height;
%    new_std=   new_std(1:floor(f_x/box_size)*box_size , 1:floor(f_y/box_size)*box_size ); 
%    new_mask = oversample(cilia.Results(bsi).ind_good_boxes,box_size);
%    mask_overlay(new_std, new_mask, [1 0 0], 0.3);
%    fig = get(groot,'CurrentFigure');
%    saveas(fig,[cilia.Filename(1:end-5),'jpg']);
%    close 'all'

    % plot frequencies
    cilia.gather_results;
    
    %% correlation function from SAVA results
    box_size_SAVA= 4;
    freq=cilia.SAVAlike.frequency_map;
    [ind_good_bins]= SAVA_good_boxes_nico(cilia,box_size_SAVA,20000);
    f_filter= (freq>10 & freq<40) & ind_good_bins ;   
    [r_bin,N_bin,cf,f_res]= FreqCorrFunc_fft(freq,f_filter(:),box_size_SAVA);
    

    
    %% correlation with DDM_boxes
%    histogram(cilia.Results(bsi).MedianFrequencyVec); %% plot histogram 
%    fig = get(groot,'CurrentFigure');
%    saveas(fig,[cilia.Filename(1:end-6),'_freq_.jpg']);
   

%    freq=cilia.Results(bsi).MedianFrequencyVec(:);
%    f_filter= freq(:)>5 & freq(:)<35;
%    [r_bin,cf,N_bin,f_res]= FreqCorrFunc_DDM(cilia,freq(:),f_filter,box_size);
    
%    plot(r_bin,cf,'o');
    CF(nf,:)=cf;
    NF(nf,:)=N_bin;
    F_filter(nf,:,:)= f_filter;
%    FREQ(nf,:)=freq(:);
    
    %save('freq.mat','freq');
end
save('correlation_fun_SAVA.mat','CF','r_bin','NF','F_filter');

end
%% check for metachronalwave

d=dir('64*.mat');
N_exp= numel(dir);
box_size=1024;

for nf=1:size(d,1)
    load(d(nf).name)
    
    boxes=[];for bb= 1:numel(cilia.Results); boxes(bb)= cilia.Results(bb).BoxSize;
    end
    bsi = find(boxes == box_size);
    cilia.gather_results;
    A(nf,:)=mean(cilia.Results(1).Box.Iqtau,2)
   % plot(A);hold on;
    
    
end

%% compare all correaltion functions

analysis_folder = {'/run/user/10704/gvfs/smb-share:server=sf3.bss.phy.private.cam.ac.uk,share=space/np451/ependymalJune/DDM_sigmoids/Analysis_18.6_NT',...
    '/run/user/10704/gvfs/smb-share:server=sf3.bss.phy.private.cam.ac.uk,share=space/np451/ependymalJune/DDM_sigmoids/Analysis_13.6_NT',...   
    '/run/user/10704/gvfs/smb-share:server=sf3.bss.phy.private.cam.ac.uk,share=space/np451/ependymalJune/DDM_sigmoids/Analysis_13.6_MT',...
    '/run/user/10704/gvfs/smb-share:server=sf3.bss.phy.private.cam.ac.uk,share=space/np451/ependymalJune/DDM_sigmoids/Analysis_18.6_MT',...
    '/run/user/10704/gvfs/smb-share:server=sf3.bss.phy.private.cam.ac.uk,share=space/np451/ependymalJune/DDM_sigmoids/Analysis_13.6_T',...
    '/run/user/10704/gvfs/smb-share:server=sf3.bss.phy.private.cam.ac.uk,share=space/np451/ependymalJune/DDM_sigmoids/Analysis_18.6_T/'};
px2mu= 0.26;
lambda=[];
for dd=1:numel(analysis_folder)
 cd(analysis_folder{dd})
 load('correlation_fun_SAVA.mat');
 hold on;
 sp= (1:size(CF,2))*box_size*px2mu;
 cc=sum((CF.*NF),1)./sum(NF,1);  %%%% weighed average with the number of data
 plot(sp,cc,'-o');hold on;
 
 %%fit with exp
ft = fittype('a*exp(-x/b)+c','independent','x');
fo = fitoptions('Method','NonlinearLeastSquares',...
               'Lower',[0,1,-1],...
               'Upper',[1,Inf,1],...
               'StartPoint',[0.9, box_size*px2mu*5, 0]);
f_res = fit(sp',cc',ft,fo);
lambda(dd)= f_res.b; 
 
cd('..')
end

xlabel('distance [um]');ylabel('correlation');xlim([5,100]);
%legend(le)
legend({'18.6_NT','13.6_NT','13.6_MT','18.6_MT','13.6_T','18.6_T'})
set(gca,'XScale', 'log', 'YScale', 'log');

figure(2);
plot(mean(reshape(lambda,[2,3]),1));


%% Load Data

ana_dir = '/home/np451/Desktop/Mouse/MAR2018/analysis_multiDDM/'
cd(ana_dir);
d=dir('*.movie');
N_exp= numel(dir);


for nf=1:size(d,1)
        filename= d(nf).name;
        cd(ana_dir);cd(filename);
        s = str2num(filename(end-7:end-6));
        m = str2num(filename(end-10:end-9));
        h = str2num(filename(end-13:end-12));
        t(nf)= h*3600 + m*60 + s; 
        
        if isempty(strfind(filename,'Pa'));
            shear(nf)=0;
        else
            shear(nf)=str2num(filename(strfind(filename,'Pa')-1)); 
        end
            
        
%        load(strcat(filename(1:end-6),'.mat')); cilia.gather_results;
%        freq=cilia.Results.MedianFrequencyVec(:);
%        save('freq.mat','freq');
%        F(jj,nf)=nanmean( cilia.Results.MedianFrequencyVec(:));        
        load('freq.mat');
        freq=freq(freq>4 & freq<16);
        F(nf)=nanmean( freq(:));
        N(nf)= numel(~isnan(freq(:)));
        std_F(nf)= nanstd(freq(:))/sqrt(N(nf));
        c{nf}=freq(:);
end    
shear= 10.^(-shear);
shear(shear==1)=0;
[t,I] =sort(t);
t=(t-min(t(:)))/60;
shear=shear(I);
F=F(I);
%%
plot(t,F,'o-');
hold on;
plot(t,shear*1000*max(F(:)));
ylim([min(F(:)),max(F(:))]);
%%

cd(path_dir);

%%% mean without weights with error:  F_t ; eF_t ; rv  (CBF, errCBF, viscosity)
F_t=[mean(mean(F(mu==0,:),1),2),mean(mean(F(mu==0.5,:),1),2),mean(mean(F(mu==1,:),1),2),mean(mean(F(mu==1.5,:),1),2),mean(mean(F(mu==2,:),1),2)];
std_t = [std(std(F(mu==0,:),[],1),[],2),std(std(F(mu==0.5,:),[],1),[],2),std(std(F(mu==1,:),[],1),[],2),std(std(F(mu==1.5,:),[],1),[],2),std(std(F(mu==2,:),[],1),[],2)];
std_N = [numel(F(mu==0,:)),numel(F(mu==0.5,:)),numel(F(mu==1,:)),numel(F(mu==1.5,:)),numel(F(mu==2,:))];
eF_t = std_t./sqrt(std_N);    
load('MC_mu.mat')
rv= M(2,:); %%% real viscosity Pa*s


%% mean with weights with error:  F_t ; eF_t ; rv  (CBF, errCBF, viscosity)
tCBF= zeros(size(mu));
MC=[0,0.5,1,1.5,2];
for k=1:numel(mu); tCBF(k)=sum( (F(k,:).*(std_F(k,:).^2))/(sum(std_F(k,:).^2)) ); teCBF(k)=sum( ((std_F(k,:).^2).*(std_F(k,:).^2))/(sqrt(N_exp)*(sum(std_F(k,:).^2))) ); end;

MC=[0,0.5,1,1.5,2];
clear CBF;clear eCBF;
for k=1:numel(MC); temp_CBF=tCBF(mu==MC(k)); temp_eCBF= teCBF(mu==MC(k))'; CBF(k)=[sum( temp_CBF.*temp_eCBF)/(sum(temp_eCBF))]; eCBF(k)=[std(temp_CBF)/sqrt(numel(temp_CBF))];end; %eCBF(k)=[sum( temp_eCBF.*temp_eCBF)/(sum(temp_eCBF))];end;

figure(); errorbar(rv,CBF,eCBF,'-o'); hold on;
title('CBF vs viscosity weighted average'); legend('exp point');
xlabel('{$\eta[Pa*s]$}','Interpreter','latex','FontSize',15);
ylabel('{$CBF [Hz]$}','Interpreter','latex','FontSize',15);
%saveas(gcf,'CBF_viscosity_lin_weighted','pdf');



%% linear 

figure(); plot(mu,F,'o');
title('CBF vs viscosity'); xlabel('viscosity [Pa*s]');ylabel('CBF [Hz]');
saveas(gcf,'CBF_allfield_viscosity.png');


figure(); plot(rv,F_t,'o');
title('CBF vs viscosity'); xlabel('viscosity [Pa*s]');ylabel('CBF [Hz]');


figure(); errorbar(rv,F_t,eF_t,'o');
title('CBF vs viscosity'); xlabel('viscosity [Pa*s]');ylabel('CBF [Hz]');

%% loglog


figure(); errorbar(rv,F_t,eF_t,'o'); hold on;
title('CBF vs viscosity'); xlabel('viscosity [Pa*s]');ylabel('CBF [Hz]');
saveas(gcf,'CBF_viscosity_lin','jpg');
p= polyfit(log(rv),log(F_t),1);
fit_f=polyval(p, log(rv));

plot(rv,exp(fit_f),'-');

legend('exp',['y=x^{',num2str(p(1),2),'}']);
set(gca,'xscale','linear','yscale','linear');
saveas(gcf,'CBF_viscosity_lin_fit','jpg');

set(gca,'xscale','log','yscale','log');
saveas(gcf,'CBF_viscosity_log_fit','jpg');



%% semilogx
figure(); errorbar(rv,F_t,eF_t,'d'); hold on;

%%%%% fit semilogx CBF = log(visc)

px= polyfit(log(rv),(F_t),1);
fit_fx=polyval(px, log(rv));
plot(rv,(fit_fx),'-');

set(gca,'xscale','log','yscale','lin');
legend('exp',['y=',num2str(px(1),2),'*log(x)']);
title('CBF vs viscosity'); xlabel('viscosity [Pa*s]');ylabel('CBF [Hz]');
saveas(gcf,'CBF_viscosity_logx','jpg');


%% correlate number of cells with frequency

for ii=1:size(c,1)
for i=1:size(c,2);
    plot(numel(c{ii,i}), nanmean(c{ii,i}),'o');hold on;
end;
end

%%
cc=zeros(size(c,1));
for ii=1:size(c,1)
for i=1:size(c,2);
    cc(ii)= cc(ii)+numel(c{ii,i}); 
end;
end

cc=cc(:,1);cc=cc(:);

cc=reshape(cc,[numel(cc)/2,2]);
c_m= mean(cc,2);
ec_m= std(cc,[],2)/2;















%% previuos
path_dir = '/home/np451/Desktop/ependymal data/6.6/NOFLOW/'
cd(path_dir);
mkdir('analysis');
d=dir('*.movie');
R=[]

for nf=10:20%size(d,1);
    
    filename= d(nf).name;
    cd(path_dir);
    cilia=DDM_Analysis(filename);
    cd('analysis');

    cilia.VariableBoxSize_Analysis([box_size]);
    save([cilia.Filename(1:end-5),'mat'],'cilia');  %%% save sruslts



    % plot good boxes with std %


    new_std = cilia.std_fs;
    new_std=   new_std(1:floor(f_x/box_size)*box_size , 1:floor(f_y/box_size)*box_size ) 
    new_mask = oversample(cilia.Results.ind_good_boxes,box_size);
    mask_overlay(new_std, new_mask, [1 0 0], 0.3);
    fig = get(groot,'CurrentFigure');
    saveas(fig,[cilia.Filename(1:end-5),'jpg']);
    close 'all'

    % plot frequencies
    cilia.gather_results;
    histogram(cilia.Results.MedianFrequencyVec); %% plot histogram 
    fig = get(groot,'CurrentFigure');
    saveas(fig,[cilia.Filename(1:end-6),'_freq_.jpg']);
    
    
    
    
    % add frequencies to R
    R=cat(1,cilia.Results.MedianFrequencyVec(:));
end


%%%% also for the other

path_dir = '/home/np451/Desktop/ependymal data/6.6/0.1mlmin/'
cd(path_dir);
mkdir('analysis');
d=dir('*.movie');
R=[]

for nf=10:20%size(d,1);
    
    filename= d(nf).name;
    cd(path_dir);
    cilia=DDM_Analysis(filename);
    cd('analysis');

    cilia.VariableBoxSize_Analysis([box_size]);
    save([cilia.Filename(1:end-5),'mat'],'cilia');  %%% save sruslts



    % plot good boxes with std %


    new_std = cilia.std_fs;
    new_std=   new_std(1:floor(f_x/box_size)*box_size , 1:floor(f_y/box_size)*box_size ) 
    new_mask = oversample(cilia.Results.ind_good_boxes,box_size);
    mask_overlay(new_std, new_mask, [1 0 0], 0.3);
    fig = get(groot,'CurrentFigure');
    saveas(fig,[cilia.Filename(1:end-5),'jpg']);
    close 'all'

    % plot frequencies
    cilia.gather_results;
    histogram(cilia.Results.MedianFrequencyVec); %% plot histogram 
    fig = get(groot,'CurrentFigure');
    saveas(fig,[cilia.Filename(1:end-6),'_freq_.jpg']);
    
    
    
    
    % add frequencies to R
    R=cat(1,cilia.Results.MedianFrequencyVec(:));
end



%%%%%%%

path_dir = '/home/np451/Desktop/ependymal data/6.6/0.5mlmin/'
cd(path_dir);
mkdir('analysis');
d=dir('*.movie');
R=[]

for nf=10:20%size(d,1);
    
    filename= d(nf).name;
    cd(path_dir);
    cilia=DDM_Analysis(filename);
    cd('analysis');

    cilia.VariableBoxSize_Analysis([box_size]);
    save([cilia.Filename(1:end-5),'mat'],'cilia');  %%% save sruslts



    % plot good boxes with std %


    new_std = cilia.std_fs;
    new_std=   new_std(1:floor(f_x/box_size)*box_size , 1:floor(f_y/box_size)*box_size ) 
    new_mask = oversample(cilia.Results.ind_good_boxes,box_size);
    mask_overlay(new_std, new_mask, [1 0 0], 0.3);
    fig = get(groot,'CurrentFigure');
    saveas(fig,[cilia.Filename(1:end-5),'jpg']);
    close 'all'

    % plot frequencies
    cilia.gather_results;
    histogram(cilia.Results.MedianFrequencyVec); %% plot histogram 
    fig = get(groot,'CurrentFigure');
    saveas(fig,[cilia.Filename(1:end-6),'_freq_.jpg']);
    
    
    
    
    % add frequencies to R
    R=cat(1,cilia.Results.MedianFrequencyVec(:));
end


%%%%%

path_dir = '/home/np451/Desktop/ependymal data/6.6/1mlmin/'
cd(path_dir);
mkdir('analysis');
d=dir('*.movie');
R=[]

for nf=10:20%size(d,1);
    
    filename= d(nf).name;
    cd(path_dir);
    cilia=DDM_Analysis(filename);
    cd('analysis');

    cilia.VariableBoxSize_Analysis([box_size]);
    save([cilia.Filename(1:end-5),'mat'],'cilia');  %%% save sruslts



    % plot good boxes with std %


    new_std = cilia.std_fs;
    new_std=   new_std(1:floor(f_x/box_size)*box_size , 1:floor(f_y/box_size)*box_size ) 
    new_mask = oversample(cilia.Results.ind_good_boxes,box_size);
    mask_overlay(new_std, new_mask, [1 0 0], 0.3);
    fig = get(groot,'CurrentFigure');
    saveas(fig,[cilia.Filename(1:end-5),'jpg']);
    close 'all'

    % plot frequencies
    cilia.gather_results;
    histogram(cilia.Results.MedianFrequencyVec); %% plot histogram 
    fig = get(groot,'CurrentFigure');
    saveas(fig,[cilia.Filename(1:end-6),'_freq_.jpg']);
    
    
    
    
    % add frequencies to R
    R=cat(1,cilia.Results.MedianFrequencyVec(:));
end


