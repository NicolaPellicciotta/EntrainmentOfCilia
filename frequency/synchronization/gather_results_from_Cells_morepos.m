%%% gather results from Cell array in plots
clear all;

data_dir1='/media/np451/Seagate Expansion Drive/29.10.18/';
Positions1={'P1','P2','P3','P4','P5','P6','P7','P8','P9','P10'}
data_dir2='/media/np451/Seagate Backup Plus Drive/DATA/1.11.18/' ;
Positions2={'P1','P2','P3','P4','P5','P6','P7','P8','P9','P10','P11','P12','P13'}
data_dir3='/media/np451/Seagate Backup Plus Drive/DATA/7.11.18/';
Positions3={'P1','P2','P3','P4','P5','P6','P7'};


data_dir4='/media/np451/Seagate Backup Plus Drive/DATA/26.11.18/';
%Positions3={'P1','P2','P3','P4','P5','P6','P7'};

cd('/media/np451/Seagate Backup Plus Drive/DATA/synchronisation_paper/')

cc=1;
for dd=[1:30];
Positions4{cc}=strcat('P',num2str(dd)); cc=cc+1;
end
cd(data_dir4);


PT= [strcat(data_dir1,Positions1), strcat(data_dir2,Positions2)...
    ,strcat(data_dir3,Positions3),strcat(data_dir4,Positions4)];

    cc=1;
for ps=1:numel(PT)
cd(PT{ps});
a=load('Cell.mat');
temp_Ncells= (size(a.Cell,2));
for nc=1:1
    Cell(cc).Rc=     a.Cell(1).Rc; 
    Cell(cc).good=   a.Cell(1).good;
    Cell(cc).F_rest= a.Cell(1).F_rest;
    Cell(cc).Pos=    PT{ps};
%    Cell(cc).Noise=a.Cell(nc).Noise;
%    Cell(cc).Centrin=a.Cell(nc).Centrin;
    Cell(cc).Ncilia=a.Cell(1).Ncilia;
%    Cell(cc).Sp_noise= a.Cell(nc).Sp_noise;

    cc=cc+1;
end
%cd(data_dir2);    
end

%%%% load calibration matrix
cd('/u/homes/np451/Documents/MATLAB/frequency/synchronization/');
Cal= load('calibration_matrix.mat');
cd(data_dir4)
%% 
%%%% find  
Ncells= size(Cell,2);

%%%%% clean up from Cell class the cells tahat are not good
good_ind=zeros([1,Ncells]);
for nc=1:Ncells;
    if Cell(nc).good==1;
        good_ind(nc)=1;
    end    
end
Cell=Cell(logical(good_ind));
Ncells= size(Cell,2);
%%%%%%

%%%%%% loading the variable Freq on Cell
Volts=[2:0.5:6];

%%% which frequency

for nc=1:Ncells
    cc=1;Periods=[];clear Freq;
 for jj=1:(numel(Cell(nc).Rc)-1); if Cell(nc).Rc{jj}.Volt == 4; 
         Periods(cc)= Cell(nc).Rc{jj}.Period; cc=cc+1;
     end;end;
%%%%% occhio un po di casini con quando prende gli infiniti!!! da P=0;

 Freq=1000./double(Periods); Freq=sort(Freq);Freq=Freq(1:end-1);
 Cell(nc).Freq=Freq;
%%%%%
end


for nc=1:Ncells
    Freq=Cell(nc).Freq;
    Syn=zeros([numel(Volts),numel(Freq)]);
%    Syn=nan([numel(Volts),numel(Freq)]);
    for jj=1:numel(Cell(nc).Rc);
        Cell(nc).Rc{jj}.Freq= 1000./double(Cell(nc).Rc{jj}.Period);
        p= Cell(nc).Rc{jj}.Period;
        if p~=0% & p<100;
        v=Cell(nc).Rc{jj}.Volt;
        f=Cell(nc).Rc{jj}.Freq;
        Syn(v==Volts,f==Freq) = Cell(nc).Rc{jj}.syn;
        end
     end
    Cell(nc).Syn=Syn;
    %%%%%% integral for find the synchonisation region epsilon
    for v=1:numel(Volts); Cell(nc).eps(v)=  trapz(Freq,Cell(nc).Syn(v,:)); end;
end


%cd(data_dir);
%%% frequency distribution
F_rest=[];
for nc=1:Ncells;
    F_rest(nc)=Cell(nc).F_rest;
end

histogram(F_rest,5);
saveas(gca,'frequency_histogram.fig');
saveas(gca,'frequency_histogram.pdf');


%% plot results with image
    s= Cell.s;
    BW=Cell.BW
    whichVolt=6;
    v_eps= (Volts==whichVolt);
    Reds= cat(3,s(:,:,1),0*s(:,:,1),0*s(:,:,1));  
    imshow(imadjust(s));
    BWT=zeros(size(s));
    for nc=1:size(Cell,2);
        if Cell(nc).good
        BWT= BWT+BW{nc};
        [maskx,masky]=find(BW{nc});
        text(masky(1)-30,maskx(1)-30,num2str(Cell(nc).eps(v_eps),2),'Color','red','FontSize',24)
        end
        end
    hold on
    h = imshow(Reds); % Save the handle; we'll need it later
    hold off        
    set(h, 'AlphaData', BWT); 
    saveas(gca,'cell_analysed_epsilon.jpg')




%% epsilon vs frequency
figure();
whichVolt=3;
v_eps= (Volts == whichVolt);
clear EPS, clear FREQS;
EPS=[];
FREQS=[];
cc=1;
for nc=1:Ncells;
    if Cell(nc).good
    plot(Cell(nc).F_rest,Cell(nc).eps(v_eps),'o');hold on;
    leg{cc}=Cell(nc).Pos; 
    EPS(cc)=Cell(nc).eps(v_eps);
    FREQS(cc)=Cell(nc).F_rest;
    cc=cc+1;
    end
end


%legend(leg)
xlabel('Frequency at rest[Hz]');
ylabel('Epsilon [Hz] synch region');
title(strcat('Epsilon vs Frequency V=',num2str(whichVolt)));
%legend(leg)
saveas(gca,'epsilon_frequency.jpg');
saveas(gca,'epsilon_frequency.fig');

nbins=5;
figure();
[histw,ehistw,vinterval]=hist_nico(FREQS,EPS,nbins);
errorbar(vinterval,histw,ehistw,'ko','MarkerSize',10,'LineWidth',1.5);
xlabel('Frequency at rest[Hz]');
ylabel('Epsilon [Hz] synch region');
%set(gca, 'YScale', 'log')

%% join figures epsilon vs frequency
parent_directory= '/media/np451/Seagate Expansion Drive/DATA/'
directory_eperiments= {'13.10.18','16.10.18','19.10.18'};
cd(parent_directory);
close all

Xdata_array= zeros([numel(directory_eperiments),50]);
Xdata_array(:)=nan;
Ydata_array= zeros([numel(directory_eperiments),50]);
Ydata_array(:)=nan;


for ff=1:numel(directory_eperiments)

cd(directory_eperiments{ff});
open 'epsilon_frequency.fig';
D=get(gca,'Children'); %get the handle of the line object
Xdata{ff}=get(D,'XData'); %get the x data
Ydata{ff}=get(D,'YData'); %get the y data

for hh=1:numel(Xdata{ff}); Xdata_array(ff,hh)= Xdata{ff}{hh};end 
for hh=1:numel(Ydata{ff}); Ydata_array(ff,hh)= Ydata{ff}{hh};end 

%plot(Xdata_array, Ydata_array,'ko','MarkerSize',7,'LineWidth',1,'MarkerFaceColor',[ff,ff,ff]./(numel(directory_eperiments)+1)); hold on;
cd(parent_directory);
end

close all

figure(1);
for ff=1:numel(Xdata)
plot(Xdata_array(ff,:), Ydata_array(ff,:),'ko','MarkerSize',7,'LineWidth',1,'MarkerFaceColor',[ff,ff,ff]./(numel(directory_eperiments))); hold on;
end
xlabel('Frequency[Hz]');
ylabel('epsilon[Hz]'); 

hold on; 
%%%%%%%
bins=[10,13,16,20,23,26,30]
n_bins= numel(bins);
for b=1:(numel(bins)-1)
Eps_mean(b)=nanmean(Ydata_array(Xdata_array(:)>bins(b) & Xdata_array(:)<bins(b+1)));
end
plot((bins(1:end-1)+bins(2:end))/2,Eps_mean,'k-');

legend({'Day6','Day9','Day12','mean'});

saveas(gca,'epsilon_frequency_total.jpg');
saveas(gca,'epsilon_frequency_total.fig');


%% epsilon vs Volts
set(0,'defaulttextinterpreter','latex')

load('/media/np451/Seagate Backup Plus Drive/DATA/synchronisation_paper/allvariables_area.mat');
lw=0.3;
mks=14

%Cal=load('/u/homes/np451/Documents/MATLAB/frequency/synchronization/calibration_matrix.mat');


N_plot=[1,5,10,20 ;5,10,20,30];
%color_scatter_array={'ko','bo','ro','mo'};
%color_median_array={'k--o','k-->','k--<','k--^'};
color_median_array={'ko','r>','g<','b^'};
facecolor_median_array={'k','r','g','b'};


clear l;cleg=1;
for kk=1:size(N_plot,2)


Flow_mean= mean(Cal.Flow_int,2);
Volts_plot=[2,3,4,5];
Volts_plot2=[4,5];

Eps=nan([Ncells,numel(Volts)]);
volt_ind=zeros(size(Volts)); %%% to find an array with the index with the right voltage
volt_ind2=zeros(size(Volts));

for vv=1:numel(Volts_plot);
    volt_ind=volt_ind+ (Volts==Volts_plot(vv));    
end
volt_ind=logical(volt_ind);

for vv=1:numel(Volts_plot2);
    volt_ind2=volt_ind2+ (Volts==Volts_plot2(vv));    
end
volt_ind2=logical(volt_ind2);

num_good=zeros(size(volt_ind));
for nc=1:Ncells
       if Cell(nc).good & (Cell(nc).Ncilia>=N_plot(1,kk) & Cell(nc).Ncilia<N_plot(2,kk)) 
            if nc>28
                Eps(nc,volt_ind)= Cell(nc).eps(volt_ind);
            else
                Eps(nc,volt_ind2)= Cell(nc).eps(volt_ind2);                                
            end
       end
end
 
num_good= sum(~isnan(Eps));
Eps_mean= nanmedian(Eps,1); Eps_std= nanstd(Eps,1);%./sqrt(num_good);
Eps_mean(Eps_mean>9.5)=nan; 

%figure()
figure(1)
errorbar(Flow_mean(volt_ind),Eps_mean(volt_ind),Eps_std(volt_ind),color_median_array{kk},'MarkerSize',mks,'LineWidth',lw,'MarkerFaceColor',facecolor_median_array{kk});hold on;
%xlabel('flow velocity [mm/s]');
%ylabel('\epslion[Hz]');
%saveas(gca,'epsilon_voltage.jpg');
%saveas(gca,'epsilon_voltage.fig');
hold on;

end

legend({'N_{c}[1,4]', 'N_{c}[5,9]', 'N_{c}[10,19]', 'N_{c}[20,25]'},'FontSize',10,'Location','northwest')
ylim([0,13])
xlim([0,2])
set(gca,'FontSize',15);
xlabel('$v_{EX} [mm/s]$','FontSize',20);
ylabel('$\epsilon $ [Hz]','FontSize',20);



 %% join figures from different experiments in one, in this case the epsilon vs voltage

parent_directory= '/media/np451/Seagate Expansion Drive/DATA/'
directory_eperiments= {'13.10.18','16.10.18','19.10.18','22.10.18'};
cd(parent_directory);
close all


load('/home/np451/Documents/MATLAB/frequency/synchronization/calibration/calibration_matrix.mat')


for ff=1:numel(directory_eperiments)
cd(directory_eperiments{ff});
open 'epsilon_voltage.fig';
D=get(gca,'Children'); %get the handle of the line object
Xdata{ff}=get(D,'XData'); %get the x data
Ydata{ff}=get(D,'YData'); %get the y data
close all
%plot(XData, YData,'ko','MarkerSize',7,'LineWidth',1,'MarkerFaceColor',[ff,ff,ff]./(numel(directory_eperiments)+1)); hold on;
cd(parent_directory);
end 
figure(1);

mean_Flow= mean(Flow_int(:,fq_int>13),2);
for ff=1:numel(Xdata)
Xdata_int=[]; for vv=1:numel(Xdata{ff}); Xdata_int(vv)= mean_Flow(volt_int==Xdata{ff}(vv));end;     
%plot(Xdata, Ydata{ff},'-ko','MarkerSize',7,'LineWidth',1,'MarkerFaceColor',[ff,ff,ff]./(numel(directory_eperiments))); hold on;
plot(Xdata_int, Ydata{ff},'-ko','MarkerSize',7,'LineWidth',1,'MarkerFaceColor',[ff,ff,ff]./(numel(directory_eperiments))); hold on;
end
xlabel('Average Flow [mm/s]');
ylabel('epsilon[Hz]'); 
legend({'Day6','Day9','Day12','Day15'})
saveas(gca,'epsilon_voltage_total.jpg');
 saveas(gca,'epsilon_voltage_total.fig');
 

 
 
 %% epsilon vs N_cilia
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
EPS=[]; NCILIA=[];
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
xlabel('$N_{c}$ number of cilia','Interpreter','latex', 'FontSize',20);
ylabel('$\epsilon [Hz]$ ','Interpreter','latex', 'FontSize',20);
title(strcat('average flow=  ',num2str(Flow_mean(v_eps),'%.1f'),'[mm/s]'),'Interpreter','latex', 'FontSize',15);
saveas(gca,strcat('average flow  ',num2str(Flow_mean(v_eps)),'.jpg'));


nbins=5;
figure(5);
[histw,ehistw,vinterval]=hist_nico(NCILIA,EPS,nbins,[1:5:30]);
histw(histw>9)=nan;  %%%% to remove from statistic the one the went out of range 
%errorbar(vinterval,histw,ehistw,color_median,'MarkerSize',mks,'LineWidth',lw,'MarkerFaceColor',color_face);hold on;
errorbar(vinterval,histw./(Flow_mean(v_eps)),ehistw./(Flow_mean(v_eps)),color_median,'MarkerSize',mks,'LineWidth',lw,'MarkerFaceColor',color_face);hold on;

xlabel('Number of cilia per cell','FontSize',20);
ylabel('$\epsilon $ [Hz]','FontSize',20);
ylim([0,20]);
set(gca,'FontSize',16);
%title(strcat('average flow ',num2str(Flow_mean(v_eps)),'[mm/s]'));
end
 legend(l,'FontSize',13);
 set(gca,'XScale', 'log', 'YScale', 'log');

%set('XScale', 'log', 'YScale', 'log');
%% flagellar noise and frequency
figure();
cc=1;
for nc=1:Ncells
plot(Cell(nc).F_rest,Cell(nc).Noise.value,'o');hold on;
leg{cc}=Cell(nc).Pos(end-2:end); ;cc=cc+1;
end
legend(leg)
xlabel('freq[Hz]'); ylabel('Noise');

%% %% epsilon vs flagellar noise
figure();
v_eps= (Volts==4);
cc=1;
for nc=1:Ncells;
    plot(Cell(nc).Noise.value,Cell(nc).eps(v_eps),'o','MarkerSize',8,'LineWidth',3,'MarkerFaceColor',[1,1,1]);hold on;
    leg{cc}=num2str(nc); ;cc=cc+1;

end
legend(leg)
xlabel('noise (std(f)/<f>2)');
ylabel('Epsilon [Hz] synch region');

%% epsilon vs spatial noise

figure();
v_eps= (Volts==4);
cc=1;
clear EPS; clear STDFREQ; clear NCILIA;

EPS=[]; STDFREQ=[]; NCILIA=[];
for nc=1:Ncells;
    F=Cell(nc).Sp_noise.F(:);
    box_size= Cell(nc).Sp_noise.box_size*0.13/40*60 %%%% boxsize in um using 60X
    
    diff_med=F-nanmedian(F); median_var= median(abs(diff_med(~isnan(diff_med))));
    std_freq= Cell(nc).Sp_noise.std_freq/Cell(nc).F_rest;
    N_boxes= numel(F(~isnan(F)))*box_size^2;
%    plot(median_var/sqrt(N_boxes) , Cell(nc).eps(v_eps),'o','MarkerSize',8,'LineWidth',3,'MarkerFaceColor',[1,1,1]);hold on;
    plot(std_freq/sqrt(N_boxes),Cell(nc).eps(v_eps),'o','MarkerSize',8,'LineWidth',3,'MarkerFaceColor',[1,1,1]);hold on;
%     plot(N_boxes,Cell(nc).eps(v_eps),'o','MarkerSize',8,'LineWidth',3,'MarkerFaceColor',[1,1,1]);hold on;

    EPS(cc)=Cell(nc).eps(v_eps);
    STDFREQ(cc)= std_freq/sqrt(N_boxes);
    NCILIA(cc)= Cell(nc).Ncilia;
%    STDFREQ(cc)=median_var/sqrt(N_boxes)
    NBOXES(cc)=N_boxes;
    leg{cc}=num2str(nc); ;cc=cc+1;

end
legend(leg)
xlabel('spatial noise [Hz]');
ylabel('Epsilon [Hz] synch region');

nbins=5;
figure();
[histw,ehistw,vinterval]=hist_nico(STDFREQ,EPS,nbins);
errorbar(vinterval,histw,ehistw,'ko','MarkerSize',7,'LineWidth',1);
xlabel('<f^2>/<f> [HZ]');
ylabel('Epsilon [Hz] synch region');

nbins=6;
figure();
[histw,ehistw,vinterval]=hist_nico(NBOXES,EPS,nbins);
errorbar(vinterval,histw,ehistw,'ko','MarkerSize',7,'LineWidth',1);
xlabel('cilia beating area [um2]');
ylabel('Epsilon [Hz] synch region');

nbins=6;
figure();
[histw,ehistw,vinterval]=hist_nico(NCILIA,NBOXES,nbins);
errorbar(vinterval,histw,ehistw,'ko','MarkerSize',7,'LineWidth',1);
ylabel('cilia beating area [um2]');
xlabel('Number of cilia per cell');




%% Centrin area and intensity vs epsilon

figure();
v_eps= (Volts==4);
cc=1;
for nc=1:Ncells;
    if isempty(Cell(nc).Centrin.movie)==0
    plot(Cell(nc).Ncilia./ Cell(nc).Centrin.size_mu,Cell(nc).eps(v_eps),'o','MarkerSize',8,'LineWidth',3,'MarkerFaceColor',[1,1,1]);hold on;
    leg{cc}=Cell(nc).Pos(end-2:end); ;cc=cc+1;
    end
end
legend(leg)
xlabel('centrin area ');
ylabel('Epsilon [Hz] synch region');


%% Number of cilia counted manually vs epsilon 


figure();
v_eps= (Volts==4);
cc=1;
for nc=1:Ncells;
    [~,arg_min]=min(abs(Cal.fq_int-Cell(nc).F_rest))
    flow_norm = Cal.Flow_int(v_eps,arg_min);
%%%%normilised with flow
%    plot(Cell(nc).Ncilia, Cell(nc).eps(v_eps)/flow_norm,'o','MarkerSize',8,'LineWidth',3,'MarkerFaceColor',[1,1,1]);hold on;   
%%%% normalised with frequency
%    plot(Cell(nc).Ncilia,Cell(nc).eps(v_eps)/Cell(nc).F_rest,'o','MarkerSize',8,'LineWidth',3,'MarkerFaceColor',[1,1,1]);hold on;
%%%% not normalised
    plot(Cell(nc).Ncilia,Cell(nc).eps(v_eps),'o','MarkerSize',8,'LineWidth',3,'MarkerFaceColor',[1,1,1]);hold on;

leg{cc}=Cell(nc).Pos(end-2:end); ;cc=cc+1;
end
legend(leg)
xlabel('N cilia ');
ylabel('Epsilon [Hz] synch region');
%%  Number of cilia and spatial frequency deviation


figure();
v_eps= (Volts==4);
cc=1;
for nc=1:Ncells;

    plot(Cell(nc).Ncilia,Cell(nc).Sp_noise.std_freq,'o','MarkerSize',8,'LineWidth',3,'MarkerFaceColor',[1,1,1]);hold on;

leg{cc}=Cell(nc).Pos(end-2:end); ;cc=cc+1;
end
legend(leg)
xlabel('N cilia ');
ylabel('Epsilon [Hz] synch region');



%% Number of cilia and Frequency
figure();
v_eps= (Volts==4);
cc=1;FREQ=[];NCILIA=[];

lw=1;
mks=14


for nc=1:Ncells;
    %plot(Cell(nc).Ncilia, Cell(nc).F_rest,'bo','LineWidth',lw,'MarkerSize',mks);hold on;
    leg{cc}=Cell(nc).Pos(end-2:end); 

    FREQ(cc)=Cell(nc).F_rest;
    NCILIA(cc)= Cell(nc).Ncilia;
%    STDFREQ(cc)=median_var/sqrt(N_boxes)
    leg{cc}=num2str(nc); ;cc=cc+1;



end
%legend(leg)

plot(NCILIA, FREQ,'bo','LineWidth',lw,'MarkerSize',mks);hold on;
set(gca,'FontSize',16);



nbins=6;
%figure();
[histw,ehistw,vinterval]=hist_nico(NCILIA(NCILIA<23),FREQ(NCILIA<23),nbins);
errorbar(vinterval,histw,ehistw,'--ko','LineWidth',lw,'MarkerSize',mks,'MarkerFaceColor','k');
xlabel('Number of cilia per cell','FontSize',20);
ylabel('CBF [Hz]','FontSize',20);
legend({'single cell','binnig average'})
ylim([10,32])

%%
figure(); 
Ncilia_sim=[1,2,3,4,6,8,10,20,30];
errorbar(vinterval(1:end-1),histw(1:end-1),ehistw(1:end-1),'--ko','MarkerSize',7,'LineWidth',1);
hold on;
%plot(aa.Ncilia_sim, aa.medfd0d4*15,'--r<','MarkerSize',7,'LineWidth',1);
%plot(aa.Ncilia_sim, aa.medfd0d6*15,'--b>','MarkerSize',7,'LineWidth',1);
%plot(aa.Ncilia_sim, aa.medfd1d2*15,'--mv','MarkerSize',7,'LineWidth',1);
legend({'exp'});
%legend({'exp','d=4','d=6','d=12'});
xlabel('Ncilia [HZ]');
ylabel('Frequency [Hz]');
xlim([0,20]);


%%
figure(); 
Ncilia_sim=[1,2,3,4,6,8,10,20,30];
errorbar(vinterval(1:end-1),histw(1:end-1),ehistw(1:end-1),'--ko','MarkerSize',7,'LineWidth',1);
hold on;
plot(Ncilia_sim, bb.fmedb1d2*15,'--r<','MarkerSize',7,'LineWidth',1);
plot(Ncilia_sim, bb.fmedb1d5*15,'--b>','MarkerSize',7,'LineWidth',1);
plot(Ncilia_sim, bb.fmedb1d8*15,'--mv','MarkerSize',7,'LineWidth',1);
legend({'exp'});
%legend({'exp','d=4','d=6','d=12'});
xlabel('Ncilia [HZ]');
ylabel('Frequency [Hz]');
xlim([0,20]);




