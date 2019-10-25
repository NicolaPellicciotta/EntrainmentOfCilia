%% This is a script for the plot of the Figure 1 version 25.10.19

clear all
set(0,'defaulttextinterpreter','latex')

load('/media/np451/Seagate Backup Plus Drive/DATA/synchronisation_paper/allvariables_area.mat');
lw=0.3;
mks=14

%Cal=load('/u/homes/np451/Documents/MATLAB/frequency/synchronization/calibration_matrix.mat');


N_plot=[1,5,10,20 ;5,10,20,30];
%color_scatter_array={'ko','bo','ro','mo'};
%color_median_array={'k--o','k-->','k--<','k--^'};
color_median_array={'k--o','r-->','g--<','b--^'};
facecolor_median_array={'k','r','g','b'};

figure(98)
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
errorbar(Flow_mean(volt_ind),Eps_mean(volt_ind),Eps_std(volt_ind),color_median_array{kk},'MarkerSize',mks,'LineWidth',lw,'MarkerFaceColor',facecolor_median_array{kk});hold on;
%xlabel('flow velocity [mm/s]');
%ylabel('\epslion[Hz]');
%saveas(gca,'epsilon_voltage.jpg');
%saveas(gca,'epsilon_voltage.fig');
hold on;

end

%legend({'N_{c}\in[1,4]', 'N_{c}\in[5,9]', 'N_{c}\in[10,19]', 'N_{c}\in[20,25]'},'FontSize',10,'Location','northwest')
legend({'N \in[1,4]', 'N \in[5,9]', 'N \in[10,19]', 'N \in[20,25]'},'FontSize',13,'Location','northwest')

ylim([0,13])
xlim([0,2])
set(gca,'FontSize',15);
%xlabel('$v_{EX} [mm/s]$','FontSize',20);
%ylabel('$\epsilon $ [Hz]','FontSize',20);
x0=0;
y0=0;
width=500;
height=400;


fig=figure(98);
%saveas(fig,'/u/homes/np451/Dropbox/Synchonisation of mammalian cilia with hydrodynamic forces/PNAS/figures_review/epsilon_vex_std.pdf')


