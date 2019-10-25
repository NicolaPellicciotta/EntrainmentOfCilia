%% Getting the Q values
%channel_factor(0.007, 0, 1, 1) % (z,y,w,h)
cd '/home/np451/Documents/MATLAB/frequency/synchronization/calibration/'

a=load('track_results_y0.mat');C0=a.C;
a=load('track_results_y250.mat');C250=a.C;
a=load('track_results_yneg250.mat');Cneg250=a.C;
a=load('track_results_with20x_0.mat');C20_0=a.C;
a=load('track_results_with20x_250.mat');C20_250=a.C;
a=load('track_results_with20x_neg250.mat');C20_neg250=a.C;




%%%%%%%% this is for the flow rate with the 40X

[u07_neg250, u14_neg250, u21_neg250] = get_velocities_fromC(Cneg250);
Q07_neg250 = u07_neg250 ./ channel_factor(0.007, 0, 1, 1);
Q14_neg250 = u14_neg250 ./ channel_factor(0.014, 0, 1, 1);
Q21_neg250 = u21_neg250 ./ channel_factor(0.021, 0, 1, 1);

[u07_250, u14_250, u21_250] = get_velocities_fromC(C250);
Q07_250 = u07_250 ./ channel_factor(0.007, 0, 1, 1);
Q14_250 = u14_250 ./ channel_factor(0.014, 0, 1, 1);
Q21_250 = u21_250 ./ channel_factor(0.021, 0, 1, 1);

[u07_0, u14_0, u21_0] = get_velocities_fromC(C0);
Q07_0 = u07_0 ./ channel_factor(0.007, 0, 1, 1);
Q14_0 = u14_0 ./ channel_factor(0.014, 0, 1, 1);
Q21_0 = u21_0 ./ channel_factor(0.021, 0, 1, 1);

Qstack = cat(3, Q07_neg250, Q14_neg250, Q21_neg250, Q07_250, Q14_250, Q21_250, Q07_0, Q14_0, Q21_0);
Qstack(Qstack==0 | Qstack<0 )=nan;
Qmeans = nanmean(Qstack, 3);

u07 =   nanmean( calibration_remove_nan(cat(3,u07_neg250,u07_0,u07_250)),3); 
u14 = nanmean( calibration_remove_nan(cat(3,u14_neg250,u14_0,u14_250)),3);
u21 = nanmean(calibration_remove_nan(cat(3,u21_neg250,u21_0,u21_250)),3);

%Qstd = nanstd(Qstack,[], 3);

%surf(Qstd) % There are outlier points that should be excluded...

%%%%%%%% this is for the flow rate with the 40X

[u07_neg250, u14_neg250, u21_neg250] = get_velocities_fromC(C20_neg250);
Q07_neg250 = u07_neg250 ./ channel_factor(0.007, 0, 1, 1);
Q14_neg250 = u14_neg250 ./ channel_factor(0.014, 0, 1, 1);
Q21_neg250 = u21_neg250 ./ channel_factor(0.021, 0, 1, 1);

[u07_250, u14_250, u21_250] = get_velocities_fromC(C20_250);
Q07_250 = u07_250 ./ channel_factor(0.007, 0, 1, 1);
Q14_250 = u14_250 ./ channel_factor(0.014, 0, 1, 1);
Q21_250 = u21_250 ./ channel_factor(0.021, 0, 1, 1);

[u07_0, u14_0, u21_0] = get_velocities_fromC(C20_0);
Q07_0 = u07_0 ./ channel_factor(0.007, 0, 1, 1);
Q14_0 = u14_0 ./ channel_factor(0.014, 0, 1, 1);
Q21_0 = u21_0 ./ channel_factor(0.021, 0, 1, 1);

Q20_stack = cat(3, Q07_neg250, Q14_neg250, Q21_neg250, Q07_250, Q14_250, Q21_250, Q07_0, Q14_0, Q21_0);
Q20_stack(Q20_stack==0 | Q20_stack==-1 )=nan;
Q20_means = nanmean(Q20_stack, 3);   %%%% it is moltipled by 2 because it is with 20X
%Q20_std = nanstd(Q20_stack,[], 3);

%%%%% replace dodgy values from the 40X with the one at 20X
Qmeans(9,1:4)=Q20_means(9,1:4);
Qmeans(8,1:2)=Q20_means(8,1:2);
Qmeans(7,1)=Q20_means(7,1);
Qmeans(6,1)=Q20_means(6,1);



%% velocities plots
Z= 0.01;%%%% in mm
figure();

Periods = [100,80,66,57,50,44,40,36,33,28];
Volts=[2,2.5,3,3.5,4,4.5,5,5.5,6];
plot_volts= [2,3,4,5];
Frequencies= 1000./Periods;
lcc=1;
Flow=Qmeans.*channel_factor(0.07, 0, 1, 1);  
Flow14=Qmeans.*channel_factor(0.014, 0, 1, 1);
%%%%% flow is defined for [2,2.5,3,3.5,4,4.5,5,5.5,6] and ten periods


for volt=1:size(Qmeans,1)
    
    if ~isempty(plot_volts(Volts(volt) == plot_volts)) ; 
    plot(Frequencies, Flow14(volt,:),'ko-','MarkerSize',7,'LineWidth',1,'MarkerFaceColor',volt./numel(Volts)*[1,1,1]);
    hold on;
    l{lcc}=num2str(Volts(volt));
    lcc=lcc+1;
end
end

legend(l)
xlim([10,35]);
xlabel('external frequency [Hz]');
ylabel('flow [mm/s]');
%% velocities plots
Z= 0.01;%%%% in mm
figure();

Periods = [100,80,66,57,50,44,40,36,33,28];
Volts=[2,2.5,3,3.5,4,4.5,5,5.5,6]
Frequencies= 1000./Periods;
lcc=1;
Flow=Qmeans.*channel_factor(0.007, 0, 1, 1);

for volt=1:size(Qmeans,1)

    plot(Frequencies, Flow(volt,:),'ko-','MarkerSize',7,'LineWidth',1,'MarkerFaceColor',volt./numel(Volts)*[1,1,1]);
    hold on;
    l{lcc}=num2str(Volts(volt));
    lcc=lcc+1;
end

legend(l)
xlim([10,35]);
xlabel('external frequency [Hz]');
ylabel('flow [mm/s]');

%% plot Flow as function of voltage

Z= 0.01;%%%% in mm
figure();

Periods = [100,80,66,57,50,44,40,36,33,28];
Volts=[2,2.5,3,3.5,4,4.5,5,5.5,6]
Frequencies= 1000./Periods;
lcc=1;
Flow=Qmeans.*channel_factor(Z, 0, 1, 1);
%ErFlow=Qstd.*channel_factor(Z, 0, 1, 1);

for fq=1:size(Qmeans,2)

%    errorbar(Volts, Flow(:,fq)./Flow(1,fq),ErFlow(:,fq),'ko-','MarkerSize',7,'LineWidth',1,'MarkerFaceColor',fq./numel(Frequencies)*[1,1,1]);
    plot(Volts, Flow(:,fq)./Flow(1,fq),'ko-','MarkerSize',7,'LineWidth',1,'MarkerFaceColor',fq./numel(Frequencies)*[1,1,1]);
    hold on;
    l{lcc}=num2str(Frequencies(fq));
    lcc=lcc+1;
end

legend(l)
xlabel('Voltage [V]');
ylabel('normalised flow with flow at V==2');
%% compare velocity at different high

figure();
p07=nanmean(u07(:)./u07(:)); p07E=nanstd(u07(:)./u07(:))%/sqrt(numel(u07));
p14=nanmean(u14(:)./u07(:)); p14E=nanstd(u14(:)./u07(:))%/sqrt(numel(u07));
p21=nanmean(u21(:)./u07(:)); p21E=nanstd(u21(:)./u07(:))%/sqrt(numel(u07));

%%%% experimental  average ratio between velocity at  
errorbar([7,14,21],[p07,p14,p21],[p07E,p14E,p21E],'ko','MarkerSize',7,'LineWidth',1);
hold on;
tp07 = channel_factor(0.007, 0, 1, 1)/channel_factor(0.007, 0, 1, 1);
tp12 = channel_factor(0.014, 0, 1, 1)/channel_factor(0.007, 0, 1, 1);
tp21 = channel_factor(0.021, 0, 1, 1)/channel_factor(0.007, 0, 1, 1);
plot([7,14,21],[tp07,tp12,tp21],'-ko','MarkerSize',7,'LineWidth',1);

legend({'exp','theory'});

%% interpolation

%%%%% linear interpolation of the data

Periods = [100,80,66,57,50,44,40,36,33,28];
Volts=[2,2.5,3,3.5,4,4.5,5,5.5,6]
Frequencies= 1000./Periods;

period_int=[28:100];
fq_int= sort(1000./period_int);  %%%% interpolated range of frequencies
[X,Y]=meshgrid(Frequencies,Volts);
V= Qmeans;
[Xq,Yq]=meshgrid(fq_int,Volts);
Vq = interp2(X,Y,V,Xq,Yq);

Z= 0.01;%%%% in mm
Flow_int=Vq.*channel_factor(0.007, 0, 1, 1);

figure();
lcc=1;  %%%% counting number for the legend
for volt=1:size(Flow_int,1)

    plot(fq_int, Flow_int(volt,:),'ko-','MarkerSize',7,'LineWidth',1,'MarkerFaceColor',volt./numel(Volts)*[1,1,1]);
    hold on;
    l{lcc}=num2str(Volts(volt));
    lcc=lcc+1;
end

legend(l)
xlim([10,35]);
xlabel('external frequency [Hz]');
ylabel('flow [mm/s]');

volt_int = Volts;
save('calibration_matrix.mat','fq_int','period_int','Flow_int','volt_int');
