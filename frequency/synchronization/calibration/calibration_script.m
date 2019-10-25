%%%% calibration script. runs through all movies to find amplitude

%data_dir='/u/homes/np451/Desktop/ependymal/sync/28.7.18/P3/';
data_dir='/media/np451/Seagate Expansion Drive/ependymal/syn/calibration/28.7.18/';
cd(data_dir);
load('track_results.mat');
%%

%%% select a mask from the movie without external flow %%%%
d=dir('*movie');
for jj=1:numel(d);
    filename=d(jj).name;
    C{jj}.filename=filename;
    Vind=strfind(filename,'V');
    Vind=strfind(filename,'V');
    Vend=strfind(filename,'_O');
    C{jj}.Volt= str2num(filename(Vind+1:Vend-1));
    Offind=strfind(filename,'O');
    C{jj}.Offset= str2num(filename(Offind+1));
    Pind=strfind(filename,'P');
    Endind=strfind(C{jj}.filename,'_Z');
    C{jj}.Period= str2num(filename(Pind+1:Endind-1));
    Zind= strfind(filename,'Z');
    C{jj}.Z= str2num(filename(Zind+1:Zind+2));
    Dind= strfind(filename,'D');
    C{jj}.D= str2num(filename(Dind+1:Dind+4));
%    C{jj}.A1=nan;
%    C{jj}.A2=nan;
    if C{jj}.Period == 50  & ( C{jj}.Volt == 3 )   ;
    %if C{jj}.Period ~=0 & C{jj}.Z == 7 & C{jj}.D == 250 & C{jj}.Volt = 2.5;
%%%%% load data and the roi
    mo=moviereader(filename);
    fs=mo.read();
    fs= double(fs)- (mean(fs,3));
    N_frames= size(fs,3);
    FR =mo.FrameRate;
    C{jj}.FR = FR;
   
   disp('press k for take a point and then press enter');
   scroll_stack(fs); 
   disp(jj)
   disp(filename)
   k = waitforbuttonpress;
   [x1,y1] = getpts; 
   k = waitforbuttonpress;
   [x2,y2] = getpts;
   k = waitforbuttonpress;
   [x3,y3] = getpts;
   k = waitforbuttonpress;
   [x4,y4] = getpts;
   C{jj}.A1= mean([pdist([x1,y1;x2,y2],'euclidean'),pdist([x3,y3;x4,y4],'euclidean')]);
   disp(num2str(C{jj}.A1))
%    k = waitforbuttonpress;
%    [x1,y1] = getpts; 
%    k = waitforbuttonpress;
%    [x2,y2] = getpts;
%    k = waitforbuttonpress;
%    [x3,y3] = getpts;
%    k = waitforbuttonpress;
%    [x4,y4] = getpts;
%    C{jj}.A2= mean([pdist([x1,y1;x2,y2],'euclidean'),pdist([x3,y3;x4,y4],'euclidean')]);
%    disp(num2str(C{jj}.A2))
   close all;
    end
end

%save('track_results.mat','C');

%%  plot defined Z and D
Volts=2:0.5:6;
Periods=[100,67,50,40,33,29];
Mat=zeros([numel(Volts),numel(Periods)]);
for jj=1:numel(C)
 if C{jj}.Period ~=0 & C{jj}.Z == 7 & C{jj}.D == 250;% & C{jj}.Volt==2;
   
     col= C{jj}.Volt/6;
     v=C{jj}.Volt;
     p=C{jj}.Period;
     fq= 1000/C{jj}.Period;
     Mat(v==Volts,p==Periods)= C{jj}.A1*0.14*fq;
 plot(fq,C{jj}.A1*0.14*fq,'ko','MarkerSize',7,'LineWidth',1,'MarkerFaceColor',[col,col,col]);hold on;

 
 
 end
end
xlabel('freq [Hz]');
ylabel('Oscillation Amplitude [um]');

figure(); 
plot(Volts,mean(Mat,2),'ko','MarkerSize',7,'LineWidth',1); 
xlabel('Volts');ylabel('average Amplitude [um]');
%% interpolate data from calibration matrix Mat (section above)
Volts=2:0.5:6;
Periods=[100,67,50,40,33,29];
int_points= 29:100;
int_Mat= zeros([numel(Volts),numel(int_points)]);
for ii=1:size(Mat,1)
    int_Mat(ii,:)= interp1(Periods,Mat(ii,:),int_points);
    
end

figure();
plot(1000./int_points,int_Mat,'-');hold on;
plot(1000./Periods,Mat,'o')

%save('calibration_matrix.mat','Mat','Periods','Volts','int_Mat','int_points');

%% plot defined V and Z (varying Distance

for jj=1:numel(C)
% if C{jj}.Period == 50 & C{jj}.Z == 7 & C{jj}.Volt == 3;% & C{jj}.Volt==2;
 if C{jj}.Period == 50 & C{jj}.Volt == 3;% & C{jj}.Volt==2;
     disp(jj);
     col= C{jj}.Z/15
 %    fq= 1000/C{jj}.Period;
     
 plot(C{jj}.D,C{jj}.A1*0.14,'ko','MarkerSize',7,'LineWidth',1,'MarkerFaceColor',[col,col,col]);hold on;

 end
end
xlabel('distance from walls [um]');
ylabel('Oscillation Amplitude [um]')
title('V3 20Hz Z07 and Z15')

%%

Volts=[2,4,6];
Periods= [29,33,40,50,67,100];
Mat= zeros([numel(Periods),numel(Volts)]);
for kk=1:numel(C);
   if isfield(C{kk},'A1')==0; C{kk}.A1=0;
   end
   if isfield(C{kk},'A2')==0; C{kk}.A2=0;
   end
   p= C{kk}.Period;
   v= C{kk}.Volt; 
   Mat1(p==Periods,v==Volts)= mean([C{kk}.A1,C{kk}.A2]);
    
end

plot(Mat1(1,:))


