%%% gather results from Cell array in plots
clear all;
%data_dir='/media/np451/Seagate Expansion Drive/ependymal/syn/29.7.18/P1';
%data_dir='/media/np451/Seagate Expansion Drive1/DATA/13.10.18/P5/';
%%%% for ALI 15 
% data_dir='/media/np451/Seagate Expansion Drive/DATA/22.10.18/';
% Positions={'P1','P2','P3'};

% %%%%% for ALI 12
% data_dir='/media/np451/Seagate Expansion Drive/DATA/19.10.18/'
% Positions={'P5','P6'};

% %%%% for ALI 9
% data_dir='/media/np451/Seagate Expansion Drive/DATA/16.10.18/'
% Positions={'P3','P4'};

% %%%%% for ALI 6
% data_dir='/media/np451/Seagate Expansion Drive/DATA/13.10.18/'
% Positions={'P4','P5','P6'};
data_dir1='/media/np451/Seagate Expansion Drive/29.10.18/';
Positions1={'P1','P2','P3','P4','P5','P6','P7','P8','P9','P10'}
data_dir2='/media/np451/Seagate Backup Plus Drive/DATA/1.11.18/' ;
Positions2={'P1','P2','P3','P4','P5','P6','P7','P8','P9','P10','P11','P12','P13'}
data_dir3='/media/np451/Seagate Backup Plus Drive/DATA/7.11.18/';
Positions3={'P1','P2','P3','P4','P5','P6','P7'};


data_dir4='/media/np451/Seagate Backup Plus Drive/DATA/26.11.18/';
%Positions3={'P1','P2','P3','P4','P5','P6','P7'};

cc=1;
for dd=[1:30];
Positions4{cc}=strcat('P',num2str(dd)); cc=cc+1;
end
cd(data_dir4);


PT= [strcat(data_dir1,Positions1), strcat(data_dir2,Positions2)...
    ,strcat(data_dir3,Positions3),strcat(data_dir4,Positions4)];


cd(data_dir);
%% for counting number of cilia in a cell
for dd=3:30numel(PT)
    disp(dd);
    cd(PT{dd});
    load('Cell.mat');

    d=dir('*P0*movie');
    mo=moviereader(d(1).name);
    fs=mo.read();

     fsbw= double(fs(:,:,2:300))- movmean(fs(:,:,2:300),100,3);
     s= std(double(fs(:,:,1:200)),[],3);
    clear fsbw2; 
     for kk=1:size(fsbw,3)
     fsbw2(:,:,kk)= wiener2(fsbw(:,:,kk),[5,5]);
     end
     k=0; clear x1;clear x2; clear y1;clear y2;
     
     scroll_stack(fsbw2)
     
   k = waitforbuttonpress;
   [x1,y1] = getpts; 
   k = waitforbuttonpress;
   [x2,y2] = getpts;

 
   A2(dd)= pdist([x1,y1;x2,y2],'euclidean');    %%%%% ampiezza
   Ncilia2(dd)= Cell.Ncilia;
   close all
   
end

%save('Amplitude_vs_Ncilia2.mat','Ncilia2','A2');

%% figure measure amplitude
imagesc(s);hold on;
plot([x1,x2],[y1,y2],'r-')

%%
load('Amplitude_vs_Ncilia.mat');
load('Amplitude_vs_Ncilia2.mat')

A=cat(1,A(:),A2(:));
Ncilia= cat(1,Ncilia(:),Ncilia2(:))



nbins=5;
ppm=0.13*40/60;

%Ncilia=Ncilia((A*ppm)>7.5);
%A=A((A*ppm)>7.5);
figure();
[histw,ehistw,vinterval]=hist_nico(Ncilia,A*ppm,nbins);
errorbar(vinterval,histw,ehistw,'-ko','MarkerSize',10,'MarkerFaceColor','k','LineWidth',1.5);hold on;
set(0,'defaulttextinterpreter','latex')
xlabel('Number of cilia per cell','FontSize',20);
ylabel('Beating Amplitude [$\mu m$]','FontSize',20);

%figure();
plot(Ncilia,A*ppm,'bo','MarkerSize',10,'LineWidth',2);
xlabel('Number of cilia per cell','FontSize',20);
ylabel('Beating Amplitude [$\mu m$]','FontSize',20);

legend('median binning average','single cell data');


%% figure centrin 


dd= 14;

    disp(dd);
    cd(strcat(data_dir,directories{dd}));
    load('Cell.mat');

    d=dir('*P0*movie');
    mo=moviereader(d(1).name);
    fs=mo.read();

     fsbw= double(fs(:,:,2:300))- movmean(fs(:,:,2:300),100,3);
     s= std(double(fs(:,:,1:200)),[],3);
     
     d=dir('*centrin*movie');
    mo=moviereader(d(1).name);
    fs=mo.read();
    C= mean(fs,3);
    C=C/max(C(:));
    figure(1)
    imshow(mat2gray(C));
    saveas(gcf,'centriolar_patch_P28.jpg');
    figure(2);
    imshow(mat2gray(s))
    saveas(gcf,'standard_deviation_P28.jpg');
     
