%% plots Arnaud tongue 
clear all;
close all;
data_dir= '/home/np451/Desktop/video to transer/rawdata'; %%% your directory
Cal=[]; 

%%%% load the results from the first step1_fft.m%%%%
cd(data_dir);
%%% this is the result of my analysis you can use this one if you like,
%%% otherwise just upload the one that you did with the script step1_fft
load('fft_results_fft_nicola.mat'); 

mkdir('plots');
     
 %%%% since we have only one cell we can put in Rc all the information of
 %%%% the first cell, this script was made for eventually have also more
 %%%% than a cell recorded in the same field ov view. 
   close all; clear Rc; nc=1; cc=1;
    for kk= nc:size(R,1):numel(R);    
    Rc{cc}=R{kk};
    cc=cc+1;
    end

%%%% get all the information about Period and Voltage right and then put in
%%%% the array 'ind' all the inex when P==0 and when Voltage Gain V<4 so
%%%% then to plot it
%%%
ind=[];
cc=1;


for jj=1:(numel(Rc)) 
        
    Pind=strfind(Rc{jj}.filename,'P');
    Rc{jj}.Period= str2num(Rc{jj}.filename(Pind+1:(end-25)));
    Vind=strfind(Rc{jj}.filename,'V');
    Vend=strfind(Rc{jj}.filename,'_O');
    Rc{jj}.Volt= str2num(Rc{jj}.filename(Vind+1:Vend-1));
    if (Rc{jj}.Period==0 ); ind(cc)=jj;cc=cc+1;end    
    Rc{jj}.Freq= 1000./Rc{jj}.Period;    %%% from the Period we get the externla frequency of the oscillatory flow 
end

%%%% since there can be multiple harmonic peaks, we give a first approximate
%%%%CBF of the cell from the analysis of the videos in the absence of external flow
%%%%plot all the spectra when P=0 e V<4 and decide the frequexy at rest of
%%%% the ciliated cell

figure(1);
cleg=1;
clear leg;
 for jj=ind;
 figure(1);plot(Rc{jj}.fq(5:end),smooth(Rc{jj}.m_pxx(5:end)),'LineWidth',2); hold on;
 %figure(2); plot(str2num(Hz{jj}),Rc{jj}.fp-str2num(Hz{jj}),'o'); hold on;
 leg{cleg}=num2str(Rc{jj}.Volt);cleg=cleg+1;
 end
 figure(1); legend(leg);
 xlabel('frequecy [Hz]'),ylabel('periodogram');
 xlim([12,35]);
  
prompt= 'what is the frequency guess ';
fq_guess = (input(prompt));

%%% populate results in the right range of frq, use the frequecy at rest
%%% of the cilium to evaluate a good range where to look for
%%% synchronisation
  fq_max= fq_guess + (7);
  fq_min  =fq_guess- (7);


  
  for jj=1:(numel(Rc))

    fq=Rc{jj}.fq;
    f_range=Rc{jj}.fq> (fq_min) & Rc{jj}.fq<(fq_max);
%    baseline= min((Rc{jj}.m_pxx(f_range)));

%%%% find the peaks frequency in the selected freq range
    [pks,locs,w,p] = findpeaks((Rc{jj}.m_pxx(f_range)),fq(f_range));
    baseline= 0;
    [~,ind_sort]= sort(pks);                                           %%%% sort peaks and get an index 
    pks=pks(ind_sort); locs=locs(ind_sort);w=w(ind_sort);p=p(ind_sort);  %%% order all the variables with the same index
    Rc{jj}.pks=pks; Rc{jj}.locs=locs; Rc{jj}.w=w; Rc{jj}.p=p;     %%%%% load results in Rc
    Rc{jj}.baseline=baseline;
    [~,where] = max(pks);                                       
      
 %%% fp1c is the frequency with higher peak. if there are no peaks, put to nan    
    if numel(locs)==0; Rc{jj}.fp1c =nan;   
        else Rc{jj}.fp1 = locs(end);
    end
%%% look for a second peak
%%% fp2c is the second frequency high second higher peak, it is a nan if it does not exixt
    if numel(locs)>1
        Rc{jj}.fp2 = locs(end-1);   
        Rc{jj}.pk2 = pks(end-1);
    else
        Rc{jj}.fp2 = 0; 
        Rc{jj}.pk2=0;
    end
  end

%%% now we need to estimate the CBF and its value in the periodogram for the cell 
%%% without external flow.   
%%% f_rest is an array with all the frequency at rest of the cell, f_peak the corresponding peak value of the fft  
%%% F_rest is the mean frequency and F_peak the average peak of the
%%% periodogram

f_rest=[];
 volt_exp=[];
 f_peak=[];
 cc=0;
 for jj=1:(numel(Rc)-1)
      if Rc{jj}.Period==0; 
          cc=cc+1; f_rest(cc)=Rc{jj}.fp1;f_peak(cc)= Rc{jj}.pks(end); 
          volt_exp(cc)= Rc{jj}.Volt;
      end
 end
 
 F_rest=mean(f_rest);  
 F_peak= mean(f_peak);
 F_pstd= std(f_peak);
 
 
 %%%% start to add result to the class Cell, this will be useful later for
%%%% analysing multiple data

Cell(nc).s=s;
Cell(nc).BW=BW; 
 Cell(nc).F_rest=mean(f_rest);
 Cell(nc).f_rest=f_rest;
 Cell(nc).f_peak=f_peak;
 Cell(nc).F_peak= mean(f_peak);
 Cell(nc).F_pstd= std(f_peak);
 Cell(nc).f_range=f_range; 
 
 %%%% load the calibaration matrix of the flow vs voltage and offset  
cd(data_dir); 
Cal= load('calibration_matrix.mat');


%%%%% now start the plot of the Arnold tounge
 for jj=1:(numel(Rc));

 v=Rc{jj}.Volt;
 p=Rc{jj}.Period;    
 
 %%%%%% find the flow using the Calibration matrix and put it in the Matrix
 %%%%%% Mat1
 if isempty(Cal) == 1
        flow= Rc{jj}.Volt;
 else   flow = Cal.Flow_int(Cal.volt_int==v , Cal.fq_int== 1000/p );  
 end

 fp1=Rc{jj}.fp1; %%% frequency peak corrected in a smaller range
 fpw=Rc{jj}.w(end);
 fp2=Rc{jj}.fp2;
 pk2= Rc{jj}.pk2;
 baseline= Rc{jj}.baseline;
  Periods=51:120;
 Volts=2:0.5:6;
 Mat1(p==Periods,v==Volts)= fp1-(1000/p);
 ff= f_rest(1);
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 %%%%%% set a value for the identification of a synchronisation event. 
 %%%%%% If there is in the singla a second peak that it is similar in 
 %%%%%% magnitude to the one of the cell in the absence of external flow,
 %%%%%% then it is not a synchronisation event. For more info read the
 %%%%%% paper. 
 %%%  This allows us to identify that signal is coming from the cell 
 %%%% and not from other sources such as possible period defocusing of the
 %%%% field of view.
 
 threshold=0.66; 
 Cell(nc).threshold=threshold;

 if p~=0  %% if there is an external flow (period differnt from 0)
     
     
%%% first condition to identify synchronisation: the frequency of the new 
%%%measured peak must be very close to the applied exernal frequency 
    peak1_cond=  abs(fp1-(1000/p))<= 0.5;
    
%%% second condition: there should not be another peak that it is similar 
%%% in magnitude and frequency to the signal of the cell in the absence of
%%% externla flow     
    peak2_cond=   abs(fp2-ff )> 3 | abs( (F_peak-pk2)/F_peak )> threshold ;  %% 20% of the inital peak 

    Rc{jj}.peak1_cond=peak1_cond;
    Rc{jj}.peak2_cond=peak2_cond;
 
    
 %%%%% plot Arnould tonge with black and white depending on the above condition 
 figure(3);
 if peak1_cond & peak2_cond

     plot(1000./p,flow,'ko','MarkerSize',7,'LineWidth',1,'MarkerFaceColor',[0,0,0]);hold on;
     Rc{jj}.syn=1;
 else
     plot(1000./p,flow,'ko','MarkerSize',7,'LineWidth',1,'MarkerFaceColor',[1,1,1]);hold on;
     Rc{jj}.syn=0;
 end
 
 %%%%% plot Arnould tonge with greyscale depending on the value Sm
 %figure(4);
 %plot(1000./p,flow,'ko','MarkerSize',7,'LineWidth',1,'MarkerFaceColor',Sm*[1,1,1]);hold on;
 
 end
 
 end
 
 cd('plots');
 
 figure(3); title(strcat('Cell ',num2str(nc),' mean freq',num2str(F_rest), 'GreyScale'));
 xlabel('External Frequency [Hz]')
 ylabel('Velocity [mm/s]');
 saveas(gca,strcat('Cell_',num2str(nc),'_BaW.jpg'))
 
 cd(data_dir);
 save('Cell.mat','Cell');

  %% you can use this part to check that everything is all right. 
%Here you plot the periodogram at different external frequency and you
%compare it with the one when there is not an external flow.

close all
% here put the range of external frequency in which you want to check the 
%periodogram, use a small range if you dont want to have a untidy plot :)

Hz_toplot=[21,25]; 
%%%% here the voltage (from 2,3,4,5)
Volt_toplot=[5];

Period_toplot= unique(floor(1000./Hz_toplot));
ind=[];
cc=1;
for jj=1:(numel(Rc)-1) 
if any((Rc{jj}.Period)< max(Period_toplot)) & (Rc{jj}.Period> min(Period_toplot)) & any((Rc{jj}.Volt)== Volt_toplot); ind(cc)=jj;cc=cc+1;end
end

for jj=1:(numel(Rc)-1) 
if (Rc{jj}.Period) == 0 & (Rc{jj}.Volt)== Volt_toplot ; ind(cc)=jj;cc=cc+1;end
end


figure(1);
cleg=1;
 for jj=ind;
 figure(1);plot(Rc{jj}.fq(5:end),(Rc{jj}.m_pxx(5:end)),'LineWidth',2); hold on;
 %figure(2); plot(str2num(Hz{jj}),Rc{jj}.fp-str2num(Hz{jj}),'o'); hold on;
 if Rc{jj}.Period==0;
     leg{cleg}='at rest'
 else leg{cleg}=num2str(1000/Rc{jj}.Period);cleg=cleg+1;
 end
 end
 figure(1); legend(leg)
 xlabel('frequency [Hz]'),ylabel('fft over px intensity');
 xlim([Hz_toplot(1)-3,Hz_toplot(end)+3]);
 