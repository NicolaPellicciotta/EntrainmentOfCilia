%%%% this script load the data that I have found using a script that
%%%% unfortunately I have lost. I am an idiot. fortuantely in the class Tam
%%%% there is almost everything that I need. In the saved class I have
%%%% found some weird things about the phase that I want to solve
%load('Tam.mat');
load('Tam2.mat');


%% figure
Periods_plot=[40,41,42,43,44,46,48];
cleg=1;
figure('Renderer', 'painters', 'Position', [10 10 900 300])
for jj=1:numel(Tam.Periods_check)
    Period= double(Tam.Periods_check(jj));Volt=Tam.Volt;
    d=dir(strcat('*V',num2str(Volt),'*P',num2str(Period),'*.movie'));
    mo=moviereader(d(1).name); fps=mo.FrameRate;
    times=Tam.Res(jj).times;
    f_clock= 1000/(Period);
    phase_clock= 2*pi*f_clock*times;
    phase_cilium= (1:numel(times))*2*pi;
    
    dphi= (phase_clock(:)-phase_cilium(:))/(2*pi);
    dphi=dphi-dphi(1);
    if any(Period==Periods_plot)
    plot(times,dphi,'-','LineWidth',1.5);hold on;
    leg{cleg}=num2str(f_clock);cleg=cleg+1;
    %leg{cleg}=num2str(Period);cleg=cleg+1;

       end
end

ylim([-6,6]);
legend(leg,'FontSize',5);
xlabel('time [s]');
ylabel('phase ')




%% make nice plot for methods paper how to recognise the phase of a cilium. the real method was more complicated than only doing a threeshold bu I dont remember

cd '/media/np451/Seagate Backup Plus Drive/DATA/26.11.18/P28'

jj=1;
d=dir('*.movie')
mo=moviereader(d(1).name);
fps=mo.FrameRate;
pxint= Tam.Res(jj).pxint;
threshold=Tam.Res(jj).threeshold;
t=(1:numel(pxint))/fps;
[pks,locs,w,p] = findpeaks(pxint,t);
Locs= locs(pks>threshold);
Pks= pks(pks>threshold);

plot(t,pxint); hold on;
plot(Locs,Pks,'ro');


%% another figure fo the paper, plot the frequency distribution of the cell at rest and when there is an external flow

load('Cell.mat');
load('BW.mat');
Rc=Cell.Rc;


%Hz_toplot=[20,22,24,26]; 
Volt_toplot=[3];
%Period_toplot= unique(floor(1000./Hz_toplot));
%Period_toplot= [37,42,44,48];  %%% specific for P_28
Period_toplot= [40,42,44,46];  %%% specific for P_28
Hz_toplot=1000./Period_toplot;
ind=[];cc=1;

for jj=1:(numel(Rc)-1) 
if any((Rc{jj}.Period)== (Period_toplot)) & (Rc{jj}.Volt)== Volt_toplot ;
    mo= moviereader(Rc{jj}.filename); fs=mo.read(); 
    [fq,m_pxx] = average_fft(fs,BW{1},mo.FrameRate,1);
    Rc{jj}.fq=fq;
    Rc{jj}.m_pxx=m_pxx;
    ind(cc)=jj;cc=cc+1;end
end


for jj=1:(numel(Rc)-1) 
if (Rc{jj}.Period) == 0 & (Rc{jj}.Volt)== 4%Volt_toplot ;
    mo= moviereader(Rc{jj}.filename); fs=mo.read(); 
    [fq,m_pxx] = average_fft(fs,BW{1},mo.FrameRate,1);
    Rc{jj}.fq=fq;
    Rc{jj}.m_pxx=m_pxx;
    s= std(double(fs(:,:,1:200)),[],3); s=s/max(s(:));
    imshow(s);
    ind(cc)=jj;cc=cc+1;end
end

 plot_y= 0.05 +(0.9/6)*(1:5);
 dy=diff(plot_y);
cleg=1;
figure('Renderer', 'painters', 'Position', [10 10 500 500])


 for kk=1:numel(ind);
 jj=ind(kk);    
 %subplot(numel(Period_toplot)+1,1,kk);
 subplot('Position',[0.1 plot_y(kk) 0.8 dy(1) ]); 
 if Rc{jj}.Period==0;
     plot(Rc{jj}.fq(5:end),(Rc{jj}.m_pxx(5:end)),'k-','LineWidth',1.5); 
 %    legend('control','Location','east')
 else
     plot(Rc{jj}.fq(5:end),(Rc{jj}.m_pxx(5:end)),'k-','LineWidth',1.5); hold on;
     vline(1000/Rc{jj}.Period,'r','')%, strcat(num2str(1000/Rc{jj}.Period),'[Hz]') );
%     legend(num2str(1000/Rc{jj}.Period),'Location','east');
 end
 max_y=max((Rc{jj}.m_pxx(5:end))) +max((Rc{jj}.m_pxx(5:end)))/5; 
 
 xlim([min(Hz_toplot)-1,max(Hz_toplot)+1]);
 ylim([0,max_y])
 end
% legend(leg)
 %xlabel('frequency [Hz]'),ylabel('fft over px intensity');


%% picture with motion map and interrogation window 
d=dir('*V2*P0*.movie');
mo=moviereader(d(1).name);
fs=mo.read();
s= std(double(fs(:,:,1:200)),[],3);
imagesc(s)
rect=getrect()
[rows, cols]= rect2sub(rect);
fs=fs(rows,cols,:);
s= std(double(fs(:,:,1:200)),[],3);
s=s/max(s(:));
imshow(s);
rect=getrect()
mask=zeros(size(s));
[rows, cols]= rect2sub(rect);
mask(rows,cols)=1;
mask_overlay(s,mask,[1,0,0],0.3);
size(s)
