data_dir='/media/np451/Seagate Backup Plus Drive/DATA/7.11.18/'  %%% experiment in DMEM P/S
cc=1;
for dd=[1:7];
directories{cc}=strcat('P',num2str(dd)); cc=cc+1;
end
cd(data_dir);
%% for counting number of cilia in a cell
for dd=1:numel(directories)
    disp(dd);
    cd(strcat(data_dir,directories{dd}));
    load('Cell.mat');

    
    d=dir('*P0*movie');
    mo=moviereader(d(1).name);
    fs=mo.read();
    fsbw= double(fs(:,:,2:300))- mean(fs(:,:,2:300),3);
    k=0;
    scroll_stack(fsbw)
    prompt= 'number of cilia';
    Ncilia = (input(prompt));
    if numel(Cell)>1;
        Cell(1).Ncilia=Ncilia;    
    else
        Cell(1).Ncilia=Ncilia;
    end
    save('Cell.mat','Cell');
    close all;
    
end
%% for addin Centrin and Noise
for dd=1:numel(directories)
    disp(dd);
    cd(strcat(data_dir,directories{dd}));
    load('Cell.mat');
%     load('BW.mat');
%     Noise=noise_fft(Cell(1).F_rest,BW{1});
%     Cell(1).Noise=Noise;
    if numel(Cell)>1;
        Cell(1).Centrin= Centrin_fun();
        saveas(gcf,'plots/centrin.pdf');
        %%%% make images    
    else
        Cell.Centrin= Centrin_fun();
        saveas(gcf,'plots/centrin.pdf');
    end
        save('Cell.mat','Cell');
        close all;
end