function [ r_bin,th_bin,cf,Ecf ] = phase_correlation_function_2d( ph,Ntheta,NR )

[height, width, times] = size(ph);
[X,Y]=meshgrid(1:width,1:height);

X=X(~isnan(ph(:,:,1)));
Y=Y(~isnan(ph(:,:,1)));

R = sqrt(((repmat(X,[1,numel(X)])) - (repmat(X,[1,numel(X)]))').^2 ...
       + ((repmat(Y,[1,numel(Y)])) - (repmat(Y,[1,numel(Y)]))').^2 );
TH= (atan2d( (repmat(Y,[1,numel(Y)])) - (repmat(Y,[1,numel(Y)]))',(repmat(X,[1,numel(X)])) - (repmat(X,[1,numel(X)]))'));    
TH(TH<0)=TH(TH<0)+180;   


Rl=(min(R(:))+1);
Ru=(max(R(:)))/2;
dr= floor((Ru-Rl)/NR);
if dr<1; dr=1; end

R_bin=Rl:dr:Ru;
dth=floor(180/Ntheta);
TH_bin=0:dth:180;
th_bin= (TH_bin(1:end-1)+TH_bin(2:end))/2;
r_bin= (R_bin(1:end-1)+R_bin(2:end))/2;

for tt=1:times
         pht=ph(:,:,tt);
         pht=pht(~isnan(ph(:,:,1)));
    
         MatP=repmat(pht,[1,numel(pht)]);
         CP= (exp(i*MatP).*conj(exp(i*MatP))');
        
         for ii=1:(numel(R_bin)-1)
             for jj=1:(numel(TH_bin)-1)
             
             slot= tril((R>=R_bin(ii) & R < (R_bin(ii+1)) ) & (TH>=TH_bin(jj) & TH < (TH_bin(jj+1)) ),-1) ;
             
             cf(ii,jj,tt)= angle(nanmean(CP(slot)));
             N(ii,jj,tt)=sum(slot(:));
             Ecf(ii,jj,tt)=  angle(sqrt( mean((CP(slot)-cf(ii,jj)).^2 )./N(ii,jj))) ;   %nansum(sqrt(MatF(slot)-F).^2)*nansum(sqrt(MatF(slot)'-F).^2);
             
  %           ecf(i)=  nanstd(CF (R>=R_bin(i) & R < (R_bin(i+1)) ))/sqrt(numel(CF (R>=R_bin(i) & R < (R_bin(i+1)) )));
 %            r_bin(ii,jj,tt)= (R_bin(ii)+R_bin(ii+1))/2; 
 %            th_bin(ii,jj,tt)= (TH_bin(jj)+TH_bin(jj+1))/2; 
             end
         end
         
end
% plot(dist_plot,av_corr)