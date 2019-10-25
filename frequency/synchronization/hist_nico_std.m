% HISTWC  Weighted histogram count given number of bins
%
% This function generates a vector of cumulative weights for data
% histogram. Equal number of bins will be considered using minimum and 
% maximum values of the data. Weights will be summed in the given bin.
%
% Usage: [histw, vinterval] = histwc(vv, ww, nbins)
%
% Arguments:
%       vv    - values as a vector
%       ww    - weights as a vector
%       nbins - number of bins
%
% Returns:
%       histw     - weighted histogram
%       vinterval - intervals used
%       
%
%
% See also: HISTC, HISTWCV
% Author:
% mehmet.suzen physics org
% BSD License
% July 2013
function [histw,ehistw,mean_bins] = hist_nico_std(vv, ww, nbins,vinterval)
  if  isempty(vinterval) & ~isempty(nbins)
        minV  = min(vv);
        maxV  = max(vv);
        vinterval = linspace(minV, maxV, nbins+1);
  elseif isempty(nbins) & ~isempty(vinterval)
        nbins=numel(vinterval)-1;
  elseif isempty(nbins) & isempty(vinterval)
        nbins=floor(numel(vv)/10); vinterval = linspace(minV, maxV, nbins+1);
  end
  
  
  [~,~,ind]=histcounts(vv,vinterval);
%  delta = (maxV(1)-minV(1))/nbins;
   good= ind~=0;
   ind=ind(good);ww=ww(good);ww=ww';ind=ind'
   histw=accumarray(ind,ww,[],@nanmedian);
   ehistw=accumarray(ind,ww,[],@nanstd);%./sqrt((accumarray(ind,ww,[],@numel)));
   N=(accumarray(ind,ww,[],@numel));
   mean_bins=(vinterval(1:end-1)+vinterval(2:end))/2;

end