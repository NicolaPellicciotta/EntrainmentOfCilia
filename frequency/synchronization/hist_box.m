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
function [Box] = hist_box(vv, ww, nbins,vinterval)
  minV  = min(vv);
  maxV  = max(vv);
  delta = (maxV(1)-minV(1))/nbins;
  g=[];
  values=[];
  if nargin<4
    vinterval = linspace(minV, maxV, nbins);
  end
  histw = zeros(nbins-1, 1);
  ehistw = zeros(nbins-1, 1);
  ones_vec=ones(size(vv));
  for i=1:(numel(vinterval)-1)
    ind = vv>= vinterval(i) & vv < vinterval(i+1);
    if ~isempty(ind)
      values=cat(2,values,ww(ind));
      g=cat(2,g,i*ones([numel(ww(ind)),1])');
    end
  end
    vinterval=vinterval(1:end-1)+delta/2;
    Box.vinterval=vinterval;
    Box.g=g;
    Box.values=values;
end