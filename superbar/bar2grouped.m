%BAR2GROUPED Utility for BARE and BAREC.

function [XX, Y, varargout] = bar2grouped(X, Y, width)

% Parameters
group_width = 0.75;

nGroups = size(Y,1);
nElePerGroup = size(Y,2);

if size(X,1)==1
    X = X';
end

subwidth = group_width / nElePerGroup;
dX = (0:nElePerGroup-1)*subwidth - group_width/2 + subwidth/2;

XX = bsxfun(@plus, X, dX);

if nargin>=3
    if nGroups>1
        width = width * group_width / nElePerGroup;
    end
    varargout = {width};
end

end
