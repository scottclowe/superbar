%BARE Plots bar graph with error bars included
%
% Inputs
%   X: location of bar. Empty for 1:size(Y,1)
%   Y: height of bars. Matrix for groups of bars [nGroups x nBarsPerGroup]
%   E: error on Y. Empty for no error bars. (def: empty)
%   width: (def: 0.8)
%   ori: orientation, 'v' or 'h' (def: 'v')
% Outputs:
%   hb: handles of bars
%   he: handles of error bars

function [hb,he] = bare(X, Y, E, width, ori)

% Input handling
if nargin<3
    E = [];
end
if nargin<4 || isempty(width)
    width = 0.8;
end
if nargin<5 || isempty(ori)
    ori = 'v';
end

% Default inputs
if isempty(X)
    X = 1:size(Y,1);
end
[X, Y, width] = bar2grouped(X, Y, width);

% Check if hold is already on
wasHeld = ishold(gca);
% If not, clear the axes and turn hold on
if ~wasHeld;
    cla;
    hold(gca,'on');
end;

if strncmp(ori,'h',1)
    % Horizontal bars & errors
    hb = barh(X, Y, width);
    if isempty(E)
        he = [];
    else
        he = ploterr(Y, X, E, []);
    end
else
    % Vertical bars & errors
    hb = bar(X, Y, width);
    if isempty(E)
        he = [];
    else
        he = ploterr(X, Y, [], E);
    end
end

% If hold was off, turn it off again
if ~wasHeld; hold(gca,'off'); end;
end
