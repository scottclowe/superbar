%BAREC Plots bar graph with error bars and each bar in a different colour
%
% Inputs
%   X: location of bar. Empty for 1:size(Y,1)
%   Y: height of bars. Matrix for groups of bars [nGroups x nBarsPerGroup]
%   E: error on Y. Empty for no error bars. (def: empty)
%   C: colour of bars (def: dark gray) [nColours x 3]
%   CE: colour of error bars (def: 70% of C)
%   width: (def: 0.8)
%   ori: orientation, 'v' or 'h' (def: 'v')
%   baseval: base value for bars (you can't set them individually from the
%            handles if you plotted more bars for some reason)
% Outputs:
%   hb: handles of bars
%   he: handles of error bars

function [hhb,hhe] = barec(X, Y, E, C, CE, width, ori, baseval)

% Input handling
if isempty(X)
    X = 1:size(Y,1);
end
if nargin<3 || isempty(E)
    E = [];
end
if nargin<4 || isempty(C)
    C = [.4 .4 .4];
end
if nargin<5 || isempty(CE)
    CE = 0.7*C;
end
if nargin<6 || isempty(width)
    width = 0.8;
end
if nargin<7 || isempty(ori)
    ori = 'v';
end
if nargin<8 || isempty(baseval)
    baseval = [];
end

if size(C,1)~=size(CE,1)
    error('Number of colours for error does not match colours for bars');
end

[X, Y, width] = bar2grouped(X, Y, width);

% Check if hold is already on
wasHeld = ishold(gca);
% If not, clear the axes and turn hold on
if ~wasHeld;
    cla;
    hold(gca,'on');
end;

nBar = numel(Y);
hhb = nan(nBar,1);
hhe = nan(nBar,2);
for i=1:nBar
    % Check which colour to use
    k = mod(i-1,size(C,1))+1;
    % Plot bar with error
    if isempty(E)
        hhb(i) = bare(X(i), Y(i), [], width, ori);
    else
        [hhb(i),hhe(i,:)] = bare(X(i), Y(i), E(i), width, ori);
        set(hhe(i,:),'Color',CE(k,:));
    end
    % Colour it in correctly
    set(hhb(i),'FaceColor',C(k,:),'EdgeColor','none');
    if ~isempty(baseval)
        set(hhb{i},'BaseValue',baseval);
    end
end

% If hold was off, turn it off again
if ~wasHeld; hold(gca,'off'); end;
end
