%SUPERBAR Plots bar graph with error bars and each bar in a different colour
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

function [hhb, hhe] = superbar(X, Y, varargin)

% Input handling
if ischar(Y)
    % Deal with omitted X input
    varargin = [{Y}, varargin];
    Y = X;
    X = 1:size(Y, 1);
end
if isempty(X)
    X = 1:size(Y, 1);
end
% Use parser for the rest of the arguments
parser = inputParser;
addParameter(parser, 'E', []);
addParameter(parser, 'C', []);
addParameter(parser, 'CE', []);
addParameter(parser, 'width', [], @isnumeric);
addParameter(parser, 'orientation', 'v');
addParameter(parser, 'baseval', 0);
% addParameter(parser, 'theme', 'light');
parse(parser, varargin{:});

input = parser.Results;

if isempty(input.width)
    % Default with 0.8 of the smallest distance between bars
    input.width = 0.8 * min(diff(sort(X(:))));
end
if isempty(input.C)
    input.C = [.4 .4 .4];
end
if isempty(input.CE)
    input.CE = 0.7 * input.C;
end

if ~ismatrix(Y)
    error('Y should have no more than 2-dimensions (%d given)', ndims(Y));
end


if size(input.C,1)~=size(input.CE,1)
    error('Number of colours for error does not match colours for bars');
end

[X, Y, input.width] = bar2grouped(X, Y, input.width);

% Check if hold is already on
wasHeld = ishold(gca);
% If not, clear the axes and turn hold on
if ~wasHeld;
    cla;
    hold(gca, 'on');
end;

nBar = numel(Y);
hhb = nan(nBar, 1);
hhe = nan(nBar, 2);
for i=1:nBar
    % Check which colour to use
    k = mod(i-1, size(input.C, 1)) + 1;
    % Plot bar with error
    if isempty(input.E)
        hhb(i) = bare(X(i), Y(i), [], input.width, input.orientation);
    else
        [hhb(i),hhe(i,:)] = bare(X(i), Y(i), input.E(i), input.width, input.orientation);
        set(hhe(i,:), 'Color', input.CE(k,:));
    end
    % Colour it in correctly
    set(hhb(i), 'FaceColor', input.C(k,:), 'EdgeColor', 'none');
    if ~isempty(input.baseval)
        set(hhb(i), 'BaseValue', input.baseval);
    end
end

% If hold was off, turn it off again
if ~wasHeld; hold(gca, 'off'); end;
end


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

function [hb, he] = bare(X, Y, E, width, ori)

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
    hold(gca, 'on');
end;

if strncmp(ori, 'h', 1)
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
if ~wasHeld; hold(gca, 'off'); end;
end


%BAR2GROUPED Utility for BARE and BAREC.
function [XX, Y, varargout] = bar2grouped(X, Y, width)

% Parameters
group_width = 0.75;

nGroups = size(Y, 1);
nElePerGroup = size(Y, 2);

if size(X, 1)==1
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