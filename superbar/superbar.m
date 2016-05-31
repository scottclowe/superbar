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

function [hhb, hhe, hht, hhl] = superbar(X, Y, varargin)

% Check number of inputs is okay
narginchk(1, Inf);

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
addParameter(parser, 'P', []);
addParameter(parser, 'width', [], @isnumeric);
addParameter(parser, 'orientation', 'v');
addParameter(parser, 'baseval', 0);
addParameter(parser, 'p_threshold', [0.05, 0.01, 0.001, 0.0001]);
addParameter(parser, 'show_gt', true);
addParameter(parser, 'p_offset', []);
addParameter(parser, 'max_dx_single', []);
addParameter(parser, 'max_dx_full', []);
addParameter(parser, 'p_line_color', [.5 .5 .5]);
addParameter(parser, 'pad_lines', true);
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
if isempty(input.p_offset)
    input.p_offset = 0.075 * max(Y(:));
end

if ~ismatrix(Y)
    error('Y should have no more than 2-dimensions (%d given)', ndims(Y));
end


if size(input.C,1)~=size(input.CE,1)
    error('Number of colours for error does not match colours for bars');
end

[X, Y, input.width] = bar2grouped(X, Y, input.width);

if isempty(input.max_dx_single)
    input.max_dx_single = input.width * 0.25;
end
if isempty(input.max_dx_full)
    input.max_dx_full = input.width * 0.75;
end

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
        [hhb(i), hhe(i,:)] = bare(X(i), Y(i), input.E(i), input.width, ...
            input.orientation);
        set(hhe(i,:), 'Color', input.CE(k,:));
    end
    % Colour it in correctly
    set(hhb(i), 'FaceColor', input.C(k,:), 'EdgeColor', 'none');
    if ~isempty(input.baseval)
        set(hhb(i), 'BaseValue', input.baseval);
    end
end
if isempty(input.P)
    % Do nothing
    hht = [];
    hhl = [];
elseif numel(input.P)==numel(X)
    % Add stars above bars
    hht = plot_p_values_single(X, Y, input.E, input.P, input.p_threshold, ...
        input.p_offset, input.show_gt, input.orientation);
    hhl= [];
elseif numel(input.P)==numel(X)^2
    % Add lines and stars between pairs of bars
    [hhl, hht] = plot_p_values_pairs(X, Y, input.E, input.P, ...
        input.p_threshold, input.p_offset, input.show_gt, ...
        input.max_dx_single, input.max_dx_full, input.pad_lines, ...
        'Color', input.p_line_color);
else
    error('Bad number of P-values');
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
        he = supererr(Y, X, E, [], 'I', width/2);
        he = he(:, 1);
    end
else
    % Vertical bars & errors
    hb = bar(X, Y, width);
    if isempty(E)
        he = [];
    else
        he = supererr(X, Y, [], E, 'I', width/2);
        he = he(:, 2);
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


%plot_p_values_single
%   Plot stars above bars to indicate which are statistically significant.
%   Can be used with bars in either horizontal or vertical direction.
function h = plot_p_values_single(X, Y, E, P, p_threshold, offset, ...
    show_gt, orientation, baseval)

% Default inputs
if nargin < 9
    baseval = 0;
end
if nargin < 8
    orientation = 'v';
end
if nargin < 7
    show_gt = true;
end

% Validate inputs
assert(numel(X)==numel(Y), 'Number of datapoints mismatch {X,Y}.');
assert(numel(X)==numel(E), 'Number of datapoints mismatch {X,E}.');
assert(numel(X)==numel(P), 'Number of datapoints mismatch {X,P}.');
assert(all(E(:) >= 0), 'Error must be a non-negative value.');
assert(offset > 0, 'Offset must be a positive value.');

% Loop over every bar
h = nan(size(X));
for i=1:numel(X)
    % Check how many stars to put in the text
    num_stars = sum(P(i) <= p_threshold);
    str = repmat('*', 1, num_stars);
    % Check whether to include a > sign too
    if show_gt && all(P(i) < p_threshold)
        str = ['>' str];
    end
    % Work out where to put the text
    x = X(i);
    y = Y(i);
    if strncmpi(orientation, 'h', 1)
        if x >= baseval
            x = x + E(i) + offset;
        else
            x = x - E(i) - offset;
        end
    else
        if y >= baseval
            y = y + E(i) + offset;
        else
            y = y - E(i) - offset;
        end
    end
    % Add the text for the stars
    h(i) = text(x, y, str, ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle');
end

end


%plot_p_values_pairs
%   Plot lines and stars to indicate pairwise comparisons and whether they
%   are significant. Only works for error bars in the Y-direction.
function [hl, ht] = plot_p_values_pairs(X, Y, E, P, p_threshold, offset, ...
    show_gt, max_dx_single, max_dx_full, pad_lines, varargin)

% Validate inputs
N = numel(X);
assert(numel(Y)==N, 'Number of datapoints mismatch {X,Y}.');
assert(numel(E)==N, 'Number of datapoints mismatch {X,E}.');
assert(numel(P)==N^2, 'Number of datapoints mismatch {X,P}.');

% Color to pad lines with
bg_color = get(gca, 'Color');

% Turn into vectors
X = X(:);
Y = Y(:);
E = E(:);
P = reshape(P, N, N);
assert(all(all(P==P' | isnan(P))), 'P must be symmetric between pairs');

% Sort by bar location
[X, IX] = sort(X);
Y = Y(IX);
E = E(IX);
P = P(IX, IX);

% Ensure P is symmetric
P = max(P, P');
% Remove lower triangle
P(logical(tril(ones(size(P))))) = NaN;

% Find the max of each pair of bars
pair_max_y = max(repmat(Y + E, 1, N), repmat(Y' + E', N, 1));
% Find the distance between the bars
pair_distance = abs(repmat(X, 1, N) - repmat(X', N, 1));
% Remove pairs which are not measured
li = isnan(P);
pair_max_y(li) = NaN;
pair_distance(li) = NaN;

% Sort by maximum value, smallest first
[~, I1] = sort(pair_max_y(:));
pair_distance_tmp = pair_distance(I1);
% Sort by pair distance
[~, I2] = sort(pair_distance_tmp);
% Combine the two mappings into a single indexing step
IS = I1(I2);
% Now we have primary sort by pair_distance and secondary sort by max value
[ISi, ISj] = ind2sub(size(pair_distance), IS);

% For each bar, check how many lines there will be
num_comp_per_bar = sum(~isnan(max(P, P')), 2);
dX_list = nan(size(P));
for i=1:numel(X)
    dX_list(1:num_comp_per_bar(i), i) = ...
        (0:(num_comp_per_bar(i)-1)) - (num_comp_per_bar(i)-1) / 2;
end
dX_each = min(max_dx_single, max_dx_full / max(num_comp_per_bar));
dX_list = dX_list * dX_each;

% Minimum value for lines over each bar
YEO = Y + E + offset / 2;
current_height = repmat(YEO(:)', N, 1) + offset / 2;

% Loop over every pair with a measurement
num_comparisons = sum(~isnan(P(:)));
hl = nan(size(P));
ht = nan(size(P));
hbl = nan(size(P));
coords = nan(4, 2, num_comparisons);
for iPair=1:num_comparisons
    % Get index of left and right pairs
    i = min(ISi(iPair), ISj(iPair));
    j = max(ISi(iPair), ISj(iPair));
    % Check which bar origin point we're up to
    il = find(~isnan(dX_list(:, i)), 1, 'last');
    jl = find(~isnan(dX_list(:, j)), 1, 'first');
    % Offset the X value to get the non-intersecting origin point
    xi = X(i) + dX_list(il, i);
    xj = X(j) + dX_list(jl, j);
    % Clear these origin points so they aren't reused
    dX_list(il, i) = NaN;
    dX_list(jl, j) = NaN;
    xx = [xi, xi, xj, xj];
    % Work out how high the line must be
    yi = YEO(i);
    yj = YEO(j);
    % It must be higher than all intermediate lines; check which these are
    intermediate_index = (il + N*(i-1)) : (jl + N*(j-1));
    % Also offset so we are higher than these lines
    y_ = max(current_height(intermediate_index)) + offset;
    yy = [yi, y_, y_, yj];
    % Update intermediates so we know the new hight above them
    current_height(intermediate_index) = y_;
    % Save the co-ordinates to plot later
    coords(:, 1, iPair) = xx;
    coords(:, 2, iPair) = yy;
end
for iPair=num_comparisons:-1:1
    % Get co-ordinates back again
    xx = coords(:, 1, iPair);
    yy = coords(:, 2, iPair);
    % Draw the line
    if pad_lines
        hbl(iPair) = line(xx, yy, 'Color', bg_color, 'LineWidth', 3);
    end
    hl(iPair) = line(xx, yy, varargin{:});
    % Check how many stars to put in the text
    num_stars = sum(P(iPair) < p_threshold);
    str = repmat('*', 1, num_stars);
    % Check whether to include a > sign too
    if show_gt && all(P(i) < p_threshold)
        str = ['>' str];
    end
    % Add the text for the stars, slightly above the middle of the line
    ht(iPair) = text(mean(xx), yy(2) + offset/4, str, ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle');
end

end