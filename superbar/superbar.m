%SUPERBAR  Plots bar graph with errors and p-values indicated
%   SUPERBAR(X,Y) draws the columns of the M-by-N matrix Y as M groups of
%   N vertical bars with location given by X and length given by Y. The
%   vector X should not have duplicate values.
%
%   SUPERBAR(Y) uses the default value of X=1:M. For vector inputs,
%   SUPERBAR(X,Y) or SUPERBAR(Y) draws LENGTH(Y) bars.
%
%   SUPERBAR(...,'E',E) plots the error, E, on each bar using SUPERERR. E
%   can be a matrix with the same number of elements as Y to specify
%   symmetric or single-directional errors, or an array with twice as many
%   elements (sized M-by-N-by-2 for instance) to specify asymmetric
%   errorbars. If E contains only 1 or 2 elements, the same symmetric or
%   asymmetric error bounds are used for every bar. Note that the ambiguous
%   case of plotting 2 bars with a/symmetric error bars should be
%   disambiguated by using a 1-by-1-by-2 array to apply the same asymmetric
%   error bounds to both bars and a 1-by-2 or 2-by-1 array for two
%   different but symmetric errors.
%
%   SUPERBAR(...,'P',P) with P the same size as Y adds stars above bars to
%   indicated whether the p-values in P are significant. Significance
%   thresholds can be set with the 'PStarThreshold' parameter (see below).
%
%   SUPERBAR(...,'P',P) with P a symmetric matrix sized (N*M)-by-(N*M) adds
%   lines between each bar indicating comparisons between them, with stars
%   above the lines which correspond to significant comparisons. Unmeasured
%   comparisons can be indicated with NaN values in P. This option is only
%   available for vertically oriented bars.
%
%   SUPERBAR(AX,...) plots into the axes with handle AX instead of GCA.
%
%   The inputs can be followed by parameter/value pairs to specify
%   additional properties, as follows.
%
%   General plot attributes
%   -----------------------
%       'Orientation' : Orientation of the bar plot. Set to 'v' for
%           vertical bars, or 'h' for horizontal bars. Default is 'v'. Note
%           that X is still the location of the bars and Y the length of
%           the bars, even if orientation 'h' is used.
%
%   Bar attributes
%   --------------
%       'BaseValue' : Base value from which bars begin. Default is 0.
%       'BarWidth' : Width of bars. Default is 80% of the minimum
%           separation between bars as specified in X.
%       'BarRelativeGroupWidth' : Relative width of bars when they are
%           grouped. Setting BarRelativeGroupWidth to 1 will have no spaces
%           between bars. Default is 0.8.
%       'BarFaceColor' : Color of the bars. Can be a colorspec string
%           (one of 'rgbymckw'), or an RGB array. In the case of an RGB
%           array, the size can be m-by-3 or m-by-n-by-3. If the input
%           contains fewer than M rows (or N columns in the later case),
%           the colors are repeated cyclically. Alternatively, can be set
%           to 'none' for transparent bar faces. Default is [.4, .4, .4].
%       'BarEdgeColor' : Color of the bars edges. For input options, see
%           'BarFaceColor'. Default is 'none'.
%
%   Errorbar attributes
%   -------------------
%       'E' : Errorbar magnitudes. Can be the same size as Y for specifying
%           symmetric or one-sided errorbars, or M-by-N-by-2 for asymmetric
%           errorbars. If E contains only a single value or two values, the
%           same symmetric or asymmetric errorbars are used for each bar.
%           Note that the ambiguous case with two bars and two errorbar
%           values should be disambiguated by ensuring E is 3-dimensional
%           when specifying asymmetric error bounds. If empty, no errorbars
%           are shown. Default is [].
%       'ErrorbarRelativeWidth' : Width of the errorbar caps, relative to
%           the bar width. Default is 0.7.
%       'ErrorbarColor' : Color of the errorbars. For input options, see
%           'BarFaceColor'. Default is the same as BarEdgeColor if it is
%           not 'none', otherwise 0.75 * BarFaceColor.
%       'ErrorbarStyle' : Shape of the errorbars to plot. Different
%           combinations allow plotting only stave, only caps, only
%           errorbars in a single direction, etc. Default is 'I', which has
%           a stave and cap in both directions always. See SUPERERR for a
%           list of possible errorbar styles. For instance, single-
%           directional errorbars can be acheived with the 'T' style.
%       'ErrorbarLineWidth' : LineWidth for errorbar lines. Default is 2.
%
%   P-value comparison attributes
%   -----------------------------
%       'P' : P-values. Can be either the same size as Y for specifying the
%           significance of each bar, or an (N*M)-by-(N*M) symmetric matrix
%           to indicate comparisons between each bar. If empty, no stars or
%           comparison lines are shown. Default is [].
%       'PStarThreshold' : Values which p-values must exceed (be smaller
%           than or equal to) to earn a star. Default is [0.05, 0.01,
%           0.001, 0.0001].
%       'PStarColor' : Color of the text for significance stars. Default is
%           [.2 .2 .2].
%       'PStarShowGT' : Whether to show a greater-than sign (>) for
%           p-values which are smaller than every value in PStarThreshold.
%           Default is true.
%       'PStarOffset' : Distance of the stars from the top of the errorbars
%           (or bars if no errorbars used). Default is 8% of the tallest
%           bar for single comparisons, or a quarter of PLineOffset for
%           paired comparisons.
%       'PLineColor' : Color of the lines indicating comparisons between
%           bars. Default is [.5 .5 .5].
%       'PLineWidth' : Width of the lines indicating comparisons between
%           bars. Default is 2.
%       'PLineOffset' : Vertical space between the comparison lines.
%           Default is 8% of the tallest bar.
%       'PLineSourceRelativeSpacing' : Maximum space between each line
%           coming from the top of a bar, relative to the width of the bar.
%           Default is 0.5.
%       'PLineSourceRelativeBreadth' : Maximum space which the lines coming
%           from each bar can collectively occupy, relative to the width of the
%           bar. Default is the same as ErrorbarRelativeWidth, if it is
%           non-zero, otherwise 0.8.
%       'PLineBacking' : Whether to pad p-value comparison lines by
%           plotting them on top of a backing line. Default is true.
%       'PLineBackingWidth' : Width of the line giving a backing color
%           behind each the comparison line. Default is 3 times PLineWidth,
%           so that the space on each side of the line is the same width as
%           the line itself.
%       'PLineBackingColor' : Color to use for the backing behind
%           comparison lines. Default is the axes background color.
%
%   [HB, HE, HPT, HPL, HPB] = SUPERBAR(...) returns handles to the
%   generated graphics objects. HB contains handles to the bars themselves,
%   in a matrix whose size matches that of Y. HE contains handles to the
%   errorbars, in a matrix whose size matches that of Y. HPT contains
%   handles to the text showing p-value siginficance levels with stars. HPL
%   contains handles to the comparison lines between bars. HPB contains
%   handles to the background behind the comparison lines.
%
%   Note that unlike BAR and BARH, bars plotted with SUPERBAR are always
%   grouped and never stacked.
%
%   See also SUPERERR, BAR, BARH.

function varargout = superbar(X, Y, varargin)

% Check number of inputs is okay
narginchk(1, Inf);

% Extend the reach of varargin
if nargin>=2
    varargin = [{Y}, varargin];
end
varargin = [{X}, varargin];

% Strip out axes input if it is there
[ax, varargin, nargs] = axescheck(varargin{:});
% Otherwise, default with the current axes
if isempty(ax)
    ax = gca;
end
% Check number of inputs is still okay
if nargs<1
    error('Must provide at least 1 input argument, in addition to axes.');
end

% Input handling for X and Y
if nargs==1 || ischar(varargin{2})
    % Deal with omitted X input
    Y = varargin{1};
    X = 1:size(Y, 1);
    % Drop the argument
    varargin = varargin(2:end);
else
    % Take X and Y out of varargin
    X = varargin{1};
    Y = varargin{2};
    % Drop these arguments
    varargin = varargin(3:end);
end
if ~ismatrix(Y)
    error('Y should have no more than 2-dimensions (%d given).', ndims(Y));
end
if size(Y, 1)==1
    Y = Y';
end
if isempty(X)
    X = (1:size(Y, 1))';
end
if size(X, 1)~=size(Y, 1) && size(X, 2)==size(Y, 1)
    X = X';
end
if size(X, 1)~=size(Y, 1)
    error('X and Y must be the same size in dimension 1 (%d and %d given)', ...
        size(X, 1), size(Y, 1));
end

% Use parser for the rest of the arguments
parser = inputParser;
% Plot attributes
addParameter(parser, 'Orientation', 'v', ...
    @ischar);
addParameter(parser, 'BaseValue', 0, ...
    @(t) (isscalar(t)) && isnumeric(t));
% Bar attributes
addParameter(parser, 'BarWidth', [], ...
    @(t) (isempty(t) || isscalar(t)) && isnumeric(t));
addParameter(parser, 'BarRelativeGroupWidth', 0.8, ...
    @(t) (isscalar(t)) && isnumeric(t));
addParameter(parser, 'BarFaceColor', [.4, .4, .4]);
addParameter(parser, 'BarEdgeColor', 'none');
% Errorbar attributes
addParameter(parser, 'E', []);
addParameter(parser, 'ErrorbarRelativeWidth', 0.7, ...
    @(t) (isscalar(t)) && isnumeric(t));
addParameter(parser, 'ErrorbarColor', []);
addParameter(parser, 'ErrorbarStyle', 'I', ...
    @ischar);
addParameter(parser, 'ErrorbarLineWidth', 2, ...
    @(t) (isscalar(t)) && isnumeric(t));
% P-value attributes
addParameter(parser, 'P', []);
addParameter(parser, 'PStarThreshold', [0.05, 0.01, 0.001, 0.0001], ...
    @isnumeric);
addParameter(parser, 'PStarColor', [.2 .2 .2]);
addParameter(parser, 'PStarShowGT', true, ...
    @isscalar);
addParameter(parser, 'PStarOffset', [], ...
    @(t) (isempty(t) || isscalar(t)) && isnumeric(t));
addParameter(parser, 'PLineColor', [.5 .5 .5]);
addParameter(parser, 'PLineWidth', 2, ...
    @(t) (isscalar(t)) && isnumeric(t));
addParameter(parser, 'PLineOffset', [], ...
    @(t) (isempty(t) || isscalar(t)) && isnumeric(t));
addParameter(parser, 'PLineSourceRelativeSpacing', 0.5, ...
    @(t) (isscalar(t)) && isnumeric(t));
addParameter(parser, 'PLineSourceRelativeBreadth', [], ...
    @(t) (isempty(t) || isscalar(t)) && isnumeric(t));
addParameter(parser, 'PLineBacking', true, ...
    @isscalar);
addParameter(parser, 'PLineBackingWidth', [], ...
    @(t) (isempty(t) || isscalar(t)) && isnumeric(t));
addParameter(parser, 'PLineBackingColor', []);
% Parse the arguments
parse(parser, varargin{:});
input = parser.Results;

% Default input arguments which inherit values from others
% Bar defaults
if isempty(input.BarWidth)
    % Default with 0.8 of the smallest distance between bars
    input.BarWidth = 0.8 * min(diff(sort(X(:))));
end
% Errorbar defaults
if isempty(input.ErrorbarColor)
    % Try taking the colour from the errorbar style
    COLOR_SPECS = 'rgbwcmyk';
    isColorSpec = ismember(COLOR_SPECS, input.ErrorbarStyle);
    if any(isColorSpec)
        idx = find(isColorSpec, 1, 'last');
        input.ErrorbarColor = input.ErrorbarStyle(idx);
    elseif ~isequal(input.BarEdgeColor, 'none')
        % Try taking the color from the bar edge
        input.ErrorbarColor = input.BarEdgeColor;
    elseif ~isequal(input.BarFaceColor, 'none')
        % Try taking the color from the bar face
        color = input.BarFaceColor;
        if ischar(color)
            % Convert string into RGB colour
            color = colorspec2rgb(color);
        end
        % Make the bar colour darker
        input.ErrorbarColor = 0.7 * color;
    else
        % Well if you want everything transparent you can have it, I
        % guess, though maybe this should be an error instead
        input.ErrorbarColor = 'none';
    end
end
% P-value defaults
if isempty(input.PLineOffset)
    % Base the offset on the maximum of the bars
    input.PLineOffset = 0.075 * max(abs(Y(:)));
end
if isempty(input.PStarOffset)
    if numel(input.P)==numel(X) || numel(input.P)==numel(Y)
        % If we're just showing the stars and no lines, base the offset on
        % the maximum of the bars
        input.PStarOffset = 0.08 * max(abs(Y(:)));
    else
        % If we're showing comparison lines, make the stars be a little
        % above the lines so its clear to which they belong
        input.PStarOffset = input.PLineOffset / 4;
    end
end
if isempty(input.PLineSourceRelativeBreadth)
    if ~isempty(input.ErrorbarRelativeWidth) && input.ErrorbarRelativeWidth>0
        % The breadth of the space the lines come from should be the same
        % as the width of the errorbars, if possible
        input.PLineSourceRelativeBreadth = input.ErrorbarRelativeWidth;
    else
        % Otherwise base it on the bar width
        input.PLineSourceRelativeBreadth = 0.8;
    end
end
if isempty(input.PLineBackingColor)
    % Color to pad lines with
    input.PLineBackingColor = get(ax, 'Color');
end
if isempty(input.PLineBackingWidth)
    input.PLineBackingWidth = 3 * input.PLineWidth;
end

% Split up bars which are composed of groups
if size(X, 2)==1
    [X, Y, input.BarWidth] = bar2grouped(X, Y, input.BarWidth, ...
        input.BarRelativeGroupWidth);
end
% Check size of X and Y match
assert(isequal(size(X), size(Y)), ...
    'Sizes of X and Y must match. Sizes were %s and %s.', ...
    mat2str(size(X)), mat2str(size(Y)));

% Fix relative widths
errorbarWidth = input.ErrorbarRelativeWidth * input.BarWidth;
PLineSourceBreadth = input.PLineSourceRelativeBreadth * input.BarWidth;
PLineSourceSpacing = input.PLineSourceRelativeSpacing * input.BarWidth;

% Fix shape of E
if numel(input.E)==1
    input.E = repmat(input.E, numel(Y), 1);
elseif numel(input.E)==2*numel(Y)
    input.E = reshape(input.E, numel(Y), 2);
elseif numel(input.E)==2 && (numel(Y)~=2 || ~ismatrix(input.E))
    input.E = repmat(input.E(:)', numel(Y), 1);
elseif numel(input.E)==numel(Y)
    input.E = input.E(:);
elseif ~isempty(input.E)
    error(...
        ['E input must contain either the same number of values as Y' ...
        ' (for symmetric errorbars), or twice as many values (for' ...
        ' asymmetric errorbars), or a single value/pair of values (to' ...
        ' use the same error for every bar).']);
end

% Extend colors to be per bar
function C = extend_colors(C)
    siz = size(Y);
    siz_ = size(C);
    if ~ischar(C)
        assert(length(siz_)<=3, 'Too many dimensions for C.');
        assert(siz_(end)==3, 'Must be RGB color in C with 3 channels.');
    end
    if length(siz_)==2
        C = permute(C, [1, 3, 2]);
    end
    siz_ = size(C);
    C = repmat(C, ceil(siz(1) / siz_(1)), ceil(siz(2) / siz_(2)));
    C = C(1:siz(1), 1:siz(2), :);
end
input.BarFaceColor = extend_colors(input.BarFaceColor);
input.BarEdgeColor = extend_colors(input.BarEdgeColor);
input.ErrorbarColor = extend_colors(input.ErrorbarColor);

% Check if hold is already on
wasHeld = ishold(ax);
% If not, clear the axes and turn hold on
if ~wasHeld;
    cla(ax);
    hold(ax, 'on');
end;

nBar = numel(Y);
hb = nan(nBar, 1);
for iBar=1:nBar
    % Get indices to tell which row and column to take colour from
    [i, j] = ind2sub(size(Y), iBar);
    % Plot bar
    if strncmpi(input.Orientation, 'h', 1)
        hb(iBar) = barh(X(iBar), Y(iBar), input.BarWidth);
    else
        hb(iBar) = bar(X(iBar), Y(iBar), input.BarWidth);
    end
    % Colour it in correctly
    set(hb(iBar), ...
        'FaceColor', input.BarFaceColor(i,j,:), ...
        'EdgeColor', input.BarEdgeColor(i,j,:), ...
        'BaseValue', input.BaseValue);
end
% Add errorbars
if isempty(input.E)
    he = [];
elseif strncmpi(input.Orientation, 'h', 1)
    % Horizontal errorbars
    he = supererr(Y, X, input.E, [], input.ErrorbarStyle, ...
        errorbarWidth, ...
        'Color', input.ErrorbarColor, ...
        'LineWidth', input.ErrorbarLineWidth);
    he = reshape(he(:, 1), size(Y));
else
    % Vertical errorbars
    he = supererr(X, Y, [], input.E, input.ErrorbarStyle, ...
        errorbarWidth, ...
        'Color', input.ErrorbarColor, ...
        'LineWidth', input.ErrorbarLineWidth);
    he = reshape(he(:, 2), size(Y));
end
% Add p-values
if isempty(input.P)
    % Do nothing
    hpt = [];
    hpl = [];
    hpb = [];
elseif numel(input.P)==numel(X)
    % Add stars above bars
    hpt = plot_p_values_single(ax, X, Y, input.E, input.P, ...
        input.Orientation, input.BaseValue, input.PStarThreshold, ...
        input.PStarOffset, input.PStarShowGT, input.PStarColor);
    hpl = [];
    hpb = [];
elseif numel(input.P)==numel(X)^2
    % Add lines and stars between pairs of bars
    [hpt, hpl, hpb] = plot_p_values_pairs(ax, X, Y, input.E, input.P, ...
        input.PStarThreshold, input.PLineOffset, input.PStarOffset, ...
        input.PStarShowGT, PLineSourceSpacing, ...
        PLineSourceBreadth, input.PLineBacking, ...
        {'Color', input.PLineColor, ...
         'LineWidth', input.PLineWidth}, ...
        {'Color', input.PLineBackingColor, ...
         'LineWidth' input.PLineBackingWidth}, ...
        {'Color', input.PStarColor});
else
    error('Bad number of P-values');
end

% If hold was off, turn it off again
if ~wasHeld; hold(ax, 'off'); end;

if nargout==0
    varargout = {};
else
    varargout = {hb, he, hpt, hpl, hpb};
end

end


%colorspec2rgb
%   Convert a color string to an RGB value.
%          b     blue
%          g     green
%          r     red
%          c     cyan
%          m     magenta
%          y     yellow
%          k     black
%          w     white
function color = colorspec2rgb(color)

% Define the lookup table
rgb = [1 0 0; 0 1 0; 0 0 1; 1 1 1; 0 1 1; 1 0 1; 1 1 0; 0 0 0];
colspec = 'rgbwcmyk';

idx = find(colspec==color(1));
if isempty(idx)
    error('colorstr2rgb:InvalidColorString', 'Unknown color string.');
end

if idx~=3 || length(color)==1,
    color = rgb(idx, :);
elseif length(color)>2,
    if strcmpi(color(1:3), 'bla')
        color = [0 0 0];
    elseif strcmpi(color(1:3), 'blu')
        color = [0 0 1];
    else
        error('colorstr2rgb:UnknownColorString', 'Unknown color string.');
    end
end

end


%bar2grouped
%   Split vector X to position all bars in a group into appropriate places.
function [X, Y, width] = bar2grouped(X, Y, width, group_width)

% Parameters
if nargin<4
    group_width = 0.75;
end

nElePerGroup = size(Y, 2);

if nElePerGroup==1
    % No need to do anything to X as the groups only contain one element
    return;
end

if size(X, 1)==1
    % Transpose X if necessary
    X = X';
end
if ~ismatrix(X) || size(X, 2)>1
    error('X must be a column vector.')
end

% Compute the offset for each bar, such that they are centred correctly
dX = width / nElePerGroup * ((0:nElePerGroup-1) - (nElePerGroup-1)/2);
% Apply the offset to each bar in X
X = bsxfun(@plus, X, dX);
% Reduce width of bars so they only take up group_width as much space, and
% divide what there is evenly between the bars per group
width = width * group_width / nElePerGroup;

end


%plot_p_values_single
%   Plot stars above bars to indicate which are statistically significant.
%   Can be used with bars in either horizontal or vertical direction.
function h = plot_p_values_single(ax, X, Y, E, P, orientation, baseval, ...
    p_threshold, offset, show_gt, text_color)

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
    h(i) = text(ax, x, y, str, ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle', ...
        'Color', text_color);
end

end


%plot_p_values_pairs
%   Plot lines and stars to indicate pairwise comparisons and whether they
%   are significant. Only works for error bars in the Y-direction.
function [ht, hl, hbl] = plot_p_values_pairs(ax, X, Y, E, P, p_threshold, ...
    offset, star_offset, show_gt, max_dx_single, max_dx_full, ...
    pad_lines, line_args, pad_args, text_args)

% Validate inputs
N = numel(X);
assert(numel(Y)==N, 'Number of datapoints mismatch {X,Y}.');
assert(numel(E)==N, 'Number of datapoints mismatch {X,E}.');
assert(numel(P)==N^2, 'Number of datapoints mismatch {X,P}.');

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
    % Check we're not failing terribly
    if isnan(P(i,j))
        error('This shouldnt be NaN!');
    end
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
        hbl(iPair) = line(ax, xx, yy, pad_args{:});
    end
    hl(iPair) = line(ax, xx, yy, line_args{:});
    % Check how many stars to put in the text
    num_stars = sum(P(ISi(iPair), ISj(iPair)) <= p_threshold);
    str = repmat('*', 1, num_stars);
    % Check whether to include a > sign too
    if show_gt && all(P(ISi(iPair), ISj(iPair)) < p_threshold)
        str = ['>' str];
    end
    % Add the text for the stars, slightly above the middle of the line
    ht(iPair) = text(ax, mean(xx), yy(2) + star_offset, str, ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle', ...
        text_args{:});
end

end