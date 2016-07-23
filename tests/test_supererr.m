% Run and advance each cell and output is generated appropriately
clear all;
figure;

%% Plotting single X errorbar
clear all;
clf;
supererr(1, 1, 0.2);
xlim([0 2]);
ylim([0 2]);
title('Plotting single X errorbar');

%% Plotting single X errorbar, empty Y
clear all;
clf;
supererr(1, 1, 0.2, []);
xlim([0 2]);
ylim([0 2]);
title('Plotting single X errorbar, empty Y');

%% Plotting single Y errorbar
clear all;
clf;
supererr(1, 1, [], 0.2);
xlim([0 2]);
ylim([0 2]);
title('Plotting single Y errorbar');

%% Plot three bars, one at a time
clear all;
clf;
supererr(1, 1, 0.2, 0.4);
hold on;
supererr(2, 3, 0.1, 0.2);
supererr(3, 2, 0.3, 0.6);
title('Plot three errorbars, one at a time');

%% Plot three errorbars, at once, same errors
clear all;
clf;
supererr([1 2 3], [1 3 2], 0.2, 0.4);
title('Plot three errorbars, at once, same errors');

%% Plot 3 errorbars, different upper and lower errors
clear all;
clf;
supererr([1 2 3], [1 3 2], [0.1 0.3], [0.2 0.4]);
title('Plot 3 errorbars, different upper and lower errors');

%% Plot 3 errorbars, at once, different symmetric errors
clear all;
clf;
supererr([1 2 3], [1 3 2], [0.1; 0.3; 0.2], [0.3; 0.9; 0.6]);
title('Plot 3 errorbars, at once, different symmetric errors');

%% Plot 3 errorbars, at once, different asymmetric errors
clear all;
clf;
supererr([1 2 3], [1 3 2], [0.1 0.2; 0.3 0.6; 0.4 0.2], ...
    [0.2 0.4; 0.9 0.6; 0.4 0.8]);
title('Plot 3 errorbars, at once, different asymmetric errors');

%% Plot matrix input
clear all;
clf;
supererr([1 4; 3 2], [7 2 3 6]', [0.1; 0.3; 0.4; 0.2], ...
    [0.4 0.2 0.6 0.8]');
title('Plot matrix input');

%% Plot 2 errorbars, different upper and lower errors
clear all;
clf;
supererr([1 2], [1 3], [0.1 0.3], [0.2 0.4]);
title('Plot 2 errorbars, different upper and lower errors');

%% Plot 2 errorbars, at once, different symmetric errors
clear all;
clf;
supererr([1 2], [1 3], [0.1; 0.3], [0.3; 0.9]);
title('Plot 2 errorbars, at once, different symmetric errors');

%% Plot into non-current axes
clear all;
clf;
ax1 = subplot(1, 2, 1);
ax2 = subplot(1, 2, 2);
supererr(ax1, [1 2 3], [1 3 2], [0.1 0.2; 0.3 0.6; 0.4 0.2], ...
    [0.2 0.4; 0.9 0.6; 0.4 0.8]);
text(.5, .5, 'current axes', 'HorizontalAlignment', 'center');
title(ax1, 'Plot should be here');
title(ax2, 'Not here');

%% Red style string
clear all;
clf;
supererr([1 2 3], [1 3 2], [0.1 0.2; 0.3 0.6; 0.4 0.2], ...
    [0.2 0.4; 0.9 0.6; 0.4 0.8], 'r');
title('Red style string');

%% One-directional T style
clear all;
clf;
supererr([1 2 3], [1 3 2], [0.2; 0.6; 0.2], ...
    [0.4; 0.6; 0.8], 'T');
title('One-directional T style');

%% Two-directional T style
clear all;
clf;
supererr([1 2 3], [1 3 2], [0.1 0.2; 0.3 0.6; 0.4 0.2], ...
    [0.2 0.4; 0.9 0.6; 0.4 0.8], 'T');
title('Two-directional T style');

%% Stave, no caps, magenta
clear all;
clf;
supererr([1 2 3], [1 3 2], [0.1 0.2; 0.3 0.6; 0.4 0.2], ...
    [0.2 0.4; 0.9 0.6; 0.4 0.8], '|m');
title('Stave, no caps, magenta');

%% Blue, caps only
clear all;
clf;
supererr([1 2 3], [1 3 2], [0.1 0.2; 0.3 0.6; 0.4 0.2], ...
    [0.2 0.4; 0.9 0.6; 0.4 0.8], 'b=');
title('Blue, caps only');

%% Different X and Y styles
clear all;
clf;
supererr([1 2 3], [1 3 2], [0.1 0.2; 0.3 0.6; 0.4 0.2], ...
    [0.2 0.4; 0.9 0.6; 0.4 0.8], 'gI', 'k_');
title('Different X and Y styles');

%% Set errorbar width
clear all;
clf;
supererr([1 2 3], [1 3 2], [0.1 0.2; 0.3 0.6; 0.4 0.2], ...
    [0.2 0.4; 0.9 0.6; 0.4 0.8], 0.2);
title('Set errorbar width');

%% Set errorbar width after styles
clear all;
clf;
supererr([1 2 3], [1 3 2], [0.1 0.2; 0.3 0.6; 0.4 0.2], ...
    [0.2 0.4; 0.9 0.6; 0.4 0.8], 'gI', 'k_', 0.2);
title('Set errorbar width after styles');

%% Linewidth
clear all;
clf;
supererr([1 2 3], [1 3 2], [0.1 0.2; 0.3 0.6; 0.4 0.2], ...
    [0.2 0.4; 0.9 0.6; 0.4 0.8], 'gI', 'b|', 'LineWidth', 4);
title('Linewidth');

%% Color input colorstr single
clear all;
clf;
supererr([1 2 3], [1 3 2], [0.1 0.2; 0.3 0.6; 0.4 0.2], ...
    [0.2 0.4; 0.9 0.6; 0.4 0.8], '|', 0.2, 'Color', 'r');
title('Color input colorstr single ''r''');

%% Color input colorstr column vector
clear all;
clf;
supererr([1 2 3], [1 3 2], [0.1 0.2; 0.3 0.6; 0.4 0.2], ...
    [0.2 0.4; 0.9 0.6; 0.4 0.8], '|', 0.2, 'Color', ['rb']');
title('Color input colorstr column vector ''rb''');

%% Color input colorstr row vector
clear all;
clf;
supererr([1 2 3], [1 3 2], [0.1 0.2; 0.3 0.6; 0.4 0.2], ...
    [0.2 0.4; 0.9 0.6; 0.4 0.8], '|', 0.2, 'Color', 'rb');
title('Color input colorstr row vector ''rb''');

%% Color input RGB single
clear all;
clf;
supererr([1 2 3], [1 3 2], [0.1 0.2; 0.3 0.6; 0.4 0.2], ...
    [0.2 0.4; 0.9 0.6; 0.4 0.8], '|', 0.2, 'Color', [.6 .6 .1]);
title('Color input RGB single');

%% Color input RGB multi
clear all;
clf;
supererr([1 2 3], [1 3 2], [0.1 0.2; 0.3 0.6; 0.4 0.2], ...
    [0.2 0.4; 0.9 0.6; 0.4 0.8], 'Color', [.6 .6 .1; .1 .6 .6]);
title('Color input RGB multi');

%% Different X and Y colors
clear all;
clf;
supererr([1 2 3], [1 3 2], [0.1 0.2; 0.3 0.6; 0.4 0.2], ...
    [0.2 0.4; 0.9 0.6; 0.4 0.8], ...
    'XColor', [.6 .6 .0], 'YColor', [.0 .6 .6]);
title('Different X and Y colors');
