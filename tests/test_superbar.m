% Run and advance each cell and output is generated appropriately
clear all;
figure;

%% Plotting single bar
clear all;
clf;
superbar(5);
xlim([0 2]);
ylim([-1 6]);
title('Plotting a single bar');

%% Plot 2 bars, one at a time
clear all;
clf;
superbar(5);
hold on;
superbar(2, 5.5);
xlim([0 3]);
ylim([-1 6]);
title('Plot 2 bars, one at a time');

%% Plot 2 bars, at once (row vector)
clear all;
clf;
superbar([5, 5.5]);
xlim([0 3]);
ylim([-1 6]);
title('Plot 2 bars, at once (row vector)');

%% Plot 2 bars, at once (column vector)
clear all;
clf;
superbar([5, 5.5]');
xlim([0 3]);
ylim([-1 6]);
title('Plot 2 bars, at once (column vector)');

%% Manually specify bar width
clear all;
clf;
superbar([11 14 13]', 'BarWidth', 0.5);
xlim([0 4]);
title('Manually specify bar width');

%% Manual X location
clear all;
clf;
X = [ 4  8 14];
Y = [11 14 13];
superbar(X, Y');
title('Manual X location');

%% Groups of bars
clear all;
clf;
Y = [11 14 13;
     12 15 14];
superbar(Y);
title('Two groups of three bars');

%% Groups of bars, manually placed
clear all;
clf;
X = [1.0 1.25 1.6;
     2.0 2.4 2.7];
Y = [11 14 13;
     12 15 14];
superbar(X, Y);
title('Groups of bars, manually placed');

%% Empty X input
clear all;
clf;
Y = [11 14 13;
     12 15 14];
superbar([], Y');
title('Two groups of three bars, empty X input');

%% Plot into non-current axes
clear all;
clf;
Y = [11 14 13;
     12 15 14];
ax1 = subplot(1, 2, 1);
ax2 = subplot(1, 2, 2);
superbar(ax1, Y);
text(.5, .5, 'current axes', 'HorizontalAlignment', 'center');
title(ax1, 'Plot should be here');
title(ax2, 'Not here');

%% Offset base value
clear all;
clf;
Y = [11 14 13;
     12 15 14];
superbar(Y', 'BaseValue', 2);
ylim([0, 18]);
title('Offset base value');

%% Hollow bars
clear all;
clf;
Y = [11 14 13;
     15 12 16];
superbar(Y, 'BarFaceColor', 'none');
title('Hollow bars');

%% Coloured bars, by group
clear all;
clf;
Y = [11 14 13;
     15 12 16];
C = [.8 .2 .2;
     .2 .2 .8];
superbar(Y, 'BarFaceColor', C);
title('Coloured bars, by group');

%% Hollow coloured bars, by group
clear all;
clf;
Y = [11 14 13;
     15 12 16];
C = [.8 .2 .2;
     .2 .2 .8];
superbar(Y, 'BarFaceColor', 'none', 'BarEdgeColor', C);
title('Hollow coloured bars, by group');

%% Implicitly hollow coloured bars, by group
clear all;
clf;
Y = [11 14 13;
     15 12 16];
C = [.8 .2 .2;
     .2 .2 .8];
superbar(Y, 'BarEdgeColor', C);
title('Implicitly hollow coloured bars, by group');

%% Coloured bars, by index in group
clear all;
clf;
Y = [11 14 13;
     15 12 16];
C = [.8 .2 .2;
     .2 .8 .2;
     .2 .2 .8];
superbar(Y, 'BarFaceColor', permute(C, [3 1 2]));
title('Coloured bars, by index in group');

%% Coloured bars, individually
clear all;
clf;
Y = [11 14 13;
     15 12 16];
C = nan(2, 3, 3);
C(1, 1, :) = [.9 .3 .3];
C(1, 2, :) = [.3 .9 .3];
C(1, 3, :) = [.3 .3 .9];
C(2, 1, :) = [.6 .1 .1];
C(2, 2, :) = [.1 .6 .1];
C(2, 3, :) = [.1 .1 .6];
superbar(Y, 'BarFaceColor', C);
title('Coloured bars, individually');

%% Coloured bars, a single char
clear all;
clf;
Y = [11 14 13;
     15 12 16];
C = 'c';
superbar(Y, 'BarFaceColor', C);
title('Coloured bars, a single char');

%% Coloured bars, a char array
clear all;
clf;
Y = [11 14 13;
     15 12 16];
C = 'krb';
superbar(Y, 'BarFaceColor', C);
title('Coloured bars, char array');

%% Coloured bars, a char array by group
clear all;
clf;
Y = [11 14 13;
     15 12 16];
C = 'gm';
superbar(Y, 'BarFaceColor', C');
title('Coloured bars, char array by group');

%% Coloured bars, with a cell array of strings
clear all;
clf;
Y = [11 14 13;
     15 12 16];
C = {'k', 'r'};
superbar(Y, 'BarFaceColor', C);
title('Coloured bars, with a cell array of strings');

%% Coloured bars, with a cell array of strings and RGB vectors
clear all;
clf;
Y = [11 14 13;
     15 12 16];
C = {'g'; [1 0 0]};
superbar(Y, 'BarFaceColor', C);
title('Coloured bars, with a cell array of strings and RGB vectors');

%% Coloured bars, with a cell array of strings, RGBs, and 'none'
clear all;
clf;
Y = [11 14 13;
     15 12 16];
C = {'c', [1 0 0], [.2 .7 .4]};
CE = {'b', [.5 0.1 0], 'none'};
h = superbar(Y, 'BarFaceColor', C, 'BarEdgeColor', CE);
set(h, 'LineWidth', 3);
title('Coloured bars, with a cell array of strings, RGBs, and ''none''');

%% Bars with errorbars
clear all;
clf;
Y = [11 14 13;
     15 12 16];
E = [ 3  4  2;
      5  2  3];
C = [.8 .2 .2;
     .2 .2 .8];
superbar(Y, 'E', E, 'BarFaceColor', C);
title('Bars with errorbars');

%% Bars with errorbars, coloured with a char array
clear all;
clf;
Y = [11 14 13;
     15 12 16];
E = [ 3  4  2;
      5  2  3];
C = 'ckr';
superbar(Y, 'E', E, 'BarFaceColor', C);
title('Bars with errorbars, coloured with a char array');

%% Coloured bars, with a cell array of strings, RGBs, and 'none'
clear all;
clf;
Y = [11 14 13;
     15 12 16];
E = [ 3  4  2;
      5  2  3];
C = {'c', [1 0 0], [.2 .7 .4]};
CE = {'b', [.5 0.1 0], 'none'};
h = superbar(Y, 'E', E, 'BarFaceColor', C, 'BarEdgeColor', CE);
set(h, 'LineWidth', 3);
title('Coloured bars, with a cell array of strings, RGBs, and ''none''');

%% Bars with T errorbar style
clear all;
clf;
Y = [11 14 13;
     15 12 16];
E = [ 3  4  2;
      5  2  3];
C = [.8 .2 .2;
     .2 .2 .8];
superbar(Y, 'E', E, 'BarFaceColor', C, 'ErrorbarStyle', 'T');
title('Bars with T errorbar style');

%% Bar with p-values
clear all;
clf;
Y = [11 14 13;
     15 12 16];
P = [0.10 0.04   0.00001;
     0.05 0.0009 0.0017 ];
superbar(Y, 'P', P);
title('Bar with p-values');

%% Bar with errorbar and p-values
clear all;
clf;
Y = [11 14 13;
     15 12 16];
E = [ 3  4  2;
      5  2  3];
P = [0.10 0.04   0.00001;
     0.05 0.0009 0.0017 ];
superbar(Y, 'E', E, 'P', P);
title('Bar with errorbar and p-values');

%% Bar with massive errorbar and p-values
clear all;
clf;
X = [];
Y = [11 14 13;
     15 12 16];
E = [ 3  4  2;
      5  2  3];
P = [0.10 0.04   0.00001;
     0.05 0.0009 0.0017 ];
superbar(Y, 'E', 100 * E, 'P', P);
title('Bar with massive errorbar and p-values');

%% Horizontal bars with errorbar and p-values
clear all;
clf;
Y = [11 14 13;
     15 12 16];
E = [ 3  4  2;
      5  2  3];
P = [0.10 0.04   0.00001;
     0.05 0.0009 0.0017 ];
superbar(Y, 'E', E, 'P', P, 'Orientation', 'h');
title('Horizontal bars with errorbar and p-values');

%% Horizontal bars with errorbar and p-values (rotated)
clear all;
clf;
Y = [11 14 13;
     15 12 16];
E = [ 3  4  2;
      5  2  3];
P = [0.10 0.04   0.00001;
     0.05 0.0009 0.0017 ];
superbar(Y, 'E', E, 'P', P, 'Orientation', 'h', 'PStarFixedOrientation', false);
title('Horizontal bars with errorbar and p-values (rotated)');

%% Many pair-wise comparison p values
clear all;
clf;

Y = [11 14 13;
     15 12 16];
E = [ 3  4  2;
      5  2  3];
P = [NaN     1     1    1    1    NaN ;
     0.05    NaN   1    1    1    1   ;
     0.04    0.02  NaN  1    1    1   ;
     0.0009  0.005 1    NaN  1    1   ;
     0.00001 0.17  0.17 0.17 NaN  1   ;
     NaN     0.17  0.17 0.17 NaN  NaN ];
for i=1:sqrt(numel(P))
    P(i, i) = NaN;
end
% Make symmetric
P = min(P, P');

superbar(Y, 'E', E, 'P', P);
title('Many pair-wise comparison p values');

%% Many pair-wise comparison p values, horizontally
clear all;
clf;

Y = [11 14 13;
     15 12 16];
E = [ 3  4  2;
      5  2  3];
P = [NaN     1     1    1    1    NaN ;
     0.05    NaN   1    1    1    1   ;
     0.04    0.02  NaN  1    1    1   ;
     0.0009  0.005 1    NaN  1    1   ;
     0.00001 0.17  0.17 0.17 NaN  1   ;
     NaN     0.17  0.17 0.17 NaN  NaN ];
for i=1:sqrt(numel(P))
    P(i, i) = NaN;
end
% Make symmetric
P = min(P, P');

superbar(Y, 'E', E, 'P', P, 'orientation', 'h');
title('Many pair-wise comparison p values, horizontally');

%% Many pair-wise comparison p values, into non-current axes
clear all;
clf;

Y = [11 14 13;
     15 12 16];
E = [ 3  4  2;
      5  2  3];
P = [NaN     1     1    1    1    NaN ;
     0.05    NaN   1    1    1    1   ;
     0.04    0.02  NaN  1    1    1   ;
     0.0009  0.005 1    NaN  1    1   ;
     0.00001 0.17  0.17 0.17 NaN  1   ;
     NaN     0.17  0.17 0.17 NaN  NaN ];
for i=1:sqrt(numel(P))
    P(i, i) = NaN;
end
% Make symmetric
P = min(P, P');

ax1 = subplot(1, 2, 1);
ax2 = subplot(1, 2, 2);
superbar(ax1, Y, 'E', E, 'P', P);
text(.5, .5, 'current axes', 'HorizontalAlignment', 'center');
title(ax1, 'Plot should be here');
title(ax2, 'Not here');

%% Many pair-wise comparison p values, no source spacing
clear all;
clf;

Y = [11 14 13;
     15 12 16];
E = [ 3  4  2;
      5  2  3];
P = [NaN     1     1    1    1    NaN ;
     0.05    NaN   1    1    1    1   ;
     0.04    0.02  NaN  1    1    1   ;
     0.0009  0.005 1    NaN  1    1   ;
     0.00001 0.17  0.17 0.17 NaN  1   ;
     NaN     0.17  0.17 0.17 NaN  NaN ];
for i=1:sqrt(numel(P))
    P(i, i) = NaN;
end
% Make symmetric
P = min(P, P');

superbar(Y, 'E', E, 'P', P, 'PLineSourceRelativeSpacing', 0);
title('Many pair-wise comparison p values, no source spacing');

%% Many pair-wise comparison p values, no breadth spacing
clear all;
clf;

Y = [11 14 13;
     15 12 16];
E = [ 3  4  2;
      5  2  3];
P = [NaN     1     1    1    1    NaN ;
     0.05    NaN   1    1    1    1   ;
     0.04    0.02  NaN  1    1    1   ;
     0.0009  0.005 1    NaN  1    1   ;
     0.00001 0.17  0.17 0.17 NaN  1   ;
     NaN     0.17  0.17 0.17 NaN  NaN ];
for i=1:sqrt(numel(P))
    P(i, i) = NaN;
end
% Make symmetric
P = min(P, P');

superbar(Y, 'E', E, 'P', P, 'PLineSourceRelativeBreadth', 0);
title('Many pair-wise comparison p values, no breadth spacing');
