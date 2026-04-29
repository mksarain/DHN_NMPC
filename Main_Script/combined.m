clear; clc; close all;

%% -------------------- Load saved simulation data --------------------
load('DHN_sim_results.mat', 'logsout');

%% -------------------- Create output folder for figures --------------------
figFolder = 'figures_exported';
if ~exist(figFolder, 'dir')
    mkdir(figFolder);
end

%% -------------------- Extract logged signals --------------------
P_ts   = logsout.getElement('P').Values;
q12_ts = logsout.getElement('q12').Values;
q23_ts = logsout.getElement('q23').Values;
q34_ts = logsout.getElement('q34').Values;
q41_ts = logsout.getElement('q41').Values;
T1_ts  = logsout.getElement('T1').Values;
T2_ts  = logsout.getElement('T2').Values;
T3_ts  = logsout.getElement('T3').Values;
T4_ts  = logsout.getElement('T4').Values;
Pd_ts  = logsout.getElement('Pd').Values;
st_ts  = logsout.getElement('nlp_status').Values;

%% -------------------- Global plot formatting --------------------
fontName = 'Times New Roman';
axisFS   = 14;
labelFS  = 16;
titleFS  = 17;
legendFS = 12;
lineW    = 2.0;

figW = 9;    % inches
figH = 10;   % taller figure for 3 stacked plots

% Journal-style color palette
c1 = [0.000 0.447 0.741];   % blue
c2 = [0.850 0.325 0.098];   % orange
c3 = [0.929 0.694 0.125];   % yellow
c4 = [0.494 0.184 0.556];   % purple
c5 = [0.466 0.674 0.188];   % green
c6 = [0.301 0.745 0.933];   % light blue
c7 = [0.635 0.078 0.184];   % dark red
c8 = [0.250 0.250 0.250];   % dark gray

% Distinct line styles
ls1 = '-';
ls2 = '--';
ls3 = ':';
ls4 = '-.';

%% -------------------- Convert signals --------------------
t_st = st_ts.Time(:)/3600;
y_st = squeeze(st_ts.Data); y_st = y_st(:);

t_T1 = T1_ts.Time(:)/3600; y_T1 = squeeze(T1_ts.Data); y_T1 = y_T1(:);
t_T2 = T2_ts.Time(:)/3600; y_T2 = squeeze(T2_ts.Data); y_T2 = y_T2(:);
t_T3 = T3_ts.Time(:)/3600; y_T3 = squeeze(T3_ts.Data); y_T3 = y_T3(:);
t_T4 = T4_ts.Time(:)/3600; y_T4 = squeeze(T4_ts.Data); y_T4 = y_T4(:);

t_q12 = q12_ts.Time(:)/3600; y_q12 = squeeze(q12_ts.Data); y_q12 = y_q12(:);
t_q23 = q23_ts.Time(:)/3600; y_q23 = squeeze(q23_ts.Data); y_q23 = y_q23(:);
t_q34 = q34_ts.Time(:)/3600; y_q34 = squeeze(q34_ts.Data); y_q34 = y_q34(:);
t_q41 = q41_ts.Time(:)/3600; y_q41 = squeeze(q41_ts.Data); y_q41 = y_q41(:);

t_P  = P_ts.Time(:)/3600;
y_P  = squeeze(P_ts.Data);  y_P  = y_P(:)/1000;

t_Pd = Pd_ts.Time(:)/3600;
y_Pd = squeeze(Pd_ts.Data); y_Pd = y_Pd(:)/1000;

%% =========================================================
%% 1) Optimization status (optional separate figure)
%% =========================================================
figure(1);
set(gcf, 'Color', 'w', 'Units', 'inches', 'Position', [1 1 9 4.5]);

plot(t_st, y_st, 'LineStyle', ls1, 'Color', c8, 'LineWidth', lineW);
grid on; box on;

ax = gca;
ax.FontName = fontName;
ax.FontSize = axisFS;
ax.LineWidth = 1.0;

xlabel('Time (h)', 'Interpreter', 'latex', 'FontSize', labelFS);
ylabel('Status code (-)', 'Interpreter', 'latex', 'FontSize', labelFS);
title('NMPC Optimization Status', 'Interpreter', 'latex', 'FontSize', titleFS);
legend({'$\mathrm{nlp\_status}$'}, ...
    'Interpreter', 'latex', 'FontSize', legendFS, 'Location', 'best');

xlim([0 24]);

exportgraphics(gcf, fullfile(figFolder,'optimization_status.png'), 'Resolution', 300);

%% =========================================================
%% 2) Combined comparison figure
%%    Top    : P and Pd
%%    Middle : all temperatures
%%    Bottom : all flow rates
%% =========================================================
figure(2);
set(gcf, 'Color', 'w', 'Units', 'inches', 'Position', [1 1 figW figH]);

tl = tiledlayout(3,1, 'TileSpacing', 'compact', 'Padding', 'compact');

%% -------------------- Top tile: Power --------------------
ax1 = nexttile;
plot(t_P,  y_P,  'LineStyle', ls1, 'Color', c1, 'LineWidth', lineW); hold on;
plot(t_Pd, y_Pd, 'LineStyle', ls2, 'Color', c7, 'LineWidth', lineW);
grid on; box on;

ax1.FontName = fontName;
ax1.FontSize = axisFS;
ax1.LineWidth = 1.0;

ylabel('Thermal power (kW)', 'Interpreter', 'latex', 'FontSize', labelFS);
title('DHN Closed-Loop Response', 'Interpreter', 'latex', 'FontSize', titleFS);
legend({'$P$', '$P_d$'}, ...
    'Interpreter', 'latex', 'FontSize', legendFS, ...
    'Location', 'best', 'Orientation', 'horizontal');

xlim([0 24]);

%% -------------------- Middle tile: Temperatures --------------------
ax2 = nexttile;
plot(t_T1, y_T1, 'LineStyle', ls1, 'Color', c1, 'LineWidth', lineW); hold on;
plot(t_T2, y_T2, 'LineStyle', ls2, 'Color', c2, 'LineWidth', lineW);
plot(t_T3, y_T3, 'LineStyle', ls3, 'Color', c5, 'LineWidth', lineW);
plot(t_T4, y_T4, 'LineStyle', ls4, 'Color', c4, 'LineWidth', lineW);
grid on; box on;

ax2.FontName = fontName;
ax2.FontSize = axisFS;
ax2.LineWidth = 1.0;

ylabel('Temperature ($^\circ$C)', 'Interpreter', 'latex', 'FontSize', labelFS);
legend({'$T_1$', '$T_2$', '$T_3$', '$T_4$'}, ...
    'Interpreter', 'latex', 'FontSize', legendFS, ...
    'Location', 'best', 'Orientation', 'horizontal');

xlim([0 24]);

%% -------------------- Bottom tile: Flow rates --------------------
ax3 = nexttile;
plot(t_q12, y_q12, 'LineStyle', ls1, 'Color', c7, 'LineWidth', lineW); hold on;
plot(t_q23, y_q23, 'LineStyle', ls2, 'Color', c6, 'LineWidth', lineW);
plot(t_q34, y_q34, 'LineStyle', ls3, 'Color', c3, 'LineWidth', lineW);
plot(t_q41, y_q41, 'LineStyle', ls4, 'Color', c8, 'LineWidth', lineW);
grid on; box on;

ax3.FontName = fontName;
ax3.FontSize = axisFS;
ax3.LineWidth = 1.0;

xlabel('Time (h)', 'Interpreter', 'latex', 'FontSize', labelFS);
ylabel('Mass flow rate (kg/s)', 'Interpreter', 'latex', 'FontSize', labelFS);
legend({'$q_{12}$', '$q_{23}$', '$q_{34}$', '$q_{41}$'}, ...
    'Interpreter', 'latex', 'FontSize', legendFS, ...
    'Location', 'best', 'Orientation', 'horizontal');

xlim([0 24]);

% Optional: link x-axes for consistent zoom/pan
linkaxes([ax1, ax2, ax3], 'x');

exportgraphics(gcf, fullfile(figFolder,'combined_dhn_response.png'), 'Resolution', 300);

disp('All figures exported to folder: figures_exported');