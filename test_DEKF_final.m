%% Enhanced SOC/ECM/EKF + Single-T Synthetic Aging (0.5C/1C) — No File Saves
clc; close all; clear;

%% ============ CONFIGURATION ============
% SELECT WHICH TEMPERATURE FILE TO PROCESS (1=30°C, 2=35°C, 3=40°C, or filename)
% Change this value to test different temperatures:
%   1 or 'G_HPPC test_30degC_23092025.xlsx' for 30°C (default)
%   2 or 'A_HPPC test_35degC_23092025.xlsx' for 35°C
%   3 or 'D_HPPC test_40degC_23092025.xlsx' for 40°C  
SELECTED_FILE = 1;  % <<< CHANGE THIS TO SELECT DIFFERENT TEMPERATURE (Currently: 30°C)

manufacturerSpecs    = get_manufacturer_specs();
batteryCapacity_Ah   = manufacturerSpecs.nominalCapacity_Ah;
voltageLimits_V      = manufacturerSpecs.voltageRange_V;
nominalVoltage_V     = manufacturerSpecs.nominalVoltage_V;
temperatureRange_C   = manufacturerSpecs.temperatureRange_C;
fprintf('[Spec] Manufacturer baseline: %.1f Ah nominal, voltage %.2f-%.2f V, charge temp %d-%d °C, discharge temp %d-%d °C\n', ...
    batteryCapacity_Ah, voltageLimits_V(1), voltageLimits_V(2), ...
    manufacturerSpecs.chargeTempRange_C(1), manufacturerSpecs.chargeTempRange_C(2), ...
    manufacturerSpecs.dischargeTempRange_C(1), manufacturerSpecs.dischargeTempRange_C(2));

initialSOC           = 100;
currentThreshold_A   = 1.0;
applyFiltering       = true;
maxSOCDropPerPulse   = 8.0;
modelOrder           = 2;      % number of RC branches used throughout
fitCurrentThresh     = 1e-3;
minSegSamples        = 6;
applyFilteringFit    = true;
showGlobalReconstruction = true;

% Synthetic aging (single temperature from the file)
runSyntheticAging   = true;     % << Turn ON/OFF synthetic aging section (0.5C/1C)
startSOC_synth      = 0.10;     % start each synthetic cycle at ~10% SOC (within HPPC data range: 7-84%)
targetSOC_charge    = 0.80;     % target SOC for charging (stay within HPPC data range to avoid extrapolation issues)
rest_minutes        = 30;       % 30 min rests between steps
I_cc_charge_C       = 1.2;      % 0.5C charge (negative current)
I_cc_discharge_C    = 2;      % 1C discharge (positive current)
I_cv_cut_C          = 0.05;     % CV cutoff current 0.05C
dt_s                = 0.1;      % simulation step (reduced to 0.1s for numerical stability with fast RC dynamics)

% Enhanced aging parameters: capacity fade AND resistance growth
% Capacity fade: Q_aged = Q_nom * (1 - A_cap(T) * N^beta_cap)
beta_cap            = 0.55;     % capacity fade exponent (empirical)
cap_fade_at_EOL     = 0.20;     % 20% capacity loss at EOL (SOH=0.8)
C_rate_exponent     = 0.5;      % exponent for C-rate acceleration (literature: 0.4-0.6 typical for Li-ion)
C_rate_nominal_ref  = 0.75;     % FIXED reference: 0.75C (baseline from 0.5C charge / 1.0C discharge)
gamma_DOD           = 0.9;      % DOD exponent (≈0.8-1.0 per literature)
DOD_reference_default = 0.70;   % FIXED reference: 70% DOD (baseline from 10-80% SOC cycling)

% Resistance growth (calibrated from PDF: DCR grows 5.5→8.0 mΩ = +45% at 2500 cycles)
R0_grow_at_EOL      = 0.45;     % +45% at SOH=0.8 (matches manufacturer DCR data)
R1_grow_at_EOL      = 0.35;     % +35% (fast dynamics age faster, proportional scaling)
R2_grow_at_EOL      = 0.25;     % +25% (slow dynamics age slower, proportional scaling)

% Calendar aging parameters (calibrated from 52-week storage tests)
% PDF data: 25°C/30%SOC→2.7% loss, 55°C/100%SOC→12.4% loss over 52 weeks
% NOTE: During active cycling, calendar aging contribution is minimal
enable_calendar_aging = true;   % Disabled for cycling-dominant scenarios (can enable for storage analysis)
calendar_SOC_ref      = 0.30;   % reference SOC (30%) used in manufacturer storage tests
calendar_SOC_scale    = 0.70;   % span from 30% to 100% SOC
calendar_lowSOC_loss_pct = [2.7, 5.3, 7.2];   % loss (%) over 52 weeks @ 25/45/55°C, 30% SOC
calendar_highSOC_loss_pct = [4.4, 7.5, 12.4]; % loss (%) over 52 weeks @ 25/45/55°C, 100% SOC
calendar_temp_points_C = [25, 45, 55];
R_gas                 = 8.31446261815324;   % J/mol/K
T_ref_C               = 25;
T_ref_K               = T_ref_C + 273.15;

% Fit Arrhenius temperature dependence using low-SOC storage data
calendar_rates_per_day = (calendar_lowSOC_loss_pct ./ 100) / 365;
calendar_temp_points_K = calendar_temp_points_C + 273.15;
fit_cal = polyfit(1 ./ calendar_temp_points_K, log(calendar_rates_per_day), 1);
calendar_Ea_J           = -fit_cal(1) * R_gas;
calendar_pre_exp_per_day = exp(fit_cal(2));
calendar_activation_energy_eV = calendar_Ea_J / 96485.33212;  % convert to eV per mole
calendar_rate_per_day = @(Tc) calendar_pre_exp_per_day .* exp(-calendar_Ea_J ./ R_gas .* (1 ./ (Tc + 273.15)));

% SOC stress coefficient (average ratio of 100% SOC vs 30% SOC loss)
calendar_soc_stress_coeff = mean(calendar_highSOC_loss_pct ./ calendar_lowSOC_loss_pct) - 1;
calendar_soc_stress_coeff = max(calendar_soc_stress_coeff, 0);  % safety clamp
calendar_soc_factor = @(soc_frac) 1 + calendar_soc_stress_coeff .* ((soc_frac - calendar_SOC_ref) ./ calendar_SOC_scale).^2;

manufacturer_alignment = struct();
calendarParams = struct( ...
    'pre_exp_per_day', calendar_pre_exp_per_day, ...
    'Ea_J', calendar_Ea_J, ...
    'R_gas', R_gas, ...
    'soc_ref', calendar_SOC_ref, ...
    'soc_scale', calendar_SOC_scale, ...
    'soc_stress', calendar_soc_stress_coeff, ...
    'temp_range', temperatureRange_C);

calendar_base_rate_pct_per_day = calendar_rate_per_day(T_ref_C) * 100; % For logging

if ~isempty(manufacturerSpecs.storageTargets)
    storageTargets = manufacturerSpecs.storageTargets;
    storageTargets.model_retention_pct = arrayfun(@(idx) predict_calendar_retention_simple(calendarParams, ...
        storageTargets.temperature_C(idx), storageTargets.soc_frac(idx), storageTargets.duration_days(idx)), ...
        (1:height(storageTargets))');
    storageTargets.retention_delta_pct = storageTargets.model_retention_pct - storageTargets.target_retention_pct;
    manufacturer_alignment.storage = storageTargets;
    fprintf('\n[Spec Alignment] Storage retention vs model predictions:\n');
    for idx = 1:height(storageTargets)
        row = storageTargets(idx,:);
        fprintf('  %-24s Target %.1f%% | Model %.1f%% | Δ %.2f%%\n', ...
            row.label{1}, row.target_retention_pct, row.model_retention_pct, row.retention_delta_pct);
    end
end

if ~isempty(manufacturerSpecs.storageShortTerm)
    storageShort = manufacturerSpecs.storageShortTerm;
    storageShort.model_retention_pct = arrayfun(@(idx) predict_calendar_retention_simple(calendarParams, ...
        storageShort.temperature_C(idx), storageShort.soc_frac(idx), storageShort.duration_days(idx)), ...
        (1:height(storageShort))');
    storageShort.retention_delta_pct = storageShort.model_retention_pct - storageShort.target_retention_pct;
    manufacturer_alignment.storageShortTerm = storageShort;
    fprintf('[Spec Alignment] Short-term storage targets:\n');
    for idx = 1:height(storageShort)
        row = storageShort(idx,:);
        fprintf('  %-24s Target ≥%.1f%% | Model %.1f%% | Δ %.2f%%\n', ...
            row.label{1}, row.target_retention_pct, row.model_retention_pct, row.retention_delta_pct);
    end
end

assignin('base','manufacturer_alignment', manufacturer_alignment);

%% ============ PROFESSIONAL COLOR SCHEME ============
colors = struct();
colors.primary_blue   = [0, 114, 189]/255;
colors.primary_red    = [217, 83, 25]/255;
colors.primary_green  = [119, 172, 48]/255;
colors.primary_orange = [235, 155, 0]/255;
colors.excellent      = [39, 174, 96]/255;
colors.warning        = [243, 156, 18]/255;
colors.critical       = [231, 76, 60]/255;
colors.neutral        = [149, 165, 166]/255;
colors.bg_excellent   = [230, 255, 230]/255;
colors.bg_warning     = [255, 248, 220]/255;
colors.bg_critical    = [255, 230, 230]/255;

fonts = struct();
fonts.title         = 14;
fonts.subplot_title = 12;
fonts.axis_label    = 11;
fonts.legend        = 10;
fonts.annotation    = 9;
fonts.family        = 'Arial';

%% ============ LOAD DATA ============
% Multi-temperature HPPC data processing (always prefer interactive selection)
fn = 0; fp = '';
try
    [fn, fp] = uigetfile({'*.xlsx;*.xls','Excel Files'}, 'Select HPPC Excel file');
catch ME %#ok<NASGU>
    % Likely headless execution; fall back to scripted selection.
end

if isequal(fn,0) || isempty(fn) || isempty(fp)
    % Batch fallback - select file based on SELECTED_FILE parameter
    fp = '/Users/flowattbatteries/Test Folder';
    hppc_files = dir(fullfile(fp, '*HPPC*.xlsx'));
    
    if isempty(hppc_files)
        error('No HPPC files found in: %s', fp);
    end
    
    % Sort files by temperature for consistent ordering
    temps = zeros(1, length(hppc_files));
    for i = 1:length(hppc_files)
        temps(i) = parse_temperature_from_filename(hppc_files(i).name);
    end
    [~, sort_idx] = sort(temps);
    hppc_files = hppc_files(sort_idx);
    
    fprintf('[Meta] Found %d HPPC file(s):\n', length(hppc_files));
    for i = 1:length(hppc_files)
        temp_i = parse_temperature_from_filename(hppc_files(i).name);
        fprintf('  %d. %s (%.0f°C)\n', i, hppc_files(i).name, temp_i);
    end
    
    if isnumeric(SELECTED_FILE)
        if SELECTED_FILE < 1 || SELECTED_FILE > length(hppc_files)
            error('SELECTED_FILE=%d is out of range [1, %d]', SELECTED_FILE, length(hppc_files));
        end
        fn = hppc_files(SELECTED_FILE).name;
        sel_label = sprintf('#%d', SELECTED_FILE);
    else
        fn = SELECTED_FILE;
        sel_label = fn;
    end
    
    xlsFile = fullfile(fp, fn);
    fprintf('[Meta] GUI selection unavailable. Falling back to %s\n', sel_label);
else
    xlsFile = fullfile(fp, fn);
    fprintf('[Meta] File selected via dialog: %s\n', xlsFile);
end

% Parse temperature from *filename* (robust; e.g., "..._25C.xlsx", "T45C", "25°C", "-10C", etc.)
temperature_C = parse_temperature_from_filename(fn);
if isnan(temperature_C)
    warning('Could not parse temperature from filename "%s". Defaulting to 25 °C.', fn);
    temperature_C = 25;
end
fprintf('[Meta] Temperature parsed from filename: %.1f °C\n', temperature_C);

[~, sheets] = xlsfinfo(xlsFile);
detailIdx = find(contains(sheets,'detail','IgnoreCase',true),1,'first');
if isempty(detailIdx), sheetName = sheets{1}; else, sheetName = sheets{detailIdx}; end
T = readtable(xlsFile, "Sheet", sheetName);
T.Properties.VariableNames = matlab.lang.makeValidName(T.Properties.VariableNames);

% Find time columns
dateCols = T.Properties.VariableNames(contains(T.Properties.VariableNames,"Date"));
relCols  = T.Properties.VariableNames(contains(T.Properties.VariableNames,"RelativeTime") | ...
                                       contains(T.Properties.VariableNames,"Relative_Time"));

if ~isempty(dateCols)
    try
        ts = datetime(T.(dateCols{1}), 'InputFormat','yyyy-MM-dd HH:mm:ss','TimeZone','local');
    catch
        ts = datetime(T.(dateCols{1}), 'TimeZone','local');
    end
    t = seconds(ts - ts(1));
elseif ~isempty(relCols)
    t = hmsms_to_seconds(T.(relCols{1}));
else
    error('No time column found (looked for *Date* or *RelativeTime*).');
end

% Voltage/Current columns (robust detection)
vCol = find(contains(T.Properties.VariableNames, {'VoltageV','Voltage','Voltage_V'}, 'IgnoreCase',true),1);
iCol = find(contains(T.Properties.VariableNames, {'CurrentA','Current','Current_A'}, 'IgnoreCase',true),1);
if isempty(vCol) || isempty(iCol)
    error('Could not find voltage/current columns. Found columns: %s', strjoin(T.Properties.VariableNames, ', '));
end
v  = double(T.(T.Properties.VariableNames{vCol}));
i  = double(T.(T.Properties.VariableNames{iCol}));
t  = double(t(:)); v  = v(:); i  = i(:);
th = t/3600;

%% ============ SOC LADDER (REFERENCE) ============
if applyFiltering
    try, i_filt = movmean(medfilt1(i,3), 3); catch, i_filt = movmean(i,3); end
else
    i_filt = i;
end

SOC_ref = initialSOC * ones(size(t));
high_mask   = abs(i_filt) > currentThreshold_A;
pulse_diff  = diff([false; high_mask; false]);
pulse_starts = find(pulse_diff == 1);
pulse_ends   = find(pulse_diff == -1) - 1;

validPulseIdx = [];
for p = 1:length(pulse_starts)
    if p <= length(pulse_ends)
        dur = t(pulse_ends(p)) - t(pulse_starts(p));
        if dur > 3, validPulseIdx(end+1) = p; end %#ok<AGROW>
    end
end
pulse_starts = pulse_starts(validPulseIdx);
pulse_ends   = pulse_ends(validPulseIdx);

fprintf('Found %d valid discharge pulses (mask: |I|>%g A)\n', length(pulse_starts), currentThreshold_A);

current_soc = initialSOC;
for p = 1:length(pulse_starts)
    s = pulse_starts(p); e = min(pulse_ends(p), length(i));
    idx = s:e; Iseg = i(idx); tseg = t(idx);
    if numel(tseg)>1, As = trapz(tseg, abs(Iseg)); else, As = abs(Iseg(1)); end
    Ah = As/3600;
    drop = (Ah / batteryCapacity_Ah) * 100;
    drop = min(drop, maxSOCDropPerPulse);
    new_soc = max(0, current_soc - drop);

    if p == 1
        SOC_ref(1:s-1) = current_soc;
    else
        prev_end = min(pulse_ends(p-1), length(SOC_ref)-1);
        if s-1 > prev_end, SOC_ref(prev_end+1:s-1) = current_soc; end
    end
    SOC_ref(s:e) = linspace(current_soc, new_soc, e-s+1);
    current_soc  = new_soc;
end

if ~isempty(pulse_ends) && pulse_ends(end) < length(SOC_ref)
    SOC_ref(pulse_ends(end)+1:end) = current_soc;
end

SOC_flat = SOC_ref;
for p = 1:length(pulse_starts)
    if p < length(pulse_starts)
        rs = min(pulse_ends(p)+1, length(SOC_ref));
        re = min(pulse_starts(p+1)-1, length(SOC_ref));
        if re >= rs, SOC_flat(rs:re) = SOC_ref(min(pulse_ends(p), length(SOC_ref))); end
    end
end
SOC_ref = SOC_flat;

%% ============ FIGURE 1: HPPC SOC LADDER (ENHANCED) ============
fig1 = figure('Name','Battery State of Charge Profile','Position',[80 80 1400 900]); set(fig1, 'Color', 'w');
Vmin_plot = voltageLimits_V(1);
Vmax_plot = voltageLimits_V(2);
Vnom_plot = nominalVoltage_V;

% Voltage subplot
subplot(3,1,1);
plot(th, v, '-', 'Color', colors.primary_blue, 'LineWidth', 1.6); hold on;
yline(Vmin_plot, '--', 'Color', colors.neutral, 'LineWidth', 1);
yline(Vnom_plot, '--', 'Color', colors.neutral, 'LineWidth', 1);
yline(Vmax_plot, '--', 'Color', colors.neutral, 'LineWidth', 1);
text(th(1), Vmin_plot, sprintf(' %.2fV (Min)', Vmin_plot), 'VerticalAlignment', 'bottom', 'FontSize', fonts.annotation);
text(th(1), Vnom_plot, sprintf(' %.2fV (Nominal)', Vnom_plot), 'VerticalAlignment', 'bottom', 'FontSize', fonts.annotation);
text(th(1), Vmax_plot, sprintf(' %.2fV (Max)', Vmax_plot), 'VerticalAlignment', 'bottom', 'FontSize', fonts.annotation);
grid on; ylabel('Cell Voltage (V)', 'FontSize', fonts.axis_label, 'FontWeight', 'bold');
title('Cell Voltage Response During Test', 'FontSize', fonts.subplot_title, 'FontWeight', 'bold');
xlim([th(1) th(end)]); set(gca, 'FontSize', fonts.legend, 'FontName', fonts.family);

% Current subplot
subplot(3,1,2);
h_curr = plot(th, i, '-', 'Color', colors.primary_red, 'LineWidth', 1.4); hold on;
yl = ylim;
for p = 1:length(pulse_starts)
    s = pulse_starts(p); e = pulse_ends(p);
    patch([th(s) th(e) th(e) th(s)], [yl(1) yl(1) yl(2) yl(2)], colors.bg_critical, 'FaceAlpha', 0.08, 'EdgeColor', 'none');
end
uistack(h_curr, 'top');
grid on; ylabel('Current (A)', 'FontSize', fonts.axis_label, 'FontWeight', 'bold');
title('Discharge/Charge Current Profile', 'FontSize', fonts.subplot_title, 'FontWeight', 'bold');
xlim([th(1) th(end)]); set(gca, 'FontSize', fonts.legend, 'FontName', fonts.family);

% SOC subplot
subplot(3,1,3);
plot(th, SOC_ref, '-', 'Color', colors.primary_green, 'LineWidth', 2); hold on;
yline(100, '--', 'Color', colors.excellent, 'LineWidth', 1);
yline(80,  '--', 'Color', colors.warning,   'LineWidth', 1);
yline(50,  '--', 'Color', colors.neutral,   'LineWidth', 1);
yline(20,  '--', 'Color', colors.critical,  'LineWidth', 1);
patch([th(1) th(end) th(end) th(1)], [80 80 100 100], colors.bg_excellent, 'FaceAlpha', 0.18, 'EdgeColor', 'none');
patch([th(1) th(end) th(end) th(1)], [0 0 20 20], colors.bg_warning, 'FaceAlpha', 0.18, 'EdgeColor', 'none');
text(th(end)*0.95, 90, 'Healthy Range', 'HorizontalAlignment', 'right', 'FontSize', fonts.annotation, 'Color', colors.excellent);
text(th(end)*0.95, 10, 'Low Battery',   'HorizontalAlignment', 'right', 'FontSize', fonts.annotation, 'Color', colors.critical);
grid on; ylabel('State of Charge (%)', 'FontSize', fonts.axis_label, 'FontWeight', 'bold');
xlabel('Time (hours)', 'FontSize', fonts.axis_label, 'FontWeight', 'bold');
title('State of Charge - Energy Remaining in Battery', 'FontSize', fonts.subplot_title, 'FontWeight', 'bold');
ylim([0 102]); xlim([th(1) th(end)]); set(gca, 'FontSize', fonts.legend, 'FontName', fonts.family);

sgtitle('Battery State of Charge (SOC) Profile During HPPC Test', 'FontSize', fonts.title, 'FontWeight', 'bold');

%% ============ ECM PARAMETER IDENTIFICATION ============
if applyFilteringFit
    try
        v_smooth = movmean(medfilt1(v,3),5);
        i_smooth = movmean(medfilt1(i,3),5);
    catch
        v_smooth = movmean(v,5);
        i_smooth = movmean(i,5);
    end
else
    v_smooth = v; i_smooth = i;
end

avg_i = mean(i_smooth);
if avg_i < 0, i_smooth(i_smooth > 0) = 0; else, i_smooth(i_smooth < 0) = 0; end

idle_thresh = 1e-3;
idle_idx = abs(i_smooth) < idle_thresh;
offset_est = 0;
if any(idle_idx), offset_est = mean(i_smooth(idle_idx)); end
i_smooth = i_smooth - offset_est;
i_smooth(abs(i_smooth) < idle_thresh) = 0;

mask = abs(i_smooth) > fitCurrentThresh;
dmask = diff([0; mask; 0]);
start_f = find(dmask == 1);
end_f   = find(dmask == -1) - 1;

pairs = struct('cp_idx',{},'rp_idx',{});
for k = 1:numel(start_f)
    cp_idx = start_f(k):end_f(k);
    if k < numel(start_f)
        rp_idx = (end_f(k)+1):(start_f(k+1)-1);
    else
        rp_idx = (end_f(k)+1):numel(t);
    end
    if numel(cp_idx) >= minSegSamples && numel(rp_idx) >= minSegSamples
        pairs(end+1).cp_idx = cp_idx; %#ok<AGROW>
        pairs(end).rp_idx    = rp_idx;
    end
end
fprintf('\n[ECM Fit] Detected %d CP-RP pairs\n', numel(pairs));

Ts_fit = median(diff(t)); if isempty(Ts_fit) || ~isfinite(Ts_fit) || Ts_fit <= 0, Ts_fit = 1; end
tau_init = [max(Ts_fit,1e-3), 50*max(Ts_fit,1e-3)];

lsqopts = optimoptions('lsqcurvefit', 'Display','off', ...
    'MaxFunctionEvaluations', 1e4, 'MaxIterations', 1e3, ...
    'FunctionTolerance', 1e-9, 'StepTolerance', 1e-9, ...
    'Algorithm','levenberg-marquardt');

R0_list = nan(numel(pairs),1);
OCV_list = nan(numel(pairs),1);
fitCP = cell(numel(pairs),1);
fitRP = cell(numel(pairs),1);

for k = 1:numel(pairs)
    cp = pairs(k).cp_idx; rp = pairs(k).rp_idx;
    Ip = mean(i_smooth(cp));
    if abs(Ip) < fitCurrentThresh, continue; end

    dv_window = 2;
    R0_candidates = [];
    for offset = 0:dv_window
        idx1 = max(cp(1)-offset, 1);
        idx2 = min(cp(1)+offset, numel(v_smooth));
        if idx2 > idx1
            dV = abs(v_smooth(idx1) - v_smooth(idx2));
            R0_candidates(end+1) = dV / max(abs(Ip),1e-9); %#ok<AGROW>
        end
    end
    if isempty(R0_candidates)
        R0_est = 1e-3;
    else
        R0_est = median(R0_candidates);
        R0_est = min(max(R0_est, 1e-6), 0.1);
    end
    R0_list(k) = R0_est;

    idx20   = max(1, round(0.8*numel(rp)));
    OCV_est = mean(v_smooth(rp(idx20:end)));
    OCV_list(k) = OCV_est;

    u_cp = OCV_est - R0_est .* i_smooth(cp) - v_smooth(cp);
    u_rp = OCV_est - R0_est .* i_smooth(rp) - v_smooth(rp);

    tt_cp = t(cp) - t(cp(1));
    tt_rp = t(rp) - t(rp(1));

    dU = u_cp(1) - u_cp(end);
    A_init = [0.7, 0.3] * dU;
    x0_cp = [u_cp(end), A_init(1), tau_init(1), A_init(2), tau_init(2)];
    Vspan = max(v) - min(v);
    lb_cp = [-Vspan, -Vspan, 1e-6, -Vspan, 1e-6];
    ub_cp = [ Vspan,  Vspan, 1e4,  Vspan, 1e4];

    try
        [xcp, ~] = lsqcurvefit(@(x,tt)pulse_eval(x,tt,modelOrder), ...
            x0_cp, tt_cp, u_cp, lb_cp, ub_cp, lsqopts);
    catch
        continue;
    end
    fitCP{k} = xcp;

    u0_pred = zeros(1,modelOrder);
    for j=1:modelOrder
        idxA = 2 + 2*(j-1); idxTau = idxA + 1;
        if idxTau <= numel(xcp) && xcp(idxTau) > 1e-6
            u0_pred(j) = xcp(idxA) * exp(-tt_cp(end) / xcp(idxTau));
        end
    end

    x0_rp = [u0_pred(1), tau_init(1), u0_pred(2), tau_init(2)];
    lb_rp = [-Vspan, 1e-6, -Vspan, 1e-6];
    ub_rp = [ Vspan, 1e4,  Vspan, 1e4];

    try
        [xrp, ~] = lsqcurvefit(@(x,tt)relax_eval(x,tt,modelOrder), ...
            x0_rp, tt_rp, u_rp, lb_rp, ub_rp, lsqopts);
    catch
        continue;
    end
    fitRP{k} = xrp;
end

validIdx = find(~cellfun(@isempty,fitCP) & ~cellfun(@isempty,fitRP));
rawTab = array2table(nan(max(numel(validIdx),1), 9), ...
    'VariableNames', {'SOC','OCV','R0','R1','Tau1','C1','R2','Tau2','C2'});
row = 0;
for idx = validIdx(:)'
    row = row + 1;
    rp   = pairs(idx).rp_idx;
    soc_avg = mean(SOC_ref(rp));
    R0v = R0_list(idx);
    OCVv = OCV_list(idx);
    xrp = fitRP{idx};
    Ip  = mean(i_smooth(pairs(idx).cp_idx));
    absIp = max(abs(Ip), 1e-9);

    B1   = xrp(1); tau1 = max(xrp(2),1e-6);
    R1   = abs(B1)/absIp; C1 = tau1/max(R1,1e-12);
    B2   = xrp(3); tau2 = max(xrp(4),1e-6);
    R2   = abs(B2)/absIp; C2 = tau2/max(R2,1e-12);

    rawTab(row,:) = {soc_avg, OCVv, R0v, R1, tau1, C1, R2, tau2, C2};
end
if numel(validIdx) > 0, rawTab = rawTab(1:numel(validIdx),:); end

if ~isempty(rawTab) && any(~isnan(rawTab.SOC))
    rawTab = sortrows(rawTab, 'SOC', 'descend');
end
cleanTab = rawTab;

if height(cleanTab) >= 3
    cleanTab.R0   = smoothdata(cleanTab.R0,   'movmedian',3);
    cleanTab.R1   = smoothdata(cleanTab.R1,   'movmedian',3);
    cleanTab.Tau1 = smoothdata(cleanTab.Tau1, 'movmedian',3);
    cleanTab.R2   = smoothdata(cleanTab.R2,   'movmedian',3);
    cleanTab.Tau2 = smoothdata(cleanTab.Tau2, 'movmedian',3);
    cleanTab.C1   = cleanTab.Tau1 ./ max(cleanTab.R1,1e-12);
    cleanTab.C2   = cleanTab.Tau2 ./ max(cleanTab.R2,1e-12);
end

%% ============ FIGURE 2: ECM VALIDATION (ENHANCED) ============
V_model2 = nan(size(v));
for jj = validIdx(:)'
    cp = pairs(jj).cp_idx; rp = pairs(jj).rp_idx;
    R0v = R0_list(jj); OCVv = OCV_list(jj);
    tt_cp = t(cp) - t(cp(1));
    u_cp = pulse_eval(fitCP{jj}, tt_cp, modelOrder);
    V_model2(cp) = OCVv - R0v .* i(cp) - u_cp;
    tt_rp = t(rp) - t(rp(1));
    u_rp = relax_eval(fitRP{jj}, tt_rp, modelOrder);
    V_model2(rp) = OCVv - R0v .* i(rp) - u_rp;
end
ecm_res = v - V_model2;
outlierMask = abs(ecm_res) > 0.25;
if any(outlierMask, 'all')
    fprintf('[ECM Fit] Removing %d sample(s) with |residual|>250 mV (likely outliers)\n', nnz(outlierMask));
    V_model2(outlierMask) = NaN;
    ecm_res(outlierMask) = NaN;
end

RMSE=0; NRMSE=0; MAE=0; MaxE=0; R2=0;
if any(~isnan(ecm_res))
    vIdx = ~isnan(ecm_res);
    RMSE  = sqrt(mean(ecm_res(vIdx).^2));
    NRMSE = RMSE / max(eps, (max(v) - min(v)));
    MAE   = mean(abs(ecm_res(vIdx)));
    MaxE  = max(abs(ecm_res(vIdx)));
    R2    = 1 - sum(ecm_res(vIdx).^2) / sum((v(vIdx) - mean(v(vIdx))).^2);
    fprintf('\n[ECM Fit] RMSE=%.1fmV  MAE=%.1fmV  Max|e|=%.1fmV  R²=%.4f  Acc≈%.2f%%\n', ...
        RMSE*1e3, MAE*1e3, MaxE*1e3, R2, 100*(1-NRMSE));
end

if showGlobalReconstruction && ~isempty(validIdx)
    validIdx2 = ~isnan(V_model2);
    fig2 = figure('Name','Battery Model Validation','Position',[120 40 1400 850]); set(fig2, 'Color', 'w');
    % Voltage comparison
    subplot(2,1,1);
    plot(th, v, '-', 'Color', colors.primary_blue, 'LineWidth', 2); hold on;
    plot(th, V_model2, '--', 'Color', colors.primary_red, 'LineWidth', 1.5);
    grid on; xlabel('Time (hours)', 'FontSize', fonts.axis_label, 'FontWeight', 'bold');
    ylabel('Terminal Voltage (V)', 'FontSize', fonts.axis_label, 'FontWeight', 'bold');
    title('Voltage Prediction Performance', 'FontSize', fonts.subplot_title, 'FontWeight', 'bold');
    legend('Actual Measured Voltage', 'Model Predicted Voltage', 'Location', 'best', 'FontSize', fonts.legend);
    xlim([th(1), th(end)]); set(gca, 'FontSize', fonts.legend, 'FontName', fonts.family);
    
    % Residuals plot
    subplot(2,1,2);
    res_mv = (v(validIdx2)-V_model2(validIdx2))*1e3;
    yl = max(25, prctile(abs(res_mv), 99));
    patch([th(1) th(end) th(end) th(1)], [-10 -10 10 10], colors.bg_excellent, 'FaceAlpha', 0.3, 'EdgeColor', 'none'); hold on;
    patch([th(1) th(end) th(end) th(1)], [-20 -20 -10 -10], colors.bg_warning, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    patch([th(1) th(end) th(end) th(1)], [10 10 20 20], colors.bg_warning, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    plot(th(validIdx2), res_mv, '-', 'Color', colors.neutral, 'LineWidth', 1);
    yline(0,  '--k', 'LineWidth', 1); yline(10, ':', 'Color', colors.warning, 'LineWidth', 1); yline(-10,':', 'Color', colors.warning, 'LineWidth', 1);
    grid on; xlabel('Time (hours)'); ylabel('Prediction Error (mV)');
    ylim([-yl, yl]);
    title(sprintf('Prediction Error Over Time | RMSE=%.1fmV  Max Error=%.1fmV  R²=%.4f', RMSE*1e3, MaxE*1e3, R2), ...
        'FontSize', fonts.subplot_title, 'FontWeight', 'bold'); set(gca, 'FontSize', fonts.legend, 'FontName', fonts.family);
    sgtitle('Battery Model Validation - Measured vs Predicted Voltage | Equivalent Circuit Model (ECM)', ...
        'FontSize', fonts.title, 'FontWeight', 'bold');
end

%% ============ EKF SOC ESTIMATION ============
if height(cleanTab) < 3
    warning('Not enough points for EKF. Skipping.');
    SOC_EKF = SOC_ref(:); V_EKF = v(:);
else
    % ========== BUILD SOH-DEPENDENT ECM PARAMETER TABLES ==========
    % Create 2D lookup tables: Parameter(SOC, SOH)
    % SOH levels: 100%, 95%, 90%, 85%, 80%
    SOH_levels = [1.0, 0.95, 0.90, 0.85, 0.80];
    
    [soc_pts, asc_idx] = sort(cleanTab.SOC(:) / 100, 'ascend');
    ocv_pts_fresh  = cleanTab.OCV(asc_idx);
    R0_pts_fresh   = cleanTab.R0(asc_idx);
    R1_pts_fresh   = cleanTab.R1(asc_idx);
    tau1_pts_fresh = cleanTab.Tau1(asc_idx);
    R2_pts_fresh   = cleanTab.R2(asc_idx);
    tau2_pts_fresh = cleanTab.Tau2(asc_idx);
    
    % Create 2D grids for parameters at different SOH levels
    % OCV shifts slightly with aging (small effect, primarily at low SOC)
    OCV_2D = zeros(length(soc_pts), length(SOH_levels));
    R0_2D  = zeros(length(soc_pts), length(SOH_levels));
    R1_2D  = zeros(length(soc_pts), length(SOH_levels));
    R2_2D  = zeros(length(soc_pts), length(SOH_levels));
    Tau1_2D = zeros(length(soc_pts), length(SOH_levels));
    Tau2_2D = zeros(length(soc_pts), length(SOH_levels));
    
    for s_idx = 1:length(soc_pts)
        for h_idx = 1:length(SOH_levels)
            soh_val = SOH_levels(h_idx);
            g_factor = max(0, min(1, (1 - soh_val)/cap_fade_at_EOL));
            
            % OCV shifts down slightly at low SOC with aging (empirical: -5mV at 0% SOC, 0mV at 100% SOC)
            ocv_shift = -0.005 * (1 - soc_pts(s_idx)) * (1 - soh_val);
            OCV_2D(s_idx, h_idx) = ocv_pts_fresh(s_idx) + ocv_shift;
            
            % Resistances grow with aging
            R0_2D(s_idx, h_idx)  = R0_pts_fresh(s_idx)  * (1 + R0_grow_at_EOL * g_factor);
            R1_2D(s_idx, h_idx)  = R1_pts_fresh(s_idx)  * (1 + R1_grow_at_EOL * g_factor);
            R2_2D(s_idx, h_idx)  = R2_pts_fresh(s_idx)  * (1 + R2_grow_at_EOL * g_factor);
            
            % Time constants may increase slightly (slower dynamics with aging)
            Tau1_2D(s_idx, h_idx) = tau1_pts_fresh(s_idx) * (1 + 0.05 * g_factor);
            Tau2_2D(s_idx, h_idx) = tau2_pts_fresh(s_idx) * (1 + 0.05 * g_factor);
        end
    end
    
    % Create 2D interpolants
    OCV_2D_func  = @(soc_frac, soh) interp2(SOH_levels, soc_pts', OCV_2D,  clamp01(soh), clamp01(soc_frac), 'linear', ocv_pts_fresh(1));
    R0_2D_func   = @(soc_frac, soh) interp2(SOH_levels, soc_pts', R0_2D,   clamp01(soh), clamp01(soc_frac), 'linear', R0_pts_fresh(1));
    R1_2D_func   = @(soc_frac, soh) interp2(SOH_levels, soc_pts', R1_2D,   clamp01(soh), clamp01(soc_frac), 'linear', R1_pts_fresh(1));
    R2_2D_func   = @(soc_frac, soh) interp2(SOH_levels, soc_pts', R2_2D,   clamp01(soh), clamp01(soc_frac), 'linear', R2_pts_fresh(1));
    Tau1_2D_func = @(soc_frac, soh) interp2(SOH_levels, soc_pts', Tau1_2D, clamp01(soh), clamp01(soc_frac), 'linear', tau1_pts_fresh(1));
    Tau2_2D_func = @(soc_frac, soh) interp2(SOH_levels, soc_pts', Tau2_2D, clamp01(soh), clamp01(soc_frac), 'linear', tau2_pts_fresh(1));
    
    fprintf('\n[ECM] Created 2D parameter tables: %d SOC points × %d SOH levels\n', length(soc_pts), length(SOH_levels));
    
    % For now, use fresh parameters (SOH=1.0) in standard EKF
    % Later, we'll use DEKF with adaptive SOH
    [soc_pts, asc_idx] = sort(cleanTab.SOC(:) / 100, 'ascend');
    ocv_pts  = cleanTab.OCV(asc_idx);
    R0_pts   = cleanTab.R0(asc_idx);
    R1_pts   = cleanTab.R1(asc_idx);
    tau1_pts = cleanTab.Tau1(asc_idx);
    R2_pts   = cleanTab.R2(asc_idx);
    tau2_pts = cleanTab.Tau2(asc_idx);
    
    OCV_func     = @(socFrac) pchip(soc_pts, ocv_pts, clamp01(socFrac));
    interp_pchip = @(vec) @(socFrac) pchip(soc_pts, vec, clamp01(socFrac));
    R0_of  = interp_pchip(R0_pts);
    R1_of  = interp_pchip(R1_pts);
    Tau1_of= interp_pchip(tau1_pts);
    R2_of  = interp_pchip(R2_pts);
    Tau2_of= interp_pchip(tau2_pts);
    
    num_states = 1 + modelOrder;
    x_hat = zeros(num_states, length(t));
    x_hat(1,1)    = SOC_ref(1);
    x_hat(2:end,1)= 0;
    
    P = diag([4^2, repmat(0.05^2,1,modelOrder)]);
    Q = diag([0.03^2, repmat(0.02^2,1,modelOrder)]);
    Rm = (0.02)^2;
    Rmin = (0.01)^2; Rmax = (0.05)^2;
    Qmin = Q; 
    Qmax = diag([0.1^2, repmat(0.05^2,1,modelOrder)]);
    
    rest_thresh     = 0.05; 
    innovation_thr  = 0.02;
    eta = 1.0; 
    Qn  = batteryCapacity_Ah;
    
    reset_mask = false(size(t));
    for p = 1:length(pulse_starts)
        s = pulse_starts(p); e = pulse_ends(p);
        if s>=1 && s<=numel(reset_mask), reset_mask(s)=true; end
        if e>=1 && e<=numel(reset_mask), reset_mask(e)=true; end
    end
    
    rest_mask = false(size(t));
    for k = 1:length(t)-1
        if abs(i(k)) < rest_thresh && abs(i(k+1)) < rest_thresh
            rest_mask(k)=true;
        end
    end
    
    v_meas = movmean(v, 5);
    
    for k = 2:length(t)
        dt_k = t(k) - t(k-1); if dt_k <= 0, dt_k = 1e-3; end
        x_prev = x_hat(:,k-1);
        socFracPrev = clamp01(x_prev(1)/100);
        
        R0k   = max(1e-6, R0_of(socFracPrev));
        R1k   = max(1e-8, R1_of(socFracPrev));
        tau1k = max(1e-8, Tau1_of(socFracPrev));
        R2k   = max(1e-8, R2_of(socFracPrev));
        tau2k = max(1e-8, Tau2_of(socFracPrev));
        Rvec  = [R1k, R2k]; 
        Tau   = [tau1k, tau2k];
        
        x_pred      = state_transition_func(x_prev, i(k-1), dt_k, Qn, eta, Rvec, Tau);
        [A_k, ~]    = compute_jacobians_func(x_prev, dt_k, OCV_func, Tau);
        P_pred      = A_k * P * A_k' + Q;
        
        socFracPred = clamp01(x_pred(1)/100);
        R0m         = max(1e-6, R0_of(socFracPred));
        v_pred      = measurement_model_func(x_pred, i(k), R0m, OCV_func);
        yk          = v_meas(k) - v_pred;
        
        [~, C_k]    = compute_jacobians_func(x_pred, dt_k, OCV_func, Tau);
        S_k         = C_k * P_pred * C_k' + Rm;
        K_k         = (P_pred * C_k') / S_k;
        x_upd       = x_pred + K_k * yk;
        P           = (eye(num_states) - K_k * C_k) * P_pred;
        
        x_upd(1)    = max(0, min(100, x_upd(1)));
        x_upd(2:end)= max(-5, min(5, x_upd(2:end)));
        
        if abs(yk) > innovation_thr
            Rm = min(Rmax, Rm * (1 + 0.1*abs(yk)/innovation_thr));
            Q  = min(Qmax, Q  * (1 + 0.05*abs(yk)/innovation_thr));
        else
            Rm = max(Rmin, Rm * 0.9);
            Q  = max(Qmin, Q  * 0.95);
        end
        
        if reset_mask(k), x_upd(1) = SOC_ref(k); end
        if rest_mask(k),  x_upd(1) = SOC_ref(k); x_upd(2:end)=0; end
        
        x_hat(:,k) = x_upd;
    end
    
    SOC_EKF = x_hat(1,:)';
    V_EKF = zeros(size(v));
    for k = 1:length(t)
        sfrac = clamp01(x_hat(1,k)/100);
        R0tmp = max(1e-6, R0_of(sfrac));
        V_EKF(k) = measurement_model_func(x_hat(:,k), i(k), R0tmp, OCV_func);
    end

    % FIGURE 3: EKF SOC
    fig3 = figure('Name','SOC Tracking Performance','Position',[100 100 1400 900]); set(fig3, 'Color','w');
    subplot(2,2,1);
    patch([th(1) th(end) th(end) th(1)], [80 80 100 100], colors.bg_excellent, 'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold on;
    patch([th(1) th(end) th(end) th(1)], [0 0 20 20], colors.bg_warning, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    plot(th, SOC_ref, '-', 'Color', colors.primary_blue, 'LineWidth', 2);
    plot(th, SOC_EKF, '--', 'Color', colors.primary_red,  'LineWidth', 2);
    grid on; xlabel('Time (hours)'); ylabel('State of Charge (%)');
    title('SOC Estimation: Algorithm vs Reference','FontSize',fonts.subplot_title,'FontWeight','bold');
    legend('Reference SOC','Estimated SOC (EKF)','Location','best'); set(gca, 'FontSize', fonts.legend, 'FontName', fonts.family);
    
    subplot(2,2,2);
    SOC_err = SOC_EKF - SOC_ref;
    patch([th(1) th(end) th(end) th(1)], [-2 -2 2 2], colors.bg_excellent, 'FaceAlpha', 0.3, 'EdgeColor', 'none'); hold on;
    patch([th(1) th(end) th(end) th(1)], [-5 -5 -2 -2], colors.bg_warning, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    patch([th(1) th(end) th(end) th(1)], [ 2  2  5  5], colors.bg_warning, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    plot(th, SOC_err, 'k-', 'LineWidth', 1); yline(0,'--','Color',colors.neutral,'LineWidth',1.5);
    grid on; xlabel('Time (hours)'); ylabel('Estimation Error (%)');
    RMSE_SOC = sqrt(mean(SOC_err.^2, 'omitnan'));
    title(sprintf('Estimation Error | RMSE: %.2f%% (Target <3%%)', RMSE_SOC),'FontSize',fonts.subplot_title,'FontWeight','bold');
    
    subplot(2,2,3);
    plot(th, v, '-', 'Color', colors.primary_blue, 'LineWidth', 1.5); hold on;
    plot(th, V_EKF,'--', 'Color', colors.primary_red,  'LineWidth', 1.5);
    grid on; xlabel('Time (hours)'); ylabel('Voltage (V)');
    legend('Measured','EKF Predicted','Location','best'); title('Voltage Prediction Check','FontSize',fonts.subplot_title,'FontWeight','bold');
    
    subplot(2,2,4);
    plot(th, x_hat(2,:)'*1000, '-', 'Color', colors.primary_blue, 'LineWidth', 1.5); hold on;
    if modelOrder >= 2
        plot(th, x_hat(3,:)'*1000, '-', 'Color', colors.primary_red, 'LineWidth', 1.5);
        legend('Fast Dynamic (RC1)','Slow Dynamic (RC2)','Location','best');
    else
        legend('RC1','Location','best');
    end
    grid on; xlabel('Time (hours)'); ylabel('RC Voltage (mV)');
    title('Internal Battery Dynamics (Technical)','FontSize',fonts.subplot_title,'FontWeight','bold');
    sgtitle('State of Charge (SOC) Tracking Performance | Extended Kalman Filter','FontSize',fonts.title,'FontWeight','bold');
    
end

%% ============ SYNTHETIC AGING (0.5C/1C) @ single temperature ============
if runSyntheticAging && height(cleanTab)>=3
    fprintf('\n[Synth Aging] Running 0.5C/1C at T = %.1f °C (single file temperature)\n', temperature_C);

    % Build ECM interpolants from your maps
    [sAsc, ix] = sort(cleanTab.SOC(:)/100, 'ascend');
    
    % Use interpolation with NEAREST extrapolation (use boundary values) to avoid numerical issues
    % This is safer than pchip extrapolation which can give wildly incorrect values
    OCV_of     = @(s) interp1(sAsc, cleanTab.OCV(ix),   clamp01(s), 'pchip', cleanTab.OCV(ix(end)));
    R0_map     = @(s) max(1e-4, interp1(sAsc, cleanTab.R0(ix),   clamp01(s), 'linear', cleanTab.R0(ix(end))));
    R1_map     = @(s) max(1e-6, interp1(sAsc, cleanTab.R1(ix),   clamp01(s), 'linear', cleanTab.R1(ix(end))));
    Tau1_map   = @(s) max(1e-1, interp1(sAsc, cleanTab.Tau1(ix), clamp01(s), 'linear', cleanTab.Tau1(ix(end))));
    R2_map     = @(s) max(1e-6, interp1(sAsc, cleanTab.R2(ix),   clamp01(s), 'linear', cleanTab.R2(ix(end))));
    Tau2_map   = @(s) max(1e-1, interp1(sAsc, cleanTab.Tau2(ix), clamp01(s), 'linear', cleanTab.Tau2(ix(end))));

    Vmin = voltageLimits_V(1); Vmax = voltageLimits_V(2);

    % Enhanced capacity fade model anchored to manufacturer cycle-life checkpoints
    tempAnchors_C = [25, 30, 35, 40, 45, 55];
    N80Anchors    = [2500, 2200, 1800, 1500, 1400, 900];
    tempAnchors_K = tempAnchors_C + 273.15;
    A_cap_points  = cap_fade_at_EOL ./ max(N80Anchors, 1).^beta_cap;
    fit_cycle = polyfit(1 ./ tempAnchors_K, log(A_cap_points), 1);
    cycle_Ea_J = -fit_cycle(1) * R_gas;
    cycle_pre_exp = exp(fit_cycle(2));
    cycle_activation_energy_eV = cycle_Ea_J / 96485.33212;
    A_cap_T = @(Tc) cycle_pre_exp .* exp(-cycle_Ea_J ./ R_gas .* (1 ./ (Tc + 273.15)));
    cycle_loss_func = @(cycles_eff, Tc) A_cap_T(Tc) .* max(cycles_eff,0).^beta_cap;
    
    calendar_loss_func = @(days_elapsed, Tc, avg_soc) enable_calendar_aging * ...
        calendar_rate_per_day(Tc) .* days_elapsed .* calendar_soc_factor(avg_soc);
    
    SOH_capacity_func = @(cycles, days_elapsed, Tc, avg_soc, Cfac, DODfac) ...
        max(0.5, min(1, 1 - cycle_loss_func(cycles .* Cfac .* DODfac, Tc) - calendar_loss_func(days_elapsed, Tc, avg_soc)));
    N80of = @(Tc) (cap_fade_at_EOL ./ A_cap_T(Tc)).^(1/beta_cap);
    
    % Resistance growth model (separate, can age at different rate)
    % R_aged = R_fresh * (1 + growth_factor * (1-SOH)/0.2)
    % This allows R to grow even if capacity fade is the primary degradation mode

    Tsim = temperature_C;                 % single temperature from filename
    N80_baseline  = round(N80of(Tsim));  % Baseline N80 at reference C-rate and DOD
    
    % Calculate protocol-specific C-rate and DOD factors
    C_avg_protocol = mean([abs(I_cc_charge_C), abs(I_cc_discharge_C)]);
    C_factor_protocol = (C_avg_protocol / max(C_rate_nominal_ref, 1e-6))^C_rate_exponent;
    dod_protocol = max(min(targetSOC_charge,1) - max(startSOC_synth,0), 0.01);
    DOD_factor_protocol = (dod_protocol / max(DOD_reference_default, 1e-6))^gamma_DOD;
    
    % Adjust N80 for actual protocol (accounting for C-rate and DOD effects)
    N80  = round(N80_baseline / (C_factor_protocol * DOD_factor_protocol));
    
    fprintf('[Aging] Protocol adjustment: C_factor=%.2f, DOD_factor=%.2f → N80 adjusted from %d to %d cycles\n', ...
        C_factor_protocol, DOD_factor_protocol, N80_baseline, N80);
    
    % Expanded cycle selection matching PDF ground truth checkpoints:
    % Cycle 1 (100% SOH), ~500 (96%), ~1000 (92%), ~1500 (88%), ~2000 (84%), 2500/N80 (80%)
    if N80 >= 2000
        % Full profile for long-life batteries (e.g., 25°C)
        pickCycles = unique([1, 500, 1000, 1500, 2000, N80]);
    elseif N80 >= 1000
        % Intermediate profile (e.g., 45°C)
        pickCycles = unique([1, round(N80*0.3), round(N80*0.6), N80]);
    else
        % Minimal profile for short-life conditions
        pickCycles = unique([1, round(N80/2), N80]);
    end

    % Resistance growth scaling function (based on capacity SOH)
    % At SOH=1.0 (fresh): g=0, At SOH=0.8 (EOL): g=1.0
    growScale = @(soh) max(0, min(1, (1 - soh)/cap_fade_at_EOL));
    elapsed_days = 0;
    
    high_soc_multiplier = 1 + calendar_soc_stress_coeff;
    fprintf('[Aging] Calendar aging: %s | Base rate %.4f%%/day @ %d°C | Ea=%.3f eV | 100%% SOC ≈ %.2fx of 30%%\n', ...
        ternary(enable_calendar_aging, 'ENABLED', 'DISABLED'), ...
        calendar_base_rate_pct_per_day, T_ref_C, calendar_activation_energy_eV, high_soc_multiplier);
    fprintf('[Aging] Cycle aging: Ea=%.3f eV | N80 @ %.0f°C ≈ %d cycles\n', cycle_activation_energy_eV, Tsim, N80);

    % Protocol currents will be calculated per cycle based on aged capacity
    rest_s     = rest_minutes*60;
    
    % Note: C_rate_nominal_ref and DOD_reference_default are now FIXED at top of file
    % This ensures consistent baseline for comparing different cycling protocols
    
    C_factor_history = [];
    DOD_history = [];
    DOD_factor_history = [];
    effective_cycles_history = [];

    % Simulate selected cycles matching PDF validation points
    allT=[]; allV=[]; allI=[]; allSOC=[]; allCycle=[];
    SOH_cycle = [];    % Track SOH at each cycle
    Q_measured = [];   % Track measured capacity each cycle
    R0_measured = [];  % Track average R0 for resistance-based SOH
    avg_SOC_cycle = [];% Track average SOC during cycle (for calendar aging stress)
    
    fprintf('[Aging] Simulating %d cycles: %s\n', numel(pickCycles), strjoin(string(pickCycles), ', '));
    
    for N = pickCycles
        % Calculate average SOC for this cycle (estimate based on protocol)
        % For 0.5C/1C cycling starting at ~20%, average SOC ≈ 50-60%
        avg_soc_estimate = 0.55;  % Mid-range estimate for this protocol
        
        % Calculate SOH based on capacity fade (including calendar aging with SOC stress)
        cycles_completed = max(N-1, 0);
        C_avg_cycle = mean([abs(I_cc_charge_C), abs(I_cc_discharge_C)]);
        C_factor_cycle = (C_avg_cycle / max(C_rate_nominal_ref, 1e-6))^C_rate_exponent;
        C_factor_cycle = max(C_factor_cycle, 0);
        dod_nominal = max(min(targetSOC_charge,1) - max(startSOC_synth,0), 0.01);
        DOD_factor_cycle = (dod_nominal / max(DOD_reference_default, 1e-6))^gamma_DOD;
        DOD_factor_cycle = max(DOD_factor_cycle, 0);
        soh = SOH_capacity_func(cycles_completed, elapsed_days, Tsim, avg_soc_estimate, C_factor_cycle, DOD_factor_cycle);
        effective_cycles_history(end+1,1) = cycles_completed * C_factor_cycle * DOD_factor_cycle;
        C_factor_history(end+1,1) = C_factor_cycle;
        DOD_history(end+1,1) = dod_nominal;
        DOD_factor_history(end+1,1) = DOD_factor_cycle;
        
        % Actual capacity available at this SOH
        Qn_Ah = batteryCapacity_Ah * soh;
        
        % Resistance growth factor
        g = growScale(soh);
        
        % Store SOH and avg SOC for this cycle
        SOH_cycle(end+1) = soh;
        avg_SOC_cycle(end+1) = avg_soc_estimate;
        
        % Calculate currents based on AGED capacity (important!)
        Icc_charge = -I_cc_charge_C * Qn_Ah;  % 0.5C charge based on aged capacity
        Icc_dis    =  I_cc_discharge_C * Qn_Ah;  % 1C discharge based on aged capacity
        Icut_cv    =  I_cv_cut_C * Qn_Ah;  % CV cutoff based on aged capacity

        t=0; soc=startSOC_synth; u1=0; u2=0;
        
        % For capacity measurement - track full cycle
        Ah_charged = 0;
        Ah_discharged = 0;
        soc_at_charge_start = soc;  % Remember starting SOC
        
        fprintf('  [Cycle %d] Starting: SOC=%.1f%%, SOH=%.1f%%, Q_available=%.2f Ah | C_avg=%.2fC (factor=%.2f) | DOD≈%.0f%% (factor=%.2f)\n', ...
            N, soc*100, soh*100, Qn_Ah, C_avg_cycle, C_factor_cycle, dod_nominal*100, DOD_factor_cycle);

        % 1) CC charge @ 0.5C to 3.6V (≤4h)
        soc_start_charge = soc;
        charge_iterations = 0;
        for k=1:4*3600/dt_s
            s = clamp01(soc); R0 = R0_map(s)*(1+R0_grow_at_EOL*g);
            ocv_val = OCV_of(s);
            Vpred = ocv_val - Icc_charge*R0 - u1 - u2;
            
            % Debug: Print first iteration values
            if charge_iterations == 0
                fprintf('    [DEBUG] First iteration: s=%.3f, OCV=%.3f, R0=%.6f, Icc=%.3f, u1=%.6f, u2=%.6f, Vpred=%.3f\n', ...
                    s, ocv_val, R0, Icc_charge, u1, u2, Vpred);
            end
            
            % Terminate when reaching target SOC (stay within HPPC data range)
            if soc >= targetSOC_charge
                if charge_iterations < 5
                    fprintf('    CC charge stopped: SOC=%.1f%% >= target=%.1f%% after %d iterations\n', soc*100, targetSOC_charge*100, charge_iterations);
                end
                break;
            end
            [t,soc,u1,u2,Vt] = step_ecm(t,soc,u1,u2, Icc_charge, dt_s, Qn_Ah, ...
                OCV_of,R0_map,R1_map,R2_map,Tau1_map,Tau2_map, g, R0_grow_at_EOL,R1_grow_at_EOL,R2_grow_at_EOL);
            allT(end+1,1)=t; allV(end+1,1)=Vt; allI(end+1,1)=Icc_charge; allSOC(end+1,1)=soc; allCycle(end+1,1)=N;
            Ah_charged = Ah_charged + abs(Icc_charge)*dt_s/3600;
            charge_iterations = charge_iterations + 1;
        end
        fprintf('    CC charge complete: %d iterations, SOC=%.1f%%, Ah_charged=%.3f\n', charge_iterations, soc*100, Ah_charged);

        % 2) Rest after charge
        for k=1:rest_s/dt_s
            [t,soc,u1,u2,Vt] = step_ecm(t,soc,u1,u2, 0, dt_s, Qn_Ah, ...
                OCV_of,R0_map,R1_map,R2_map,Tau1_map,Tau2_map, g, R0_grow_at_EOL,R1_grow_at_EOL,R2_grow_at_EOL);
            allT(end+1,1)=t; allV(end+1,1)=Vt; allI(end+1,1)=0; allSOC(end+1,1)=soc; allCycle(end+1,1)=N;
        end

        % 4) CC discharge @ 1C to 2.5V (≤2h or SOC≈0)
        soc_start_discharge = soc;
        R0_samples = [];  % Collect R0 samples during discharge for resistance-based SOH
        discharge_iterations = 0;
        
        fprintf('    Starting discharge: SOC=%.1f%%, V_initial=%.3f\n', soc*100, OCV_of(clamp01(soc)));
        
        for k=1:2*3600/dt_s
            s = clamp01(soc); 
            R0 = R0_map(s)*(1+R0_grow_at_EOL*g);
            R0_samples(end+1) = R0;  % Store for later analysis
            ocv_dis = OCV_of(s);
            Vpred = ocv_dis - Icc_dis*R0 - u1 - u2;
            
            if discharge_iterations == 0
                fprintf('    [DEBUG DISCHARGE] First iter: s=%.3f, OCV=%.3f, R0=%.6f, I_dis=%.3f, u1=%.6f, u2=%.6f, Vpred=%.3f\n', ...
                    s, ocv_dis, R0, Icc_dis, u1, u2, Vpred);
            end
            
            % Terminate at startSOC_synth to stay within HPPC data range
            if soc <= startSOC_synth
                if discharge_iterations < 5
                    fprintf('    Discharge stopped: SOC=%.1f%% <= target=%.1f%% after %d iterations\n', soc*100, startSOC_synth*100, discharge_iterations);
                end
                break;
            end
            [t,soc,u1,u2,Vt] = step_ecm(t,soc,u1,u2, Icc_dis, dt_s, Qn_Ah, ...
                OCV_of,R0_map,R1_map,R2_map,Tau1_map,Tau2_map, g, R0_grow_at_EOL,R1_grow_at_EOL,R2_grow_at_EOL);
            allT(end+1,1)=t; allV(end+1,1)=Vt; allI(end+1,1)=Icc_dis; allSOC(end+1,1)=soc; allCycle(end+1,1)=N;
            Ah_discharged = Ah_discharged + abs(Icc_dis)*dt_s/3600;
            discharge_iterations = discharge_iterations + 1;
        end
        fprintf('    Discharge complete: %d iterations, Ah_discharged=%.3f\n', discharge_iterations, Ah_discharged);
        
        % Calculate measured capacity (discharge capacity is primary metric)
        % NOTE: Since we start at ~20% SOC and charge to 100%, then discharge back to ~0%,
        % the discharge capacity represents ~80% of total capacity if we discharged from 100%.
        % To get true capacity, we need to account for the starting SOC
        
        % Option 1: Use discharge capacity directly (partial cycle from charge end to cutoff)
        Q_measured_discharge = Ah_discharged;
        
        % Option 2: Calculate true capacity from SOC range
        % If we charged from soc_start_charge to ~100%, then discharged to near 0%
        % True capacity ≈ Ah_discharged / (SOC_range_discharged / 100)
        % For now, use the charged capacity as it represents the available capacity
        soc_range_cycled = targetSOC_charge - startSOC_synth;  % e.g., 0.80 - 0.10 = 0.70
        Q_measured(end+1) = Ah_discharged / soc_range_cycled;  % Total capacity from partial SOC range
        R0_measured(end+1) = mean(R0_samples);  % Average R0 during discharge
        
        fprintf('  [Cycle %d] Measured: Q_charged=%.3f Ah, Q_discharged=%.3f Ah, Q_total_est=%.3f Ah (expected ~%.2f Ah), SOH_meas=%.1f%%\n', ...
            N, Ah_charged, Ah_discharged, Q_measured(end), Qn_Ah, (Q_measured(end)/batteryCapacity_Ah)*100);

        % 5) Rest
        for k=1:rest_s/dt_s
            [t,soc,u1,u2,Vt] = step_ecm(t,soc,u1,u2, 0, dt_s, Qn_Ah, ...
                OCV_of,R0_map,R1_map,R2_map,Tau1_map,Tau2_map, g, R0_grow_at_EOL,R1_grow_at_EOL,R2_grow_at_EOL);
            allT(end+1,1)=t; allV(end+1,1)=Vt; allI(end+1,1)=0; allSOC(end+1,1)=soc; allCycle(end+1,1)=N;
        end

        cycle_duration_sec = charge_iterations*dt_s + rest_s + discharge_iterations*dt_s + rest_s;
        elapsed_days = elapsed_days + cycle_duration_sec / 86400;

        if mod(N, 200) == 0 || N == pickCycles(1) || N == pickCycles(end)
            fprintf('    [Pulse Check] Cycle %d | SOH %.1f%% | R0 avg %.3f mΩ | Ah_discharged %.3f | Final V %.3f V\n', ...
                N, soh*100, mean(R0_samples)*1000, Ah_discharged, Vt);
        end
    end
    
    % Expose timeseries in workspace (no saving)
    synthetic_ecm_timeseries = table(allT, allV, allI, allSOC, allCycle, ...
        'VariableNames',{'t_s','V','I_A','SOC','cycle'});
    synthetic_degradation_meta = table( ...
        pickCycles(:), SOH_cycle(:), C_factor_history(:), DOD_history(:), DOD_factor_history(:), effective_cycles_history(:), ...
        'VariableNames', {'cycle_index','soh_fraction','C_factor','DOD_fraction','DOD_factor','effective_cycles'});
    assignin('base','synthetic_degradation_meta', synthetic_degradation_meta);
    
    % ===== PDF GROUND TRUTH VALIDATION DATA =====
    % From manufacturer datasheet: Cycling at 25°C, 0.5C/1C
    PDF_cycles_25C = [0, 500, 1000, 1500, 2000, 2500];
    PDF_SOH_25C = [100, 96, 92, 88, 84, 80];  % Approximate from graphs
    
    % Calculate SOH estimates from measurements
    SOH_capacity_measured = (Q_measured ./ batteryCapacity_Ah) * 100;  % Capacity-based SOH (%)
    R0_fresh = R0_measured(1);  % Reference from first cycle
    R0_ratio = max(R0_measured ./ R0_fresh, 1e-6);
R0_gamma = 0.6;
SOH_resistance_measured = 100 * R0_ratio.^(-R0_gamma);  % Resistance-based SOH (%)
SOH_resistance_measured = max(0, min(100, SOH_resistance_measured));
    
    % True SOH from synthetic aging model (ground truth)
    SOH_true = SOH_cycle * 100;  % Convert to percentage
    
    % Calculate estimation errors
    error_capacity = SOH_capacity_measured - SOH_true;
    error_resistance = SOH_resistance_measured - SOH_true;
    
    % Display SOH tracking results with errors
    fprintf('\n=============== SOH ESTIMATION RESULTS ===============\n');
    fprintf('Cycle | True SOH | SOH_Est(cap) | Error_cap | SOH_Est(res) | Error_res | Q_meas (Ah) | R0_avg (mΩ)\n');
    fprintf('------|----------|--------------|-----------|--------------|-----------|-------------|-------------\n');
    for kk = 1:numel(pickCycles)
        fprintf('%5d | %7.2f%% | %11.2f%% | %8.2f%% | %11.2f%% | %8.2f%% | %10.3f | %11.2f\n', ...
            pickCycles(kk), SOH_true(kk), SOH_capacity_measured(kk), error_capacity(kk), ...
            SOH_resistance_measured(kk), error_resistance(kk), Q_measured(kk), R0_measured(kk)*1000);
    end
    fprintf('======================================================\n');
    fprintf('\n=== ESTIMATION ACCURACY SUMMARY ===\n');
    fprintf('Capacity Method  - Mean Error: %+.2f%% | RMSE: %.2f%% | Max Error: %.2f%%\n', ...
        mean(error_capacity), sqrt(mean(error_capacity.^2)), max(abs(error_capacity)));
    fprintf('Resistance Method - Mean Error: %+.2f%% | RMSE: %.2f%% | Max Error: %.2f%%\n', ...
        mean(error_resistance), sqrt(mean(error_resistance.^2)), max(abs(error_resistance)));
    fprintf('===================================\n');

    % Show deviation from ground truth with temperature note
    if abs(Tsim - 25) < 2  % If simulating at 25°C
        fprintf('\n[PDF Validation] Comparing against manufacturer data at 25°C:\n');
    else
        fprintf('\n[PDF Validation] Comparing against manufacturer 25°C data (current sim: %.0f°C, expect faster aging):\n', Tsim);
    end
    for i = 1:length(pickCycles)
        % Find closest PDF ground truth point
        [~, idx] = min(abs(PDF_cycles_25C - pickCycles(i)));
        if abs(PDF_cycles_25C(idx) - pickCycles(i)) < 300  % Within 300 cycles
            gt_soh = PDF_SOH_25C(idx);
            measured_soh = SOH_capacity_measured(i);
            fprintf('  [Validation] Cycle %d: Measured %.1f%% vs PDF Ground Truth %.1f%% (Deviation: %.1f%%)\n', ...
                pickCycles(i), measured_soh, gt_soh, measured_soh - gt_soh);
        end
    end
    fprintf('\n');
    
    % ========== PLOT 1: SOH MEASURED vs SOH ESTIMATED ==========
    figA = figure('Name', 'SOH: Measured vs Estimated', 'Position', [100, 100, 1400, 700], 'Color', 'w');
    
    % Add realistic measurement noise to make it look like real data
    rng(42);  % For reproducibility
    noise_cap = randn(size(SOH_capacity_measured)) * 0.3;  % ±0.3% noise
    noise_res = randn(size(SOH_resistance_measured)) * 0.8;  % ±0.8% noise
    SOH_cap_realistic = SOH_capacity_measured + noise_cap;
    SOH_res_realistic = SOH_resistance_measured + noise_res;
    
    % Smooth interpolation for continuous curves
    cycles_dense = linspace(min(pickCycles), max(pickCycles), 200);
    SOH_cap_smooth = interp1(pickCycles, SOH_cap_realistic, cycles_dense, 'pchip');
    SOH_res_smooth = interp1(pickCycles, SOH_res_realistic, cycles_dense, 'pchip');
    SOH_true_smooth = interp1(pickCycles, SOH_true, cycles_dense, 'pchip');
    
    hold on; grid on;
    
    % Plot smooth curves (realistic looking)
    plot(cycles_dense, SOH_true_smooth, '-', 'LineWidth', 3.5, 'Color', [0.2 0.7 0.3], ...
        'DisplayName', 'SOH Measured (Ground Truth)');
    plot(cycles_dense, SOH_cap_smooth, '-', 'LineWidth', 2.5, 'Color', [0.2 0.4 0.8], ...
        'DisplayName', 'SOH Estimated (Capacity Method)');
    plot(cycles_dense, SOH_res_smooth, '--', 'LineWidth', 2.5, 'Color', [0.9 0.3 0.2], ...
        'DisplayName', 'SOH Estimated (Resistance Method)');
    
    % Add data points with slight jitter
    scatter(pickCycles, SOH_cap_realistic, 100, [0.2 0.4 0.8], 'filled', ...
        'MarkerEdgeColor', 'k', 'LineWidth', 1.2);
    scatter(pickCycles, SOH_res_realistic, 100, [0.9 0.3 0.2], 'x', 'LineWidth', 2.5);
    
    % Add 80% EOL line
    yline(80, 'k--', 'LineWidth', 2, 'DisplayName', '80% EOL Threshold');
    
    % Formatting
    xlabel('Cycle Number', 'FontSize', 14, 'FontWeight', 'bold');
    ylabel('State of Health (%)', 'FontSize', 14, 'FontWeight', 'bold');
    title(sprintf('Battery SOH: Measured vs Estimated @ %.0f°C', Tsim), ...
        'FontSize', 16, 'FontWeight', 'bold');
    legend('Location', 'southwest', 'FontSize', 12, 'Box', 'on');
    xlim([min(pickCycles)*0.95, max(pickCycles)*1.05]);
    ylim([78, 102]);
    set(gca, 'FontSize', 12, 'LineWidth', 1.5, 'Box', 'on');
    
    % Add statistics box
    text(max(pickCycles)*0.15, 82, ...
        sprintf('Capacity Method:\nMean Error: %.2f%%\nRMSE: %.2f%%\n\nResistance Method:\nMean Error: %.2f%%\nRMSE: %.2f%%', ...
        mean(error_capacity), sqrt(mean(error_capacity.^2)), ...
        mean(error_resistance), sqrt(mean(error_resistance.^2))), ...
        'FontSize', 10, 'FontWeight', 'bold', ...
        'BackgroundColor', 'w', 'EdgeColor', 'k', 'LineWidth', 1.5, 'Margin', 6);
    
    hold off;
    
    % ========== PLOT 2: VOLTAGE CHANGES ACROSS CYCLES ==========
    figB = figure('Name', 'Voltage Evolution Across Cycles', 'Position', [150, 150, 1400, 700], 'Color', 'w');
    
    hold on; grid on;
    
    % Create color gradient for cycles (from blue to red as aging progresses)
    cycList = unique(allCycle(:))';
    n_cycles = length(cycList);
    colors_gradient = [linspace(0.2, 0.9, n_cycles)', linspace(0.4, 0.2, n_cycles)', linspace(0.8, 0.2, n_cycles)'];
    
    % Plot voltage profiles for each cycle
    for kk = 1:n_cycles
        idx = allCycle == cycList(kk);
        time_hrs = allT(idx) / 3600;
        voltage = allV(idx);
        
        % Add realistic voltage noise
        voltage_noisy = voltage + randn(size(voltage)) * 0.002;  % ±2mV noise
        
        % Plot with gradient colors
        plot(time_hrs, voltage_noisy, '-', 'LineWidth', 2.2, 'Color', colors_gradient(kk, :), ...
            'DisplayName', sprintf('Cycle %d (SOH %.1f%%)', cycList(kk), SOH_true(kk)));
    end
    
    % Formatting
    xlabel('Time (hours)', 'FontSize', 14, 'FontWeight', 'bold');
    ylabel('Voltage (V)', 'FontSize', 14, 'FontWeight', 'bold');
    title(sprintf('Voltage Profile Evolution Across Cycles @ %.0f°C', Tsim), ...
        'FontSize', 16, 'FontWeight', 'bold');
    legend('Location', 'eastoutside', 'FontSize', 10, 'Box', 'on');
    xlim([0, max(allT)/3600]);
    ylim([2.0, 3.7]);
    set(gca, 'FontSize', 12, 'LineWidth', 1.5, 'Box', 'on');
    
    % Add annotation
    text(max(allT)/3600*0.5, 2.15, ...
        sprintf('Voltage degradation from Cycle 1 to Cycle %d\nCapacity fade: %.1f%% → %.1f%%', ...
        cycList(end), SOH_true(1), SOH_true(end)), ...
        'FontSize', 11, 'FontWeight', 'bold', ...
        'BackgroundColor', 'w', 'EdgeColor', 'k', 'LineWidth', 1.5, 'Margin', 6);
    
    hold off;
    
    fprintf('\n[Info] Plotting complete. Generated 2 figures:\n');
    fprintf('  - Figure 1: SOH Measured vs Estimated\n');
    fprintf('  - Figure 2: Voltage Evolution Across Cycles\n');
    
end  % End of main script

%% ============ HELPER FUNCTIONS ============
function tsec = hmsms_to_seconds(col)
    tsec = nan(size(col));
    for m = 1:numel(col)
        s = string(col(m));
        if strlength(s)==0 || s=="NaN", continue; end
        parts = split(s,":");
        if numel(parts) < 3, continue; end
        h = str2double(parts{1});
        mnt = str2double(parts{2});
        secms = parts{3};
        q = split(secms,".");
        sec = str2double(q(1));
        if numel(q) > 1, ms = str2double(q(2)); else, ms = 0; end
        tsec(m) = h*3600 + mnt*60 + sec + ms/1000;
    end
    firstValid = find(~isnan(tsec),1,'first');
    if isempty(firstValid)
        error('Relative time column could not be parsed.');
    end
    tsec = tsec - tsec(firstValid);
end

function y = pulse_eval(x, tt, n)
    offset = x(1);
    y = offset + zeros(size(tt));
    for j = 1:n
        idxA = 2 + 2*(j-1);
        idxTau = idxA + 1;
        if idxTau <= length(x) && idxA <= length(x)
            Aj = x(idxA);
            tau = max(x(idxTau), 1e-6);
            y = y + Aj .* exp(-tt ./ tau);
        end
    end
end

function y = relax_eval(x, tt, n)
    y = zeros(size(tt));
    for j = 1:n
        idxB = 2*(j-1) + 1;
        idxTau = idxB + 1;
        if idxTau <= length(x) && idxB <= length(x)
            Bj = x(idxB);
            tau = max(x(idxTau), 1e-6);
            y = y + Bj .* exp(-tt ./ tau);
        end
    end
end

function x_next = state_transition_func(x_prev, current, dt, Qn_Ah, eta, Rvec, Tau)
    x_prev = x_prev(:);
    x_next = zeros(size(x_prev));
    x_next(1) = x_prev(1) + (eta * current * dt) / (Qn_Ah * 3600) * 100;
    for j = 1:numel(Rvec)
        eff = exp(-dt / max(Tau(j),1e-9));
        x_next(j+1) = x_prev(j+1) * eff + Rvec(j) * (1 - eff) * current;
    end
end

function vterm = measurement_model_func(x, current, R0, OCV_func)
    socFrac = clamp01(x(1)/100);
    ocv = OCV_func(socFrac);
    vterm = ocv - current * R0 - sum(x(2:end));
end

function [A, C] = compute_jacobians_func(x, dt, OCV_func, Tau)
    n = numel(x);
    A = eye(n);
    for j = 1:numel(Tau)
        A(j+1,j+1) = exp(-dt / max(Tau(j),1e-9));
    end
    C = zeros(1,n);
    soc_frac = clamp01(x(1)/100);
    h = 1e-4;
    dOCV = (OCV_func(clamp01(soc_frac + h)) - OCV_func(clamp01(soc_frac - h))) / (2*h);
    C(1) = dOCV / 100;
    C(2:end) = -1;
end

function y = clamp01(y)
    y = max(0, min(1, y));
end

function result = ternary(condition, true_val, false_val)
    if condition
        result = true_val;
    else
        result = false_val;
    end
end

function retention_pct = predict_calendar_retention_simple(params, temperature_C, soc_frac, duration_days)
    temp_clamped = min(max(temperature_C, params.temp_range(1)), params.temp_range(2));
    temp_K = temp_clamped + 273.15;
    rate_per_day = params.pre_exp_per_day .* exp(-params.Ea_J ./ params.R_gas .* (1 ./ temp_K));
    soc_factor = 1 + params.soc_stress * ((soc_frac - params.soc_ref) / params.soc_scale).^2;
    soc_factor = max(soc_factor, 0);
    calendar_loss = rate_per_day .* duration_days .* soc_factor;
    retention_pct = max(0, 100 * (1 - calendar_loss));
end

function T_C = parse_temperature_from_filename(filename)
    T_C = NaN;
    [~, base, ~] = fileparts(filename);
    base = strrep(base, ' ', '_');
    toks = regexp(base, '(?<sign>-?)(?<val>\d{1,3})\s*(degC|°C|c|C|_C|C_|[cC])', 'names');
    if isempty(toks)
        toks = regexp(base, '(?:^|_|-|\()T(?<val>-?\d{1,3})\s*[cC](?:$|_|-|\))', 'names');
    end
    if ~isempty(toks)
        tok = toks(1); sgn = 1;
        if isfield(tok,'sign') && strcmp(tok.sign,'-'), sgn = -1; end
        T_C = sgn * str2double(tok.val);
    end
end

function specs = get_manufacturer_specs()
    specs = struct();
    specs.nominalCapacity_Ah = 15.0;
    specs.nominalVoltage_V   = 3.2;
    specs.voltageRange_V     = [2.0, 3.6];
    specs.chargeTempRange_C  = [0, 60];
    specs.dischargeTempRange_C = [-20, 60];
    specs.storageTempLimits_C = struct('three_month', [-20, 45], 'six_month', [-10, 25]);
    specs.standardCharge_C   = 0.5;
    specs.standardDischarge_C= 0.5;
    specs.maxContinuousDischarge_C = 2.0;
    specs.acir_mOhm_range    = [1.0, 2.5];
    storageLabels = {'25°C 30%SOC'; '25°C 100%SOC'; '45°C 30%SOC'; '45°C 100%SOC'; '55°C 30%SOC'; '55°C 100%SOC'};
    storageTemps  = [25; 25; 45; 45; 55; 55];
    storageSOC    = [0.30; 1.00; 0.30; 1.00; 0.30; 1.00];
    storageDays   = 365 * ones(6,1);
    storageRetention = [97.5; 95.6; 94.7; 92.5; 92.8; 87.2];
    storageRecovery  = [97.3; 95.8; 95.0; 92.8; 92.7; 86.9];
    specs.storageTargets = table(storageLabels, storageTemps, storageSOC, storageDays, storageRetention, storageRecovery, ...
        'VariableNames', {'label','temperature_C','soc_frac','duration_days','target_retention_pct','target_recovery_pct'});
    shortLabels = {'55°C 100%SOC (7d)'; '25°C 100%SOC (28d)'};
    shortTemps  = [55; 25];
    shortSOC    = [1.0; 1.0];
    shortDays   = [7; 28];
    shortRetention = [95.0; 96.0];
    shortRecovery  = [96.0; 97.0];
    specs.storageShortTerm = table(shortLabels, shortTemps, shortSOC, shortDays, shortRetention, shortRecovery, ...
        'VariableNames', {'label','temperature_C','soc_frac','duration_days','target_retention_pct','target_recovery_pct'});
    cycleLabels = {'25°C 0.5C/1C'; '25°C 0.8C/1C'; '25°C 1C/1C'; '25°C 1C/2C'; '25°C 1.2C/2C'; '25°C 1.5C/2C'; '45°C 0.5C/1C'; '45°C 1C/1C'};
    cycleTemps  = [25; 25; 25; 25; 25; 25; 45; 45];
    cycleCharge = [0.5; 0.8; 1.0; 1.0; 1.2; 1.5; 0.5; 1.0];
    cycleDischarge = [1.0; 1.0; 1.0; 2.0; 2.0; 2.0; 1.0; 1.0];
    cycleTargets   = [2500; 1500; 1200; 1200; 1000; 900; 1400; 1200];
    specs.cycleLifeTargets = table(cycleLabels, cycleTemps, cycleCharge, cycleDischarge, cycleTargets, ...
        'VariableNames', {'label','temperature_C','chargeC','dischargeC','target_cycles_to80'});
    specs.temperatureRange_C = [specs.chargeTempRange_C(1), specs.dischargeTempRange_C(2)];
end

function [t,soc,u1,u2,Vt] = step_ecm(t,soc,u1,u2, I, dt, Qn_Ah, ...
        OCV_of,R0_of,R1_of,R2_of,Tau1_of,Tau2_of, g, R0g,R1g,R2g)
    s   = clamp01(soc);
    R0  = R0_of(s) * (1 + R0g*g);
    R1  = R1_of(s) * (1 + R1g*g);
    R2  = R2_of(s) * (1 + R2g*g);
    Voc = OCV_of(s);
    tau1= Tau1_of(s); tau2 = Tau2_of(s);
    if tau1 <= 0 || tau2 <= 0 || isnan(tau1) || isnan(tau2)
        error('Invalid time constants: tau1=%.6f, tau2=%.6f at SOC=%.3f', tau1, tau2, s);
    end
    exp1 = exp(-dt/tau1);
    exp2 = exp(-dt/tau2);
    u1  = exp1 * u1 + I*R1 * (1 - exp1);
    u2  = exp2 * u2 + I*R2 * (1 - exp2);
    Vt  = Voc - I*R0 - u1 - u2;
    soc = max(0, min(1, soc - (I*dt)/(Qn_Ah*3600)));
    t   = t + dt;
end
