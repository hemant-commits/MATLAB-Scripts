% Capacity_OCV_Simple.m
% Clean, client-ready plots for 9 Excel logs (30/35/40 °C × 3 files)
% - Multi-select Excel files (uigetfile)
% - Sheet autodetect ("Detail", "Statis")
% - Robust time handling (RelTime resets per step -> cumulative across test)
% - Capacity% from per-step coulomb counting (Statis capacity if present)
% - Simple styling:
%     Color = Temperature (30/35/40 °C)
%     Style = solid (Discharge), dashed (Charge)
%     Width = high-rate (bold), low-rate (thin)
% - Separate OCV-like (~0.03C) plots (Charge/Discharge)
% - Time-domain plots with side-by-side High-rate vs Low-rate
% - NO file saving

clear; clc; close all;

%% Select Excel files
[files, path] = uigetfile('*.xlsx', 'Select Excel file(s)', 'MultiSelect','on');
if isequal(files,0), disp('No files selected.'); return; end
if ischar(files), files = {files}; end

%% Constants / Style
nominalCap   = 15;               % Ah
I_high_thr   = 5.0;              % A  -> "high-rate"
lowI_bandA   = [0.35 0.65];      % A  -> ~0.03C windows
tempsOrder   = [30 35 40];       % legend order

% Temperature colors
tempCol = containers.Map('KeyType','double','ValueType','any');
tempCol(30) = [0.26 0.54 0.96];  % blue
tempCol(35) = [1.00 0.58 0.00];  % orange
tempCol(40) = [0.89 0.10 0.11];  % red
defaultCol = [0.35 0.35 0.35];

LW_high = 2.2; LW_low = 1.3;     % bold vs thin

%% Aggregators
traceDis = {}; traceChg = {};           % Capacity% traces
ocvChg   = {}; ocvDis  = {};            % ~0.03C Capacity% traces
fileTab  = struct('D',{},'tempC',{});   % keep D table & temp for time plots

%% Loop files
for f = 1:numel(files)
    fname = files{f};
    fpath = fullfile(path, fname);
    fprintf('[%d/%d] %s\n', f, numel(files), fname);

    % Temperature from filename
    tok = regexp(fname, '(\d+(?:\.\d+)?)\s*(?:degC|°C|C)', 'tokens', 'once', 'ignorecase');
    if ~isempty(tok), tempC = str2double(tok{1}); else, tempC = NaN; end

    % Detect sheets
    try, [~, sheets] = xlsfinfo(fpath); catch, sheets = {}; end
    detSheet = 'Detail_96_1_4'; staSheet = 'Statis_96_1_4';

    try
        rawD = readtable(fpath, 'Sheet', detSheet, 'ReadVariableNames', false);
    catch
        cand = sheets(contains(lower(sheets),'detail'));
        if isempty(cand), error('Detail sheet not found in %s', fname); end
        detSheet = cand{1};
        rawD = readtable(fpath, 'Sheet', detSheet, 'ReadVariableNames', false);
    end

    rawS = table();
    try
        rawS = readtable(fpath, 'Sheet', staSheet, 'ReadVariableNames', false);
    catch
        candS = sheets(contains(lower(sheets),'statis'));
        if ~isempty(candS)
            staSheet = candS{1};
            rawS = readtable(fpath, 'Sheet', staSheet, 'ReadVariableNames', false);
        end
    end

    % Standardize Detail columns and numeric conversion
    rawD.Properties.VariableNames = {'RecordNum','State','Jump','Cycle','Steps', ...
        'Current_A','Voltage_V','Capacity_Ah','Energy_Wh','RelTime_str','Date_str'};
    D = rawD;

    numCols = {'RecordNum','Jump','Cycle','Steps','Current_A','Voltage_V','Capacity_Ah','Energy_Wh'};
    for c = 1:numel(numCols)
        cn = numCols{c};
        if ~ismember(cn, D.Properties.VariableNames), continue; end
        v = D.(cn);
        if iscell(v) || isstring(v), v = str2double(string(v)); end
        if ~isnumeric(v), v = double(v); end
        v(~isfinite(v)) = 0;
        D.(cn) = v;
    end

    % Keep Cycle 1
    if ismember('Cycle', D.Properties.VariableNames)
        D = D(D.Cycle==1, :);
    end
    if isempty(D), warning('No Cycle 1 rows in %s', fname); continue; end

    % Parse RelTime -> seconds (robust)
    relStr = string(D.RelTime_str);
    try
        dtmp = duration(relStr, 'InputFormat','hh:mm:ss.SSS');
    catch
        try
            dtmp = duration(relStr, 'InputFormat','hh:mm:ss');
        catch
            % manual hh:mm:ss(.fff) parsing
            mm = regexp(relStr, '^\s*(\d+):(\d+):(\d+(?:[.,]\d+)?)', 'tokens', 'once');
            svec = zeros(size(relStr));
            for k=1:numel(relStr)
                if ~isempty(mm{k})
                    h = str2double(mm{k}{1});
                    m = str2double(mm{k}{2});
                    s = str2double(strrep(mm{k}{3},',','.'));
                    svec(k) = h*3600 + m*60 + s;
                end
            end
            dtmp = seconds(svec);
        end
    end
    D.Time_s_rel = seconds(dtmp);
    D.Time_s_rel(~isfinite(D.Time_s_rel)) = 0;

    % Build cumulative time across steps (RelTime resets each step)
    steps = double(D.Steps);
    rel   = D.Time_s_rel;
    tCum  = zeros(size(rel));
    if ~isempty(rel)
        [uSteps,~,ix] = unique(steps,'stable');
        stepStart = accumarray(ix, rel, [], @(v) v(1));
        stepEnd   = accumarray(ix, rel, [], @(v) v(end));
        stepDur   = max(stepEnd - stepStart, 0);              % seconds
        stepOff   = [0; cumsum(stepDur(1:end-1))];            % seconds per unique step
        offByRow  = stepOff(ix);
        tCum      = offByRow + (rel - stepStart(ix));
    end
    D.Time_s_cum = tCum;  % seconds

    % Build step capacity map from Statis if available
    stepCaps = containers.Map('KeyType','double','ValueType','double');
    if ~isempty(rawS) && width(rawS) >= 19
        rawS.Properties.VariableNames = {'Channel','Cycle','Steps','OriginalStep','State', ...
            'OnsetV','EndV','StartI','EndI','Capacity_Ah','ContTime','RelTime','Date', ...
            'NetDChgCap','ChgCap','DChgCap','NetDChgEnergy','ChgEnergy','DChgEnergy'};
        if iscell(rawS.Cycle), rawS.Cycle = str2double(rawS.Cycle); end
        if iscell(rawS.Steps), rawS.Steps = str2double(rawS.Steps); end
        S1 = rawS(rawS.Cycle==1, :);
        caps = S1.Capacity_Ah; if iscell(caps), caps = str2double(caps); end
        for j = 1:height(S1)
            sNum = S1.Steps(j);
            if ~isnan(sNum) && ~isnan(caps(j)), stepCaps(sNum) = abs(caps(j)); end
        end
    end

    % Build Capacity% traces per step (Discharge / Charge)
    st = string(D.State);
    isD = contains(st,'DChg');
    isC = contains(st,'Chg') & ~isD;
    Iabs = abs(D.Current_A);

    % Discharge steps
    uD = unique(D.Steps(isD));
    for k = 1:numel(uD)
        idx = isD & (D.Steps == uD(k));
        if ~any(idx), continue; end
        dt_h_rel = [0; diff(D.Time_s_rel(idx))]/3600;
        q_step   = cumsum(abs(D.Current_A(idx)).*dt_h_rel);   % Ah
        if ~isempty(stepCaps) && isKey(stepCaps, uD(k)), q_tot = stepCaps(uD(k));
        else, q_tot = max(q_step); end
        if ~(q_tot > 0), q_tot = 1; end
        tr.tempC  = tempC;
        tr.capPct = (q_step / q_tot) * 100;
        tr.V      = D.Voltage_V(idx);
        tr.isHigh = any(Iabs(idx) > I_high_thr);
        traceDis{end+1} = tr; %#ok<SAGROW>
    end

    % Charge steps
    uC = unique(D.Steps(isC));
    for k = 1:numel(uC)
        idx = isC & (D.Steps == uC(k));
        if ~any(idx), continue; end
        dt_h_rel = [0; diff(D.Time_s_rel(idx))]/3600;
        q_step   = cumsum(abs(D.Current_A(idx)).*dt_h_rel);   % Ah
        if ~isempty(stepCaps) && isKey(stepCaps, uC(k)), q_tot = stepCaps(uC(k));
        else, q_tot = max(q_step); end
        if ~(q_tot > 0), q_tot = 1; end
        tr.tempC  = tempC;
        tr.capPct = (q_step / q_tot) * 100;
        tr.V      = D.Voltage_V(idx);
        tr.isHigh = any(Iabs(idx) > I_high_thr);
        traceChg{end+1} = tr; %#ok<SAGROW>
    end

    % OCV-like (~0.03C) traces (Discharge & Charge)
    lowI = Iabs >= lowI_bandA(1) & Iabs <= lowI_bandA(2);

    uDlo = unique(D.Steps(isD & lowI));
    for k = 1:numel(uDlo)
        idx = isD & lowI & (D.Steps == uDlo(k));
        if ~any(idx), continue; end
        dt_h_rel = [0; diff(D.Time_s_rel(idx))]/3600;
        q_step   = cumsum(abs(D.Current_A(idx)).*dt_h_rel);
        if ~isempty(stepCaps) && isKey(stepCaps, uDlo(k)), q_tot = stepCaps(uDlo(k));
        else, q_tot = max(q_step); end
        if ~(q_tot > 0), q_tot = 1; end
        tr.tempC  = tempC;
        tr.capPct = (q_step / q_tot) * 100;
        tr.V      = D.Voltage_V(idx);
        ocvDis{end+1} = tr; %#ok<SAGROW>
    end

    uClo = unique(D.Steps(isC & lowI));
    for k = 1:numel(uClo)
        idx = isC & lowI & (D.Steps == uClo(k));
        if ~any(idx), continue; end
        dt_h_rel = [0; diff(D.Time_s_rel(idx))]/3600;
        q_step   = cumsum(abs(D.Current_A(idx)).*dt_h_rel);
        if ~isempty(stepCaps) && isKey(stepCaps, uClo(k)), q_tot = stepCaps(uClo(k));
        else, q_tot = max(q_step); end
        if ~(q_tot > 0), q_tot = 1; end
        tr.tempC  = tempC;
        tr.capPct = (q_step / q_tot) * 100;
        tr.V      = D.Voltage_V(idx);
        ocvChg{end+1} = tr; %#ok<SAGROW>
    end

    % Keep table for time-domain plots
    fileTab(end+1).D = D; %#ok<SAGROW>
    fileTab(end).tempC = tempC;
end

%% ---------- FIGURES (simple, consistent) ----------

% Helper inline pick color
pickCol = @(t) ( ~isnan(t) && isKey(tempCol,t) ) .* tempCol(t) + ...
               ( isnan(t) || ~isKey(tempCol,t) ) .* defaultCol;

%% 1) DISCHARGE: Voltage vs Capacity%
figure('Name','Discharge: Voltage vs Capacity%','Color','w','Position',[60 80 1250 720]);
ax = axes; hold on; grid on; grid minor; box on;
for i = 1:numel(traceDis)
    tr = traceDis{i};
    c  = pickCol(tr.tempC);
    lw = LW_high; if ~tr.isHigh, lw = LW_low; end
    plot(tr.capPct, tr.V, '-', 'Color', c, 'LineWidth', lw);
end
xlabel('Capacity (%)'); ylabel('Voltage (V)');
title('Discharge: Voltage vs Capacity (%)');
subtitle('Color = Temp (30/35/40 °C) | Solid = Discharge | Width: highC (bold), lowC (thin)');
xlim([0 100]); ylim([1.8 3.7]);

% Legend proxies: per Temp, high/low
hLeg = []; labLeg = {};
for t = tempsOrder
    hasT = any(arrayfun(@(k) isequal(traceDis{k}.tempC,t), 1:numel(traceDis)));
    if ~hasT, continue; end
    c = tempCol(t);
    % highC
    hasHigh = any(arrayfun(@(k) isequal(traceDis{k}.tempC,t) && traceDis{k}.isHigh, 1:numel(traceDis)));
    if hasHigh
        hLeg(end+1) = plot(nan,nan,'-','Color',c,'LineWidth',LW_high); %#ok<SAGROW>
        labLeg{end+1} = sprintf('%g^\\circC — highC (DChg)', t);       %#ok<SAGROW>
    end
    % lowC
    hasLow = any(arrayfun(@(k) isequal(traceDis{k}.tempC,t) && ~traceDis{k}.isHigh, 1:numel(traceDis)));
    if hasLow
        hLeg(end+1) = plot(nan,nan,'-','Color',c,'LineWidth',LW_low);  %#ok<SAGROW>
        labLeg{end+1} = sprintf('%g^\\circC — lowC (DChg)', t);        %#ok<SAGROW>
    end
end
legend(hLeg, labLeg, 'Location','eastoutside');

%% 2) CHARGE: Voltage vs Capacity%
figure('Name','Charge: Voltage vs Capacity%','Color','w','Position',[60 80 1250 720]);
ax = axes; hold on; grid on; grid minor; box on;
for i = 1:numel(traceChg)
    tr = traceChg{i};
    c  = pickCol(tr.tempC);
    lw = LW_high; if ~tr.isHigh, lw = LW_low; end
    plot(tr.capPct, tr.V, '--', 'Color', c, 'LineWidth', lw);
end
xlabel('Capacity (%)'); ylabel('Voltage (V)');
title('Charge: Voltage vs Capacity (%)');
subtitle('Color = Temp (30/35/40 °C) | Dashed = Charge | Width: highC (bold), lowC (thin)');
xlim([0 100]); ylim([1.8 3.7]);

% Legend proxies: per Temp, high/low
hLeg = []; labLeg = {};
for t = tempsOrder
    hasT = any(arrayfun(@(k) isequal(traceChg{k}.tempC,t), 1:numel(traceChg)));
    if ~hasT, continue; end
    c = tempCol(t);
    % highC
    hasHigh = any(arrayfun(@(k) isequal(traceChg{k}.tempC,t) && traceChg{k}.isHigh, 1:numel(traceChg)));
    if hasHigh
        hLeg(end+1) = plot(nan,nan,'--','Color',c,'LineWidth',LW_high); %#ok<SAGROW>
        labLeg{end+1} = sprintf('%g^\\circC — highC (Chg)', t);         %#ok<SAGROW>
    end
    % lowC
    hasLow = any(arrayfun(@(k) isequal(traceChg{k}.tempC,t) && ~traceChg{k}.isHigh, 1:numel(traceChg)));
    if hasLow
        hLeg(end+1) = plot(nan,nan,'--','Color',c,'LineWidth',LW_low);  %#ok<SAGROW>
        labLeg{end+1} = sprintf('%g^\\circC — lowC (Chg)', t);          %#ok<SAGROW>
    end
end
legend(hLeg, labLeg, 'Location','eastoutside');

%% 3) OCV-like (~0.03C) CHARGE
figure('Name','OCV-like Charge (~0.03C)','Color','w','Position',[60 80 1200 700]);
axes; hold on; grid on; grid minor; box on;
for i = 1:numel(ocvChg)
    tr = ocvChg{i};
    c  = pickCol(tr.tempC);
    plot(tr.capPct, tr.V, '--', 'Color', c, 'LineWidth', 2.0);
end
xlabel('Capacity (%)'); ylabel('Voltage (V)');
title('OCV-like Charge (~0.03C)'); subtitle('Color = Temp (30/35/40 °C)');
xlim([0 100]); ylim([1.8 3.7]);
% legend: temps only
hLeg = []; labLeg = {};
for t = tempsOrder
    if any(arrayfun(@(k) isequal(ocvChg{k}.tempC,t), 1:numel(ocvChg)))
        hLeg(end+1) = plot(nan,nan,'--','Color',tempCol(t),'LineWidth',2); %#ok<SAGROW>
        labLeg{end+1} = sprintf('%g^\\circC', t);                           %#ok<SAGROW>
    end
end
legend(hLeg, labLeg, 'Location','eastoutside');

%% 4) OCV-like (~0.03C) DISCHARGE
figure('Name','OCV-like Discharge (~0.03C)','Color','w','Position',[60 80 1200 700]);
axes; hold on; grid on; grid minor; box on;
for i = 1:numel(ocvDis)
    tr = ocvDis{i};
    c  = pickCol(tr.tempC);
    plot(tr.capPct, tr.V, '-', 'Color', c, 'LineWidth', 2.0);
end
xlabel('Capacity (%)'); ylabel('Voltage (V)');
title('OCV-like Discharge (~0.03C)'); subtitle('Color = Temp (30/35/40 °C)');
xlim([0 100]); ylim([1.8 3.7]);
% legend: temps only
hLeg = []; labLeg = {};
for t = tempsOrder
    if any(arrayfun(@(k) isequal(ocvDis{k}.tempC,t), 1:numel(ocvDis)))
        hLeg(end+1) = plot(nan,nan,'-','Color',tempCol(t),'LineWidth',2); %#ok<SAGROW>
        labLeg{end+1} = sprintf('%g^\\circC', t);                          %#ok<SAGROW>
    end
end
legend(hLeg, labLeg, 'Location','eastoutside');

%% 5) Voltage vs Time — DISCHARGE (High-rate | Low-rate)
figure('Name','Voltage vs Time — Discharge','Color','w','Position',[60 80 1400 700]);
for p = 1:2
    subplot(1,2,p); hold on; grid on; grid minor; box on;
    for f = 1:numel(fileTab)
        D = fileTab(f).D; tC = fileTab(f).tempC;
        st = string(D.State); isD = contains(st,'DChg');
        if p==1, idx = isD & (abs(D.Current_A) > I_high_thr); ttl = 'Discharge — High-rate only';
        else,    idx = isD & (abs(D.Current_A) >= lowI_bandA(1) & abs(D.Current_A) <= lowI_bandA(2));
                 ttl = 'Discharge — Low-rate (~0.03C) only';
        end
        if any(idx)
            c = pickCol(tC);
            plot(D.Time_s_cum(idx)/3600, D.Voltage_V(idx), '-', 'Color', c, 'LineWidth', 1.8);
        end
    end
    xlabel('Time (h)'); ylabel('Voltage (V)'); ylim([1.8 3.7]); title(ttl);
end
% legend (temps only)
axes('Position',[0 0 1 1],'Visible','off'); hold on;
hLeg = []; labLeg = {};
for t = tempsOrder
    hLeg(end+1) = plot(nan,nan,'-','Color',tempCol(t),'LineWidth',2); %#ok<AGROW>
    labLeg{end+1} = sprintf('%g^\\circC', t);                           %#ok<AGROW>
end
legend(hLeg, labLeg, 'Position',[0.90 0.35 0.08 0.25]);  % right edge

%% 6) Voltage vs Time — CHARGE (High-rate | Low-rate)
figure('Name','Voltage vs Time — Charge','Color','w','Position',[60 80 1400 700]);
for p = 1:2
    subplot(1,2,p); hold on; grid on; grid minor; box on;
    for f = 1:numel(fileTab)
        D = fileTab(f).D; tC = fileTab(f).tempC;
        st = string(D.State); isC = contains(st,'Chg') & ~contains(st,'DChg');
        if p==1, idx = isC & (abs(D.Current_A) > I_high_thr); ttl = 'Charge — High-rate only';
        else,    idx = isC & (abs(D.Current_A) >= lowI_bandA(1) & abs(D.Current_A) <= lowI_bandA(2));
                 ttl = 'Charge — Low-rate (~0.03C) only';
        end
        if any(idx)
            c = pickCol(tC);
            plot(D.Time_s_cum(idx)/3600, D.Voltage_V(idx), '--', 'Color', c, 'LineWidth', 1.8);
        end
    end
    xlabel('Time (h)'); ylabel('Voltage (V)'); ylim([1.8 3.7]); title(ttl);
end
axes('Position',[0 0 1 1],'Visible','off'); hold on;
hLeg = []; labLeg = {};
for t = tempsOrder
    hLeg(end+1) = plot(nan,nan,'--','Color',tempCol(t),'LineWidth',2); %#ok<AGROW>
    labLeg{end+1} = sprintf('%g^\\circC', t);                            %#ok<AGROW>
end
legend(hLeg, labLeg, 'Position',[0.90 0.35 0.08 0.25]);

%% 7) Current vs Time — DISCHARGE (High-rate | Low-rate)
figure('Name','Current vs Time — Discharge','Color','w','Position',[60 80 1400 700]);
for p = 1:2
    subplot(1,2,p); hold on; grid on; grid minor; box on;
    for f = 1:numel(fileTab)
        D = fileTab(f).D; tC = fileTab(f).tempC;
        st = string(D.State); isD = contains(st,'DChg');
        if p==1, idx = isD & (abs(D.Current_A) > I_high_thr); ttl = 'Discharge — High-rate only'; yl = [-16 2];
        else,    idx = isD & (abs(D.Current_A) >= lowI_bandA(1) & abs(D.Current_A) <= lowI_bandA(2));
                 ttl = 'Discharge — Low-rate (~0.03C) only'; yl = [-1 1];
        end
        if any(idx)
            c = pickCol(tC);
            plot(D.Time_s_cum(idx)/3600, D.Current_A(idx), '-', 'Color', c, 'LineWidth', 1.8);
        end
    end
    xlabel('Time (h)'); ylabel('Current (A)'); ylim(yl); title(ttl);
end
axes('Position',[0 0 1 1],'Visible','off'); hold on;
hLeg = []; labLeg = {};
for t = tempsOrder
    hLeg(end+1) = plot(nan,nan,'-','Color',tempCol(t),'LineWidth',2); %#ok<AGROW>
    labLeg{end+1} = sprintf('%g^\\circC', t);                           %#ok<AGROW>
end
legend(hLeg, labLeg, 'Position',[0.90 0.35 0.08 0.25]);

%% 8) Current vs Time — CHARGE (High-rate | Low-rate)
figure('Name','Current vs Time — Charge','Color','w','Position',[60 80 1400 700]);
for p = 1:2
    subplot(1,2,p); hold on; grid on; grid minor; box on;
    for f = 1:numel(fileTab)
        D = fileTab(f).D; tC = fileTab(f).tempC;
        st = string(D.State); isC = contains(st,'Chg') & ~contains(st,'DChg');
        if p==1, idx = isC & (abs(D.Current_A) > I_high_thr); ttl = 'Charge — High-rate only'; yl = [0 16];
        else,    idx = isC & (abs(D.Current_A) >= lowI_bandA(1) & abs(D.Current_A) <= lowI_bandA(2));
                 ttl = 'Charge — Low-rate (~0.03C) only'; yl = [0 1];
        end
        if any(idx)
            c = pickCol(tC);
            plot(D.Time_s_cum(idx)/3600, D.Current_A(idx), '--', 'Color', c, 'LineWidth', 1.8);
        end
    end
    xlabel('Time (h)'); ylabel('Current (A)'); ylim(yl); title(ttl);
end
axes('Position',[0 0 1 1],'Visible','off'); hold on;
hLeg = []; labLeg = {};
for t = tempsOrder
    hLeg(end+1) = plot(nan,nan,'--','Color',tempCol(t),'LineWidth',2); %#ok<AGROW>
    labLeg{end+1} = sprintf('%g^\\circC', t);                            %#ok<AGROW>
end
legend(hLeg, labLeg, 'Position',[0.90 0.35 0.08 0.25]);

disp('All figures generated.');
