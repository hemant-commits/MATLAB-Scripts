% plot_capacity_and_ocv_by_temperature.m
% Select one or more HPPC/OCV Excel files. The script:
%   • extracts the "Detail" sheet
%   • computes SOC for each data set (properly anchored to 0-100%)
%   • parses the temperature from the filename (e.g. "…_30degC_…", "…_15°C_…")
%   • shows one combined figure for all capacity-test charge/discharge curves
%   • shows a second combined figure for all OCV charge/discharge curves

clear; clc;

[fileNames, filePath] = uigetfile( ...
    {'*.xlsx;*.xls', 'Excel files (*.xlsx, *.xls)'}, ...
    'Select HPPC/OCV data files', 'MultiSelect', 'on');

if isequal(fileNames, 0)
    disp('No file selected. Exiting.');
    return;
end
if ischar(fileNames), fileNames = {fileNames}; end

% C-rate mapping for different temperatures
cRateMap = containers.Map();
cRateMap('30') = struct('charge', 0.8, 'discharge', 1.0);
cRateMap('35') = struct('charge', 0.5, 'discharge', 1.0);
cRateMap('40') = struct('charge', 0.33, 'discharge', 0.5);

% Containers for combined plotting
capChargeAll    = struct('label', {}, 'temp', {}, 'cRate', {}, 'soc', {}, 'volt', {});
capDischargeAll = struct('label', {}, 'temp', {}, 'cRate', {}, 'soc', {}, 'volt', {});
ocvChargeAll    = struct('label', {}, 'temp', {}, 'soc', {}, 'volt', {});
ocvDischargeAll = struct('label', {}, 'temp', {}, 'soc', {}, 'volt', {});
ocvWholeData    = struct('label', {}, 'temp', {}, 'time', {}, 'volt', {});
capWholeData    = struct('label', {}, 'temp', {}, 'chgCRate', {}, 'dchgCRate', {}, 'time', {}, 'volt', {});

for f = 1:numel(fileNames)
    fullFile = fullfile(filePath, fileNames{f});
    fprintf('\nProcessing %s\n', fullFile);

    % Identify the "Detail" sheet
    sheets = sheetnames(fullFile);
    idxDetail = find(contains(lower(sheets), 'detail'), 1, 'first');
    if isempty(idxDetail)
        warning('Skipping %s: no sheet containing the word \"Detail\".', fileNames{f});
        continue;
    end
    detailSheet = sheets{idxDetail};
    fprintf('Using sheet: %s\n', detailSheet);

    T = readtable(fullFile, 'Sheet', detailSheet);
    varNames = T.Properties.VariableNames;

    % Locate relevant columns (handles MATLAB's sanitized headings)
    findCol = @(token) varNames{find(contains(varNames, token, 'IgnoreCase', true), 1, 'first')};
    colState    = findCol('State');
    colSteps    = findCol('Steps');
    colVoltage  = findCol('Voltage');
    colCapacity = findCol('Capacity');

    requiredCols = {colState, colSteps, colVoltage, colCapacity};
    if any(cellfun(@isempty, requiredCols))
        warning('Skipping %s: missing expected columns.', fileNames{f});
        continue;
    end

    % Optional date column (for sorting)
    dateColIdx = find(contains(varNames, 'Date', 'IgnoreCase', true), 1, 'first');
    if ~isempty(dateColIdx)
        try
            T.DateTime = datetime(T.(varNames{dateColIdx}), 'InputFormat', 'yyyy-MM-dd HH:mm:ss');
        catch
            % fall back to MATLAB default parsing
            T.DateTime = datetime(T.(varNames{dateColIdx}));
        end
    else
        T.DateTime = (0:height(T)-1).';
    end
    T = sortrows(T, 'DateTime');

    State    = string(T.(colState));
    Steps    = T.(colSteps);
    Voltage  = T.(colVoltage);
    Capacity = T.(colCapacity);
    DateTime = T.DateTime;  % Keep datetime for time-based plotting

    % Determine reference capacity (largest excursion among main steps)
    refIdx = ismember(Steps, [3 5 9 11]);
    refCapacityAh = max(Capacity(refIdx));
    if refCapacityAh <= 0 || isnan(refCapacityAh)
        warning('Skipping %s: reference capacity <= 0.', fileNames{f});
        continue;
    end

    % Build SOC profile across steps with proper 0-100% anchoring
    soc = nan(height(T), 1);
    uniqueSteps = unique(Steps, 'stable');

    for k = 1:numel(uniqueSteps)
        stepVal = uniqueSteps(k);
        idxStep = Steps == stepVal;
        stateName = State(find(idxStep, 1, 'first'));

        if contains(stateName, 'DChg', 'IgnoreCase', true)
            % Discharge: SOC decreases from 100% to 0%
            deltaCap = Capacity(idxStep);
            % Start from 100% and subtract capacity
            soc(idxStep) = 100 - (deltaCap ./ refCapacityAh) * 100;
            
        elseif contains(stateName, 'Chg', 'IgnoreCase', true)
            % Charge: SOC increases from 0% to 100%
            deltaCap = Capacity(idxStep);
            % Start from 0% and add capacity
            soc(idxStep) = (deltaCap ./ refCapacityAh) * 100;
            
        else
            % Rest or other state: maintain constant SOC
            if k > 1
                % Use the last SOC value from previous step
                prevStepIdx = Steps == uniqueSteps(k-1);
                lastSOC = soc(find(prevStepIdx, 1, 'last'));
                if ~isnan(lastSOC)
                    soc(idxStep) = lastSOC;
                else
                    soc(idxStep) = 50; % default fallback
                end
            else
                soc(idxStep) = 50; % default for first step if it's a rest
            end
        end

        % Clamp SOC to [0, 100]
        soc(idxStep) = max(min(soc(idxStep), 100), 0);
    end

    % Pull temperature from filename (e.g., "…_30degC…", "…_15°C…")
    baseName = erase(fileNames{f}, {'.xlsx', '.xls'});
    lowerName = lower(baseName);
    tempTok = regexp(lowerName, '(-?\d+(\.\d+)?)\s*(deg|°)?c', 'tokens', 'once');
    if isempty(tempTok)
        tempLabel = 'Unknown';
        tempNum = 'Unknown';
    else
        tempLabel = tempTok{1};
        tempNum = tempLabel;
        if tempLabel(1) ~= '-'  % make sure there's no leading plus sign
            tempLabel = strtrim(tempLabel);
        end
    end
    
    % Get C-rates for this temperature
    if isKey(cRateMap, tempNum)
        cRates = cRateMap(tempNum);
        chgCRate = cRates.charge;
        dchgCRate = cRates.discharge;
    else
        chgCRate = NaN;
        dchgCRate = NaN;
    end
    
    label = sprintf('%s°C', tempLabel);

    % Helper to snag step data
    grab = @(stepList) struct( ...
        'soc',  soc(ismember(Steps, stepList)), ...
        'volt', Voltage(ismember(Steps, stepList)));
    
    % Helper for time-based data
    grabTime = @(stepList) struct( ...
        'time', DateTime(ismember(Steps, stepList)), ...
        'volt', Voltage(ismember(Steps, stepList)));

    capCharge    = grab(3);
    capDischarge = grab(5);
    ocvCharge    = grab(9);
    ocvDischarge = grab(11);
    
    % Get all OCV data (steps 9 and 11) for voltage vs time plot
    ocvWhole = grabTime([9 11]);
    % Get all Capacity data (steps 3 and 5) for voltage vs time plot
    capWhole = grabTime([3 5]);

    capChargeAll(end+1)    = struct('label', label, 'temp', str2double(tempNum), ...
        'cRate', chgCRate, 'soc', capCharge.soc, 'volt', capCharge.volt);
    capDischargeAll(end+1) = struct('label', label, 'temp', str2double(tempNum), ...
        'cRate', dchgCRate, 'soc', capDischarge.soc, 'volt', capDischarge.volt);
    ocvChargeAll(end+1)    = struct('label', label, 'temp', str2double(tempNum), ...
        'soc', ocvCharge.soc, 'volt', ocvCharge.volt);
    ocvDischargeAll(end+1) = struct('label', label, 'temp', str2double(tempNum), ...
        'soc', ocvDischarge.soc, 'volt', ocvDischarge.volt);
    ocvWholeData(end+1)    = struct('label', label, 'temp', str2double(tempNum), ...
        'time', ocvWhole.time, 'volt', ocvWhole.volt);
    capWholeData(end+1)    = struct('label', label, 'temp', str2double(tempNum), ...
        'chgCRate', chgCRate, 'dchgCRate', dchgCRate, 'time', capWhole.time, 'volt', capWhole.volt);
end

if isempty(capChargeAll)
    error('No valid files processed—nothing to plot.');
end

%% Combined capacity figure
figure('Name', 'Capacity Test - Voltage vs SOC at Different Temperatures', ...
    'Position', [100 100 1400 600], 'Color', 'w');

% Define professional color scheme - one color per temperature
colorMap = containers.Map();
colorMap('30') = [0.0000 0.4470 0.7410];  % Blue
colorMap('35') = [0.8500 0.3250 0.0980];  % Red-orange
colorMap('40') = [0.4660 0.6740 0.1880];  % Green
colorMap('25') = [0.9290 0.6940 0.1250];  % Yellow
colorMap('20') = [0.4940 0.1840 0.5560];  % Purple
colorMap('45') = [0.3010 0.7450 0.9330];  % Cyan
% Default color for unknown temperatures
defaultColor = [0.5 0.5 0.5];  % Gray

% Sort by temperature for consistent plotting
[~, sortIdx] = sort([capChargeAll.temp]);
capChargeAll = capChargeAll(sortIdx);
[~, sortIdx] = sort([capDischargeAll.temp]);
capDischargeAll = capDischargeAll(sortIdx);

% Charge subplot
subplot(1,2,1); hold on; grid on; box on;
for k = 1:numel(capChargeAll)
    tempStr = num2str(capChargeAll(k).temp);
    if isKey(colorMap, tempStr)
        plotColor = colorMap(tempStr);
    else
        plotColor = defaultColor;
    end
    
    if ~isnan(capChargeAll(k).cRate)
        legendStr = sprintf('%s (%.2gC)', capChargeAll(k).label, capChargeAll(k).cRate);
    else
        legendStr = capChargeAll(k).label;
    end
    plot(capChargeAll(k).soc, capChargeAll(k).volt, 'LineWidth', 2, ...
        'Color', plotColor, 'DisplayName', legendStr);
end
xlabel('State of Charge (%)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Voltage (V)', 'FontSize', 12, 'FontWeight', 'bold');
title('Charge Profile', 'FontSize', 14, 'FontWeight', 'bold');
legend('Location', 'southeast', 'FontSize', 10);
xlim([0 100]);
set(gca, 'FontSize', 11, 'LineWidth', 1.5);

% Discharge subplot
subplot(1,2,2); hold on; grid on; box on;
for k = 1:numel(capDischargeAll)
    tempStr = num2str(capDischargeAll(k).temp);
    if isKey(colorMap, tempStr)
        plotColor = colorMap(tempStr);
    else
        plotColor = defaultColor;
    end
    
    if ~isnan(capDischargeAll(k).cRate)
        legendStr = sprintf('%s (%.2gC)', capDischargeAll(k).label, capDischargeAll(k).cRate);
    else
        legendStr = capDischargeAll(k).label;
    end
    plot(capDischargeAll(k).soc, capDischargeAll(k).volt, 'LineWidth', 2, ...
        'Color', plotColor, 'DisplayName', legendStr);
end
set(gca, 'XDir', 'reverse');
xlabel('State of Charge (%)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Voltage (V)', 'FontSize', 12, 'FontWeight', 'bold');
title('Discharge Profile', 'FontSize', 14, 'FontWeight', 'bold');
legend('Location', 'southwest', 'FontSize', 10);
xlim([0 100]);
set(gca, 'FontSize', 11, 'LineWidth', 1.5);

% Add main title
sgtitle('Capacity Test: Voltage vs SOC at Different Temperatures and C-Rates', ...
    'FontSize', 16, 'FontWeight', 'bold');

%% Combined OCV figure
% OCV C-rate (C/30)
ocvCRate = 1/30; % C/30 rate

figure('Name', 'OCV Test - Voltage vs SOC at Different Temperatures', ...
    'Position', [150 150 1400 600], 'Color', 'w');

% Sort by temperature for consistent plotting
[~, sortIdx] = sort([ocvChargeAll.temp]);
ocvChargeAll = ocvChargeAll(sortIdx);
[~, sortIdx] = sort([ocvDischargeAll.temp]);
ocvDischargeAll = ocvDischargeAll(sortIdx);

% Charge subplot
subplot(1,2,1); hold on; grid on; box on;
for k = 1:numel(ocvChargeAll)
    tempStr = num2str(ocvChargeAll(k).temp);
    if isKey(colorMap, tempStr)
        plotColor = colorMap(tempStr);
    else
        plotColor = defaultColor;
    end
    
    plot(ocvChargeAll(k).soc, ocvChargeAll(k).volt, 'LineWidth', 2, ...
        'Color', plotColor, 'DisplayName', ocvChargeAll(k).label);
end
xlabel('State of Charge (%)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Voltage (V)', 'FontSize', 12, 'FontWeight', 'bold');
title('OCV Charge Profile', 'FontSize', 14, 'FontWeight', 'bold');
legend('Location', 'southeast', 'FontSize', 10);
xlim([0 100]);
set(gca, 'FontSize', 11, 'LineWidth', 1.5);

% Discharge subplot
subplot(1,2,2); hold on; grid on; box on;
for k = 1:numel(ocvDischargeAll)
    tempStr = num2str(ocvDischargeAll(k).temp);
    if isKey(colorMap, tempStr)
        plotColor = colorMap(tempStr);
    else
        plotColor = defaultColor;
    end
    
    plot(ocvDischargeAll(k).soc, ocvDischargeAll(k).volt, 'LineWidth', 2, ...
        'Color', plotColor, 'DisplayName', ocvDischargeAll(k).label);
end
set(gca, 'XDir', 'reverse');
xlabel('State of Charge (%)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Voltage (V)', 'FontSize', 12, 'FontWeight', 'bold');
title('OCV Discharge Profile', 'FontSize', 14, 'FontWeight', 'bold');
legend('Location', 'southwest', 'FontSize', 10);
xlim([0 100]);
set(gca, 'FontSize', 11, 'LineWidth', 1.5);

% Add main title with C-rate
sgtitle(sprintf('OCV Test at C/30 (%.4gC): Voltage vs SOC at Different Temperatures', ocvCRate), ...
    'FontSize', 16, 'FontWeight', 'bold');

%% OCV Voltage vs Time figure
figure('Name', 'OCV Test - Voltage vs Time', ...
    'Position', [200 200 1200 700], 'Color', 'w');

% Sort by temperature for consistent plotting
[~, sortIdx] = sort([ocvWholeData.temp]);
ocvWholeData = ocvWholeData(sortIdx);

hold on; grid on; box on;

for k = 1:numel(ocvWholeData)
    tempStr = num2str(ocvWholeData(k).temp);
    if isKey(colorMap, tempStr)
        plotColor = colorMap(tempStr);
    else
        plotColor = defaultColor;
    end
    
    % Convert time to relative time in hours from start of OCV test
    timeData = ocvWholeData(k).time;
    if ~isempty(timeData) && isdatetime(timeData)
        relativeTime = hours(timeData - timeData(1));
    else
        relativeTime = (1:length(ocvWholeData(k).volt))';
    end
    
    plot(relativeTime, ocvWholeData(k).volt, 'LineWidth', 2, ...
        'Color', plotColor, 'DisplayName', ocvWholeData(k).label);
end

xlabel('Time (hours)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Voltage (V)', 'FontSize', 12, 'FontWeight', 'bold');
title('OCV Test: Voltage vs Time', 'FontSize', 14, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 10);
set(gca, 'FontSize', 11, 'LineWidth', 1.5);

% Add subtitle with C-rate info
subtitle(sprintf('C-Rate: C/30 (%.4gC)', ocvCRate), 'FontSize', 12);

%% Capacity Test Voltage vs Time figure
figure('Name', 'Capacity Test - Voltage vs Time', ...
    'Position', [250 250 1200 700], 'Color', 'w');

% Sort by temperature for consistent plotting
[~, sortIdx] = sort([capWholeData.temp]);
capWholeData = capWholeData(sortIdx);

hold on; grid on; box on;

for k = 1:numel(capWholeData)
    tempStr = num2str(capWholeData(k).temp);
    if isKey(colorMap, tempStr)
        plotColor = colorMap(tempStr);
    else
        plotColor = defaultColor;
    end
    
    % Convert time to relative time in hours from start of capacity test
    timeData = capWholeData(k).time;
    if ~isempty(timeData) && isdatetime(timeData)
        relativeTime = hours(timeData - timeData(1));
    else
        relativeTime = (1:length(capWholeData(k).volt))';
    end
    
    % Create legend with both charge and discharge C-rates
    if ~isnan(capWholeData(k).chgCRate) && ~isnan(capWholeData(k).dchgCRate)
        legendStr = sprintf('%s (Chg: %.2gC, DChg: %.2gC)', ...
            capWholeData(k).label, capWholeData(k).chgCRate, capWholeData(k).dchgCRate);
    else
        legendStr = capWholeData(k).label;
    end
    
    plot(relativeTime, capWholeData(k).volt, 'LineWidth', 2, ...
        'Color', plotColor, 'DisplayName', legendStr);
end

xlabel('Time (hours)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Voltage (V)', 'FontSize', 12, 'FontWeight', 'bold');
title('Capacity Test: Voltage vs Time', 'FontSize', 14, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 10);
set(gca, 'FontSize', 11, 'LineWidth', 1.5);

% Add subtitle
subtitle('Complete Charge and Discharge Cycles at Different Temperatures', 'FontSize', 12);

fprintf('\nPlotted combined capacity and OCV curves for %d data set(s).\n', numel(capChargeAll));
