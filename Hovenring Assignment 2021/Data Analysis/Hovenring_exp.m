%% Clear Commands
clear
close all
clc

%% Check and Load Data
folder = what('Experimental Data');
% user = input(' Name of the input file *.inp (without extension) for the structure to be analyzed = ' , 's');
user = 'Cable13'; % temporary
file_i = strcat(user, '.mat');

try
    
    load(strcat(folder.path, '\', file_i));
    
catch load_err    % error handling
    
    if strcmp(load_err.identifier, 'MATLAB:load:couldNotReadFile')
        
        error('ERROR: invalid file name')
        
    end
    
end

%% Plot Time History
window_beg = 15;  % beginning of linear behaviour
window_end = 35;  % end of appreciable TMD influence

figure('Name', 'Time History', 'NumberTitle', 'off')
plot(time, acc)
rectangle('Edgecolor', 'g', 'LineWidth', 2, 'Position', ...
    [window_beg, min(acc), window_end - window_beg, max(acc) - min(acc)])    % rectangle
grid on
title('Time History')
xlabel('Time [s]')
ylabel('Acceleration [m/s^2]')

%% Peaks Extraction
[pks_acc, idx] = findpeaks(acc - mean(acc));  % peaks in complete signal (we subtract the mean)
pks_time = zeros(1, length(pks_acc));    % initiate empty array

for ii = 1 : length(idx)    % peaks' time extraction in complete signal
    
   pks_time(ii) = time(idx(ii)); 
    
end

every_other = 4;    % number of values 'jumped' for evaluation

%% Dimensionless Damping Calculation
log_decs = zeros(1, (length(pks_acc) - every_other));   % initiate logarithmic decrements array

for ii = 1 : (length(pks_acc) - every_other)    % calculate logarithmic decrements
    
    log_decs(ii) = (1 / every_other) * (pks_acc(ii) / pks_acc(ii + every_other));
    
end

adim_damp = log_decs ./ (2 * pi);    % dimensionless damping array
adim_damp_time = pks_time(1 : (end - every_other));   % dimensionless damping time array

npt = 10;   % number of points considered in moving average
adim_damp_avg = movmean(adim_damp, npt);    % moving average of dimensionless damping

%% Oscillation Period Calculation
osc_period = zeros(1, (length(pks_acc) - every_other)); % initiate oscillation period array

for ii = 1 : (length(pks_acc) - every_other)
    
    osc_period(ii) = (1 / every_other) * (pks_time(ii + every_other) - pks_time(ii));
    
end

eigenfreq = 1 ./ osc_period;    % array of first eigenfrequencies
eigenfreq_time = pks_time(1 : (end - every_other));    % time array of first eigenfrequencies   

%% Plot Extracted Data
figure('Name', 'Non-dimensional Damping', 'NumberTitle', 'off') % dimensionless damping plot
scatter(adim_damp_time, adim_damp, 30, 'r')
grid on
title('Non-dimensional Damping')
xlabel('Time [s]')
ylabel('h [-]')

figure('Name', 'Non-dimensional Damping (Averaged)', 'NumberTitle', 'off')  % dimensionless damping plot (averaged)
scatter(adim_damp_time, adim_damp_avg, 30, 'r')
title('Non-dimensional Damping (Averaged)')
xlabel('Time [s]')
ylabel('h [-]')

figure('Name', 'First Eigenfrequencies', 'NumberTitle', 'off')  % dimensionless damping plot (averaged)
scatter(eigenfreq_time, eigenfreq, 30, 'b')
title('Eigenfrequency')
xlabel('Time [s]')
ylabel('Eigenfrequency [Hz]')

%% Parameters Calculation in Linear Zone
lin_beg = fsamp * window_beg;   % beginning of window in array index
lin_end = fsamp * window_end;   % end of window in array index
lin_time = time(lin_beg : lin_end) - window_beg;  % time axis (begins from 0)
lin_acc = acc(lin_beg : lin_end) - mean(acc);    % acc axis (we subtract the mean)

[lin_pks_acc, lin_idx] = findpeaks(lin_acc);    % data analysis function returns peaks value and peaks index
lin_pks_time = zeros(1, length(lin_pks_acc));    % initiate empty array of peaks' time

for ii = 1 : length(lin_idx)    % peak's time extraction in linear zone
    
    lin_pks_time(ii) = lin_time(lin_idx(ii));
    
end

lin_n_param = floor(length(lin_pks_acc) / every_other);
lin_log_decs = (1 / every_other) * log(lin_pks_acc(1 : every_other : ((lin_n_param - 1) ...
    * every_other)) ./ lin_pks_acc((1 + every_other) : every_other : ...
    (lin_n_param * every_other)));   % calculate logarithmic decrements

lin_adim_damp = lin_log_decs ./ (2 * pi);    % dimensionless damping array

lin_npt = 5;   % number of points considered in moving average
lin_adim_damp_avg = zeros(1, (length(lin_adim_damp) - npt));    % initiate moving average array

for ii = 1 : (length(lin_adim_damp) - lin_npt)    % moving average calculation
    
    lin_adim_damp_avg(ii) = (1 / lin_npt) * sum(lin_adim_damp(ii : (ii + lin_npt)));
    
end

disp('The dimensionless damping evaluated in the linear zone is:')
disp(mean(adim_damp_avg))








figure('Name', 'Linear Window', 'NumberTitle', 'off')
hold on
plot(lin_time, lin_acc)
scatter(lin_pks_time, lin_pks_acc)   % peaks' plot
yline(0, '-g')  % horizontal green line
grid on
title('Linear Behaviour')
xlabel('Time [s]')
ylabel('Acceleration [m/s^2]')
