%% Clear Commands
clear
close all
clc

%% Check and Load Data
folder = what('Experimental Data');
user_exp = input(' Name of the variables file *.mat (without extension) for the time history to be analyzed = ' , 's');
file_i = strcat(user_exp, '.mat');

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
plot_err = 1;

figure('Name', 'Time History', 'NumberTitle', 'off')
plot(time, acc)
rectangle('Edgecolor', 'r', 'LineWidth', 2, 'Position', ...
    [0, min(acc), window_beg - (plot_err / 2), max(acc) - min(acc)])    % red rectangle
rectangle('Edgecolor', 'g', 'LineWidth', 2, 'Position', ...
    [window_beg + (plot_err / 2), min(acc), window_end - window_beg - (plot_err / 2), max(acc) - min(acc)])    % green rectangle
rectangle('Edgecolor', 'b', 'LineWidth', 2, 'Position', ...
    [window_end + plot_err, min(acc), max(time) - window_end - plot_err, max(acc) - min(acc)])    % blue rectangle
grid on
title('Time History')
xlabel('Time [s]')
ylabel('Acceleration [m/s^2]')

%% Peaks Extraction
[pks_acc, idx] = findpeaks(acc - mean(acc));  % peaks in complete signal (we subtract the mean)
pks_time = time(idx);   % peaks time axis

%% Main Loop
npt = 20;   % number of points considered in moving average
every_other = 8;    % number of values 'jumped' for evaluation

adim_damp_avg_values = zeros(1, every_other);   % initiate array of dimensionless averaged damping
adim_damp_values = zeros(1, every_other);   % initiate array of dimensionless  damping
eigenfreq_values = zeros(1, every_other);   % initiate array of first eigenfrequencies

for jj = 1 : every_other    % main loop that iterates over different value of "every_other"
    
    % Dimensionless Damping Calculation
    log_decs = zeros(1, (length(pks_acc) - jj));   % initiate logarithmic decrements array
    
    for ii = 1 : (length(pks_acc) - jj)    % calculate logarithmic decrements
        
        log_decs(ii) = (1 / jj) * log(pks_acc(ii) / pks_acc(ii + jj));
        
    end
    
    adim_damp = log_decs ./ (2 * pi);    % dimensionless damping array
    adim_damp_avg = movmean(adim_damp, npt);    % moving average of dimensionless damping
    adim_damp_time = pks_time(1 : (end - jj));   % dimensionless damping time array
    
    % Oscillation Period Calculation
    osc_period = zeros(1, (length(pks_acc) - jj)); % initiate oscillation period array
    
    for ii = 1 : (length(pks_acc) - jj)
        
        osc_period(ii) = (1 / jj) * (pks_time(ii + jj) - pks_time(ii));
        
    end
    
    eigenfreq = 1 ./ osc_period;    % array of first eigenfrequencies
    eigenfreq_time = pks_time(1 : (end - jj));    % time array of first eigenfrequencies
    
    % Final Parameters Calculation
    adim_damp_avg_values(jj) = mean(adim_damp_avg(adim_damp_time >= 15 & adim_damp_time <= 35));    % fill array of dimensionless damping (averaged)
    adim_damp_values(jj) = mean(adim_damp(adim_damp_time >= 15 & adim_damp_time <= 35));    % fill array of dimensionless damping
    eigenfreq_values(jj) = mean(eigenfreq(eigenfreq_time >= 15 & eigenfreq_time <= 35));    % fill array of eigenfrequencies
    
end

%% Plot Example of Extracted Data
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

figure('Name', 'First Eigenfrequencies', 'NumberTitle', 'off')  % eigenfrequencies plot
scatter(eigenfreq_time, eigenfreq, 30, 'b')
title('Eigenfrequency')
xlabel('Time [s]')
ylabel('Eigenfrequency [Hz]')

%% Display Final Parameters Vectors
disp('Dimensionless Dampings (averaged):')  % averaged dimensionless damping prints
disp(adim_damp_avg_values)
disp(mean(adim_damp_avg_values))
disp('Dimensionless Dampings:') % dimensionless damping prints
disp(adim_damp_values)
disp(mean(adim_damp_values))
disp('First Eigenfrequencies:') % first eigenfrequencies prints
disp(eigenfreq_values)
disp(mean(eigenfreq_values))

save(user_exp, 'adim_damp_avg_values', 'adim_damp_values', 'eigenfreq_values')  % save variables to .mat file

%% Plot of Linear Zone Time History
lin_beg = fsamp * window_beg;   % beginning of window in array index
lin_end = fsamp * window_end;   % end of window in array index
lin_time = time(lin_beg : lin_end) - window_beg;  % time axis (begins from 0)
lin_acc = acc(lin_beg : lin_end) - mean(acc);    % acc axis (we subtract the mean)

[lin_pks_acc, lin_idx] = findpeaks(lin_acc);    % data analysis function returns peaks value and peaks index
lin_pks_time = lin_time(lin_idx);    % linear zone peaks time axis

figure('Name', 'Linear Window', 'NumberTitle', 'off')
hold on
plot(lin_time, lin_acc)
scatter(lin_pks_time, lin_pks_acc)   % peaks plot
yline(0, '-g')  % horizontal green line
grid on
title('Linear Behaviour')
xlabel('Time [s]')
ylabel('Acceleration [m/s^2]')
