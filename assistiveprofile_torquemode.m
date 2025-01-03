% Inputs
start_time = 55;          % Start time as a percentage of the gait cycle
peak_time = start_time + 12.5;           % Peak time as a percentage of the gait cycle
end_time = start_time + 35;            % End time as a percentage of the gait cycle
peak_torque = 15;         % Peak torque value

% Ensure that times are within valid bounds
start_time = max(1, min(100, start_time));
peak_time = max(start_time, min(100, peak_time));  % Ensure peak_time is after start_time
end_time = max(peak_time, min(100, end_time));     % Ensure end_time is after peak_time

% Generate rising segment (start_time to peak_time)
rising_duration = peak_time - start_time + 1;      % Duration of rising segment
rising_ang = linspace(0, pi/2, rising_duration);   % Rising segment angle
rising_torque = sin(rising_ang) * peak_torque;     % Rising torque profile

% Generate falling segment (peak_time to end_time)
falling_duration = end_time - peak_time + 1;       % Duration of falling segment
falling_ang = linspace(pi/2, pi, falling_duration);% Falling segment angle
falling_torque = sin(falling_ang) * peak_torque;   % Falling torque profile

% Combine torque profile
torqueAssistive = zeros(1, 100);
torqueAssistive(start_time:peak_time) = rising_torque;
torqueAssistive(peak_time:end_time) = falling_torque;

% Plot torque profile for visualization
figure;
plot(1:100, torqueAssistive, 'LineWidth', 2);
xlabel('Gait Cycle (%)');
ylabel('Torque (Nm)');
title('Assistive Torque Profile');
grid on;
