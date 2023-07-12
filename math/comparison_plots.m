function rmse_ekf = comparison_plots(simData,ekfData,meanMotion)
    [pitch, roll, yaw] = quat2angle(simData(2:5, 1:end)', 'YXZ');
    
    pitch_true = rad2deg(pitch);
    roll_true = rad2deg(roll);
    yaw_true = rad2deg(yaw);
    
    len = length(pitch);
    startIndex = round(len * 3 / 4);

    simData(6:8,1:end) = simData(6:8,1:end) *(180/pi); % [rad/sec] to [deg/sec]
    ekfData(6:8,1:end) = ekfData(6:8,1:end) *(180/pi);
    
    [pitch_est, roll_est, yaw_est] = quat2angle(ekfData(2:5, 1:end)', 'YXZ');
    
    pitch_est = rad2deg(pitch_est);
    roll_est = rad2deg(roll_est);
    yaw_est = rad2deg(yaw_est);
    
    
    rmseRoll = sqrt(immse(roll_est(startIndex:end), roll_true(startIndex:end)));
    rmsePitch = sqrt(immse(pitch_est(startIndex:end), pitch_true(startIndex:end)));
    rmseYaw = sqrt(immse(yaw_est(startIndex:end), yaw_true(startIndex:end)));
    rmseW1 = sqrt(immse(simData(6, startIndex:end), ekfData(6,startIndex:end)));
    rmseW2 = sqrt(immse(simData(7, startIndex:end) - meanMotion, ekfData(7,startIndex:end)));
    rmseW3 = sqrt(immse(simData(8, startIndex:end), ekfData(8,startIndex:end)));
    
    red = [203/255, 37/255, 37/255];
    green = [138/255, 181/255, 73/255];
    blue = [33/255, 144/255, 209/255];

    timeInHours = simData(1, :) / 3600;
    
    %% Euler Angles
    figure;
    plot(timeInHours, pitch_true, 'Color', red, 'LineWidth', 2)
    hold on
    plot(timeInHours, pitch_est, 'Color', green, 'LineWidth', 2)
    grid on
    xlabel('Time in hours')
    ylabel('Pitch Angle (\theta), [deg]')
    lgd = legend('True','Estimated',Location='southeast');
    lgd.Title.String = ['rmse_\theta = ' num2str(rmsePitch) 'deg'];
    title('Comparison of True & Estimated Pitch Angle');

    ax2 = axes('Position',[.6 .6 .3 .3]);
    box on;
    plot(timeInHours(timeInHours<=1), pitch_true(timeInHours<=1), 'Color', red, 'LineWidth', 2);
    hold on;
    plot(timeInHours(timeInHours<=1), pitch_est(timeInHours<=1), 'Color', green, 'LineWidth', 2);
    grid on;
    

    figure;
    plot(timeInHours, roll_true, 'Color', red, 'LineWidth', 2)
    hold on
    plot(timeInHours, roll_est, 'Color', green, 'LineWidth', 2)
    grid on
    xlabel('Time in hours')
    ylabel('Roll Angle (\phi), [deg]')
    lgd = legend('True','Estimated',Location='southeast');
    lgd.Title.String = ['rmse_\phi = ' num2str(rmseRoll) 'deg'];
    title('Comparison of True & Estimated Roll Angle');
    
    ax2 = axes('Position',[.6 .6 .3 .3]);
    box on;
    plot(timeInHours(timeInHours<=1), roll_true(timeInHours<=1), 'Color', red, 'LineWidth', 2);hold on;
    plot(timeInHours(timeInHours<=1), roll_est(timeInHours<=1), 'Color', green, 'LineWidth', 2); grid on;
    
    figure;
    plot(timeInHours, yaw_true, 'Color', red, 'LineWidth', 2)
    hold on
    plot(timeInHours, yaw_est, 'Color', green, 'LineWidth', 2)
    grid on
    xlabel('Time in hours')
    ylabel('Yaw Angle (\psi), [deg]')
    lgd = legend('True','Estimated',Location='southeast');
    lgd.Title.String = ['rmse_\psi = ' num2str(rmseYaw) 'deg'];
    title('Comparison of True & Estimated Yaw Angle');

    ax2 = axes('Position',[.6 .6 .3 .3]);
    box on;
    plot(timeInHours(timeInHours<=1), yaw_true(timeInHours<=1), 'Color', red, 'LineWidth', 2);hold on;
    plot(timeInHours(timeInHours<=1), yaw_est(timeInHours<=1), 'Color', green, 'LineWidth', 2); grid on;
    
    %% Angular Velocity
    figure;
    plot(timeInHours, simData(6,1:end), 'Color', red, 'LineWidth', 2)
    hold on
    plot(timeInHours, ekfData(6,1:end), 'Color', green, 'LineWidth', 2)
    grid on
    xlabel('Time in hours')
    ylabel('\omega_1, [deg/sec]')
    lgd = legend('True','Estimated',Location='southeast');
    lgd.Title.String = ['rmse_\omega_1 = ' num2str(rmseW1) 'deg/s'];
    title('Comparison of True & Estimated Angular Velocity (\omega_1)');
    ax2 = axes('Position',[.6 .6 .3 .3]);
    box on;
    plot(timeInHours(timeInHours<=1), simData(6,timeInHours<=1), 'Color', red, 'LineWidth', 2);hold on;
    plot(timeInHours(timeInHours<=1), ekfData(6,timeInHours<=1), 'Color', green, 'LineWidth', 2); grid on;
    
    figure;
    plot(timeInHours, simData(7,1:end), 'Color', red, 'LineWidth', 2)
    hold on
    plot(timeInHours, ekfData(7,1:end), 'Color', green, 'LineWidth', 2)
    grid on
    xlabel('Time in hours')
    ylabel('\omega_2, [deg/sec]')
    lgd = legend('True','Estimated',Location='southeast');
    lgd.Title.String = ['rmse_\omega_2 = ' num2str(rmseW2) 'deg/s'];
    title('Comparison of True & Estimated Angular Velocity (\omega_2)');
    ax2 = axes('Position',[.6 .6 .3 .3]);
    box on;
    plot(timeInHours(timeInHours<=1), simData(7,timeInHours<=1), 'Color', red, 'LineWidth', 2);hold on;
    plot(timeInHours(timeInHours<=1), ekfData(7,timeInHours<=1), 'Color', green, 'LineWidth', 2); grid on;
    
    figure;
    plot(timeInHours, simData(8,1:end), 'Color', red, 'LineWidth', 2)
    hold on
    plot(timeInHours, ekfData(8,1:end), 'Color', green, 'LineWidth', 2)
    grid on
    xlabel('Time in hours')
    ylabel('\omega_3, [deg/sec]')
    lgd = legend('True','Estimated',Location='southeast');
    lgd.Title.String = ['rmse_\omega_3 = ' num2str(rmseW2) 'deg/s'];
    title('Comparison of True & Estimated Angular Velocity (\omega_3)');
    ax2 = axes('Position',[.6 .6 .3 .3]);
    box on;
    plot(timeInHours(timeInHours<=1), simData(8,timeInHours<=1), 'Color', red, 'LineWidth', 2);hold on;
    plot(timeInHours(timeInHours<=1), ekfData(8,timeInHours<=1), 'Color', green, 'LineWidth', 2); grid on;

    figure;
    plot(timeInHours, ekfData(9,1:end),'Color', red, 'LineWidth', 2)
    hold on
    plot(timeInHours, ekfData(10,1:end), 'Color', green, 'LineWidth', 2)
    hold on
    plot(timeInHours, ekfData(11,1:end), 'Color', blue, 'LineWidth', 2)
    grid on
    xlabel('Time in hours')
    ylabel('Magnetometer Bias [T]')
    legend('b_x','b_y','b_z');
    title('Magnetometer Bias Estimation')

    figure;
    plot(timeInHours, ekfData(12,1:end),'Color', red, 'LineWidth', 2)
    hold on
    plot(timeInHours, ekfData(13,1:end), 'Color', green, 'LineWidth', 2)
    hold on
    plot(timeInHours, ekfData(14,1:end), 'Color', blue, 'LineWidth', 2)
    grid on
    xlabel('Time in hours')
    ylabel('Gyroscope Bias [rad/s]')
    legend('b_x','b_y','b_z');
    title('Gyroscope Bias Estimation')

    rmse_ekf = [rmseRoll rmsePitch rmseYaw rmseW1 rmseW2 rmseW3];

end
