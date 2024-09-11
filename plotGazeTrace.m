function plotGazeTrace(CurrentTrial, TrialEnter, TrialExit, dx, dy, fractal_win, fixation_x, fixation_y, fixation_w, fixation_h)
    % Plot the gaze trace during the main part of the trial

    % -------------------------------------------------------------------------
    %  Made By Mohamad Hosein Faramarzi - 2024
    %-------------------------------------------------------------------------- 

    % Define the 8 points of 8 regions in a matrix
    regions = [768, 408;
              662, 662;
              408, 768;
              152, 662;
              48,  408;
              153, 153;
              407, 48;
              662, 153];
    % Find indices where time is between TrialEnter and TrialExit
    indices_in_range = find(CurrentTrial.time >= TrialEnter & CurrentTrial.time <= TrialExit);

    % Extract the gx and gy values in that time range
    gx_in_range = CurrentTrial.gx(indices_in_range);
    gy_in_range = CurrentTrial.gy(indices_in_range);
    time_in_range = CurrentTrial.time(indices_in_range);

    % Normalize time to use as color values (between 0 and 1)
    time_normalized = (time_in_range - min(time_in_range)) / (max(time_in_range) - min(time_in_range));

    % Use a smoother colormap (e.g., parula or viridis)
    colors = parula(length(gx_in_range));

    % Plot the gaze trace with color gradient
    hold on;

    % Adjust the gaze positions (for calibration purposes)
    gx_in_range = gx_in_range - dx;
    gy_in_range = gy_in_range - dy;
        % Define the fixation window rectangle and plot it
    hold on;
    rectangle('Position', [fixation_x, fixation_y, fixation_w, fixation_h], 'FaceColor', 'r', 'LineWidth', 2, 'LineStyle', '-');
    axis equal;


    % Set axis labels
    xlabel('X Position', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('Y Position', 'FontSize', 12, 'FontWeight', 'bold');
    title('Fixation Window', 'FontSize', 14, 'FontWeight', 'bold');
    grid on;
    hold on;

    % Plot each segment of the gaze trace with corresponding color and larger markers
    scatter(gx_in_range, gy_in_range, 50, colors, 'filled', 'MarkerEdgeColor', 'k');

    % Adjust axis limits dynamically
    xlim([min(gx_in_range) - 50, max(gx_in_range) + 50]);
    ylim([min(gy_in_range) - 50, max(gy_in_range) + 50]);

    % Set axis labels and title
    xlabel('Gaze X Position', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('Gaze Y Position', 'FontSize', 12, 'FontWeight', 'bold');
    title('Gaze Trace from TrialEnter to TrialExit', 'FontSize', 14, 'FontWeight', 'bold');

    % Add colorbar to show the time progression with a label
    c = colorbar;
    c.Label.String = 'Normalized Time';
    c.Label.FontSize = 12;

    % Add gridlines for better readability
    grid on;

    % Enhance figure aesthetics
    set(gca, 'FontSize', 10, 'FontWeight', 'bold');
    axis equal; % Keep the aspect ratio equal

    % Draw the fractal window with improved styling
    if(CurrentTrial.TrialType==0)
        x1 = regions(CurrentTrial.FractalRegions + 1, 1);
        y1 = regions(CurrentTrial.FractalRegions + 1, 2);
        if strcmp(CurrentTrial.FractalType, 'Good')
            rectangle('Position', [x1, y1, fractal_win, fractal_win], 'EdgeColor', 'g', 'LineWidth', 2, 'LineStyle', '--');
        else
            rectangle('Position', [x1, y1, fractal_win, fractal_win], 'EdgeColor', 'r', 'LineWidth', 2, 'LineStyle', '--');
        end
    elseif(CurrentTrial.TrialType==1)
        x1 = regions(CurrentTrial.FractalRegions(1)+1,1);
        y1 = regions(CurrentTrial.FractalRegions(1)+1,2);
        x2 = regions(CurrentTrial.FractalRegions(2)+1,1);
        y2 = regions(CurrentTrial.FractalRegions(2)+1,2);
        rectangle('Position', [x1, y1, fractal_win, fractal_win], 'EdgeColor', 'g', 'LineWidth', 2, 'LineStyle', '--');
        rectangle('Position', [x2, y2, fractal_win, fractal_win], 'EdgeColor', 'r', 'LineWidth', 2, 'LineStyle', '--');

    end

    xlim([0, 1000]);
    ylim([0, 1000]);


    hold off;
end
