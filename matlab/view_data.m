% ---------------------------------------------------------------------------------
% BER plot
% ---------------------------------------------------------------------------------

figure
grid on  % Ensure grid is on
set(gca, 'YScale', 'log')  % Set the Y-axis to logarithmic scale
set(gca, 'YMinorGrid', 'on')  % Enable minor grid lines on Y-axis
set(gca, 'XMinorGrid', 'on')  % Enable minor grid lines on X-axis
set(gca, 'GridLineStyle', '-')  % Solid lines for grids
set(gca, 'MinorGridLineStyle', ':')  % Dotted lines for minor grids
set(gca, 'XMinorTick', 'on')  % Enable minor ticks on X-axis
set(gca, 'YMinorTick', 'on')  % Enable minor ticks on Y-axis
set(gca, 'FontSize', 12)  % Increase font size for better readability
title('Bit Error Rate vs. Signal-to-Noise Ratio')
xlabel('SNR (dB)')
ylabel('Bit Error Rate (BER)')

% Define a set of regular, easily distinguishable colors
colors = {[0, 0.4470, 0.7410], [0.8500, 0.3250, 0.0980], [0.9290, 0.6940, 0.1250], ...
          [0.4940, 0.1840, 0.5560], [0.4660, 0.6740, 0.1880], [0.3010, 0.7450, 0.9330], ...
          [0.6350, 0.0780, 0.1840]};
markers = {'o', 's', 'd', '^', '>', 'v', '<', 'p'};  % Different markers
lineStyles = {'-', '--', '-.', ':'};  % Different line styles

for i = 1:nb_algorithms
    mean_errors = zeros(1, length(snr));
    for j = 1:length(snr)
        mean_errors(j) = sum(squeeze(errors(i,j,:,:)), "all") / (Kb * cpd_L * frames);
    end

    % Plot each algorithm with a unique style
    semilogy(snr, mean_errors, 'LineWidth', 1.5, 'Color', colors{mod(i-1, length(colors)) + 1}, ...
             'Marker', markers{mod(i-1, length(markers)) + 1}, ...
             'LineStyle', lineStyles{mod(i-1, length(lineStyles)) + 1});
    hold on
end

legendInfo = cell(1, nb_algorithms);
for i = 1:nb_algorithms
    if i == 1
        legendInfo{i} = 'No decoding';
    end
    if i == 2
        legendInfo{i} = 'SGDBF, \epsilon = 0.001';
    end
    if i == 3
        legendInfo{i} = 'SGDBF, \epsilon = 0.1'; %'ISGDBF';
    end
    if i == 4
        legendInfo{i} = 'SGDBF, \epsilon = 1'; %'LC-HWBF, \theta_{att} = 0.9';
    end
    if i == 5
        legendInfo{i} = 'SGDBF, \epsilon = 2';
    end
    if i == 6
        legendInfo{i} = 'SGDBF, \epsilon = 5';
    end
    if i == 7
        legendInfo{i} = 'SGDBF, \epsilon = 10';
    end
    if i == 8
        legendInfo{i} = 'SGDBF, \epsilon = 0.5';
    end
    if i == 9
        legendInfo{i} = 'SGDBF, \epsilon = 0.9';
    end
    if i == 10
        legendInfo{i} = 'SGDBF, \epsilon = 0';
    end
    
    % print the total number of errors for each algorithm
    total_errors = sum(sum(squeeze(errors(i,:,:,:)),2),1);
    disp(['Total number of errors for algorithm ' num2str(i) ' is ' num2str(total_errors)]);
end

legend(legendInfo, 'Location', 'southwest')  % Position the legend in the bottom left corner
set(gca, 'YScale', 'log')  % Set the Y-axis to logarithmic scale
set(gca, 'GridLineStyle', '-')  % Solid lines for main grids
set(gca, 'MinorGridLineStyle', ':')  % Dotted lines for minor grids
title('Bit Error Rate vs. Signal-to-Noise Ratio')
xlabel('SNR (dB)')
ylabel('Bit Error Rate (BER)')
grid on 
hold off



% ---------------------------------------------------------------------------------
% BER plot
% ---------------------------------------------------------------------------------
%
%figure
%grid on  % Ensure grid is on
%set(gca, 'YScale', 'log')  % Set the Y-axis to logarithmic scale
%set(gca, 'YMinorGrid', 'on')  % Enable minor grid lines on Y-axis
%set(gca, 'XMinorGrid', 'on')  % Enable minor grid lines on X-axis
%set(gca, 'GridLineStyle', '-')  % Solid lines for grids
%set(gca, 'MinorGridLineStyle', ':')  % Dotted lines for minor grids
%set(gca, 'XMinorTick', 'on')  % Enable minor ticks on X-axis
%set(gca, 'YMinorTick', 'on')  % Enable minor ticks on Y-axis
%set(gca, 'FontSize', 12)  % Increase font size for better readability
%title('Bit Error Rate vs. Signal-to-Noise Ratio')
%xlabel('SNR (dB)')
%ylabel('Bit Error Rate (BER)')
%
%% Define a set of more distinguishable colors
%colors = {[0, 0.4470, 0.7410], [0.8500, 0.3250, 0.0980], [0.9290, 0.6940, 0.1250], ...
%          [0.4940, 0.1840, 0.5560], [0.4660, 0.6740, 0.1880], [0.3010, 0.7450, 0.9330], ...
%          [0.6350, 0.0780, 0.1840], [0.8500, 0.3250, 0.0980], [0.4940, 0.1840, 0.5560]};  % More distinct colors
%markers = {'o', 's', 's', 's', 's', 'v', 'v', 'v', 'v'};  % Same markers for SGDBF and MGDBF
%lineStyles = {'-', '-', ':', ':', ':', '-', ':', ':', ':'};  % Different line styles
%
%for i = 1:nb_algorithms
%    mean_errors = zeros(1, length(snr));
%    for j = 1:length(snr)
%        mean_errors(j) = sum(squeeze(errors(i,j,:,:)), "all") / (Kb * cpd_L * frames);
%    end
%
%    % Plot each algorithm with a unique style
%    semilogy(snr, mean_errors, 'LineWidth', 1.5, 'Color', colors{i}, ...
%             'Marker', markers{i}, ...
%             'LineStyle', lineStyles{i});
%    hold on
%end
%
%legendInfo = cell(1, nb_algorithms);
%for i = 1:nb_algorithms
%    if i == 1
%        legendInfo{i} = 'SGDBF';
%    end
%    if i == 2
%        legendInfo{i} = 'MGDBF, \epsilon = 0.001';
%    end
%    if i == 3
%        legendInfo{i} = 'SGDBF, \epsilon = 0.1';
%    end
%    if i == 4
%        legendInfo{i} = 'SGDBF, \epsilon = 1';
%    end
%    if i == 5
%        legendInfo{i} = 'SGDBF, \epsilon = 2';
%    end
%    if i == 6
%        legendInfo{i} = 'SGDBF, \epsilon = 5';
%    end
%    if i == 7
%        legendInfo{i} = 'SGDBF, \epsilon = 10';
%    end
%    if i == 8
%        legendInfo{i} = 'SGDBF, \epsilon = 0.5';
%    end
%    if i == 9
%        legendInfo{i} = 'SGDBF, \epsilon = 0.9';
%    end
%
%    % print the total number of errors for each algorithm
%    total_errors = sum(sum(squeeze(errors(i,:,:,:)),2),1);
%    disp(['Total number of errors for algorithm ' num2str(i) ' is ' num2str(total_errors)]);
%end
%
%legend(legendInfo, 'Location', 'southwest')  % Position the legend in the bottom left corner
%set(gca, 'YScale', 'log')  % Set the Y-axis to logarithmic scale
%set(gca, 'GridLineStyle', '-')  % Solid lines for main grids
%set(gca, 'MinorGridLineStyle', ':')  % Dotted lines for minor grids
%title('Bit Error Rate vs. Signal-to-Noise Ratio')
%xlabel('SNR (dB)')
%ylabel('Bit Error Rate (BER)')
%grid on 
%hold off




% ---------------------------------------------------------------------------------
% BER plot
% ---------------------------------------------------------------------------------
%figure
%grid on  % Ensure grid is on
%set(gca, 'YScale', 'log')  % Set the Y-axis to logarithmic scale
%set(gca, 'YMinorGrid', 'on')  % Enable minor grid lines on Y-axis
%set(gca, 'XMinorGrid', 'on')  % Enable minor grid lines on X-axis
%set(gca, 'GridLineStyle', '-')  % Solid lines for grids
%set(gca, 'MinorGridLineStyle', ':')  % Dotted lines for minor grids
%set(gca, 'XMinorTick', 'on')  % Enable minor ticks on X-axis
%set(gca, 'YMinorTick', 'on')  % Enable minor ticks on Y-axis
%set(gca, 'FontSize', 12)  % Increase font size for better readability
%title('Bit Error Rate vs. Signal-to-Noise Ratio')
%xlabel('SNR (dB)')
%ylabel('Bit Error Rate (BER)')
%
%% Define a set of more distinguishable colors
%colors = {[0, 0.4470, 0.7410], [0.8500, 0.3250, 0.0980], [0.9290, 0.6940, 0.1250], ...
%          [0.4940, 0.1840, 0.5560], [0.4660, 0.6740, 0.1880], [0.3010, 0.7450, 0.9330], ...
%          [0.6350, 0.0780, 0.1840], [0.75, 0, 0.75], [0, 0.5, 0], [0.75, 0.75, 0]} ;  % Unique colors for each curve
%markers = {'o', 'o', 's', 's', '^', 's', 'v', 'v', 'v'};  % Same markers for SGDBF and MGDBF
%lineStyles = {'-', ':', '--', '-.', '-', '-.', ':', ':', ':'};  % Different line styles
%
%for i = 1:nb_algorithms
%    mean_errors = zeros(1, length(snr));
%    for j = 1:length(snr)
%        mean_errors(j) = sum(squeeze(errors(i,j,:,:)), "all") / (Kb * cpd_L * frames);
%    end
%
%    % Plot each algorithm with a unique style
%    semilogy(snr, mean_errors, 'LineWidth', 1.5, 'Color', colors{i}, ...
%             'Marker', markers{i}, ...
%             'LineStyle', lineStyles{i});
%    hold on
%end
%
%legendInfo = cell(1, nb_algorithms);
%for i = 1:nb_algorithms
%   %if i == 1
%   %    legendInfo{i} = 'MGDBF, \alpha = 1, \beta = 1';
%   %end
%   %if i == 2
%   %    legendInfo{i} = 'MGDBF, \alpha = 2.5, \beta = 1.1';
%   %end
%   %if i == 3
%   %    legendInfo{i} = 'MGDBF, \alpha = 3, \beta = 1.5';
%   %end
%   %if i == 4
%   %    legendInfo{i} = 'MGDBF, \alpha = 1.9, \beta = 1';
%   %end
%   %if i == 5 
%   %    legendInfo{i} = 'MGDBF, \alpha = 1.5, \beta = 1.5';
%   %end
%    if i == 1
%        legendInfo{i} = 'No decoding';
%    elseif i == 2
%        legendInfo{i} = 'LC-HWBF, \theta_{att} = 0.9';
%    elseif i == 3
%        legendInfo{i} = 'SGDBF';
%    elseif i == 4 
%        legendInfo{i} = 'ISGDBF, \epsilon = 0.001';
%    end
%    if i == 5 
%        legendInfo{i} = 'MGDBF, \theta = -0.5,';
%    end
%    if i == 6 
%        legendInfo{i} = 'IMGDBF, \theta = -0.5, \alpha = 1, \beta = 1';
%    end
%    if i == 7 
%        legendInfo{i} = 'IMGDBF, \theta = -0.5, \alpha = 2.5, \beta = 1.1';
%    end
%    if i == 8
%        legendInfo{i} = 'IMGDBF, \theta = -0.5, \alpha = 1.9, \beta = 0.9';
%    end
%    if i == 9
%        legendInfo{i} = 'IMGDBF, \theta = -0.5, \alpha = 1.9, \beta = 0.9';
%    end
%
%    % print the total number of errors for each algorithm
%    total_errors = sum(sum(squeeze(errors(i,:,:,:)),2),1);
%    disp(['Total number of errors for algorithm ' num2str(i) ' is ' num2str(total_errors)]);
%end
%
%legend(legendInfo, 'Location', 'southwest')  % Position the legend in the bottom left corner
%set(gca, 'YScale', 'log')  % Set the Y-axis to logarithmic scale
%set(gca, 'GridLineStyle', '-')  % Solid lines for main grids
%set(gca, 'MinorGridLineStyle', ':')  % Dotted lines for minor grids
%title('Bit Error Rate vs. Signal-to-Noise Ratio')
%xlabel('SNR (dB)')
%ylabel('Bit Error Rate (BER)')
%grid on 
%hold off






%%% Plotting the number of flips
%figure;
%hold on;
%
%% Define the line styles, colors, and markers
%lineStyles = {'-', ':', '--', '-.', '-', '-.', ':', ':', ':'};  % Different line styles
%colors = {[0, 0.4470, 0.7410], [0.8500, 0.3250, 0.0980], [0.9290, 0.6940, 0.1250], ...
%          [0.4940, 0.1840, 0.5560], [0.4660, 0.6740, 0.1880], [0.3010, 0.7450, 0.9330], ...
%          [0.6350, 0.0780, 0.1840], [0.75, 0, 0.75], [0, 0.5, 0], [0.75, 0.75, 0]};  % Unique colors for each curve
%markers = {'o', 's', 'd', '^', '>', 'v', '<', 'p', 'h', 'x'};  % Different markers
%
%for i = 1:nb_algorithms
%    mean_num_flips = mean(squeeze(num_flips_all(i, :, :)), 2);
%    % Plot each algorithm with a unique style
%    plot(snr, mean_num_flips, 'LineWidth', 1.5, 'Color', colors{i}, 'LineStyle', lineStyles{i}); %num_flips_all(i, :)
%    %plot(0:1:99, mean_num_flips, 'LineWidth', 1.5, 'Color', colors{i}, 'LineStyle', lineStyles{i}); %num_flips_all(i, :)
%end
%
%% Create legend information
%legendInfo = cell(1, nb_algorithms);
%for i = 1:nb_algorithms
%    if i == 1
%        legendInfo{i} = 'Single GDBF'; % Single 'LC-HWBF, \theta = 0.9';
%    elseif i == 2
%        legendInfo{i} = 'Multi GDBF'; %'SGDBF, \epsilon = 0.001';
%    elseif i == 3
%        legendInfo{i} = 'Multi GDBF';
%    elseif i == 4
%        legendInfo{i} = 'MGDBF, \theta = -0.5';
%    elseif i == 5 
%        legendInfo{i} = 'IMGDBF, \theta = 0.5';
%    elseif i == 6 
%        legendInfo{i} = 'MGDBF, \theta = 0.01';
%    elseif i == 7 
%        legendInfo{i} = 'MGDBF, \theta = 0.1';
%    elseif i == 8 
%        legendInfo{i} = 'MGDBF, \theta = 0.5';
%    elseif i == 9
%        legendInfo{i} = 'MGDBF, \theta = 0.9';
%    end
%end
%
%% Plot settings
%xlabel('SNR');
%ylabel('Number of Iterations');
%%xlabel('SNR (dB)')
%%ylabel('Average number of bits flipped')
%title('Number of Iterations vs. SNR');
%legend(legendInfo, 'Location', 'southwest');  % Position the legend in the bottom left corner
%grid on;
%hold off;
%
%