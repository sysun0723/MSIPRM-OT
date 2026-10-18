sizes = {'32×32', '64×64', '128×128', '256×256', '512×512', '1024×1024', '2048×2048'};
solvelinear_time = [0.0026 0.012 0.049 0.21 0.83 3.6 16.3];
per_mille = [16.016 3.968 0.993 0.244 0.061 0.016 0.004];  
 
figure('Position', [100, 100, 500, 500]);
 
yyaxis left 
p1 = plot(1:length(solvelinear_time), solvelinear_time, 'b-o', ...
    'LineWidth', 1.5, ...
    'MarkerSize', 5, ...
    'MarkerFaceColor', 'b');
 
ylabel('Linear System Solving Time (s)', 'FontSize', 12, 'FontWeight', 'bold');

ylim([-1 17]);
 

yyaxis right 
p2 = plot(1:length(per_mille), per_mille, 'r-^', ...
    'LineWidth', 1.5, ...
    'MarkerSize', 5, ...
    'MarkerFaceColor', 'r');
 
ylabel('Filling Rate (‰)', 'FontSize', 12, 'Color', 'r', 'FontWeight', 'bold');
ax = gca;
ax.YAxis(2).Color = 'r';
ylim([-1 17]); 
 
xtick_pos = 1:length(sizes);
set(gca, 'XTick', xtick_pos, 'XTickLabel', sizes, 'FontSize', 10, 'FontWeight', 'bold');
xlim([0.5 length(sizes)+0.5]);
xlabel('Resolution', 'FontSize', 12, 'FontWeight', 'bold');
%title('Performance Metrics by Resolution', 'FontSize', 14);
grid on;
 

% for i = 1:length(solvelinear_time)
%     text(xtick_pos(i), solvelinear_time(i), ...
%         sprintf('%.3f s', solvelinear_time(i)), ...
%         'HorizontalAlignment', 'center', ...
%         'VerticalAlignment', 'bottom', ...
%         'FontSize', 8);
% end 
%  

% for i = 1:length(per_mille)
%     text(xtick_pos(i), per_mille(i), ...
%         sprintf('%.3f‰', per_mille(i)), ...
%         'HorizontalAlignment', 'left', ...
%         'VerticalAlignment', 'middle', ...
%         'Color', 'r', ...
%         'FontSize', 8);
% end 
ax_pos = get(gca, 'Position');
 

legend_pos = [0.4, 0.8, 0.3, 0.1];  
legend([p1, p2], ...
    {'Linear System Solving Time (s)', 'Filling Rate (‰)'}, ...
    'Position', legend_pos, ...
    'FontSize', 10, ...
    'FontWeight', 'bold'); 
set(gca, 'GridColor', [0.3 0.3 0.3], 'GridAlpha', 0.5);  
set(gcf, 'Color', 'w');