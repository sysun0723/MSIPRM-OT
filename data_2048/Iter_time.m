figure('Position', [100, 100, 500, 500]); 
data1 = [ 
    26 23 25; % 32×32 
    36 34 34; % 64×64 
    37 51  0; % 128×128 
    57  0  0; % 256×256 
    83  0  0; % 512×512 
    39  0  0; % 1024×1024 
    59  0  0  % 2048×2048 
]; 
data2=[
0.152 0.1097 0.12264;
0.5305 0.53671 0.53772;
2.4595 3.4945 0;
15.4515 0 0;
89.6669 0 0;
181.1988 0 0;
1186.1366 0 0;
];


sizes = {'32×32', '64×64', '128×128', '256×256', '512×512', '1024×1024', '2048×2048'};
colors = [0.2 0.6 0.8; 0.8 0.4 0.2; 0.4 0.8 0.4]; 
 
group_gap = 0.8; 
bar_width = 0.3; 
 

x_pos = 1:group_gap:(1 + (length(sizes)-1)*group_gap);
x_left = x_pos - bar_width/2;  
x_right = x_pos + bar_width/2; 
 

yyaxis left;
b1 = bar(x_left, data1, 'stacked', 'BarWidth', bar_width); 
for i = 1:3 
    b1(i).FaceColor = colors(i,:); 
    b1(i).EdgeColor = 'none'; 
end 
ylabel('Number of IPRM Iterations', 'FontSize', 12, 'FontWeight', 'bold'); 
ylim([0 max(sum(data1,2))*1.2]);

yyaxis right;
b2 = bar(x_right, data2, 'stacked', 'BarWidth', bar_width, 'FaceAlpha', 0.3); 
for i = 1:3 
    b2(i).FaceColor = colors(i,:); 
    b2(i).EdgeColor = 'none'; 
end 
set(gca, 'YScale', 'log', 'FontWeight', 'bold');
ylabel('IPRM Time (s)', 'FontSize', 12, 'FontWeight', 'bold'); 
ylim([10^-1, 10^4]);
%set(gca, 'YTick', 10.^(-1:4));
set(gca, 'YTickLabel', {'10^{-1}', '10^{0}', '10^{1}', '10^{2}', '10^{3}', '10^{4}'}, 'FontWeight', 'bold');

set(gca, 'XTick', x_pos, 'XTickLabel', sizes, 'FontSize', 10, 'FontWeight', 'bold'); 
xlim([x_pos(1)-0.5, x_pos(end)+0.5]); 
xlabel('Resolution', 'FontSize', 12); 
 

grid on;
set(gca, 'GridColor', [0.3 0.3 0.3], 'GridAlpha', 0.5, 'YMinorGrid', 'off', 'FontWeight', 'bold');
set(gcf, 'Color', 'w');

legend_entries = [b1(1), b2(1),b1(2), b2(2), b1(3),  b2(3)];
 
legend_labels = {
    ' ', ' ', ' ', ' ', '   IPRM Iterations','   IPRM Time(s)' 
};

legend(legend_entries, legend_labels, ...
       'Location', 'northeast', ...
       'NumColumns', 3);