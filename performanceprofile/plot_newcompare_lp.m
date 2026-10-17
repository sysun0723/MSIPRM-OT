clc; clear;
rng(3)
scale= 256;

fileID1 = fopen(['C:\Users\Napoleon\Desktop\IPRM4OT\plotport\solvingtime\data-', num2str(scale), '\iprm.txt'], 'r');
fileID2 = fopen(['C:\Users\Napoleon\Desktop\IPRM4OT\plotport\solvingtime\data-', num2str(scale), '\HOT.txt'], 'r');
fileID3 = fopen(['C:\Users\Napoleon\Desktop\IPRM4OT\plotport\solvingtime\data-', num2str(scale), '\CPLEXNet.txt'], 'r');
fileID4 = fopen(['C:\Users\Napoleon\Desktop\IPRM4OT\plotport\solvingtime\data-', num2str(scale), '\SparseSinkhorn.txt'], 'r');
%fileID5 = fopen(['C:\Users\Napoleon\Desktop\IPRM4OT\plotport\solvingtime\data-', num2str(scale), '\lemon.txt'], 'r');
A1_t = fscanf(fileID1,'%f');
A2_t = fscanf(fileID2,'%f');
A3_t = fscanf(fileID3,'%f');
A4_t = fscanf(fileID4,'%f');
%A5_t = fscanf(fileID5,'%f');


numsolver = 4; 

if size(A1_t,1)<360
pickproblem=(1:size(A1_t,1))';
else
random_remainders = randsample(0:44, 10, false); 
remainders = (mod(1:size(A1_t,1), 45))';
pickproblem = find(ismember(remainders, random_remainders));
end

num_of_probs = size(pickproblem,1);
c_t = zeros(num_of_probs,numsolver);


c_t(:,1) = A1_t(pickproblem,1);      %2 norm
c_t(:,2) = A2_t(pickproblem,1);   
c_t(:,3) = A3_t(pickproblem,1);   
c_t(:,4) = A4_t(pickproblem,1); 
%c_t(:,5) = A5_t(pickproblem,1);
%c_t(:,6) = A6_t(1:num_of_probs,1); 
tol= 10;
time_gen=log(min(c_t,40000)+tol);
tt=sum(time_gen,1)./num_of_probs;
tt=exp(tt)-tol;
re_tt=tt./(min(tt));

% disp('Arithmetic mean:')
% fprintf('iprm: %fs HOT: %fs  cplex-ns: %fs  Sp-Sinkhorn: %fs \nlemon: %fs', tt(1), tt(2), tt(3), tt(4), tt(5));
% disp('Geometric mean:')
% fprintf('iprm: %f  HOT: %f   cplex-ns: %f   Sp-Sinkhorn: %f  \nlemon: %f   ', re_tt(1), re_tt(2), re_tt(3), re_tt(4), re_tt(5));

disp('Arithmetic mean:')
fprintf('iprm: %fs HOT: %fs  cplex-ns: %fs  Sp-Sinkhorn: %fs\n', tt(1), tt(2), tt(3), tt(4));
disp('Geometric mean:')
fprintf('iprm: %f  HOT: %f   cplex-ns: %f   Sp-Sinkhorn: %f\n', re_tt(1), re_tt(2), re_tt(3), re_tt(4));

% disp('Arithmetic mean:')
% fprintf('iprm: %fs HOT: %fs  cplex-ns: %fs ', tt(1), tt(2), tt(3));
% disp('Geometric mean:')
% fprintf('iprm: %f  HOT: %f   cplex-ns: %f  ', re_tt(1), re_tt(2), re_tt(3));

c_min_t = zeros(num_of_probs,1);

for i = 1:num_of_probs
    c_min_t(i) = min(c_t(i,:));
end

r_t = zeros(num_of_probs,numsolver);

for j = 1:numsolver
    for i = 1:num_of_probs
        r_t(i,j) = c_t(i,j)/c_min_t(i);
       % r_it(i,j) = c_it(i,j)/c_min_it(i);
       % r_itall(i,j) = c_itall(i,j)/c_min_itall(i);
    end
end
N = 1000;
split = 2000;
t_t = linspace(1,N,split);
perf_func_t = zeros(split,numsolver);
for k = 1:split
    for j = 1:numsolver
        perf_func_t(k,j) = sum(r_t(:,j) <= t_t(k))/num_of_probs;
    end
end
figure(1);
%plot(t_t,perf_func_t(:,1),'m--',t_t,perf_func_t(:,2),'b.-',t_t,perf_func_t(:,3),'k:','LineWidth',1.5);
%plot(t_t,perf_func_t(:,1),'m--',t_t,perf_func_t(:,2),'b.-','LineWidth',1.5);
% semilogx(t_t,perf_func_t(:,1),'m--',t_t,perf_func_t(:,2),'c-',t_t,perf_func_t(:,3),'r-',t_t,perf_func_t(:,4),'k:',...
%     t_t,perf_func_t(:,5),'b.-',t_t,perf_func_t(:,6),'y.-','LineWidth',1.8);
%plot(t_t,perf_func_t(:,1),'m--',t_t,perf_func_t(:,2),'b.-',t_t,perf_func_t(:,3),'k:',t_t,perf_func_t(:,4),'c-','LineWidth',1.5);
%semilogx(t_t,perf_func_t(:,1),'m--',t_t,perf_func_t(:,2),'k:','LineWidth',1.5);
 % plot(t_t,perf_func_t(:,1),'m--',t_t,perf_func_t(:,2),'b.-');
axis([1 N 0 1 ]);
%lgd = legend('IP-PMM','IPRM','OSQP','Location','southeast');
%lgd = legend('IP-PMM','IPRM','Location','southeast');

% semilogx(t_t,perf_func_t(:,1), 'LineStyle' , '--' , 'LineWidth', 2);
% hold on;
% semilogx(t_t,perf_func_t(:,2), 'LineStyle',':', 'LineWidth', 2);
% hold on;
% semilogx(t_t,perf_func_t(:,3), 'LineStyle',':', 'LineWidth', 2);
% hold on;
% semilogx(t_t,perf_func_t(:,4), 'LineStyle',':', 'LineWidth', 2);
% hold on;
% semilogx(t_t,perf_func_t(:,5), 'LineStyle',':', 'LineWidth', 2);
% hold on;
% semilogx(t_t,perf_func_t(:,6), 'LineStyle',':', 'LineWidth', 2);
% hold on;
% semilogx(t_t,perf_func_t(:,7), 'LineStyle',':', 'LineWidth', 2);
% hold on;

% Step 5: 绘制图形
hold on;
grid on;
fig.Position = [1000,600,400,300];
stairs(t_t, perf_func_t(:,1), 'LineWidth', 2,'LineStyle','-.');
stairs(t_t, perf_func_t(:,2), 'LineWidth', 2,'LineStyle','-.');
stairs(t_t, perf_func_t(:,3), 'LineWidth', 2,'LineStyle','-.');
stairs(t_t, perf_func_t(:,4), 'LineWidth', 2,'LineStyle','-.');
%stairs(t_t, perf_func_t(:,5), 'LineWidth', 2,'LineStyle','-.');
% stairs(t_t, perf_func_t(:,6), 'LineWidth', 2,'LineStyle','-.');
% stairs(t_t, perf_func_t(:,7), 'LineWidth', 2,'LineStyle','-.');
set(gca,'XScale','log');


lgd = legend('MSIPRM','HOT','MCPLEX-NS','Sp-Sinkhorn','LEMON-NS','Location','southeast');
set(lgd,'FontName','Times New Roman','FontSize',10,'FontWeight','normal')
xlabel({'\tau'},'FontSize',20,'FontWeight','bold');
ylabel({'\rho_{s}'},'FontSize',20,'FontWeight','bold');
%saveas(gcf, 'port_time', 'png');
%title('time');

