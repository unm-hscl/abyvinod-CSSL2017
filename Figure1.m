clear
close all

%%%%% Implementation of FTBU %%%%%%%%%%%
usePatternSearch=1;
useFmincon=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Generate Figure 1 (Choose the following options)');
code_flow_flag = input('1-Compute via the scripts for FTBU and DP | 2-Load the matfile\n');
if code_flow_flag == 2
    load('Figure1_curseOfDim_full_10_5_20points_n1to40_reqdOnly');
    run_DP_script = 0;
elseif code_flow_flag == 1
    table_1_generation = 0;
    run_DP_script = 1;
    scriptForChainDI
else
    error('Invalid option');
end

timeDP=cat(1, elapsed_time_DP{:});
maxvalue=max(timeDP)*1.2;
% logDP=log10(timeDP);
stateDP=1:length(timeDP);

% timeFT=cat(1, elapsed_time_FT{:})/no_of_testing_initial_states;
% logFT=log(timeFT);
% stateFT=1:length(timeFT);

timeFT_FM=cat(1, elapsed_time_FT_FM{:})/no_of_testing_initial_states;
% logFT_FM=log10(timeFT_FM);
stateFT_FM=1:length(timeFT_FM);

timeFT_PS=cat(1, elapsed_time_FT_PS{:})/no_of_testing_initial_states;
% logFT_PS=log10(timeFT_PS);
stateFT_PS=1:length(timeFT_PS);

figure(1)
clf
hold on
% plot(stateDP,logDP,'bx--');
% plot(stateFT_FM,logFT_FM,'ro--');
% plot(stateFT_PS,logFT_PS,'kd--');
% ylabel('Runtime ($10^x$ seconds)','interpreter','latex')
% axis([0 41 -1 5])
plot(stateDP,timeDP,'bx--');
plot(stateFT_FM,timeFT_FM,'ro--');
plot(stateFT_PS,timeFT_PS,'kd--');
ylabel('Runtime (seconds)')
set(gca,'YScale','log');
set(gca,'YScale','log');
axis([0 41 0.1 maxvalue])
% l=legend('$\log(\mathrm{DPBDA\ runtime})$','$\log(\mathrm{FTBU\ runtime})$ (fmincon)','$\log(\mathrm{FTBU\ runtime})$ (patternsearch)');%,'$\log(n)+2$');
l=legend('DPBDA','FTBU (fmincon)','FTBU (patternsearch)');%,'$\log(n)+2$');
set(l,'interpreter','latex');
box on
grid on
set(gca,'FontSize',20)
xlabel('State dimension (n)')
set(gca,'YTick',10.^[-1:log10(maxvalue)])
% 
% figure(2)
% clf
% subplot(2,1,1);
% hold on
% plot(1:4,log(timeDP),'bx--');
% plot(1:20,log(timeFT),'ro--');
% box on
% grid on
% subplot(2,1,2);
% hold on
% plot(log(1:4),log(timeDP),'bx--');
% plot(log(1:20),log(timeFT),'ro--');
% box on
% grid on
