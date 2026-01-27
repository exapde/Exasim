digitize_figs_new
drawstr = '-';
load('results.mat'); iplt = i_av;
viz_av_cont


viz_av_cont% load('results_k_cv_bou2.mat');iplt = ®i_av;

figure(105); hold on; scatter(exp_Q_I(:,1), exp_Q_I(:,2), 150, 'square','k','filled')
legend(["Noncat","Supercat","$\gamma=1.0$"],"Interpreter","latex")
ylabel("$q_w$","Interpreter","latex")
xlabel("$\theta$", "Interpreter","latex")
set(gca,'FontSize',16);
grid on;
%%
figure(100); hold on; scatter(exp_P_I(:,1), exp_P_I(:,2), 150, 'square','k','filled')
legend(["Noncat","Supercat","$\gamma=1.0$"],"Interpreter","latex")
ylabel("$p_w$","Interpreter","latex")
xlabel("$\theta$", "Interpreter","latex")
set(gca,'FontSize',16);
grid on;