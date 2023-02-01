%y = load ("I:\Data\capturas_elf_cal\2016\04\ .mat");
y = (load("c_2013_11_16_04.dat"))';
start_plot = 111634;
end_plot = 114876;
y_select = y([start_plot:end_plot]);
t = [1:length(y_select)] * 1/187;
figure(1)
clf;


hold on;    
start_qburst = 112742;
end_qburst = 113064;
x_q_burst = [start_qburst - start_plot:end_qburst - start_plot];
t_q_burst = t(x_q_burst);
y_q_burst = y_select(x_q_burst);
plot(t,y_select,'Color', [0.3 0.8 0.6],'LineWidth', 0.5);
plot(t_q_burst,y_q_burst,'Color',[0.99 0.1 0.3],'LineWidth', 0.6)


grid on

xlim([t(1),t(end)]);
xlabel("Seconds (s)");
ylabel("Magnetic Field Amplitude (pT)");
title("ELF register");
legend(["Background ELF Noise","Q-Burst"]);
save_fig('ultrawide', "Intro_test_1", "png");
