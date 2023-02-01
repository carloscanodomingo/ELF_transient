close all    
for i = 1:4
        current_data = (load("c_2013_11_16_0" + num2str(i) + ".dat"))';
        figure(i)
        subplot(2,1,1)
        plot(current_data);
        subplot(2,1,2)
        [values, f] = SR_peak_welch(current_data, 187,1870, 0);
        plot(f, values)
        sgtitle("Capture " + num2str(i))
        save_fig('wide', "IMG/new_captures" + num2str(i), "pdf");
    end

if (exist("data","var") == 0)
    data = (load("c_2013_11_16_04.dat"))';
end

figure()
fs = 187;
ts = 1/187;
subplot(2,1,1)
total_time = ts * length(data);
tx = linspace(0,length(data) * ts, length(data));
plot(tx,data);
set(gca,'xticklabel',[])
set(gca,'yticklabel',[])

subplot(2,1,2)


[values, f] = SR_peak_welch(data, 187, 1870,0);

values = smoothdata(values,'sgolay',15);
values = smoothdata(values,'sgolay',15);

data_f_power = 10 .^ (values / 10);
plot_y = sqrt(data_f_power);


semilogy(f, plot_y,'LineWidth',2);
set(gca,'xticklabel',[])
set(gca,'yticklabel',[])
save_fig('wide', "IMG/comparison", "png");
figure();

semilogy(f, plot_y,'LineWidth',2, "Color", 'k');
save_fig('wide', "IMG/welch_complete", "png");
figure();



ylim([1e3 1e6]);
lower_band_f = f(f < 15);
lower_band_y = plot_y(f < 15);

semilogy(lower_band_f, lower_band_y,'LineWidth',2,'Color', 'k');

ylabel("Power/frequency (dB/Hz)");
xlabel("Frequency (Hz)");
title("FFT 1st SR")
save_fig('narrow', "IMG/welch_narrow", "png");
figure();
hold off



load('AM_SR_ray_v.mat')


rays = AM_SR_ray_vector.select_ray(AM_SR_ray_v);

Qburst = AM_SR_ray_vector.separate_qburst(rays);

select_time = 298;

[Time_min,index_time_min] = min(abs([AM_SR_ray_v.start_time] - select_time));

plots = AM_SR_ray_v(index_time_min).plot();

save_fig('normal', "IMG/ray", "png");
    title([]);
    ylabel([]);
    xlabel([]);
    set(gca,'xticklabel',[])
    set(gca,'yticklabel',[])
    lg1 = legend(gca);
    set(lg1, "Visible", 0)
    set(plots(1), "Color", [70, 130, 180]/255);
    set(plots(1), "LineStyle", ":");
    set(plots(1), "LineWidth", 3);
    set(plots(3), "LineWidth", 2.5);
    set(plots(3), "Color", [20, 20, 80]/255);
save_fig('normal', "IMG/ray_visio", "png");





figure()
non_rays = AM_SR_ray_vector.select_non_ray(AM_SR_ray_v);



figure()
subplot(4,2,1)
plts = non_rays(2).plot();
remove_ray_visio(plts);
set(plts(3), "Color", [180, 20, 50]/255);

subplot(4,2,2)
plts = Qburst(2).plot();
remove_ray_visio(plts);
set(plts(3), "Color", [20, 180, 50]/255);

subplot(4,2,3)
plts = non_rays(67).plot();
remove_ray_visio(plts);
set(plts(3), "Color", [180, 20, 50]/255);

subplot(4,2,4)
plts = Qburst(5).plot();
remove_ray_visio(plts);
set(plts(3), "Color", [20, 180, 50]/255);

subplot(4,2,5)
plts = non_rays(180).plot();
remove_ray_visio(plts);
set(plts(3), "Color", [180, 20, 50]/255);

subplot(4,2,6)
plts = Qburst(8).plot();
remove_ray_visio(plts);
set(plts(3), "Color", [20, 180, 50]/255);

subplot(4,2,7)
plts = non_rays(360).plot();
remove_ray_visio(plts);
set(plts(3), "Color", [180, 20, 50]/255);

subplot(4,2,8)
plts = Qburst(11).plot();
remove_ray_visio(plts);
set(plts(3), "Color", [20, 180, 50]/255);
save_fig("narrow", "IMG/reject", "png")



Qburst(11).instant_freq(7.8)
emd(Qburst(13).data_raw,'MaxNumIMF',4,'Interpolation', 'spline' ,'SiftRelativeTolerance' , 0.01);
save_fig('normal', "IMG/emd", "png")
axes = get(gcf, 'children');


for i =2:length(axes)
    set(axes(i), 'title', []);
    set(axes(i), 'xlabel', []);
    set(axes(i), 'ylabel', []);
    set(axes(i), 'yticklabel', []);
    set(axes(i), 'xticklabel', []);
end
save_fig('normal', "IMG/emd_visio", "png")


imf = emd(Qburst(13).data_raw,'MaxNumIMF',4,'Interpolation', 'spline' ,'SiftRelativeTolerance' , 0.01);
hht(imf, 187/2)
save_fig('normal', "IMG/hht", "png")
lg1 = legend(gca);
set(lg1, "Visible", 0);
title([]);
ylabel([]);
xlabel([]);
set(gca,'xticklabel',[])
set(gca,'yticklabel',[])
cb = colorbar(gca);
set(cb, "Visible", 0);
save_fig('normal', "IMG/hht_visio", "png")


%Comparative

Qburst_present = Qburst(13);

plt = Qburst(13).plot_psd_in_band();

for i=1:length(plt)
    set(plt(i), "LineWidth", 3);
    set(plt(i), "Color", [70, 130, 180]/255);
    subplot(2,1,i)
    xticklabels([])
    yticklabels([])
    title([])
end
save_fig("normal","IMG/psds", "png");

function remove_ray_visio(plt)
    title([]);
    ylabel([]);
    xlabel([]);
    set(gca,'xticklabel',[])
    set(gca,'yticklabel',[])
    lg1 = legend(gca);
    set(lg1, "Visible", 0)
    set(plt(1), "Color", [70, 130, 180]/255);
    set(plt(1), "LineStyle", ":");
    set(plt(1), "LineWidth", 3);
    set(plt(3), "LineWidth", 2.5);
    set(plt(3), "Color", [20, 20, 80]/255);
end