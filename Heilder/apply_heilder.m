% Program to filter and apply fit_heilder function

if (exist("data","var") == 0)
    data = (load("c_2013_11_16_04.dat"))';
end


figure(1)



subplot(4,1,1)
hold off
fs = 187;
ts = 1/187;

total_time = ts * length(data);
tx = linspace(0,length(data) * ts, length(data));
tx_lim = [210, 220];
selected = select_chunk(tx, tx_lim);

plot(tx(selected), data(selected), "b:")
title("Data original");
ylabel("AU");
xlabel("Time (s)");

subplot(4,1,2)

L = length(data);
f = fs*(0:(L/2))/L;
Y = fft(data);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

[P1_v, P1_f] = SR_peak_welch(data,fs, SR_config.window, SR_config.overlap);
plot(P1_f, P1_v);
ylim([40 120]);
title("F transform original");
ylabel("dB");
xlabel("Time (s)");


data_preprocess = preprocess_band_pass(data, tx, 1, 16);

subplot(4,1,3)
plot(tx(selected), data_preprocess(selected), "r");

title("Signal prefiltered [1 - 16]Hz");
ylabel("dB");
xlabel("Time (s)");

subplot(4,1,4)
[P2_v, P2_f] = SR_peak_welch(data_preprocess,fs, SR_config.window, SR_config.overlap);
plot(P2_f, P2_v);
ylim([40 120]);
title("F transform prefiltered [1 - 16]Hz");
ylabel("dB");
xlabel("Time (s)");

sig = data_preprocess;
 hold off
%First method:

figure(2)
subplot(2,1,1)
Fs = 187;
numSamples = length(tx);
DownsampleFactor = 15;
frameSize = 10 * DownsampleFactor;

sigsq = 2 * sig .* sig;


rp = 1;           % Passband ripple in dB 
rs = 60;          % Stopband ripple in dB
fs = 187;        % Sampling frequency
f = [4 6];    % Cutoff frequencies
a = [1 0];        % Desired amplitudes

dev = [(10^(rp/20)-1)/(10^(rp/20)+1) 10^(-rs/20)]; 

[n,fo,ao,w] = firpmord(f,a,dev,fs);

b = firpm(n,fo,ao,w);


sigsqdownlp = filter(b, 1, sigsq);

shifted_data = delayseq(sigsqdownlp'  , -n/2);

output_1 = abs(sqrt(shifted_data'));


plot(tx(selected), data(selected), "b:")
hold on
title("Data Method 1 ");
ylabel("AU");
xlabel("Time (s)");
plot(tx(selected), output_1(selected), "r");

L = length(output_1);
f = fs*(0:(L/2))/L;

Y = fft(output_1);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
subplot(2,1,2)
hold off
plot(f, P1, "r");

%% Second Method
figure(3)

N = 60; % FILTER ORDER
b = firpm(N,[0.01 .95],[1 1],'hilbert');

sighil = filter(b,1,sig);

sighildelay = delayseq(sighil'  , -N/2);

sigprocess = abs(complex (0,sighildelay'));

subplot(3,2,1);
plot(tx(selected), sigprocess(selected));

title("Hilbert Signal");
ylabel("AU");
xlabel("Time (s)");

subplot(3,2,2);
[Hilbert_v_1, Hilbert_f_1] = SR_peak_welch(sigprocess,fs, SR_config.window, SR_config.overlap);
plot(Hilbert_f_1, Hilbert_v_1);
ylim([40 120]);
title("FFT");
ylabel("dB");
xlabel("Frequency (Hz)");

DownsampleFactor = 2;

signal_downsample = downsample(sigprocess,DownsampleFactor);
tx_downsampled = linspace(0, total_time, length(signal_downsample));
select_downsampled = select_chunk(tx_downsampled, tx_lim);

subplot(3,2,3);
plot(tx_downsampled(select_downsampled), signal_downsample(select_downsampled));

title("Hilbert Signal Downsampled");
ylabel("AU");
xlabel("Time (s)");
fs_downsample = linspace(1,fs/ DownsampleFactor, length(signal_downsample));

subplot(3,2,4);
[Hilbert_v_2, Hilbert_f_2] = SR_peak_welch(signal_downsample,fs/DownsampleFactor, SR_config.window, SR_config.overlap);
plot(Hilbert_f_2, Hilbert_v_2);
title("FFT");
ylabel("dB");
xlabel("Frequency (Hz)");
ylim([40 120]);



rp = 1;           % Passband ripple in dB 
rs = 60;          % Stopband ripple in dB
%fs = 187/DownsampleFactor ;        % Sampling frequency
f = [4.5 6];    % Cutoff frequencies
a = [1 0];        % Desired amplitudes

dev = [(10^(rp/20)-1)/(10^(rp/20)+1) 10^(-rs/20)]; 

[n,fo,ao,w] = firpmord(f,a,dev,fs/DownsampleFactor);

b = firpm(n,fo,ao,w);

sigfilt = filter(b, 1, signal_downsample);

output_2 = delayseq(sigfilt'  , -n/2);
subplot(3,2,5);
plot(tx_downsampled(select_downsampled), output_2(select_downsampled),":", "LineWidth", 2.5);

title("Hilbert Signal filtered");
ylabel("AU");
xlabel("Time (s)");


subplot(3,2,6);
[Hilbert_v_3, Hilbert_f_3] = SR_peak_welch(output_2,fs/DownsampleFactor, SR_config.window, SR_config.overlap);
plot(Hilbert_f_3, Hilbert_v_3);
title("FFT");
ylabel("dB");
xlabel("Frequency (Hz)");
ylim([40 120]);
hold off
save_fig('large', "IMG/hilbert_process", "png");


for i=1:6
    subplot(3,2,i)
    title([]);
    ylabel([]);
    xlabel([]);
    set(gca,'xticklabel',[])
    set(gca,'yticklabel',[])
end
subplot(3,2,1)
ylabel("Hilbert Tansform");
subplot(3,2,3)
ylabel("Downsampling");
subplot(3,2,5)
ylabel("Envelope Extract");
sgtitle([]);
save_fig('normal', "IMG/hilbert_process_visio", "png");
%% Find peak
figure(4)


select_min_peaks = tx_lim;%[210 220];



select_down = select_chunk(tx_downsampled, select_min_peaks);
select = select_chunk(tx, select_min_peaks);
Fs = 1/mean(diff(tx_downsampled)); % Average sample rate

values_find_peak = output_2(select_down);

tx_find_peak = tx_downsampled(select_down);

%findpeaks(values_find_peak,Fs,'MinPeakProminence',0.25e6,'MinPeakDistance',0.3);
[peaks, locs, width] = findpeaks(values_find_peak,tx_find_peak,'MinPeakProminence',0.1e6,'MinPeakDistance',0.3);


hold on

for ii = 1:length(peaks)
    text(locs(ii),peaks(ii) + 0.14 * max(peaks),num2str(ii),'Color','k', 'FontSize', 15)
end
ptl2 = plot( tx(select)  , data(select), "Color", [0.3 0.2 0.8], 'DisplayName','Raw Data');
ptl2 = plot(tx_find_peak, values_find_peak, "Color", [0.3 0.8 0.6], "LineWidth", 4, 'DisplayName','Envelope Hilbert');
%plot(tx_find_peak, values_find_peak, "c:");
minIndices = islocalmin(values_find_peak,'MinProminence',0.06e6,'FlatSelection','first','MinSeparation',0.1);


% Plot local minima
h1 = plot(tx_find_peak((minIndices)) ,values_find_peak(minIndices),'v','Color',[50 40 230]/255,...
    'MarkerFaceColor',[30 190 230]/255,'DisplayName','Possible separation points');
%title(['Number of extrema: ' num2str(nnz(minIndices))])

% Plot peaks maxima
h1 = plot(locs,peaks,'^','Color',[170 40 32]/255,...
    'MarkerFaceColor',[237 177 32]/255,'DisplayName','Peaks','MarkerSize',6, "LineWidth", 2);
%title(['Number of extrema: ' num2str(nnz(minIndices))])

tx_min = tx_find_peak((minIndices));
lg = legend
title("Transient Peak Detector");
ylabel("AU");
xlabel("Time (s)");
grid on

    title([]);
    ylabel([]);
    xlabel([]);
    set(gca,'xticklabel',[])
    set(gca,'yticklabel',[])
set(lg,'visible','off');
set(ptl2, "Color", [70, 130, 180]/255);
save_fig('normal', "IMG/segmentation_visio", "png");



hold off
%% Heilder 
tx_selected = tx_find_peak;
save_data = zeros(8,length(peaks));
minimum_value_left = (1/7.8);
minimum_value_right = (1/7.8);

for i = 1:length(peaks)
    
    
    distant_left = locs(i) - tx_min;
    % Discard negative: Theses minimums are in the other side
    
    distant_left(distant_left < minimum_value_left) = Inf; 
    
    if ~all(distant_left == Inf)
        min_distant_left = min(distant_left);
        idx = find(distant_left == min_distant_left);
        min_left_t = tx_min(idx);
    else 
        min_left_t = tx_selected(1);
    end
    
    distant_right = tx_min - locs(i);
    % Discard negative: Theses minimums are in the other side
    
    distant_right(distant_right < minimum_value_right) = Inf;
    
    if ~all(distant_right == Inf)
        min_distant_right = min(distant_right);
        idx = find(distant_right == min_distant_right);
        min_right_t = tx_min(idx);
    else 
        min_right_t = tx_selected(end);
    end
    
    envelope_select = select_chunk(tx_find_peak,[min_left_t min_right_t]);
    tx_to_fit = tx(envelope_select) ;
    
    values_to_fit = values_find_peak(envelope_select);
    [fit, gof, fitparam] = fit_heilder(tx_to_fit, values_to_fit');  
    
    save_data(1,i) = gof.rmse;
    save_data(2,i) = fitparam.t1;
    save_data(3,i) = fitparam.t2;
    save_data(4,i) = fitparam.nh;
    save_data(5,i) = fitparam.Io;
    save_data(6,i) = peaks(i);
    save_data(7,i) = fitparam.T1;
    save_data(8,i) = fitparam.T2;
    
    
end
selected = select_ray(save_data);

function selected = select_ray(save_data)
minimum_value_left = (1/7.8);
minimum_value_right = (1/7.8) * 2;
    selected = ones(1,size(save_data,2));
    % FALL TIME HAS TO BE LONGER THAN RISE TIME
    for i = 1:size(save_data,2)
        selected(i) = selected(i) & (save_data(3,i) > save_data(2,i));
        selected(i) = selected(i) & (save_data(2,i) <  minimum_value_left);
        selected(i) = selected(i) & (save_data(2,i) <  minimum_value_right);
    end
end


function selected  = select_chunk(x, x_lim)

    selected = (x > x_lim(1) & x < x_lim(2));
    
end