% Program to filter and apply fit_heilder function

    data = (load("c_2013_11_16_04.dat"))';

close all;
figure(1)



subplot(3,1,1)
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

subplot(3,1,2)

L = length(data);
f = fs*(0:(L/2))/L;
Y = fft(data);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

plot(f, P1);


data_preprocess = preprocess_2hz(data, tx, 1, 20);

subplot(3,1,3)
plot(tx(selected), data_preprocess(selected), "r");


sig = data_preprocess;

%First method:

figure(2)
subplot(2,1,1)
hold off
Fs = 187;
numSamples = length(tx);
DownsampleFactor = 15;
frameSize = 10 * DownsampleFactor;

sigsq = 2 * sig .* sig;


rp = 1;           % Passband ripple in dB 
rs = 60;          % Stopband ripple in dB
fs = 187;        % Sampling frequency
f = [5 6];    % Cutoff frequencies
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
subplot(2,1,1)
hold off
N = 60; % FILTER ORDER
b = firpm(N,[0.01 .95],[1 1],'hilbert');

sighil = filter(b,1,sig);

sighildelay = delayseq(sighil'  , -N/2);

sigprocess = abs(complex (0,sighildelay'));

DownsampleFactor = 2;

signal_downsample = downsample(sigprocess,DownsampleFactor);

rp = 1;           % Passband ripple in dB 
rs = 60;          % Stopband ripple in dB
fs = 187/DownsampleFactor ;        % Sampling frequency
f = [5 7];    % Cutoff frequencies
a = [1 0];        % Desired amplitudes

dev = [(10^(rp/20)-1)/(10^(rp/20)+1) 10^(-rs/20)]; 

[n,fo,ao,w] = firpmord(f,a,dev,fs);

b = firpm(n,fo,ao,w);

sigfilt = filter(b, 1, signal_downsample);

output_2 = delayseq(sigfilt'  , -n/2);

plot(tx(selected), data(selected), "b:")
hold on
title("Data Method 2");
ylabel("AU");
xlabel("Time (s)");


tx_downsampled = linspace(0, total_time, length(output_2));

select_downsampled = select_chunk(tx_downsampled, tx_lim);



plot(tx_downsampled(select_downsampled), output_2(select_downsampled), "r");
subplot(2,1,2)

L = length(output_2);
f = fs*(0:(L/2))/L;

Y = fft(output_2);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

plot(f, P1, "r");

%% Find peak
figure(4)
hold off

select_min_peaks = [210 220];



select_down = select_chunk(tx_downsampled, select_min_peaks);
select = select_chunk(tx, select_min_peaks);
Fs = 1/mean(diff(tx_downsampled)); % Average sample rate

values_find_peak = output_2(select_down);

tx_find_peak = tx_downsampled(select_down);

%findpeaks(values_find_peak,Fs,'MinPeakProminence',0.25e6,'MinPeakDistance',0.3);
[peaks, locs, width] = findpeaks(values_find_peak,tx_find_peak,'MinPeakProminence',0.1e6,'MinPeakDistance',0.4);


% Plot peaks maxima
h1 = plot(locs,peaks,'v','Color',[0 177 32]/255,...
    'MarkerFaceColor',[237 177 32]/255,'DisplayName','Local minima');
%title(['Number of extrema: ' num2str(nnz(minIndices))])

hold on

for ii = 1:length(peaks)
    text(locs(ii),peaks(ii) + 0.05 * max(peaks),num2str(ii),'Color','r', 'FontSize', 16)
end
plot( tx(select)  , data(select), "m");
plot(tx_find_peak, values_find_peak, "c");
%plot(tx_find_peak, values_find_peak, "c:");
minIndices = islocalmin(values_find_peak,'MinProminence',0.1e6,'FlatSelection','first','MinSeparation',0.1);


% Plot local minima
h1 = plot(tx_find_peak((minIndices)) ,values_find_peak(minIndices),'v','Color',[237 177 32]/255,...
    'MarkerFaceColor',[237 177 32]/255,'DisplayName','Local minima');
%title(['Number of extrema: ' num2str(nnz(minIndices))])

tx_min = tx_find_peak((minIndices));
legend

%% Heilder 
tx_selected = tx_find_peak;
save_data = zeros(5,length(peaks));
minimum_value_left = (1/7.8);
minimum_value_right = (1/7.8) * 2;
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
    [fit, gof, fitparam] = fit_heilder(tx_to_fit, values_to_fit);  
    title("PLOT NUMBER: " + i);
    save_data(1,i) = gof.rmse;
    save_data(2,i) = fitparam.t1;
    save_data(3,i) = fitparam.t2;
    save_data(4,i) = fitparam.nh;
    save_data(5,i) = fitparam.Io;
    save_data(6,i) = peaks(i);
    
    
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