function [AM_SR_ray_vector] = find_ray(tx, data)
    %FIND_RAY Summary of this function goes here
    %   Detailed explanation goes here

    %Sampling Frequency 
    Fs = 187;
    % Period
    Ts = 1/187;
    
    % Minimal separation rise part
    minimum_value_left = 0;%(1/7.8);
    % Minimal separation fall part
    minimum_value_right = 0;%(1/7.8) * 2;

    total_time = Ts * length(data);

    band_pass_lower_limit = 1; % 1Hz
    band_pass_upper_limit = 20; % 1Hz
    sig = preprocess_band_pass(data, tx, band_pass_lower_limit, band_pass_upper_limit);





    %% Hilbert Method
    N = 60; % FILTER ORDER

    % Calculate b coefficients for the hilbert filter
    b = firpm(N,[0.01 .95],[1 1],'hilbert');

    % Apply hilbert filter
    sig_hil = filter(b,1,sig);

    % Correct delay produced by the hilber filter
    sig_hil_delay = delayseq(sig_hil'  , -N/2);

    %  Magnitude of the complex num with the hilbert signal in the imaginary part 
    sig_process = abs(complex (0,sig_hil_delay'));

    % Downsample in order to reduce the ...
    DownsampleFactor = 2;

    signal_downsample = downsample(sig_process,DownsampleFactor);

    % to compare
    raw_downsample = downsample(data,DownsampleFactor);

     % LOW PASS FILTER
    rp = 1;           % Passband ripple in dB 
    rs = 60;          % Stopband ripple in dB
    fs = 187/DownsampleFactor ;        % Sampling frequency
    f = [4.5 6];    % Cutoff frequencies
    a = [1 0];        % Desired amplitudes

    %  ripple of the low pass filter
    dev = [(10^(rp/20)-1)/(10^(rp/20)+1) 10^(-rs/20)]; 

    [n,fo,ao,w] = firpmord(f,a,dev,fs);

    b = firpm(n,fo,ao,w);

    sig_filt = filter(b, 1, signal_downsample);

    output = delayseq(sig_filt'  , -n/2);

    tx_downsampled = linspace(0, total_time, length(output));


    %% Find peak


    values_find_peak = output;

    tx_find_peak = tx_downsampled;

    %findpeaks(values_find_peak,Fs,'MinPeakProminence',0.25e6,'MinPeakDistance',0.3);

    % FIND MAX POINTS
    [peaks, locs, width] = findpeaks(values_find_peak,tx_find_peak,'MinPeakProminence',0.1e6,'MinPeakDistance',0.3);

    % FIND MIN POINTS
    minIndices = islocalmin(values_find_peak,'MinProminence',0.06e6,'FlatSelection','first','MinSeparation',0.1);
    tx_mins = tx_find_peak((minIndices));


    %% Heilder 
    tx_selected = tx_find_peak;
    save_data = zeros(6,length(peaks));

 for i = 1:length(peaks)


        distant_left = locs(i) - tx_mins;
        % Discard negative: Theses minimums are in the other side

        distant_left(distant_left < minimum_value_left) = Inf; 

        if ~all(distant_left == Inf)
            min_distant_left = min(distant_left);
            idx = find(distant_left == min_distant_left);
            min_left_t = tx_mins(idx);
        else 
            min_left_t = tx_selected(1);
        end

        distant_right = tx_mins - locs(i);
        % Discard negative: Theses minimums are in the other side

        distant_right(distant_right < minimum_value_right) = Inf;

        if ~all(distant_right == Inf)
            min_distant_right = min(distant_right);
            idx = find(distant_right == min_distant_right);
            min_right_t = tx_mins(idx);
        else 
            min_right_t = tx_selected(end);
        end
        
        minimuns{i} = [min_left_t min_right_t];

        
        % Downsample in order to reduce the ...
 end
 index_minimuns_selected = {};
 current_index_miniums = minimuns{i};
  for index_minimuns = 2:(length(minimuns))
      if(minimuns{index_minimuns}(1) ~= current_index_miniums(1) && ...
        minimuns{index_minimuns}(2) ~= current_index_miniums(2) )  
        index_minimuns_selected{end+1} = minimuns{index_minimuns};
        current_index_miniums =  minimuns{index_minimuns};      
      end
  end
 for index_minimus = 1:length(index_minimuns_selected)
        envelope_select = select_chunk(tx_find_peak,index_minimuns_selected{index_minimus});
        tx_to_fit{index_minimus} = tx_find_peak(envelope_select) ;
        values_to_fit{index_minimus} = values_find_peak(envelope_select);
        raw_signal{index_minimus} = raw_downsample(envelope_select);
 end 
parfor index = 1:length(tx_to_fit)
    try 
        AM_SR_ray_vector_temp{index} = AM_SR_ray(  tx_to_fit{index}, values_to_fit{index}, raw_signal{index});
    catch ME
        AM_SR_ray_vector_temp{index} = [];
        disp(ME.message);
    end
end
for index = 1:length(AM_SR_ray_vector_temp)
    AM_SR_ray_vector(index) = AM_SR_ray_vector_temp{index};
end
 
    %% AUX FUNCTION
    function selected = select_ray(save_data)
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
end
