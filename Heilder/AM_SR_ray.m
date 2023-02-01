classdef AM_SR_ray < handle
    %RAY_SR Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        start_time;
        duration;
        tx;
        data_envelope;
        data_raw;
        data_fit;
        t_rise;
        t_fall;
        nh;
        Io;
        peak_value;
        peak_location;
        rmse;
        rsquare;
        select;
        fit;
        T1;
        T2;
        RiseTime;
        FWHMT;
        instant_freq_carrier;
        instant_freq_mod;
        
        PSD_in_band;
        PSD_total;
        PSD_rest;
        
    end
    properties (Constant)
        donwsample_rate = 2;
    end
    methods
        function obj = AM_SR_ray(tx, signal, raw_signal)
            obj.start_time = tx(1);
            obj.duration = tx(end) - tx(1);
            obj.data_envelope = signal;
            obj.data_raw = raw_signal;
            obj.tx = tx;
            
            
            [obj.fit, gof, fitparam, obj.data_fit] = fit_heilder(tx, signal');  


            obj.t_rise = fitparam.t1;
            obj.t_fall = fitparam.t2;
            obj.nh = fitparam.nh;
            obj.Io = fitparam.Io;
            obj.T1 = fitparam.T1;
            obj.T2 = fitparam.T2;
            obj.RiseTime = fitparam.RiseTime;
            obj.FWHMT = fitparam.FWHMT;
            obj.peak_value = max(signal);
            obj.peak_location = tx((signal == obj.peak_value));
            obj.rmse = gof.rmse;
            obj.rsquare = gof.rsquare;
            
            obj.PSD_in_band = obj.get_SPC([6, 9]);
            obj.PSD_total = obj.get_SPC([12 100]);
            obj.PSD_rest = obj.get_SPC([0 100]) - obj.PSD_in_band;
            obj.instant_freq_carrier = obj.instant_freq(7.8);
            obj.instant_freq_mod = obj.instant_freq(1);
            obj.select = is_ray(obj);
        end
        function select = is_ray(obj)
            Ts = (1 / 7.8);
            
            select = 1;

            select = select & (obj.rmse < 2);
            select = select & (abs(obj.FWHMT / obj.RiseTime)  > 2);
            select = select & (abs(obj.T2 / obj.T1)  > 1);
            select = select & (obj.RiseTime <  Ts * 3);
            select = select & (obj.FWHMT <  Ts * 20);
            select = select & (obj.instant_freq_mod <  4);
            select = select & (obj.peak_value < 10e6); %Problem 
            select = select & (obj.instant_freq_carrier >  6 & obj.instant_freq_carrier < 10);
        end
       function select = is_qburst(obj)
            select = (obj.peak_value >= 2e6);
        end
        function plots =  plot(obj)
            
            plots(1) = plot(obj.tx, obj.data_envelope,"Color", [0.3 0.8 0.6],'LineWidth', 2,'DisplayName','Envelope Hilbert');
            
            hold on
            plots(2) = plot(obj.tx, obj.data_raw,"--",'Color', [0.3 0.2 0.8],'LineWidth', 1.5,'DisplayName','Raw Signal');
            
            plots(3) = plot(obj.tx, obj.data_fit,'Color', 'k','LineWidth', 1,'DisplayName','Heilder Fit');
            ylabel("AU");
            xlabel("seconds (s)");
            xlim([min(obj.tx), max(obj.tx)]);
            legend();
            hold off
            
        end
        function instant_freq = instant_freq(obj, centre_point)
            Ts = 1 / 7.8;
            lower_limit = (obj.peak_location - 1 * Ts);
            upper_limit = (obj.peak_location + 4 * Ts);
            MaxNumIMF = 5;
            %instant_freq_vector = zeros(1,MaxNumIMF);
            
            select_tx = (obj.tx > lower_limit) & (obj.tx < upper_limit);
            imf = emd(obj.data_raw(select_tx),'MaxNumIMF',MaxNumIMF,'Interpolation', 'spline' ,'SiftRelativeTolerance' , 0.01);
            temp_instant = mean(instfreq(imf,187/2,'Method','hilbert'));
            %instant_freq_vector(1:length(temp_instant)) = temp_instant;
            [value,index] = min(abs(temp_instant - centre_point));
            instant_freq = temp_instant(index);
        end
        
        function new_ray = new(obj)
            vector_time = linspace(obj.start_time, obj.start_time + obj.duration, length(obj.data_envelope));
            new_ray = AM_SR_ray(vector_time, obj.data_envelope, obj.data_raw);
        end
        function psd_power = get_SPC(obj, limits)
            if (obj.peak_value > 2e6)
                is_qburst = 1;
            end
            start_point = 1;
            if (obj.duration < obj.T1 + obj.T2)
                end_point = length(obj.data_raw);
            else
                final_time = obj.T1 + obj.T2 + obj.start_time;
                [~,index_min] = min(abs(obj.tx - final_time));
                end_point = index_min;
            end
            data_to_psd = obj.data_raw(start_point:end_point);
            [values, f] = periodogram(data_to_psd,rectwin(length(data_to_psd)),length(data_to_psd),187/2);
            psd_power = sum(values(f > limits(1) & f < limits(2)))/ length(data_to_psd);%fft(data_to_psd,SR_config.fs,[0 10])/length(data_to_psd);
        end
        function plts = plot_psd_in_band(obj)
            subplot(2,1,1)
            f = [6,9];
            fs = 187/2;
            y_in_band = bandpass(obj.data_raw, f, fs);
            plts(1) = plot(obj.tx, y_in_band);
            xlim([min(obj.tx), max(obj.tx)])
            title("Band Pass 7.8Hz")
            
            subplot(2,1,2)
            y_out_band = bandstop(obj.data_raw, f, fs);
            plts(2) = plot(obj.tx, y_out_band);
            xlim([min(obj.tx), max(obj.tx)])
            title("Band Stop 7.8Hz")
        end

    end
end

