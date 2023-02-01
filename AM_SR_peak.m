classdef AM_SR_peak < handle
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        data_original;
        data_envelope;
        data_hilbert;
        start_time;
        duration_simple;
        duration_p2;
        duration_lorentz;
        max_value;
        max_time;
        increment_value;
        x_data;
        y_data;
        ftpoly2;
        ft_lorentz_left;
        ft_lorentz_right;
        duration_fit_exponential;
        fall_times;
        fall_fit;
        rise_times;
        rise_fit;
        middle_point;
        exp_rise_values;
        exp_fall_values;
    end
    
    properties (Access = private)

        time_vector;
        max_time_deviation;
    end
    properties( Constant)
        fs = 187;
        ts = 1/AM_SR_peak.fs;
        estimate_floor_value = 0;
        min_left_position = 1;
        min_right_position = 2;
        max_samples_deviation = 45;
        maximum_lorentz_width = 0.8;
        
        
    end
    
    methods
        function obj = AM_SR_peak(data_envelope,start_sample,data_hilbert, data_original)
            %UNTITLED2 Construct an instance of this class
            %   Detailed explanation goes here
            obj.data_envelope = data_envelope * 1e-6;
            
            obj.start_time = (start_sample - 1) * AM_SR_peak.ts;
            
            obj.duration_simple = length(data_envelope) * AM_SR_peak.ts;
            
           obj.increment_value = max(obj.data_envelope) - obj.data_envelope(1);
        
           obj.max_value = max(obj.data_envelope);
           
           obj.time_vector = [1:(length(obj.data_envelope) )] * AM_SR_peak.ts;
           
           obj.max_time = obj.time_vector(find(obj.data_envelope == obj.max_value));
           
           obj.max_time_deviation = AM_SR_peak.max_samples_deviation * AM_SR_peak.ts;
           
           obj.get_x_data_y_data();
           
           obj.get_poly2();
           
           obj.get_lorentz();
           
           obj.data_hilbert = data_hilbert;
           
           obj.data_original = data_original;
           
           %obj.find_abrupt();
           
           
            
        end
        
        function get_x_data_y_data(obj)
            
                if(obj.max_time - obj.max_time_deviation < 0)
                    left_time = 0;
                else
                    left_time = obj.max_time - obj.max_time_deviation;
                end
                
                if((obj.max_time + obj.max_time_deviation) > (length(obj.data_envelope) * AM_SR_peak.ts))
                    right_time = length(obj.data_envelope) * AM_SR_peak.ts;
                else
                   right_time =  obj.max_time + obj.max_time_deviation;
                end
                
                obj.x_data = obj.time_vector(and((obj.time_vector > left_time), (obj.time_vector < right_time)));
                obj.y_data = obj.data_envelope(and((obj.time_vector > left_time), (obj.time_vector < right_time)));
        end
        
        function get_poly2(obj)
  
                obj.ftpoly2 = polyfit(obj.x_data,obj.y_data,2);
                obj.duration_poly2();
 
        end
        
        function get_lorentz(obj)
                %%Split x_data and y_data
                x_data_left = obj.x_data(1:find(obj.x_data == obj.max_time));
                y_data_left = obj.y_data(1:find(obj.x_data == obj.max_time));
                
                obj.ft_lorentz_left = lorentz_estimation(x_data_left,y_data_left, obj.max_time);
                left_width = obj.ft_lorentz_left.C1;
                
                x_data_right = obj.x_data(find(obj.x_data == obj.max_time):length(obj.x_data));
                y_data_right = obj.y_data(find(obj.x_data == obj.max_time):length(obj.x_data));
                obj.ft_lorentz_right = lorentz_estimation(x_data_right,y_data_right, obj.max_time);
                right_width = obj.ft_lorentz_right.C1;
                
                if ((left_width > AM_SR_peak.maximum_lorentz_width) || (right_width > AM_SR_peak.maximum_lorentz_width))
                    obj.duration_lorentz = 2 * min([abs(right_width),abs(left_width)]);
                else
                   obj.duration_lorentz = abs(right_width) + abs(left_width);
                end
                
                
        end
        
        function duration_poly2(obj)
                estimated_limits = roots(obj.ftpoly2) ;
                obj.duration_p2 = abs((max(estimated_limits) - min(estimated_limits)))  ;
        end
       function duration_half_power_poly2(obj)
                estimated_limits = roots(p) ;
                obj.duration_p2 = (max(estimated_limits) - min(estimated_limits))  ;
       end
       function duration_half_power_lorentz(obj)

       end
        
       function plot_envelope(obj)
           plot(obj.x_data, obj.y_data);
       end
       function plot_select_segment(obj)
           plot(obj.time_vector,obj.data_original,'LineWidth', 1);
           hold on;
           plot(obj.time_vector,obj.data_hilbert,'LineWidth', 2);
           hold off;
       end
       function vector = get_time_vector(obj)
           vector = obj.time_vector;
       end
       
       function find_abrupt(obj)
           
           [test_peak, ~] = envelope(obj.data_original,90,'analytic');
           %findchangepts(test_peak,'MaxNumChanges',1,'Statistic','linear');
           
           hold on
           %plot(obj.data_original,'b:','LineWidth', 1);

           
           [ipt, error] = findchangepts(test_peak,'MaxNumChanges',1,'Statistic','linear');
           half_window = 25;
           if  (((isempty(ipt) == 0) && ((ipt - half_window) >= 1) && ( (ipt + half_window) <= length(test_peak))))
               
               start_point = (ipt - half_window);
               pre_rise_window = test_peak((ipt - half_window):(ipt+half_window));
               max_rise_sample = find(pre_rise_window == max(pre_rise_window));
               rise_window = pre_rise_window(1:max_rise_sample);
               t_rise_window = [1:length(rise_window)] * AM_SR_peak.ts;
               [obj.rise_fit,rise_gof] = fit_exp(t_rise_window,rise_window,length(rise_window)* (AM_SR_peak.ts),rise_window(end));
               obj.exp_rise_values = obj.rise_fit(t_rise_window);
               obj.rise_times = obj.start_time + (start_point * AM_SR_peak.ts) + t_rise_window;
               obj.middle_point = (max_rise_sample * AM_SR_peak.ts)  + obj.start_time ;
               


               fall_window = test_peak(max_rise_sample:end);
               t_fall_window = [1:length(fall_window)] * AM_SR_peak.ts;
               [obj.fall_fit,fall_gof] = fit_exp(t_fall_window,fall_window,t_fall_window(1),fall_window(1));
               obj.exp_fall_values = obj.fall_fit(t_fall_window);
               %t_plot_fall = [middle_point:(length(test_peak))];
               %values_fit = fall_fitter(1:length(fall_window));
               %plot(t_plot_fall, values_fit,'Color',	'#77AC30', 'LineWidth', 3);
               obj.fall_times = obj.middle_point + (start_point * AM_SR_peak.ts)  +  t_fall_window;
               obj.duration_fit_exponential = (length(rise_window) + length(fall_window)) * (AM_SR_peak.ts);
           else
               obj.duration_fit_exponential = 0;
           end
       end
       function property_fit_heilder = fit_heilder(obj)
           t_length = (max(obj.time_vector) - min(obj.time_vector));
           t_start = min(obj.time_vector);
           d_max = max(obj.data_envelope);
           t_normalized = (obj.time_vector - t_start) / t_length;
           d_normalized = obj.data_envelope / d_max;
           [fitresult, gof] = fit_heilder(t_normalized, d_normalized);
           property_fit_heilder = fitresult;
       end
       function plot(obj,option)
           if strcmp("fit2", option)
               plot(obj.fall_times,(obj.exp_fall_values)',obj.rise_times,(obj.exp_rise_values'));
           end
       end
    end
end

