classdef AM_SR_ray_vector
    %AM_SR_RAY_VECTOR Summary of this class goes here
    %   Detailed explanation goes here
   
    properties(Constant)
        bins = 50;
    end
    methods (Static)
        function histogram(AM_SR_ray_array)
            subplot(3,3,1);
            selected_rays = AM_SR_ray_vector.select_ray(AM_SR_ray_array);
            
            hist([selected_rays.peak_value],AM_SR_ray_vector.bins);
            title("Peak Values");
            subplot(3,3,2);
            hist(1./[selected_rays.FWHMT],AM_SR_ray_vector.bins);
            title("Freq (1 / FWHMT)");
            subplot(3,3,3);
            hist([selected_rays.RiseTime],AM_SR_ray_vector.bins);
            title("Rise Time");
            subplot(3,3,4);
            hist([selected_rays.FWHMT],AM_SR_ray_vector.bins);
            title("FWHMT");
            subplot(3,3,5);
            hist([selected_rays.t_fall]./[selected_rays.t_rise], AM_SR_ray_vector.bins);
           
            title("Fall Time / Rise Time")
            subplot(3,3,6);
            hist([selected_rays.T1],AM_SR_ray_vector.bins);
            title("Time UP: 0 to Max");
            subplot(3,3,7);
            hist([selected_rays.T2],AM_SR_ray_vector.bins);
            title("Time DOWN: max to 10%");
            subplot(3,3,8);
            hist([selected_rays.T2]./[selected_rays.T1], AM_SR_ray_vector.bins);
           title(" Fall DOWN / Rise UP");
             subplot(3,3,9);
            hist([selected_rays.T2] + [selected_rays.T1], AM_SR_ray_vector.bins);
           title("Time DOWN + Time UP");
           
        end
        
        function scatter(AM_SR_ray_array)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            selected_rays = AM_SR_ray_vector.select_ray(AM_SR_ray_array);
            subplot(3,2,1);
            scatter([selected_rays.peak_value],[selected_rays.RiseTime]);
            title("Peak Value vs Time Rise");
            subplot(3,2,2);
            scatter([selected_rays.peak_value],[selected_rays.FWHMT]);
            title("Peak Value vs FWHMT");
            subplot(3,2,3);
            scatter([selected_rays.T1],[selected_rays.T2]);
            title("Time Rise vs Time Fall");
            subplot(3,2,4);
            scatter([selected_rays.peak_value],[selected_rays.T2]./[selected_rays.T1]);
            title("Peak Value vs Relation between Time fall and Rise");
            figure();
            scatter([selected_rays.peak_value],[selected_rays.SPC]);
            
            title("Peak Value vs PSD");
        end
        function scatter_outliers_PSD(AM_SR_ray_array_raw)
            AM_SR_ray_array = AM_SR_ray_vector.select_ray(AM_SR_ray_array_raw);
            
            PSDs = AM_SR_ray_vector.get_PSD(AM_SR_ray_array);
            TF = isoutlier(PSDs, 'grubbs');
            display(sum(TF));
            AM_SR_ray_array_outliers = AM_SR_ray_array(TF);
            AM_SR_ray_array_not_outliers = AM_SR_ray_array(~TF);
            std_not_outlier = std(PSDs(~TF));
            norm_PSDs = PSDs - mean(PSDs);
            
            %Scatter not outliers
            scatter([AM_SR_ray_array_not_outliers.peak_value], norm_PSDs(~TF)./std_not_outlier);
            hold on
            
            %Extreme Outliers
            outliers_psds = norm_PSDs(TF);
            
            
            TF = (norm_PSDs(TF)./std_not_outlier > 100);
            display("Extreme Outliers: " + sum(TF));
            
            AM_SR_ray_array_extreme_outliers = AM_SR_ray_array_outliers(TF);
            AM_SR_ray_array_rest_outliers = AM_SR_ray_array_outliers(~TF);
            
            scatter([AM_SR_ray_array_rest_outliers.peak_value], outliers_psds(~TF)./std_not_outlier,'d','MarkerEdgeColor',[.5 .2 .2]);
            
            %Scatter outliers
            scatter([AM_SR_ray_array_extreme_outliers.peak_value], outliers_psds(TF)./std_not_outlier,'x','MarkerEdgeColor',[1 .2 .2]);
           hold off 
        end
    
        function selected_rays = select_ray(AM_SR_ray_array)
                select = arrayfun(@(x) x.is_ray(),AM_SR_ray_array);
                selected_rays = AM_SR_ray_array(select);
        end
        function selected_non_rays = select_non_ray(AM_SR_ray_array)
                select = arrayfun(@(x) x.is_ray(),AM_SR_ray_array);
                selected_non_rays = AM_SR_ray_array(~select);
        end
        function  selected_qburst = separate_qburst(AM_SR_ray_array)
              select = arrayfun(@(x) x.is_qburst(),AM_SR_ray_array);
              selected_qburst = AM_SR_ray_array(select);
        end
        function values = get_PSD(AM_SR_ray_array)
            values  = arrayfun(@(x) x.get_SPC(),AM_SR_ray_array);
        end
        function instant_freq = get_instant_freq(AM_SR_ray_array, centre_point)
            freq_carrier = 7.8;
            freq_mod = 1.5;
            
            carrier_instant = arrayfun(@(x) x.instant_freq(freq_carrier),AM_SR_ray_array);
            mod_instant = arrayfun(@(x) x.instant_freq(freq_mod),AM_SR_ray_array);
            figure(1)
            scatter(carrier_instant, mod_instant);
            figure(2)
            scatter([AM_SR_ray_array.peak_value], mod_instant);
            figure(3)
            scatter([AM_SR_ray_array.peak_value], carrier_instant);
            instant_freq = carrier_instant;
            
            psd = AM_SR_ray_vector.get_PSD(AM_SR_ray_array);
            figure(4);
            scatter(psd, carrier_instant);
            figure(5);
            scatter(psd, mod_instant);
        end
        function export_r(AM_SR_ray_array)
            
            delete R/r_am_sr.mat;
            
            data.start_time = [AM_SR_ray_array.start_time]';
            
            data.duration = [AM_SR_ray_array.duration]';

            data.t_rise = [AM_SR_ray_array.t_rise]';

            data.t_fall = [AM_SR_ray_array.t_fall]';

            data.nh = [AM_SR_ray_array.nh]';

            data.Io = [AM_SR_ray_array.Io]';

            data.peak_value = [AM_SR_ray_array.peak_value]';

            data.peak_location = [AM_SR_ray_array.peak_location]';

            data.rmse = [AM_SR_ray_array.rmse]';

            data.rsquare = [AM_SR_ray_array.rsquare]';

            data.T1 = [AM_SR_ray_array.T1]';

            data.T2 = [AM_SR_ray_array.T2]';

            data.RiseTime = [AM_SR_ray_array.RiseTime]';

            data.FWHMT = [AM_SR_ray_array.FWHMT]';

            data.instant_freq_carrier = [AM_SR_ray_array.instant_freq_carrier]';

            data.instant_freq_mod =  [AM_SR_ray_array.instant_freq_mod]';

            data.PSD_in_band = [AM_SR_ray_array.PSD_in_band]';
            
            data.PSD_total = [AM_SR_ray_array.PSD_total]';
            
            data.PSD_rest = [AM_SR_ray_array.PSD_rest]';
            
            save R/r_am_sr.mat data;
        end
    end
end

