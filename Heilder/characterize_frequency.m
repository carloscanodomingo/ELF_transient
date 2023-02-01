frequency = [];
for i = 2016:2017
    year_data = SR_data.load_year(i,"NS",SR_config.current_version);
    for index_month = 1:12
        current_month = year_data{index_month};
        select = SR_peak_process_array.select_not_noisy(current_month, "NS", 3);
        not_noity = current_month(select);
        
        current_month_lorentz = [not_noity.SR_classification];
        current_month_freqs = [current_month_lorentz.FITLORENTZ_first_lorentz];
        
        frequency = [frequency, current_month_freqs];
    end

end