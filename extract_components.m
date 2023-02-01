function [SC,f] = extract_components(raw_data, start_frequency,end_frequency, BW)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    validateattributes(start_frequency,{'numeric'},{'>=',0.1,'<=',end_frequency});
    %validateattributes(start_frequency - floor(start_frequency), {'numeric'}, {'=',0});

    fs = 187;
    SC_components = (end_frequency - start_frequency) * (1/BW);
    SC = zeros(1,SC_components);



    for i=1:SC_components
        f_start = start_frequency + (BW * (i - 1));
        f_end = f_start + BW;
        data_filtered = band_pass_filter(raw_data, fs, f_start, f_end);
        SC(i) = rms(data_filtered).^2;
        f(i) = (f_start + f_end) / 2;
    end
end
