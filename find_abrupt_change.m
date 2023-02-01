
if (exist("data","var") == 0)
    data = (load("c_2013_11_16_04.dat"))';
end

fs = 187;
ts = 1/187;
tx = linspace(0,length(data) * ts, length(data));

N = length(data);
data_filter = preprocess(data,tx);

fl3 = 300;
[up3,lo3] = envelope(data_filter,fl3,'analytic');

step = 187*2;
semi_step = step / 2;
n_chunck = round(N/(step)) - 1;

factor_max_num_changes = 2;



for i = 1:n_chunck
    start_point = (step * (i - 1) + 1);
    end_point = (step * i);
    envelope_chunck = up3(start_point:end_point);
    max_num_changes = round((length(envelope_chunck) * factor_max_num_changes) / fs);
    findchangepts(envelope_chunck,'MaxNumChanges', max_num_changes,'Statistic','linear');
    [ipt,residual] = findchangepts(envelope_chunck,'MaxNumChanges', max_num_changes,'Statistic','linear');
    
    
    
    
end
data = up3;