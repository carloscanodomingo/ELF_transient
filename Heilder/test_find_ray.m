



data = load("c_2013_11_16_04.dat");

%data = repmat(data);
fs = 187;
ts = 1/187;

total_time = ts * length(data);
tx = linspace(0,length(data) * ts, length(data));
tx_lim = [210, 220.5];

selected = select_chunk(tx, tx_lim);

tx_test = tx;%(selected);
data_test = data';%(selected);

[AM_SR_ray_v] = find_ray(tx_test, data_test);

    function selected  = select_chunk(x, x_lim)

        selected = (x > x_lim(1) & x < x_lim(2));

    end
