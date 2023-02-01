
x = data_filter;
t = tx;
fl1 = 12;
[up1,lo1] = envelope(x,fl1,'analytic');
fl2 = 30;
[up2,lo2] = envelope(x,fl2,'analytic');
fl3 = 300;
[up3,lo3] = envelope(x,fl3,'analytic');
param_small = {'Color',[0.9 0.4 0.1],'Linewidth',2};
param_large = {'Color',[0 0.4 0],'Linewidth',2};
param_3 = {'Color',[0 0.9 0.3],'Linewidth',2};

plot(t,x)
hold on
%p1 = plot(t,up1,param_small{:});
%plot(t,lo1,param_small{:});
p2 = plot(t,up2,param_large{:});
plot(t,lo2,param_large{:});
p3 = plot(t,up3,param_3{:});
plot(t,lo3,param_3{:});
hold off

%legend([p1 p2 p3],'fl = 12','fl = 30', 'fl = 50')
ylim([-9e6,9e6]);
xlim([294,302]);
title('Analytic Envelope')