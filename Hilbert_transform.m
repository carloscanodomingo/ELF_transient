load train

findchangepts(y,'MaxNumChanges',10,'Statistic','rms')
rng('default')

lr = 20;

mns = [0 1 4 -5 2 0 1];
nm = length(mns);

vrs = [1 4 6 1 3];
nv = length(vrs);

v = randn(1,lr*nm*nv)/2;

f = reshape(repmat(mns,lr*nv,1),1,lr*nm*nv);
y = reshape(repmat(vrs,lr*nm,1),1,lr*nm*nv);

t = v.*y+f;

subplot(2,2,1)
plot(v)
title('Original')
xlim([0 700])

subplot(2,2,2)
plot([f;v+f]')
title('Means')
xlim([0 700])

subplot(2,2,3)
plot([y;v.*y]')
title('Variances')
xlim([0 700])

subplot(2,2,4)
plot(t)
title('Final')
xlim([0 700])