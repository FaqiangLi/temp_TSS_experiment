function TSS_distance_histogram(theta,inp)

xp = inp.xp;
ix_JJ2C = inp.ix2;
NT = inp.NT;
J = inp.J;

[~,~,P_pred] = TSS_quantities(theta,inp);
P_obs = inp.storepairvisited;

dist = xp(:,:,:,2);
dist = reshape(dist,NT,J*J);
dist = dist(:,ix_JJ2C);
dist_pred = sum(P_pred.*dist,2);
dist_obs = sum(P_obs.*dist,2);

stores = 1 + xp(:,:,:,1);
stores = reshape(stores,NT,J*J);
stores = stores(:,ix_JJ2C);
stores_pred = sum(P_pred.*stores,2);
% fix rounding error:
ix_below = stores_pred<1;
ix_above = stores_pred>2;
stores_pred(ix_below) = 1;
stores_pred(ix_above) = 2;
stores_obs = sum(P_obs.*stores,2);

subplot(2,2,1)
h1 = histogram(dist_pred);
h1.Normalization = 'probability';
h1.BinWidth = 1;
titlestring = {'Distance travelled (predicted)',};
xaxisstring = 'Kilometres';
title(titlestring)
xlabel(xaxisstring)
axis([0 100 0 0.095])

subplot(2,2,2)
h2 = histogram(dist_obs);
h2.Normalization = 'probability';
h2.BinWidth = 1;
titlestring = {'Distance travelled (observed)',};
xaxisstring = 'Kilometres';
title(titlestring)
xlabel(xaxisstring)
axis([0 100 0 0.095])

subplot(2,2,3)
h3 = histogram(stores_pred);
h3.Normalization = 'probability';
h3.BinWidth = .05;
titlestring = {'No. stores visited (predicted)',};
xaxisstring = 'stores';
title(titlestring)
xlabel(xaxisstring)
axis([1 2 0 0.8])

subplot(2,2,4)
h4 = histogram(stores_obs);
h4.Normalization = 'probability';
h4.BinWidth = .05;
titlestring = {'No. stores visited (observed)',};
xaxisstring = 'stores';
title(titlestring)
xlabel(xaxisstring)
axis([1 2 0 0.8])

filename = ['Output/' 'Distance_histogram','_',datestr(now, 'ddmm_yyyy_HHMM')];
print(filename,'-djpeg')
print(filename,'-deps')
print(filename,'-dpdf')


