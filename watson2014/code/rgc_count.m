%compute comulative counts as a function of eccentricity
%integrate density at each eccentricity r multiplied by 2pir to account for increasing area
clear all;  
dataLoad_all;
x_nameIdx = {'xi_deg_rgcd','xn_deg_rgcd','xt_deg_rgcd','xs_deg_rgcd'}; 
y_nameIdx = {'yi_deg_rgcd','yn_deg_rgcd','yt_deg_rgcd','ys_deg_rgcd'}; 
%Since the number of data points provided is not equal across all retinal locations
%we need to account for this:
for ii = 1:length(x_nameIdx)
    temp = round(eval(x_nameIdx{ii}));
    eccIdx(ii) = find(temp==19);
    clear temp;
end
%if each foveal cone drives exactly two midget retinal ganglion cells
%then we must also scale the final output by 2*1.12 (the 1.12 - Drasdo et al., 2007)
%integrate at each specified ecc(r)*2pir = data*2pir^2
for jj = 1:length(y_nameIdx)
    y_data = eval(y_nameIdx{jj});
    x_data = eval(x_nameIdx{jj});
    pointInt = cumsum(y_data(1:eccIdx(jj)).*(2*pi*x_data(1:eccIdx(jj)).^2))/10000;
    pointInt = interp(pointInt,10);
    locInt(:,jj) = pointInt*(2*1.12); %scaling factor taken from paper
    clear pointInt
end






