function [sse,fit] = fit_displacementZoneNasalTemporal(params,y,ecc,type,test)
    %this function fits values outside the displacement zone for Nasal and
    %Temporal formula
    %minimizes sse
    % had to fix 3 parameters because fmincon was getting lost otherwise...
    % Not sure why...
    if test == 1
        if type == 2
            gamma_t = 0.77754;
        elseif type == 1
            gamma_t = 0.91565;
        end
        mu = params(2);
        if type == 2
            beta = 1.746;
        elseif type == 1
            beta = 2.4598;
        end
        delta = params(4); %gain
        if type == 2
            alpha = 2.4607;
        elseif type == 1
            alpha = 1.8938;
        end
    else
        % in case everything should be free params
        gamma_t = params(1);
        mu = params(2);
        beta = params(3);
        delta = params(4);
        alpha = params(5);
    end
    %
    fit = delta*((gamma_t.*exp(-((ecc-mu)/beta).^gamma_t)).*(((ecc-mu)/beta).^(alpha*gamma_t-1))/(beta.*gamma(alpha)));
    sse = sum(y-real(fit))^2;
    
    %% paper params ...for testing
%     gamma_t = 0.91;
%     mu = -0.09;
%     beta = 2.45; %2.2 works (2-3) range
%     delta = 14.9; %gain
%     alpha = 1.89;   
