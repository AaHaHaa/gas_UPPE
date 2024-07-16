function ZdYd = calc_ZdYd(n_gas, n_out, mode)
%CALC_ZDYD Summary of this function goes here
%   Detailed explanation goes here

relative_n_out = complex(n_out./n_gas);

switch mode
    case 'TE'
        Zd = 1./sqrt(relative_n_out.^2-1); % impedance
        ZdYd = Zd;
    case 'TM'
        Yd = relative_n_out.^2./sqrt(relative_n_out.^2-1); % admittance
        ZdYd = Yd;
    case 'EH'
        % 1/2*(Zd+Yd)
        ZdYd = 1/2*(relative_n_out.^2+1)./sqrt(relative_n_out.^2-1);
end

end