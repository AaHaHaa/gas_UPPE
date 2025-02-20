function gas = Raman_T2( gas,eta )
%RAMAN_T2 It calculates the Raman dephasing time of several materials.
%   Input:
%       gas: a structure containing
%           gas.material
%           gas.pressure
%           gas.temperature
%       eta: gas density (amagat)
%
%   Output:
%       gas

switch gas.material
    case 'H2'
        % T2
        gas.H2.R.T2 = 1e6/(pi*(6.15/eta+114*eta)); % ps
        gas.H2.V.T2 = 1e6/(pi*(309/eta*(gas.temperature/298)^0.92+(51.8+0.152*(gas.temperature-298)+4.85e-4*(gas.temperature-298)^2)*eta)); % ps
    case {'N2','O2','air'}
        if ismember(gas.material,{'N2','air'})
            % T2
            gas.N2.R.T2 = 1e6/(pi*3570*eta); % ps
            gas.N2.V.T2 = 1e6/(pi*22.5); % ps
        end
        if ismember(gas.material,{'O2','air'})
            % T2
            gas.O2.R.T2 = 1e6/(pi*1701*eta); % ps
            gas.O2.V.T2 = 1e6/(pi*54); % ps
        end
    case 'CH4'
        gas.CH4.V.T2 = 1e6/(pi*(8220+384*eta)); % ps
end

end