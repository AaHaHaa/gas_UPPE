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

c = 299792458; % m/s

num_gas = length(gas.material);
for gas_i = 1:num_gas
    switch gas.material{gas_i}
        case 'H2'
            % From Herring et al., "Temperature and density dependence of the linewidths and line shifts of the rotational Raman lines in N2 and H2," Phys. Rev. A 34, 1944-1951 (1986)
            gas.H2.R.T2 = 1e6/(pi*(6.15/eta(gas_i)+114*eta(gas_i))); % ps
            % From Bischel and Dyer, "Temperature dependence of the Raman linewidth and line shift for the Q(1) and Q(0) transitions in normal and para-H_2," Phys. Rev. A 33, 3113-3123 (1986)
            gas.H2.V.T2 = 1e6/(pi*(309/eta(gas_i)*(gas.temperature/298)^0.92+(51.8+0.152*(gas.temperature-298)+4.85e-4*(gas.temperature-298)^2)*eta(gas_i))); % ps
        case 'D2'
            % From Hanson and Poirier, "Stimulated rotational Raman conversion in H_2, D_2, and HD," IEEE J. Quantum Electron. 29, 2342-345 (1993)
            % It contains only the 134 factor.
            % From H2's vib.T2 relation, it seems that 51.8 --> 120 is accompanied by 309 --> 101 with approximately the same factor, possibly due to the inverse eta relation.
            % From this speculation, I scale D2's rot.T2 accordingly.
            % H2's 114 --> D2's 134 leads to H2's 6.15 --> D2's 6.15*(134/114)=7.2, so I put 7 below.
            gas.D2.R.T2 = 1e6/(pi*(7/eta(gas_i)+134*eta(gas_i))); % ps
            gas.D2.V.T2 = 1e6/(pi*(101/eta(gas_i)+120*eta(gas_i))); % ps
        case {'N2','O2','air'}
            if ismember(gas.material{gas_i},{'N2','air'})
                % From Table 1 of
                % Miller et al., "Communication: Time-domain measurement of high-pressure N2 and O2 self-broadened linewidths using hybrid femtosecond/picosecond coherent anti-Stokes Raman scattering," J. Chem. Phys. 135(20), 201104 (2011)
                % 0.098 cm^(-1)/atm is converted to 2938 MHz/atm (room temperature) = 3207 MHz/amagat
                %
                % Other reference:
                % Fig.9 of Lempert et al., "Stimulated Raman scattering and coherent anti-Stokes Raman spectroscopy in high-pressure oxygen," J. Opt. Soc. Am. B 7(5), 715-721 (1990)
                gas.N2.R.T2 = 1e6/(pi*3207*eta(gas_i)); % ps
                % From Fig.9 of
                % Jammu et al., "PRESSURE BROADENING OF THE ROTATIONAL RAMAN LINES OF SOME SIMPLE GASES," Can. J. Phys. 44, 797-814 (1966)
                % 0.1 cm^(-1)/atm is converted to 2998 MHz/atm (room temperature) = 3273 MHz/amagat
                gas.N2.V.T2 = 1e6/(pi*3273*eta(gas_i)); % ps

                % MEG model for finding linewidth/T2 and strength ratio
                if eta(gas_i) ~= 0
                    try
                        if eta(gas_i) < 0.05 % to avoid ignoring Raman when eta<0.05 whose rounded_eta=0
                            rounded_eta = 0.1;
                        else
                            rounded_eta = round(gather(eta(gas_i)),1); % MATLAB round() has some problem with gpuArray
                        end
                        eta_idx = find(abs(gas.N2.V.T2_data.eta - rounded_eta) < 0.05,1);
                    catch % if there is no T2 data, compute it
                        eta_idx = [];
                    end
                    if ~isempty(eta_idx)
                        gas.N2.V.T2 = gas.N2.V.T2_data.T2(eta_idx,:);
                    else
                        % % From Rahn and Palmer, "Studies of nitrogen self-broadening at high temperature with inverse Raman spectroscopy," J. Opt. Soc. Am. B 3(9), 1164-1169 (1986)
                        gas.N2.V.MEG = struct('alpha', 0.0231*1e2*c*1e-12*2*pi,... % 2*pi*THz/atm
                                               'beta', 1.67,...
                                               'delta', 1.21,...
                                               'a', 1.5,...
                                               'm', 0.1487); % from Rahn et al., "Comparison of rotationally inelastic collision models for Q-branch Raman spectra of N2," Chem. Phys. Lett. 133(6), 513-516 (1987)
                        gas.N2.V.T2 = find_Raman_linewidth(gas,gas.material{gas_i},eta(gas_i),gas.N2.R,gas.N2.V);
                    end
                end
            end
            if ismember(gas.material{gas_i},{'O2','air'})
                % From Table 2 of
                % Miller et al., "Communication: Time-domain measurement of high-pressure N2 and O2 self-broadened linewidths using hybrid femtosecond/picosecond coherent anti-Stokes Raman scattering," J. Chem. Phys. 135(20), 201104 (2011)
                % 0.107 cm^(-1)/atm is converted to 3207 MHz/atm (room temperature) = 3500 MHz/amagat
                gas.O2.R.T2 = 1e6/(pi*3500*eta(gas_i)); % ps
                % From Fig.9 of
                % Fletcher et al., "High resolution vibrational Raman spectrum of oxygen," J. Raman Spectrosc. 2, 3-14 (1974)
                % 0.146 cm^(-1)/atm is converted to 4377 MHz/atm (room temperature) = 4778 MHz/amagat
                gas.O2.V.T2 = 1e6/(pi*4778*eta(gas_i)); % ps

                % MEG model for finding linewidth/T2 and strength ratio
                if eta(gas_i) ~= 0
                    try
                        if eta(gas_i) < 0.05 % to avoid ignoring Raman when eta<0.05 whose rounded_eta=0
                            rounded_eta = 0.1;
                        else
                            rounded_eta = round(gather(eta(gas_i)),1); % MATLAB round() has some problem with gpuArray
                        end
                        eta_idx = find(abs(gas.N2.V.T2_data.eta - rounded_eta) < 0.05,1);
                    catch % if there is no T2 data, compute it
                        eta_idx = [];
                    end
                    if ~isempty(eta_idx)
                        gas.O2.V.T2 = gas.N2.V.T2_data.T2(eta_idx,:);
                    else
                        % From Lempert et al., "Stimulated Raman scattering and coherent anti-Stokes Raman spectroscopy in high-pressure oxygen," J. Opt. Soc. Am. B 7, 715-721 (1990)
                        gas.O2.V.MEG = struct('alpha', 0.01528*1e2*c*1e-12*2*pi,... % 2*pi*THz/atm
                                              'beta', 1.362,...
                                              'delta', 1.786,...
                                              'a', 0.50,...
                                              'm', 0.1487); % assumed the same as N2; from Rahn et al., "Comparison of rotationally inelastic collision models for Q-branch Raman spectra of N2," Chem. Phys. Lett. 133(6), 513-516 (1987)
                        gas.O2.V.T2 = find_Raman_linewidth(gas,gas.material{gas_i},eta(gas_i),gas.O2.R,gas.O2.V);
                    end
                end
            end
        case 'N2O' % we assume that N2O follows the same relation as N2
            gas.(gas.material{gas_i}).R.T2 = 1e6/(pi*3570*eta(gas_i)); % ps
        case 'CO2'
            gas.(gas.material{gas_i}).R.T2 = 1e6/(pi*3900*eta(gas_i)); % ps
        case 'CH4'
            gas.CH4.V.T2 = 1e6/(pi*(8220+384*eta(gas_i))); % ps
    end
end

end