function [u07, u14, u21] = get_velocities_fromC(C)
% INPUT: 'C' track results for a location
% OUTPUT: [u07, u14, u21] the velocities at different heights in mm/s
volts = ["2_", "2.5", "3_", "3.5", "4_", "4.5", "5_", "5.5", "6_"];
volts_double = [2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6];
Periods = [100,80,66,57,50,44,40,36,33,28];
Periods_string = ["100","80","66","57","50","44","40","36","33","28"];
freqs = 1000./Periods;

u07 = zeros(numel(volts), numel(freqs));
u14 = zeros(numel(volts), numel(freqs));
u21 = zeros(numel(volts), numel(freqs));
for ii=1:numel(volts)
    v = strcat('V', volts(ii));
    for jj=1:numel(freqs)
        p = strcat('P', Periods_string(jj));
        for kk=1:numel(C)
            if contains(C{kk}.filename, v) && contains(C{kk}.filename, p)
                if contains(C{kk}.filename, 'Z14')
                    %Include 2*pi from the formula u_max = A * angular_freq
                    %= A * frequency * 2pi
                    % A1 is in pixels, 0.14um is the respective distance
                    % per pixel.
                    u14(ii,jj) = C{kk}.A1*0.14*freqs(jj)*2;  
                elseif contains(C{kk}.filename, 'Z07')
                    u07(ii,jj) = C{kk}.A1*0.14*freqs(jj)*2;
                elseif contains(C{kk}.filename, 'Z21')
                    u21(ii,jj) = C{kk}.A1*0.14*freqs(jj)*2;
                end
            end
        end
    end
end

% Set units of u to mm/s (before are um/s)
u07 = u07./1000;
u14 = u14./1000;
u21 = u21./1000;

end