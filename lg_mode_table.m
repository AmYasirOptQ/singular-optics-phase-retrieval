function [Lin, Pin] = lg_mode_table(num)
% LG_MODE_TABLE  Build the table of LG mode indices for the half-integer
%                j hierarchy used in the phase retrieval paper.
%
%   [Lin, Pin] = lg_mode_table(num)
%
%   Input
%     num : highest level index; j runs from 0 to num/2 in steps of 1/2.
%           num = 5 gives j in {0, 1/2, 1, 3/2, 2, 5/2} -> 21 modes.
%
%   Outputs
%     Lin : azimuthal indices l = 2m   (1 x num0 vector)
%     Pin : radial indices    p = j-|m| (1 x num0 vector)
%
%   Physics
%     The hierarchy follows the SU(2) representation structure.
%     Each level j contributes (2j+1) modes with m = -j, ..., +j.
%     The LG indices are:
%
%       l = 2m          (azimuthal; can be fractional for half-integer j)
%       p = j - |m|     (radial;    always a non-negative integer)
%
%     Total number of modes: num0 = (num+1)(num+2)/2.

Min = 0;
Jin = 0;

for n = 1:num
    j   = 0.5*n;
    Min = [Min,  j:-1:-j]; 
    Jin = [Jin,  j*ones(1, n+1)]; 
end

Lin = 2*Min;
Pin = Jin - abs(Min);

end
