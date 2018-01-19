function [value,isterminal,direction] = event_retrograde(t,Y)

% Stops the function when ATPm reaches 0.5.
% This condition can be changed via the 'value' parameter

value = Y(2)-.5;  
isterminal = 1; % Stop the fcn
direction = 0;

end