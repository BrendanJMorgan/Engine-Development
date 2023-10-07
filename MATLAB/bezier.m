function curve = bezier(points) % Generates a quartic bezier curve for four inputted control points. 
	% points = [startx, starty; startcontrolx, startcontroly; endcontrolx, endcontroly; endx, endy]
						
	binomialCoefficients = [1, 3, 3, 1];
   
    bezierPoints = []; % Initialize an empty array to store the Bezier curve points
    
    for u = 0:0.001:1 % Loop through parameter 'u' from 0 to 1 to generate points on the Bezier curve
        bernsteinBasis = zeros(1, 4); % Initialize an array to store the Bernstein basis values
        for d = 1:4 % Compute the Bernstein basis for each control point for the current value of 'u'
            bernsteinBasis(d) = binomialCoefficients(d) * ((1 - u)^(4 - d)) * (u^(d - 1));
		end
        bezierPoints = [bezierPoints; bernsteinBasis]; % Append the newly calculated Bernstein basis to the bezierPoints array
    end
	
    curve = bezierPoints * points;     % Curve points are equivalent to a weighted sum of control points

	% % Calculate Tangents at u=0 and u=1
    % T0 = 4*(points(2,:) - points(1,:));
    % T1 = 4*(points(4,:) - points(3,:));
	% 
    % % Normalize tangents for better control (optional)
    % T0 = T0 / norm(T0);
    % T1 = T1 / norm(T1);
	% 
    % % Adjust the first and last points on the curve
    % curve(1,:) = curve(1,:) + 0.001 * T0;
    % curve(end,:) = curve(end,:) + 0.001 * T1;
end