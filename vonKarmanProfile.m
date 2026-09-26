function [x, y] = vonKarmanProfile(noseLength, baseDiameter, options)
    % vonKarmanProfile Coordinates of the von Karman (LD-Haack) nose cone profile.
    %
    %   [x, y] = vonKarmanProfile(noseLength, baseDiameter) returns the axial
    %   position x, measured from the tip, and the local radius y of the Haack
    %   series profile with C = 0. This is the profile of minimum linearized
    %   wave drag for a given length and base diameter.
    %
    %   [x, y] = vonKarmanProfile(..., NumPoints=n) sets the number of points
    %   (default 401). Points are equally spaced in the Haack parameter
    %   theta = acos(1 - 2*x/noseLength), which clusters them near the tip
    %   where the curvature is largest.
    %
    %   Inputs
    %       noseLength   - nose cone length, same units as baseDiameter
    %       baseDiameter - diameter at the base of the nose cone
    %
    %   Outputs
    %       x - row vector, 0 at the tip and noseLength at the base
    %       y - row vector, 0 at the tip and baseDiameter/2 at the base
    %
    %   Profile equations (Haack series):
    %       theta = acos(1 - 2*x/L)
    %       y = (R/sqrt(pi))*sqrt(theta - sin(2*theta)/2 + C*sin(theta)^3), C = 0

    arguments
        noseLength (1, 1) double {mustBePositive}
        baseDiameter (1, 1) double {mustBePositive}
        options.NumPoints (1, 1) double {mustBeInteger, mustBeGreaterThanOrEqual(options.NumPoints, 3)} = 401
    end

    baseRadius = baseDiameter/2;
    theta = linspace(0, pi, options.NumPoints);

    % haack series with C = 0; the max guards against round-off at the tip
    x = (noseLength/2)*(1 - cos(theta));
    y = (baseRadius/sqrt(pi))*sqrt(max(theta - sin(2*theta)/2, 0));
end
