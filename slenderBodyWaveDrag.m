function [dragCoefficient, fourierCoefficients] = slenderBodyWaveDrag(x, radius, options)
    % slenderBodyWaveDrag Linearized supersonic wave drag of a slender body of revolution.
    %
    %   cd = slenderBodyWaveDrag(x, radius) evaluates the Karman-Moore (slender
    %   body) wave drag of a body with cross-section area S(x) = pi*radius.^2
    %   and returns the drag coefficient referenced to the maximum cross-section
    %   area, which is the base area for a nose cone. The area slope must vanish
    %   at both ends, S'(0) = S'(L) = 0: a pointed tip, and either a closed tail
    %   or a base tangent to a following cylinder. The Haack series satisfies
    %   this. The linearized result does not depend on Mach number.
    %
    %   [cd, an] = slenderBodyWaveDrag(...) also returns the Fourier
    %   coefficients A_n of the area slope, S'(x) = sum(A_n*sin(n*theta)) with
    %   x = (L/2)*(1 - cos(theta)), from which D/q = (pi/4)*sum(n*A_n^2).
    %
    %   slenderBodyWaveDrag(..., NumTerms=n) sets the number of Fourier terms
    %   (default 20).
    %   slenderBodyWaveDrag(..., ReferenceArea=a) sets the reference area for
    %   the drag coefficient (default pi*max(radius)^2).
    %
    %   Inputs
    %       x      - axial positions from the tip, monotonically increasing
    %       radius - local radius at each x
    %
    %   Method
    %       A_n = (2/pi)*integral(S'(x(theta))*sin(n*theta), theta = 0..pi)
    %       With S'(x) = (dS/dtheta)/(dx/dtheta) = (2/L)*(dS/dtheta)/sin(theta):
    %       A_n = (4/(pi*L))*integral((dS/dtheta)*U_{n-1}(cos(theta)), 0..pi)
    %       where U_{n-1}(cos(theta)) = sin(n*theta)/sin(theta) is a Chebyshev
    %       polynomial of the second kind, evaluated by recurrence so that no
    %       division by sin(theta) occurs at the tip or the base.

    arguments
        x double {mustBeVector, mustBeNonnegative}
        radius double {mustBeVector, mustBeNonnegative}
        options.NumTerms (1, 1) double {mustBeInteger, mustBePositive} = 20
        options.ReferenceArea (1, 1) double {mustBePositive} = pi*max(radius)^2
    end

    x = reshape(x, 1, []);
    radius = reshape(radius, 1, []);
    if numel(x) ~= numel(radius)
        error("slenderBodyWaveDrag:sizeMismatch", ...
            "x and radius must have the same number of elements. Received %d and %d.", ...
            numel(x), numel(radius));
    end

    bodyLength = x(end) - x(1);

    % map x onto the haack parameter: theta = 0 at the tip, pi at the base
    theta = acos(1 - 2*(x - x(1))/bodyLength);
    area = pi*radius.^2;
    areaSlopeTheta = gradient(area, theta);

    % fourier coefficients of S'(x), using the chebyshev recurrence for sin(n*theta)/sin(theta)
    cosTheta = cos(theta);
    chebyshevPrevious = zeros(size(theta));     % U_{-1}
    chebyshevCurrent = ones(size(theta));       % U_0
    fourierCoefficients = zeros(1, options.NumTerms);
    for n = 1:options.NumTerms
        fourierCoefficients(n) = (4/(pi*bodyLength))*trapz(theta, areaSlopeTheta.*chebyshevCurrent);
        chebyshevNext = 2*cosTheta.*chebyshevCurrent - chebyshevPrevious;
        chebyshevPrevious = chebyshevCurrent;
        chebyshevCurrent = chebyshevNext;
    end

    % karman-moore wave drag D/q = (pi/4)*sum(n*A_n^2)
    termIndex = 1:options.NumTerms;
    dragArea = (pi/4)*sum(termIndex.*fourierCoefficients.^2);
    dragCoefficient = dragArea/options.ReferenceArea;
end
