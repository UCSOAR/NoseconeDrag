function [totalCd, pressureCd, frictionCd, speedOfSound] = vonKarmanDragCoefficient(mach, finenessRatio, options)
% vonKarmanDragCoefficient Zero-lift drag coefficient of a von Karman (LD-Haack) nose cone.
%
%   totalCd = vonKarmanDragCoefficient(mach, finenessRatio) returns the drag
%   coefficient of a von Karman nose cone mounted on a cylindrical body, referenced
%   to the base area pi*D^2/4. mach is the free-stream Mach number and finenessRatio
%   is L/D. totalCd has one row per fineness ratio and one column per Mach number,
%   so surf(mach, finenessRatio, totalCd) plots it directly.
%
%   [totalCd, pressureCd, frictionCd, speedOfSound] = vonKarmanDragCoefficient(...)
%   also returns the pressure (wave) and skin friction parts and the speed of
%   sound in m/s, so velocity = mach*speedOfSound.
%
%   name-value arguments
%       BaseDiameter     - nose base diameter in m (default 0.1)
%       SurfaceRoughness - surface roughness height in m (default 60e-6, regular paint)
%       Temperature      - static air temperature in K (default 288.15, sea level)
%       Pressure         - static air pressure in Pa (default 101325, sea level)
%
%   method (semi-empirical component build-up used by OpenRocket)
%       pressure drag: free-flight data for a von Karman nose with L/D = 3 between
%       Mach 0.9 and 3 (NASA TR-R-100), rescaled to other fineness ratios f with
%           cd = cdStagnation*(cd3/cdStagnation)^(log(f + 1)/log(4))
%       where cdStagnation is the drag of a flat-faced cylinder (f = 0). the
%       pressure drag is zero below Mach 0.9 and held constant above Mach 3.
%       friction drag: fully turbulent flat-plate skin friction with compressibility
%       and roughness limits, times the wetted-to-base area ratio and the body
%       correction 1 + 1/(2*L/D).
%       base drag is excluded because the nose is assumed to sit on a body tube.
%
%   sources
%       https://github.com/openrocket/openrocket
%       core/src/main/java/info/openrocket/core/aerodynamics/barrowman/SymmetricComponentCalc.java
%       core/src/main/java/info/openrocket/core/aerodynamics/BarrowmanCalculator.java
%
%   see also plotVonKarmanDragSurface

    arguments
        mach (1, :) double {mustBeNonempty, mustBeFinite, mustBePositive}
        finenessRatio (:, 1) double {mustBeNonempty, mustBeFinite, mustBePositive}
        options.BaseDiameter (1, 1) double {mustBeFinite, mustBePositive} = 0.1
        options.SurfaceRoughness (1, 1) double {mustBeFinite, mustBeNonnegative} = 60e-6
        options.Temperature (1, 1) double {mustBeFinite, mustBePositive} = 288.15
        options.Pressure (1, 1) double {mustBeFinite, mustBePositive} = 101325
    end

    % dry-air properties from the ideal gas law and sutherland's viscosity law
    gasConstant = 287.053;              % J/(kg*K)
    heatCapacityRatio = 1.4;
    sutherlandConstant = 1.458e-6;      % kg/(m*s*K^0.5)
    sutherlandTemperature = 110.4;      % K
    temperature = options.Temperature;

    density = options.Pressure/(gasConstant*temperature);
    speedOfSound = sqrt(heatCapacityRatio*gasConstant*temperature);
    dynamicViscosity = sutherlandConstant*temperature^1.5/(temperature + sutherlandTemperature);
    kinematicViscosity = dynamicViscosity/density;

    % flight conditions; rows follow finenessRatio, columns follow mach
    velocity = mach*speedOfSound;
    noseLength = finenessRatio*options.BaseDiameter;
    reynolds = noseLength.*velocity/kinematicViscosity;

    pressureCd = nosePressureDrag(mach, finenessRatio);

    % body correction evaluated with the fineness ratio L/D
    bodyCorrection = 1 + 1./(2*finenessRatio);
    skinFriction = turbulentSkinFriction(mach, reynolds, noseLength, options.SurfaceRoughness);
    frictionCd = skinFriction.*(wettedAreaRatio(finenessRatio).*bodyCorrection);

    totalCd = pressureCd + frictionCd;
end

function pressureCd = nosePressureDrag(mach, finenessRatio)
% nosePressureDrag Pressure drag of the von Karman nose, referenced to the base area.

    % free-flight data at L/D = 3 from NASA TR-R-100, as tabulated in OpenRocket
    dataMach = [0.9, 0.95, 1.0, 1.05, 1.1, 1.2, 1.4, 1.6, 2.0, 3.0];
    dataCd = [0, 0.010, 0.027, 0.055, 0.070, 0.081, 0.095, 0.097, 0.091, 0.083];
    dataFinenessRatio = 3;

    % geometric blend between the flat-faced cylinder (L/D = 0) and the L/D = 3 data
    stagnationCd = stagnationDragCoefficient(dataMach);
    blendExponent = log(finenessRatio + 1)/log(dataFinenessRatio + 1);
    scaledCd = stagnationCd.*(dataCd./stagnationCd).^blendExponent;

    % linear in mach between data points, clamped to the end values outside them
    clampedMach = min(max(mach(:), dataMach(1)), dataMach(end));
    pressureCd = interp1(dataMach, scaledCd.', clampedMach).';
end

function stagnationCd = stagnationDragCoefficient(mach)
% stagnationDragCoefficient Pressure drag of a flat-faced cylinder from its stagnation pressure.

    stagnationFactor = 0.85;

    % stagnation-to-free-stream pressure ratio, subsonic and supersonic fits
    pressureRatio = 1 + mach.^2/4 + mach.^4/40;
    isSupersonic = mach > 1;
    supersonicMach = mach(isSupersonic);
    pressureRatio(isSupersonic) = 1.84 - 0.76./supersonicMach.^2 + 0.166./supersonicMach.^4 ...
        + 0.035./supersonicMach.^6;

    stagnationCd = stagnationFactor*pressureRatio;
end

function frictionCoefficient = turbulentSkinFriction(mach, reynolds, noseLength, roughness)
% turbulentSkinFriction Turbulent skin friction coefficient with compressibility and roughness limits.

    lowReynolds = 1.0e4;
    lowReynoldsCf = 1.48e-2;
    roughnessCoefficient = 0.032;
    subsonicCorrection = @(m) 1 - 0.1*m.^2;

    % transonic blending weight: 0 at mach 0.9 and below, 1 at mach 1.1 and above
    lowerMach = 0.9;
    upperMach = 1.1;
    weight = min(max((mach - lowerMach)/(upperMach - lowerMach), 0), 1);

    % smooth turbulent flat plate, constant below re = 1e4
    smoothCf = 1./(1.50*log(reynolds) - 5.6).^2;
    smoothCf(reynolds < lowReynolds) = lowReynoldsCf;
    compressibility = (1 - weight).*subsonicCorrection(mach) + weight./(1 + 0.15*mach.^2).^0.58;
    smoothCf = smoothCf.*compressibility;

    % roughness-limited value; its transonic blend uses the corrections at mach 0.9 and 1.1
    roughnessCorrection = (1 - weight).*subsonicCorrection(min(mach, lowerMach)) ...
        + weight./(1 + 0.18*max(mach, upperMach).^2);
    roughnessScale = (roughness./noseLength).^0.2;
    roughCf = roughnessCoefficient*roughnessScale.*roughnessCorrection;

    % the larger of the smooth and roughness-limited values governs
    frictionCoefficient = max(smoothCf, roughCf);
end

function areaRatio = wettedAreaRatio(finenessRatio)
% wettedAreaRatio Wetted area of the von Karman nose divided by its base area.
%   with x = (L/2)*(1 - cos(theta)) and y = (R/sqrt(pi))*sqrt(g), g = theta - sin(2*theta)/2,
%   the surface of revolution gives
%   S_wet/S_base = 2*integral(sin(theta)*sqrt(f^2*g/pi + sin(theta)^2/pi^2), 0, pi)
%   for f = L/D. the integrand is smooth at both ends, so trapz converges quickly.

    numPoints = 2001;
    theta = linspace(0, pi, numPoints);
    haackTerm = max(theta - sin(2*theta)/2, 0);
    integrand = sin(theta).*sqrt(finenessRatio.^2.*haackTerm/pi + sin(theta).^2/pi^2);
    areaRatio = 2*trapz(theta, integrand, 2);
end
