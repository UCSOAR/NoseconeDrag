function air = standardAtmosphere(altitude)
    % standardAtmosphere Air properties of the 1976 standard atmosphere up to 71 km.
    %
    %   air = standardAtmosphere(altitude) returns a struct of air properties at
    %   geometric altitudes above mean sea level, in metres. Each field has the
    %   size of altitude:
    %       Temperature        - K
    %       Pressure           - Pa
    %       Density            - kg/m^3
    %       SpeedOfSound       - m/s
    %       DynamicViscosity   - Pa*s
    %       KinematicViscosity - m^2/s
    %
    %   Method
    %       Geometric altitude z becomes geopotential altitude h = r0*z/(r0 + z)
    %       with r0 = 6356766 m. Temperature is linear in h within each layer and
    %       pressure follows the barometric formula
    %           L ~= 0   P = Pb*(Tb/(Tb + L*(h - hb)))^(g0*M/(R*L))
    %           L == 0   P = Pb*exp(-g0*M*(h - hb)/(R*Tb))
    %       with g0 = 9.80665 m/s^2, R = 8314.32 J/(kmol*K) and M = 28.9644 kg/kmol.
    %       Density is P*M/(R*T) and the speed of sound is sqrt(gamma*R*T/M) with
    %       gamma = 1.4. Dynamic viscosity uses the linear fit from OpenRocket's
    %       AtmosphericConditions.java, mu = 3.7291e-6 + 4.9944e-8*T.

    arguments
        altitude double {mustBeNonnegative, mustBeLessThanOrEqual(altitude, 71000)}
    end

    % layer table: base geopotential altitude, lapse rate, base temperature, base pressure
    baseAltitude = [0, 11000, 20000, 32000, 47000, 51000];                    % m
    lapseRate = [-6.5, 0, 1.0, 2.8, 0, -2.8]/1000;                            % K/m
    baseTemperature = [288.15, 216.65, 216.65, 228.65, 270.65, 270.65];       % K
    basePressure = [101325, 22632.1, 5474.89, 868.019, 110.9063, 66.9389];    % Pa

    earthRadius = 6356766;              % m, for geopotential altitude
    standardGravity = 9.80665;          % m/s^2
    gasConstant = 8314.32;              % J/(kmol*K)
    molarMass = 28.9644;                % kg/kmol
    heatCapacityRatio = 1.4;
    viscosityIntercept = 3.7291e-6;     % Pa*s
    viscositySlope = 4.9944e-8;         % Pa*s/K

    geopotential = earthRadius*altitude./(earthRadius + altitude);
    layer = discretize(geopotential, [baseAltitude, Inf]);

    % temperature and pressure, one layer at a time
    temperature = zeros(size(altitude));
    pressure = zeros(size(altitude));
    hydrostaticConstant = standardGravity*molarMass/gasConstant;              % K/m
    for k = 1:numel(baseAltitude)
        inLayer = layer == k;
        height = geopotential(inLayer) - baseAltitude(k);
        temperature(inLayer) = baseTemperature(k) + lapseRate(k)*height;
        if abs(lapseRate(k)) < eps
            pressure(inLayer) = basePressure(k)*exp(-hydrostaticConstant*height/baseTemperature(k));
        else
            exponent = hydrostaticConstant/lapseRate(k);
            pressure(inLayer) = basePressure(k)*(baseTemperature(k)./temperature(inLayer)).^exponent;
        end
    end

    density = pressure*molarMass./(gasConstant*temperature);
    dynamicViscosity = viscosityIntercept + viscositySlope*temperature;

    air = struct(Temperature=temperature, Pressure=pressure, Density=density, ...
        SpeedOfSound=sqrt(heatCapacityRatio*gasConstant*temperature/molarMass), ...
        DynamicViscosity=dynamicViscosity, KinematicViscosity=dynamicViscosity./density);
end
