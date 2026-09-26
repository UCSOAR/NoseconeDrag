function [dragCoefficient, details] = rocketFrictionDrag(mach, altitude, vehicleLength, bodyDiameter, ...
        bodyWettedArea, options)
    % rocketFrictionDrag Skin friction drag coefficient of a finned rocket (OpenRocket model).
    %
    %   cd = rocketFrictionDrag(mach, altitude, vehicleLength, bodyDiameter, bodyWettedArea)
    %   returns the skin friction drag coefficient of a rocket, referenced to the
    %   body cross-section pi*bodyDiameter^2/4, at the given Mach numbers and
    %   geometric altitudes in the standard atmosphere. It evaluates equation
    %   3.85 of the OpenRocket technical documentation,
    %       CD = Cf*((1 + 1/(2*fB))*Awet,body + (1 + 2*t/cbar)*Awet,fins)/Aref
    %   where fB = vehicleLength/bodyDiameter, Awet,fins counts both sides of
    %   every fin, and Cf comes from skinFrictionCoefficient at the Reynolds
    %   number of the whole vehicle, Re = V*vehicleLength/nu, as in OpenRocket.
    %
    %   cd = rocketFrictionDrag(..., Name=Value) adds fins and sets the finish:
    %       FinCount        - number of fins (default 0)
    %       FinPlanformArea - one-sided planform area of one fin, m^2
    %       FinMeanChord    - mean aerodynamic chord of one fin, m
    %       FinThickness    - fin thickness, m (default 0)
    %       Roughness       - roughness height of the body, m (default 60e-6)
    %       FinRoughness    - roughness height of the fins, m (default 60e-6)
    %   The 60 um default is OpenRocket's default finish, "regular paint".
    %
    %   [cd, details] = rocketFrictionDrag(...) also returns a struct with the
    %   body and fin parts of cd (BodyDrag, FinDrag), the vehicle Reynolds number
    %   (Reynolds) and the friction coefficients of the body and the fins
    %   (BodyFrictionCoefficient, FinFrictionCoefficient).
    %
    %   Inputs
    %       mach           - free-stream Mach number
    %       altitude       - geometric altitude above mean sea level, m
    %       vehicleLength  - nose tip to aft end, m
    %       bodyDiameter   - maximum body diameter, m
    %       bodyWettedArea - wetted area of the nose cone, body tubes and transitions, m^2
    %
    %   mach and altitude must have compatible sizes. A row of Mach numbers and a
    %   column of altitudes gives a grid with one row per altitude.

    arguments
        mach double {mustBeNonnegative}
        altitude double {mustBeNonnegative}
        vehicleLength (1, 1) double {mustBePositive}
        bodyDiameter (1, 1) double {mustBePositive}
        bodyWettedArea (1, 1) double {mustBePositive}
        options.FinCount (1, 1) double {mustBeInteger, mustBeNonnegative} = 0
        options.FinPlanformArea (1, 1) double {mustBeNonnegative} = 0
        options.FinMeanChord (1, 1) double {mustBeNonnegative} = 0
        options.FinThickness (1, 1) double {mustBeNonnegative} = 0
        options.Roughness (1, 1) double {mustBeNonnegative} = 60e-6
        options.FinRoughness (1, 1) double {mustBeNonnegative} = 60e-6
    end

    hasFins = options.FinCount > 0;
    if hasFins && ~(options.FinPlanformArea > 0 && options.FinMeanChord > 0)
        error("rocketFrictionDrag:missingFinGeometry", ...
            "FinPlanformArea and FinMeanChord must be positive when FinCount is %d. " + ...
            "Set both, or use FinCount=0 for a body without fins.", options.FinCount);
    end

    % flight condition
    air = standardAtmosphere(altitude);
    velocity = mach.*air.SpeedOfSound;
    reynolds = velocity*vehicleLength./air.KinematicViscosity;

    % friction coefficients, with roughness relative to the vehicle length as in OpenRocket
    bodyFrictionCoefficient = skinFrictionCoefficient(reynolds, mach, ...
        RelativeRoughness=options.Roughness/vehicleLength);
    finFrictionCoefficient = skinFrictionCoefficient(reynolds, mach, ...
        RelativeRoughness=options.FinRoughness/vehicleLength);

    % equation 3.85: cylindrical body correction and fin thickness correction
    referenceArea = pi*bodyDiameter^2/4;
    finenessRatio = vehicleLength/bodyDiameter;
    bodyDrag = (1 + 1/(2*finenessRatio))*bodyFrictionCoefficient*bodyWettedArea/referenceArea;
    finDrag = zeros(size(bodyDrag));
    if hasFins
        finWettedArea = 2*options.FinCount*options.FinPlanformArea;
        thicknessCorrection = 1 + 2*options.FinThickness/options.FinMeanChord;
        finDrag = thicknessCorrection*finFrictionCoefficient*finWettedArea/referenceArea;
    end

    dragCoefficient = bodyDrag + finDrag;
    details = struct(BodyDrag=bodyDrag, FinDrag=finDrag, Reynolds=reynolds, ...
        BodyFrictionCoefficient=bodyFrictionCoefficient, FinFrictionCoefficient=finFrictionCoefficient);
end
