function frictionCoefficient = skinFrictionCoefficient(reynolds, mach, options)
    % skinFrictionCoefficient Average turbulent skin friction coefficient (OpenRocket model).
    %
    %   cf = skinFrictionCoefficient(reynolds, mach) returns the average skin
    %   friction coefficient, referenced to wetted area, of a fully turbulent
    %   boundary layer at the given length-based Reynolds number and free-stream
    %   Mach number. The formulas are those of the OpenRocket technical
    %   documentation, section 3.4.2 (equations 3.78 to 3.84), applied the way
    %   OpenRocket's BarrowmanDragCalculator.java applies them:
    %       smooth turbulent   cf = 1/(1.50*ln(Re) - 5.6)^2, set to 1.48e-2 below Re = 1e4
    %       roughness limited  cf = 0.032*(Rs/L)^0.2
    %   with the compressibility corrections
    %       M < 0.9            cf*(1 - 0.1*M^2)
    %       M > 1.1, turbulent cf/(1 + 0.15*M^2)^0.58
    %       M > 1.1, rough     cf/(1 + 0.18*M^2)
    %   interpolated linearly in Mach between 0.9 and 1.1. The result is the
    %   larger of the corrected turbulent and roughness-limited values. The
    %   technical documentation instead switches to the roughness-limited value
    %   above Re = 51*(Rs/L)^-1.039, which gives a slightly lower value for
    %   subsonic Reynolds numbers just below that threshold.
    %
    %   cf = skinFrictionCoefficient(..., RelativeRoughness=r) sets the surface
    %   roughness height divided by the reference length, Rs/L (default 0).
    %
    %   reynolds and mach must have compatible sizes.

    arguments
        reynolds double {mustBeNonnegative}
        mach double {mustBeNonnegative}
        options.RelativeRoughness (1, 1) double {mustBeNonnegative} = 0
    end

    % expand both inputs to a common size
    reynolds = reynolds.*ones(size(mach));
    mach = mach.*ones(size(reynolds));

    % smooth turbulent flat plate, constant below Re = 1e4
    minimumReynolds = 1e4;
    lowReynoldsValue = 1.48e-2;
    turbulent = 1./(1.50*log(max(reynolds, minimumReynolds)) - 5.6).^2;
    turbulent(reynolds < minimumReynolds) = lowReynoldsValue;

    % weight of the supersonic correction: 0 below Mach 0.9, 1 above Mach 1.1
    subsonicLimit = 0.9;
    supersonicLimit = 1.1;
    supersonicWeight = min(max((mach - subsonicLimit)/(supersonicLimit - subsonicLimit), 0), 1);

    % the turbulent correction blends both formulas evaluated at the local Mach number
    subsonicFactor = 1 - 0.1*mach.^2;
    turbulentSupersonic = 1./(1 + 0.15*mach.^2).^0.58;
    turbulentFactor = (1 - supersonicWeight).*subsonicFactor + supersonicWeight.*turbulentSupersonic;

    % the roughness correction blends the values at the band edges, as OpenRocket does
    roughSubsonic = 1 - 0.1*min(mach, subsonicLimit).^2;
    roughSupersonic = 1./(1 + 0.18*max(mach, supersonicLimit).^2);
    roughFactor = (1 - supersonicWeight).*roughSubsonic + supersonicWeight.*roughSupersonic;

    % the roughness-limited value applies wherever it exceeds the turbulent value
    roughnessLimited = 0.032*options.RelativeRoughness^0.2*roughFactor;
    frictionCoefficient = max(turbulent.*turbulentFactor, roughnessLimited);
end
