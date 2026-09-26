%[text] # Von Karman nose cone: drag versus fineness ratio
%[text] Supersonic wave drag and skin friction of the LD-Haack (von Karman) nose profile as a function of fineness ratio $L/d$, and the trade-off that sets the optimum length. Wave drag comes from linearized slender-body theory and is checked numerically from the profile geometry. Skin friction uses the OpenRocket (Barrowman) turbulent flat-plate model. Requires `vonKarmanProfile.m`, `slenderBodyWaveDrag.m` and `skinFrictionCoefficient.m` on the path.

%%
%[text] ## Design point
%[text] Change these values to match the vehicle. Reynolds number and relative roughness use the nose length as reference length, so the friction estimate is for the nose alone.
baseDiameter = 0.15;                        % m, nose base diameter
flightMach = 2.0;                           % design Mach number for the study
speedOfSound = 340;                         % m/s
kinematicViscosity = 1.5e-5;                % m^2/s, air near sea level
surfaceRoughness = 20e-6;                   % m, painted surface (OpenRocket table 3.2)
finenessRatios = linspace(2, 10, 161);      % L/d values to sweep
unitReynolds = flightMach*speedOfSound/kinematicViscosity   % 1/m

%%
%[text] ## Profile
%[text] The Haack series radius is $y = \\frac{R}{\\sqrt{\\pi}}\\sqrt{\\theta - \\frac{\\sin 2\\theta}{2} + C \\sin^3\\theta}$ with $\\theta = \\arccos(1 - 2x/L)$. $C = 0$ gives the von Karman (LD-Haack) profile, the shape of minimum wave drag for a given length and base diameter. Fineness ratio only stretches the same curve along the axis.
profileRatios = [2, 3, 5];
profileColors = ["#86b6ef", "#3987e5", "#184f95"];   % one-hue ordinal ramp, light to dark
hold on
for k = 1:numel(profileRatios)
    [xProfile, yProfile] = vonKarmanProfile(profileRatios(k)*baseDiameter, baseDiameter);
    plot(xProfile/baseDiameter, yProfile/baseDiameter, Color=profileColors(k), LineWidth=1.5, ...
        DisplayName="L/d = " + profileRatios(k))
    text(profileRatios(k), 0.56, "L/d = " + profileRatios(k), HorizontalAlignment="center", ...
        Color="#0b0b0b", FontSize=9)
end
hold off
axis equal
xlim([0, 5.6])
ylim([0, 0.7])
box off
xlabel("x / d")
ylabel("y / d")
title("Von Karman profile normalized by base diameter")
legend(Location="southeast", Box="off")

%%
%[text] ## Wave drag from slender-body theory
%[text] Karman-Moore linearized theory gives the wave drag of a slender body of revolution from its area distribution $S(x)$ as $D = -\\frac{\\rho U^2}{4\\pi}\\int_0^L\\int_0^L S''(x_1) S''(x_2) \\ln|x_1 - x_2| \\, dx_1 dx_2$. Writing $S'(x) = \\sum_n A_n \\sin n\\theta$ turns this into $D/q = \\frac{\\pi}{4}\\sum_n n A_n^2$, and the base-area constraint fixes $A_1 = 4 S_B/(\\pi L)$. The von Karman profile keeps only $A_1$, so its wave drag coefficient referenced to base area is $C_{D,w} = \\frac{4 S_B}{\\pi L^2} = (d/L)^2$, independent of Mach number within the theory. The check below evaluates the Fourier coefficients numerically from the profile coordinates alone.
checkRatio = 3;
[xCheck, yCheck] = vonKarmanProfile(checkRatio*baseDiameter, baseDiameter, NumPoints=2001);
[waveDragNumeric, fourierCoefficients] = slenderBodyWaveDrag(xCheck, yCheck);
higherHarmonics = max(abs(fourierCoefficients(2:end)))/fourierCoefficients(1);
waveDragCheck = table((1/checkRatio)^2, waveDragNumeric, higherHarmonics, ...
    VariableNames=["Analytic", "Numeric", "HigherHarmonicsRelativeToA1"])

%%
%[text] The wave drag falls with the inverse square of the fineness ratio: doubling $L/d$ cuts it by a factor of four.
waveDrag = 1./finenessRatios.^2;

%%
%[text] ## Skin friction
%[text] Friction drag referenced to base area is $C_{D,f} = C_f S_{wet}/S_B$. The wetted area is integrated from the profile, and the average friction coefficient $C_f$ uses the OpenRocket turbulent flat-plate model with its supersonic compressibility correction and the roughness-limited branch. The wetted area grows almost linearly with $L/d$, so friction rises linearly while wave drag falls quadratically.
noseLengths = finenessRatios*baseDiameter;
reynoldsNumbers = unitReynolds*noseLengths;
wettedAreaRatio = zeros(size(finenessRatios));      % S_wet / S_B
frictionCoefficients = zeros(size(finenessRatios));
baseArea = pi*baseDiameter^2/4;
for k = 1:numel(finenessRatios)
    [xProfile, yProfile] = vonKarmanProfile(noseLengths(k), baseDiameter, NumPoints=1001);
    arcLength = [0, cumsum(hypot(diff(xProfile), diff(yProfile)))];
    wettedAreaRatio(k) = trapz(arcLength, 2*pi*yProfile)/baseArea;
    frictionCoefficients(k) = skinFrictionCoefficient(reynoldsNumbers(k), flightMach, ...
        RelativeRoughness=surfaceRoughness/noseLengths(k));
end
frictionDrag = frictionCoefficients.*wettedAreaRatio;
wettedAreaPerFineness = mean(wettedAreaRatio./finenessRatios)   % S_wet/S_B per unit L/d

%%
%[text] ## Total drag and optimum fineness ratio
%[text] The total is the sum of the two terms. The minimum sits where the friction penalty of extra length balances the wave drag saving.
totalDrag = waveDrag + frictionDrag;
[minimumDrag, minimumIndex] = min(totalDrag);
optimumFineness = finenessRatios(minimumIndex)

%%
seriesColors = ["#2a78d6", "#eb6834", "#1baf7a"];   % total, wave, friction
plot(finenessRatios, totalDrag, Color=seriesColors(1), LineWidth=1.5, DisplayName="Total")
hold on
plot(finenessRatios, waveDrag, Color=seriesColors(2), LineWidth=1.5, DisplayName="Wave drag (d/L)^2")
plot(finenessRatios, frictionDrag, Color=seriesColors(3), LineWidth=1.5, DisplayName="Skin friction")
plot(optimumFineness, minimumDrag, "o", MarkerSize=8, MarkerFaceColor=seriesColors(1), ...
    MarkerEdgeColor="#fcfcfb", HandleVisibility="off")
text(optimumFineness, minimumDrag + 0.02, "minimum at L/d = " + string(round(optimumFineness, 1)), ...
    HorizontalAlignment="center", Color="#0b0b0b", FontSize=9)
hold off
box off
grid on
xlabel("Fineness ratio L/d")
ylabel("C_D referenced to base area")
title("Von Karman nose drag at Mach " + flightMach)
legend(Location="northeast", Box="off")

%%
%[text] Values at integer fineness ratios. Below $L/d \\approx 4$ the wave drag dominates and each added caliber of length pays for itself; above the optimum the curve is flat, so the choice is driven by mass and structure rather than drag.
tableRatios = (2:10).';
[~, tableIndex] = min(abs(finenessRatios - tableRatios), [], 2);
dragTable = table(tableRatios, waveDrag(tableIndex).', frictionDrag(tableIndex).', totalDrag(tableIndex).', ...
    VariableNames=["FinenessRatio", "WaveDrag", "FrictionDrag", "TotalDrag"])

%%
%[text] ## Comparison with the OpenRocket data-based model
%[text] OpenRocket estimates nose pressure drag from Stoney's free-flight data at fineness ratio 3 and extrapolates in fineness ratio with $C_D = C_0 \\, (C_3/C_0)^{\\log_4(f_N + 1)}$, where $C_0 = 0.85 \\, q_{stag}/q$ is the blunt-cylinder value. The von Karman data curve below is reproduced from the OpenRocket source. The theory line is the same at every Mach number; the data-based curves fall slowly with Mach and sit below the theory near $L/d = 3$.
vonKarmanMach = [0.9, 0.95, 1.0, 1.05, 1.1, 1.2, 1.4, 1.6, 2.0, 3.0];
vonKarmanDragAtThree = [0, 0.010, 0.027, 0.055, 0.070, 0.081, 0.095, 0.097, 0.091, 0.083];
comparisonMach = [1.5, 2.0, 2.5];
stagnationRatio = 1.84 - 0.76./comparisonMach.^2 + 0.166./comparisonMach.^4 + 0.035./comparisonMach.^6;
bluntDrag = 0.85*stagnationRatio;
dragAtThree = interp1(vonKarmanMach, vonKarmanDragAtThree, comparisonMach);
finenessExponent = log(finenessRatios + 1)/log(4);
empiricalWaveDrag = bluntDrag.'.*(dragAtThree.'./bluntDrag.').^finenessExponent;   % rows: Mach

machColors = ["#ef9668", "#d1521f", "#8a3410"];     % one-hue ordinal ramp, light to dark
plot(finenessRatios, waveDrag, Color="#2a78d6", LineWidth=1.5, DisplayName="Slender-body theory, any Mach")
hold on
for k = 1:numel(comparisonMach)
    plot(finenessRatios, empiricalWaveDrag(k, :), Color=machColors(k), LineWidth=1.5, ...
        DisplayName="OpenRocket model, Mach " + comparisonMach(k))
end
hold off
box off
grid on
xlabel("Fineness ratio L/d")
ylabel("Wave drag C_D referenced to base area")
title("Von Karman wave drag: theory versus data-based model")
legend(Location="northeast", Box="off")

%%
%[text] Wave drag at integer fineness ratios: theory in the second column, the data-based model at Mach 1.5, 2.0 and 2.5 in the last three.
comparisonTable = table(tableRatios, waveDrag(tableIndex).', empiricalWaveDrag(:, tableIndex).', ...
    VariableNames=["FinenessRatio", "Theory", "OpenRocketModel"])

%%
%[text] Optimum fineness ratio with each wave drag model at the design Mach number, using the same friction estimate.
designStagnation = 1.84 - 0.76/flightMach^2 + 0.166/flightMach^4 + 0.035/flightMach^6;
designBluntDrag = 0.85*designStagnation;
designDragAtThree = interp1(vonKarmanMach, vonKarmanDragAtThree, flightMach);
designEmpiricalWaveDrag = designBluntDrag*(designDragAtThree/designBluntDrag).^finenessExponent;
totalDragEmpirical = designEmpiricalWaveDrag + frictionDrag;
[minimumDragEmpirical, minimumIndexEmpirical] = min(totalDragEmpirical);
optimumTable = table(["Slender-body theory"; "OpenRocket data model"], ...
    [optimumFineness; finenessRatios(minimumIndexEmpirical)], [minimumDrag; minimumDragEmpirical], ...
    VariableNames=["WaveDragModel", "OptimumFineness", "MinimumTotalDrag"])

%%
%[text] ## Assumptions and limits
%[text] - Nose alone. Body tube, fins and base drag are not included. Extra nose length on a fixed body also adds mass and moves the centre of gravity forward, which this study does not price. \
%[text] - Slender-body theory needs a slender nose and Mach numbers away from 1. It has no transonic drag rise and no Mach dependence, and it over-predicts the free-flight data near $L/d = 3$ by roughly 15 to 25 percent at Mach 1.6 to 3. \
%[text] - Friction assumes a fully turbulent boundary layer from the tip, with the nose length as reference length for Reynolds number and relative roughness. Surface finish matters: with the 20 um roughness above the roughness-limited branch is active over most of the sweep. \
%[text] - If total vehicle length is fixed, lengthening the nose shortens the body tube and the friction penalty nearly vanishes, so the drag optimum moves to higher $L/d$ than shown here. \

%%
%[text] ## References
%[text] - Haack series profile equation: [Wikipedia, Nose cone design](https://en.wikipedia.org/wiki/Nose_cone_design) \
%[text] - Karman-Moore wave drag integral: [Wikipedia, Karman-Moore theory](https://en.wikipedia.org/wiki/K%C3%A1rm%C3%A1n%E2%80%93Moore_theory) and [Wikipedia, Sears-Haack body](https://en.wikipedia.org/wiki/Sears%E2%80%93Haack_body) \
%[text] - Von Karman ogive result $C_D (l/d)^2 = 1$ and Sears-Haack $9\\pi^2/8$: [W. H. Mason, Curiosity 16](https://archive.aoe.vt.edu/mason/Mason_f/C16K-OandS-H.pdf) \
%[text] - Skin friction, compressibility corrections, roughness and the nose cone drag model: [OpenRocket technical documentation v13.05, sections 3.4.2 and 3.4.3, appendix B](https://openrocket.sourceforge.net/techdoc.pdf) \
%[text] - Von Karman fineness-ratio-3 data curve: [OpenRocket source, SymmetricComponentCalc.java](https://raw.githubusercontent.com/openrocket/openrocket/unstable/core/src/main/java/info/openrocket/core/aerodynamics/barrowman/SymmetricComponentCalc.java) \
%[text] - Free-flight data compilation behind the OpenRocket curves: [Stoney, NASA TR R-100](https://ntrs.nasa.gov/citations/19630004995) \

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline"}
%---
