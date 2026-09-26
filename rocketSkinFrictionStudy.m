%[text] # Skin friction drag of Athena
%[text] Skin friction drag coefficient of a finned rocket over Mach number and altitude, using the method in OpenRocket. A turbulent flat-plate friction coefficient at the Reynolds number of the whole vehicle is corrected for compressibility, raised to the roughness-limited value where the surface finish controls, then scaled from wetted area to the body cross-section. Requires `skinFrictionCoefficient.m`, `standardAtmosphere.m`, `rocketFrictionDrag.m` and `vonKarmanProfile.m` on the path.

%%
%[text] ## Vehicle
%[text] The Athena technical baseline puts the diameter close to 10 in, allows three or four fins and calls for an automotive-grade painted finish, which the 20 µm OpenRocket smooth paint preset stands in for. Fin planform and thickness come from the Aegis/Athena fin drawing in the sandwich panel fins CDR. The baseline leaves the overall length to the MDO, so the vehicle and nose lengths below are placeholders. Replace them with the OpenRocket model's values before quoting any result.
inchToMeter = 0.0254;
bodyDiameter = 10*inchToMeter;              % m, technical baseline
vehicleLength = 6.0;                        % m, nose tip to aft end, placeholder
noseLength = 5*bodyDiameter;                % m, von Karman nose, placeholder
finCount = 3;                               % baseline allows 3 or 4
finRootChord = 17.717*inchToMeter;          % m, outline including the edge inserts
finTipChord = 10.402*inchToMeter;           % m
finSpan = 5.512*inchToMeter;                % m
finThickness = 0.39*inchToMeter;            % m, average after the tip-to-tip layup
surfaceRoughness = 20e-6;                   % m, OpenRocket smooth paint preset
fineness = vehicleLength/bodyDiameter

%%
%[text] ## Wetted areas
%[text] The body wetted area is the nose cone surface, integrated from the von Karman profile, plus the cylinder behind it. Each trapezoidal fin has a one-sided planform area $A_{fin} = s(C_r + C_t)/2$ and a mean aerodynamic chord $\\bar{c} = \\frac{1}{A_{fin}}\\int_0^s c^2(y)\\,dy = \\frac{2}{3}\\frac{C_r^2 + C_r C_t + C_t^2}{C_r + C_t}$. The fin wetted area counts both sides of every fin. OpenRocket references drag coefficients to the body cross-section $A_{ref} = \\pi d^2/4$, and so does this script. Areas are in m².
[xNose, yNose] = vonKarmanProfile(noseLength, bodyDiameter, NumPoints=1001);
noseArcLength = [0, cumsum(hypot(diff(xNose), diff(yNose)))];
noseWettedArea = trapz(noseArcLength, 2*pi*yNose);
tubeWettedArea = pi*bodyDiameter*(vehicleLength - noseLength);
bodyWettedArea = noseWettedArea + tubeWettedArea;
finPlanformArea = finSpan*(finRootChord + finTipChord)/2;
finMeanChord = (2/3)*(finRootChord^2 + finRootChord*finTipChord + finTipChord^2)/(finRootChord + finTipChord);
finWettedArea = 2*finCount*finPlanformArea;
referenceArea = pi*bodyDiameter^2/4;
wettedAreas = [noseWettedArea; tubeWettedArea; finWettedArea];
areaTable = table(["Nose cone"; "Body tube"; "Fins"], wettedAreas, wettedAreas/referenceArea, ...
    VariableNames=["Surface", "WettedArea", "RatioToReferenceArea"])

%%
%[text] ## Worked example at one flight condition
%[text] Mach 2.5 at 6 km. The altitude is a placeholder for where Athena reaches peak Mach. OpenRocket bases the Reynolds number on the length of the whole vehicle, $Re = VL/\\nu$, for the body and the fins alike.
designMach = 2.5;
designAltitude = 6000;                      % m, placeholder
air = standardAtmosphere(designAltitude);
velocity = designMach*air.SpeedOfSound;
reynolds = velocity*vehicleLength/air.KinematicViscosity

%%
%[text] There are two candidate friction coefficients. The smooth turbulent flat plate gives $C_f = 1/(1.50 \\ln Re - 5.6)^2$. A rough wall gives $C_f = 0.032 (R_s/L)^{0.2}$, independent of Reynolds number. Above Mach 1.1 compressibility divides the first by $(1 + 0.15M^2)^{0.58}$ and the second by $1 + 0.18M^2$. OpenRocket uses whichever corrected value is larger.
turbulentCoefficient = 1/(1.50*log(reynolds) - 5.6)^2/(1 + 0.15*designMach^2)^0.58;
roughCoefficient = 0.032*(surfaceRoughness/vehicleLength)^0.2/(1 + 0.18*designMach^2);
frictionCoefficient = max(turbulentCoefficient, roughCoefficient);
coefficientTable = table(turbulentCoefficient, roughCoefficient, frictionCoefficient, ...
    VariableNames=["Turbulent", "RoughnessLimited", "Used"])

%%
%[text] Equation 3.85 of the OpenRocket technical documentation converts $C_f$ to a drag coefficient on the reference area, $C_{D,f} = C_f \\left[\\left(1 + \\frac{1}{2f_B}\\right) A_{wet,body} + \\left(1 + \\frac{2t}{\\bar{c}}\\right) A_{wet,fins}\\right] / A_{ref}$. The term with the fineness ratio $f_B = L/d$ corrects the body for its cylindrical shape, and the term with $t/\\bar{c}$ corrects the fins for their thickness. The drag force is $D = \\frac{1}{2}\\rho V^2 C_{D,f} A_{ref}$.
bodyDragExample = frictionCoefficient*(1 + 1/(2*fineness))*bodyWettedArea/referenceArea;
finDragExample = frictionCoefficient*(1 + 2*finThickness/finMeanChord)*finWettedArea/referenceArea;
dragCoefficients = [bodyDragExample; finDragExample; bodyDragExample + finDragExample];
dynamicPressure = 0.5*air.Density*velocity^2;
exampleTable = table(["Body"; "Fins"; "Total"], dragCoefficients, dynamicPressure*referenceArea*dragCoefficients, ...
    VariableNames=["Part", "FrictionDragCoefficient", "DragForceNewtons"])

%%
%[text] The same condition through `rocketFrictionDrag`, which the rest of the script uses. The relative difference from the step-by-step result should be at round-off level.
frictionDrag = @(mach, altitude, roughness) rocketFrictionDrag(mach, altitude, vehicleLength, bodyDiameter, ...
    bodyWettedArea, FinCount=finCount, FinPlanformArea=finPlanformArea, FinMeanChord=finMeanChord, ...
    FinThickness=finThickness, Roughness=roughness, FinRoughness=roughness);
[designDrag, designDetails] = frictionDrag(designMach, designAltitude, surfaceRoughness);
relativeDifference = designDrag/(bodyDragExample + finDragExample) - 1

%%
%[text] ## Friction drag over the flight envelope
%[text] Friction drag coefficient against Mach number at four altitudes with the 20 µm finish. The coefficient falls with Mach number because of the compressibility corrections. Near the ground the finish sets the friction coefficient over almost the whole Mach range. At 10 km it does so between about Mach 0.4 and 1.6, where the 0 and 10 km curves lie on top of each other. At 20 and 30 km the Reynolds number is low enough that the smooth turbulent value governs at every Mach number, and the coefficient is higher.
machNumbers = linspace(0.1, 2.5, 241);
altitudes = [0; 10; 20; 30]*1000;                                   % m
envelopeDrag = frictionDrag(machNumbers, altitudes, surfaceRoughness);  % one row per altitude
altitudeColors = ["#86b6ef", "#3987e5", "#1c5cab", "#0d366b"];      % one-hue ordinal ramp, light to dark
hold on
for k = 1:numel(altitudes)
    plot(machNumbers, envelopeDrag(k, :), Color=altitudeColors(k), LineWidth=1.5, ...
        DisplayName=string(altitudes(k)/1000) + " km")
end
hold off
box off
grid on
xlim([0, 2.5])
ylim([0, 0.6])
xlabel("Mach number")
ylabel("C_{D,f} referenced to body cross-section")
title("Athena skin friction drag coefficient, 20 \mum finish")
legend(Location="northeast", Box="off")

%%
%[text] Values at selected Mach numbers, one column per altitude.
tableMach = [0.5; 1.0; 1.5; 2.0; 2.5];
envelopeTable = array2table(interp1(machNumbers, envelopeDrag.', tableMach), ...
    VariableNames=string(altitudes.'/1000) + " km");
envelopeTable = addvars(envelopeTable, tableMach, Before=1, NewVariableNames="Mach")

%%
%[text] ## Surface finish
%[text] Friction drag coefficient against roughness height at the worked-example condition, with the OpenRocket finish presets marked. On the flat part the smooth turbulent value governs and a smoother finish gains nothing. Past the threshold the roughness-limited value takes over and the coefficient grows as $R_s^{0.2}$.
roughnessHeights = logspace(log10(0.5e-6), log10(500e-6), 200);     % m
finishDrag = zeros(size(roughnessHeights));
for k = 1:numel(roughnessHeights)
    finishDrag(k) = frictionDrag(designMach, designAltitude, roughnessHeights(k));
end
presetNames = ["Polished"; "Optimum paint"; "Smooth paint"; "Regular paint"; "Unfinished"; "Rough"];
presetHeights = [2; 5; 20; 60; 150; 500]*1e-6;                      % m, OpenRocket finish presets
presetDrag = arrayfun(@(height) frictionDrag(designMach, designAltitude, height), presetHeights);
semilogx(roughnessHeights*1e6, finishDrag, Color="#2a78d6", LineWidth=1.5)
hold on
plot(presetHeights*1e6, presetDrag, "o", MarkerSize=8, MarkerFaceColor="#2a78d6", MarkerEdgeColor="#fcfcfb")

% labels sit below and right of their markers, clear of the rising curve; the last one sits above and left
labelShift = 1.2;                                                   % factor along the log axis
labelGap = 0.004;
text(presetHeights(1:end-1)*1e6*labelShift, presetDrag(1:end-1) - labelGap, presetNames(1:end-1), ...
    HorizontalAlignment="left", VerticalAlignment="top", Color="#0b0b0b", FontSize=9)
text(presetHeights(end)*1e6/labelShift, presetDrag(end) + labelGap, presetNames(end), ...
    HorizontalAlignment="right", VerticalAlignment="bottom", Color="#0b0b0b", FontSize=9)
hold off
box off
grid on
xlim([0.5, 500])
ylim([0, 0.25])
xlabel("Roughness height (\mum)")
ylabel("C_{D,f} referenced to body cross-section")
title("Friction drag against surface finish, Mach 2.5 at 6 km")

%%
%[text] Friction drag coefficient at each preset, with roughness in µm and the change from smooth paint in percent. OpenRocket's default finish is regular paint.
finishTable = table(presetNames, presetHeights*1e6, presetDrag, 100*(presetDrag/designDrag - 1), ...
    VariableNames=["Finish", "Roughness", "CDf", "ChangePercent"])

%%
%[text] The threshold is the roughness at which the two candidate values are equal. Below it, sanding and polishing no longer reduce friction in this model. The threshold moves with the flight condition. Below Mach 0.9 both candidate values carry the same compressibility factor, so it cancels. Thresholds are in µm.
thresholdRoughness = vehicleLength*(turbulentCoefficient*(1 + 0.18*designMach^2)/0.032)^5;
seaLevel = standardAtmosphere(0);
lowSpeedMach = 0.5;
lowSpeedReynolds = lowSpeedMach*seaLevel.SpeedOfSound*vehicleLength/seaLevel.KinematicViscosity;
lowSpeedTurbulent = 1/(1.50*log(lowSpeedReynolds) - 5.6)^2;
lowSpeedThreshold = vehicleLength*(lowSpeedTurbulent/0.032)^5;
thresholdTable = table(["Mach 2.5 at 6 km"; "Mach 0.5 at sea level"], 1e6*[thresholdRoughness; lowSpeedThreshold], ...
    VariableNames=["Condition", "ThresholdRoughness"])

%%
%[text] ## Sensitivity to the fin Reynolds length
%[text] OpenRocket evaluates the fins at the Reynolds number of the whole vehicle. The boundary layer on a fin starts at its own leading edge, so basing the fin Reynolds number and relative roughness on the fin mean chord is a reasonable alternative. That choice is engineering judgement, not part of the OpenRocket method. At the worked-example condition it raises the fin friction coefficient by about 65 percent and the total friction drag by about 4 percent.
finReynolds = velocity*finMeanChord/air.KinematicViscosity;
finCoefficientChord = skinFrictionCoefficient(finReynolds, designMach, RelativeRoughness=surfaceRoughness/finMeanChord);
finDragChord = finCoefficientChord*(1 + 2*finThickness/finMeanChord)*finWettedArea/referenceArea;
lengthTable = table(["Vehicle length (OpenRocket)"; "Fin mean chord"], [reynolds; finReynolds], ...
    [designDetails.FinFrictionCoefficient; finCoefficientChord], [designDetails.FinDrag; finDragChord], ...
    designDetails.BodyDrag + [designDetails.FinDrag; finDragChord], ...
    VariableNames=["Length", "FinRe", "FinCf", "FinCDf", "TotalCDf"])

%%
%[text] ## Friction drag along a simulated flight
%[text] `rocketFrictionDrag` accepts arrays, so a trajectory works directly. Export Mach number and altitude above sea level from an OpenRocket simulation and pass the two columns to `frictionDrag` with the chosen roughness. OpenRocket can also export its own friction drag coefficient, Cdf, for comparison. Expect small differences, because OpenRocket computes the speed of sound from a linear fit in temperature and integrates the modelled nose, transition and fin shapes.

%%
%[text] ## Assumptions and limits
%[text] - The boundary layer is turbulent from the nose tip, as in OpenRocket.
%[text] - Every component uses the Reynolds number of the whole vehicle. The fin sensitivity section shows the effect of that choice.
%[text] - Body and fins have the same roughness. Rail buttons, fin fillets and the edge of the tip-to-tip layup are not modelled.
%[text] - Vehicle length, nose length, fin count and the altitude of peak Mach are placeholders. The body term scales with wetted area, so the length placeholder matters most.
%[text] - Only skin friction is included. Wave drag, base drag and fin pressure drag are separate terms.
%[text] - Air properties follow the standard atmosphere with no launch-site temperature offset. \

%%
%[text] ## References
%[text] - Skin friction method, equations 3.12, 3.30 and 3.76 to 3.85, and roughness table 3.2: [OpenRocket technical documentation v13.05](https://openrocket.sourceforge.net/techdoc.pdf)
%[text] - Transonic interpolation, the larger-value rule, and body and fin terms: [BarrowmanDragCalculator.java](https://raw.githubusercontent.com/openrocket/openrocket/unstable/core/src/main/java/info/openrocket/core/aerodynamics/BarrowmanDragCalculator.java), [FinSetCalc.java](https://raw.githubusercontent.com/openrocket/openrocket/unstable/core/src/main/java/info/openrocket/core/aerodynamics/barrowman/FinSetCalc.java), [SymmetricComponentCalc.java](https://raw.githubusercontent.com/openrocket/openrocket/unstable/core/src/main/java/info/openrocket/core/aerodynamics/barrowman/SymmetricComponentCalc.java)
%[text] - Finish presets and the default finish: [ExternalComponent.java](https://raw.githubusercontent.com/openrocket/openrocket/unstable/core/src/main/java/info/openrocket/core/rocketcomponent/ExternalComponent.java)
%[text] - Viscosity fit, speed of sound fit and heat capacity ratio in OpenRocket: [AtmosphericConditions.java](https://raw.githubusercontent.com/openrocket/openrocket/unstable/core/src/main/java/info/openrocket/core/models/atmosphere/AtmosphericConditions.java)
%[text] - Exported friction drag coefficient and altitude above sea level: [FlightDataType.java](https://raw.githubusercontent.com/openrocket/openrocket/unstable/core/src/main/java/info/openrocket/core/simulation/FlightDataType.java)
%[text] - Atmosphere layers and geopotential altitude: [Wikipedia, International Standard Atmosphere](https://en.wikipedia.org/wiki/International_Standard_Atmosphere)
%[text] - Barometric formula and its constants: [Wikipedia, Barometric formula](https://en.wikipedia.org/wiki/Barometric_formula)
%[text] - Speed of sound: [NASA Glenn, Speed of Sound](https://www.grc.nasa.gov/www/BGH/sound.html)
%[text] - Athena diameter, fin count and finish: Aegis Technical Baseline.pdf. Fin planform and thickness: Fins CDR 2025.pptx. Peak Mach 2.5: Nose Cone Profiles.pdf. \

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline"}
%---
