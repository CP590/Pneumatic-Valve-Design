k = 2.8e3; %C1125-085-1500S spring constant, Nm^-1
P1 = 1e5; %1 Bar pressure, Pa
P10 = 10e5; %10 Bar pressure, Pa
DBore1 = 40e-3; %Outer tube bore diameter, m
DShaft1 = 40e-3; %Piston head and inner tube shaft diameter, m
DBore2 = 16e-3; %Inner tube bore diameter, m
DBore3 = 0;
DShaft2 = 16e-3; %Servo & spool rod diameter, m
IDORing1 = 34e-3; %Piston head O-ring inner diameter, m
CSORing1 = 3e-3; %Piston head O-ring cross-sectional diameter, m
IDORing2 = 11.6e-3; %Spool rod O-ring inner diameter, m
CSORing2 = 2.4e-3; %Spool rod O-ring cross-sectional diameter, m
IDORing3 = 16e-3; %Spool rod-spool piston head O-ring inner diameter, m
CSORing3 = 1.5e-3; %Spool rod-spool piston head O-ring cross-sectional diameter, m
lbf = 4.448; %lbf to N conversion factor
in  = 25.4e-3; %in to m conversion factor

fprintf('Note that all inputted data does not need converting into metres and should be inputted based off of the units displayed in brackets \n');
pause(1);
prompt = 'Input outer tube H9 bore tolerance (0 to 62 microns) \n';
TolBore1 = input(prompt)*1e-6;
prompt2 = 'Input servo piston head e9 shaft tolerance (-50 to -112 microns) \n';
TolShaft1 = input(prompt2)*1e-6;
prompt3 = 'Input inner tube H7 bore tolerance (0 to 18 microns) \n';
TolBore2 = input(prompt3)*1e-6;
prompt4 = 'Input servo/spool rod g6 shaft tolerance (-6 to -17 microns) \n';
TolShaft2 = input(prompt4)*1e-6;
prompt5 = 'Input servo piston head groove depth (2.7 to 2.77 millimetres) \n';
tGroove1 = input(prompt5)*1e-3;
prompt6 = 'Input servo piston head O-Ring CS tolerance (-0.1 to 0.1 millimetres) \n';
TolCSORing1 = input(prompt6)*1e-3;
prompt7 = 'Input spool rod groove depth (2.13 to 2.2 millimetres) \n';
tGroove2 = input(prompt7)*1e-3;
prompt8 = 'Input servo spool rod O-Ring CS tolerance (-0.08 to 0.08 millimetres) \n';
TolCSORing2 = input(prompt8)*1e-3;
prompt9 = 'Input spool rod-spool piston head chamfer (2.08 to 2.2 millimetres) \n';
tGroove3 = input(prompt9)*1e-3;
prompt10 = 'Input spool rod-spool piston head O-Ring CS tolerance (-0.08 to 0.08 millimetres) \n';
TolCSORing3 = input(prompt10)*1e-3;

%AGroove2 = pi*((DBore2 +TolBore)^2 - ((DShaft2 +TolShaft2) - 2*tGroove2)^2)/4;
%%AGroove3 = pi*((DShaft2 +TolShaft2 + 2*tGroove3)^2 - (DShaft2 + TolShaft2)^2)/4;
%AGroove3 = pi*((20.525e-3 + TolBore2)^2 - (12.9877e-3 + TolShaft2)^2)/4;

LpORing1 = pi*(DBore1 + TolBore1);
ApORing1 = pi*((DBore1 + TolBore1)^2 - (DShaft1 + TolShaft1 - 2*tGroove1)^2)/4;
LpORing2 = pi*(DBore2 + TolBore2);
ApORing2 = pi*((DBore2 + TolBore2)^2 - (DShaft2 + TolShaft2 - 2*tGroove2)^2)/4;
%LpORing3 = pi*(DBore2 + TolBore2 + tGroove3);
LpORing3 = pi*(DBore2 + TolBore2 + 4*tGroove3/3);
%ApORing3 = pi*((DBore2 + TolBore2 + tGroove3)^2 - (DShaft2 + TolShaft2)^2)/4;
ApORing3 = pi*((DBore2 + 2*tGroove3 + TolBore2)^2 - (DShaft2 + TolShaft2)^2)/4;
SealCom1 = abs(((0.5*((DBore1 +TolBore1) - ((DShaft1 + TolShaft1 - 2*tGroove1)))/(CSORing1 + TolCSORing1))-1)*100);
SealCom2 = abs(((0.5*((DBore2 +TolBore2) - ((DShaft2 +TolShaft2 - 2*tGroove2)))/(CSORing2 + TolCSORing2))-1)*100);
%SealCom3 = abs(((0.5*((DBore2 +TolBore2 + tGroove3) - ((DShaft2 + TolShaft2)))/(CSORing3 + TolCSORing3))-1)*100);
%SealCom3 = abs(((0.5*((DBore2 + TolBore2) - (DShaft2 + TolShaft2 - 4*tGroove3/3))/(CSORing3 + TolCSORing3))-1)*100);

delPL2 = sym('delPL2');
FNC1 = P1*((ApORing3 - ApORing2));
FNC10 = P10*((ApORing3 - ApORing2));
eqnPreLoad2 = FNC10 == k*delPL2;
delLPreLoad2 = vpasolve(eqnPreLoad2, delPL2); %Spring 2 preload extension solved

SealCom = sym('SealCom');
grad1 = 1.5*(lbf/in)/22.5;
eqnfc = grad1*SealCom;
fh = 10*lbf/in^2;
fc1 = vpa(subs(eqnfc, SealCom, SealCom1));
fc2 = vpa(subs(eqnfc, SealCom, SealCom2));
%fc3 = vpa(subs(eqnfc, SealCom, SealCom3));

Ffriction1 = vpa(fc1*LpORing1 + fh*ApORing1);
Ffriction2 = vpa(fc2*LpORing2 + fh*ApORing2);
%Ffriction3 = vpa(fc3*LpORing3 + fh*ApORing3);

FServo = P1*pi*(DBore1 + TolBore1)^2/4;
FPreLoad2 = k*delLPreLoad2;

delL1 = sym('delL1');
delLPreLoad1 = 2e-3;
FPreLoad1 = k*delLPreLoad1;
eqnSpring1 = FServo - Ffriction1 == k*(delL1 + delLPreLoad1);
delL1F1 = solve(eqnSpring1, delL1); %Distance servo rod can be pushed until equilibrium, m
FNet1 = sym('FNet1');
delLSpring1 = 4e-3; %Set distance between servo rod and spool rod, m
eqnTotalNet1 = FServo - Ffriction1- k*delLSpring1 == FNet1;
FNet1 = solve(eqnTotalNet1, FNet1); %Net force when servo rod hits spool rod, N
delL2 = sym('delL2');
eqnSpring2 = (FServo + FNC1 - (Ffriction1 + Ffriction2) - FPreLoad2 - FPreLoad1 - k*(delLSpring1 + delL2)) == k*delL2;
delL1F2 = solve(eqnSpring2, delL2); %Distance spool rod can be pushed until equilibrium, m

fprintf('Percentage compression for servo piston head O-ring = %3.3f \n', SealCom1);
vpa(SealCom1,5);
fprintf('Percentage compression for spool rod O-ring = %3.3f \n', SealCom2);
vpa(SealCom2,5);
fprintf('fc value for servo piston head O-ring = %3.3f Nm^-1 \n', fc1);
vpa(fc1,5);
fprintf('fc value for spool rod O-ring = %3.3f Nm^-1 \n', fc2);
vpa(fc2,5);
fprintf('fh value used for both O-rings = %5.0f Nm^-2 \n', fh);
vpa(fh,5);
fprintf('Lp value for servo piston head O-ring = %0.5f m  \n', LpORing1);
vpa(LpORing1,5);
fprintf('Ap value for servo piston head O-ring = %1.4e m^2 \n', ApORing1);
vpa(ApORing1,5);
fprintf('Lp value for spool rod O-ring = %0.5f m \n', LpORing2);
vpa(LpORing2,5);
fprintf('Ap value for spool rod O-ring = %1.4e m^2 \n', ApORing2);
vpa(ApORing2,5);
fprintf('Preload compression of spool spring = %1.4e m \n', delLPreLoad2);
vpa(delLPreLoad2,5);
fprintf('Max stroke displacement of servo piston = %1.4e m \n', delL1F1);
vpa(delL1F1,5);
fprintf('Net force on servo rod connection with spool rod = %3.3f N \n', FNet1);
vpa(FNet1,5);
fprintf('Max stroke displacement of spool piston = %1.4e m \n', delL1F2);
vpa(delL1F2,5);

promptc = 'Continue (Y/N)? : ';
if input(promptc, 's') == 'Y'
ValveForcesInteractive
else
end