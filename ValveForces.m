%Program to determne forces in the valve system over the course of its
%cycle

k = 2.8e3; %C1125-085-1500S spring constant, Nm^-1
P1 = 1e5; %1 Bar pressure, Pa
DBore1 = 40e-3; %Outer tube bore diameter, m
TolBore1 = [62 0]*1e-6; %H9 bore tolerance, m
DShaft1 = 40e-3; %Inner tube shaft diameter, m
TolShaft1 = [-50 -112]*1e-6; %e9 shaft tolerance, m
DBore2 = 16e-3; %Inner tube bore diameter, m
TolBore2 = [18 0]*1e-6; %H7 bore tolerance, m
DShaft2 = 16e-3; %Inner tube bore diameter, m
TolShaft2 = [-17 -6]*1e-6; %g6 shaft tolerance, m
tGroove1Max = 2.77e-3; %Servo piston head max groove depth, m
tGroove1Min = 2.7e-3; %Servo piston head min groove depth, m
tGroove2Max = 2.2e-3; %Spool rod max groove depth, m
tGroove2Min = 2.13e-3; %Spool rod min groove depth, m
tGroove3Max = 3.42e-3; %Spool rod-spool piston head max chamfer, m
tGroove3Min = 3.3e-3; %Spool rod-spool piston head min chamfer, m

AGroove2Max = pi*((DShaft2 +TolShaft2).^2 - ((DShaft2 +TolShaft2) -2*tGroove2Max).^2)/4;
AGroove2Min = pi*((DShaft2 +TolShaft2).^2 - ((DShaft2 +TolShaft2) -2*tGroove2Min).^2)/4;
AGroove3Max = pi*((DShaft2 +TolShaft2 + 2*tGroove3Max).^2 - (DShaft2 +TolShaft2).^2)/4;
AGroove3Min = pi*((DShaft2 +TolShaft2 + 2*tGroove3Min).^2 - (DShaft2 +TolShaft2).^2)/4;

delPL2 = sym('delPL2');
eqnPreLoad2 = P1*(max(AGroove3Max) - min(AGroove2Min)) == k*delPL2;
delLPreLoad2 = vpasolve(eqnPreLoad2, delPL2)

IDORing1 = 34e-3;
CSORing1 = 3e-3;
TolCSORing1 = [0.1 -0.1]*1e-3;
LpORing1 = pi*max(DBore1 + TolBore1);
ApORing1Max = pi*(max(DBore1 + TolBore1)^2 - min(DShaft1 + TolShaft1 - 2*tGroove1Max)^2)/4;
ApORing1Min = pi*(max(DBore1 + TolBore1)^2 - min(DShaft1 + TolShaft1 - 2*tGroove1Min)^2)/4;
IDORing2 = 11.6e-3;
CSORing2 = 2.4e-3;
TolCSORing2 = [0.08 -0.08]*1e-3;
LpORing2 = pi*max(DBore2 + TolBore2);
ApORing2Max = pi*(max(DBore2 + TolBore2)^2 - min(DShaft2 + TolShaft2 - 2*tGroove2Max)^2)/4;
ApORing2Min = pi*(max(DBore2 + TolBore2)^2 - min(DShaft2 + TolShaft2 - 2*tGroove2Min)^2)/4;

SealCom1Max = abs(((0.5*(min(DBore1 +TolBore1) - max((DShaft1 +TolShaft1 - 2*tGroove1Min)))/max(CSORing1 + TolCSORing1))-1)*100);
SealCom1Min = abs(((0.5*(max(DBore1 +TolBore1) - min((DShaft1 +TolShaft1 - 2*tGroove1Max)))/min(CSORing1 + TolCSORing1))-1)*100);
SealCom2Max = abs(((0.5*(min(DBore2 +TolBore2) - max((DShaft2 +TolShaft2 - 2*tGroove2Min)))/max(CSORing2 + TolCSORing2))-1)*100);
SealCom2Min = abs(((0.5*(max(DBore2 +TolBore2) - min((DShaft2 +TolShaft2 - 2*tGroove2Max)))/min(CSORing2 + TolCSORing2))-1)*100);
SealComArray1 = [SealCom1Max SealCom1Min];
SealComArray2 = [SealCom2Max SealCom2Min];

lbf = 4.448;
in  = 25.4e-3;

SealCom = sym('SealCom');
grad1 = 1.5*(lbf/in)/22.5;
eqnfc = grad1*SealCom;
fh = 10*lbf/in^2;
fc1 = vpa(subs(eqnfc, SealCom, SealComArray1));
fc2 = vpa(subs(eqnfc, SealCom, SealComArray2));

Ffriction1 = vpa(fc1*LpORing1 + fh*ApORing1Max);
Ffriction2 = vpa(fc2*LpORing2 + fh*ApORing2Max);

FServo = P1*pi*min(DBore1 + TolBore1)^2/4;
FNC1Max = P1*min((AGroove3Max - AGroove2Max));
FNC1Min = P1*min((AGroove3Min - AGroove2Min));
FPreLoad2 = k*delLPreLoad2;

delL1 = sym('delL1');
%eqnFinal1 = (FNC1Max + FServo - FPreLoad2 - max(Ffriction1 + Ffriction2)) == k*delL1;
%eqnFinal2 = (FNC1Min + FServo - FPreLoad2 - max(Ffriction1 + Ffriction2)) == k*delL1;
eqnSpring1 = (FServo - max(Ffriction1)) == k*delL1;
delL1F1 = solve(eqnSpring1, delL1)
%delL1F2 = solve(eqnFinal2, delL1)
FNet1 = sym('FNet1');
delLSpring1 = 10e-3;
eqnTotalNet1 = FServo - max(Ffriction1) - k*delLSpring1 == FNet1;
FNet1 = solve(eqnTotalNet1, FNet1)
delL2 = sym('delL2');
eqnSpring2 = (FServo + FNC1Min - max(Ffriction1 + Ffriction2) - FPreLoad2) == k*delL2;
delL1F2 = solve(eqnSpring2, delL2)

fprintf('Max percentage compression for 2.7mm groove = \n');
vpa(SealCom1Max,5)
fprintf('Min percentage compression for 2.77mm groove = \n');
vpa(SealCom1Min,4)
fprintf('Min bore diameter used for 2.7mm groove = \n');
vpa((min(DBore1 +TolBore1)),6)
fprintf('mm \n');
fprintf('Max bore diameter used for 2.77mm groove = \n');
vpa(max(DBore1 +TolBore1),6)
fprintf ('mm \n');
fprintf('fc used for 2.7mm groove = \n');
vpa(fc1(1),6)
fprintf ('Nm^-1 \n');
fprintf('fc used for 2.77mm groove = \n');
vpa(fc1(2),6)
fprintf ('Nm^-1 \n');
fprintf('fh used for both grooves = \n');
fh
fprintf ('Nm^-2 \n');
fprintf('Lp used for both grooves = \n');
LpORing1
fprintf ('m \n');
fprintf('Ap used for 2.7mm groove = \n');
ApORing1Min
fprintf ('m^2 \n');
fprintf('Ap used for 2.77mm groove = \n');
ApORing1Max
fprintf ('m^2 \n');
