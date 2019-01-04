#!/usr/bin/env wolframscript

(*scriptin' stuff*)

Needs["DifferentialEquations`NDSolveUtilities`"];

SetDirectory["/home/storluffarn/phd-stuff/research/graphene/code/"];

(*reference values*)

(* single slip*)
(*
c11=4.0*^-20;
c22=3.0;
\[Eta]n=5.0*^12;

k = 2.0; (*Li: 0.8*)
v = 1; (*Li: *2*)
v2=0.5;
c1=c11; (*Li: ~9*^-17?!*)
c2=c22;
c3=14*^20; (*3.5*^1 for q^2*)
c4=7.0; (*period is ~99.5?*)
c5=0.5c11; (* 2.5*^-21 *)
c6=0.5c22; (* 1 *)
c7=1.0; (*period is ~7.1072?*)
a = 2.5*^-10; (*Li: 2.5*^-10*)
mx = 1*^-23;
mq = 5*^-22;
\[Eta]x=1.25\[Eta]n;
\[Eta]q=\[Eta]n;
\[Zeta]x=0*^-2; 
\[Zeta]q=75.0*^21;
h = 1*^-15;
tmax=5000000h;
*)

(* vel dep *)
(*
c11=1.3*^-20;
c22=0.2;
\[Eta]n=3.0*^12;

k = 1.0; (*Li: 0.8*)
v = 1; (*Li: *2*)
v2=0.5;
c1=c11; (*Li: ~9*^-17?!*)
c2=c22;
c3=8*^18; (*1.0*^0 for q^2*)
c4=1.0; (*period is ~99.5?*)
c5=0.0c11; (* 2.5*^-21 *)
c6=0.0c22; (* 1 *)
c7=1.0; (*period is ~7.1072?*)
a = 2.5*^-10; (*Li: 2.5*^-10*)
mx = 1*^-23;
mq = 5*^-22;
\[Eta]x=1.25\[Eta]n;
\[Eta]q=\[Eta]n;
\[Zeta]x=0*^-2; 
\[Zeta]q=0.0*^21;
h = 1*^-15;
tmax=7500000h;
*)

(*parameter values*)

c11=4.0*^-20;
c22=3.0;
\[Eta]n=5.0*^12;

k = 2.0; (*Li: 0.8*)
v = 1; (*Li: *2*)
v2=0.5;
c1=c11; (*Li: ~9*^-17?!*)
c2=c22;
(* c3=14*^20; (*3.5*^1 for q^2*) *)
c3 = ToExpression@StringReplace[$ScriptCommandLine[[2]],"e"->"*^"];
c4=7.0; (*period is ~99.5?*)
c5=0.5c11; (* 2.5*^-21 *)
c6=0.5c22; (* 1 *)
c7=1.0; (*period is ~7.1072?*)
a = 2.5*^-10; (*Li: 2.5*^-10*)
mx = 1*^-23;
mq = 5*^-22;
\[Eta]x=1.25\[Eta]n;
\[Eta]q=\[Eta]n;
\[Zeta]x=0*^-2; 
\[Zeta]q=75.0*^21;
h = 1*^-15;
tmax=5000000h;

(*potential to be used*)

anpotstring=
"1/2k(x[t]-v t)^2+(c1+c2 q[t]^2)(1-Cos[2Pi(x[t]-c4 q[t])/a])+(c5+c6 q[t]^2)(1-Cos[c7 2Pi(x[t])/a])+c3 q[t]^4";
anpot=ToExpression[anpotstring];
	
(*diff solver, stiffness switching is good for numerically striff stiff problems*)
s=NDSolve[{-mx x''[t]-D[anpot,x[t]]-mx \[Eta]x (x'[t])==0,-mq q''[t]-D[anpot,q[t]]-mq \[Eta]q(q'[t])-mq \[Zeta]q q[t]==0,x[0]==0,x'[0]==0,q[0]==0,q'[0]==0},{x[t],q[t]},{t,0,tmax},MaxSteps->100000];

(* find extrema *)
extrema=Reap[k(v t-x[t]/.First[NDSolve[{-mx x''[t]-D[anpot,x[t]]-mx \[Eta]x (x'[t])==0,-mq q''[t]-D[anpot,q[t]]-mq \[Eta]q(q'[t])-mq \[Zeta]q q[t]==0,x[0]==0,x'[0]==0,q[0]==0,q'[0]==0},{x,q},{t,0,tmax},Method->{"EventLocator","Event"->D[k(v t-x[t]),t],"EventAction":>Sow[{t,k(v t-x[t])}]},MaxSteps->100000]])][[2,1]];

extrema={v #[[1]],#[[2]]}&/@extrema;

maxfric = Max[extrema];

meanfric = (extrema[[Length[extrema]]][[2]] + extrema[[Length[extrema]-1]][[2]])/2;

Print[meanfric]

meanfric>>>"fric_out.dat"

