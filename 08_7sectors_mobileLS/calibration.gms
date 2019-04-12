* ======================================================
* 7 Sectors
* Calibration file
* ======================================================

$INCLUDE data.gms

*++++++++++++++++++++++++++++++++++++++++++++++++++++++++
* Initializing
*++++++++++++++++++++++++++++++++++++++++++++++++++++++++

* Flag for Endogenous Labor Supply
PARAMETER
FlEndo		Flag (0) Exogenous (1) Endogenous Labor Supply
;


* Model Parameters
SigInter	= 1;
SigIntra	= 3;
sigCon(g)	= 2.5;
sigInv		= 3;
sigLdem(s)	= 1.5;
SigGov		= 2.5;
sigX(s)		= 2.5;
Rint0   	= 1.40;

IF(CARD(g) EQ 4,
	EP(g)   = 1 + 0.35*ORD(g) - 0.09*(ORD(g)**2);	
	EP(gr)  = 0;
	delta   = 0.2;
    );

Leis0(gf,q)	= 0.35;
Leis0("g2",q)	= 0.20;
Leis0("g3",q)	= 0.90;

* Specifying Earnings Profiles for each qualification 
EPQ(g,"q1")	= EP(g)*1.2;
EPQ(g,"q2")	= EP(g)*0.8;
DISPLAY EPQ;

Lab0(q)   =	(SUM(s,LQ0(s,q)) / SUM(g,PopQ0(q,g)*EPQ(g,q)));

DISPLAY Lab0;

BeqR(g) 	= 0;
InhR(g) 	= 0;
rho0    	= 0.80;
AlDemQ(s,q)	= LQ0(s,q)/Ldem0(s);
LsupQ0(q) 	= (SUM(g, PopQ0(q,g)*Lab0(q)*EPQ(g,q)));
Lsup0		= SUM((q),LsupQ0(q));


*++++++++++++++++++++++++++++++++++++++++++++++++++++++++
* Calibration
*++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Y0(s)		= YK0(s) + Ldem0(s) + X0(s);
I0		= SUM(s,IS0(s));
K0              = I0 / (delta + GPop0 - 1);
AlK(s)$Y0(s)    = YK0(s) / Y0(s);
AlX(s)$Y0(s)	= X0(s) / Y0(s);
Ldem0(s)	= (1-AlK(s)-AlX(s))*Y0(s);
C0		= SUM(s, CS0(s));
Con0(g,q)	= C0 / SUM(gg,Pop0(gg));
DIFFC0		= C0 - SUM((g,q),PopQ0(q,g)*Con0(g,q));
	DISPLAY DIFFC0
	ABORT$(ABS(DIFFC0) GT 1.E-10) "Consumption not balanced";


*++++++++++++++++++++++++++++++++++++++++++++++++++++++++
* Equations
*++++++++++++++++++++++++++++++++++++++++++++++++++++++++

EQUATIONS
    HBudg0Eq1(g,q)    	HH Budget Constraint
    HBudg0Eq2(g,q)    	HH Budget Constraint last generation
    GammaEq(g,q)	Consumption intensity parameter calibration 
    Beq0Eq(g,q)       	Bequests
    Inh0Eq(g,q) 	Inheritances
    VV0vEq(g,q)		Definition of V	
    Con0Eq(g,q)       	Intertemporal Consumption (Euler Equation)
    C0Eq                Aggregate Consumption
    Rint0VEq            Balance of interest and rental rates
    Asset0Eq            Asset Market Equilibrium
    GBudg0Eq		Government Budget Balance
    OBJEq               Objective Function
;

* HH Budget Constraint
HBudg0Eq1(g+1,q)..
    (1+CTxR0)*Con0v(g,q) + B0v(g+1,q)
    =E=
    ((1-WTxR0)*Lab0(q)*EPQ(g,q)) +
    ((Rint0v-KTxR0*(Rint0v-1))*B0v(g,q))$(ORD(g) GT 1)
    + Inh0v(g,q) - Beq0v(g,q)
    ;

* HH Budget Constraint last generation
HBudg0Eq2(gl,q)..
    (1+CTxR0)*Con0v(gl,q)
    =E=
    (1-WTxR0)*Lab0(q)*EPQ(gl,q) +
    (Rint0v-KTxR0*(Rint0v-1))*B0v(gl,q) +
    Inh0v(gl,q) - Beq0v(gl,q)
    ;

* Bequests
Beq0Eq(gl,q)..
    Beq0v(gl,q)
    =E=
    BeqR(gl)*Con0v(gl,q)
    ;

* Inheritance
Inh0Eq(gw,q)..
    PopQ0(q,gw)*Inh0v(gw,q)
    =E=
    InhR(gw)*SUM(gr, PopQ0(q,gr)*Beq0v(gr,q))
    ;

* Intertemporal Consumption (Euler Equation)
Con0Eq(g+1,q)..
    Con0v(g+1,q)/Con0v(g,q)
    =E=
    ((Rint0v-KTxR0*(Rint0v-1))/rho0v)**SigInter 
    ;

* Aggregate Consumption
C0Eq..
    SUM((g,q), PopQ0(q,g)*Con0v(g,q))
    =E=
    C0
    ;

* Balance of interest and rental rates
Rint0vEq..
    Rint0v - 1
    =E=
    Rent0v - delta
    ;

* Asset Market Equiliubrium
Asset0Eq..
    SUM((g,q), PopQ0(q,g+1)*B0v(g+1,q))
    =E=
    K0 + Bond0v
    ;

* Government Budget Balance
GBudg0Eq..
    Gpop0*Bond0v +
    (SUM((g,q),PopQ0(q,g)*(WTxR0*Lab0(q)*EPQ(g,q)+
    CTxR0*Con0v(g,q)+KTxR0*(Rint0v-1)*B0v(g,q))))
    =E=
    G0 + Rint0v*Bond0v
    ;

* Objective Equation
ObjEq..
    Obj
    =E=
    0
    ;

*++++++++++++++++++++++++++++++++++++++++++++++++++++++++
* Model Definition
*++++++++++++++++++++++++++++++++++++++++++++++++++++++++

MODEL OLG0 /
      HBudg0Eq1
      HBudg0Eq2
      Beq0Eq
      Inh0Eq
      Con0Eq
      C0Eq
      Rint0vEq
      Asset0Eq
      GBudg0Eq
      ObjEq
      /
;

* Treat fixed variables as constants
OLG0.HOLDFIXED=1;


*++++++++++++++++++++++++++++++++++++++++++++++++++++++++
* Initial Values and Bounds
*++++++++++++++++++++++++++++++++++++++++++++++++++++++++

B0v.L(g,q)    	=       K0/SUM(gg$(ORD(gg) GT 1), Pop0(gg));
B0v.FX(gf,q)  	=       0;
Rint0v.L        =       Rint0;
Rint0v.LO       =       0.75;
Rint0v.UP       =       2.5;
Rent0v.L        =       Rint0v.L-(1-delta);
Rent0v.LO       =       0.01;
rho0v.L         =       rho0;
rho0v.LO        =       0.5*rho0v.L;
rho0v.UP        =       2.5;
Con0v.L(g,q)  	=       Con0(g,q);
Con0v.LO(g,q)	=	0.01;
Beq0v.L(gr,q) 	=       BeqR(gr)*Con0(gr,q);
Beq0v.FX(g,q)$(ORD(g) LT CARD(g)) =       BeqR(g);
Inh0v.L(gw,q) 	=       InhR(gw)*SUM(gr,Pop0(gr)*BeqR(gr)*Con0(gr,q))/Pop0(gw);
Inh0v.FX(gr,q)	=       InhR(gr);
VV0v.L(g,q)	=	1.5;
gammav.L(g,q)	=	1;
gammav.LO(g,q)	=	0;
Bond0v.L	=	Bond0;

*++++++++++++++++++++++++++++++++++++++++++++++++++++++++
* Options
*++++++++++++++++++++++++++++++++++++++++++++++++++++++++

OPTION DECIMALS = 6;
OPTIONS SOLPRINT=ON, LIMROW=0, LIMCOL=0, ITERLIM=1000;
OPTION NLP=CONOPT;

*++++++++++++++++++++++++++++++++++++++++++++++++++++++++
* Solve
*++++++++++++++++++++++++++++++++++++++++++++++++++++++++

SOLVE OLG0 MAXIMIZING OBJ USING NLP;
ABORT$(OLG0.ModelStat GT 2)"Check Calibration: No Feasible Solution for ---> Model OLG0";

*++++++++++++++++++++++++++++++++++++++++++++++++++++++++
* Save Results
*++++++++++++++++++++++++++++++++++++++++++++++++++++++++

B0(g,q)       	=       B0v.L(g,q);
Con0(g,q)     	=       Con0v.L(g,q);
rho0            =       rho0v.L;
Beq0(g,q)     	=       Beq0v.L(g,q);
Inh0(g,q)     	=       Inh0v.L(g,q);
Rent0           =       Rent0v.L;
Rint0           =       Rint0v.L;
Kdem0(s)        =       YK0(s)/Rent0;
Ksc0            =       SUM(s,Kdem0(s))/K0;
A0(s)           =       Y0(s) / (Kdem0(s)**AlK(s) * X0(s)**AlX(s) * Ldem0(s)**(1-AlK(s)-AlX(s)));
Gamma(g,q)	=	Gammav.L(g,q);
VV0(g,q)	=	VV0v.L(g,q);
Bond0		=	Bond0v.L;
AlInvS(s)$I0	=	IS0(s)/I0;
AlGovS(s)$G0	=	GS0(s)/G0;
AlXS(s,ss)	=	Input0(s,ss)/X0(ss);
	TESTJ		=	SUM(s,AlInvS(s))-1;
	ABORT$(ABS(TESTJ) GT 1.E-7) "InvS(s): shares do not sum to one";
	

$ONTEXT
*++++++++++++++++++++++++++++++++++++++++++++++++++++++++
* Equilibrium Tests
*++++++++++++++++++++++++++++++++++++++++++++++++++++++++

PARAMETER
    Walras0
;
Walras0         =       INF;
Walras0         =       Y0 - SUM(g, Pop0(g)*Con0(g)) - I0;
Walras0$(ABS(Walras0) LT 1.E-7) =       0;
DISPLAY Walras0;

PARAMETER
    WAssetBal0
;
WAssetBal0      =       INF;
WAssetBal0      =       SUM(g, Pop0(g+1)*B0(g+1)) - K0;
WAssetBal0$(ABS(WAssetBal0) LT 1.E-7)   =       0;
DISPLAY WAssetBal0;
$OFFTEXT


* Updating sectoral structure of consumption baskets by generations

$CALL GDXXRW.EXE ALCONSC.xlsx par=ALCONSC rng=sheet1!A1:H8
$GDXIN ALCONSC.gdx
$LOAD ALCONSC
DISPLAY ALCONSC;

AlConS(s,g,q)	=	ALCONSC(s,g)/SUM(ss,ALCONSC(ss,g));
ConS0(s,g,q)	=	AlConS(s,g,q)*Con0(g,q);
display AlConS;
EQUATIONS
    ConS0Eq1(g,q)
    ConS0Eq2(s)
    ObjConSEq
;

ConS0Eq1(g,q)..
    Con0(g,q)
    =E=
    SUM(s,ConS0v(s,g,q))
    ;

ConS0Eq2(s)..
    CS0(s)
    =E=
    SUM((g,q),PopQ0(q,g)*ConS0v(s,g,q))
    ;

ObjConSEq..
    OBJ
    =E=
    SUM((s,g,q),SQR(ConS0v(s,g,q)-ConS0(s,g,q)))
    ;

MODEL ModConS0 /
      ConS0Eq1
      ConS0Eq2
      ObjConSEq
      /
;

ConS0v.L(s,g,q)	=	ConS0(s,g,q);
ConS0v.LO(s,g,q)	=	EPS;
ConS0v.UP(s,g,q)	=	Con0(g,q);

SOLVE ModConS0 MINIMIZING OBJ USING NLP;
ABORT$(ModConS0.ModelStat GT 2)"Check Calibration: No Feasible Solution for ---> Model ModConS0";

ConS0(s,g,q)	=	ConS0v.L(s,g,q);
AlConS(s,g,q)	=	ConS0(s,g,q)/Con0(g,q);
	TESTJG(g,q)	=	SUM((s),AlConS(s,g,q))-1;
display TESTJG,AlConS;
	LOOP((g,q),
	ABORT$(ABS(TESTJG(G,Q)) GT 1.E-7) "ConS0(s,g): Shares do not sum to one";
	);

DISPLAY AlConS;



EXECUTE_UNLOAD 'calibration.gdx';

*EXECUTE '=gdx2xls calibration.gdx';
