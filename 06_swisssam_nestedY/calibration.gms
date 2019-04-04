* ======================================================
* 7 Sectors
* Imperfect Substitutability between Age groups
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
SigInter	= 6;
SigIntra	= 6;
sigCon(g)	= 2.5;
sigInv		= 5;
sigLdem(s)	= 3.5;
SigGov		= 2.5;
sigX(s)		= 5;
sigAge(s)	= 6;
sigVA(s)	= 1.5;	
Rint0   	= 1.4;

IF(CARD(g) EQ 4,
	EP(g)   = 1 + 0.35*ORD(g) - 0.09*(ORD(g)**2);	
	EP(gr)  = 0;
	delta   = 1;
    );

Leis0(gf,q,e)	= 0.35;
Leis0("g2",q,e)	= 0.20;
Leis0("g3",q,e)	= 0.90;
LeisS0(g,q,"s1")=Leis0(g,q,"e1");
LeisS0(g,q,"s2")=Leis0(g,q,"e2");
LeisS0(g,q,"s3")=Leis0(g,q,"e3");
LeisS0(g,q,"s4")=Leis0(g,q,"e4");
LeisS0(g,q,"s5")=Leis0(g,q,"e5");
LeisS0(g,q,"s6")=Leis0(g,q,"e6");
LeisS0(g,q,"s7")=Leis0(g,q,"e7");

* Specifying Earnings Profiles for each qualification 
EPQ(g,"q1")	= EP(g)*1.2;
EPQ(g,"q2")	= EP(g)*0.8;
DISPLAY EPQ;

LabS0(q,s,gw)  	= LQA0(s,q,gw) / (PopQS0(q,s,gw)*EPQ(gw,q));

Lab0(q,g,"e1")	= LabS0(q,"s1",g);
Lab0(q,g,"e2")	= LabS0(q,"s2",g);
Lab0(q,g,"e3")	= LabS0(q,"s3",g);
Lab0(q,g,"e4")	= LabS0(q,"s4",g);
Lab0(q,g,"e5")	= LabS0(q,"s5",g);
Lab0(q,g,"e6")	= LabS0(q,"s6",g);
Lab0(q,g,"e7")	= LabS0(q,"s7",g);
Lab0(q,g,"e8")	= LabS0(q,"s8",g);
DISPLAY Lab0;

BeqR(g) 	= 0;
InhR(g) 	= 0;
rho0    	= 0.70;
AlDemQ(s,q)	= LQ0(s,q)/Ldem0(s);
AlAge(s,q,g)	= LQA0(s,q,g)/LQ0(s,q);
LsupEQ0(e,q,g) 	= PopQE0(q,e,g)*Lab0(q,g,e)*EPQ(g,q);
Lsup0		= SUM((e,q,g),LsupEQ0(e,q,g));


*++++++++++++++++++++++++++++++++++++++++++++++++++++++++
* Calibration
*++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Y0(s)		= Va0(s) + X0(s);
I0		= SUM(s,IS0(s));
K0              = I0 / (delta + GPop0 - 1);
AlK(s)$Y0(s)    = YK0(s) / VA0(s);
AlX(s)$Y0(s)	= X0(s) / Y0(s);
Ldem0(s)	= (1-AlK(s))*VA0(s);
C0		= SUM(s, CS0(s));
Con0(g,q,e)	= C0 / SUM(gg,Pop0(gg));
DIFFC0		= C0 - SUM((g,q,e),PopQE0(q,e,g)*Con0(g,q,e));
	DISPLAY DIFFC0
	ABORT$(ABS(DIFFC0) GT 1.E-10) "Consumption not balanced";


*++++++++++++++++++++++++++++++++++++++++++++++++++++++++
* Equations
*++++++++++++++++++++++++++++++++++++++++++++++++++++++++

EQUATIONS
    HBudg0Eq1(g,q,e)    HH Budget Constraint
    HBudg0Eq2(g,q,e)    HH Budget Constraint last generation
    GammaEq(g,q,e)	Consumption intensity parameter calibration 
    Beq0Eq(g,q,e)       Bequests
    Inh0Eq(g,q,e)       Inheritances
    VV0vEq(g,q,e)	Definition of V	
    Con0Eq(g,q,e)       Intertemporal Consumption (Euler Equation)
    C0Eq                Aggregate Consumption
    Rint0VEq            Balance of interest and rental rates
    Asset0Eq            Asset Market Equilibrium
    GBudg0Eq		Government Budget Balance
    OBJEq               Objective Function
;

* HH Budget Constraint
HBudg0Eq1(g+1,q,e)..
    (1+CTxR0)*Con0v(g,q,e) + B0v(g+1,q,e)
    =E=
    ((1-WTxR0(q))*Lab0(q,g,e)*EPQ(g,q)) +
    ((Rint0v-KTxR0*(Rint0v-1))*B0v(g,q,e))$(ORD(g) GT 1)
    + Inh0v(g,q,e) - Beq0v(g,q,e)
    ;

* HH Budget Constraint last generation
HBudg0Eq2(gl,q,e)..
    (1+CTxR0)*Con0v(gl,q,e)
    =E=
    (1-WTxR0(q))*Lab0(q,gl,e)*EPQ(gl,q) +
    (Rint0v-KTxR0*(Rint0v-1))*B0v(gl,q,e) +
    Inh0v(gl,q,e) - Beq0v(gl,q,e)
    ;

* Bequests
Beq0Eq(gl,q,e)..
    Beq0v(gl,q,e)
    =E=
    BeqR(gl)*Con0v(gl,q,e)
    ;

* Inheritance
Inh0Eq(gw,q,e)..
    PopQE0(q,e,gw)*Inh0v(gw,q,e)
    =E=
    InhR(gw)*SUM(gr, PopQE0(q,e,gr)*Beq0v(gr,q,e))
    ;

* Intertemporal Consumption (Euler Equation)
Con0Eq(g+1,q,e)..
    Con0v(g+1,q,e)/Con0v(g,q,e)
    =E=
    ((Rint0v-KTxR0*(Rint0v-1))/rho0v)**SigInter 
    ;

* Aggregate Consumption
C0Eq..
    SUM((g,q,e), PopQE0(q,e,g)*Con0v(g,q,e))
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
    SUM((g,q,e), PopQE0(q,e,g+1)*B0v(g+1,q,e))
    =E=
    K0 + Bond0v
    ;

* Government Budget Balance
GBudg0Eq..
    Gpop0*Bond0v +
    (SUM((g,q,e),PopQE0(q,e,g)*(WTxR0(q)*Lab0(q,g,e)*EPQ(g,q)+
    CTxR0*Con0v(g,q,e)+KTxR0*(Rint0v-1)*B0v(g,q,e))))
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

B0v.L(g,q,e)    =       K0/SUM(gg$(ORD(gg) GT 1), Pop0(gg));
B0v.FX(gf,q,e)  =       0;
Rint0v.L        =       Rint0;
Rint0v.LO       =       .1;
Rint0v.UP       =       2.5;
Rent0v.L        =       Rint0v.L-(1-delta);
Rent0v.LO       =       .1;
rho0v.L         =       rho0;
rho0v.LO        =       0.1*rho0v.L;
rho0v.UP        =       2.5;
Con0v.L(g,q,e)  =       Con0(g,q,e);
Con0v.LO(g,q,e)	=	0.01;
Beq0v.L(gr,q,e) =       BeqR(gr)*Con0(gr,q,e);
Beq0v.FX(g,q,e)$(ORD(g) LT CARD(g)) =       BeqR(g);
Inh0v.L(gw,q,e) =       InhR(gw)*SUM(gr,Pop0(gr)*BeqR(gr)*Con0(gr,q,e))/Pop0(gw);
Inh0v.FX(gr,q,e)=       InhR(gr);
VV0v.L(g,q,e)	=	1.5;
gammav.L(g,q,e)	=	1;
gammav.LO(g,q,e)=	0;
Bond0v.L	=	Bond0;
*Bond0v.fx 	= 	0;

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

B0(g,q,e)       =       B0v.L(g,q,e);
Con0(g,q,e)     =       Con0v.L(g,q,e);
rho0            =       rho0v.L;
Beq0(g,q,e)     =       Beq0v.L(g,q,e);
Inh0(g,q,e)     =       Inh0v.L(g,q,e);
Rent0           =       Rent0v.L;
Rint0           =       Rint0v.L;
Kdem0(s)        =       YK0(s)/Rent0;
Ksc0            =       SUM(s,Kdem0(s))/K0;
A0(s)           =       Y0(s) / (((1-AlX(s))*VA0(s)**((SigX(s)-1)/SigX(s)) + AlX(s)*X0(s)**((SigX(s)-1)/SigX(s)))**(SigX(s)/(SigX(s)-1)));
Gamma(g,q,e)	=	Gammav.L(g,q,e);
VV0(g,q,e)	=	VV0v.L(g,q,e);
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

$CALL GDXXRW.EXE ALCONSC.xlsx par=ALCONSC rng=sheet1!A1:E9
$GDXIN ALCONSC.gdx
$LOAD ALCONSC
DISPLAY ALCONSC;

AlConS(s,g,q,e)	=	ALCONSC(s,g)/SUM(ss,ALCONSC(ss,g));
ConS0(s,g,q,e)	=	AlConS(s,g,q,e)*Con0(g,q,e);
display AlConS;
EQUATIONS
    ConS0Eq1(g,q,e)
    ConS0Eq2(s)
    ObjConSEq
;

ConS0Eq1(g,q,e)..
    Con0(g,q,e)
    =E=
    SUM(s,ConS0v(s,g,q,e))
    ;

ConS0Eq2(s)..
    CS0(s)
    =E=
    SUM((g,q,e),PopQE0(q,e,g)*ConS0v(s,g,q,e))
    ;

ObjConSEq..
    OBJ
    =E=
    SUM((s,g,q,e),SQR(ConS0v(s,g,q,e)-ConS0(s,g,q,e)))
    ;

MODEL ModConS0 /
      ConS0Eq1
      ConS0Eq2
      ObjConSEq
      /
;

ConS0v.L(s,g,q,e)	=	ConS0(s,g,q,e);
ConS0v.LO(s,g,q,e)	=	EPS;
ConS0v.UP(s,g,q,e)	=	Con0(g,q,e);

SOLVE ModConS0 MINIMIZING OBJ USING NLP;
ABORT$(ModConS0.ModelStat GT 2)"Check Calibration: No Feasible Solution for ---> Model ModConS0";

ConS0(s,g,q,e)	=	ConS0v.L(s,g,q,e);
AlConS(s,g,q,e)	=	ConS0(s,g,q,e)/Con0(g,q,e);
	TESTJG(g,q,e)	=	SUM((s),AlConS(s,g,q,e))-1;
display TESTJG,AlConS;
	LOOP((g,q,e),
	ABORT$(ABS(TESTJG(G,Q,E)) GT 1.E-7) "ConS0(s,g): Shares do not sum to one";
	);

DISPLAY AlConS;



EXECUTE_UNLOAD 'calibration.gdx';

EXECUTE '=gdx2xls calibration.gdx';
