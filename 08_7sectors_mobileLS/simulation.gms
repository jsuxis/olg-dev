* ========================================================
* 7 Sectors
* Simulation file
* ========================================================

$INCLUDE calibration.gms
$OFFSYMXREF OFFSYMLIST OFFUELLIST OFFUELXREF

 PARAMETERS
FlWalras	Flag (1) for using Walras' Law
FlDemog1     	Flag (1) for Demographic shock 1
FlDemog2     	Flag (1) for Demographic shock 2
FlProdu     	Flag (1) for Productivity shock
FlQual		Flag (1) for Qualification Shock
;

* --- Choose Flags

FlWalras	=	0;
FlDemog1	= 	0;
FlDemog2     	= 	0;
FlProdu     	= 	1;
FlQual		=	0;

*++++++++++++++++++++++++++++++++++++++++++++++++++++++++
* Equations
*++++++++++++++++++++++++++++++++++++++++++++++++++++++++

EQUATIONS
    QEq(s,t)            Production Function
    KdemEq(s,t)         Capital Demand
    LdemEq(s,t)         Labor Demand
    XdemEq(s,t)		Input Demand
    HBudgEq1(q,t,g)   	HH Budget Constraint
    HBudgEq2(q,t,g)   	HH Budget Constraint last generation
    BeqEq(q,t,g) 	Bequests
    InhEq(q,t,g)      	Inheritance
    ConEq(q,t,g)      	Euler Equation
    ConEqSS(q,t,g)    	Steady State Consumption profile
    LeisEq(q,t,g)	Intratemporal FOC
    MuEq(q,t,g)		Mu is only positive if not working (Kuhn-Tucker)
    VVEq(q,t,g)		Definition of V
    LsupEq(q,t)       	Labor Supply Balance
    RintEq(t)           Balance of interest and rental rates
    PConEq(q,t,g)	FOC of Intratemporal (sectoral) consumption allocation
    ConSEq(q,s,t,g)	FOC of Intratemporal (sectoral) consumption allocation
    PInvEq(t)		FOC of Investment origin allocation
    InvSEq(s,t)		FOC of Investment origin allocation
    WLdemEq(s,t)	FOC of qualification specific labor demand allocation
    LdemQEq(s,q,t)	FOC of qualification specific labor demand allocation
    XSdemEq(s,t)	FOC of Input allocation
    XdemSEq(s,ss,t)	FOC of Input allocation
    KStockEq(t)         Capital Stock Dynamics
    KStockEqSS(t)       Capital Stock in Steady State
    GBudgEq1(t)		Government Budget Balance 1
    GBudgEqSS(t)	Government Budget in Steady State
    GBudgEq2(t)		Government Budget Balance 2 (Constant Debt)
    PGovEq(t)		FOC of Government Consumption Allocation
    GovSEq(s,t)		FOC of Government Consumption Allocation
    PEq(s,t)            Goods Market Equilibrium
    WageEq(q,t)       	Labor Market Equilibrium
    RentEq(t)           Capital Market Equilibrium
    AssetEq(t)          Asset Market Equilibrium
    GovGEq(t)		Government Spending Growth
;


* Production Function (in logs)
QEq(s,t)        $(ORD(t) GT CARD(tp) AND ORD(t) LE CARD(tp)+CARD(tm))..
    LOG(Y(s,t))
    =E=
    LOG(A(s,t)) + AlK(s)*LOG(Kdem(s,t)) + AlX(s)*LOG(Xdem(s,t)) + (1-AlK(s)-AlX(s))*LOG(Ldem(s,t)) 
    ;

* Capital Demand (FOC1)
KdemEq(s,t)     $(ORD(t) GT CARD(tp) AND ORD(t) LE CARD(tp)+CARD(tm))..
    LOG(Rent(t)) + LOG(Kdem(s,t))
    =E=
    LOG(Alk(s)) + LOG(Y(s,t)) + LOG(P(s,t))
    ;

* Labor Demand (FOC2)
LdemEq(s,t)     $(ORD(t) GT CARD(tp) AND ORD(t) LE CARD(tp)+CARD(tm))..
    LOG(W(t)) + LOG(Ldem(s,t))
    =E=
    LOG((1-AlK(s)-AlX(s))) + LOG(Y(s,t)) + LOG(P(s,t))
    ;

* Input Demand (FOC3)
XdemEq(s,t)	$(ORD(t) GT CARD(tp) AND ORD(t) LE CARD(tp)+CARD(tm))..
    LOG(PX(s,t)) + LOG(Xdem(s,t))
    =E=
    LOG(AlX(s)) + LOG(Y(s,t)) + LOG(P(s,t))
    ;

* HH Budget Constraint
HBudgEq1(q,t+1,g+1)	$(ORD(t) GT CARD(tp) AND ORD(t) LT CARD(tp) + CARD(tm))..
    PCon(q,t,g)*(1+CTxR(t))*Con(q,t,g) + B(q,t+1,g+1)
    =E=
    ((1-WTxR(t))*wage(q,t)*Lab(q,t)*EPQ(g,q))+
    ((Rint(t)-KTxR(t)*(Rint(t)-1))*B(q,t,g))$(ORD(g) GT 1) +
    Inh(q,t,g) - Beq(q,t,g)
    ;

* HH Budget Constraint last generation
HBudgEq2(q,t,gl)	$(ORD(t) GT CARD(tp) AND ORD(t) LE CARD(tp)+CARD(tm))..
    PCon(q,t,gl)*(1+CTxR(t))*Con(q,t,gl)
    =E=
    ((1-WTxR(t))*wage(q,t)*Lab(q,t)*EPQ(gl,q)) +
    ((Rint(t)-KTxR(t)*(Rint(t)-1))*B(q,t,gl)) +
    Inh(q,t,gl) - Beq(q,t,gl)
    ;
* Bequests
BeqEq(q,t,g)      $(ORD(t) GT CARD(tp) AND ORD(t) LE CARD(tp)+CARD(tm) AND BeqR(g) NE 0)..
    Beq(q,t,g)
    =E=
    BeqR(g)*PCon(q,t,g)*Con(q,t,g)
    ;

* Inheritance
InhEq(q,t,g)     $(ORD(t) GT CARD(tp) AND ORD(t) LE CARD(tp)+CARD(tm) AND InhR(g) NE 0)..
    PopQ(q,t,g)*Inh(q,t,g)
    =E=
    InhR(g)*SUM(gr, PopQ(q,t,g)*Beq(q,t,g))
    ;

* Euler Equation
ConEq(q,t+1,g+1)  $(ORD(t) GT CARD(tp) AND ORD(t) LE CARD(tp)+CARD(tm)-1)..
    Con(q,t+1,g+1)/Con(q,t,g)
    =E=
    ((((Rint(t+1)-KTxR(t+1)*(Rint(t+1)-1))*PCon(q,t,g) / (rho*PCon(q,t+1,g+1)))*((1+CTxR(t))/(1+CTxR(t+1))))**SigInter)
    ;

* Steady State Consumption Profile
ConEqSS(q,t+1,g)  $((ORD(t) EQ (CARD(tp) + CARD(tm)-1) AND ORD(g) LT CARD(g)))..
    Con(q,t+1,g)
    =E=
    Con(q,t,g)
    ;

* FOC of Intratemporal (sectoral) consumption allocation
PConEq(q,t,g)	$(ORD(t) GT CARD(tp) AND ORD(t) LE CARD(tp)+CARD(tm))..
    PCon(q,t,g)**(1-SigCon(g))
    =E=
    SUM(s$(AlConS(s,g,q) GT 1.E-13), AlConS(s,g,q) * P(s,t)**(1-SigCon(g)))
    ;

* FOC of Intratemporal (sectoral) consumption allocation
ConSEq(q,s,t,g)	$(ORD(t) GT CARD(tp) AND ORD(t) LE CARD(tp)+CARD(tm)
		AND AlConS(s,g,q) GT 1.E-13)..
    ConS(s,q,t,g)
    =E=
    AlConS(s,g,q)*(PCon(q,t,g)/P(s,t))**SigCon(g)*Con(q,t,g)
    ;

* FOC of Investment origin allocation
PInvEq(t)	$(ORD(t) GT CARD(tp) AND ORD(t) LE CARD(tp)+CARD(tm))..
    PI(t)**(1-SigInv)
    =E=
    SUM(s$(AlInvS(s) GT 1.E-13), AlInvS(s)*P(s,t)**(1-SigInv))
    ;

* FOC of Investment origin allocation
InvSEq(s,t)	$(ORD(t) GT CARD(tp) AND ORD(t) LE CARD(tp)+CARD(tm)
		AND AlInvS(s) GT 1.E-13)..
    IS(s,t)
    =E=
    AlInvS(s)*(PI(t)/P(s,t))**SigInv * I(t)
    ;

* FOC of qualification specific labor demand allocation
WLdemEq(s,t)	$(ORD(t) GT CARD(tp) AND ORD(t) LE CARD(tp)+CARD(tm))..
    W(t)**(1-sigLdem(s))
    =E=
    SUM(q$(AlDemQ(s,q) GT 1.E-13), AlDemQ(s,q) * wage(q,t)**(1-sigLdem(s)))
    ;

* FOC of qualification specific labor demand allocation
LdemQEq(s,q,t)	$(ORD(t) GT CARD(tp) AND ORD(t) LE CARD(tp)+CARD(tm)
		AND AlDemQ(s,q) GT 1.E-13)..
    LQ(s,q,t)
    =E=
    AlDemQ(s,q)*(W(t)/wage(q,t))**sigLdem(s)*Ldem(s,t)
    ;
    
* FOC of sectoral inputs
XSdemEq(s,t)	$(ORD(t) GT CARD(tp) AND ORD(t) LE CARD(tp)+CARD(tm))..
    PX(s,t)**(1-sigX(s))
    =E=
    SUM(ss, AlXS(ss,s) * P(ss,t)**(1-sigX(s)))
    ;

* FOC of sectoral inputs
XdemSEq(s,ss,t)	$(ORD(t) GT CARD(tp) AND ORD(t) LE CARD(tp)+CARD(tm)
		AND AlXS(s,ss) GT 1.E-13)..
    Input(ss,s,t)
    =E=
    AlXS(ss,s)*(PX(s,t)/P(ss,t))**sigX(s)*Xdem(s,t)
    ;

* Labor Supply Balance
LsupEq(q,t)       $(ORD(t) GT CARD(tp) AND ORD(t) LE CARD(tp)+CARD(tm))..
    LsupQ(q,t)
    =E=
    (SUM(g, PopQ(q,t,g)*Lab(q,t)*EPQ(g,q)))
    ;

* Balance of interest and rental rates
RintEq(t)       $(ORD(t) GT CARD(tp) AND ORD(t) LE CARD(tp)+CARD(tm))..
    Rint(t)
    =E=
    (Rent(t)+1-delta*PI(t))/PI(t-1)
    ;

* Capital Stock Dynamics
KstockEq(t+1)   $(ORD(t) GT CARD(tp) AND ORD(t) LE CARD(tp) + CARD(tm)-1)..
    K(t+1)
    =E=
    I(t)+(1-delta)*K(t)
    ;

* Capital Stock in Steady State
KstockEqSS(t)   $((ORD(t) EQ CARD(tp)+CARD(tm)))..
    I(t)
    =E=
    ((GPop(t)-1) + delta)*K(t)
    ;

* Government Budget Constraint 1
GBudgEq1(t+1)	$(ORD(t) GT CARD(tp) AND ORD(t) LT CARD(tp)+CARD(tm))..
    PGov(t)*Bond(t+1)+(SUM((q,g),PopQ(q,t,g)*(
      WTxR(t)*wage(q,t)*Lab(q,t)*EPQ(g,q)
    + CTxR(t)*PCon(q,t,g)*Con(q,t,g)
    + KTxR(t)*(Rint(t)-1)*B(q,t,g)))) 
    =E=
    PGov(t)*Gov(t)+(Rint(t-1)*(PGov(t)/PGov(t-1)))*PGov(t-1)*Bond(t)
    ;

* Government Budget in Steady State
GBudgEqSS(t)	$(ORD(t) EQ CARD(tp)+CARD(tm))..
    GPop(t)*PGov(t)*Bond(t)+(SUM((q,g),PopQ(q,t,g)*(
      WTxR(t)*wage(q,t)*Lab(q,t)*EPQ(g,q)
    + CTxR(t)*PCon(q,t,g)*Con(q,t,g)
    + KTxR(t)*(Rint(t)-1)*B(q,t,g))))
    =E=
    PGov(t)*Gov(t)+Rint(t-1)*PGov(t)*Bond(t)
    ;

* Government Budget Constraint 2
GBudgEq2(t+1)	$(ORD(t) GT CARD(tp) AND ORD(t) LT CARD(tp)+CARD(tm))..
    PGov(t)*Bond(t+1)
    =E=
    GPop(t)*PGov(t-1)*Bond(t)
    ;

* FOC of Government consumption allocation
PGovEq(t)	$(ORD(t) GT CARD(tp) AND ORD(t) LE CARD(tp)+CARD(tm))..
    PGov(t)**(1-SigGov)
    =E=
    SUM(s$(AlGovS(s) GT 1.E-13), AlGovS(s) * P(s,t)**(1-SigGov))
    ;

* FOC of Government consumption allocation
GovSEq(s,t)	$(ORD(t) GT CARD(tp) AND ORD(t) LE CARD(tp)+CARD(tm)
		AND AlGovS(s) GT 1.E-13)..
    GovS(s,t)
    =E=
    AlGovS(s)*(PGov(t)/P(s,t))**SigGov * Gov(t)
    ;

* Goods market equilibrium
PEq(s,t)        $(ORD(t) GT CARD(tp) AND ORD(t) LE CARD(tp)+CARD(tm))..
    Y(s,t)
    =E=
    SUM((q,g), PopQ(q,t,g)*ConS(s,q,t,g)) + IS(s,t) + GovS(s,t)
    + SUM(ss,Input(s,ss,t))
    ;

* Labor market equilibrium
WageEq(q,t)       $(ORD(t) GT CARD(tp) AND ORD(t) LE CARD(tp)+CARD(tm))..
    LsupQ(q,t)
    =E=
    SUM(s,LQ(s,q,t))
    ;
        
* Capital Market Equilibrium
RentEq(t)       $(ORD(t) GT CARD(tp) AND ORD(t) LE CARD(tp)+CARD(tm))..
    Ksc0*K(t)
    =E=
    SUM(s,Kdem(s,t))
    ;

* Asset Market Equilibrium
 AssetEq(t+1)    $(ORD(t) GT CARD(tp) AND ORD(t) LE CARD(tp) + CARD(tm)-1 AND (FlWalras EQ 0))..
    SUM((g,q), PopQ(q,t+1,g+1)*B(q,t+1,g+1))
    =E=
    PI(t)*K(t+1) + PGov(t)*Bond(t+1)
    ;

*++++++++++++++++++++++++++++++++++++++++++++++++++++++++
* Model Definition
*++++++++++++++++++++++++++++++++++++++++++++++++++++++++

MODEL OLG /
      QEq
      KdemEq
      LdemEq
      XdemEq
      HBudgEq1
      HBudgEq2
      BeqEq
      InhEq
      ConEq
      ConEqSS
      LsupEq
      PConEq
      ConSEq
      PInvEq
      InvSEq
      WLdemEq
      LdemQEq
      RintEq
      KstockEq
      KstockEqSS
      GBudgEq1
      GBudgEq2
      GBudgEqSS
      PGovEq	
      GovSEq
      XSdemEq
      XdemSEq
      PEq
      WageEq
      RentEq
      AssetEq
      ObjEq
      /
;

* Treat fixed variables as constants
OLG.WORKSPACE = 1000;          
OLG.HOLDFIXED=1;


*++++++++++++++++++++++++++++++++++++++++++++++++++++++++
* Initialisation
*++++++++++++++++++++++++++++++++++++++++++++++++++++++++

* Population
NN(q,t)                           	=      	1;
PopQ(q,t,g) $(ORD(t) EQ CARD(g))  	=	PopQ0(q,g);

LOOP(t $(ORD(t) GT CARD(tp)),
    PopQ(q,t+1,g)                  	=       NN(q,t)*PopQ(q,t,g)
);

LOOP(tp,
    PopQ(q,tp,g)	                =       SUM(t $(ORD(t) EQ CARD(g)), PopQ(q,t,g))
);

TPop(t) $(ORD(t) LE CARD(t))    =       SUM((q,g), PopQ(q,t,g));

LOOP(t $(ORD(t) LE CARD(t)),
    GPop(t)                     =       TPop(t+1)/TPop(t)
);
Gov(t)				=	G0;
LOOP(t $(ORD(t) GT CARD(tp)),
    Gov(t+1)			=	GPop(t)*Gov(t)
);

* Parameters
A(s,t)          =       A0(s);
rho		=       rho0;
Lab(q,t)      	=       Lab0(q);
CTxR(t)		=	CTxR0;
KTxR(t)		=	KTxR0;


$INCLUDE bounds.gms

*++++++++++++++++++++++++++++++++++++++++++++++++++++++++
* Options
*++++++++++++++++++++++++++++++++++++++++++++++++++++++++

OPTIONS SOLPRINT=ON, LIMROW=0, LIMCOL=0, ITERLIM=1000, RESLIM=500;
OPTION NLP=CONOPT;

*++++++++++++++++++++++++++++++++++++++++++++++++++++++++
* Solve
*++++++++++++++++++++++++++++++++++++++++++++++++++++++++

SOLVE OLG MAXIMIZING OBJ USING NLP;

EXECUTE_UNLOAD 'results_base.gdx'
*EXECUTE '=gdx2xls results_base.gdx'
;


*++++++++++++++++++++++++++++++++++++++++++++++++++++++++
* Shocks - Demographic
*++++++++++++++++++++++++++++++++++++++++++++++++++++++++

* 1) Transitory increase in birth rate

$ontext

IF(FlDemog1 EQ 1,
NN(q,e,t) $(ORD(t) GT CARD(tp) AND ORD(t) LE (CARD(tp)+6))
                                =       1.01*NN(q,e,t);

LOOP(t $(ORD(t) GT CARD(tp)),
    PopQE(q,e,t+1,gf)           =       NN(q,e,t)*PopQE(q,e,t,gf)
);

LOOP(t,
    PopQE(q,e,t+1,g+1)          =       PopQE(q,e,t,g)
);

TPop(t) $(ORD(t) LE CARD(t))    =       SUM((q,g,e),PopQE(q,e,t,g));

LOOP(t $(ORD(t) LE CARD(t)),
    GPop(t)                     =       TPop(t+1)/TPop(t) )
;
LOOP(t $(ORD(t) GT CARD(tp)),
    Gov(t+1)			=	GPop(t)*Gov(t)
);

OPTION NLP=CONOPT;
SOLVE OLG MAXIMIZING OBJ USING NLP;

*CS(s,t) = SUM(g,Pop(t,g)*ConS.L(s,t,g));

EXECUTE_UNLOAD 'results_dem1.gdx'
*EXECUTE '=gdx2xls results_dem1.gdx'
);

$INCLUDE bounds.gms
NN(q,e,t)                           	=      	1;
PopQE(q,e,t,g) $(ORD(t) EQ CARD(g))  	=	PopQE0(q,e,g);

LOOP(t $(ORD(t) GT CARD(tp)),
    PopQE(q,e,t+1,g)                  	=       NN(q,e,t)*PopQE(q,e,t,g)
);

LOOP(tp,
    PopQE(q,e,tp,g)	                =       SUM(t $(ORD(t) EQ CARD(g)), PopQE(q,e,t,g))
);

TPop(t) $(ORD(t) LE CARD(t))    =       SUM((q,g,e), PopQE(q,e,t,g));

LOOP(t $(ORD(t) LE CARD(t)),
    GPop(t)                     =       TPop(t+1)/TPop(t)
);
Gov(t)				=	G0;
LOOP(t $(ORD(t) GT CARD(tp)),
    Gov(t+1)			=	GPop(t)*Gov(t)
);

OPTION NLP=CONOPT;
SOLVE OLG MAXIMIZING OBJ USING NLP;
$offtext
$ontext
IF(FlDemog2 EQ 1,
* 1) Baby Boom and Baby Bust
NN(q,e,t) $(ORD(t) GT CARD(tp) AND ORD(t) LE CARD(tp)+3)    = 1.02*NN(q,e,t);
NN(q,e,t) $(ORD(t) GT CARD(tp)+3 AND ORD(t) LE CARD(tp)+6)  = 0.98*NN(q,e,t);

LOOP(t $(ORD(t) GT CARD(tp)),
    PopQE(q,e,t+1,gf)           =       NN(q,e,t)*PopQE(q,e,t,gf)
);

LOOP(t,
    PopQE(q,e,t+1,g+1)          =       PopQE(q,e,t,g)
);

TPop(t) $(ORD(t) LE CARD(t))    =       SUM((q,g,e),PopQE(q,e,t,g));

LOOP(t $(ORD(t) LE CARD(t)),
    GPop(t)                     =       TPop(t+1)/TPop(t) )
;
LOOP(t $(ORD(t) GT CARD(tp)),
    Gov(t+1)			=	GPop(t)*Gov(t)
);

OPTION NLP=CONOPT;
SOLVE OLG MAXIMIZING OBJ USING NLP;

EXECUTE_UNLOAD 'results_dem2.gdx'
*EXECUTE '=gdx2xls results_dem2.gdx'
);

$INCLUDE bounds.gms
NN(q,e,t)                           	=      	1;
PopQE(q,e,t,g) $(ORD(t) EQ CARD(g))  	=	PopQE0(q,e,g);

LOOP(t $(ORD(t) GT CARD(tp)),
    PopQE(q,e,t+1,g)                  	=       NN(q,e,t)*PopQE(q,e,t,g)
);

LOOP(tp,
    PopQE(q,e,tp,g)	                =       SUM(t $(ORD(t) EQ CARD(g)), PopQE(q,e,t,g))
);

TPop(t) $(ORD(t) LE CARD(t))    =       SUM((q,g,e), PopQE(q,e,t,g));

LOOP(t $(ORD(t) LE CARD(t)),
    GPop(t)                     =       TPop(t+1)/TPop(t)
);
Gov(t)				=	G0;
LOOP(t $(ORD(t) GT CARD(tp)),
    Gov(t+1)			=	GPop(t)*Gov(t)
);

OPTION NLP=CONOPT;
SOLVE OLG MAXIMIZING OBJ USING NLP;


$offtext
* Prductivity Shock
IF(FlProdu EQ 1,
    A("s1",t)	$(ORD(t) GT CARD(tp)+1) = 1.01*A("s1",t);

    SOLVE OLG USING NLP MINIMIZING OBJ;       

*    CS(s,t) = SUM(g,Pop(t,g)*ConS.L(s,t,g));
    EXECUTE_UNLOAD 'results_prod.gdx'
*	EXECUTE '=gdx2xls results_prod.gdx'

);

