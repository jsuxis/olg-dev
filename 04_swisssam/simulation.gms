* ========================================================
* 7 Sectors
* Imperfect Substitutability between Age groups
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
FlDemog1	= 	1;
FlDemog2     	= 	0;
FlProdu     	= 	0;
FlQual		=	0;

*++++++++++++++++++++++++++++++++++++++++++++++++++++++++
* Equations
*++++++++++++++++++++++++++++++++++++++++++++++++++++++++

EQUATIONS
    QEq(s,t)            Production Function
    KdemEq(s,t)         Capital Demand
    LdemEq(s,t)         Labor Demand
    XdemEq(s,t)		Input Demand
    HBudgEq1(q,e,t,g)   HH Budget Constraint
    HBudgEq2(q,e,t,g)   HH Budget Constraint last generation
    BeqEq(q,e,t,g) 	Bequests
    InhEq(q,e,t,g)      Inheritance
    ConEq(q,e,t,g)      Euler Equation
    ConEqSS(q,e,t,g)    Steady State Consumption profile
    LeisEq(q,e,t,g)	Intratemporal FOC
    MuEq(q,e,t,g)	Mu is only positive if not working (Kuhn-Tucker)
    VVEq(q,e,t,g)	Definition of V
    LsupEq(e,q,g,t)       Labor Supply Balance
    RintEq(t)           Balance of interest and rental rates
    PConEq(q,e,t,g)	FOC of Intratemporal (sectoral) consumption allocation
    ConSEq(q,e,s,t,g)	FOC of Intratemporal (sectoral) consumption allocation
    PInvEq(t)		FOC of Investment origin allocation
    InvSEq(s,t)		FOC of Investment origin allocation
    WLdemEq(s,t)	FOC of qualification specific labor demand allocation
    LdemQEq(s,q,t)	FOC of qualification specific labor demand allocation
    WAdemEq(s,q,t)	FOC of age specific labor demand allocation
    AdemQEq(s,q,g,t)	FOC of age specific labor demand allocation
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
    WageEq(s,q,g,t)       Labor Market Equilibrium
    WageEq2(s,q,g,t)
    WageEq3(s,q,g,t)      
    WageEq4(s,q,g,t)
    WageEq5(s,q,g,t)
    WageEq6(s,q,g,t)
    WageEq7(s,q,g,t)
    WageEq8(s,q,g,t)	
    Wage2Eq1(s,q,g,t)
    Wage2Eq2(s,q,g,t)
    Wage2Eq3(s,q,g,t)
    Wage2Eq4(s,q,g,t)
    Wage2Eq5(s,q,g,t)
    Wage2Eq6(s,q,g,t)
    Wage2Eq7(s,q,g,t)
    Wage2Eq8(s,q,g,t)
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
    LOG(W(s,t)) + LOG(Ldem(s,t))
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
HBudgEq1(q,e,t+1,g+1)	$(ORD(t) GT CARD(tp) AND ORD(t) LT CARD(tp) + CARD(tm))..
    PCon(q,e,t,g)*(1+CTxR(t))*Con(q,e,t,g) + B(q,e,t+1,g+1)
    =E=
    ((1-WTxR(q,t))*wageAE(e,q,g,t)*Lab(q,g,e,t)*EPQ(g,q))+
    ((Rint(t)-KTxR(t)*(Rint(t)-1))*B(q,e,t,g))$(ORD(g) GT 1) +
    Inh(q,e,t,g) - Beq(q,e,t,g)
    ;

* HH Budget Constraint last generation
HBudgEq2(q,e,t,gl)	$(ORD(t) GT CARD(tp) AND ORD(t) LE CARD(tp)+CARD(tm))..
    PCon(q,e,t,gl)*(1+CTxR(t))*Con(q,e,t,gl)
    =E=
    ((1-WTxR(q,t))*wageAE(e,q,gl,t)*Lab(q,gl,e,t)*EPQ(gl,q)) +
    ((Rint(t)-KTxR(t)*(Rint(t)-1))*B(q,e,t,gl)) +
    Inh(q,e,t,gl) - Beq(q,e,t,gl)
    ;
* Bequests
BeqEq(q,e,t,g)      $(ORD(t) GT CARD(tp) AND ORD(t) LE CARD(tp)+CARD(tm) AND BeqR(g) NE 0)..
    Beq(q,e,t,g)
    =E=
    BeqR(g)*PCon(q,e,t,g)*Con(q,e,t,g)
    ;

* Inheritance
InhEq(q,e,t,g)     $(ORD(t) GT CARD(tp) AND ORD(t) LE CARD(tp)+CARD(tm) AND InhR(g) NE 0)..
    PopQE(q,e,t,g)*Inh(q,e,t,g)
    =E=
    InhR(g)*SUM(gr, PopQE(q,e,t,g)*Beq(q,e,t,g))
    ;

* Euler Equation
ConEq(q,e,t+1,g+1)  $(ORD(t) GT CARD(tp) AND ORD(t) LE CARD(tp)+CARD(tm)-1)..
    Con(q,e,t+1,g+1)/Con(q,e,t,g)
    =E=
    ((((Rint(t+1)-KTxR(t+1)*(Rint(t+1)-1))*PCon(q,e,t,g) / (rho*PCon(q,e,t+1,g+1)))*((1+CTxR(t))/(1+CTxR(t+1))))**SigInter)
    ;

* Steady State Consumption Profile
ConEqSS(q,e,t+1,g)  $((ORD(t) EQ (CARD(tp) + CARD(tm)-1) AND ORD(g) LT CARD(g)))..
    Con(q,e,t+1,g)
    =E=
    Con(q,e,t,g)
    ;

* FOC of Intratemporal (sectoral) consumption allocation
PConEq(q,e,t,g)	$(ORD(t) GT CARD(tp) AND ORD(t) LE CARD(tp)+CARD(tm))..
    PCon(q,e,t,g)**(1-SigCon(g))
    =E=
    SUM(s$(AlConS(s,g,q,e) GT 1.E-13), AlConS(s,g,q,e) * P(s,t)**(1-SigCon(g)))
    ;

* FOC of Intratemporal (sectoral) consumption allocation
ConSEq(q,e,s,t,g)	$(ORD(t) GT CARD(tp) AND ORD(t) LE CARD(tp)+CARD(tm)
		AND AlConS(s,g,q,e) GT 1.E-13)..
    ConS(s,q,e,t,g)
    =E=
    AlConS(s,g,q,e)*(PCon(q,e,t,g)/P(s,t))**SigCon(g)*Con(q,e,t,g)
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
    W(s,t)**(1-sigLdem(s))
    =E=
    SUM(q$(AlDemQ(s,q) GT 1.E-13), AlDemQ(s,q) * wage(s,q,t)**(1-sigLdem(s)))
    ;

* FOC of qualification specific labor demand allocation
LdemQEq(s,q,t)	$(ORD(t) GT CARD(tp) AND ORD(t) LE CARD(tp)+CARD(tm)
		AND AlDemQ(s,q) GT 1.E-13)..
    LQ(s,q,t)
    =E=
    AlDemQ(s,q)*(W(s,t)/wage(s,q,t))**sigLdem(s)*Ldem(s,t)
    ;

* FOC of age specific labor demand allocation
WAdemEq(s,q,t)	$(ORD(t) GT CARD(tp) AND ORD(t) LE CARD(tp)+CARD(tm))..
    wage(s,q,t)**(1-sigAge(s))
    =E=
    SUM(g$(AlAge(s,q,g) GT 1.E-13), AlAge(s,q,g) * wageA(s,q,g,t)**(1-sigAge(s)))
    ;

* FOC of age specific labor demand allocation
AdemQEq(s,q,g,t)	$(ORD(t) GT CARD(tp) AND ORD(t) LE CARD(tp)+CARD(tm) AND AlAge(s,q,g) GT 1.E-13)..
    LQA(s,q,g,t)
    =E=
    AlAge(s,q,g)*(wage(s,q,t)/wageA(s,q,g,t))**sigAge(s)*LQ(s,q,t)
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
LsupEq(e,q,g,t)       $(ORD(t) GT CARD(tp) AND ORD(t) LE CARD(tp)+CARD(tm))..
    LsupQE(e,q,g,t)
    =E=
    PopQE(q,e,t,g)*Lab(q,g,e,t)*EPQ(g,q)
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
    PGov(t)*Bond(t+1)+(SUM((q,e,g),PopQE(q,e,t,g)*(
      WTxR(q,t)*wageAE(e,q,g,t)*Lab(q,g,e,t)*EPQ(g,q)
    + CTxR(t)*PCon(q,e,t,g)*Con(q,e,t,g)
    + KTxR(t)*(Rint(t)-1)*B(q,e,t,g)))) 
    =E=
    PGov(t)*Gov(t)+(Rint(t-1)*(PGov(t)/PGov(t-1)))*PGov(t-1)*Bond(t)
    ;

* Government Budget in Steady State
GBudgEqSS(t)	$(ORD(t) EQ CARD(tp)+CARD(tm))..
    GPop(t)*PGov(t)*Bond(t)+(SUM((q,e,g),PopQE(q,e,t,g)*(
      WTxR(q,t)*wageAE(e,q,g,t)*Lab(q,g,e,t)*EPQ(g,q)
    + CTxR(t)*PCon(q,e,t,g)*Con(q,e,t,g)
    + KTxR(t)*(Rint(t)-1)*B(q,e,t,g))))
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
    SUM((q,e,g), PopQE(q,e,t,g)*ConS(s,q,e,t,g)) + IS(s,t) + GovS(s,t)
    + SUM(ss,Input(s,ss,t))
    ;

* Labor market equilibrium
WageEq(s,q,g,t)       $(ORD(t) GT CARD(tp) AND ORD(t) LE CARD(tp)+CARD(tm))..
    LsupQE("e1",q,g,t)
    =E=
    LQA("s1",q,g,t)
    ;

* Labor market equilibrium
WageEq2(s,q,g,t)       $(ORD(t) GT CARD(tp) AND ORD(t) LE CARD(tp)+CARD(tm))..
    LsupQE("e2",q,g,t)
    =E=
    LQA("s2",q,g,t)
    ;

* Labor market equilibrium
WageEq3(s,q,g,t)       $(ORD(t) GT CARD(tp) AND ORD(t) LE CARD(tp)+CARD(tm))..
    LsupQE("e3",q,g,t)
    =E=
    LQA("s3",q,g,t)
    ;

* Labor market equilibrium
WageEq4(s,q,g,t)       $(ORD(t) GT CARD(tp) AND ORD(t) LE CARD(tp)+CARD(tm))..
    LsupQE("e4",q,g,t)
    =E=
    LQA("s4",q,g,t)
    ;

* Labor market equilibrium
WageEq5(s,q,g,t)       $(ORD(t) GT CARD(tp) AND ORD(t) LE CARD(tp)+CARD(tm))..
    LsupQE("e5",q,g,t)
    =E=
    LQA("s5",q,g,t)
    ;

* Labor market equilibrium
WageEq6(s,q,g,t)       $(ORD(t) GT CARD(tp) AND ORD(t) LE CARD(tp)+CARD(tm))..
    LsupQE("e6",q,g,t)
    =E=
    LQA("s6",q,g,t)
    ;

* Labor market equilibrium
WageEq7(s,q,g,t)       $(ORD(t) GT CARD(tp) AND ORD(t) LE CARD(tp)+CARD(tm))..
    LsupQE("e7",q,g,t)
    =E=
    LQA("s7",q,g,t)
    ;

* Labor market equilibrium
WageEq8(s,q,g,t)       $(ORD(t) GT CARD(tp) AND ORD(t) LE CARD(tp)+CARD(tm))..
    LsupQE("e8",q,g,t)
    =E=
    LQA("s8",q,g,t)
    ;

* Labor market equilibrium
Wage2Eq1(s,q,g,t)       $(ORD(t) GT CARD(tp) AND ORD(t) LE CARD(tp)+CARD(tm))..
    wageAE("e1",q,g,t)
    =E=
    wageA("s1",q,g,t)
    ;

* Labor market equilibrium
Wage2Eq2(s,q,g,t)       $(ORD(t) GT CARD(tp) AND ORD(t) LE CARD(tp)+CARD(tm))..
    wageAE("e2",q,g,t)
    =E=
    wageA("s2",q,g,t)
    ;

* Labor market equilibrium
Wage2Eq3(s,q,g,t)       $(ORD(t) GT CARD(tp) AND ORD(t) LE CARD(tp)+CARD(tm))..
    wageAE("e3",q,g,t)
    =E=
    wageA("s3",q,g,t)
    ;
    
* Labor market equilibrium
Wage2Eq4(s,q,g,t)       $(ORD(t) GT CARD(tp) AND ORD(t) LE CARD(tp)+CARD(tm))..
    wageAE("e4",q,g,t)
    =E=
    wageA("s4",q,g,t)
    ;

* Labor market equilibrium
Wage2Eq5(s,q,g,t)       $(ORD(t) GT CARD(tp) AND ORD(t) LE CARD(tp)+CARD(tm))..
    wageAE("e5",q,g,t)
    =E=
    wageA("s5",q,g,t)
    ;

* Labor market equilibrium
Wage2Eq6(s,q,g,t)       $(ORD(t) GT CARD(tp) AND ORD(t) LE CARD(tp)+CARD(tm))..
    wageAE("e6",q,g,t)
    =E=
    wageA("s6",q,g,t)
    ;

* Labor market equilibrium
Wage2Eq7(s,q,g,t)       $(ORD(t) GT CARD(tp) AND ORD(t) LE CARD(tp)+CARD(tm))..
    wageAE("e7",q,g,t)
    =E=
    wageA("s7",q,g,t)
    ;

* Labor market equilibrium
Wage2Eq8(s,q,g,t)       $(ORD(t) GT CARD(tp) AND ORD(t) LE CARD(tp)+CARD(tm))..
    wageAE("e8",q,g,t)
    =E=
    wageA("s8",q,g,t)
    ;
    
* Capital Market Equilibrium
RentEq(t)       $(ORD(t) GT CARD(tp) AND ORD(t) LE CARD(tp)+CARD(tm))..
    Ksc0*K(t)
    =E=
    SUM(s,Kdem(s,t))
    ;

* Asset Market Equilibrium
 AssetEq(t+1)    $(ORD(t) GT CARD(tp) AND ORD(t) LE CARD(tp) + CARD(tm)-1 AND (FlWalras EQ 0))..
    SUM((g,q,e), PopQE(q,e,t+1,g+1)*B(q,e,t+1,g+1))
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
      WAdemEq
      AdemQEq
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
      WageEq2
      WageEq3
      WageEq4
      WageEq5
      WageEq6
      WageEq7
      WageEq8
      Wage2Eq1
      Wage2Eq2
      Wage2Eq3
      Wage2Eq4
      Wage2Eq5
      Wage2Eq6
      Wage2Eq7
      Wage2Eq8
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

* Parameters
A(s,t)          =       A0(s);
rho		=       rho0;
Lab(q,g,e,t)    =       Lab0(q,g,e);
CTxR(t)		=	CTxR0;
KTxR(t)		=	KTxR0;


*++++++++++++++++++++++++++++++++++++++++++++++++++++++++
* Initial Values and Bounds
*++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Y.L(s,t)        	=       Y0(s);
Y.FX(s,tp)		=	Y0(s);
Y.LO(s,t)		=	1.E-13;
P.L(s,t)		=	1;
P.FX("s1",t)		=	1;
P.FX(s,tp)		=	1;
P.LO(s,t)		=	1.E-13;
Kdem.L(s,t)     	=       Kdem0(s);
Kdem.FX(s,tp)		=	Kdem0(s);
Kdem.LO(s,t)		=	1.E-13;
Ldem.L(s,t)     	=       Ldem0(s);
Ldem.FX(s,tp)		=	Ldem0(s);
Ldem.LO(s,t)		=	1.E-13;
LQ.L(s,q,t)		=	LQ0(s,q);
LQ.FX(s,q,tp)		=	LQ0(s,q);
LQA.L(s,q,g,t)		=	LQA0(s,q,g);
LQA.FX(s,q,g,tp)	=	LQA0(s,q,g);
LQA.LO(s,q,g,t)		=	1.E-13;
B.L(q,e,t,g)      	=       B0(g,q,e);
B.LO(q,e,t,g)		=	0;
B.FX(q,e,tmf,g)   	=       B0(g,q,e);
B.FX(q,e,t,gf)		=	0;
B.FX(q,e,tp,g)		=	B0(g,q,e);
Beq.FX(q,e,t,g)   	=       Beq0(g,q,e);
Inh.FX(q,e,t,g)   	=       Inh0(g,q,e);
Leis.L(q,e,t,g)		=	Leis0(g,q,e);
Leis.LO(q,e,t,g)	=	0;
Leis.UP(q,e,t,g)	=	1;
Leis.FX(q,e,tp,g)	=	Leis0(g,q,e);
Mu.L(q,e,t,g)		=	0;
Mu.UP(q,e,t,g)		=	10;
Mu.LO(q,e,t,g)		=	0;
VV.L(q,e,t,g)		=	VV0(g,q,e);
Con.L(q,e,t,g)    	=       Con0(g,q,e);
Con.LO(q,e,t,g)		=	1.E-13;
Con.FX(q,e,tp,g)	=	Con0(g,q,e);
ConS.L(s,q,e,t,g)	=	ConS0(s,g,q,e);
ConS.FX(s,q,e,tp,g)	=	ConS0(s,g,q,e);
PCon.L(q,e,t,g)		=	1;
PCon.FX(q,e,tp,g)	=	1;
I.L(t)          	=       I0;
I.FX(tp)		=	I0;
I.LO(t)			=	0;
PI.L(t)			=	1;
PI.FX(tp)		=	1;
PI.LO(t)		=	0;
PGov.L(t)		=	1;
PGov.FX(tp)		=	1;
PGov.LO(t)		=	0;
GovS.L(s,t)		=	GS0(s);
GovS.FX(s,tp)		=	GS0(s);
IS.L(s,t)		=	IS0(s);
IS.FX(s,tp)		=	IS0(s);
IS.LO(s,t)		=	0;
K.L(t)          	=       K0;
K.LO(t)			=	1.E-13;
K.FX(tp)       		=       K0;
K.FX(tmf)       	=       K0;
Rent.L(t)       	=       Rent0;
Rent.FX(tp)     	=       Rent0;
Rint.L(t)       	=       Rint0;
Rint.FX(tp)     	=       Rint0;
W.L(s,t)          	=       1;
W.FX(s,tp)        	=       1;
wage.L(s,q,t)		=	1;
wage.LO(s,q,t)		=	1.E-13;
wage.FX(s,q,tp)		=	1;
wageE.L(e,q,t)		=	1;
wageE.LO(e,q,t)		=	1.E-13;
wageE.FX(e,q,tp)	=	1;
wageA.L(s,q,g,t)	=	1;
wageA.LO(s,q,g,t)	=	1.E-13;
wageA.FX(s,q,g,tp)	=	1;
wageAE.L(e,q,g,t)	=	1;
wageAE.LO(e,q,g,t)	=	1.E-13;
wageAE.FX(e,q,g,tp)	=	1;
LsupQE.L(e,q,g,t)	=	LsupEQ0(e,q,g);
LsupQE.FX(e,q,g,tp)	=	LsupEQ0(e,q,g);
LsupQE.LO(e,q,g,t)	=	1.E-13;
Bond.L(t)		=	Bond0;
Bond.FX(tp)		=	Bond0;
Bond.FX(tmf)		=	Bond0;
WTxR.L(q,t)		=	WTxR0(q);
WTxR.FX(q,tp)		=	WTxR0(q);
Input.L(s,ss,t)		=	Input0(s,ss);
Input.FX(s,ss,tp)	=	Input0(s,ss);
Input.LO(s,ss,t)	=	0;
PX.L(s,t)		=	1;
PX.FX(s,tp)		=	1;
PX.LO(s,t)		=	0;
Xdem.L(s,t)		=	X0(s);
Xdem.FX(s,tp)		=	X0(s);
Xdem.LO(s,t)		=	1.E-13;
Obj.FX          	=       0;




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
EXECUTE '=gdx2xls results_base.gdx'
;


*++++++++++++++++++++++++++++++++++++++++++++++++++++++++
* Shocks - Demographic
*++++++++++++++++++++++++++++++++++++++++++++++++++++++++

* 1) Transitory increase in birth rate



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
EXECUTE '=gdx2xls results_dem1.gdx'
);

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
EXECUTE '=gdx2xls results_dem2.gdx'
);

* Prductivity Shock
IF(FlProdu EQ 1,
    A("s1",t)	$(ORD(t) GT CARD(tp) AND ORD(t) LE CARD(tp)+6) = 1.01*A("s1",t);


    SOLVE OLG USING NLP MINIMIZING OBJ;       

*    CS(s,t) = SUM(g,Pop(t,g)*ConS.L(s,t,g));
    EXECUTE_UNLOAD 'results_prod.gdx'
	EXECUTE '=gdx2xls results_prod.gdx'

);


* Increase in the share of High Skill Workers
IF(FlQual EQ 1,

NN("q1",e,t) $(ORD(t) GT CARD(tp) AND ORD(t) LE (CARD(tp)+4))
                                =       1.02*NN("q1",e,t);
NN("q2",e,t) $(ORD(t) GT CARD(tp) AND ORD(t) LE (CARD(tp)+4))
                                =       0.98*NN("q2",e,t);				

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
OPTION NLP=CONOPT;
SOLVE OLG MAXIMIZING OBJ USING NLP;

EXECUTE_UNLOAD 'results_qual.gdx'
EXECUTE '=gdx2xls results_qual.gdx'
);
