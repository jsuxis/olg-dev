* ========================================================
* 7 Sectors
* Imperfect Substitutability between Age groups
* Declarations File
* ========================================================

* Generations
SET
    g           Generations                     /g1 * g4/
    gf(g)       First Generation
    gw(g)       Working Generations             /g1 * g3/
    gr(g)       Retired Generations             /g4/
    gl(g)       Last Generation
;

        gf(g) = YES$(ORD(g) EQ 1);
        gl(g) = YES$(ORD(g) EQ CARD(g));
        ALIAS(g,gg);

* Sectors
SET
    s		Sectors		/s1*s7/
;
	ALIAS (S,SS),(SS,SSD),(S,SD);

* Sectoral Employment / Sectoral Affiliation
SET
    e		Sectoral Employment	/e1*e7/
;
	ALIAS (e,ee);

* Qualifications
SET
    q		Qualifications	/q1*q2/
;

	ALIAS (q,qq);
	
* Time
SET
    t           Total Time Horizon              /t1 * t50/
    tp(t)       Periods of Previously Born      /t1 * t3/
    tm(t)       Model Horizon
    tmf(t)      First Period of Model Horizon
;
        tm(t) = YES$(ORD(t) GT CARD(tp) AND ORD(t) LT CARD(t));
        tmf(t) = YES$(ORD(t) EQ CARD(tp)+1);

* Parameters for Data Input
PARAMETER
    C0          	Initial Aggregate Consumption
    I0          	Initial Investment
    Ldem0(s)    	Initial Labor Demand
    Input0(s,ss)	Initial Sectoral Input Matrix
    YK0(s)      	Initial Capital
    Y0(s)       	Initial Output
    Tpop0       	Initial Total Population
    TPopQ0(q)		Initial Population by Qualification
    TPopQE0(q,e)	Initial Population by Qualification and Sector		
    GPop0       	Initial Gross Population Growth Rate
    Pop0(g)     	Initial Population by Generation
    PopQ0(g,q)		Initial Population by Generation and Qualification
    PopQE0(q,e,g)	Initial Population by Generation Sector-Affil. and Qualification    
    PopQS0(q,s,g)	Initial Population by Generation Sector-Empl. and Qualification
    SAM(*,*)		Social Accounting Matrix
    CS0(s)		Initial Sectoral Consumption
    IS0(s)		Initial Sectoral Investment
    G0			Initial Government Spending
    GS0(s)		Initial Sectoral Government Spending
    GdpBondR		GDP to Bond Ratio (Public Debt Ratio)
    Bond0		Public Debt (Scaled due to long periods)
    CTxR0		Initial Consumption Tax Rate
    WTxR0		Initial Income Tax Rate
    KTxR0		Initial Capital Tax Rate
;

* Parameters for Calibration
PARAMETER
    A0(s)       	Initial Scaling Factor
    AlK(s)      	Capital Share (alpha)
    AlX(s)		Input Share
    X0(s)		Total Input
    sigma       	Sigma (Intertemporal Elasticity of Substitution)
    sigCon(g)		Elasticity of Substitution between consumption goods
    SigInter		Intertemporal Elasticity of Substitution
    SigIntra		Intratemporal Elasticity of Substitution
    sigAge(s)
    AlConS(s,g,q,e)	Share parameter for consumption
    AlConSC(s,g) 	Share of spending by sectors and generation
    eta0(s,q)		Initial Sector affinity
    eta(s,q,t)		Sector affinity
    rho0        	Initial Discount Rate
    B0(g,q,e)     	Initial Bond (Lend)
    Con0(g,q,e)   	Initial Consumption Demand
    ConS0(s,g,q,e)	Initial Sectoral Consumption Demand
    Beq0(g,q,e)   	Initial Bequest
    Inh0(g,q,e)   	Initial Inheritance
    BeqR(g)     	Bequest Rate
    InhR(g)     	Inheritance Rate
    EP(g)       	Earnings Profile (base)
    EPQ(g,q)		Earnings Profile for each qualification
    Leis0(g,q,e)	Initial Leisure
    LeisS0(g,q,s)	Initial Leisure
    Lab0(q,g,e)        	Initial Physical (non-effective) Labor Supply (hours)
    LabS0(q,s,g)	Sectoral Initial Physical (non-effective) Labor Supply (hours)
    Lsup0       	Initial Labor Supply
    LsupQ0(q)		Initial Labor Supply by qualification
    LsupEQ0(e,q,g) 	Initial Labor Supply by qualification and sector
    delta       	Depreciation Rate
    K0          	Initial Capital Stock
    Rent0       	Initial Rental Rate of Capital
    Rint0       	Initial Interest Rate
    Kdem0       	Initial Capital Demand (check)
    Ksc0        	Share of Capital Demand on total Stock of capital (Capital Share?)
    CTR0 		Initial Contribution Rate to Public Pension System
    SigInv		Elasticity of substitution between Investment goods
    AlInvS(s)		Share Parameter for Sectoral Investment
    LQ0(s,q)		Initial Labor Demand by Qualification
    LQA0(s,q,g)		Initial Labor Demand by Qualification and Age
    AlDemQ(s,q)		Share Parameter for labor demand by qualification
    sigLdem(s)		Elasticity of Substitution between labor types
    AlXS(s,ss)		Share Parameter for sectoral inputs
    sigX(s)		Elasticity of Substitution between inputs
    AlAge(s,q,g)		Share Parameter for Age groups in production
    SigAge(s)		Elasticity of Substitution between age groups
    Gamma(g,q,e)	Consumption Intensity Parameter
    VV0(g,q,e)		Initial V
    AlGovS(s)		Share Parameter for Sectoral Government Spending
    DIFFC0		Test Parameter
    TESTJ		Test Parameter
    TESTJG(g,q,e)	Test Parameter
    TESTDEM(s)		Test Parameter
;

* Skipped Parameters: Z0, PopGR0

* Postcalibration Parameters
PARAMETER
    LABEP(g)    Effective Units of Labor per generation
    SUMLABEP    Effective Units of Labor in total
;

* Variables for Calibration
VARIABLES
    Beq0v(g,q,e)
    Inh0v(g,q,e)
    Con0v(g,q,e)
    B0v(g,q,e)
    rho0v
    Rent0v
    Rint0v
    Obj
    ConS0v(s,g,q,e)
    LQ0v(s,q)
    gammav(g,q,e)
    VV0v(g,q,e)
    Bond0v
;

* Parameters for Simulation
PARAMETER
    A(s,t)              Scaling factor for C-D Production function
    NN(q,e,t)           Birth Rate
    GPop(t)             Gross Population Growth Rate
    TPop(t)             Total Population
    Pop(t,g)            Population by Generation
    PopQE(q,e,t,g)	Household size by Generation
    Lab(q,g,e,t)        Physical (non-effective) Labor Supply (hours)
    rho			Discount Rate
    rhoGen(t,g)         Generation Specific Discount Rate
    RkG(g)              ???
    TB(t)               Trade Balance
    WTB(t)              World(?) Trade Balance
    Walras(t)           Walras Condition
    WAssetBal(t)        Asset Balance Condition
    CS(s,t)		Total Sectoral Comsumption
    CTxR(t)		Consumption Tax Rate
    KTxR(t)		Capital Tax Rate (interest income)
    Gov(t)		Government Spending
    SigGov		Government Spending Elasticity
;

RkG(g)  =       ORD(g);

* Variables for Simulation
VARIABLE
    Y(s,t)      	Output
    P(s,t)		Price of Consumption goods
    Kdem(s,t)   	Capital Demand
    Ldem(s,t)   	Labor Demand
    LQ(s,q,t)		Labor Demand by Qualification
    LQA(s,q,g,t)	Labor Demand by Qualification and Age
    Xdem(s,t)		Sectoral Input Demand
    Input(s,ss,t)	Sectoral Input Matrix
    Con(q,e,t,g)	Consumption per generation
    PCon(q,e,t,g)	Price Index of total consumption per generation
    ConS(s,q,e,t,g)	Sectoral consumption per generations
    B(q,e,t,g)      	Bond (Lend)
    Beq(q,e,t,g)    	Bequests
    Inh(q,e,t,g)    	Inheritance
    Leis(q,e,t,g)	Leisure
    VV(q,e,t,g)		V
    mu(q,e,t,g)		Mu
    LsupQE(e,q,g,t)    	Labor Supply by sector and qualification
    I(t)        	Investments
    PI(t)		Price of Investment Good
    IS(s,t)		Investment using goods produced by sector s
    K(t)        	Capital Stock
    W(s,t)      	Composite wage index for effective units of labor
    PX(s,t)		Composite Input Demand Price Index
    wage(s,q,t)		Skill specific wage
    wageE(e,q,t)	Skill specific wage (earned)
    wageA(s,q,g,t)	Age specific wage
    wageAE(e,q,g,t)	Age specific wage (earned)
    Rent(t)     	Rental Rate of Capital
    Rint(t)     	Interest Rate
    Bond(t)		Government Debt
    WTxR(t)		Wage Tax Rate
    GovS(s,t)		Sectoral Government Spending
    PGov(t)		Price of Government Spending
;

EXECUTE_UNLOAD 'declarations.gdx';


