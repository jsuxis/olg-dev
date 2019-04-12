*++++++++++++++++++++++++++++++++++++++++++++++++++++++++
* Initial Values and Bounds
*++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Y.L(s,t)        	=       Y0(s);
Y.FX(s,tp)		=	Y0(s);
Y.LO(s,t)		=	1.E-13;
VA.L(s,t)        	=       VA0(s);
VA.FX(s,tp)		=	VA0(s);
VA.LO(s,t)		=	1.E-13;
PVA.L(s,t)		=	1;
PVA.FX(s,tp)		=	1;
PVA.LO(s,t)		=	1.E-13;
P.L(s,t)		=	1;
P.FX("s1",t)		=	1;
P.FX(s,tp)		=	1;
Kdem.L(s,t)     	=       Kdem0(s);
Kdem.FX(s,tp)		=	Kdem0(s);
Kdem.LO(s,t)		=	1.E-13;
Ldem.L(s,t)     	=       Ldem0(s);
Ldem.FX(s,tp)		=	Ldem0(s);
Ldem.LO(s,t)		=	1.E-13;
LQ.L(s,q,t)		=	LQ0(s,q);
LQ.FX(s,q,tp)		=	LQ0(s,q);
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
ConS.LO(s,q,e,t,g)	=	0;
PCon.L(q,e,t,g)		=	1;
PCon.FX(q,e,tp,g)	=	1;
PCon.LO(q,e,t,g)	=	1.E-13;
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
K.L(t)          	=       K0;
K.LO(t)			=	0;		
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
LsupQE.L(e,q,t)		=	LsupEQ0(e,q);
LsupQE.FX(e,q,tp)	=	LsupEQ0(e,q);
Bond.L(t)		=	Bond0;
Bond.FX(tp)		=	Bond0;
Bond.FX(tmf)		=	Bond0;
WTxR.L(t)		=	WTxR0;
WTxR.FX(tp)		=	WTxR0;
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
