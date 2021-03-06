* ========================================================
* 7 Sectors
* Data File
* ========================================================

$INCLUDE declarations.gms

$CALL GDXXRW.EXE SAM.xlsx par=SAM rng=sheet1!A13:O28
$GDXIN SAM.gdx
$LOAD SAM
DISPLAY SAM;

CS0("s1")   	  	= SAM("firm1","consumption");
CS0("s2")    		= SAM("firm2","consumption");
CS0("s3")		= SAM("firm3","consumption");
CS0("s4")		= SAM("firm4","consumption");
CS0("s5")		= SAM("firm5","consumption");
CS0("s6")		= SAM("firm6","consumption");
CS0("s7")		= SAM("firm7","consumption");
IS0("s1")     	 	= SAM("firm1","investment");
IS0("s2")		= SAM("firm2","investment");
IS0("s3")		= SAM("firm3","investment");
IS0("s4")		= SAM("firm4","investment");
IS0("s5")		= SAM("firm5","investment");
IS0("s6")		= SAM("firm6","investment");
IS0("s7")		= SAM("firm7","investment");
Input0("s1","s1")	= SAM("firm1","firm1");
Input0("s1","s2")	= SAM("firm1","firm2");
Input0("s1","s3")	= SAM("firm1","firm3");
Input0("s1","s4")	= SAM("firm1","firm4");
Input0("s1","s5")	= SAM("firm1","firm5");
Input0("s1","s6")	= SAM("firm1","firm6");
Input0("s1","s7")	= SAM("firm1","firm7");
Input0("s2","s1")	= SAM("firm2","firm1");
Input0("s2","s2")	= SAM("firm2","firm2");
Input0("s2","s3")	= SAM("firm2","firm3");
Input0("s2","s4")	= SAM("firm2","firm4");
Input0("s2","s5")	= SAM("firm2","firm5");
Input0("s2","s6")	= SAM("firm2","firm6");
Input0("s2","s7")	= SAM("firm2","firm7");
Input0("s3","s1")	= SAM("firm3","firm1");
Input0("s3","s2")	= SAM("firm3","firm2");
Input0("s3","s3")	= SAM("firm3","firm3");
Input0("s3","s4")	= SAM("firm3","firm4");
Input0("s3","s5")	= SAM("firm3","firm5");
Input0("s3","s6")	= SAM("firm3","firm6");
Input0("s3","s7")	= SAM("firm3","firm7");
Input0("s4","s1")	= SAM("firm4","firm1");
Input0("s4","s2")	= SAM("firm4","firm2");
Input0("s4","s3")	= SAM("firm4","firm3");
Input0("s4","s4")	= SAM("firm4","firm4");
Input0("s4","s5")	= SAM("firm4","firm5");
Input0("s4","s6")	= SAM("firm4","firm6");
Input0("s4","s7")	= SAM("firm4","firm7");
Input0("s5","s1")	= SAM("firm5","firm1");
Input0("s5","s2")	= SAM("firm5","firm2");
Input0("s5","s3")	= SAM("firm5","firm3");
Input0("s5","s4")	= SAM("firm5","firm4");
Input0("s5","s5")	= SAM("firm5","firm5");
Input0("s5","s6")	= SAM("firm5","firm6");
Input0("s5","s7")	= SAM("firm5","firm7");
Input0("s6","s1")	= SAM("firm6","firm1");
Input0("s6","s2")	= SAM("firm6","firm2");
Input0("s6","s3")	= SAM("firm6","firm3");
Input0("s6","s4")	= SAM("firm6","firm4");
Input0("s6","s5")	= SAM("firm6","firm5");
Input0("s6","s6")	= SAM("firm6","firm6");
Input0("s6","s7")	= SAM("firm6","firm7");
Input0("s7","s1")	= SAM("firm7","firm1");
Input0("s7","s2")	= SAM("firm7","firm2");
Input0("s7","s3")	= SAM("firm7","firm3");
Input0("s7","s4")	= SAM("firm7","firm4");
Input0("s7","s5")	= SAM("firm7","firm5");
Input0("s7","s6")	= SAM("firm7","firm6");
Input0("s7","s7")	= SAM("firm7","firm7");
X0(s)			= SUM(ss,Input0(ss,s));	
YK0("s1")	 	= SAM("capital","firm1");
YK0("s2")	 	= SAM("capital","firm2");
YK0("s3")		= SAM("capital","firm3");
YK0("s4")		= SAM("capital","firm4");
YK0("s5")		= SAM("capital","firm5");
YK0("s6")		= SAM("capital","firm6");
YK0("s7")		= SAM("capital","firm7");
LQ0("s1","q1")		= SAM("labor1","firm1");
LQ0("s1","q2")		= SAM("labor2","firm1");
LQ0("s2","q1")		= SAM("labor1","firm2");
LQ0("s2","q2")		= SAM("labor2","firm2");
LQ0("s3","q1")		= SAM("labor1","firm3");
LQ0("s3","q2")		= SAM("labor2","firm3");
LQ0("s4","q1")		= SAM("labor1","firm4");
LQ0("s4","q2")		= SAM("labor2","firm4");
LQ0("s5","q1")		= SAM("labor1","firm5");
LQ0("s5","q2")		= SAM("labor2","firm5");
LQ0("s6","q1")		= SAM("labor1","firm6");
LQ0("s6","q2")		= SAM("labor2","firm6");
LQ0("s7","q1")		= SAM("labor1","firm7");
LQ0("s7","q2")		= SAM("labor2","firm7");
Ldem0("s1")		= SUM(q,LQ0("s1",q));
Ldem0("s2")		= SUM(q,LQ0("s2",q));
Ldem0("s3")		= SUM(q,LQ0("s3",q));
Ldem0("s4")		= SUM(q,LQ0("s4",q));
Ldem0("s5")		= SUM(q,LQ0("s5",q));
Ldem0("s6")		= SUM(q,LQ0("s6",q));
Ldem0("s7")		= SUM(q,LQ0("s7",q));
Y0(s)    		= Ldem0(s) + YK0(s) + X0(s);
GS0("s1")		= SAM("firm1","government");
GS0("s2")		= SAM("firm2","government");
GS0("s3")		= SAM("firm3","government");
GS0("s4")		= SAM("firm4","government");
GS0("s5")		= SAM("firm5","government");
GS0("s6")		= SAM("firm6","government");
GS0("s7")		= SAM("firm7","government");
G0			= SUM(s, GS0(s));
CTxR0			= SAM("tax","consumption")/SUM(s,CS0(s));
KTxR0 			= SAM("tax","capital")/SUM(s,YK0(s));
WTxR0			= SAM("tax","labor1")/SUM(s,LQ0(s,"q1"));


GdpBondR		= 0.5;
Bond0			= GdpBondR*SUM(s,Y0(s))/20;
TPop0   		= 150;
TPopQ0("q1")		= TPop0*0.4;
TPopQ0("q2")		= Tpop0*0.6;
Gpop0   		= 1;
Pop0(g) 		= Tpop0/CARD(g);
PopQ0(q,g)		= TPopQ0(q)/CARD(g);


DISPLAY CTxR0,KTxR0,WTxR0;
EXECUTE_UNLOAD 'data.gdx';
