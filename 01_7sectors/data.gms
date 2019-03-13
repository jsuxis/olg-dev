
* =================================================================================
* October 2017   --- Simple Overlapping Generations Model ---
*
* This version
*        - 3 generations
*        - 1 region
*        - 1 sector (S)
*        - 2 factors of production (capital and labour)
*        - Dynamics vs steady state effects

* Data file
* =================================================================================

$INCLUDE declarations.gms

$CALL GDXXRW.EXE SAM.xlsx par=SAM rng=sheet1!A13:K24
$GDXIN SAM.gdx
$LOAD SAM
DISPLAY SAM;

CS0("s1")   	  	= SAM("firm1","consumption");
CS0("s2")    		= SAM("firm2","consumption");
CS0("s3")		= SAM("firm3","consumption");
IS0("s1")     	 	= SAM("firm1","investment");
IS0("s2")		= SAM("firm2","investment");
IS0("s3")		= SAM("firm3","investment");
Input0("s1","s1")	= SAM("firm1","firm1");
Input0("s1","s2")	= SAM("firm1","firm2");
Input0("s1","s3")	= SAM("firm1","firm3");
Input0("s2","s1")	= SAM("firm2","firm1");
Input0("s2","s2")	= SAM("firm2","firm2");
Input0("s2","s3")	= SAM("firm2","firm3");
Input0("s3","s1")	= SAM("firm3","firm1");
Input0("s3","s2")	= SAM("firm3","firm2");
Input0("s3","s3")	= SAM("firm3","firm3");
X0(s)			= SUM(ss,Input0(ss,s));	
YK0("s1")	 	= SAM("capital","firm1");
YK0("s2")	 	= SAM("capital","firm2");
YK0("s3")		= SAM("capital","firm3");
LQ0("s1","q1")		= SAM("labor1","firm1");
LQ0("s1","q2")		= SAM("labor2","firm1");
LQ0("s2","q1")		= SAM("labor1","firm2");
LQ0("s2","q2")		= SAM("labor2","firm2");
LQ0("s3","q1")		= SAM("labor1","firm3");
LQ0("s3","q2")		= SAM("labor2","firm3");
Ldem0("s1")		= SUM(q,LQ0("s1",q));
Ldem0("s2")		= SUM(q,LQ0("s2",q));
Ldem0("s3")		= SUM(q,LQ0("s3",q));
Y0(s)    		= Ldem0(s) + YK0(s) + X0(s);
GS0("s1")		= SAM("firm1","government");
GS0("s2")		= SAM("firm2","government");
GS0("s3")		= SAM("firm3","government");
G0			= SUM(s, GS0(s));
CTxR0			= SAM("tax","consumption")/SUM(s,CS0(s));
KTxR0 			= SAM("tax","capital")/SUM(s,YK0(s));
WTxR0			= SAM("tax","labor1")/SUM(s,LQ0(s,"q1"));


GdpBondR		= 0.5;
Bond0			= GdpBondR*SUM(s,Y0(s))/20;
TPop0   		= 150;
TPopQ0("q1")		= TPop0*0.4;
TPopQ0("q2")		= Tpop0*0.6;
* Assuming all sector-populations are the same size, can be adjusted in the final model:
TPopQE0(q,e)		= TpopQ0(q)/CARD(e);
Gpop0   		= 1;
Pop0(g) 		= Tpop0/CARD(g);
PopQ0(g,q)		= TPopQ0(q)/CARD(g);
PopQE0(q,e,g)		= TPopQE0(q,e)/CARD(g);
PopQS0(q,"s1",g)	= PopQE0(q,"e1",g);
PopQS0(q,"s2",g)	= PopQE0(q,"e2",g);
PopQS0(q,"s3",g)	= PopQE0(q,"e3",g);


DISPLAY CTxR0,KTxR0,WTxR0;
EXECUTE_UNLOAD 'data.gdx';
