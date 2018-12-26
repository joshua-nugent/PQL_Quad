*Note: 'dsn' = number of replications of the data set;

options nonotes;

%macro gendata(sigbsq = 1, beta0 = -2, beta1 = 0.5, n = 40, p = 100, dsn = 5);
data sim;
  sigbsq = &sigbsq; n = &n; p = &p; dsn = &dsn;
  do beta0 = &beta0;
	  do beta1 = &beta1;
    	do dsn = 1 to &dsn;
        	do cluster = 1 to n;
            	x1 = (cluster lt (n+1)/2);
            	randint = normal(0) * sqrt(sigbsq);
            	do obs = 1 to p;
                	linpred = beta0 + beta1*x1 + randint;
                	expit = exp(linpred)/(1 + exp(linpred));
                	y = (uniform(0) lt expit);
                	output;
            	end;
        	end;
    	end;
  	end;
  end;
run;
%mend gendata;

/* for checking */
/* %gendata(dsn = 2, n=4, p = 2, beta0= %str(-2,-1), beta1 = %str(0.5,2)); */

%macro pql;
ods select none;
ods output parameterestimates = pqlestimates covparms = cpbin;
title "pql";
proc glimmix
	data = sim;
/*	pconv = .00001;*/
	/* Convergence is mabe too loose with this value? */
	/*  pconv = .00001; /* Convergence is good with this value! */
	/*	pconv = .000001; /*Convergence is poor with this value */
		nloptions
			maxiter = 200
			technique = DBLDOG;
by beta0 beta1 dsn;
class obs;
model y = x1 / solution dist=bin;
random int / subject = cluster;
run;
ods select all;
%mend pql;


%macro laplace;
ods select none; /* Note "ods select none" -- suppresses all output being printed */
ods output parameterestimates = laplaceestimates;
title "laplace";
proc glimmix data = sim method = laplace;
by beta0 beta1 dsn;
class obs;
model y = x1 /solution dist=bin;
random int / subject = cluster;
run;
ods select all;
%mend laplace;


%macro quad(qpoints = 4);
ods select none;
ods output parameterestimates = quadestimates&qpoints covparms = cpbin_quad&qpoints;
title "quad &qpoints";
proc glimmix data = sim method = quad(qpoints= &qpoints);
	by beta0 beta1 dsn;
	class obs;
	model y = x1 /solution dist=bin;
	random int / subject = cluster;
run;
ods select all;
%mend quad;

%macro runsim(sigbsq = 1, beta0 = %str(-2,-1), beta1 = %str(1.25,1.5), n = 40, p = 100, dsn = 5);
ods select none;

%gendata(sigbsq = &sigbsq, beta0 = &beta0, beta1 = &beta1, 
                                          n = &n, p = &p, dsn = &dsn);
%pql;   /*%laplace;*/
%quad(qpoints=4);
%icc;
%iccq;

title "sigbsq = &sigbsq, beta0 = &beta0, beta1 = &beta1, n = &n, p = &p, dsn = &dsn";
data results; 
set pqlestimates (in = pql)
/*    laplaceestimates (in = laplace)*/
    quadestimates4 (in = quad4);
/*if laplace then source = "laplace";*/
if pql then source = "pql";
if quad4 then source = "quad4";
%mend;


%macro icc;
ods select none;
ods output covparms = cp;
title "ICC";
proc glimmix
    data = sim;
    by beta0 beta1 dsn;
    class obs;
    model y = x1 / solution dist=normal;
    random int / subject = cluster;
run;
ods select all;
%mend icc;

%macro iccq;
ods select none;
ods output covparms = cp_quad4;
title "ICC Quad4";
proc glimmix
    data = sim method = quad(qpoints=4);
    by beta0 beta1 dsn;
    class obs;
    model y = x1 / solution dist=normal;
    random int / subject = cluster;
run;
ods select all;
%mend iccq;


/* Key step ... Run the whole thing */
%runsim(							
	dsn = 2000, p = 25, n = 100, sigbsq = 1,	/* n is number of clusters, p is people per cluster */
	beta0 = %str(-4.59511985,-3.891820298,-3.47609869,-2.944438979,-2.197224577,-1.386294361,-0.8472978604,-0.4054651081),
/*  beta0 = %str(-4.59511985,-3.891820298,-3.47609869,-2.944438979,-2.197224577,-1.386294361,-0.8472978604,-0.4054651081),*/
	beta1 = %str(-0.6931471806,-0.2876820725,-0.1053605157,0.0953101798)/*,0.2851789422,0.4054651081,0.6931471806)
/*	beta1 = %str(-0.6931471806,-0.2876820725,-0.1053605157,0.0953101798,0.2851789422,0.4054651081,0.6931471806)*/
);


proc means data = results (where = (effect = "x1"));
class source beta1 beta0;
var estimate;
run;


data results2;
set results;
Odds_Ratio = exp(estimate);
run;

proc means data = results2 (where = (effect = "x1"));
class source beta1 beta0;
var Odds_Ratio;
run;

proc means data = results (where = (effect = "x1"));
class source beta1 beta0;
var stderr;
run;

proc sort data = cp; by dsn beta0 beta1; run;

/*proc print data = cp; run;*/

data icc_pql_normal;
set cp;
by dsn beta0 beta1;
retain sb2;
if first.beta1 then sb2 = estimate;
icc_pql_normal = sb2 / (sb2 + estimate);
if last.beta1 then output;
run;

proc means data = icc_pql_normal;
class beta1 beta0;
var icc_pql_normal;
run;

/*proc print data = cp_quad4; run;*/
proc sort data = cp_quad4; by dsn beta0 beta1; run;

data icc_quad4_normal;
set cp_quad4;
by dsn beta0 beta1;
retain sb2;
if first.beta1 then sb2 = estimate;
icc_quad4_normal = sb2 / (sb2 + estimate);
if last.beta1 then output;
run;
proc means data = icc_quad4_normal;
class beta1 beta0;
var icc_quad4_normal;
run;

/*proc print data = icc; run;*/

data icc_bin_pi_sq_pql;
set cpbin (where = (covparm = "Intercept"));
icc_pi_sq_pql = (estimate) / (estimate + 3.29);
run;
data icc_bin_pi_sq_quad4;
set cpbin_quad4 (where = (covparm = "Intercept"));
icc_pi_sq_quad4 = (estimate) / (estimate + 3.29);
run;

proc means data = icc_bin_pi_sq_pql;
class beta1 beta0;
var icc_pi_sq_pql;
run;

proc means data = icc_bin_pi_sq_quad4;
class beta1 beta0;
var icc_pi_sq_quad4;
run;

/* Want the histograms? */
/*proc univariate data = results;
var stderr;
histogram stderr;
run;*/

/*proc univariate data = results (where = (effect = "x1"));
var estimate;
histogram estimate;
run;*/
