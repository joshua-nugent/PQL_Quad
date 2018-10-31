*Note: unclear on what 'dsn' variable is -
*number of replications of the data set?
*I set it to be 10 just for the test runs here
*but not sure if that's correct.;

%macro gendata(sigbsq = 1, beta0 = -2, beta1 = 0.5, n = 40, p = 100, dsn = 5);

data sim;
    sigbsq = &sigbsq; n = &n; p = &p; dsn = &dsn;
/*    beta0 = &beta0; beta1 = &beta1; */
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



/* Josh-- check the log-- note that some may fail to converge and we may need 
    to do something about this */
%macro pql;
ods select none;
ods output parameterestimates = pqlestimates;
title "pql";
proc glimmix data = sim;
by beta0 beta1 dsn;
class obs;
nloptions maxiter = 100;
model y = x1 / solution dist=bin;
random int / subject = cluster;
run;
ods select all;
%mend pql;

/*%pql;*/


/* josh-- note ods select none-- suppresses all output being printed */
%macro laplace;
ods select none;
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


%macro quad(qpoints = &qpoints);
ods select none;
ods output parameterestimates = quadestimates&qpoints;
title "quad &qpoints";
proc glimmix data = sim method = quad(qpoints= &qpoints);
by beta0 beta1 dsn;
class obs;
model y = x1 /solution dist=bin;
random int / subject = cluster;
run;
ods select all;
%mend quad;

%macro runsim(sigbsq = 1, beta0 = %str(-2,-1), beta1 = %str(0.5,0.75), n = 40, p = 100, dsn = 5);
%gendata(sigbsq = &sigbsq, beta0 = &beta0, beta1 = &beta1, 
                                          n = &n, p = &p, dsn = &dsn);
%pql;
%laplace;
%quad(qpoints=4);
/*%quad(qpoints=25);*/
title "sigbsq = &sigbsq, beta0 = &beta0, beta1 = &beta1, n = &n, p = &p, dsn = &dsn";
data results; 
set pqlestimates (in = pql)
    laplaceestimates (in = laplace)
    quadestimates4 (in = quad4);
/*    quadestimates25 (in = quad25);*/
if laplace then source = "laplace";
if pql then source = "pql";
if quad4 then source = "quad4";
if quad25 then source = "quad25";
%mend;

%runsim (dsn = 10, p = 75, n = 40, beta0=%str(-2,-1), beta1=%str(1.25,1.5));


/*proc print data = results (obs = 10); run;*/

/* added below to do the odds ratio */
data results2;
set results;
oddsratio = exp(estimate);
run;

proc means data = results2 (where = (effect = "x1"));
class  source beta0 beta1;
var estimate oddsratio;
run;




