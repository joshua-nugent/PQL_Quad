*Note: unclear on what 'dsn' variable is -
*number of replications of the data set?
*I set it to be 10 just for the test runs here
*but not sure if that's correct.;

%macro gendata(sigbsq = 1, beta0 = -2, beta1 = 0.5, n = 40, p = 100, dsn = 5);

data sim;
    sigbsq = &sigbsq; beta0 = &beta0; beta1 = &beta1; n = &n; p = &p; dsn = &dsn;
    do dsn = 1 to dsn;
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
run;
%mend gendata;

%gendata(dsn = 20, n=4, p = 10);

/* Josh-- check the log-- note that some may fail to converge and we may need 
    to do something about this */
%macro pql;
ods select none;
ods output parameterestimates = pqlestimates;
title "pql";
proc glimmix data = sim;
by dsn;
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
by dsn;
class obs;
model y = x1 /solution dist=bin;
random int / subject = cluster;
run;
ods select all;
%mend laplace;


%macro quad(qpoints = 25);
ods select none;
ods output parameterestimates = quadestimates&qpoints;
title "quad 25";
proc glimmix data = sim method = quad(qpoints= &qpoints);
by dsn;
class obs;
model y = x1 /solution dist=bin;
random int / subject = cluster;
run;
ods select all;
%mend quad;

%macro runsim(sigbsq = 1, beta0 = -2, beta1 = 0.5, n = 40, p = 100, dsn = 5);
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

%runsim (dsn = 100, p = 80, n = 30);

/*proc print data = results (obs = 10); run;*/

proc means data = results (where = (effect = "x1"));
class  source;
var estimate;
run;




