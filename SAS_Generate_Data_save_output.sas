*Note: unclear on what 'dsn' variable is -
*number of replications of the data set?
*I set it to be 10 just for the test runs here
*but not sure if that's correct.;


data sim;
    sigbsq = 1; beta0 = -2; beta1 = 0.5; n = 40; p = 100; dsn = 5;
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

/* Josh-- ods output are keywords; on the LHS is another keyword that describes
          some piece of the output.  On the RHS is a name we give to the 
          data set we want to save. */ 
ods output parameterestimates = pqlestimates;
title "pql";
proc glimmix data = sim;
by dsn;
class obs;
nloptions maxiter = 100;
model y = x1 / solution dist=bin;
random int / subject = cluster;
run;


ods output parameterestimates = laplaceestimates;
title "laplace";
proc glimmix data = sim method = laplace;
by dsn;
class obs;
model y = x1 /solution dist=bin;
random int / subject = cluster;
run;


ods output parameterestimates = quadestimates;
title "quad 25";
proc glimmix data = sim method = quad(qpoints=21);
by dsn;
class obs;
model y = x1 /solution dist=bin;
random int / subject = cluster;
run;

proc print data = quadestimates (where = (effect = "x1")); run;

