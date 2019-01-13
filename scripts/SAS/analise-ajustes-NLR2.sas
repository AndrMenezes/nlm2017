proc delete data = _all_; run;
%include "C:\Users\User\Dropbox\4° Série\Modelos Lineares Generalizados\trabalho\scripts\tablatex.sas";

proc import out = dados replace dbms = tab
			datafile = 'C:\Users\User\Dropbox\4° Série\Modelos Não Lineares\trabalhos\trabalho2\goiaba.txt';
			getnames = yes;
run;

***************************** Modelos Homocedásticos *****************************;

/*Modelo 4-Gompertz - ei ~ Normal(0, sigma) */
proc nlmixed data = dados tech = tr  df = 99999;
	parms  theta1 = 190, theta2 = 18, theta3 = 0.09, theta4 = 107, sigma = 1;
	bounds sigma > 0;
	mu	 	= 	theta1 - (theta1 - theta2) * exp(-exp(theta3 * (daa - theta4) ));
	sigma2i =	sigma**2;
	ll		=	- 0.5 * log(2*constant("pi")) - 0.5 * log(sigma2i) - 0.5 / (sigma2i) * (peso - mu)**2;
*	ll		=	logpdf('normal', peso, mu, sigma * daa ** lambda);
	model peso ~ general(ll);
	ods output ParameterEstimates = mhm1(drop = df alpha gradient tValue Probt lower upper
			   rename = (estimate = estmhm1 StandardError = stdmhm1));
	predict mu out = pred1(drop = DF tValue alpha StdErrPred Probt rename = (Pred = predhm1 lower = lowerhm1 upper = upperhm1));
run;

/*Modelo 4-Gompertz com ei ~ Gumbel(0, sigma) */
proc nlmixed data = dados tech = tr  df = 99999;
	parms  theta1 = 190, theta2 = 18, theta3 = 0.09, theta4 = 107, sigma = 20;
	bounds sigma > 0;
	mu	 	= 	theta1 - (theta1 - theta2) * exp(-exp(theta3 * (daa - theta4) ));
	Ex	 	= 	mu + sigma * constant('EULER');
	zedd	=	(peso - mu) / sigma;
	ll		=	-zedd - exp(-zedd) - log(sigma);
	model peso ~ general(ll);
	ods output ParameterEstimates = mhm2(drop = df alpha gradient tValue Probt lower upper
			   rename = (estimate = estmhm2 StandardError = stdmhm2));
	predict Ex out = pred2(keep = pred lower upper rename = (Pred = predhm2 lower = lowerhm2 upper = upperhm2));
run;

***************************** Modelos Heterocedásticos *****************************;

/*Modelo 4-Gompertz com ei ~ Normal(0, sigma) e função variância potência*/
proc nlmixed data = dados tech = quanew update=bfgs df = 99999;
	parms  theta1 = 190, theta2 = 18, theta3 = 0.09, theta4 = 107, sigma = 1, lambda = 1;
	bounds sigma > 0, lambda > 0;
	mu	 	= 	theta1 - (theta1 - theta2) * exp(-exp(theta3 * (daa - theta4) ));
	sigmai  =	sigma * (daa ** lambda);
	ll		=	- 0.5 * log(2*constant("pi")) - 0.5 * log(sigmai**2) - 0.5 * ( (peso - mu) / sigmai )**2;
	model peso ~ general(ll);
	ods output ParameterEstimates = mht1(drop = df alpha gradient tValue Probt lower upper
			   rename = (estimate = estmht1 StandardError = stdmht1));
	predict mu out = pred3(keep = pred lower upper rename = (Pred = predht1 lower = lowerht1 upper = upperht1));
run;

/*Modelo 4-Gompertz com ei ~ Gumbel(0, sigma) e função variância potência*/
proc nlmixed data = dados tech = tr  df = 99999;
	parms  theta1 = 190, theta2 = 18, theta3 = 0.09, theta4 = 107, sigma = 20, lambda = 1;
	bounds sigma > 0, lambda >0;
	sigmai  =	sigma * (daa ** lambda);
	mu	 	= 	theta1 - (theta1 - theta2) * exp(-exp(theta3 * (daa - theta4) ));
	Ex	 	= 	mu + sigmai * constant('EULER');
	zedd	=	(peso - mu) / sigmai;
	ll		=	-zedd - exp(-zedd) - log(sigmai);
	model peso ~ general(ll);
	ods output ParameterEstimates = mht2(drop = df alpha gradient tValue Probt lower upper
			   rename = (estimate = estmht2 StandardError = stdmht2));
	predict Ex out = pred4(keep = pred lower upper rename = (Pred = predht2 lower = lowerht2 upper = upperht2));
run;

/*Modelo 4-Gompertz com ei ~ Normal(0, sigma) e função variância exponencial*/
proc nlmixed data = dados tech = quanew update=bfgs df = 99999;
	parms  theta1 = 190, theta2 = 18, theta3 = 0.09, theta4 = 107, sigma = 1, lambda = 1;
	bounds sigma > 0, lambda > 0;
	mu	 	= 	theta1 - (theta1 - theta2) * exp(-exp(theta3 * (daa - theta4) ));
	sigmai  =	sigma * exp(daa * lambda);
	ll		=	- 0.5 * log(2*constant("pi")) - log(sigmai) - 0.5 * ( (peso - mu) / sigmai )**2;
	model peso ~ general(ll);
	ods output ParameterEstimates = mht3(drop = df alpha gradient tValue Probt lower upper
			   rename = (estimate = estmht3 StandardError = stdmht3));
	predict mu out = pred5(keep = pred lower upper rename = (Pred = predht3 lower = lowerht3 upper = upperht3));
run;

/*Modelo 4-Gompertz com ei ~ Gumbel(0, sigma) e função variância exponencial*/
proc nlmixed data = dados tech = tr  df = 99999;
	parms  theta1 = 190, theta2 = 18, theta3 = 0.09, theta4 = 107, sigma = 20, lambda = 1;
	bounds sigma > 0, lambda >0;
	sigmai  =	sigma * exp(daa * lambda);
	mu	 	= 	theta1 - (theta1 - theta2) * exp(-exp(theta3 * (daa - theta4) ));
	Ex	 	= 	mu + sigmai * constant('EULER');
	zedd	=	(peso - mu) / sigmai;
	ll		=	-zedd - exp(-zedd) - log(sigmai);
	model peso ~ general(ll);
	ods output ParameterEstimates = mht4(drop = df alpha gradient tValue Probt lower upper
			   rename = (estimate = estmht4 StandardError = stdmht4));
	predict Ex out = pred6(keep = pred lower upper rename = (Pred = predht4 lower = lowerht4 upper = upperht4));
run;


data est1;
	merge mhm1 mhm2;
run;
data est2;
	merge mht1 mht2 mht3 mht4;
run;

%uplatex(x = 'est1.tex');
proc print data = est1 noobs;
	format _numeric_ comma10.4;
run;
%downlatex;
%uplatex(x = 'est2.tex');
proc print data = est2 noobs;
	format _numeric_ comma10.4;
run;
%downlatex;

data pred;
	merge pred1 pred2 pred3 pred4 pred5 pred6;
	drop coleta rep long trans volume;
run;
proc export outfile = 'pred.txt' replace dbms = tab; run;

data preditos;
	set pred;
	keep peso daa predhm1 predhm2 predht1 predht2 predht3 predht4;
run;

proc export data = preditos outfile = 'ypred.txt' replace dbms = tab; run;















