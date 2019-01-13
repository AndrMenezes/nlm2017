proc delete data = _all_; run;
%include "C:\Users\User\Dropbox\4° Série\Modelos Lineares Generalizados\trabalho\scripts\tablatex.sas";

proc import out = dados replace dbms = tab
			datafile = 'C:\Users\User\Dropbox\4° Série\Modelos Não Lineares\trabalhos\trabalho2\ndados.txt';
			getnames = yes;
run;

***************************** Modelos Homocedásticos *****************************;

/*Modelo 4-Gompertz - ei ~ Normal(0, sigma) */
proc nlmixed data = dados tech = tr  df = 99999;
	parms  theta1 = 190, theta2 = 18, theta3 = 0.09, theta4 = 107, sigma = 1;
	bounds sigma > 0;
	mu	 	= 	theta1 - (theta1 - theta2) * exp(-exp(theta3 * (daa - theta4) ));
	sigma2i =	sigma**2;
	if dev = 1 then do;
		ll		=	- 0.5 * log(2*constant("pi")) - 0.5 * log(sigma2i) - 0.5 / (sigma2i) * (peso - mu)**2;
	end; else ll = 0;
	model peso ~ general(ll);
	predict mu out = pred1(drop = DF tValue alpha StdErrPred Probt rename = (Pred = predhm1 lower = lowerhm1 upper = upperhm1));
run;

/*Modelo 4-Gompertz com ei ~ Gumbel(0, sigma) */
proc nlmixed data = dados tech = tr  df = 99999;
	parms  theta1 = 190, theta2 = 18, theta3 = 0.09, theta4 = 107, sigma = 20;
	bounds sigma > 0;
	mu	 	= 	theta1 - (theta1 - theta2) * exp(-exp(theta3 * (daa - theta4) ));
	Ex	 	= 	mu + sigma * constant('EULER');
	zedd	=	(peso - mu) / sigma;
	if dev = 1 then do;
		ll		=	-zedd - exp(-zedd) - log(sigma);
	end; else ll = 0;
	model peso ~ general(ll);
	predict Ex out = pred2(keep = pred lower upper rename = (Pred = predhm2 lower = lowerhm2 upper = upperhm2));
run;

***************************** Modelos Heterocedásticos *****************************;

/*Modelo 4-Gompertz com ei ~ Normal(0, sigma) e função variância potência*/
proc nlmixed data = dados tech = quanew update=bfgs df = 99999;
	parms  theta1 = 190, theta2 = 18, theta3 = 0.09, theta4 = 107, sigma = 1, lambda = 1;
	bounds sigma > 0, lambda > 0;
	mu	 	= 	theta1 - (theta1 - theta2) * exp(-exp(theta3 * (daa - theta4) ));
	sigmai  =	sigma * (daa ** lambda);
	if dev = 1 then do;
		ll		=	- 0.5 * log(2*constant("pi")) - 0.5 * log(sigmai**2) - 0.5 * ( (peso - mu) / sigmai )**2;
	end; else ll = 0;
	model peso ~ general(ll);
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
	if dev = 1 then do;
		ll		=	-zedd - exp(-zedd) - log(sigmai);
	end; else ll = 0;
	model peso ~ general(ll);
	predict Ex out = pred4(keep = pred lower upper rename = (Pred = predht2 lower = lowerht2 upper = upperht2));
run;

/*Modelo 4-Gompertz com ei ~ Normal(0, sigma) e função variância exponencial*/
proc nlmixed data = dados tech = quanew update=bfgs df = 99999;
	parms  theta1 = 190, theta2 = 18, theta3 = 0.09, theta4 = 107, sigma = 1, lambda = 1;
	bounds sigma > 0, lambda > 0;
	mu	 	= 	theta1 - (theta1 - theta2) * exp(-exp(theta3 * (daa - theta4) ));
	sigmai  =	sigma * exp(daa * lambda);
	if dev = 1 then do;
		ll		=	- 0.5 * log(2*constant("pi")) - log(sigmai) - 0.5 * ( (peso - mu) / sigmai )**2;
	end; else ll = 0;
	model peso ~ general(ll);
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
	if dev = 1 then do;
		ll		=	-zedd - exp(-zedd) - log(sigmai);
	end; else ll = 0;
	model peso ~ general(ll);
	predict Ex out = pred6(keep = pred lower upper rename = (Pred = predht4 lower = lowerht4 upper = upperht4));
run;

data pred;
	merge pred1 pred2 pred3 pred4 pred5 pred6;
run;
proc export outfile = 'predicao.txt' replace dbms = tab; run;















