%include "C:\Users\User\Dropbox\4° Série\Modelos Lineares Generalizados\trabalho\scripts\tablatex.sas";
proc delete data = _all_; run;
proc import out = elisa replace
			datafile = 'C:\Users\User\Dropbox\4° Série\Modelos Não Lineares\trabalhos\trabalho1\elisa.txt'
			dbms = tab;
			getnames = yes;
run;

proc sgplot data = elisa;
	scatter x = logd y = od / group = month;
run;

/*Estimativas dos parâmetros */
proc nlmixed data = elisa df = 999999 cov tech = quanew;
	parms theta1_may  = 0.05, theta2_may  = 1.9, theta3_may  = 2.5, theta4_may  = 3.2,
		  theta1_june = 0.05, theta2_june = 1.9, theta3_june = 2.5, theta4_june = 3.2, sigma = 1;
	bounds sigma > 0;
	if(month = 'may') then do;
		mean = theta1_may + (theta2_may - theta1_may) / (1 + exp(theta3_may * (logd - theta4_may)));
	end; else do;
		mean = theta1_june + (theta2_june - theta1_june) / (1 + exp(theta3_june * (logd - theta4_june)));
	end;
	ll 		=	- 0.5 * log(2*constant("pi")) - 0.5 * log(sigma**2) - 0.5 / sigma**2 * (od - mean)**2;
	model od ~ general(ll);
	ods output   ParameterEstimates=estimativas(drop = df tValue Probt Alpha Gradient);
	estimate 'beta' theta4_may - theta4_june;
	estimate 'rho'  10**(-(theta4_may - theta4_june));
	predict mean out = media;
run;
%uplatex(x = 'est.tex');
proc print data = estimativas noobs;
format _numeric_ comma10.4 ; run;
%downlatex;

proc nlin data = elisa;
	parms theta1_may  = 0.05, theta2_may  = 1.9, theta3_may  = 2.5, theta4_may  = 3.2,
		  theta1_june = 0.05, theta2_june = 1.9, theta3_june = 2.5, theta4_june = 3.2;
	if(month = 'may') then do;
		mean = theta1_may + (theta2_may - theta1_may) / (1 + exp(theta3_may * (logd - theta4_may)));
	end; else do;
		mean = theta1_june + (theta2_june - theta1_june) / (1 + exp(theta3_june * (logd - theta4_june)));
	end;
	model od = mean;
run;


*Testando a hipótese ->	H_0: theta_may = theta_june vs.	H_1: theta_may =! theta_june;

* Modelo completo;
proc nlmixed data = elisa df=99999;
	parms theta1_may  = 0.05, theta2_may  = 1.9, theta3_may  = 2.5, theta4_may  = 3.2,
		  theta1_june = 0.05, theta2_june = 1.9, theta3_june = 2.5, theta4_june = 3.2, sigma = 1;
	bounds sigma > 0;
	if(month = 'may') then do;
		mean = theta1_may + (theta2_may - theta1_may) / (1 + exp(theta3_may * (logd - theta4_may)));
	end; else do;
		mean = theta1_june + (theta2_june - theta1_june) / (1 + exp(theta3_june * (logd - theta4_june)));
	end;
	ll 		 =	- 0.5 * log(2*constant("pi")) - 0.5 * log(sigma**2) - 0.5 / sigma**2 * (od - mean)**2;
	model od ~ general(ll);
	contrast 'theta_may = theta_june' theta1_may - theta1_june, theta2_may - theta2_june, theta3_may - theta3_june;
	ods output Contrasts 	= Wald(drop = label DenDF rename = (NumDF = df FValue = est ProbF = pvalue));
	ods output FitStatistics = H1(rename = (value = llh1 ) where = (descr = '-2 Log Likelihood'));
quit;

* Modelo restrito;
proc nlmixed data = elisa df=99999;
	parms theta1  = 0.05, theta2  = 1.9, theta3  = 2.5, theta4_may  = 3.2, theta4_june = 3.2, sigma = 1;
	bounds sigma >0;
	if(month = 'may') then do;
		mean = theta1 + (theta2 - theta1) / (1 + exp(theta3 * (logd - theta4_may)));
	end; else do;
		mean = theta1 + (theta2 - theta1) / (1 + exp(theta3 * (logd - theta4_june)));
	end;
	ll 		 =	- 0.5 * log(2*constant("pi")) - 0.5 * log(sigma**2) - 0.5 / sigma**2 * (od - mean)**2;
	model od ~ general(ll);
	ods output FitStatistics = H0(rename = (value = llh0 ) where = (descr = '-2 Log Likelihood'));
run;

* Ajuste dos resultados;
data LR(drop = Descr llH0 llH1);
	merge H0 H1;
	df		= 3;
	Est 	= llH0 - llH1;
	pvalue	= 1 - cdf('CHISQUARE', Est, df);
run;	

data wald;
	set wald;
	est = df * est;
run;

data testes;
	set Wald LR;
run;

* Exibindo resultados;
proc print data = testes noobs;
	format _numeric_ comma10.4;
run;

*Testando a hipótese
	H_0: beta =  0
	H_1: beta =! 0;
proc nlmixed data = elisa df=99999;
	parms theta1_may  = 0.05, theta2_may  = 1.9, theta3_may  = 2.5, theta4_may  = 3.2,
		  theta1_june = 0.05, theta2_june = 1.9, theta3_june = 2.5, theta4_june = 3.2, sigma = 1;
	bounds sigma > 0;
	if(month = 'may') then do;
		mean = theta1_may + (theta2_may - theta1_may) / (1 + exp(theta3_may * (logd - theta4_may)));
	end; else do;
		mean = theta1_june + (theta2_june - theta1_june) / (1 + exp(theta3_june * (logd - theta4_june)));
	end;
	ll 		 =	- 0.5 * log(2*constant("pi")) - 0.5 * log(sigma**2) - 0.5 / sigma**2 * (od - mean)**2;
	model od ~ general(ll);
	estimate 'beta' theta4_may - theta4_june;
	contrast 'beta = 0' theta4_may - theta4_june;
	ods output FitStatistics = H1(rename = (value = llh1 ) where = (descr = '-2 Log Likelihood'));
quit;

proc nlmixed data = elisa df=99999;
	parms theta1_may  = 0.05, theta2_may  = 1.9, theta3_may  = 2.5, theta4 = 3.2,
		  theta1_june = 0.05, theta2_june = 1.9, theta3_june = 2.5, sigma = 1;
	bounds sigma > 0;
	if(month = 'may') then do;
		mean = theta1_may + (theta2_may - theta1_may) / (1 + exp(theta3_may * (logd - theta4)));
	end; else do;
		mean = theta1_june + (theta2_june - theta1_june) / (1 + exp(theta3_june * (logd - theta4)));
	end;
	ll 		 =	- 0.5 * log(2*constant("pi")) - 0.5 * log(sigma**2) - 0.5 / sigma**2 * (od - mean)**2;
	model od ~ general(ll);
	ods output FitStatistics = H0(rename = (value = llh0) where = (descr = '-2 Log Likelihood'));
quit;
data LR(drop = Descr llH0 llH1);
	merge H0 H1;
	df		= 1;
	Est 	= llH0 - llH1;
	pvalue	= 1 - cdf('CHISQUARE', Est, df);
run;
proc print data = lr noobs;
	format _numeric_ comma10.4;
run;
