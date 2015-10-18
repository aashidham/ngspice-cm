#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define CURRENT 0
#define PREV 1	

//double max ( double a, double b ) { return a > b ? a : b; }
//double min ( double a, double b ) { return a < b ? a : b; }

#define min(a,b) fmin(a,b)
#define max(a,b) fmax(a,b)


void cm_mig(ARGS)
{
	double vk = -90.0e-003;
	double vna = 55e-003;
	double vleak = -65e-3;
	double area = PARAM(area);
	double v1 = INPUT(cap);
	double cm = 1e-6*area;
	double rleak=3333.33333333/area;
	double gna = 32e-3*area;
	double gk_dr = 10e-3*area;
	double gk_a = 48e-3*area;
	
	double prev_volt;
	double* m;
	double* h;
	double* i;
	double* n;
	double* n2;
	double* l;
	
	
	if(INIT)
	{
		cm_analog_alloc(0,sizeof(double));
		cm_analog_alloc(1,sizeof(double));
		cm_analog_alloc(2,sizeof(double));
		cm_analog_alloc(3,sizeof(double));
		cm_analog_alloc(4,sizeof(double));
		cm_analog_alloc(5,sizeof(double));
		cm_analog_alloc(6,sizeof(double));
		v1 = -65e-003;
		prev_volt = v1;
	}
	else
	{
		prev_volt = *(double *) cm_analog_get_ptr(0,PREV);
	}

	*(double *) cm_analog_get_ptr(0,CURRENT) = v1;
	m = (double *) cm_analog_get_ptr(1,CURRENT);
	h = (double *) cm_analog_get_ptr(2,CURRENT);
	i = (double *) cm_analog_get_ptr(3,CURRENT);
	n = (double *) cm_analog_get_ptr(4,CURRENT);
	n2 = (double *) cm_analog_get_ptr(5,CURRENT);
	l = (double *) cm_analog_get_ptr(6,CURRENT);
	
	if(ANALYSIS==TRANSIENT)
	{
		double timestep = T(0)-T(1);
		double dv_dt = (v1-prev_volt)/timestep;
		printf("%17.15f %17.15f %17.15f\n",prev_volt,v1,timestep);
		double m_old = *(double *) cm_analog_get_ptr(1,PREV);
		double h_old = *(double *) cm_analog_get_ptr(2,PREV);
		double i_old = *(double *) cm_analog_get_ptr(3,PREV);
		double n_old = *(double *) cm_analog_get_ptr(4,PREV);
		double n2_old = *(double *) cm_analog_get_ptr(5,PREV);
		double l_old = *(double *) cm_analog_get_ptr(6,PREV);
		
		//Na channel
		double alpha_m = (0.4*(v1*1e003 + 30))/(1-exp(-1*(v1*1e003+30)/7.2));
		double beta_m = (0.124*(v1*1e003+30))/(exp((v1*1e003+30)/7.2)-1);
		double m_tau = 0.5/(alpha_m+beta_m);
		double m_inf = alpha_m/(alpha_m+beta_m);
		
		double alpha_h = (0.03*(v1*1e003 + 45))/(1-exp(-1*(v1*1e003+45)/1.5));
		double beta_h = (0.01*(v1*1e003+45))/(exp((v1*1e003+45)/1.5)-1);
		double h_tau = 0.5/(alpha_h+beta_h);
		double h_inf = 1/(1+exp((v1*1e003+50)/4));
		
		double alpha_i = exp(0.45*(v1*1e003+60));
		double beta_i = exp(0.09*(v1*1e003+60));
		double i_tau = max((3e4*beta_i/(1+alpha_i)),10);
		double i_inf = (1+0.8*exp((v1*1e003+58)/2))/(1+exp((v1*1e003+58)/2));
		
		double dm_dt = 1e3*(m_inf - m_old) / m_tau;
		*m = dm_dt*timestep + m_old;
		
		double dh_dt = 1e3*(h_inf - h_old) / h_tau;
		*h = dh_dt*timestep + h_old;
		
		double di_dt = 1e3*(i_inf - i_old) / i_tau;
		*i = di_dt*timestep + i_old;
		
		//K(Dr) channel
		double alpha_n = exp(-0.11*(v1*1e003-13));
		double beta_n = exp(-0.08*(v1*1e003-13));
		double n_tau = max(2,(50*beta_n/(1+alpha_n)));
		double n_inf = 1/(1+alpha_n);
		
		double dn_dt = 1e3*(n_inf - n_old) / n_tau;
		*n = dn_dt*timestep + n_old;
		
		//K(A) channel
		double alpha_n2 = exp(-0.038*(1.5+(1/(1+exp((v1*1e003+40)/5))))*(v1*1e003-11));
		double beta_n2 = exp(-0.038*(0.825+(1/(1+exp((v1*1e003+40)/5))))*(v1*1e003-11));
		double n2_tau = max(0.1,(4*beta_n2/(1+alpha_n2)));
		double n2_inf = 1/(1+alpha_n2);
		
		double alpha_l = exp(0.11*(v1*1e003+56));
		double l_tau = max(2,(0.26*(v1*1e003+50)));
		double l_inf = 1/(1+alpha_l);
		
		double dn2_dt = 1e3*(n2_inf - n2_old) / n2_tau;
		*n2 = dn2_dt*timestep + n2_old;
		
		double dl_dt = 1e3*(l_inf - l_old) / l_tau;
		*l = dl_dt*timestep + l_old;
		
		double m_val = *m;
		double h_val = *h;
		double i_val = *i;
		double n_val = *n;
		double n2_val = *n2;
		double l_val = *l;
		
		OUTPUT(cap) = gna*h_val*m_val*m_val*m_val*i_val*(v1-vna) + gk_dr*n_val*(v1-vk) + gk_a*n2_val*l_val*(v1-vk) + (v1-vleak)/rleak + cm*dv_dt;
		
	}
}