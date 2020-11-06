/* Created by Language version: 7.5.0 */
/* NOT VECTORIZED */
#define NRN_VECTORIZED 0
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "scoplib_ansi.h"
#undef PI
#define nil 0
#include "md1redef.h"
#include "section.h"
#include "nrniv_mf.h"
#include "md2redef.h"
 
#if METHOD3
extern int _method3;
#endif

#if !NRNGPU
#undef exp
#define exp hoc_Exp
extern double hoc_Exp(double);
#endif
 
#define nrn_init _nrn_init__naf
#define _nrn_initial _nrn_initial__naf
#define nrn_cur _nrn_cur__naf
#define _nrn_current _nrn_current__naf
#define nrn_jacob _nrn_jacob__naf
#define nrn_state _nrn_state__naf
#define _net_receive _net_receive__naf 
#define _f_rates _f_rates__naf 
#define rates rates__naf 
#define state state__naf 
 
#define _threadargscomma_ /**/
#define _threadargsprotocomma_ /**/
#define _threadargs_ /**/
#define _threadargsproto_ /**/
 	/*SUPPRESS 761*/
	/*SUPPRESS 762*/
	/*SUPPRESS 763*/
	/*SUPPRESS 765*/
	 extern double *getarg();
 static double *_p; static Datum *_ppvar;
 
#define t nrn_threads->_t
#define dt nrn_threads->_dt
#define gnabar _p[0]
#define mshift _p[1]
#define hshift _p[2]
#define ina _p[3]
#define m _p[4]
#define h _p[5]
#define Dm _p[6]
#define Dh _p[7]
#define ena _p[8]
#define gna _p[9]
#define minf _p[10]
#define hinf _p[11]
#define _g _p[12]
#define _ion_ena	*_ppvar[0]._pval
#define _ion_ina	*_ppvar[1]._pval
#define _ion_dinadv	*_ppvar[2]._pval
#define mu	*_ppvar[3]._pval
#define _p_mu	_ppvar[3]._pval
 
#if MAC
#if !defined(v)
#define v _mlhv
#endif
#if !defined(h)
#define h _mlhh
#endif
#endif
 
#if defined(__cplusplus)
extern "C" {
#endif
 static int hoc_nrnpointerindex =  3;
 /* external NEURON variables */
 /* declaration of user functions */
 static void _hoc_rates(void);
 static void _hoc_table_tauhnaf(void);
 static void _hoc_tauhnaf(void);
 static void _hoc_table_taumnaf(void);
 static void _hoc_taumnaf(void);
 static int _mechtype;
extern void _nrn_cacheloop_reg(int, int);
extern void hoc_register_prop_size(int, int, int);
extern void hoc_register_limits(int, HocParmLimits*);
extern void hoc_register_units(int, HocParmUnits*);
extern void nrn_promote(Prop*, int, int);
extern Memb_func* memb_func;
 extern void _nrn_setdata_reg(int, void(*)(Prop*));
 static void _setdata(Prop* _prop) {
 _p = _prop->param; _ppvar = _prop->dparam;
 }
 static void _hoc_setdata() {
 Prop *_prop, *hoc_getdata_range(int);
 _prop = hoc_getdata_range(_mechtype);
   _setdata(_prop);
 hoc_retpushx(1.);
}
 /* connect user functions to hoc names */
 static VoidFunc hoc_intfunc[] = {
 "setdata_naf", _hoc_setdata,
 "rates_naf", _hoc_rates,
 "table_tauhnaf_naf", _hoc_table_tauhnaf,
 "tauhnaf_naf", _hoc_tauhnaf,
 "table_taumnaf_naf", _hoc_table_taumnaf,
 "taumnaf_naf", _hoc_taumnaf,
 0, 0
};
#define table_tauhnaf table_tauhnaf_naf
#define tauhnaf tauhnaf_naf
#define table_taumnaf table_taumnaf_naf
#define taumnaf taumnaf_naf
 extern double table_tauhnaf( );
 extern double tauhnaf( double );
 extern double table_taumnaf( );
 extern double taumnaf( double );
 /* declare global and static user variables */
#define hqfact hqfact_naf
 double hqfact = 3;
#define hslope hslope_naf
 double hslope = 10.7;
#define hvhalf hvhalf_naf
 double hvhalf = -62.9;
#define mqfact mqfact_naf
 double mqfact = 3;
#define mslope mslope_naf
 double mslope = -11.8;
#define mvhalf mvhalf_naf
 double mvhalf = -23.9;
#define usetable usetable_naf
 double usetable = 1;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 "usetable_naf", 0, 1,
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "mvhalf_naf", "mV",
 "mslope_naf", "mV",
 "hvhalf_naf", "mV",
 "hslope_naf", "mV",
 "gnabar_naf", "S/cm2",
 "mshift_naf", "mV",
 "hshift_naf", "mV",
 "ina_naf", "mA/cm2",
 "mu_naf", "1",
 0,0
};
 static double delta_t = 0.01;
 static double h0 = 0;
 static double m0 = 0;
 static double v = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 "mvhalf_naf", &mvhalf_naf,
 "mslope_naf", &mslope_naf,
 "hvhalf_naf", &hvhalf_naf,
 "hslope_naf", &hslope_naf,
 "mqfact_naf", &mqfact_naf,
 "hqfact_naf", &hqfact_naf,
 "usetable_naf", &usetable_naf,
 0,0
};
 static DoubVec hoc_vdoub[] = {
 0,0,0
};
 static double _sav_indep;
 static void nrn_alloc(Prop*);
static void  nrn_init(_NrnThread*, _Memb_list*, int);
static void nrn_state(_NrnThread*, _Memb_list*, int);
 static void nrn_cur(_NrnThread*, _Memb_list*, int);
static void  nrn_jacob(_NrnThread*, _Memb_list*, int);
 
static int _ode_count(int);
static void _ode_map(int, double**, double**, double*, Datum*, double*, int);
static void _ode_spec(_NrnThread*, _Memb_list*, int);
static void _ode_matsol(_NrnThread*, _Memb_list*, int);
 
#define _cvode_ieq _ppvar[4]._i
 static void _ode_matsol_instance1(_threadargsproto_);
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "7.5.0",
"naf",
 "gnabar_naf",
 "mshift_naf",
 "hshift_naf",
 0,
 "ina_naf",
 0,
 "m_naf",
 "h_naf",
 0,
 "mu_naf",
 0};
 static Symbol* _na_sym;
 
extern Prop* need_memb(Symbol*);

static void nrn_alloc(Prop* _prop) {
	Prop *prop_ion;
	double *_p; Datum *_ppvar;
 	_p = nrn_prop_data_alloc(_mechtype, 13, _prop);
 	/*initialize range parameters*/
 	gnabar = 1.5;
 	mshift = 0;
 	hshift = 0;
 	_prop->param = _p;
 	_prop->param_size = 13;
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 5, _prop);
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 prop_ion = need_memb(_na_sym);
 nrn_promote(prop_ion, 0, 1);
 	_ppvar[0]._pval = &prop_ion->param[0]; /* ena */
 	_ppvar[1]._pval = &prop_ion->param[3]; /* ina */
 	_ppvar[2]._pval = &prop_ion->param[4]; /* _ion_dinadv */
 
}
 static void _initlists();
  /* some states have an absolute tolerance */
 static Symbol** _atollist;
 static HocStateTolerance _hoc_state_tol[] = {
 0,0
};
 static void _update_ion_pointer(Datum*);
 extern Symbol* hoc_lookup(const char*);
extern void _nrn_thread_reg(int, int, void(*)(Datum*));
extern void _nrn_thread_table_reg(int, void(*)(double*, Datum*, Datum*, _NrnThread*, int));
extern void hoc_register_tolerance(int, HocStateTolerance*, Symbol***);
extern void _cvode_abstol( Symbol**, double*, int);

 void _naf_reg() {
	int _vectorized = 0;
  _initlists();
 	ion_reg("na", -10000.);
 	_na_sym = hoc_lookup("na_ion");
 	register_mech(_mechanism, nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init, hoc_nrnpointerindex, 0);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_setdata_reg(_mechtype, _setdata);
     _nrn_thread_reg(_mechtype, 2, _update_ion_pointer);
  hoc_register_prop_size(_mechtype, 13, 5);
  hoc_register_dparam_semantics(_mechtype, 0, "na_ion");
  hoc_register_dparam_semantics(_mechtype, 1, "na_ion");
  hoc_register_dparam_semantics(_mechtype, 2, "na_ion");
  hoc_register_dparam_semantics(_mechtype, 3, "pointer");
  hoc_register_dparam_semantics(_mechtype, 4, "cvodeieq");
 	hoc_register_cvode(_mechtype, _ode_count, _ode_map, _ode_spec, _ode_matsol);
 	hoc_register_tolerance(_mechtype, _hoc_state_tol, &_atollist);
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 naf /export/home1/users/bssn/serna/HBP/SpineModel/mod/x86_64/naf.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
 static double *_t_minf;
 static double *_t_hinf;
static int _reset;
static char *modelname = "Nahc - Fast sodium current for nucleus accumbens (from hippocampal pyramidal cell)";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
static int _f_rates(double);
static int rates(double);
 
static int _ode_spec1(_threadargsproto_);
/*static int _ode_matsol1(_threadargsproto_);*/
 static void _n_rates(double);
 static int _slist1[2], _dlist1[2];
 static int state(_threadargsproto_);
 
static void* _ptable_tauhnaf = (void*)0;
 
static void* _ptable_taumnaf = (void*)0;
 
double taumnaf (  double _lv ) {
 double _arg[1];
 _arg[0] = _lv;
 return hoc_func_table(_ptable_taumnaf, 1, _arg);
 }
/*  }
  */
 
static void _hoc_taumnaf(void) {
  double _r;
   _r =  taumnaf (  *getarg(1) );
 hoc_retpushx(_r);
}
 double table_taumnaf ( ) {
	hoc_spec_table(&_ptable_taumnaf, 1);
	return 0.;
}
 
static void _hoc_table_taumnaf(void) {
  double _r;
   _r =  table_taumnaf (  );
 hoc_retpushx(_r);
}
 
double tauhnaf (  double _lv ) {
 double _arg[1];
 _arg[0] = _lv;
 return hoc_func_table(_ptable_tauhnaf, 1, _arg);
 }
/*  }
  */
 
static void _hoc_tauhnaf(void) {
  double _r;
   _r =  tauhnaf (  *getarg(1) );
 hoc_retpushx(_r);
}
 double table_tauhnaf ( ) {
	hoc_spec_table(&_ptable_tauhnaf, 1);
	return 0.;
}
 
static void _hoc_table_tauhnaf(void) {
  double _r;
   _r =  table_tauhnaf (  );
 hoc_retpushx(_r);
}
 
/*CVODE*/
 static int _ode_spec1 () {_reset=0;
 {
   rates ( _threadargscomma_ v ) ;
   Dm = ( minf - m ) / ( taumnaf ( _threadargscomma_ v ) / mqfact ) ;
   Dh = ( hinf - h ) / ( tauhnaf ( _threadargscomma_ v ) / hqfact ) ;
   }
 return _reset;
}
 static int _ode_matsol1 () {
 rates ( _threadargscomma_ v ) ;
 Dm = Dm  / (1. - dt*( ( ( ( - 1.0 ) ) ) / ( taumnaf ( _threadargscomma_ v ) / mqfact ) )) ;
 Dh = Dh  / (1. - dt*( ( ( ( - 1.0 ) ) ) / ( tauhnaf ( _threadargscomma_ v ) / hqfact ) )) ;
  return 0;
}
 /*END CVODE*/
 static int state () {_reset=0;
 {
   rates ( _threadargscomma_ v ) ;
    m = m + (1. - exp(dt*(( ( ( - 1.0 ) ) ) / ( taumnaf ( _threadargscomma_ v ) / mqfact ))))*(- ( ( ( minf ) ) / ( taumnaf ( _threadargscomma_ v ) / mqfact ) ) / ( ( ( ( - 1.0 ) ) ) / ( taumnaf ( _threadargscomma_ v ) / mqfact ) ) - m) ;
    h = h + (1. - exp(dt*(( ( ( - 1.0 ) ) ) / ( tauhnaf ( _threadargscomma_ v ) / hqfact ))))*(- ( ( ( hinf ) ) / ( tauhnaf ( _threadargscomma_ v ) / hqfact ) ) / ( ( ( ( - 1.0 ) ) ) / ( tauhnaf ( _threadargscomma_ v ) / hqfact ) ) - h) ;
   }
  return 0;
}
 static double _mfac_rates, _tmin_rates;
 static void _check_rates();
 static void _check_rates() {
  static int _maktable=1; int _i, _j, _ix = 0;
  double _xi, _tmax;
  static double _sav_mshift;
  static double _sav_hshift;
  static double _sav_mslope;
  static double _sav_hslope;
  if (!usetable) {return;}
  if (_sav_mshift != mshift) { _maktable = 1;}
  if (_sav_hshift != hshift) { _maktable = 1;}
  if (_sav_mslope != mslope) { _maktable = 1;}
  if (_sav_hslope != hslope) { _maktable = 1;}
  if (_maktable) { double _x, _dx; _maktable=0;
   _tmin_rates =  - 200.0 ;
   _tmax =  200.0 ;
   _dx = (_tmax - _tmin_rates)/201.; _mfac_rates = 1./_dx;
   for (_i=0, _x=_tmin_rates; _i < 202; _x += _dx, _i++) {
    _f_rates(_x);
    _t_minf[_i] = minf;
    _t_hinf[_i] = hinf;
   }
   _sav_mshift = mshift;
   _sav_hshift = hshift;
   _sav_mslope = mslope;
   _sav_hslope = hslope;
  }
 }

 static int rates(double _lv){ _check_rates();
 _n_rates(_lv);
 return 0;
 }

 static void _n_rates(double _lv){ int _i, _j;
 double _xi, _theta;
 if (!usetable) {
 _f_rates(_lv); return; 
}
 _xi = _mfac_rates * (_lv - _tmin_rates);
 if (isnan(_xi)) {
  minf = _xi;
  hinf = _xi;
  return;
 }
 if (_xi <= 0.) {
 minf = _t_minf[0];
 hinf = _t_hinf[0];
 return; }
 if (_xi >= 201.) {
 minf = _t_minf[201];
 hinf = _t_hinf[201];
 return; }
 _i = (int) _xi;
 _theta = _xi - (double)_i;
 minf = _t_minf[_i] + _theta*(_t_minf[_i+1] - _t_minf[_i]);
 hinf = _t_hinf[_i] + _theta*(_t_hinf[_i+1] - _t_hinf[_i]);
 }

 
static int  _f_rates (  double _lv ) {
   minf = 1.0 / ( 1.0 + exp ( ( _lv - mvhalf - mshift ) / mslope ) ) ;
   hinf = 1.0 / ( 1.0 + exp ( ( _lv - hvhalf - hshift ) / hslope ) ) ;
    return 0; }
 
static void _hoc_rates(void) {
  double _r;
    _r = 1.;
 rates (  *getarg(1) );
 hoc_retpushx(_r);
}
 
static int _ode_count(int _type){ return 2;}
 
static void _ode_spec(_NrnThread* _nt, _Memb_list* _ml, int _type) {
   Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
  ena = _ion_ena;
     _ode_spec1 ();
  }}
 
static void _ode_map(int _ieq, double** _pv, double** _pvdot, double* _pp, Datum* _ppd, double* _atol, int _type) { 
 	int _i; _p = _pp; _ppvar = _ppd;
	_cvode_ieq = _ieq;
	for (_i=0; _i < 2; ++_i) {
		_pv[_i] = _pp + _slist1[_i];  _pvdot[_i] = _pp + _dlist1[_i];
		_cvode_abstol(_atollist, _atol, _i);
	}
 }
 
static void _ode_matsol_instance1(_threadargsproto_) {
 _ode_matsol1 ();
 }
 
static void _ode_matsol(_NrnThread* _nt, _Memb_list* _ml, int _type) {
   Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
  ena = _ion_ena;
 _ode_matsol_instance1(_threadargs_);
 }}
 extern void nrn_update_ion_pointer(Symbol*, Datum*, int, int);
 static void _update_ion_pointer(Datum* _ppvar) {
   nrn_update_ion_pointer(_na_sym, _ppvar, 0, 0);
   nrn_update_ion_pointer(_na_sym, _ppvar, 1, 3);
   nrn_update_ion_pointer(_na_sym, _ppvar, 2, 4);
 }

static void initmodel() {
  int _i; double _save;_ninits++;
 _save = t;
 t = 0.0;
{
  h = h0;
  m = m0;
 {
   rates ( _threadargscomma_ v ) ;
   m = minf ;
   h = hinf ;
   }
  _sav_indep = t; t = _save;

}
}

static void nrn_init(_NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; double _v; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 v = _v;
  ena = _ion_ena;
 initmodel();
 }}

static double _nrn_current(double _v){double _current=0.;v=_v;{ {
   gna = gnabar * m * m * m * h ;
   ina = gna * ( v - ena ) * ( 1.0 - ( mu - 1.0 ) * 0.05 ) ;
   }
 _current += ina;

} return _current;
}

static void nrn_cur(_NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; int* _ni; double _rhs, _v; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
  ena = _ion_ena;
 _g = _nrn_current(_v + .001);
 	{ double _dina;
  _dina = ina;
 _rhs = _nrn_current(_v);
  _ion_dinadv += (_dina - ina)/.001 ;
 	}
 _g = (_g - _rhs)/.001;
  _ion_ina += ina ;
#if CACHEVEC
  if (use_cachevec) {
	VEC_RHS(_ni[_iml]) -= _rhs;
  }else
#endif
  {
	NODERHS(_nd) -= _rhs;
  }
 
}}

static void nrn_jacob(_NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml];
#if CACHEVEC
  if (use_cachevec) {
	VEC_D(_ni[_iml]) += _g;
  }else
#endif
  {
     _nd = _ml->_nodelist[_iml];
	NODED(_nd) += _g;
  }
 
}}

static void nrn_state(_NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; double _v = 0.0; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
 _nd = _ml->_nodelist[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 v=_v;
{
  ena = _ion_ena;
 { error =  state();
 if(error){fprintf(stderr,"at line 54 in file naf.mod:\n        SOLVE state METHOD cnexp\n"); nrn_complain(_p); abort_run(error);}
 } }}

}

static void terminal(){}

static void _initlists() {
 int _i; static int _first = 1;
  if (!_first) return;
 _slist1[0] = &(m) - _p;  _dlist1[0] = &(Dm) - _p;
 _slist1[1] = &(h) - _p;  _dlist1[1] = &(Dh) - _p;
   _t_minf = makevector(202*sizeof(double));
   _t_hinf = makevector(202*sizeof(double));
_first = 0;
}
