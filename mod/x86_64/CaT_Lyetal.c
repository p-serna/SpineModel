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
 
#define nrn_init _nrn_init__itL
#define _nrn_initial _nrn_initial__itL
#define nrn_cur _nrn_cur__itL
#define _nrn_current _nrn_current__itL
#define nrn_jacob _nrn_jacob__itL
#define nrn_state _nrn_state__itL
#define _net_receive _net_receive__itL 
#define _f_trates _f_trates__itL 
#define rates rates__itL 
#define states states__itL 
#define trates trates__itL 
 
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
#define gbar _p[0]
#define gca _p[1]
#define minf _p[2]
#define hinf _p[3]
#define mtau _p[4]
#define htau _p[5]
#define m _p[6]
#define h _p[7]
#define ica _p[8]
#define eca _p[9]
#define tadj _p[10]
#define Dm _p[11]
#define Dh _p[12]
#define _g _p[13]
#define _ion_eca	*_ppvar[0]._pval
#define _ion_ica	*_ppvar[1]._pval
#define _ion_dicadv	*_ppvar[2]._pval
 
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
 static int hoc_nrnpointerindex =  -1;
 /* external NEURON variables */
 extern double celsius;
 /* declaration of user functions */
 static void _hoc_rates(void);
 static void _hoc_states(void);
 static void _hoc_trates(void);
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
 "setdata_itL", _hoc_setdata,
 "rates_itL", _hoc_rates,
 "states_itL", _hoc_states,
 "trates_itL", _hoc_trates,
 0, 0
};
 /* declare global and static user variables */
#define ahh3 ahh3_itL
 double ahh3 = 46.34;
#define ahh2 ahh2_itL
 double ahh2 = 3.683;
#define ahh ahh_itL
 double ahh = 3.1;
#define ah3 ah3_itL
 double ah3 = 134;
#define ah2 ah2_itL
 double ah2 = 177;
#define ah ah_itL
 double ah = 7.6;
#define am am_itL
 double am = 1.84;
#define cai cai_itL
 double cai = 0;
#define cao cao_itL
 double cao = 2.5;
#define usetable usetable_itL
 double usetable = 1;
#define vth vth_itL
 double vth = -60;
#define vhh2 vhh2_itL
 double vhh2 = 70.6;
#define vhh1 vhh1_itL
 double vhh1 = 37.9;
#define vh2 vh2_itL
 double vh2 = 99.5;
#define vh1 vh1_itL
 double vh1 = 56.6;
#define vm2 vm2_itL
 double vm2 = 72;
#define vm1 vm1_itL
 double vm1 = 27;
#define vwh vwh_itL
 double vwh = 6;
#define vwm vwm_itL
 double vwm = 5;
#define v12h v12h_itL
 double v12h = 80;
#define v12m v12m_itL
 double v12m = 47;
#define vmax vmax_itL
 double vmax = 100;
#define vmin vmin_itL
 double vmin = -120;
#define vshift vshift_itL
 double vshift = 0;
#define whh2 whh2_itL
 double whh2 = 8.6;
#define whh1 whh1_itL
 double whh1 = 4.7;
#define wh2 wh2_itL
 double wh2 = 5.6;
#define wh1 wh1_itL
 double wh1 = 6.3;
#define wm2 wm2_itL
 double wm2 = 20;
#define wm1 wm1_itL
 double wm1 = 8;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 "usetable_itL", 0, 1,
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "vshift_itL", "mV",
 "cao_itL", "mM",
 "cai_itL", "mM",
 "vmin_itL", "mV",
 "vmax_itL", "mV",
 "v12m_itL", "mV",
 "v12h_itL", "mV",
 "vwm_itL", "mV",
 "vwh_itL", "mV",
 "am_itL", "mV",
 "ah_itL", "mV",
 "ah2_itL", "mV",
 "ah3_itL", "mV",
 "vm1_itL", "mV",
 "vm2_itL", "mV",
 "vh1_itL", "mV",
 "vh2_itL", "mV",
 "wm1_itL", "mV",
 "wm2_itL", "mV",
 "wh1_itL", "mV",
 "wh2_itL", "mV",
 "ahh_itL", "ms",
 "ahh2_itL", "ms",
 "ahh3_itL", "ms",
 "vhh1_itL", "mV",
 "vhh2_itL", "mV",
 "whh1_itL", "mV",
 "whh2_itL", "mV",
 "vth_itL", "mV",
 "gbar_itL", "mho/cm2",
 "gca_itL", "pS/um2",
 "mtau_itL", "ms",
 "htau_itL", "ms",
 0,0
};
 static double delta_t = 1;
 static double h0 = 0;
 static double m0 = 0;
 static double v = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 "vshift_itL", &vshift_itL,
 "cao_itL", &cao_itL,
 "cai_itL", &cai_itL,
 "vmin_itL", &vmin_itL,
 "vmax_itL", &vmax_itL,
 "v12m_itL", &v12m_itL,
 "v12h_itL", &v12h_itL,
 "vwm_itL", &vwm_itL,
 "vwh_itL", &vwh_itL,
 "am_itL", &am_itL,
 "ah_itL", &ah_itL,
 "ah2_itL", &ah2_itL,
 "ah3_itL", &ah3_itL,
 "vm1_itL", &vm1_itL,
 "vm2_itL", &vm2_itL,
 "vh1_itL", &vh1_itL,
 "vh2_itL", &vh2_itL,
 "wm1_itL", &wm1_itL,
 "wm2_itL", &wm2_itL,
 "wh1_itL", &wh1_itL,
 "wh2_itL", &wh2_itL,
 "ahh_itL", &ahh_itL,
 "ahh2_itL", &ahh2_itL,
 "ahh3_itL", &ahh3_itL,
 "vhh1_itL", &vhh1_itL,
 "vhh2_itL", &vhh2_itL,
 "whh1_itL", &whh1_itL,
 "whh2_itL", &whh2_itL,
 "vth_itL", &vth_itL,
 "usetable_itL", &usetable_itL,
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
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "7.5.0",
"itL",
 "gbar_itL",
 0,
 "gca_itL",
 "minf_itL",
 "hinf_itL",
 "mtau_itL",
 "htau_itL",
 0,
 "m_itL",
 "h_itL",
 0,
 0};
 static Symbol* _ca_sym;
 
extern Prop* need_memb(Symbol*);

static void nrn_alloc(Prop* _prop) {
	Prop *prop_ion;
	double *_p; Datum *_ppvar;
 	_p = nrn_prop_data_alloc(_mechtype, 14, _prop);
 	/*initialize range parameters*/
 	gbar = 0.25974;
 	_prop->param = _p;
 	_prop->param_size = 14;
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 3, _prop);
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 prop_ion = need_memb(_ca_sym);
 nrn_promote(prop_ion, 0, 1);
 	_ppvar[0]._pval = &prop_ion->param[0]; /* eca */
 	_ppvar[1]._pval = &prop_ion->param[3]; /* ica */
 	_ppvar[2]._pval = &prop_ion->param[4]; /* _ion_dicadv */
 
}
 static void _initlists();
 static void _update_ion_pointer(Datum*);
 extern Symbol* hoc_lookup(const char*);
extern void _nrn_thread_reg(int, int, void(*)(Datum*));
extern void _nrn_thread_table_reg(int, void(*)(double*, Datum*, Datum*, _NrnThread*, int));
extern void hoc_register_tolerance(int, HocStateTolerance*, Symbol***);
extern void _cvode_abstol( Symbol**, double*, int);

 void _CaT_Lyetal_reg() {
	int _vectorized = 0;
  _initlists();
 	ion_reg("ca", -10000.);
 	_ca_sym = hoc_lookup("ca_ion");
 	register_mech(_mechanism, nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init, hoc_nrnpointerindex, 0);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_setdata_reg(_mechtype, _setdata);
     _nrn_thread_reg(_mechtype, 2, _update_ion_pointer);
  hoc_register_prop_size(_mechtype, 14, 3);
  hoc_register_dparam_semantics(_mechtype, 0, "ca_ion");
  hoc_register_dparam_semantics(_mechtype, 1, "ca_ion");
  hoc_register_dparam_semantics(_mechtype, 2, "ca_ion");
 	hoc_register_cvode(_mechtype, _ode_count, 0, 0, 0);
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 itL /home/pablo/trabajo/Paris/HBP_project/DendriticSpine_N/mod/x86_64/CaT_Lyetal.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
 static double FARADAY = 96485.3;
 static double R = 8.3145;
 static double PI = 3.14159;
 static double _zmexp , _zhexp ;
 static double *_t_minf;
 static double *_t__zmexp;
 static double *_t_hinf;
 static double *_t__zhexp;
static int _reset;
static char *modelname = "";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
static int _f_trates(double);
static int rates(double);
static int states();
static int trates(double);
 static void _n_trates(double);
 
static int  states (  ) {
   trates ( _threadargscomma_ v + vshift ) ;
   m = m + _zmexp * ( minf - m ) ;
   h = h + _zhexp * ( hinf - h ) ;
   
/*VERBATIM*/
	return 0;
  return 0; }
 
static void _hoc_states(void) {
  double _r;
   _r = 1.;
 states (  );
 hoc_retpushx(_r);
}
 static double _mfac_trates, _tmin_trates;
 static void _check_trates();
 static void _check_trates() {
  static int _maktable=1; int _i, _j, _ix = 0;
  double _xi, _tmax;
  static double _sav_dt;
  if (!usetable) {return;}
  if (_sav_dt != dt) { _maktable = 1;}
  if (_maktable) { double _x, _dx; _maktable=0;
   _tmin_trates =  vmin ;
   _tmax =  vmax ;
   _dx = (_tmax - _tmin_trates)/199.; _mfac_trates = 1./_dx;
   for (_i=0, _x=_tmin_trates; _i < 200; _x += _dx, _i++) {
    _f_trates(_x);
    _t_minf[_i] = minf;
    _t__zmexp[_i] = _zmexp;
    _t_hinf[_i] = hinf;
    _t__zhexp[_i] = _zhexp;
   }
   _sav_dt = dt;
  }
 }

 static int trates(double _lv){ _check_trates();
 _n_trates(_lv);
 return 0;
 }

 static void _n_trates(double _lv){ int _i, _j;
 double _xi, _theta;
 if (!usetable) {
 _f_trates(_lv); return; 
}
 _xi = _mfac_trates * (_lv - _tmin_trates);
 if (isnan(_xi)) {
  minf = _xi;
  _zmexp = _xi;
  hinf = _xi;
  _zhexp = _xi;
  return;
 }
 if (_xi <= 0.) {
 minf = _t_minf[0];
 _zmexp = _t__zmexp[0];
 hinf = _t_hinf[0];
 _zhexp = _t__zhexp[0];
 return; }
 if (_xi >= 199.) {
 minf = _t_minf[199];
 _zmexp = _t__zmexp[199];
 hinf = _t_hinf[199];
 _zhexp = _t__zhexp[199];
 return; }
 _i = (int) _xi;
 _theta = _xi - (double)_i;
 minf = _t_minf[_i] + _theta*(_t_minf[_i+1] - _t_minf[_i]);
 _zmexp = _t__zmexp[_i] + _theta*(_t__zmexp[_i+1] - _t__zmexp[_i]);
 hinf = _t_hinf[_i] + _theta*(_t_hinf[_i+1] - _t_hinf[_i]);
 _zhexp = _t__zhexp[_i] + _theta*(_t__zhexp[_i+1] - _t__zhexp[_i]);
 }

 
static int  _f_trates (  double _lv ) {
   double _ltinc ;
 rates ( _threadargscomma_ _lv ) ;
   _ltinc = - dt ;
   _zmexp = 1.0 - exp ( _ltinc / mtau ) ;
   _zhexp = 1.0 - exp ( _ltinc / htau ) ;
    return 0; }
 
static void _hoc_trates(void) {
  double _r;
    _r = 1.;
 trates (  *getarg(1) );
 hoc_retpushx(_r);
}
 
static int  rates (  double _lv_ ) {
   double _la , _lb ;
 minf = 1.0 / ( 1.0 + exp ( - ( _lv_ + v12m ) / vwm ) ) ;
   hinf = 1.0 / ( 1.0 + exp ( ( _lv_ + v12h ) / vwh ) ) ;
   mtau = 1.0 / ( am / ( 1.0 + exp ( - ( _lv_ + vm1 ) / wm1 ) + exp ( - ( _lv_ + vm2 ) / wm2 ) ) ) ;
   htau = 1.0 / ( ah + ah2 / ( 1.0 + exp ( - ( _lv_ + vh1 ) / wh1 ) ) + ah3 / ( 1.0 + exp ( - ( _lv_ + vh2 ) / wh2 ) ) ) ;
   if ( v > vth ) {
     htau = ( ahh + ahh2 / ( 1.0 + exp ( - ( _lv_ + vhh1 ) / whh1 ) ) + ahh3 / ( 1.0 + exp ( ( _lv_ + vhh2 ) / whh2 ) ) ) ;
     }
    return 0; }
 
static void _hoc_rates(void) {
  double _r;
   _r = 1.;
 rates (  *getarg(1) );
 hoc_retpushx(_r);
}
 
static int _ode_count(int _type){ hoc_execerror("itL", "cannot be used with CVODE"); return 0;}
 extern void nrn_update_ion_pointer(Symbol*, Datum*, int, int);
 static void _update_ion_pointer(Datum* _ppvar) {
   nrn_update_ion_pointer(_ca_sym, _ppvar, 0, 0);
   nrn_update_ion_pointer(_ca_sym, _ppvar, 1, 3);
   nrn_update_ion_pointer(_ca_sym, _ppvar, 2, 4);
 }

static void initmodel() {
  int _i; double _save;_ninits++;
 _save = t;
 t = 0.0;
{
  h = h0;
  m = m0;
 {
   trates ( _threadargscomma_ v + vshift ) ;
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
  eca = _ion_eca;
 initmodel();
 }}

static double _nrn_current(double _v){double _current=0.;v=_v;{ {
   gca = gbar * m * m * m * h ;
   ica = gca * ( v - eca ) ;
   }
 _current += ica;

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
  eca = _ion_eca;
 _g = _nrn_current(_v + .001);
 	{ double _dica;
  _dica = ica;
 _rhs = _nrn_current(_v);
  _ion_dicadv += (_dica - ica)/.001 ;
 	}
 _g = (_g - _rhs)/.001;
  _ion_ica += ica ;
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
  eca = _ion_eca;
 { error =  states();
 if(error){fprintf(stderr,"at line 90 in file CaT_Lyetal.mod:\n        SOLVE states\n"); nrn_complain(_p); abort_run(error);}
 } }}

}

static void terminal(){}

static void _initlists() {
 int _i; static int _first = 1;
  if (!_first) return;
   _t_minf = makevector(200*sizeof(double));
   _t__zmexp = makevector(200*sizeof(double));
   _t_hinf = makevector(200*sizeof(double));
   _t__zhexp = makevector(200*sizeof(double));
_first = 0;
}
