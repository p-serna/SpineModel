/* Created by Language version: 7.5.0 */
/* VECTORIZED */
#define NRN_VECTORIZED 1
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
 
#define nrn_init _nrn_init__cat
#define _nrn_initial _nrn_initial__cat
#define nrn_cur _nrn_cur__cat
#define _nrn_current _nrn_current__cat
#define nrn_jacob _nrn_jacob__cat
#define nrn_state _nrn_state__cat
#define _net_receive _net_receive__cat 
#define _f_settables _f_settables__cat 
#define settables settables__cat 
#define states states__cat 
 
#define _threadargscomma_ _p, _ppvar, _thread, _nt,
#define _threadargsprotocomma_ double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt,
#define _threadargs_ _p, _ppvar, _thread, _nt
#define _threadargsproto_ double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt
 	/*SUPPRESS 761*/
	/*SUPPRESS 762*/
	/*SUPPRESS 763*/
	/*SUPPRESS 765*/
	 extern double *getarg();
 /* Thread safe. No static _p or _ppvar. */
 
#define t _nt->_t
#define dt _nt->_dt
#define pcatbar _p[0]
#define ical _p[1]
#define m _p[2]
#define h _p[3]
#define ecal _p[4]
#define cali _p[5]
#define calo _p[6]
#define minf _p[7]
#define hinf _p[8]
#define Dm _p[9]
#define Dh _p[10]
#define v _p[11]
#define _g _p[12]
#define _ion_cali	*_ppvar[0]._pval
#define _ion_calo	*_ppvar[1]._pval
#define _ion_ical	*_ppvar[2]._pval
#define _ion_dicaldv	*_ppvar[3]._pval
 
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
 static Datum* _extcall_thread;
 static Prop* _extcall_prop;
 /* external NEURON variables */
 extern double celsius;
 /* declaration of user functions */
 static void _hoc_efun(void);
 static void _hoc_ghk(void);
 static void _hoc_htau(void);
 static void _hoc_mtau(void);
 static void _hoc_settables(void);
 static void _hoc_table_htau(void);
 static void _hoc_table_mtau(void);
 static int _mechtype;
extern void _nrn_cacheloop_reg(int, int);
extern void hoc_register_prop_size(int, int, int);
extern void hoc_register_limits(int, HocParmLimits*);
extern void hoc_register_units(int, HocParmUnits*);
extern void nrn_promote(Prop*, int, int);
extern Memb_func* memb_func;
 extern void _nrn_setdata_reg(int, void(*)(Prop*));
 static void _setdata(Prop* _prop) {
 _extcall_prop = _prop;
 }
 static void _hoc_setdata() {
 Prop *_prop, *hoc_getdata_range(int);
 _prop = hoc_getdata_range(_mechtype);
   _setdata(_prop);
 hoc_retpushx(1.);
}
 /* connect user functions to hoc names */
 static VoidFunc hoc_intfunc[] = {
 "setdata_cat", _hoc_setdata,
 "efun_cat", _hoc_efun,
 "ghk_cat", _hoc_ghk,
 "htau_cat", _hoc_htau,
 "mtau_cat", _hoc_mtau,
 "settables_cat", _hoc_settables,
 "table_htau_cat", _hoc_table_htau,
 "table_mtau_cat", _hoc_table_mtau,
 0, 0
};
#define efun efun_cat
#define ghk ghk_cat
#define htau htau_cat
#define mtau mtau_cat
#define table_htau table_htau_cat
#define table_mtau table_mtau_cat
 extern double efun( _threadargsprotocomma_ double );
 extern double ghk( _threadargsprotocomma_ double , double , double );
 extern double htau( _threadargsprotocomma_ double );
 extern double mtau( _threadargsprotocomma_ double );
 extern double table_htau( );
 extern double table_mtau( );
 
static void _check_settables(double*, Datum*, Datum*, _NrnThread*); 
static void _check_table_thread(double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt, int _type) {
   _check_settables(_p, _ppvar, _thread, _nt);
 }
 /* declare global and static user variables */
#define hshift hshift_cat
 double hshift = 0;
#define hslope hslope_cat
 double hslope = 6.7;
#define hvhalf hvhalf_cat
 double hvhalf = -80;
#define mshift mshift_cat
 double mshift = 0;
#define mslope mslope_cat
 double mslope = -6.53;
#define mvhalf mvhalf_cat
 double mvhalf = -51.73;
#define qfact qfact_cat
 double qfact = 3;
#define usetable usetable_cat
 double usetable = 1;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 "usetable_cat", 0, 1,
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "mvhalf_cat", "mV",
 "mslope_cat", "mV",
 "mshift_cat", "mV",
 "hvhalf_cat", "mV",
 "hslope_cat", "mV",
 "hshift_cat", "mV",
 "pcatbar_cat", "cm/s",
 "ical_cat", "mA/cm2",
 0,0
};
 static double delta_t = 0.01;
 static double h0 = 0;
 static double m0 = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 "mvhalf_cat", &mvhalf_cat,
 "mslope_cat", &mslope_cat,
 "mshift_cat", &mshift_cat,
 "hvhalf_cat", &hvhalf_cat,
 "hslope_cat", &hslope_cat,
 "hshift_cat", &hshift_cat,
 "qfact_cat", &qfact_cat,
 "usetable_cat", &usetable_cat,
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
"cat",
 "pcatbar_cat",
 0,
 "ical_cat",
 0,
 "m_cat",
 "h_cat",
 0,
 0};
 static Symbol* _cal_sym;
 
extern Prop* need_memb(Symbol*);

static void nrn_alloc(Prop* _prop) {
	Prop *prop_ion;
	double *_p; Datum *_ppvar;
 	_p = nrn_prop_data_alloc(_mechtype, 13, _prop);
 	/*initialize range parameters*/
 	pcatbar = 7.6e-07;
 	_prop->param = _p;
 	_prop->param_size = 13;
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 5, _prop);
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 prop_ion = need_memb(_cal_sym);
 nrn_promote(prop_ion, 1, 0);
 	_ppvar[0]._pval = &prop_ion->param[1]; /* cali */
 	_ppvar[1]._pval = &prop_ion->param[2]; /* calo */
 	_ppvar[2]._pval = &prop_ion->param[3]; /* ical */
 	_ppvar[3]._pval = &prop_ion->param[4]; /* _ion_dicaldv */
 
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

 void _cat_reg() {
	int _vectorized = 1;
  _initlists();
 	ion_reg("cal", 2.0);
 	_cal_sym = hoc_lookup("cal_ion");
 	register_mech(_mechanism, nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init, hoc_nrnpointerindex, 1);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_setdata_reg(_mechtype, _setdata);
     _nrn_thread_reg(_mechtype, 2, _update_ion_pointer);
     _nrn_thread_table_reg(_mechtype, _check_table_thread);
  hoc_register_prop_size(_mechtype, 13, 5);
  hoc_register_dparam_semantics(_mechtype, 0, "cal_ion");
  hoc_register_dparam_semantics(_mechtype, 1, "cal_ion");
  hoc_register_dparam_semantics(_mechtype, 2, "cal_ion");
  hoc_register_dparam_semantics(_mechtype, 3, "cal_ion");
  hoc_register_dparam_semantics(_mechtype, 4, "cvodeieq");
 	hoc_register_cvode(_mechtype, _ode_count, _ode_map, _ode_spec, _ode_matsol);
 	hoc_register_tolerance(_mechtype, _hoc_state_tol, &_atollist);
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 cat /export/home1/users/bssn/serna/HBP/SpineModel/mod/x86_64/cat.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
 static double FARADAY = 96485.3;
 static double R = 8.3145;
 static double *_t_minf;
 static double *_t_hinf;
static int _reset;
static char *modelname = "T-type calcium channel for nucleus accumbens neuron ";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
static int _f_settables(_threadargsprotocomma_ double);
static int settables(_threadargsprotocomma_ double);
 
static int _ode_spec1(_threadargsproto_);
/*static int _ode_matsol1(_threadargsproto_);*/
 static void _n_settables(_threadargsprotocomma_ double _lv);
 
static void* _ptable_htau = (void*)0;
 
static void* _ptable_mtau = (void*)0;
 static int _slist1[2], _dlist1[2];
 static int states(_threadargsproto_);
 
/*CVODE*/
 static int _ode_spec1 (double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {int _reset = 0; {
   settables ( _threadargscomma_ v ) ;
   Dm = ( minf - m ) / ( mtau ( _threadargscomma_ v ) / qfact ) ;
   Dh = ( hinf - h ) / ( htau ( _threadargscomma_ v ) / qfact ) ;
   }
 return _reset;
}
 static int _ode_matsol1 (double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {
 settables ( _threadargscomma_ v ) ;
 Dm = Dm  / (1. - dt*( ( ( ( - 1.0 ) ) ) / ( mtau ( _threadargscomma_ v ) / qfact ) )) ;
 Dh = Dh  / (1. - dt*( ( ( ( - 1.0 ) ) ) / ( htau ( _threadargscomma_ v ) / qfact ) )) ;
  return 0;
}
 /*END CVODE*/
 static int states (double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) { {
   settables ( _threadargscomma_ v ) ;
    m = m + (1. - exp(dt*(( ( ( - 1.0 ) ) ) / ( mtau ( _threadargscomma_ v ) / qfact ))))*(- ( ( ( minf ) ) / ( mtau ( _threadargscomma_ v ) / qfact ) ) / ( ( ( ( - 1.0 ) ) ) / ( mtau ( _threadargscomma_ v ) / qfact ) ) - m) ;
    h = h + (1. - exp(dt*(( ( ( - 1.0 ) ) ) / ( htau ( _threadargscomma_ v ) / qfact ))))*(- ( ( ( hinf ) ) / ( htau ( _threadargscomma_ v ) / qfact ) ) / ( ( ( ( - 1.0 ) ) ) / ( htau ( _threadargscomma_ v ) / qfact ) ) - h) ;
   }
  return 0;
}
 
double mtau ( _threadargsprotocomma_ double _lv ) {
 double _arg[1];
 _arg[0] = _lv;
 return hoc_func_table(_ptable_mtau, 1, _arg);
 }
/*  }
  */
 
static void _hoc_mtau(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r =  mtau ( _p, _ppvar, _thread, _nt, *getarg(1) );
 hoc_retpushx(_r);
}
 double table_mtau ( ) {
	hoc_spec_table(&_ptable_mtau, 1);
	return 0.;
}
 
static void _hoc_table_mtau(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r =  table_mtau (  );
 hoc_retpushx(_r);
}
 
double htau ( _threadargsprotocomma_ double _lv ) {
 double _arg[1];
 _arg[0] = _lv;
 return hoc_func_table(_ptable_htau, 1, _arg);
 }
/*  }
  */
 
static void _hoc_htau(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r =  htau ( _p, _ppvar, _thread, _nt, *getarg(1) );
 hoc_retpushx(_r);
}
 double table_htau ( ) {
	hoc_spec_table(&_ptable_htau, 1);
	return 0.;
}
 
static void _hoc_table_htau(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r =  table_htau (  );
 hoc_retpushx(_r);
}
 static double _mfac_settables, _tmin_settables;
  static void _check_settables(double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {
  static int _maktable=1; int _i, _j, _ix = 0;
  double _xi, _tmax;
  static double _sav_mshift;
  static double _sav_hshift;
  if (!usetable) {return;}
  if (_sav_mshift != mshift) { _maktable = 1;}
  if (_sav_hshift != hshift) { _maktable = 1;}
  if (_maktable) { double _x, _dx; _maktable=0;
   _tmin_settables =  - 100.0 ;
   _tmax =  100.0 ;
   _dx = (_tmax - _tmin_settables)/201.; _mfac_settables = 1./_dx;
   for (_i=0, _x=_tmin_settables; _i < 202; _x += _dx, _i++) {
    _f_settables(_p, _ppvar, _thread, _nt, _x);
    _t_minf[_i] = minf;
    _t_hinf[_i] = hinf;
   }
   _sav_mshift = mshift;
   _sav_hshift = hshift;
  }
 }

 static int settables(double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt, double _lv) { 
#if 0
_check_settables(_p, _ppvar, _thread, _nt);
#endif
 _n_settables(_p, _ppvar, _thread, _nt, _lv);
 return 0;
 }

 static void _n_settables(double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt, double _lv){ int _i, _j;
 double _xi, _theta;
 if (!usetable) {
 _f_settables(_p, _ppvar, _thread, _nt, _lv); return; 
}
 _xi = _mfac_settables * (_lv - _tmin_settables);
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

 
static int  _f_settables ( _threadargsprotocomma_ double _lv ) {
   minf = 1.0 / ( 1.0 + exp ( ( _lv - mvhalf - mshift ) / mslope ) ) ;
   hinf = 1.0 / ( 1.0 + exp ( ( _lv - hvhalf - hshift ) / hslope ) ) ;
    return 0; }
 
static void _hoc_settables(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 
#if 1
 _check_settables(_p, _ppvar, _thread, _nt);
#endif
 _r = 1.;
 settables ( _p, _ppvar, _thread, _nt, *getarg(1) );
 hoc_retpushx(_r);
}
 
double ghk ( _threadargsprotocomma_ double _lv , double _lci , double _lco ) {
   double _lghk;
 double _lz , _leci , _leco ;
 _lz = ( 1e-3 ) * 2.0 * FARADAY * _lv / ( R * ( celsius + 273.15 ) ) ;
   _leco = _lco * efun ( _threadargscomma_ _lz ) ;
   _leci = _lci * efun ( _threadargscomma_ - _lz ) ;
   _lghk = ( .001 ) * 2.0 * FARADAY * ( _leci - _leco ) ;
   
return _lghk;
 }
 
static void _hoc_ghk(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r =  ghk ( _p, _ppvar, _thread, _nt, *getarg(1) , *getarg(2) , *getarg(3) );
 hoc_retpushx(_r);
}
 
double efun ( _threadargsprotocomma_ double _lz ) {
   double _lefun;
 if ( fabs ( _lz ) < 1e-4 ) {
     _lefun = 1.0 - _lz / 2.0 ;
     }
   else {
     _lefun = _lz / ( exp ( _lz ) - 1.0 ) ;
     }
   
return _lefun;
 }
 
static void _hoc_efun(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r =  efun ( _p, _ppvar, _thread, _nt, *getarg(1) );
 hoc_retpushx(_r);
}
 
static int _ode_count(int _type){ return 2;}
 
static void _ode_spec(_NrnThread* _nt, _Memb_list* _ml, int _type) {
   double* _p; Datum* _ppvar; Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
  cali = _ion_cali;
  calo = _ion_calo;
     _ode_spec1 (_p, _ppvar, _thread, _nt);
  }}
 
static void _ode_map(int _ieq, double** _pv, double** _pvdot, double* _pp, Datum* _ppd, double* _atol, int _type) { 
	double* _p; Datum* _ppvar;
 	int _i; _p = _pp; _ppvar = _ppd;
	_cvode_ieq = _ieq;
	for (_i=0; _i < 2; ++_i) {
		_pv[_i] = _pp + _slist1[_i];  _pvdot[_i] = _pp + _dlist1[_i];
		_cvode_abstol(_atollist, _atol, _i);
	}
 }
 
static void _ode_matsol_instance1(_threadargsproto_) {
 _ode_matsol1 (_p, _ppvar, _thread, _nt);
 }
 
static void _ode_matsol(_NrnThread* _nt, _Memb_list* _ml, int _type) {
   double* _p; Datum* _ppvar; Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
  cali = _ion_cali;
  calo = _ion_calo;
 _ode_matsol_instance1(_threadargs_);
 }}
 extern void nrn_update_ion_pointer(Symbol*, Datum*, int, int);
 static void _update_ion_pointer(Datum* _ppvar) {
   nrn_update_ion_pointer(_cal_sym, _ppvar, 0, 1);
   nrn_update_ion_pointer(_cal_sym, _ppvar, 1, 2);
   nrn_update_ion_pointer(_cal_sym, _ppvar, 2, 3);
   nrn_update_ion_pointer(_cal_sym, _ppvar, 3, 4);
 }

static void initmodel(double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {
  int _i; double _save;{
  h = h0;
  m = m0;
 {
   settables ( _threadargscomma_ v ) ;
   m = minf ;
   h = 0.76 ;
   }
 
}
}

static void nrn_init(_NrnThread* _nt, _Memb_list* _ml, int _type){
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; double _v; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];

#if 0
 _check_settables(_p, _ppvar, _thread, _nt);
#endif
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
  cali = _ion_cali;
  calo = _ion_calo;
 initmodel(_p, _ppvar, _thread, _nt);
 }
}

static double _nrn_current(double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt, double _v){double _current=0.;v=_v;{ {
   ical = ghk ( _threadargscomma_ v , cali , calo ) * pcatbar * m * m * m * h ;
   }
 _current += ical;

} return _current;
}

static void nrn_cur(_NrnThread* _nt, _Memb_list* _ml, int _type) {
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; int* _ni; double _rhs, _v; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
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
  cali = _ion_cali;
  calo = _ion_calo;
 _g = _nrn_current(_p, _ppvar, _thread, _nt, _v + .001);
 	{ double _dical;
  _dical = ical;
 _rhs = _nrn_current(_p, _ppvar, _thread, _nt, _v);
  _ion_dicaldv += (_dical - ical)/.001 ;
 	}
 _g = (_g - _rhs)/.001;
  _ion_ical += ical ;
#if CACHEVEC
  if (use_cachevec) {
	VEC_RHS(_ni[_iml]) -= _rhs;
  }else
#endif
  {
	NODERHS(_nd) -= _rhs;
  }
 
}
 
}

static void nrn_jacob(_NrnThread* _nt, _Memb_list* _ml, int _type) {
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
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
 
}
 
}

static void nrn_state(_NrnThread* _nt, _Memb_list* _ml, int _type) {
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; double _v = 0.0; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
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
  cali = _ion_cali;
  calo = _ion_calo;
 {   states(_p, _ppvar, _thread, _nt);
  } }}

}

static void terminal(){}

static void _initlists(){
 double _x; double* _p = &_x;
 int _i; static int _first = 1;
  if (!_first) return;
 _slist1[0] = &(m) - _p;  _dlist1[0] = &(Dm) - _p;
 _slist1[1] = &(h) - _p;  _dlist1[1] = &(Dh) - _p;
   _t_minf = makevector(202*sizeof(double));
   _t_hinf = makevector(202*sizeof(double));
_first = 0;
}

#if defined(__cplusplus)
} /* extern "C" */
#endif
