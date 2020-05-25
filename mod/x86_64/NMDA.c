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
 
#define nrn_init _nrn_init__NMDA
#define _nrn_initial _nrn_initial__NMDA
#define nrn_cur _nrn_cur__NMDA
#define _nrn_current _nrn_current__NMDA
#define nrn_jacob _nrn_jacob__NMDA
#define nrn_state _nrn_state__NMDA
#define _net_receive _net_receive__NMDA 
#define _f_mgblock _f_mgblock__NMDA 
#define betadyn betadyn__NMDA 
#define mgblock mgblock__NMDA 
 
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
#define tau_r _p[1]
#define tau_d _p[2]
#define qfact _p[3]
#define mg _p[4]
#define ca_ratio _p[5]
#define itmp _p[6]
#define i _p[7]
#define ical _p[8]
#define t1 _p[9]
#define spkcnt _p[10]
#define countflag _p[11]
#define scale _p[12]
#define y1 _p[13]
#define y2 _p[14]
#define g _p[15]
#define y1_add _p[16]
#define y1_loc _p[17]
#define B _p[18]
#define Dy1 _p[19]
#define Dy2 _p[20]
#define _g _p[21]
#define _tsav _p[22]
#define _nd_area  *_ppvar[0]._pval
#define _ion_ical	*_ppvar[2]._pval
#define _ion_dicaldv	*_ppvar[3]._pval
#define mu	*_ppvar[4]._pval
#define _p_mu	_ppvar[4]._pval
 
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
 static int hoc_nrnpointerindex =  4;
 /* external NEURON variables */
 /* declaration of user functions */
 static double _hoc_mgblock();
 static int _mechtype;
extern void _nrn_cacheloop_reg(int, int);
extern void hoc_register_prop_size(int, int, int);
extern void hoc_register_limits(int, HocParmLimits*);
extern void hoc_register_units(int, HocParmUnits*);
extern void nrn_promote(Prop*, int, int);
extern Memb_func* memb_func;
 extern Prop* nrn_point_prop_;
 static int _pointtype;
 static void* _hoc_create_pnt(_ho) Object* _ho; { void* create_point_process();
 return create_point_process(_pointtype, _ho);
}
 static void _hoc_destroy_pnt();
 static double _hoc_loc_pnt(_vptr) void* _vptr; {double loc_point_process();
 return loc_point_process(_pointtype, _vptr);
}
 static double _hoc_has_loc(_vptr) void* _vptr; {double has_loc_point();
 return has_loc_point(_vptr);
}
 static double _hoc_get_loc_pnt(_vptr)void* _vptr; {
 double get_loc_point_process(); return (get_loc_point_process(_vptr));
}
 extern void _nrn_setdata_reg(int, void(*)(Prop*));
 static void _setdata(Prop* _prop) {
 _p = _prop->param; _ppvar = _prop->dparam;
 }
 static void _hoc_setdata(void* _vptr) { Prop* _prop;
 _prop = ((Point_process*)_vptr)->_prop;
   _setdata(_prop);
 }
 /* connect user functions to hoc names */
 static VoidFunc hoc_intfunc[] = {
 0,0
};
 static Member_func _member_func[] = {
 "loc", _hoc_loc_pnt,
 "has_loc", _hoc_has_loc,
 "get_loc", _hoc_get_loc_pnt,
 "mgblock", _hoc_mgblock,
 0, 0
};
 /* declare global and static user variables */
#define Erev Erev_NMDA
 double Erev = 0;
#define saturation saturation_NMDA
 double saturation = 7;
#define usetable usetable_NMDA
 double usetable = 1;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 "usetable_NMDA", 0, 1,
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "Erev_NMDA", "mV",
 "gbar", "umho",
 "tau_r", "ms",
 "tau_d", "ms",
 "mg", "mM",
 "y1", "/ms",
 "itmp", "nA",
 "i", "nA",
 "ical", "nA",
 "t1", "ms",
 "mu", "1",
 0,0
};
 static double delta_t = 0.01;
 static double v = 0;
 static double y20 = 0;
 static double y10 = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 "Erev_NMDA", &Erev_NMDA,
 "saturation_NMDA", &saturation_NMDA,
 "usetable_NMDA", &usetable_NMDA,
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
 static void _hoc_destroy_pnt(_vptr) void* _vptr; {
   destroy_point_process(_vptr);
}
 
static int _ode_count(int);
static void _ode_map(int, double**, double**, double*, Datum*, double*, int);
static void _ode_spec(_NrnThread*, _Memb_list*, int);
static void _ode_matsol(_NrnThread*, _Memb_list*, int);
 
#define _cvode_ieq _ppvar[5]._i
 static void _ode_matsol_instance1(_threadargsproto_);
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "7.5.0",
"NMDA",
 "gbar",
 "tau_r",
 "tau_d",
 "qfact",
 "mg",
 "ca_ratio",
 0,
 "itmp",
 "i",
 "ical",
 "t1",
 "spkcnt",
 "countflag",
 "scale",
 0,
 "y1",
 "y2",
 0,
 "mu",
 0};
 static Symbol* _cal_sym;
 
extern Prop* need_memb(Symbol*);

static void nrn_alloc(Prop* _prop) {
	Prop *prop_ion;
	double *_p; Datum *_ppvar;
  if (nrn_point_prop_) {
	_prop->_alloc_seq = nrn_point_prop_->_alloc_seq;
	_p = nrn_point_prop_->param;
	_ppvar = nrn_point_prop_->dparam;
 }else{
 	_p = nrn_prop_data_alloc(_mechtype, 23, _prop);
 	/*initialize range parameters*/
 	gbar = 0.000122;
 	tau_r = 5.63;
 	tau_d = 320;
 	qfact = 2;
 	mg = 1;
 	ca_ratio = 0.01;
  }
 	_prop->param = _p;
 	_prop->param_size = 23;
  if (!nrn_point_prop_) {
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 6, _prop);
  }
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 prop_ion = need_memb(_cal_sym);
 	_ppvar[2]._pval = &prop_ion->param[3]; /* ical */
 	_ppvar[3]._pval = &prop_ion->param[4]; /* _ion_dicaldv */
 
}
 static void _initlists();
  /* some states have an absolute tolerance */
 static Symbol** _atollist;
 static HocStateTolerance _hoc_state_tol[] = {
 0,0
};
 static void _net_receive(Point_process*, double*, double);
 static void _update_ion_pointer(Datum*);
 extern Symbol* hoc_lookup(const char*);
extern void _nrn_thread_reg(int, int, void(*)(Datum*));
extern void _nrn_thread_table_reg(int, void(*)(double*, Datum*, Datum*, _NrnThread*, int));
extern void hoc_register_tolerance(int, HocStateTolerance*, Symbol***);
extern void _cvode_abstol( Symbol**, double*, int);

 void _NMDA_reg() {
	int _vectorized = 0;
  _initlists();
 	ion_reg("cal", 2.0);
 	_cal_sym = hoc_lookup("cal_ion");
 	_pointtype = point_register_mech(_mechanism,
	 nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init,
	 hoc_nrnpointerindex, 0,
	 _hoc_create_pnt, _hoc_destroy_pnt, _member_func);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_setdata_reg(_mechtype, _setdata);
     _nrn_thread_reg(_mechtype, 2, _update_ion_pointer);
  hoc_register_prop_size(_mechtype, 23, 6);
  hoc_register_dparam_semantics(_mechtype, 0, "area");
  hoc_register_dparam_semantics(_mechtype, 1, "pntproc");
  hoc_register_dparam_semantics(_mechtype, 2, "cal_ion");
  hoc_register_dparam_semantics(_mechtype, 3, "cal_ion");
  hoc_register_dparam_semantics(_mechtype, 4, "pointer");
  hoc_register_dparam_semantics(_mechtype, 5, "cvodeieq");
 	hoc_register_cvode(_mechtype, _ode_count, _ode_map, _ode_spec, _ode_matsol);
 	hoc_register_tolerance(_mechtype, _hoc_state_tol, &_atollist);
 pnt_receive[_mechtype] = _net_receive;
 pnt_receive_size[_mechtype] = 2;
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 NMDA /export/home1/users/bssn/serna/HBP/SpineModel/mod/x86_64/NMDA.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
 static double *_t_B;
static int _reset;
static char *modelname = "NMDA synapse for nucleus accumbens model";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
static int _f_mgblock(double);
static int mgblock(double);
 
static int _ode_spec1(_threadargsproto_);
/*static int _ode_matsol1(_threadargsproto_);*/
 static void _n_mgblock(double);
 static int _slist1[2], _dlist1[2];
 static int betadyn(_threadargsproto_);
 
/*CVODE*/
 static int _ode_spec1 () {_reset=0;
 {
   Dy1 = - y1 / ( tau_d / qfact ) ;
   Dy2 = y1 - y2 / ( tau_r / qfact ) ;
   }
 return _reset;
}
 static int _ode_matsol1 () {
 Dy1 = Dy1  / (1. - dt*( ( - 1.0 ) / ( tau_d / qfact ) )) ;
 Dy2 = Dy2  / (1. - dt*( ( - ( 1.0 ) / ( tau_r / qfact ) ) )) ;
  return 0;
}
 /*END CVODE*/
 static int betadyn () {_reset=0;
 {
    y1 = y1 + (1. - exp(dt*(( - 1.0 ) / ( tau_d / qfact ))))*(- ( 0.0 ) / ( ( - 1.0 ) / ( tau_d / qfact ) ) - y1) ;
    y2 = y2 + (1. - exp(dt*(( - ( 1.0 ) / ( tau_r / qfact ) ))))*(- ( y1 ) / ( ( - ( 1.0 ) / ( tau_r / qfact ) ) ) - y2) ;
   }
  return 0;
}
 
static void _net_receive (_pnt, _args, _lflag) Point_process* _pnt; double* _args; double _lflag; 
{    _p = _pnt->_prop->param; _ppvar = _pnt->_prop->dparam;
  if (_tsav > t){ extern char* hoc_object_name(); hoc_execerror(hoc_object_name(_pnt->ob), ":Event arrived out of order. Must call ParallelContext.set_maxstep AFTER assigning minimum NetCon.delay");}
 _tsav = t; {
   _args[1] = _args[1] * exp ( - ( t - t1 ) / ( tau_d / qfact ) ) ;
   y1_add = ( 1.0 - _args[1] / saturation ) ;
   _args[1] = _args[1] + y1_add ;
     if (nrn_netrec_state_adjust && !cvode_active_){
    /* discon state adjustment for cnexp case (rate uses no local variable) */
    double __state = y1;
    double __primary = (y1 + y1_add) - __state;
     __primary += ( 1. - exp( 0.5*dt*( ( - 1.0 ) / ( tau_d / qfact ) ) ) )*( - ( 0.0 ) / ( ( - 1.0 ) / ( tau_d / qfact ) ) - __primary );
    y1 += __primary;
  } else {
 y1 = y1 + y1_add ;
     }
 t1 = t ;
   spkcnt = spkcnt + 1.0 ;
   scale = _args[0] * ( ( mu - 1.0 ) * 0.3 + 1.0 ) ;
   } }
 static double _mfac_mgblock, _tmin_mgblock;
 static void _check_mgblock();
 static void _check_mgblock() {
  static int _maktable=1; int _i, _j, _ix = 0;
  double _xi, _tmax;
  static double _sav_mg;
  if (!usetable) {return;}
  if (_sav_mg != mg) { _maktable = 1;}
  if (_maktable) { double _x, _dx; _maktable=0;
   _tmin_mgblock =  - 100.0 ;
   _tmax =  100.0 ;
   _dx = (_tmax - _tmin_mgblock)/201.; _mfac_mgblock = 1./_dx;
   for (_i=0, _x=_tmin_mgblock; _i < 202; _x += _dx, _i++) {
    _f_mgblock(_x);
    _t_B[_i] = B;
   }
   _sav_mg = mg;
  }
 }

 static int mgblock(double _lv){ _check_mgblock();
 _n_mgblock(_lv);
 return 0;
 }

 static void _n_mgblock(double _lv){ int _i, _j;
 double _xi, _theta;
 if (!usetable) {
 _f_mgblock(_lv); return; 
}
 _xi = _mfac_mgblock * (_lv - _tmin_mgblock);
 if (isnan(_xi)) {
  B = _xi;
  return;
 }
 if (_xi <= 0.) {
 B = _t_B[0];
 return; }
 if (_xi >= 201.) {
 B = _t_B[201];
 return; }
 _i = (int) _xi;
 _theta = _xi - (double)_i;
 B = _t_B[_i] + _theta*(_t_B[_i+1] - _t_B[_i]);
 }

 
static int  _f_mgblock (  double _lv ) {
   B = 1.0 / ( 1.0 + exp ( 0.062 * - _lv ) * ( mg / 3.57 ) ) ;
    return 0; }
 
static double _hoc_mgblock(void* _vptr) {
 double _r;
    _hoc_setdata(_vptr);
  _r = 1.;
 mgblock (  *getarg(1) );
 return(_r);
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
 _ode_matsol_instance1(_threadargs_);
 }}
 extern void nrn_update_ion_pointer(Symbol*, Datum*, int, int);
 static void _update_ion_pointer(Datum* _ppvar) {
   nrn_update_ion_pointer(_cal_sym, _ppvar, 2, 3);
   nrn_update_ion_pointer(_cal_sym, _ppvar, 3, 4);
 }

static void initmodel() {
  int _i; double _save;_ninits++;
 _save = t;
 t = 0.0;
{
  y2 = y20;
  y1 = y10;
 {
   y1_add = 0.0 ;
   B = mgblock ( _threadargscomma_ v ) ;
   scale = 1.0 ;
   spkcnt = 0.0 ;
   countflag = 0.0 ;
   t1 = 0.0 ;
   y1_loc = 0.0 ;
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
 _tsav = -1e20;
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
 initmodel();
 }}

static double _nrn_current(double _v){double _current=0.;v=_v;{ {
   mgblock ( _threadargscomma_ v ) ;
   g = gbar * y2 ;
   itmp = scale * g * B * ( v - Erev ) ;
   i = ( 1.0 - ca_ratio ) * itmp ;
   ical = ca_ratio * itmp ;
   }
 _current += i;
 _current += ical;

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
 _g = _nrn_current(_v + .001);
 	{ double _dical;
  _dical = ical;
 _rhs = _nrn_current(_v);
  _ion_dicaldv += (_dical - ical)/.001 * 1.e2/ (_nd_area);
 	}
 _g = (_g - _rhs)/.001;
  _ion_ical += ical * 1.e2/ (_nd_area);
 _g *=  1.e2/(_nd_area);
 _rhs *= 1.e2/(_nd_area);
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
 { error =  betadyn();
 if(error){fprintf(stderr,"at line 71 in file NMDA.mod:\n	SOLVE betadyn METHOD cnexp\n"); nrn_complain(_p); abort_run(error);}
 } }}

}

static void terminal(){}

static void _initlists() {
 int _i; static int _first = 1;
  if (!_first) return;
 _slist1[0] = &(y1) - _p;  _dlist1[0] = &(Dy1) - _p;
 _slist1[1] = &(y2) - _p;  _dlist1[1] = &(Dy2) - _p;
   _t_B = makevector(202*sizeof(double));
_first = 0;
}
