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
 
#define nrn_init _nrn_init__ER
#define _nrn_initial _nrn_initial__ER
#define nrn_cur _nrn_cur__ER
#define _nrn_current _nrn_current__ER
#define nrn_jacob _nrn_jacob__ER
#define nrn_state _nrn_state__ER
#define _net_receive _net_receive__ER 
#define states states__ER 
 
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
#define ical _p[0]
#define ip3ip _p[1]
#define caer _p[2]
#define Jip3h _p[3]
#define Jip3m _p[4]
#define Jip3z _p[5]
#define ERrelz _p[6]
#define ERfilz _p[7]
#define ERlekz _p[8]
#define ier _p[9]
#define cali _p[10]
#define Dcali _p[11]
#define Dcaer _p[12]
#define DJip3h _p[13]
#define DJip3m _p[14]
#define DJip3z _p[15]
#define DERrelz _p[16]
#define DERfilz _p[17]
#define DERlekz _p[18]
#define Dier _p[19]
#define _g _p[20]
#define _ion_cali	*_ppvar[0]._pval
#define _ion_ical	*_ppvar[1]._pval
#define _ion_dicaldv	*_ppvar[2]._pval
#define ip3im	*_ppvar[3]._pval
#define _p_ip3im	_ppvar[3]._pval
#define ip3id	*_ppvar[4]._pval
#define _p_ip3id	_ppvar[4]._pval
#define diam	*_ppvar[5]._pval
 
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
 static void _hoc_Jip3Q(void);
 static void _hoc_Jip3tm(void);
 static void _hoc_Jip3minf(void);
 static void _hoc_Jip3th(void);
 static void _hoc_Jip3hinf(void);
 static void _hoc_Jip3(void);
 static void _hoc_erlek(void);
 static void _hoc_erfil(void);
 static void _hoc_errel(void);
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
 "setdata_ER", _hoc_setdata,
 "Jip3Q_ER", _hoc_Jip3Q,
 "Jip3tm_ER", _hoc_Jip3tm,
 "Jip3minf_ER", _hoc_Jip3minf,
 "Jip3th_ER", _hoc_Jip3th,
 "Jip3hinf_ER", _hoc_Jip3hinf,
 "Jip3_ER", _hoc_Jip3,
 "erlek_ER", _hoc_erlek,
 "erfil_ER", _hoc_erfil,
 "errel_ER", _hoc_errel,
 0, 0
};
#define Jip3Q Jip3Q_ER
#define Jip3tm Jip3tm_ER
#define Jip3minf Jip3minf_ER
#define Jip3th Jip3th_ER
#define Jip3hinf Jip3hinf_ER
#define Jip3 Jip3_ER
#define erlek erlek_ER
#define erfil erfil_ER
#define errel errel_ER
 extern double Jip3Q( double , double , double , double );
 extern double Jip3tm( double , double , double );
 extern double Jip3minf( double , double , double , double );
 extern double Jip3th( double , double , double , double , double , double );
 extern double Jip3hinf( double , double , double , double , double );
 extern double Jip3( double , double , double , double , double , double , double , double );
 extern double erlek( double , double , double );
 extern double erfil( double , double , double , double );
 extern double errel( double , double , double , double );
 /* declare global and static user variables */
#define Jip3m0 Jip3m0_ER
 double Jip3m0 = 0.2;
#define Jip3h0 Jip3h0_ER
 double Jip3h0 = 0.2;
#define Vip3 Vip3_ER
 double Vip3 = 1e-09;
#define ainh ainh_ER
 double ainh = 0.0002;
#define aip3 aip3_ER
 double aip3 = 420000;
#define bip3 bip3_ER
 double bip3 = 4.1;
#define caer0 caer0_ER
 double caer0 = 0.2;
#define ddis ddis_ER
 double ddis = 0.00094;
#define dip3 dip3_ER
 double dip3 = 0.00013;
#define dinh dinh_ER
 double dinh = 0.00105;
#define dact dact_ER
 double dact = 8.2e-05;
#define fer fer_ER
 double fer = 0.0025;
#define kerlek kerlek_ER
 double kerlek = 6.15e-14;
#define kerfilb kerfilb_ER
 double kerfilb = 0.0002;
#define kerfila kerfila_ER
 double kerfila = 7.5e-13;
#define kerrel kerrel_ER
 double kerrel = 3e-12;
#define kerm kerm_ER
 double kerm = 0.0002;
#define rhover rhover_ER
 double rhover = 0.15;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "caer0_ER", "mM",
 "kerm_ER", "mM",
 "kerrel_ER", "/s",
 "kerfila_ER", "mM/s",
 "kerfilb_ER", "mM",
 "kerlek_ER", "/s",
 "dact_ER", "mM",
 "dinh_ER", "mM",
 "dip3_ER", "mM",
 "ddis_ER", "mM",
 "aip3_ER", "1/mM/s",
 "bip3_ER", "1/s",
 "ainh_ER", "mM/s",
 "caer_ER", "mM",
 "Jip3h_ER", "1",
 "Jip3m_ER", "1",
 "Jip3z_ER", "mM/ms",
 "ERrelz_ER", "mM/ms",
 "ERfilz_ER", "mM/ms",
 "ERlekz_ER", "mM/ms",
 "ier_ER", "mA/cm2",
 "ical_ER", "mA/cm2",
 "ip3ip_ER", "mM",
 "ip3im_ER", "mM",
 "ip3id_ER", "mM",
 0,0
};
 static double ERlekz0 = 0;
 static double ERfilz0 = 0;
 static double ERrelz0 = 0;
 static double Jip3z0 = 0;
 static double cali0 = 0;
 static double delta_t = 0.01;
 static double ier0 = 0;
 static double v = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 "caer0_ER", &caer0_ER,
 "Jip3h0_ER", &Jip3h0_ER,
 "Jip3m0_ER", &Jip3m0_ER,
 "fer_ER", &fer_ER,
 "kerm_ER", &kerm_ER,
 "kerrel_ER", &kerrel_ER,
 "kerfila_ER", &kerfila_ER,
 "kerfilb_ER", &kerfilb_ER,
 "kerlek_ER", &kerlek_ER,
 "rhover_ER", &rhover_ER,
 "Vip3_ER", &Vip3_ER,
 "dact_ER", &dact_ER,
 "dinh_ER", &dinh_ER,
 "dip3_ER", &dip3_ER,
 "ddis_ER", &ddis_ER,
 "aip3_ER", &aip3_ER,
 "bip3_ER", &bip3_ER,
 "ainh_ER", &ainh_ER,
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
 
#define _cvode_ieq _ppvar[6]._i
 static void _ode_matsol_instance1(_threadargsproto_);
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "7.5.0",
"ER",
 0,
 "ical_ER",
 "ip3ip_ER",
 0,
 "caer_ER",
 "Jip3h_ER",
 "Jip3m_ER",
 "Jip3z_ER",
 "ERrelz_ER",
 "ERfilz_ER",
 "ERlekz_ER",
 "ier_ER",
 0,
 "ip3im_ER",
 "ip3id_ER",
 0};
 static Symbol* _morphology_sym;
 static Symbol* _cal_sym;
 
extern Prop* need_memb(Symbol*);

static void nrn_alloc(Prop* _prop) {
	Prop *prop_ion;
	double *_p; Datum *_ppvar;
 	_p = nrn_prop_data_alloc(_mechtype, 21, _prop);
 	/*initialize range parameters*/
 	_prop->param = _p;
 	_prop->param_size = 21;
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 7, _prop);
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 prop_ion = need_memb(_morphology_sym);
 	_ppvar[5]._pval = &prop_ion->param[0]; /* diam */
 prop_ion = need_memb(_cal_sym);
 nrn_promote(prop_ion, 1, 0);
 	_ppvar[0]._pval = &prop_ion->param[1]; /* cali */
 	_ppvar[1]._pval = &prop_ion->param[3]; /* ical */
 	_ppvar[2]._pval = &prop_ion->param[4]; /* _ion_dicaldv */
 
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

 void _ER_reg() {
	int _vectorized = 0;
  _initlists();
 	ion_reg("cal", -10000.);
 	_morphology_sym = hoc_lookup("morphology");
 	_cal_sym = hoc_lookup("cal_ion");
 	register_mech(_mechanism, nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init, hoc_nrnpointerindex, 0);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_setdata_reg(_mechtype, _setdata);
     _nrn_thread_reg(_mechtype, 2, _update_ion_pointer);
  hoc_register_prop_size(_mechtype, 21, 7);
  hoc_register_dparam_semantics(_mechtype, 0, "cal_ion");
  hoc_register_dparam_semantics(_mechtype, 1, "cal_ion");
  hoc_register_dparam_semantics(_mechtype, 2, "cal_ion");
  hoc_register_dparam_semantics(_mechtype, 3, "pointer");
  hoc_register_dparam_semantics(_mechtype, 4, "pointer");
  hoc_register_dparam_semantics(_mechtype, 6, "cvodeieq");
  hoc_register_dparam_semantics(_mechtype, 5, "diam");
 	hoc_register_cvode(_mechtype, _ode_count, _ode_map, _ode_spec, _ode_matsol);
 	hoc_register_tolerance(_mechtype, _hoc_state_tol, &_atollist);
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 ER /export/home1/users/bssn/serna/HBP/SpineModel/mod/x86_64/ER.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
 static double FARADAY = 96485.3;
 static double PI = 3.14159;
static int _reset;
static char *modelname = "Endoplasmic Reticulum";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
 static int _deriv1_advance = 0;
 
static int _ode_spec1(_threadargsproto_);
/*static int _ode_matsol1(_threadargsproto_);*/
 static int _slist2[3]; static double _dlist2[3];
 static double _savstate1[3], *_temp1 = _savstate1;
 static int _slist1[3], _dlist1[3];
 static int states(_threadargsproto_);
 
/*CVODE*/
 static int _ode_spec1 () {_reset=0;
 {
   Dcaer = - ( 0.001 ) * ( Jip3 ( _threadargscomma_ cali , caer , ip3ip , Vip3 , dact , dinh , dip3 , ddis ) + errel ( _threadargscomma_ cali , caer , kerrel , kerm ) - erfil ( _threadargscomma_ cali , caer , kerfila , kerfilb ) + erlek ( _threadargscomma_ cali , caer , kerlek ) ) / ( rhover / fer ) ;
   DJip3h = ( Jip3hinf ( _threadargscomma_ ip3ip , cali , dinh , dip3 , ddis ) - Jip3h ) / Jip3th ( _threadargscomma_ ip3ip , ainh , cali , dinh , dip3 , ddis ) ;
   DJip3m = ( Jip3minf ( _threadargscomma_ ip3ip , cali , dip3 , dact ) - Jip3m ) / Jip3tm ( _threadargscomma_ cali , bip3 , aip3 ) ;
   }
 return _reset;
}
 static int _ode_matsol1 () {
 Dcaer = Dcaer  / (1. - dt*( (( - ( 0.001 ) * ( Jip3 ( _threadargscomma_ cali , ( caer  + .001) , ip3ip , Vip3 , dact , dinh , dip3 , ddis ) + errel ( _threadargscomma_ cali , ( caer  + .001) , kerrel , kerm ) - erfil ( _threadargscomma_ cali , ( caer  + .001) , kerfila , kerfilb ) + erlek ( _threadargscomma_ cali , ( caer  + .001) , kerlek ) ) / ( rhover / fer ) ) - ( - ( 0.001 ) * ( Jip3 ( _threadargscomma_ cali , caer , ip3ip , Vip3 , dact , dinh , dip3 , ddis ) + errel ( _threadargscomma_ cali , caer , kerrel , kerm ) - erfil ( _threadargscomma_ cali , caer , kerfila , kerfilb ) + erlek ( _threadargscomma_ cali , caer , kerlek ) ) / ( rhover / fer )  )) / .001 )) ;
 DJip3h = DJip3h  / (1. - dt*( ( ( ( - 1.0 ) ) ) / Jip3th ( _threadargscomma_ ip3ip , ainh , cali , dinh , dip3 , ddis ) )) ;
 DJip3m = DJip3m  / (1. - dt*( ( ( ( - 1.0 ) ) ) / Jip3tm ( _threadargscomma_ cali , bip3 , aip3 ) )) ;
  return 0;
}
 /*END CVODE*/
 
static int states () {_reset=0;
 { static int _recurse = 0;
 int _counte = -1;
 if (!_recurse) {
 _recurse = 1;
 {int _id; for(_id=0; _id < 3; _id++) { _savstate1[_id] = _p[_slist1[_id]];}}
 error = newton(3,_slist2, _p, states, _dlist2);
 _recurse = 0; if(error) {abort_run(error);}}
 {
   Dcaer = - ( 0.001 ) * ( Jip3 ( _threadargscomma_ cali , caer , ip3ip , Vip3 , dact , dinh , dip3 , ddis ) + errel ( _threadargscomma_ cali , caer , kerrel , kerm ) - erfil ( _threadargscomma_ cali , caer , kerfila , kerfilb ) + erlek ( _threadargscomma_ cali , caer , kerlek ) ) / ( rhover / fer ) ;
   DJip3h = ( Jip3hinf ( _threadargscomma_ ip3ip , cali , dinh , dip3 , ddis ) - Jip3h ) / Jip3th ( _threadargscomma_ ip3ip , ainh , cali , dinh , dip3 , ddis ) ;
   DJip3m = ( Jip3minf ( _threadargscomma_ ip3ip , cali , dip3 , dact ) - Jip3m ) / Jip3tm ( _threadargscomma_ cali , bip3 , aip3 ) ;
   {int _id; for(_id=0; _id < 3; _id++) {
if (_deriv1_advance) {
 _dlist2[++_counte] = _p[_dlist1[_id]] - (_p[_slist1[_id]] - _savstate1[_id])/dt;
 }else{
_dlist2[++_counte] = _p[_slist1[_id]] - _savstate1[_id];}}}
 } }
 return _reset;}
 
double errel (  double _lcali , double _lcaer , double _lkerrel , double _lkerm ) {
   double _lerrel;
 _lerrel = _lkerrel * pow ( ( _lcali / ( _lcali + _lkerm ) ) , 1.0 ) * ( _lcaer - _lcali ) ;
   
return _lerrel;
 }
 
static void _hoc_errel(void) {
  double _r;
   _r =  errel (  *getarg(1) , *getarg(2) , *getarg(3) , *getarg(4) );
 hoc_retpushx(_r);
}
 
double erfil (  double _lcali , double _lcaer , double _lkerfila , double _lkerfilb ) {
   double _lerfil;
 _lerfil = _lkerfila * pow ( _lcali / ( 1.0 ) , 2.0 ) / ( pow ( _lcali / ( 1.0 ) , 2.0 ) + pow ( _lkerfilb / ( 1.0 ) , 2.0 ) ) ;
   
return _lerfil;
 }
 
static void _hoc_erfil(void) {
  double _r;
   _r =  erfil (  *getarg(1) , *getarg(2) , *getarg(3) , *getarg(4) );
 hoc_retpushx(_r);
}
 
double erlek (  double _lcali , double _lcaer , double _lkerlek ) {
   double _lerlek;
 _lerlek = _lkerlek * ( _lcaer - _lcali ) ;
   
return _lerlek;
 }
 
static void _hoc_erlek(void) {
  double _r;
   _r =  erlek (  *getarg(1) , *getarg(2) , *getarg(3) );
 hoc_retpushx(_r);
}
 
double Jip3 (  double _lcali , double _lcaer , double _lip3ip , double _lVip3 , double _ldact , double _ldinh , double _ldip3 , double _lddis ) {
   double _lJip3;
 _lJip3 = _lVip3 * pow ( Jip3m , 3.0 ) * pow ( Jip3h , 3.0 ) * ( _lcaer - _lcali ) ;
   
return _lJip3;
 }
 
static void _hoc_Jip3(void) {
  double _r;
   _r =  Jip3 (  *getarg(1) , *getarg(2) , *getarg(3) , *getarg(4) , *getarg(5) , *getarg(6) , *getarg(7) , *getarg(8) );
 hoc_retpushx(_r);
}
 
double Jip3minf (  double _lip3ip , double _lcali , double _ldip3 , double _ldact ) {
   double _lJip3minf;
 _lJip3minf = ( _lip3ip / ( _lip3ip + _ldip3 ) ) * ( _lcali / ( _lcali + _ldact ) ) ;
   
return _lJip3minf;
 }
 
static void _hoc_Jip3minf(void) {
  double _r;
   _r =  Jip3minf (  *getarg(1) , *getarg(2) , *getarg(3) , *getarg(4) );
 hoc_retpushx(_r);
}
 
double Jip3tm (  double _lcali , double _lbip3 , double _laip3 ) {
   double _lJip3tm;
 _lJip3tm = 1.0 / ( _lbip3 + _laip3 * _lcali ) ;
   
return _lJip3tm;
 }
 
static void _hoc_Jip3tm(void) {
  double _r;
   _r =  Jip3tm (  *getarg(1) , *getarg(2) , *getarg(3) );
 hoc_retpushx(_r);
}
 
double Jip3hinf (  double _lip3ip , double _lcali , double _ldinh , double _ldip3 , double _lddis ) {
   double _lJip3hinf;
 _lJip3hinf = Jip3Q ( _threadargscomma_ _lip3ip , _ldinh , _ldip3 , _lddis ) / ( Jip3Q ( _threadargscomma_ _lip3ip , _ldinh , _ldip3 , _lddis ) + _lcali ) ;
   
return _lJip3hinf;
 }
 
static void _hoc_Jip3hinf(void) {
  double _r;
   _r =  Jip3hinf (  *getarg(1) , *getarg(2) , *getarg(3) , *getarg(4) , *getarg(5) );
 hoc_retpushx(_r);
}
 
double Jip3th (  double _lip3ip , double _lainh , double _lcali , double _ldinh , double _ldip3 , double _lddis ) {
   double _lJip3th;
 _lJip3th = 1.0 / ( _lainh * Jip3Q ( _threadargscomma_ _lip3ip , _ldinh , _ldip3 , _lddis ) + _lcali ) ;
   
return _lJip3th;
 }
 
static void _hoc_Jip3th(void) {
  double _r;
   _r =  Jip3th (  *getarg(1) , *getarg(2) , *getarg(3) , *getarg(4) , *getarg(5) , *getarg(6) );
 hoc_retpushx(_r);
}
 
double Jip3Q (  double _lip3ip , double _ldinh , double _ldip3 , double _lddis ) {
   double _lJip3Q;
 _lJip3Q = _ldinh * ( ( _lip3ip + _ldip3 ) / ( _lip3ip + _lddis ) ) ;
   
return _lJip3Q;
 }
 
static void _hoc_Jip3Q(void) {
  double _r;
   _r =  Jip3Q (  *getarg(1) , *getarg(2) , *getarg(3) , *getarg(4) );
 hoc_retpushx(_r);
}
 
static int _ode_count(int _type){ return 3;}
 
static void _ode_spec(_NrnThread* _nt, _Memb_list* _ml, int _type) {
   Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
  cali = _ion_cali;
     _ode_spec1 ();
  }}
 
static void _ode_map(int _ieq, double** _pv, double** _pvdot, double* _pp, Datum* _ppd, double* _atol, int _type) { 
 	int _i; _p = _pp; _ppvar = _ppd;
	_cvode_ieq = _ieq;
	for (_i=0; _i < 3; ++_i) {
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
  cali = _ion_cali;
 _ode_matsol_instance1(_threadargs_);
 }}
 extern void nrn_update_ion_pointer(Symbol*, Datum*, int, int);
 static void _update_ion_pointer(Datum* _ppvar) {
   nrn_update_ion_pointer(_cal_sym, _ppvar, 0, 1);
   nrn_update_ion_pointer(_cal_sym, _ppvar, 1, 3);
   nrn_update_ion_pointer(_cal_sym, _ppvar, 2, 4);
 }

static void initmodel() {
  int _i; double _save;_ninits++;
 _save = t;
 t = 0.0;
{
  ERlekz = ERlekz0;
  ERfilz = ERfilz0;
  ERrelz = ERrelz0;
  Jip3z = Jip3z0;
  Jip3m = Jip3m0;
  Jip3h = Jip3h0;
  caer = caer0;
  ier = ier0;
 {
   caer = caer0 ;
   Jip3h = Jip3h0 ;
   Jip3m = Jip3m0 ;
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
  cali = _ion_cali;
 initmodel();
 }}

static double _nrn_current(double _v){double _current=0.;v=_v;{ {
   ip3ip = ip3im + ip3id ;
   ical = - 2.0 * FARADAY * ( 2e8 * Jip3 ( _threadargscomma_ cali , caer , ip3ip , Vip3 , dact , dinh , dip3 , ddis ) + errel ( _threadargscomma_ cali , caer , kerrel , kerm ) - erfil ( _threadargscomma_ cali , caer , kerfila , kerfilb ) + erlek ( _threadargscomma_ cali , caer , kerlek ) ) * diam * ( 1e-7 ) / 4.0 ;
   ier = ical ;
   Jip3z = Jip3 ( _threadargscomma_ cali , caer , ip3ip , Vip3 , dact , dinh , dip3 , ddis ) ;
   ERrelz = errel ( _threadargscomma_ cali , caer , kerrel , kerm ) ;
   ERfilz = erfil ( _threadargscomma_ cali , caer , kerfila , kerfilb ) ;
   ERlekz = erlek ( _threadargscomma_ cali , caer , kerlek ) ;
   }
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
  cali = _ion_cali;
 _g = _nrn_current(_v + .001);
 	{ double _dical;
  _dical = ical;
 _rhs = _nrn_current(_v);
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
double _dtsav = dt;
if (secondorder) { dt *= 0.5; }
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
  cali = _ion_cali;
 { error = _deriv1_advance = 1;
 derivimplicit(_ninits, 3, _slist1, _dlist1, _p, &t, dt, states, &_temp1);
_deriv1_advance = 0;
 if(error){fprintf(stderr,"at line 83 in file ER.mod:\n	  \n"); nrn_complain(_p); abort_run(error);}
    if (secondorder) {
    int _i;
    for (_i = 0; _i < 3; ++_i) {
      _p[_slist1[_i]] += dt*_p[_dlist1[_i]];
    }}
 } }}
 dt = _dtsav;
}

static void terminal(){}

static void _initlists() {
 int _i; static int _first = 1;
  if (!_first) return;
 _slist1[0] = &(caer) - _p;  _dlist1[0] = &(Dcaer) - _p;
 _slist1[1] = &(Jip3h) - _p;  _dlist1[1] = &(DJip3h) - _p;
 _slist1[2] = &(Jip3m) - _p;  _dlist1[2] = &(DJip3m) - _p;
 _slist2[0] = &(Jip3m) - _p;
 _slist2[1] = &(Jip3h) - _p;
 _slist2[2] = &(caer) - _p;
_first = 0;
}
