#include <stdio.h>
#include "hocdec.h"
extern int nrnmpi_myid;
extern int nrn_nobanner_;

extern void _AMPA_reg(void);
extern void _bkkca_reg(void);
extern void _cad_reg(void);
extern void _cadyn_reg(void);
extern void _caL13_reg(void);
extern void _caL13PS_reg(void);
extern void _caldyn_reg(void);
extern void _caL_reg(void);
extern void _caLPS_reg(void);
extern void _can_reg(void);
extern void _canPS_reg(void);
extern void _caq_reg(void);
extern void _caqPS_reg(void);
extern void _car_reg(void);
extern void _carPS_reg(void);
extern void _cat_reg(void);
extern void _damsg_reg(void);
extern void _ER_reg(void);
extern void _exp2synNMDA_reg(void);
extern void _GABA_reg(void);
extern void _kaf_reg(void);
extern void _kas_reg(void);
extern void _kir_reg(void);
extern void _krp_reg(void);
extern void _MGLU_reg(void);
extern void _naf_reg(void);
extern void _nap_reg(void);
extern void _NMDA_reg(void);

void modl_reg(){
  if (!nrn_nobanner_) if (nrnmpi_myid < 1) {
    fprintf(stderr, "Additional mechanisms from files\n");

    fprintf(stderr," AMPA.mod");
    fprintf(stderr," bkkca.mod");
    fprintf(stderr," cad.mod");
    fprintf(stderr," cadyn.mod");
    fprintf(stderr," caL13.mod");
    fprintf(stderr," caL13PS.mod");
    fprintf(stderr," caldyn.mod");
    fprintf(stderr," caL.mod");
    fprintf(stderr," caLPS.mod");
    fprintf(stderr," can.mod");
    fprintf(stderr," canPS.mod");
    fprintf(stderr," caq.mod");
    fprintf(stderr," caqPS.mod");
    fprintf(stderr," car.mod");
    fprintf(stderr," carPS.mod");
    fprintf(stderr," cat.mod");
    fprintf(stderr," damsg.mod");
    fprintf(stderr," ER.mod");
    fprintf(stderr," exp2synNMDA.mod");
    fprintf(stderr," GABA.mod");
    fprintf(stderr," kaf.mod");
    fprintf(stderr," kas.mod");
    fprintf(stderr," kir.mod");
    fprintf(stderr," krp.mod");
    fprintf(stderr," MGLU.mod");
    fprintf(stderr," naf.mod");
    fprintf(stderr," nap.mod");
    fprintf(stderr," NMDA.mod");
    fprintf(stderr, "\n");
  }
  _AMPA_reg();
  _bkkca_reg();
  _cad_reg();
  _cadyn_reg();
  _caL13_reg();
  _caL13PS_reg();
  _caldyn_reg();
  _caL_reg();
  _caLPS_reg();
  _can_reg();
  _canPS_reg();
  _caq_reg();
  _caqPS_reg();
  _car_reg();
  _carPS_reg();
  _cat_reg();
  _damsg_reg();
  _ER_reg();
  _exp2synNMDA_reg();
  _GABA_reg();
  _kaf_reg();
  _kas_reg();
  _kir_reg();
  _krp_reg();
  _MGLU_reg();
  _naf_reg();
  _nap_reg();
  _NMDA_reg();
}
