import numpy as na
from dengo.write_cvode_solver import create_initial_conditions
from dengo.chemistry_constants import tiny, kboltz, mh

NCELLS = 1024
rho = 1e17
X = na.logspace(-9,-3, NCELLS)
all = na.ones(NCELLS)

values = dict(HI    = None, # We apply conservation here
             HII   = X * all,
             HeI   = tiny * all,
             HeII  = tiny * all,
             HeIII = tiny * all,
             HM    = tiny * all,
             de    = X * all,
             H2I   = 0.99 * all,
             H2II  = tiny * all)

values["HI"] = (1.0 - values["H2I"] - values["HII"])

for v in values: values[v] *= rho

Temperature = na.logspace(2,na.log10(7000), NCELLS)
number_density = (values["HI"] + values["HII"] + values["de"] + values["HM"]) + \
                 (values["HeI"] + values["HeII"] + values["HeIII"])*0.25 + \
                 (values["H2I"] + values["H2II"]) * 0.5

# This is all very approximate ... and flat-out wrong at high H2 fraction.

values['ge'] = ( (Temperature * rho * kboltz)
               / (number_density * mh * (5.0/3.0 - 1))) # gamma ~ 5/3
values['rho'] = rho * all

tdyn = na.sqrt(3.0*na.pi / (6.67e-8*1.67e-24*rho))
create_initial_conditions(values, "primordial", tdyn*10.0)
