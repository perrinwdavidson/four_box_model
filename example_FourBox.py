# =========================================================
# example_FourBox
# ---------------------------------------------------------
# purpose :: to run a toggweiler and sarmiento (185)-type
#            4 box coupled-ocean atmosphere model.
# author :: perrin w. davidson
# date :: 5.12.22
# =========================================================
# configure -----------------------------------------------
# import libraries ::
from FourBox import FourBox  # for model class
import numpy as np  # for array functionality

# set conversions ::
y2s = 60 ** 2 * 24 * 365
c2p = 106

# run model -----------------------------------------------
# set model ::
model = FourBox()

# set initial values ::
j = [1.5 / y2s / c2p, 0.15 / y2s / c2p]  # phosphate flux [mol P m-2 s-1]
f_hd = 50E6  # high latitude exchange [m+3 s-1]

# solve phosphate concentrations ::
c_po4 = model.solve_po4(j, f_hd)

# solve alkalinity concentrations ::
c_alk = model.solve_alk(j, f_hd)

# solve co2 concentrations ::
c_co2 = model.solve_co2(j, f_hd)

# make ranges of data ::
j_range = np.linspace(0, 2.5 / y2s / c2p, 100)
f_hd_range = np.linspace(0, 100E6, 100)

# run experiment ::
exp_data = model.contour_data(j_range, f_hd_range)

# end example script --------------------------------------
# =========================================================
