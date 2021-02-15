---gapfilled.MDS.usmms.R  Gapfilling variables with using MDS method.

---ETpm.calculation.daily.usmms.R  calculation for daily ET at US-MMS from 1999 to 2019.

---Variables in outputs.daily.usmms.csv
Tem   (deg C): Air temperature
LE    (W m-2): Latent heat flux
ustar (m s-1): Friction velocity
Hs    (W m-2): Sensible heat flux
U     (m s-1): Wind speed
Rn    (W m-2): Net radiation
VPD   (KPa): Vapor Pressure Deficit
SWC   (%): Soil water content (volumetric), range 0-100
ETlv_mmday (mm day-1) ET converted directly from LE
ga    (m s-1) aerodynamic conductance
Gs    (m s-1) surface conductance
ETpm_mmday (mm day-1) ET calculated with Penman-monteith equation
