# 
# Example usage for the FreeTTES model
import FreeTTES_model as model
import numpy as np
from scipy.interpolate import interp1d

import time

start_time = time.time()
# original profile

start_profil = {
    2.0  : 27.63,
    6.0  : 28.39,
    10.0 : 28.39,
    14.0 : 28.39,
    18.0 : 28.39,
    22.0 : 28.39,
    26.0 : 28.39,
    30.0 : 28.42,
    34.0 : 31.07,
    38.0 : 44.13
}
def GJ_to_MWh(GJ):
    return GJ / 3.6
def m3h_to_kgs(m3h, t):
    return model.__Modell_Stoffwerte("rho", t) * m3h / 3600

m_charge = m3h_to_kgs(14,90)
m_discharge = m3h_to_kgs(14,30)
if __name__ == "__main__":
    t = 0
    while t <= 8760:
        if 0 <= t < 2920: # charge
            result = model.main(
                t=t,
                dt=3600,
                m_VL=m_charge,
                m_RL=-m_charge,
                T_Zustrom=90,
                T_amb=10.0,
                eingabe_volumen=False,
                zustand_uebernehmen=True,
                zustand=start_profil.copy()
            )
            # results are stored in result          
        elif 2920 <= t < 3650: # idle
            result = model.main(
                t=t,
                dt=3600,
                m_VL=0,
                m_RL=0,
                T_Zustrom=90,
                T_amb=10.0,
                eingabe_volumen=False,
                zustand_uebernehmen=False,
                zustand={}
            )
        elif 3650 <= t <= 6570: # discharge
            result = model.main(
                t=t,
                dt=3600,
                m_VL=-m_discharge,
                m_RL=m_discharge,
                T_Zustrom=30,
                T_amb=10.0,
                eingabe_volumen=False,
                zustand_uebernehmen=False,
                zustand={}
            )
        else: # idle
            result = model.main(
                t=t,
                dt=3600,
                m_VL=0,
                m_RL=0,
                T_Zustrom=30,
                T_amb=10.0,
                eingabe_volumen=False,
                zustand_uebernehmen=False,
                zustand={}
            )
        Speicherzustand = result["speicherzustand"]
        print(f"t={t}, T_Austritt={result['T_Austritt']} E_nutz={GJ_to_MWh(result['E_nutz'])}")
        z_grid = np.array(sorted(Speicherzustand.keys()), dtype=float)
        T_grid = np.array([Speicherzustand[z][0] for z in z_grid], dtype=float)

        # interpolator
        f = interp1d(z_grid, T_grid, kind="linear", bounds_error=True)

        # interpolate at start_profil heights
        start_profil_interpolated = {
            float(z): float(f(z))
            for z in start_profil.keys()
        }
        print(start_profil_interpolated)
        t += 1
    

print("Simulation time: %.2f s" % (time.time() - start_time))