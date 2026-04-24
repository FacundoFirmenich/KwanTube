# richardson_fragment.py
import pickle, numpy as np

def load(nc, nk=1):
    with open(f"pade_refined_ckpt_NC{nc}_Nk{nk}.pkl","rb") as f: 
        return pickle.load(f)

print("="*80)
print("  EXTRAPOLACION DE RICHARDSON (Fragmento 6DPU, Nk=1)")
print("="*80)

# Cargamos los ultimos 2 puntos para el calculo
NC_A, NC_B = 7, 8
A = load(NC_A)
B = load(NC_B)

# d_N = rho(8) - rho(7)
# r = d_8 / d_7  (medido de la tabla previa)
d7 = 5.45e-3
d8 = 2.09e-3
r = d8 / d7

print(f"Puntos usados: NC={NC_A} y NC={NC_B}")
print(f"Salto d_8: {d8:.2e}")
print(f"Ratio r (medido): {r:.3f}")

# Estimador de error total en NC=8: epsilon = d8 / (1-r)
epsilon_8 = d8 / (1 - r)
print(f"Error de truncacion estimado en NC=8: {epsilon_8:.2e} ({epsilon_8*100:.2f}%)")

# Estimador de error total en NC=7: epsilon = d7 / (1-r_prev)
# r_prev = d7 / d6 = 5.45e-3 / 1.29e-2 = 0.422
r_prev = 0.422
epsilon_7 = d7 / (1 - r_prev)
print(f"Error de truncacion estimado en NC=7: {epsilon_7:.2e} ({epsilon_7*100:.2f}%)")

# Proyeccion NC=9
d9 = d8 * r
epsilon_9 = d9 / (1 - r)
print(f"Error de truncacion proyectado en NC=9: {epsilon_9:.2e} ({epsilon_9*100:.2f}%)")

print("\n" + "-"*40)
print("  ESTADO ASINTOTICO (Richardson)")
print("-"*40)
# rho_inf = rho(8) + d8 * r / (1-r)
# Usaremos las poblaciones del sitio 0 como ejemplo
pop8 = np.array([s.full()[0,0].real for s in B["rho_t"]])
diff_pop = np.array([b.full()[0,0].real - a.full()[0,0].real for a,b in zip(A["rho_t"], B["rho_t"])])

pop_inf = pop8 + diff_pop * r / (1-r)

print(f"Poblacion sitio 0 final (t=500fs):")
print(f"  NC=7:      {A['rho_t'][-1].full()[0,0].real:.6f}")
print(f"  NC=8:      {B['rho_t'][-1].full()[0,0].real:.6f}")
print(f"  Asintota:  {pop_inf[-1]:.6f}")
print(f"  Correccion total desde NC=8: {pop_inf[-1] - pop8[-1]:.2e}")
