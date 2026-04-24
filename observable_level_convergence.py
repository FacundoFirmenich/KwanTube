# observable_level_convergence.py
import pickle, glob, numpy as np

def load(nc, nk):
    with open(f"pade_ckpt_NC{nc}_Nk{nk}.pkl","rb") as f: 
        return pickle.load(f)

def populations(d):
    # P_i(t) para cada sitio i
    # rho_t contiene objetos Qobj de QuTiP
    return np.array([[s.full()[i,i].real for s in d["rho_t"]]
                     for i in range(d["rho_t"][0].shape[0])])

def max_coh(d):
    # max_{i<j} max_t |rho_ij|
    N = d["rho_t"][0].shape[0]
    out = np.zeros((N,N))
    for s in d["rho_t"]:
        M = np.abs(s.full())
        out = np.maximum(out, M)
    return out

print("="*80)
print(f"{'Comparativa':<15} | {'dPop':<10} | {'dCoh':<10} | {'dFrob':<10}")
print("-"*80)

for (nc1,nc2) in [(3,4),(4,5)]:
    for nk in (1,2,3,4):
        try:
            A, B = load(nc1,nk), load(nc2,nk)
        except FileNotFoundError:
            continue
        
        PA, PB = populations(A), populations(B)
        # Diferencia maxima en poblaciones a traves de todo el tiempo y todos los sitios
        dP   = np.max(np.abs(PA - PB))
        
        # Diferencia maxima en la envolvente de coherencias
        dC   = np.max(np.abs(max_coh(A) - max_coh(B)))
        
        # Diferencia maxima en norma de Frobenius (la metrica original)
        dF   = max(np.linalg.norm((a-b).full(),'fro')
                   for a,b in zip(A["rho_t"],B["rho_t"]))
        
        label = f"NC {nc1}->{nc2} Nk={nk}"
        print(f"{label:<15} | {dP:<10.2e} | {dC:<10.2e} | {dF:<10.2e}")
