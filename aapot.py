import numpy as np
import pandas as pd
import scipy.sparse as ss
import scipy.sparse.linalg as ssl
from scipy.special import spherical_jn, spherical_yn, gamma
from mpmath import coulombf, coulombg
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d


ali_bodmer_params = {'d2': {'mu_a':0.475, 'V_a':130, 'mu_r':0.7, 'V_r':320},
                     'e2': {'mu_a':0.5, 'V_a':150, 'mu_r':0.8, 'V_r':640}}   

def ali_bodmer_e2(r):
    return ali_bodmer(r, **ali_bodmer_params['e2'])

def ali_bodmer_d2(r):
    return ali_bodmer(r, **ali_bodmer_params['d2'])

def fg(r, A, sig, r0): #gaussian function
    return m / hbar**2 * A * np.exp(-(r - r0)**2 / (2 * sig**2))

def ali_bodmer_oo(r):
    f1 = fg(r, A=-3.0, sig=0.1, r0=2.6)
    f2 = fg(r, A=2.0, sig=0.1, r0=4.4)
    return ali_bodmer_d2(r) + f1 + f2

def ali_bodmer_av(r):
    return 0.5 * (ali_bodmer_d2(r) + ali_bodmer_e2(r))

def ali_bodmer(r, mu_a, V_a, mu_r, V_r): #r=separation in fm
    V = V_r * np.exp(-mu_r**2 * r**2) - V_a * np.exp(-mu_a**2 * r**2)
    return m / hbar**2 * V 

def centrifugal_potential(r): #r=separation in fm
    return l_orb * (l_orb + 1) / (2. * r**2)

def coulomb_potential(r): #r=separation in fm
    return Z1 * Z2 * alpha * m * c / (hbar * r)
        
def coulomb_phase_shift(E): #E=energy in MeV
    eta = alpha * Z1 * Z2 * np.sqrt(m * c**2 / (2 * E))  #sommerfeld parameter
    return np.angle(gamma(l_orb + 1 + 1j * eta), deg=True)
    
def interaction_hamiltonian(r):
    return np.ones(len(r))
    
def bound_state_wf(r):
    sigma = sigma_bound
    return np.exp(-r**2 / (2 * sigma**2))

def three_body_potential(r):
    r0 = 3.5 + 0.5 * r
    r12 = r
    r13 = np.sqrt(r0**2 + (r/2)**2)
    r23 = r13
    rho2 = 4 / 3. * (r12**2 + r13**2 + r23**2)
    S = -20 #92 #MeV
    b = 6 #fm
    V_3b = m / hbar**2 * S * np.exp(-rho2 / b**2)
    V_c = 2 * coulomb_potential(r13)
    V_2b = 2 * ali_bodmer_d2(r13)
    return V_3b + V_c + V_2b
    
def ali_bodmer_d2_3b(r):
    return ali_bodmer_d2(r) + three_body_potential(r) 




# ---- USER INPUT ----- #
l_orb = 2
alpha = 1./137.
Z1 = 2
Z2 = 2
c = 3.00e23 #fm/s
m = 0.5 * 3.727e3 / c**2  # MeV/c^2 (reduced mass)
hbar = 6.58e-22 #MeV s^-1
L = 80. #fm
n_eig = 40
n_bins = [600] #6000
potentials = {'d2': ali_bodmer_d2, 'e2': ali_bodmer_e2}
E_min = 0.1
E_max = 12.0
data_file = 'exp_data.csv'
show_comp_points = False
return_eigenvectors = True
sigma_bound = 3.0 #fm
# --------------------- #



# --- computation --- #

results = {'energy':[], 'phase_shift':[], 'label':[], 'rate':[]}

for name, nuclear_potential in potentials.items():
    print(name)
    
    for N in n_bins:
        print(f'N={N}')

        dr = L / (N + 1.)
        r = dr * np.arange(1, N+1)

        # potential operator
        V = ss.spdiags(nuclear_potential(r), [0], N, N) 

        # kinetic energy operator
        o = np.ones(N, dtype=complex)
        data = np.stack([o, -2*o, o], axis=0)
        diags = np.array([-1, 0, 1])  
        T = -0.5 / np.power(dr,2) * ss.spdiags(data, diags, N, N) 
        
        # centrifugal operator
        V_l = ss.spdiags(centrifugal_potential(r), [0], N, N) 

        # coulomb operator
        V_c = ss.spdiags(coulomb_potential(r), [0], N, N) 

        # Hamiltonian
        H = T + V + V_l + V_c

        # https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.linalg.eigsh.html#scipy.sparse.linalg.eigsh

        # convert to different sparse format needed for further calculation
        H = H.tocsc()

        # solve eigenvalue problem
        if return_eigenvectors:
            eig_val, eig_vec = ssl.eigsh(A=H, k=n_eig, which='SM', return_eigenvectors=True)
        else:
            eig_val = ssl.eigsh(A=H, k=n_eig, which='SM', return_eigenvectors=False)
            eig_vec = None

        # convert to MeV
        E = np.power(hbar, 2) / m * eig_val
        E = np.sort(E)
            
        # select scattering states
        idx = np.nonzero(E >= 0)[0] 

        # phase shift
        k = np.sqrt(2. * m * E[idx]) / hbar
        eta = alpha * Z1 * Z2 * np.sqrt(m * c**2 / (2 * E[idx]))  #sommerfeld parameter
        Fl, Gl = [], []
        for i in idx:
            k_i = np.sqrt(2. * m * E[i]) / hbar
            eta_i = alpha * Z1 * Z2 * np.sqrt(m * c**2 / (2 * E[i]))  #sommerfeld parameter
            f = float(coulombf(l=l_orb, eta=eta_i, z=k_i*L))
            g = float(coulombg(l=l_orb, eta=eta_i, z=k_i*L))
            Fl.append(f)
            Gl.append(g)
        
        kcotd = -k * np.array(Gl) / np.array(Fl)   
        delta = np.arctan(k / kcotd) * 180./np.pi
        
        # ensure phase shift is monotonously increasing    
        delta_new = [delta[0]]
        for d in delta[1:]:
            dy = d - delta_new[-1]
            if np.abs(dy) > 90: d -= np.round(dy / 90.) * 90.
            delta_new.append(d)
        
        results['energy'].append(E[idx])
        results['phase_shift'].append(delta_new)
        results['label'].append(f'{name} (N={N})')
        
        if return_eigenvectors:
            # eigen vector norm
            norm = np.real(np.sum(np.conjugate(eig_vec[:,idx]) * eig_vec[:,idx], axis=0) * dr)
            
            # Fermi's golden rule
            H_int = interaction_hamiltonian(r)[:,np.newaxis]
            wf_final = bound_state_wf(r)[:,np.newaxis]      
            rate = np.abs(np.sum(np.conjugate(wf_final) * H_int * eig_vec[:,idx], axis=0) * dr) ** 2 / norm

            results['rate'].append(rate)


# --- exp phase shift --- #
df = pd.read_csv(data_file)
x_exp = 0.5 * df['E_lab (MeV)'].values
y_exp = df['phase shift (deg)']
e_exp = df['uncertainty (deg)']


# compute deviations
results['chi2dof'] = list()
for x,y in zip(results['energy'], results['phase_shift']):
    f = interp1d(x, y, kind='cubic', bounds_error=False, fill_value=0)    
    chi2dof = np.sum((f(x_exp) - y_exp)**2 / e_exp**2) / len(x_exp)
    results['chi2dof'].append(chi2dof)

    
# --- plotting --- #

lstyles = ['-','--',':','.-']

nrows = 3 if return_eigenvectors else 2

fig, ax = plt.subplots(nrows=nrows, figsize=(6,5+2*(nrows-1)), sharex=True)

fmt = '-o' if show_comp_points else '-'

# plot phase shift
for i in range(len(results['energy'])):
    x = results['energy'][i]
    y = results['phase_shift'][i]
    chi2dof = results['chi2dof'][i]
    l = results['label'][i] + r' [$\chi^2/N={0:.2f}$]'.format(chi2dof)
    ls = lstyles[i%len(lstyles)]
    ax[1].plot(x, y, fmt, label=l, linestyle=ls)

# superimpose exp data
ax[1].errorbar(x_exp, y_exp, yerr=e_exp, fmt='o', label='exp.', color='black', capsize=3, markersize=3)

ax[1].set_ylabel('phase shift (degree)')    
if nrows==2: 
    ax[1].set_xlabel('C.M. energy (MeV)')    


# plot potential
x = np.linspace(E_min, E_max, 100)

i = 0
for name, nuclear_potential in potentials.items():
    y = nuclear_potential(x) + coulomb_potential(x) + centrifugal_potential(x)
    y *= hbar**2 / m
    ls = lstyles[i%len(lstyles)]
    ax[0].plot(x, y, label=name, linestyle=ls)
    ax[0].set_ylabel('Potential (MeV)')
    ax[0].set_xlabel('r (fm)')    
    ax[0].set_xlim(0, E_max)    
    ax[0].set_ylim(-10.0, 12.0)    
    i += 1


ax[0].legend()
ax[1].legend()

ax[0].grid(linestyle='dotted')
ax[1].grid(linestyle='dotted')


# plot capture rate
if return_eigenvectors:

    i = 0
    for x,y,l in zip(results['energy'], results['rate'], results['label']):
        f = interp1d(x, y, kind='cubic')
        _x = np.linspace(np.min(x), np.max(x), 1000)
        _y = f(_x)
        i0 = np.argmax(_y)        
        x0 = _x[i0] #peak position
        y0 = _y[i0]
        i1 = np.sum(_y[:i0] <= 0.5*y0) - 1
        i2 = len(_x) - np.sum(_y[i0:] <= 0.5*y0) - 1
        fwhm = _x[i2] - _x[i1]
        l += r' [E={0:.2f}; $\Gamma={1:.2f}$]'.format(x0, fwhm)
        ls = lstyles[i%len(lstyles)]
        ax[2].plot(_x, _y, fmt, label=l, linestyle=ls)
        i += 1

    ax[2].set_ylabel('Capture rate (a.u.)')    
    ax[2].set_xlabel('C.M. energy (MeV)')    
    ax[2].grid(linestyle='dotted')
    ax[2].legend()


plt.show()





