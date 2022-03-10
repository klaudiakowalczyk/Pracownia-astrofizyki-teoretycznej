import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as const
import astropy.constants as aconst

def ruch(nsteps,r,v):
    dim = len(r)
    ri = np.sqrt(r[0][0]**2+r[1][0]**2)
    for k in range(dim):
        a[k][0] = -G*M*r[k][0]/(ri**3)
        v[k][1] = v[k][0]+dt/2.*a[k][0]
        r[k][1] = r[k][0]+dt*v[k][1]
        t[1] = t[0]+dt

    for i in range(1,nsteps):
        ri = np.sqrt(r[0][i]**2+r[1][i]**2)

        for k in range(dim):
            a[k][i] = -G*M*r[k][i]/(ri**3)
            v[k][i+1] = v[k][i]+dt*a[k][i]
            r[k][i+1] = r[k][i]+dt*v[k][i+1]

        t[i+1] = t[i]+dt
    return t, r, v, a


params = {'font.size': 18,
          'figure.figsize': (12.,8.),
          'axes.labelsize': 18,
          'xtick.labelsize': 18,
          'ytick.labelsize': 18,
          'text.latex.preamble': r"\usepackage{amsmath}"
         }
plt.rcParams.update(params)

# ----------------- przypadek z książki -----------------

G = 1.
M = 1. # żeby G*M = 1

nsteps = 23
dt = 0.1

dim = 2
t = np.zeros(nsteps+1)
r = np.zeros((dim,nsteps+1))
v = np.zeros((dim,nsteps+1))
a = np.zeros((dim,nsteps+1))

t[0] = 0.0
r[0][0] = 0.5
r[1][0] = 0.
v[0][0] = 0.
v[1][0] = 1.63

t_out, r_out, v_out, a_out = ruch(nsteps,r,v)

f1 = open("tabela.dat","w")

for i in range(nsteps):
    r = np.sqrt(r_out[0][i]**2+r_out[1][i]**2)
    f1.write(f"{t_out[i]:.3f} & {r_out[0][i]:.3f} & & {a_out[0][i]:.3f} & {r_out[1][i]:.3f} & & {a_out[1][i]:.3f} & {r:.3f} & {1/(r**3):.3f}\\\\\n & & {v_out[0][i]:.3f} & & & {v_out[1][i]:.3f} & & & \\\\\n")

f1.close()

# ----------------- wykresy -----------------

plt.figure(1,tight_layout=True)

plt.scatter(r_out[0], r_out[1], s=20, c='k')
plt.scatter(r_out[0][::5], r_out[1][::5], s=20, c='tab:green')
plt.scatter(0, 0, s=300, c='yellow')
plt.text(0, 0.03, 'Słońce', ha='center', transform=plt.gca().transData)
plt.gca().axis('equal')
plt.xlabel(r"$x$")
plt.ylabel(r"$y$")

plt.savefig('ruch.pdf')


# ----------------- przypadek z książki dla większej liczby kroków -----------------

nsteps = 10000
dt = 0.01

dim = 2
t = np.zeros(nsteps+1)
r = np.zeros((dim,nsteps+1))
v = np.zeros((dim,nsteps+1))
a = np.zeros((dim,nsteps+1))

t[0] = 0.0
r[0][0] = 0.5
r[1][0] = 0.
v[0][0] = 0.
v[1][0] = 1.63

t_out, r_out, v_out, a_out = ruch(nsteps,r,v)

# ----------------- wykresy -----------------

plt.figure(2,figsize=(12,11),tight_layout=True)

plt.plot(r_out[0], r_out[1], lw=0.2, c='k')
plt.scatter(0, 0, s=300, c='yellow')
plt.gca().axis('equal')
plt.xlabel(r"$x$")
plt.ylabel(r"$y$")

plt.savefig('ruch2.pdf')


# ----------------- inne warunki początkowe -----------------

# 6.67430 x 10-11 m3 kg-1 s-2
pc = aconst.pc.value # m
m2pc = 1/pc
M_sun = aconst.M_sun.value # kg
kg2Ms = 1/M_sun
Myr = const.Julian_year*10**6 # s
s2Myr = 1/Myr
G = 6.67430*10**(-11)*(m2pc)**3*(kg2Ms)**(-1)*(s2Myr)**(-2) # pc3 M_s-1 Myr-2
M = 10**12 # M_sun

time = 3000 # Myr
dt = 1 # Myr
nsteps = int(time/dt)

dim = 2
r = np.zeros((dim,nsteps+1))
v = np.zeros((dim,nsteps+1))
a = np.zeros((dim,nsteps+1))
vi = np.zeros(dim)

r[0][0] = 8000. # pc
r[1][0] = 0.
r0 = np.sqrt(r[0][0]**2+r[1][0]**2)
v[0][0] = 0.
v[1][0] = np.sqrt(G*M/r0)

t_out, r_out, v_out, a_out = ruch(nsteps,r,v)

# ----------------- wykresy -----------------

plt.figure(3,figsize=(12,11),tight_layout=True)

plt.plot(r_out[0], r_out[1], lw=0.2, c='k')
plt.scatter(0, 0, s=300, c='k')
plt.gca().axis('equal')
plt.xlabel(r"$x$")
plt.ylabel(r"$y$")

plt.savefig('ruch3.pdf')


# -------------------------------------------

G = 1.
M = 1. # żeby G*M = 1

dt_list = [0.01, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5]
dim = 2

x_text = [0.095, 0.575, 0.095, 0.575, 0.095, 0.575, 0.095, 0.575]
y_text = [0.962, 0.962, 0.718, 0.718, 0.475, 0.475, 0.233, 0.233]

fig, ax = plt.subplots(nrows=4, ncols=2, sharex=True, figsize=(18,24))
plt.subplots_adjust(left=0.08, right=0.96, bottom=0.05, top=0.98)

k = 0
for i in range(4):
    for j in range(2):
        if k == 0:
            nsteps = 4000
        else:
            nsteps = 2000
        dt = dt_list[k]

        t = np.zeros(nsteps+1)
        r = np.zeros((dim,nsteps+1))
        v = np.zeros((dim,nsteps+1))
        a = np.zeros((dim,nsteps+1))

        t[0] = 0.0
        r[0][0] = 1.0
        r[1][0] = 0.0
        v[0][0] = 0.0
        v[1][0] = 1.3

        t_out, r_out, v_out, a_out = ruch(nsteps,r,v)

        ax[i,j].scatter(r_out[0], r_out[1], s=10, c='k')
        ax[i,j].scatter(0, 0, s=300, c='yellow')

        plt.figtext(x_text[k], y_text[k], f'$dt=${dt:.2f}')

        ax[i,j].axis('equal')
        ax[i,j].set_xlabel(r"$x$")
        ax[i,j].set_ylabel(r"$y$")

        k += 1

plt.savefig('ruch4.pdf')



# -------------------------------------------

G = 1.
M = 1. # żeby G*M = 1

dt_list = [0.01, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5]
dim = 2
x_start = [1.0, 2.0, 1.0, 1.5, 1.8, 1.1, 1.2, 1.0]
y_start = [0.0, 0.5, 1.0, 0.2, 0.8, -0.1, 0.4, 0.0]
vx_start = [-0.5, 0.0, 0.0, -0.5, -0.5, 0.0, 0.0, 0.5]
vy_start = [1.2, 0.8, 1.1, 0.8, 0.5, 1.3, 1.0, 1.0]

x_text = [0.095, 0.575, 0.095, 0.575, 0.095, 0.575, 0.095, 0.575]
y_text = [0.962, 0.962, 0.718, 0.718, 0.475, 0.475, 0.233, 0.233]

fig, ax = plt.subplots(nrows=4, ncols=2, figsize=(18,24))
plt.subplots_adjust(left=0.08, right=0.96, bottom=0.05, top=0.98)

k = 0
for i in range(4):
    for j in range(2):
        if k == 0:
            nsteps = 4000
        else:
            nsteps = 2000
        dt = 0.1

        t = np.zeros(nsteps+1)
        r = np.zeros((dim,nsteps+1))
        v = np.zeros((dim,nsteps+1))
        a = np.zeros((dim,nsteps+1))

        t[0] = 0.0
        r[0][0] = x_start[k]
        r[1][0] = y_start[k]
        v[0][0] = vx_start[k]
        v[1][0] = vy_start[k]

        t_out, r_out, v_out, a_out = ruch(nsteps,r,v)

        ax[i,j].scatter(r_out[0], r_out[1], s=10, c='k')
        ax[i,j].scatter(0, 0, s=300, c='yellow')

        plt.figtext(x_text[k], y_text[k]-0.03, f'$x_0=${r[0][0]:.2f}\n$y_0=${r[1][0]:.2f}\n$v_x=${v[0][0]:.2f}\n$v_y=${v[1][0]:.2f}')

        ax[i,j].axis('equal')
        ax[i,j].set_xlabel(r"$x$")
        ax[i,j].set_ylabel(r"$y$")

        k += 1

plt.savefig('ruch5.pdf')
