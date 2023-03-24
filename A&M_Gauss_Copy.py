import numpy as np
import sympy as sy  # Using sympy for arbitrary precision floating-point arithmetic in SCFT loop
import matplotlib.pyplot as plt
import pyvista as pyv  # for 3D isosurface plotting
from overlap import findoverlap  # Wrapped Fortran module containing function to calculate overlap matrix.
from laplace import findlaplace  # Wrapped Fortran module containing function to calculate laplace matrix.
from gammatensor import findgamma  # Wrapped Fortran module containing function to calculate gamma tensor.
from numpy.linalg import solve
from scipy.linalg import eigh  # Only eigh/eig from scipy has generalized eigenvalue problem solver.
from scipy.special import gamma
from time import time
from os.path import normpath
import matplotlib.style as mplstyle
mplstyle.use('fast')

np.set_printoptions(precision=16, edgeitems=10, linewidth=125, infstr='oo', sign=' ')

start = time()


# Real Spherical Harmonics with Radial Gaussians as Basis Functions.

def Gauss(u, c, l, m):

    rad = (2.0*(2.0*c)**(l+1.5)/gamma(l+1.5))**0.5*u[0]**l*sy.exp(-c*u[0]**2)
    if m > 0:
        return rad*sy.re(sy.Ynm(l, m, u[1], u[2]).expand(func=True))*2.0**0.5
    elif m == 0:
        return rad*sy.Ynm(l, 0, u[1], u[2]).expand(func=True)
    elif m < 0:
        return rad*sy.im(sy.Ynm(l, abs(m), u[1], u[2]).expand(func=True))*2.0**0.5*(-1.0)**m


# Function to calculate the electron densities and free energies self-consistently.

def GaussDensity(u, pvec, n, Ns, l, m, l_num, b, db, itermax, tol, realint, realint1, realint2, g0inv, mix, anders, andtol, angmom, atoms, rand, switch, ind_ang, ind_sh, step_mat):

    gtime = time()  # Start time of SCFT algorithm.
    N = np.sum(Ns)  # Electron number.
    Nsk = np.size(Ns)  # Electron pair number.
    Lsk = np.size(l)  # Number of angular momentum values.
    Msk = [np.size(m[i]) for i in range(Lsk)]  # Number of m values for each angular momentum value.
    ver_field = FA[FAval][-1]  # Data label for which self-interaction correction/exchange field is being used.
    vers = np.array(['New', ''])  # Data labels for angular and spherical field components.
    vers_bool = np.array([angmom, True])  # Boolean indicating whether angular components are being considered or not
    if N > 46:  # Show more digits of beta in the input/output file name when going to higher atomic number.
        beta_prec = np.around(b+db[0], decimals=3)
        beta_prec_sav = np.around(b, decimals=3)
    else:
        beta_prec = np.around(b+db[0], decimals=1)
        beta_prec_sav = np.around(b, decimals=1)

    ver = ''.join(vers[vers_bool])
    Nbase = np.sum(n)  # Total number of basis functions.
    wb_init = np.array([np.zeros(i) for i in n])  # Initial field values.
    wb_append_shape = np.array([np.zeros(n[-1]+db[1][-1])])  # Initial field values to add-on for the next pair.
    temp_index = n
    temp_andtol = andtol
    temp_it = 0
    cvec = np.array([np.logspace(pvec[i, 0], pvec[i, 1], num=l_num[i]) for i in range(Lsk)])  # Gaussian exponents for each angular momentum value.
    cvt = np.hstack(np.array([np.tile(cvec[i], Msk[i]) for i in range(Lsk)]))  # Gaussian exponents repeated according to the number of m values corresponding to each angular momentum value.
    Gaussfunc = [Gauss(u, cvt[i], lvt[i], mvt[i]) for i in range(temp_index[0])]  # List of symbolically defined basis functions.
    S = findoverlap(cvt, lvt, mvt, temp_index[0])  # Overlap matrix.
    L = findlaplace(cvt, lvt, mvt, temp_index[0])  # Laplace Matrix.
    Gamma = findgamma(cvt, lvt, mvt, temp_index[0])  # Gamma Tensor.
    w_c = 4.0*np.pi*N*solve(L, sy.lambdify(u, Gaussfunc)(0, 0, 0))  # Spectral components of the electron-nucleus field.
    w_cshell = np.tile(w_c, (Nsk, 1))  # Spectral components of the electron-nucleus field for each pair.

    if N > 1:
        if current_atom:
            a0val = atoms[N]
            a1val1 = n
            a1val2 = n
            a1val3 = n
            gval1 = np.around(g0inv+g0inv_step, decimals=2)
        else:
            a0val = atoms[N-1]
            a1val1 = n[0:Nsk-1]+db[1][0:Nsk-1]
            a1val2 = np.append(n, n[0])+np.append(db[1], db[1][0])
            a1val3 = n + db[1]
            gval1 = np.around(g0inv+g0inv_step, decimals=2)

        if np.any(Ns == 1) and np.all(N != np.array([24, 29, 41, 42, 44, 45, 46, 78, 79, 110])):
            data = np.load(normpath('C:/Users/Owner/PycharmProjects/Research/Publication Data/Spherical Averaging/{a0}{a1}beta{a2}{a3}{a4}g0inv{a5}.npz').format(a0=a0val, a1=a1val1, a2=beta_prec, a3=ver, a4=ver_field, a5=gval1))
        elif N == 46:
            data = np.load(normpath('C:/Users/Owner/PycharmProjects/Research/Publication Data/Spherical Averaging/{a0}{a1}beta{a2}{a3}{a4}g0inv{a5}.npz').format(a0=a0val, a1=a1val2, a2=beta_prec, a3=ver, a4=ver_field, a5=np.around(g0inv+g0inv_step, decimals=2)))
        else:
            data = np.load(normpath('C:/Users/Owner/PycharmProjects/Research/Publication Data/Spherical Averaging/{a0}{a1}beta{a2}{a3}{a4}g0inv{a5}.npz').format(a0=a0val, a1=a1val3, a2=beta_prec, a3=ver, a4=ver_field, a5=np.around(g0inv+g0inv_step, decimals=2)))

        wb = np.nan_to_num(data['wb'])  # Initial field guess using previous atom output.

        prev_rdensshell = np.nan_to_num(data['rdensshell'])  # Initial pair density guess using previous atom output.
        prev_rdenstot = np.sum(prev_rdensshell, axis=0)  # Initial density guess using previous atom output.
    else:
        wb = wb_init  # Initial field guess of all zeroes for hydrogen.
        prev_rdenstot = 0

    iter = 0
    dev1_path = 0
    dev2_path = 0
    dev1_pathd = 0
    dev2_pathd = 0
    wbshape = np.shape(wb)
    wold = np.zeros((anders+1, Nbase))  # Initialize field output histories.

    if wbshape[0] != Nsk:  # add a set of coefficients for a new shell if importing from a file containing elements from the previous column of the periodic table.
        if N != 46:
            wb = np.append(wb, wb_append_shape, axis=0)
        else:
            wb = wb[:Nsk]
        wbshape = np.shape(wb)
    if np.any(np.array(db[2]) > 0): # pad each field coefficient array with a number of zeros equal to the difference between the current number of basis functions and the number contained in the imported file.
        if Lsk > 1:
            wb_ind_m = np.append([slice(0, Msk[0], None)], [slice(Msk[i], Msk[i]+Msk[i+1], None) for i in range(Lsk-1)], axis=0)
            wb_ind = np.append([slice(0, l_num[0]-db[2][0], None)], np.hstack([[slice(l_num[i]-db[2][i]+j*(l_num[i+1]-db[2][i+1]), l_num[i]-db[2][i]+l_num[i+1]-db[2][i+1]+j*(l_num[i+1]-db[2][i+1]), None) for j in range(Msk[i+1])] for i in range(Lsk-1)]), axis=0)
            wb_ind_idx = np.append([[0, l_num[0]-db[2][0]]], np.hstack([[[l_num[i]-db[2][i]+j*(l_num[i+1]-db[2][i+1]), l_num[i]-db[2][i]+l_num[i+1]-db[2][i+1]+j*(l_num[i+1]-db[2][i+1])] for j in range(Msk[i+1])] for i in range(Lsk-1)]), axis=0)
            db_ind = np.hstack([[l[j] for i in range(np.size(wb_ind[wb_ind_m[j]]))] for j in range(Lsk)])
            wb_rand = [[switch[1]*((max(wb[k, wb_ind[i]])-min(wb[k, wb_ind[i]]))*rng.random((db[2][db_ind[i]], )) + min(wb[k, wb_ind[i]])) for i in range(np.size(wb_ind))] for k in range(wbshape[0])]
            wb_idx = [[rng.integers(low=wb_ind_idx[k, 0], high=wb_ind_idx[k, 1], size=db[2][db_ind[i]]) for i in range(np.size(wb_ind))] for k in range(wbshape[0])]
            wb = np.array([np.insert(wb[i], np.hstack(wb_idx[i]), np.hstack(wb_rand[i])) for i in range(wbshape[0])])
        else:
            wb_ind_m = np.array([slice(0, Msk[0], None)])
            wb_ind = np.array([slice(0, l_num[0]-db[2][0], None)])
            wb_ind_idx = np.array([[0, l_num[0]-db[2][0]]])
            db_ind = np.hstack([[l[j] for i in range(np.size(wb_ind[wb_ind_m[j]]))] for j in range(Lsk)])
            wb_rand = [[switch[1]*((max(wb[k, wb_ind[i]])-min(wb[k, wb_ind[i]]))*rng.random((db[2][db_ind[i]], )) + min(wb[k, wb_ind[i]])) for i in range(np.size(wb_ind))] for k in range(wbshape[0])]
            wb_idx = [[rng.integers(low=wb_ind_idx[k, 0], high=wb_ind_idx[k, 1], size=db[2][db_ind[i]]) for i in range(np.size(wb_ind))] for k in range(np.size(wb_ind_idx))]
            wb = np.array([np.insert(wb[i], np.hstack(wb_idx[i]), np.hstack(wb_rand[i])) for i in range(wbshape[0])])

    wtot = np.ravel(wb)
    dev = np.zeros((anders+1, Nbase))  # Initialize deviation function histories.
    devtot = np.full(anders+1, 10000.0)  # Initialize deviation total histories.
    realints = np.linspace(0, 1.5*rc, realint)  # Radial grid values.
    realints1 = np.linspace(0, np.pi, realint1)  # Theta grid values.
    realints2 = np.linspace(0, 2.0*np.pi, realint2)  # Phi grid values.
    ri, ri1, ri2 = np.meshgrid(realints, realints1, realints2, indexing='ij')  # For outer numerical integration.
    ri_1, ri1_1 = np.meshgrid(realints, realints1, indexing='ij')  # For second numerical integration.

    # Generate uniform stochastic field using 'rand' values for endpoints.
    fake_ind = [slice(0, Msk[0]*l_num[0], None)] + [slice(Msk[i]*l_num[i], Msk[i]*l_num[i]+Msk[i+1]*l_num[i+1], None) for i in range(Lsk-1)]
    rng_wold_temp = switch[0]*((rand[1]-rand[0])*rng.random((temp_index[0], )) + rand[0])
    for k in range(Lsk):  # Suppress desired angular components.
        rng_wold_temp[fake_ind[k]] = ind_ang[k]*rng_wold_temp[fake_ind[k]]

    rng_wold = np.ravel([ind_sh[j]*rng_wold_temp for j in range(Nsk)])  # Suppress desired pair components.
    wtot = wtot + rng_wold  # Add stochastic field to initial field guess to help find other solutions.

    iter_and = 0
    wold_temp = np.array([wtot[i*temp_index[i]:(i+1)*temp_index[i]] for i in range(Nsk)])
    rdensshell = np.array([np.zeros(temp_index[i]) for i in range(Nsk)])
    rdensshell_q = np.array([np.zeros((temp_index[i], temp_index[i])) for i in range(Nsk)])
    qshell = [0 for i in range(Nsk)]
    while (devtot[0] > tol) and (iter <= itermax):  # SCFT loop to calculate the electron density.

        for i in range(Nsk):  # Loop for each pair.
            A = 0.5*L-np.einsum('ijk,i -> jk', Gamma, wold_temp[i])  # Compute matrix on RHS of propagator differential equation.
            val, vec = eigh(A, S)  # Find eigenvalues and normalized eigenvectors.
            dval = np.array([sy.exp(b*val[j]) for j in range(temp_index[0])], dtype='object')  # Symbolic eigenvalues of the propagator matrix.
            q = np.sum(dval)  # Symbolic single-pair partition function.
            dvalq = np.array(dval/q).astype(np.float64)  # Compute using symbolic expressions to avoid overflow from large beta.
            if iter == 0:
                qp_path = np.einsum_path('ij,j,kj -> ik', vec, dvalq, vec, optimize='optimal')[0]
            qp = np.einsum('ij,j,kj -> ik', vec, dvalq, vec, optimize=qp_path)  # Propagator matrix divided by single-pair partition function.
            rdensshell_q[i] = Ns[i]*qp
            rdensshell[i] = solve(S, np.einsum('ijk,ij -> k', Gamma, Ns[i]*qp))  # Solve for electron density components.
            qshell[i] = q

        rdenstot = np.sum(rdensshell, axis=0)  # Total electron density.
        w_pshell = np.array([g0inv*(rdenstot-rdensshell[i]) for i in range(Nsk)])  # Pauli field per pair.
        w_eeshell = np.tile(-4.0*np.pi*solve(L, np.einsum('ij,i->j', S, rdenstot)), (Nsk, 1))  # Electron-electron field per pair.
        w_xcshell = np.array([-4.0*np.pi*FA[FAval][1][i]*solve(L, np.einsum('ij,i->j', S, FA[FAval][2]*rdenstot+FA[FAval][3]*rdensshell[i]-FA[FAval][0]*prev_rdenstot)) for i in range(Nsk)])  # Self-interaction correction/exchange field per pair.
        wshell = w_cshell + w_eeshell + w_xcshell + w_pshell  # Total field per pair.

        wnew = np.ravel(wshell)  # Combine pair field components into one array.
        wold = np.roll(wold, 1, axis=0)  # Store old output fields.
        wold[0] = wtot  # Update field histories.
        dev = np.roll(dev, 1, axis=0)  # Store old deviation functions.
        devtot = np.roll(devtot, 1, axis=0)  # Store old spectral convergence value.
        dev[0] = wnew-wold[0]  # New deviation functions.
        devshell = wshell - wold_temp
        if devtot[0] > tol2vals[A_nuc[0]-1]:  # Coarse convergence.
            andtol = temp_andtol
            temp_it = iter
            if iter == 0:
                dev1_path = np.einsum_path('i,j,ij -> ', devshell[0], devshell[0], S, optimize='optimal')[0]
                dev2_path = np.einsum_path('i,j,ij -> ', wshell[0], wshell[0], S, optimize='optimal')[0]
            devtot_upper = np.abs(np.array([np.einsum('i,j,ij -> ', devshell[h], devshell[h], S, optimize=dev1_path) for h in range(Nsk)]))
            devtot_lower = np.abs(np.array([np.einsum('i,j,ij -> ', wshell[h], wshell[h], S, optimize=dev2_path) for h in range(Nsk)]))
        else:  # Fine convergence, using the density as a weighting factor.
            andtol = tol  # Turn off Anderson acceleration.
            if iter == temp_it+1:
                dev1_pathd = np.einsum_path('i,j,k,ijk -> ', devshell[0], devshell[0], rdensshell[0], Gamma, optimize='optimal')[0]
                dev2_pathd = np.einsum_path('i,j,k,ijk -> ', wshell[0], wshell[0], rdensshell[0], Gamma, optimize='optimal')[0]
            devtot_upper = np.abs(np.array([np.einsum('i,j,k,ijk -> ', devshell[h], devshell[h], rdensshell[h], Gamma, optimize=dev1_pathd) for h in range(Nsk)]))
            devtot_lower = np.abs(np.array([np.einsum('i,j,k,ijk -> ', wshell[h], wshell[h], rdensshell[h], Gamma, optimize=dev2_pathd) for h in range(Nsk)]))
        devtot[0] = np.sqrt(np.sum(devtot_upper)/np.sum(devtot_lower))  # Spectral convergence criterion.
        perdevtot = abs((devtot[0]-devtot[1])/devtot[0])  # Percent deviation between iterations.

        if iter_and < 20:
            if (perdevtot < andtol and devtot[0] < devtol) and iter > anders:  # Decide whether to do an Anderson step.
                wmix = wold[0]
                avecmix = anderson(dev, anders)
                wmix = wmix+np.sum((wold[1:anders+1] - wold[0])*np.vstack(avecmix), axis=0)
                print('Anderson Step')
                iter_and += 1
            else:
                wmix = mix*wnew+(1.0-mix)*wtot  # Simple Picard mixing.
                iter_and = 0
        else:
            wmix = mix*wnew+(1.0-mix)*wtot  # Simple Picard mixing.
            iter_and -= 1

        print('Iteration: {a0} | Spectral Convergence: {a1}'.format(a0=iter, a1=devtot[0]))
        if iter == int(0.4*itermax):  # Decrease tolerance values if iterations have gone on for too long.
            tol = 10*tol
            tol2vals[A_nuc[0]-1] = 10*tol2vals[A_nuc[0]-1]

        # Update spectral field components.
        wtot = wmix + step_mat[iter]*rng.permutation(rng_wold)
        wold_temp = np.array([wtot[i*temp_index[i]:(i+1)*temp_index[i]] for i in range(Nsk)])

        iter += 1

    # Free energy contributions from each interaction per pair.
    F_ke_pair = -np.array([Ns[j]*sy.log(qshell[j])/b for j in range(Nsk)]).astype(np.float64)
    F_ee_pair = U_ee_pair-np.array([np.einsum('k,j,jk -> ', rdensshell[i], w_eeshell[i], S, optimize='greedy') for i in range(Nsk)])
    F_xc_pair = U_xc_pair-np.array([np.einsum('k,j,jk -> ', rdensshell[i], w_xcshell[i], S, optimize='greedy') for i in range(Nsk)])
    F_p_pair = U_p_pair-np.array([np.einsum('k,j,jk -> ', rdensshell[i], w_pshell[i], S, optimize='greedy') for i in range(Nsk)])

    woldfunc = np.sum(np.array([np.sum(wtot[temp_index[i]*i:(i+1)*temp_index[i]]*Gaussfunc) for i in range(Nsk)]))
    wshellfunc = [np.sum(wnew[temp_index[i]*i:(i+1)*temp_index[i]]*Gaussfunc) for i in range(Nsk)]

    # Potential energies for each interaction per pair
    U_c_pair = np.array([np.einsum('k,j,jk -> ', rdensshell[i], w_cshell[i], S, optimize='greedy') for i in range(Nsk)])
    U_p_pair = 0.5*np.array([np.einsum('k,j,jk -> ', rdensshell[i], w_pshell[i], S, optimize='greedy') for i in range(Nsk)])
    U_ee_pair = 0.5*np.array([np.einsum('k,j,jk -> ', rdensshell[i], w_eeshell[i], S, optimize='greedy') for i in range(Nsk)])
    U_xc_pair = 0.5*np.array([np.einsum('k,j,jk -> ', rdensshell[i], w_xcshell[i], S, optimize='greedy') for i in range(Nsk)])

    # Lambda functions for density and field functions.
    temps = [np.sum(rdensshell[i]*Gaussfunc) for i in range(Nsk)]
    temp = sum(temps)
    realn = sy.lambdify(u, temp)
    realns = sy.lambdify(u, temps)
    wnewfunc = sum(wshellfunc)
    realwnew = sy.lambdify(u, wnewfunc)
    realwold = sy.lambdify(u, woldfunc)
    realwshell = sy.lambdify(u, wshellfunc)

    enum = np.sum(rdenstot[0:l_num[0]]*(2.0*np.pi/cvec[0])**0.75)  # Total analytic electron number.
    enums = np.array([np.sum(rdensshell[i][0:l_num[0]]*(2.0*np.pi/cvec[0])**0.75) for i in range(Nsk)])  # Shell analytic electron number.

    print('Spectral Electron Number: {}'.format(enum))
    print('Spectral Electron Shell Numbers: {}'.format(repr(enums)))

    if angmom:
        realnval1 = realn(ri, ri1, ri2)
        realnvals1 = realns(ri, ri1, ri2)
        enum = np.trapz(np.trapz(np.trapz(np.sin(ri1)*ri**2*realnval1, ri2), ri1_1), realints)  # Numerically-integrated total electron number.
        enums = np.array([np.trapz(np.trapz(np.trapz(np.sin(ri1)*ri**2*realnvals1[i], ri2), ri1_1), realints) for i in range(Nsk)])  # Numerically-integrated electron number per pair.
        realconv = np.trapz(np.trapz(np.trapz(np.sin(ri1)*ri**2*np.abs(realwold(ri, ri1, ri2)-realwnew(ri, ri1, ri2))**2, ri2), ri1_1), realints)
        realconv_low = np.trapz(np.trapz(np.trapz(np.sin(ri1)*ri**2*(realwold(ri, ri1, ri2))**2, ri2), ri1_1), realints)
        position_entropies1 = np.array([np.trapz(np.trapz(np.trapz(np.sin(ri1)*ri**2*np.abs(realnvals1[i])*np.log(np.abs(realnvals1[i])), ri2), ri1_1), realints) for i in range(Nsk)])
    else:
        realnval = realn(realints, realints1, realints2)
        realnvals = realns(realints, realints1, realints2)
        enum = 4.0*np.pi*(np.trapz(realints**2*realnval, realints))  # Numerically-integrated total electron number.
        enums = 4.0*np.pi*np.array([(np.trapz(realints**2*realnvals[i], realints)) for i in range(Nsk)])  # Numerically-integrated electron number per pair.
        realwdiffval = abs(realwnew(realints, realints1, realints2)-realwold(realints, realints1, realints2))**2
        realconv_low = 4.0*np.pi*(np.trapz(realints**2*realwold(realints, realints1, realints2)**2, realints))
        realconv = 4.0*np.pi*(np.trapz(realints**2*realwdiffval, realints))
        position_entropies1 = 4.0*np.pi*np.array([np.trapz(realints**2*np.abs(realnvals[i])*np.log(np.abs(realnvals[i])), realints) for i in range(Nsk)])
    realconv = np.sqrt(realconv/realconv_low)  # Real space convergence value.

    prop_en_pair = -np.array([Ns[j]*sy.log(qshell[j]/Ns[j]) for j in range(Nsk)]).astype(np.float64)-position_entropies1  # Quantum kinetic entropy minus density-field integral per pair.
    elec_en = Ns*np.log(Ns)  # Electron number entropy per pair.
    kin_en_R = -(np.array(sum([Ns[j]*sy.log(qshell[j]/Ns[j]) for j in range(Nsk)])).astype(np.float64)+np.sum(position_entropies1))/b+2.0*(F_xc+F_p+F_ee)-U_c  # Quantum kinetic energy.
    kin_en_R_pair = prop_en_pair/b-2.0*(U_xc_pair+U_p_pair+U_ee_pair)-U_c_pair  # Quantum kinetic energy per pair.
    kin_en = -0.5*np.einsum('ij,ji -> ', np.sum(rdensshell_q, axis=0), L)  # Total kinetic energy.

    print('Kinetic Energy: {}'.format(kin_en))
    print('Configurational Entropy: {}'.format(b*kin_en_R))
    print('Configurational Entropy per Pair: {}'.format(repr(b*kin_en_R_pair)))
    print('Translational Entropy per Pair: {}'.format(repr(position_entropies1-elec_en)))
    print('Translational Entropy: {}'.format(np.sum(position_entropies1-elec_en)))

    deci = 5
    F_pair = -(U_ee_pair+U_xc_pair+U_p_pair)+F_ke_pair+U_nuc_pair
    U_pair = U_ee_pair+U_xc_pair+U_p_pair+U_c_pair+U_nuc_pair
    print('Free Energy Terms:\n Kinetic: {a0} \n Electron-Electron: {a1} \n Self-Interaction: {a2} \n Pauli: {a3} \n Total: {a5}'.format(a0=repr(F_ke_pair), a1=repr(F_ee_pair), a2=repr(F_xc_pair), a3=repr(F_p_pair), a5=repr(F_pair)))
    print('Potential Energy Terms:\n Electron-Nucleus: {a0} \n Electron-Electron: {a1} \n Self-Interaction: {a2} \n Pauli: {a3} \n Total: {a5}'.format(a0=repr(U_c_pair), a1=repr(U_ee_pair), a2=repr(U_xc_pair), a3=repr(U_p_pair), a5=repr(U_pair)))

    entropy = b*(kin_en+U_c-F_ee-F_xc-F_p-free)  # Total entropy.
    print('Entropy: {}'.format(entropy))

    if np.all(np.isfinite([np.sum(F_pair), enum])) and np.all(np.isfinite(enums)):
        np.savez(normpath('C:/Users/Owner/PycharmProjects/Research/Publication Data/Spherical Averaging/{a0}{a1}beta{a2}{a3}{a4}g0inv{a5}.npz').format(a0=atom, a1=n, a2=beta_prec_sav, a3=ver, a4=ver_field, a5=np.around(g0inv, decimals=2)), rdensshell=rdensshell, wb=wold_temp)  # Save results to a file.

    gend = time()
    print('Gaussian Computation Time: {}'.format(gend - gtime))  # Total time taken for only SCFT algorithm to run.
    for k in range(Nsk):
        print('Shell {a0} Electron Number: {a1}'.format(a0=k+1, a1=enums[k]))

    print('Atomic Shell Configuration: {a0}'.format(a0=repr(Ns)))
    print('Free Energy (N={a5}, beta={a0}, g0inv={a6}, l_num={a1}, c=[{a2}, {a3}]) : {a4}'.format(a0=np.around(b, decimals=3), a1=[l_num[i] for i in range(Lsk)], a2=[cvec[i][0] for i in range(Lsk)], a3=[cvec[i][-1] for i in range(Lsk)], a4=free, a5=N, a6=np.around(g0inv, decimals=2)))
    print('Convergence: {}'.format(devtot[0]))
    print('Real Space Convergence: {}'.format(realconv))
    print('Electron Number: {}'.format(enum))

    return rdensshell, realns, realn, realwshell, realwnew, realwold, devtot, free, enum, enums, cvec, realconv, temp, cvt


def anderson(dev, anders):  # Anderson matrix and vector definitions.

    umat = np.zeros((anders, anders))
    vvec = np.zeros(anders)

    for m1 in range(anders):
        for m2 in range(anders):
            umat[m1, m2] = np.sum((dev[0]-dev[m1+1])*(dev[0]-dev[m2+1]))

        vvec[m1] = np.sum((dev[0]-dev[m1+1])*dev[0])

    avec = solve(umat, vvec)

    return avec


angmom = False  # Controls whether you want angular dependence in calculations and plotting.
auto = False  # Controls whether you want to plot things or not.
pairs = True  # Controls whether the pairs are populated according to a maximum of 2 or according to known atomic shell occupancies (e.g. 2, 8, 18, etc.).
current_atom = False  # Toggle on to use input data from a previous run for the current atom.

r = sy.symbols('r:3', positive=True, real=True, seq=True)  # Symbolic spherical coordinates for use with sympy.

db_2 = [0]  # Change in atomic orbital basis number.
# Set physical parameters.
A_nuc = np.array([6])  # Atomic numbers for each atom in the molecule.
Nshell = np.array([2, 2, 2])  # Electron arrangement for each Pauli pair.
Nsk = np.size(Nshell)  # Number of Pauli pairs.
N = np.sum(Nshell)  # Total number of electrons in the atom/molecule.
Atomnum = [6, 6]  # Range of atomic numbers to iterate through in the plotting loop.
l = np.array([0])  # Angular momentum numbers for spherical harmonics.
m = [[j for j in range(-i, i+1)] for i in l]  # …and a list of lists for their corresponding m values.
Lsk = np.size(l)  # Total number of atomic orbitals.
Msk = [np.size(m[i]) for i in range(Lsk)]  # Total number of m values.
l_num = np.array([175])+np.array(db_2)  # Number of basis functions for each atomic orbital.
base_num = np.array([np.sum(Msk*l_num) for j in range(Nsk)])  # Number of basis functions for each Pauli pair.
lvt = np.hstack(np.array([np.tile(l[i], l_num[i]*Msk[i]) for i in range(Lsk)]))  # Stacked angular momentum numbers for spherical harmonics corresponding to the total number of basis functions.
mvt = np.hstack(np.array([np.hstack([np.tile(m[i][j], l_num[i]) for j in range(Msk[i])]) for i in range(Lsk)]))  # Stacked m values corresponding to the total number of basis functions.
FA = np.array([[1, -np.ones(Nsk), 1, 0, 'Fukui'], [0, -1/Nshell, 0, 1, 'FAshell'], [0, -(1/N)*np.ones(Nsk), 1, 0, 'FA']])  # The first element is the Fukui function exchange-correlation field, the second is Fermi-Amaldi shell, and the third is standard Fermi-Amaldi.
FAval = 1  #  Decides which exchange-correlation function to use.
betaval = np.array([600.0])  # Beta (1/k_B*T) values.
db = [0.0, [-np.sum([Msk[i]*db_2[i] for i in range(Lsk)]) for j in range(Nsk)], db_2]  # When importing data for SCFT run, these values compensate for any change in beta, Pauli pair basis number, or atomic orbital basis number respectively.
rc = 10.0  # Radial plot length cutoff.
rscale = 1*10**(0)  # scale factor for radial coordinate.
cmin = [-16]  # Minimum values in each atomic orbital for argument of Gaussian exponent.
cmax = [12]  # Maximum values in each atomic orbital for argument of Gaussian exponent.
itermax = 600  # Maximum number of self-consistent iterations.
tol = 5.0*10**(-8)  # Spectral Convergence tolerance.
tol2 = 5.0*10**(-7)  # Spectral tolerance to start using the fine convergence criterion.
tolvals = tol*np.array([1, 1, 10, 10, 10, 10] + [100 for i in range(7, 13)] + [1000 for i in range(13, 87)] + [10000 for i in range(87, 124)])  # Custom spectral convergence tolerance values for each element.
tol2vals = tol2*np.array([1, 1, 10, 10, 10, 10] + [100 for i in range(7, 13)] + [1000 for i in range(13, 87)] + [10000 for i in range(87, 124)])  # Custom fine spectral convergence tolerance values for each element.
realint = 100000#250  # Number of radial grid points for numerical integration.
realint1 = 75  # Number of theta grid points for numerical integration.
realint2 = 100  # Number of phi grid points for numerical integration.
mix = 0.3  # Picard iteration mixing parameter. 0.1-0.5 is the typical range considered.
andersval = 10  # Guess + Anderson histories count.
andtol = 10**(-1)  # When to try an Anderson step.
devtol = 0.1  # Deviation tolerance to try an Anderson step.
rand = [-1, 1]  # Range of stochastic field values.
rand1 = [0, 10**10]  # Range of seed value for random number generator of stochastic field.
isorandnum = 400  # Number of isosurface values to plot.
isopeak = 0.4  # Isosurface lower bound multiplier corresponding to percentage of the peak density value.
seedval = 474747  # Seed for random number generator that determines the seed for the random number generator of the stochastic field.
rng1 = np.random.default_rng(seed=seedval)  # Random number generator that determines the seed for the random number generator of the stochastic field.
seedval1 = rng1.integers(low=rand1[0], high=rand1[1], size=1)  # Seed value for random number generator of stochastic field.
rng = np.random.default_rng(seed=seedval1)  # Random number generator for stochastic field.
switch = [0.0, 1]  # First element controls stochastic field strength.
ind_sh = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]  # Controls which Pauli pairs the perturbation fields are added to.
ind_ang = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]  # Controls which atomic orbitals the perturbation fields are added to.
step = 0.05  # Step size decrease for each iteration of perturbation fields.
step_mat = np.zeros(itermax+1)  # Step size array with as many entries as the maximum number of iterations.
step_mat[0:int(1/step)] = np.flip(np.linspace(0, 1-step, num=int(1/step)))  # Add the step size decrease values from 1 - (step size) to 0.

g0inv = 6.0  # Density of states parameter for Pauli potential.
g0inv_step = 4.0  # Difference between the current g0inv value and the one from the imported initial field data.

# Set plotting parameters.
dentype = 0  # Controls which quantity is plotted. A series of if statements below display what quantity each number represents.
exptype = 'Even-Tempered'  # A plotting title label to specify how the exponents were generated.
atoms = {1: 'Hydrogen', 2: 'Helium', 3: 'Lithium', 4: 'Beryllium', 5: 'Boron', 6: 'Carbon', 7: 'Nitrogen', 8: 'Oxygen',
         9: 'Fluorine', 10: 'Neon', 11: 'Sodium', 12: 'Magnesium', 13: 'Aluminum', 14: 'Silicon', 15: 'Phosphorus',
         16: 'Sulfur', 17: 'Chlorine', 18: 'Argon', 19: 'Potassium', 20: 'Calcium', 21: 'Scandium', 22: 'Titanium',
         23: 'Vanadium', 24: 'Chromium', 25: 'Manganese', 26: 'Iron', 27: 'Cobalt', 28: 'Nickel', 29: 'Copper',
         30: 'Zinc', 31: 'Gallium', 32: 'Germanium', 33: 'Arsenic', 34: 'Selenium', 35: 'Bromine', 36: 'Krypton',
         37: 'Rubidium', 38: 'Strontium', 39: 'Yttrium', 40: 'Zirconium', 41: 'Niobium', 42: 'Molybdenum', 43: 'Technetium',
         44: 'Ruthenium', 45: 'Rhodium', 46: 'Palladium', 47: 'Silver', 48: 'Cadmium', 49: 'Indium', 50: 'Tin',
         51: 'Antimony', 52: 'Tellurium', 53: 'Iodine', 54: 'Xenon', 55: 'Cesium', 56: 'Barium', 57: 'Lanthanum',
         58: 'Cerium', 59: 'Praseodymium', 60: 'Neodymium', 61: 'Promethium', 62: 'Samarium', 63: 'Europium', 64: 'Gadolinium',
         65: 'Terbium', 66: 'Dysprosium', 67: 'Holmium', 68: 'Erbium', 69: 'Thulium', 70: 'Ytterbium', 71: 'Lutetium',
         72: 'Hafnium', 73: 'Tantalum', 74: 'Tungsten', 75: 'Rhenium', 76: 'Osmium', 77: 'Iridium', 78: 'Platinum', 79: 'Gold',
         80: 'Mercury', 81: 'Thallium', 82: 'Lead', 83: 'Bismuth', 84: 'Polonium', 85: 'Astatine', 86: 'Radon', 87: 'Francium',
         88: 'Radium', 89: 'Actinium', 90: 'Thorium', 91: 'Protactinium', 92: 'Uranium', 93: 'Neptunium', 94: 'Plutonium',
         95: 'Americium', 96: 'Curium', 97: 'Berkelium', 98: 'Californium', 99: 'Einsteinium', 100: 'Fermium',
         101: 'Mendelevium', 102: 'Nobelium', 103: 'Lawrencium', 104: 'Rutherfordium', 105: 'Dubnium', 106: 'Seaborgium',
         107: 'Bohrium', 108: 'Hassium', 109: 'Meitnerium', 110: 'Darmstadtium', 111: 'Roentgenium', 112: 'Copernicium',
         113: 'Nihonium', 114: 'Flerovium', 115: 'Moscovium', 116: 'Livermorium', 117: 'Tennessine', 118: 'Oganesson',
         119: 'Ununennium', 120: 'Unbinilium', 121: 'Unbiunium', 122: 'Unbibium', 123: 'Unbitrium'}
atom = atoms[N]  # A plotting title label to specify which atom is being plotted.
# Fixed values for plotting the electron density with respect to two of three variables only.
pcoord = 'nrtp'  # Specifies which type of angular plot is to be displayed. 'nrtp' is for radial only plots, 'rtp' is for 3D isosurface plots, 'rp' is for fixed theta plots, 'rt' is for fixed phi plots, and 'tp' is for fixed radius plots.
rr = np.linspace(10**(-6), 10, 100000)  # Radial values.
tt = np.pi*np.linspace(0, 1, 400)  # Polar angular values.
pp = np.pi*np.linspace(0, 2, 400)  # Azimuthal angular values.

# Plotting labels for each type of quantity.
if dentype == 0:
    denstr = 'Density'
    z_denstr = r'$n(r)$'
elif dentype == 1:
    denstr = 'Radial Electron Density'
    z_denstr = r'$4\pi r^2n(r)$'
elif dentype == 2:
    denstr = 'External Field'
    z_denstr = r'External Field $w(r)$'
elif dentype == 3:
    denstr = 'Radial External Field'
    z_denstr = r'$4\pi r^2 w(r)$'
elif dentype == 4:
    denstr = 'Radial Density Distribution'
    z_denstr = r'$D(r)$'

#figtitle = '{a0} for {a1} Using Gaussian Basis Functions'.format(a0=denstr, a1=atom)
figtitle = '{a0} for {a1}'.format(a0=denstr, a1=atom)
#figtitle = '{a0} for {a1} to {a2} Using Gaussian Basis Functions'.format(a0=denstr, a1=atoms[Atomnum[0]], a2=atoms[Atomnum[1]])

# More plotting parameters for surface plots of the functions above holding one variable fixed.
if angmom and (pcoord != 'rtp') and (dentype != 4):

    if pcoord == 'rt':
        constvals = np.pi*np.linspace(0, 2, num=4)
        radtype = r'$\phi$'
        rvec1, rvec2 = np.meshgrid(rr, tt)
    elif pcoord == 'rp':
        constvals = np.pi*np.linspace(0.1, 0.9, num=4)
        radtype = r'$\theta$'
        rvec1, rvec3 = np.meshgrid(rr, pp)
    elif pcoord == 'tp':
        constvals = np.linspace(0.1, 2, num=4)
        radtype = 'r'
        rvec2, rvec3 = np.meshgrid(tt, pp)

    #fig = plt.figure(figsize=(8, 6))
    #ax = plt.axes(projection='3d')
    #ax = plt.axes()
    #ax.set_title(figtitle)
    #ax.set_xlabel('x')
    #coord_label = 'y'
    #ax.set_ylabel('{}'.format(coord_label))
    #ax.set_zlabel(z_denstr)

    fig, axs = plt.subplots(Nsk, len(constvals), constrained_layout=True)
    #fig, axs = plt.subplots(np.shape(constvals)[0], np.shape(constvals)[1], constrained_layout=True)
    #fig.suptitle(figtitle)

# Isosurface plot grid coordinate values.
elif angmom and pcoord == 'rtp':
    rvec1, rvec2, rvec3 = np.meshgrid(rr, tt, pp)

    X = rvec1*np.sin(rvec2)*np.cos(rvec3)
    Y = rvec1*np.sin(rvec2)*np.sin(rvec3)
    Z = rvec1*np.cos(rvec2)
else:
    fig, ax = plt.subplots()

    rvec1, rvec2, rvec3 = rr, np.pi/2, 0

    ax.set_title(figtitle)
    ax.set_xlabel('r (a.u.)')
    ax.set_ylabel(z_denstr)

iterval = 0

for anum in range(Atomnum[0], Atomnum[1]+1):

    c = np.array([[cmin[i], cmax[i]] for i in range(Lsk)])  # Array of min and max Gaussian exponent pairs for each angular momentum value.
    rdensshell, realns, realn, realwshellval, realwnewval, realwoldval, devtot, fe, enum, enums, cvec, realconv, dfunc, cvt = GaussDensity(r, c, base_num, Nshell, l, m, l_num, betaval[0], db, itermax, tolvals[A_nuc[0]-1], realint, realint1, realint2, g0inv, mix, andersval, andtol, angmom, atoms, rand, switch, ind_ang, ind_sh, step_mat)

    if np.any(pcoord == np.array(['rt', 'rp', 'tp'])):
        cust_ind = [0, 1, 2]
        for k in range(Nsk):
            for i in range(len(constvals[:])):
                if pcoord == 'rt':
                    rvec3 = constvals[i]
                elif pcoord == 'rp':
                    rvec2 = constvals[i]
                elif pcoord == 'tp':
                    rvec1 = constvals[i]
                coordfunc1, coordfunc2 = rvec1*np.sin(rvec2)*np.cos(rvec3), rvec1*np.sin(rvec2)*np.sin(rvec3)
                realnval = realn(rvec1, rvec2, rvec3)
                realnvals = realns(rvec1, rvec2, rvec3)
                plotfunc1 = realnvals[k]
                ax = axs[k, i]
                axs[0, i].set_title('{a0}={a1}'.format(a0=radtype, a1=np.around(constvals[i], decimals=2)))
                axs[k, 0].set_ylabel('Pair {a0}'.format(a0=k+1))
                surf = ax.pcolormesh(coordfunc1, coordfunc2, plotfunc1, shading='gouraud', cmap='viridis', antialiased=True)
                fig.colorbar(surf, ax=ax)
        plt.show()
    elif auto:
        pass
    else:
        realnval = realn(rvec1, rvec2, rvec3)
        realnvals = realns(rvec1, rvec2, rvec3)
        realwshellvals = realwshellval(rvec1, rvec2, rvec3)

    if pairs:
        if np.any(A_nuc[0] == np.arange(2, 19, step=2)):
            Nshell = np.concatenate((Nshell, np.array([1])))
            base_num = np.concatenate((base_num, np.array([base_num[0]])))
        else:
            Nshell[-1] = Nshell[-1] + 1
    else:
        if np.any(A_nuc[0] == np.array([2, 10, 18, 36, 46, 54, 86, 118])):
            Nshell = np.concatenate((Nshell, np.array([1])))
            base_num = np.concatenate((base_num, np.array([base_num[0]])))
        elif A_nuc[0] == 45:
            Nshell = Nshell[:Nsk-1]
            Nshell[-1] = Nshell[-1] + 2
            db[1] = db[1][0:Nsk-1]
            base_num = base_num[0:Nsk-1]
        elif A_nuc[0] == 1:
            Nshell[0] = Nshell[0] + 1
        elif np.any(A_nuc[0] == np.array([i for i in range(3, 10)])):
            Nshell[1] = Nshell[1] + 1
        elif np.any(A_nuc[0] == np.array([i for i in range(11, 18)] + [20, 21, 22, 25, 26, 27])):
            Nshell[2] = Nshell[2] + 1
        elif np.any(A_nuc[0] == np.array([19, 24] + [i for i in range(29, 36)] + [38, 39, 41, 44, 57, 59, 60, 61, 62, 65, 66, 67, 68, 69])):
            Nshell[3] = Nshell[3] + 1
        elif np.any(A_nuc[0] == np.array([37, 42] + [i for i in range(47, 54)] + [56, 63] + [i for i in range(70, 77)] + [78, 91, 92, 94, 97, 98, 99, 100, 101, 120, 121, 122])):
            Nshell[4] = Nshell[4] + 1
        elif np.any(A_nuc[0] == np.array([55] + [i for i in range(79, 86)] + [88, 89, 95] + [i for i in range(102, 112)])):
            Nshell[5] = Nshell[5] + 1
        elif np.any(A_nuc[0] == np.array([87] + [i for i in range(112, 118)])):
            Nshell[6] = Nshell[6] + 1
        elif A_nuc[0] == 119:
            Nshell[7] = Nshell[7] + 1
        elif np.any(A_nuc[0] == np.array([23, 28])):
            Nshell[3] = Nshell[3] - 1
            Nshell[2] = Nshell[2] + 2
        elif np.any(A_nuc[0] == np.array([40, 43, 58, 64])):
            Nshell[4] = Nshell[4] - 1
            Nshell[3] = Nshell[3] + 2
        elif np.any(A_nuc[0] == np.array([77, 90, 93, 96])):
            Nshell[5] = Nshell[5] - 1
            Nshell[4] = Nshell[4] + 2

    if not auto:
        if dentype == 0:
            plotfunc = realnval
            textstr = r'$n(r)$'
        elif dentype == 1:
            plotfunc = 4.0*np.pi*rvec1**2*realnval
            textstr = r'$4\pi r^2n(r)$'
        elif dentype == 2:
            plotfunc = realwnewval
            textstr = r'$w(r)$'
        elif dentype == 3:
            plotfunc = 4.0*np.pi*rvec1**2*realwnewval
            textstr = r'$4\pi r^2 w(r)$'
        elif dentype == 4:
            plotfunc = (2**11*np.pi)**0.25*rvec1**2*sy.lambdify(r, np.sum(np.array([np.sum(rdensshell, axis=0)[i]*cvec[0][i]**0.75*sy.exp(-cvec[0][i]*r[0]**2) for i in range(l_num[0])])))(rvec1, rvec2, rvec3)
            textstr = r'$D(r)$'

        if pcoord != 'rtp':
            if angmom and (dentype != 4):
                ax.set_xlim((-rc, rc))
                ax.set_ylim((-rc, rc))
            else:
                ax.set_xlim((0, rc))
                #ax.plot(rr, plotfunc, label='{a0} with {a1} Exponents and {a2} Basis Functions'.format(a0=textstr, a1=exptype, a2=base_num))
                ax.plot(rr, plotfunc, label='{a0}'.format(a0=atom))

            for i in range(Nsk):
                if dentype == 0:
                    plotfunc1 = realnvals[i]
                    textstr1 = r'$n_{}(r)$'.format(i+1)
                elif dentype == 1:
                    plotfunc1 = 4.0*np.pi*rvec1**2*realnvals[i]
                    textstr1 = r'$4\pi r^2n_{}(r)$'.format(i+1)
                elif dentype == 2:
                    plotfunc1 = realwshellvals[i]
                    textstr1 = r'$w_{}(r)$'.format(i+1)
                elif dentype == 3:
                    plotfunc1 = 4.0*np.pi*rvec1**2*realwshellvals[i]
                    textstr1 = r'$4\pi r^2w_{}(r)$'.format(i+1)
                elif dentype == 4:
                    plotfunc1 = (2**11*np.pi)**0.25*rvec1**2*sy.lambdify(r, np.sum(np.array([rdensshell[i][j]*cvec[0][j]**0.75*sy.exp(-cvec[0][j]*r[0]**2) for j in range(l_num[0])])))(rvec1, rvec2, rvec3)
                    textstr1 = r'$D_{}(r)$'.format(i+1)

                if angmom and (dentype != 4):
                    ax.set_xlim((-rc, rc))
                    ax.set_ylim((-rc, rc))
                else:
                    ax.set_xlim((0, rr[-1]))
                    #ax.plot(rr, plotfunc1, label='{a0} with {a1} exponents and {a2} basis functions'.format(a0=textstr1, a1=exptype, a2=base_num[i]))
                    #ax.plot(rr, plotfunc1, label='{a0}'.format(a0=textstr1))

            if N == 1:
                if dentype == 0:
                    plotfunc_known = np.exp(-2.0*rvec1)/np.pi
                elif np.any(dentype == np.array([1, 4])):
                    plotfunc_known = 4.0*rvec1**2*np.exp(-2.0*rvec1)
                elif dentype == 2:
                    plotfunc_known = -1/rvec1
                elif dentype == 3:
                    plotfunc_known = -4.0*np.pi*rvec1

                if angmom and (dentype != 4):
                    ax.set_xlim((-rc, rc))
                    ax.set_ylim((-rc, rc))
                else:
                    ax.set_xlim((0, rc))
                    #ax.plot(rr, plotfunc_known, label='Known Solution')

        elif angmom and pcoord == 'rtp':

            Isosurf = np.append(np.linspace(np.min(plotfunc), isopeak*np.max(plotfunc), int(0.1*isorandnum)), np.linspace(isopeak*np.max(plotfunc), np.max(plotfunc), isorandnum))
            grid = pyv.StructuredGrid(X, Y, Z)  # Make an empty uniform grid of the right size.
            grid.point_data["Isosurfaces"] = plotfunc.ravel(order='F')  # Populate the grid with my data and add a label.
            contours = grid.contour(Isosurf)#.scale([10.0, 10.0, 10.0])  # Choose which isosurfaces to display.
            pyv.set_plot_theme('document')  # Set plot background to white.
            plotname = 'tot_plot'
            myplot = pyv.Plotter()  # Open an empty 3D plot.
            myplot.add_axes()
            myplot.add_title(figtitle)
            myplot.add_mesh(contours, opacity=0.5, clim=[np.min(plotfunc), np.max(plotfunc)], name=plotname)  # Plot isosurfaces.
            myplot.show_bounds(bounds=[-rc, rc, -rc, rc, -rc, rc], location='outer', grid='front', xlabel='x', ylabel='y', zlabel='z', minor_ticks=True)
            myplot.show()

            for i in range(Nsk):
                if dentype == 0:
                    plotfunc1 = realnvals[i]
                    textstr1 = r'$n_{}(r)$'.format(i+1)
                elif dentype == 1:
                    plotfunc1 = 4.0*np.pi*rvec1**2*realnvals[i]
                    textstr1 = r'$4\pi r^2n_{}(r)$'.format(i+1)
                elif dentype == 2:
                    plotfunc1 = realwshellval[i]
                    textstr1 = r'$w_{}(r)$'.format(i+1)
                elif dentype == 3:
                    plotfunc1 = 4.0*np.pi*rvec1**2*realwshellval[i]
                    textstr1 = r'$4\pi r^2w_{}(r)$'.format(i+1)

                Isosurf1 = np.append(np.linspace(np.min(plotfunc1), isopeak*np.max(plotfunc1), int(0.1*isorandnum)), np.linspace(isopeak*np.max(plotfunc1), np.max(plotfunc1), isorandnum))
                grid1 = pyv.StructuredGrid(X, Y, Z)  # Make an empty uniform grid of the right size.
                grid1.point_data["Isosurfaces"] = plotfunc1.ravel(order='F')  # Populate the grid with my data and add a label.
                contours1 = grid1.contour(Isosurf1)  # Choose which isosurfaces to display.
                pyv.set_plot_theme('document')  # Set plot background to white.
                name = 'Shell {}'.format(i)
                myplot1 = pyv.Plotter()  # Open an empty 3D plot.
                myplot1.add_axes()
                myplot1.add_title('Shell {a0} {a1} for {a2} Using Gaussian Basis Functions'.format(a0=i+1, a1=denstr, a2=atom))
                myplot1.add_mesh(contours1, opacity=0.5, clim=[np.min(plotfunc1), np.max(plotfunc1)], name=name)  # Plot isosurfaces.
                myplot1.show_bounds(bounds=[-rc, rc, -rc, rc, -rc, rc], location='outer', grid='front', xlabel='x', ylabel='y', zlabel='z', minor_ticks=True)
                myplot1.show()

            if N == 1:
                if dentype == 0:
                    plotfunc_known = np.exp(-2.0*rvec1)/np.pi
                elif dentype == 1:
                    plotfunc_known = 4.0*rvec1**2*np.exp(-2.0*rvec1)
                elif dentype == 2:
                    plotfunc_known = -1/rvec1
                elif dentype == 3:
                    plotfunc_known = -4.0*np.pi*rvec1

                Isosurf2 = np.append(np.linspace(np.min(plotfunc_known), isopeak*np.max(plotfunc_known), int(0.1*isorandnum)), np.linspace(isopeak*np.max(plotfunc_known), np.max(plotfunc_known), isorandnum))
                grid2 = pyv.StructuredGrid(X, Y, Z)  # Make an empty uniform grid of the right size.
                grid2.point_data["Isosurfaces"] = plotfunc_known.ravel(order='F')  # Populate grid with my data and add a label.
                contours2 = grid2.contour(Isosurf2)  # Choose which isosurfaces to display.
                pyv.set_plot_theme('document')  # Set plot background to white.
                myplot2 = pyv.Plotter()  # Open an empty 3D plot.
                myplot2.add_axes()
                myplot2.add_title(figtitle)
                myplot2.add_mesh(contours2, opacity=0.5, clim=[np.min(plotfunc_known), np.max(plotfunc_known)])  # Plot isosurfaces.
                myplot2.show_bounds(bounds=[-rc, rc, -rc, rc, -rc, rc], location='outer', grid='front', xlabel='x', ylabel='y', zlabel='z', minor_ticks=True)
                myplot2.show()

        end = time()
        print('Computation Time: {}'.format(end - start))

        if pcoord != 'rtp':
            ax.legend(loc='best')
            plt.show()
    else:
        pass

    Nsk = np.size(Nshell)
    N = np.sum(Nshell)
    FA[FAval][1] = -1/Nshell
    db[1] = [db[1][0] for j in range(Nsk)]
    db[0] = 0.0
    atom = atoms[N]
    A_nuc[0] = A_nuc[0]+1

if pcoord != 'rtp':
    ax.legend(loc='best')
    plt.show()
