#!/usr/bin/env python3

import numpy as np
import scipy.special
pi = np.pi
from numpy import log, exp, sqrt, sin, cos
import argparse
import h5py
import fio

# Note to self: setting nu_inner=nu_outer=0 gives QUARTER-MASS untapered Mestel disc!
class TaperedMestel:
    "Tapered Mestel disc a la Fouvry et al (2015) assuming vcirc=1"
    def __init__(self, q=None, sigma=None, Rinner=1.0, Router=11.5, nu_inner=4, nu_outer=5):
        if q is not None and sigma is None: sigma = 1/sqrt(1+q)
        elif q is None and sigma is not None: q = 1/sigma**2-1
        else:
            raise Exception("You need to set ONE of 'q' or 'sigma' in TaperedMestel!")
        self.q, self.sigma = q, sigma
        self.Li, self.Lo = Rinner, Router # assumes V_c=1
        self.nui = nu_inner
        self.nuo = nu_outer
        self.norm = 1.0/(2**(q/2)*sqrt(pi)*sigma**(q+2)*scipy.special.gamma((1+q)/2)*2*pi)

    def F(self,R,vr,vphi):
        q, nui, nuo = self.q, self.nui, self.nuo
        L = R*vphi
        E = log(R)+(vr**2+vphi**2)/2
        return self.norm*L**q*exp(-E/self.sigma**2)*\
            L**nui/(self.Li**nui+L**nui)*self.Lo**nuo/(self.Lo**nuo+L**nuo)
    def density(self,R,nv=1000,vmax=5.0):
        q, nui, nuo = self.q, self.nui, self.nuo
        dv = vmax/nv
        vphis = np.arange(0,nv)*dv
        Ls = R*vphis
        return self.norm*(
            Ls**q*Ls**nui/(self.Li**nui+Ls**nui)*self.Lo**nuo/(self.Lo**nuo+Ls**nuo)
            *exp(-vphis**2/(2*self.sigma**2)) ).sum()*dv*\
            sqrt(2*pi)*self.sigma*exp(-log(R)/self.sigma**2)
    def drawv(self,R):
        "Draw and return (vR,vphi)"
        q, nui, nuo = self.q, self.nui, self.nuo
        vR = np.random.normal(0.0,self.sigma)
        while True:
            vphi = sqrt(2*np.random.gamma(shape=(1+q)/2,scale=1))*self.sigma
            L = R*vphi
            fref = L**nui/(self.Li**nui+L**nui)*self.Lo**nuo/(self.Lo**nuo+L**nuo)
            f = np.random.uniform(0,1)
            if f<fref: break
        return (vR,vphi)
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-N', type=int,default=1000000)
    parser.add_argument('-q', type=float,default=11.4)
    parser.add_argument('-seed', type=int,default=0)
    #parser.add_argument('--fractive', type=float,default=1.0)
    parser.add_argument('output')
    args = parser.parse_args()
    N = args.N
    q = args.q
    seed = args.seed
    #fractive = args.fractive
    fnam = args.output
 
    m1 = TaperedMestel(q=q)
    log10Rs = np.arange(-2,2,0.001)
    Rs = 10**log10Rs
    rho1s = np.array([ m1.density(R) for R in Rs])

    # Calculate total mass
    dlogR = log(Rs[1])-log(Rs[0])
    M = 2*pi*(Rs**2*rho1s).sum()*dlogR
    #M *= fractive
    print("Total mass =",M)

    # Now let's generate the N-body realisation.
    # Crude is good: sample directly from these rho1s.
    ps = rho1s*Rs**2
    ps/=ps.sum() # Now ps[i] is probability of finding star at Rs[i]
    xs,vs = np.zeros((N,2)), np.zeros((N,2))
    # Fixing the random seed
    np.random.seed(seed)
    print("Sampling...")
    for n in range(N):
        print(f"\r {n+1} / {N}",end="")
        R = np.random.choice(Rs,p=ps)
        phi = np.random.uniform(0,2*pi)
        xs[n] = R*cos(phi), R*sin(phi)
        vR,vphi = m1.drawv(R)
        vs[n] = vR*cos(phi)-vphi*sin(phi), vR*sin(phi)+vphi*cos(phi)
    print("") # newline
    ms = np.ones(N)*M/N
    xs = np.array(xs)
    vs = np.array(vs)

    fio.dump_ptle(ms,xs,vs,fnam,"/",'w')
