#!/home/sthalabard/anaconda3/bin/python
#!/impa/home/f/simon.thalabard/anaconda3/bin/python3

#!/usr/bin/env python3

# -*- coding: utf-8 -*-

#--------------------------------#
from math import *
import numpy as np
import scipy as scp
import scipy.optimize 
from scipy import fftpack as fft
from matplotlib import pyplot as plt
from matplotlib.pyplot import *
from matplotlib import animation
from time import time as ctime
import os,glob

from .  import Shell4KH_boost as boost

#Simplified dynamics, without forcing

class GOY:
        def __init__(self,nu=1,kdiss=0,hyper=1,d=2.,\
                    nmin=-12, nmax=10,Lambda=2,dt=1e-3,dtype='complex128'):

            #Parameters of the shell Model
            self.nu=nu      #Viscosity
            self.d=d        #Dimension
            self.nmin=nmin  #minimal shell number
            self.nmax=nmax  #Maximal shell number
            self.Lambda=Lambda #Intershell ratio

            self.n=nmax-nmin+1
            self.dt=dt
            self.hyper=hyper #hyper viscosity
            self.k=np.arange(nmin,nmax+1);

            if d==2:self.gamma=Lambda**2
            else:self.gamma=-Lambda 

            self.Lambdas=self.Lambda**(1.*self.k)
            self.Gammas=self.gamma**(1.*self.k)

            self.kdiss=kdiss
            self.dissip= self.nu*self.Lambda**((self.k-self.kdiss)*(2*self.hyper))

            self.t=0

## Velocity Fields
            self.u=np.zeros(self.n,dtype=dtype)          #First Model

            self.v=np.zeros(self.n,dtype=dtype)          #Second Model

#TO use the boost through f2py                                                            #
            self.dtype=dtype
            self.u_xy=np.zeros((self.n,2),dtype=self.u[0].real.dtype)
            self.v_xy=np.zeros((self.n,2),dtype=self.u[0].real.dtype)


#            self.E=0. #!sum |a|^2
#            self.Z=0.# sum \gamma |a|^2 

                #STAT
#            self.S2=np.zeros_like(self.Gammas)
#            self.S3=np.zeros_like(self.Gammas)
#            self.S4=np.zeros_like(self.Gammas)

#            self.update_laws()

        def display(self):
            print('nu =',self.nu)
            print('kdiss =',self.kdiss)
            print('hyper =',self.hyper)
            print('d =',self.d)
            print('nmin =',self.nmin)
            print('nmax =',self.nmax)
            print('Lambda =',self.Lambda)
            print('dt =',self.dt)
            return None
            
        def update_z2xy(self):
            self.u_xy[:,0],self.u_xy[:,1]=self.u.real[:],self.u.imag[:]
            self.v_xy[:,0],self.v_xy[:,1]=self.v.real[:],self.v.imag[:]
            return None

        def update_xy2z(self):
            self.u[:]=self.u_xy[:,0]+1j*self.u_xy[:,1]
            self.v[:]=self.v_xy[:,0]+1j*self.v_xy[:,1]
            return None
        ######### RK ROUTINES
        def rhs_nl(self,u):
            dudt=np.zeros_like(u);
            ustar=u.conjugate()
            dudt[:-2]=ustar[1:-1]*ustar[2:] *self.Lambda**2*self.gamma#a(k+1)*a(k+2)
            dudt[1:-1]+=ustar[:-2]*ustar[2:]*(-self.Lambda)*(1+self.gamma) #a(k-1)*a(k+1)
            dudt[2:]+=ustar[:-2]*ustar[1:-1] #a(k-2)a(k-1)
            dudt[:]*=self.Lambdas[:]
            return dudt

        def update_RK_nl(self):
            h=self.dt
            y,yc=self.u.copy(),self.v.copy()
            y1,y1c=self.rhs_nl(y),self.rhs_nl(yc)
            y2,y2c=self.rhs_nl(y+0.5*h*y1),self.rhs_nl(yc+0.5*h*y1c)
            self.u=self.u+ h*y2
            self.v=self.v+ h*y2c
            return None

        def update_diss(self):
            self.u  = self.u*np.exp(-self.dissip*self.dt)
            self.v  = self.v*np.exp(-self.dissip*self.dt)
            return None

        def update(self):
            self.update_RK_nl()
            self.update_diss()
            self.t=self.t+self.dt
            return None

        def integrate_until(self,tol=1e9,Tmax=1e9,every=0.1,verbose=False):
            tmp=0
#           cond1 = lambda x: (np.abs(x)**2).sum() <tol
            f = lambda x: (np.abs(x)**2).sum()


            print(self.t,f(self.u-self.v),self.energy_umv())

            while (f(self.u-self.v)<tol) and (self.t<Tmax): 
                self.update()
                if verbose:
                    tmp=tmp+self.dt
                    if tmp>every :
                        print('t : %0.2f \t E_U : %0.2f \t E_V : %0.2f \t CE :%0.2e ' %(self.t, self.energy_u(),self.energy_u(),self.energy_umv()))
                        tmp=0

            print(self.t,f(self.u-self.v),self.energy_umv())

            return None

        def integrate_until_boost(self,tol=1e9,Tmax=10,verbose=False):
            nt=int(np.floor((Tmax-self.t)/self.dt))
            nt=max(nt,0)
            self.update_z2xy()
            dum1,dum2,nstep=boost.update(u=self.u_xy,v=self.v_xy,dt=self.dt,nt=nt,tol=tol,\
                         dissip=self.dissip,gam=self.gamma,lam=self.Lambda,\
                         etas=self.Lambdas)
            if verbose: print(nt,nstep)
            
            self.update_xy2z()
            self.t=self.t+nstep*self.dt

#            print(boost.update.__doc__)
            return None

        ########## INITIALIZATION ROUTINES
        def init_KH(self,a=0,theta=pi/4):
            self.u[:]=0
            self.v[:]=0
            self.u[self.k%3==0]=self.Lambdas[self.k%3==0]**(-a)*(np.cos(theta)+1j*np.sin(theta))
            self.v[:]=self.u[:]
            self.update_z2xy()
            return None

        def perturb(self,kappa=0.1,a=0,seed=None):
            if seed is not None: np.random.seed(seed)
            
            self.u[:]+=self.Lambdas**(-0.)*(self.Lambdas/self.Lambdas[-1])**a\
                        *np.sqrt(kappa/2.)*(1j*np.random.randn(self.n)+np.random.randn(self.n))
            self.v[:]+=self.Lambdas**(-0.)*(self.Lambdas/self.Lambdas[-1])**a\
                        *np.sqrt(kappa/2.)*(1j*np.random.randn(self.n)+np.random.randn(self.n))
            self.update_z2xy()
            return None

        def perturb_smooth(self,kmax=None,kappa=0.1,a=0,seed=None):
            if seed is not None: np.random.seed(seed)
            
            self.u[:]+=self.Lambdas**(-0.5)*(self.Lambdas/self.Lambdas[self.k==kmax])**a\
                        *np.sqrt(kappa/2.)*(1j*np.random.randn(self.n)+np.random.randn(self.n))
            self.v[:]+=self.Lambdas**(-0.5)*(self.Lambdas/self.Lambdas[self.k==kmax])**a\
                        *np.sqrt(kappa/2.)*(1j*np.random.randn(self.n)+np.random.randn(self.n))
            self.update_z2xy()
            return None

        def regularize(self,kmax=None,tol=1e-30):
            if kmax is None: 
                damp=0
            else:
                self.t=-0.5*np.log(tol)/self.dissip[self.k==kmax]
                damp=self.dissip*self.t;
            self.u=self.u*np.exp(-damp)
            self.v=self.v*np.exp(-damp)
            self.update_z2xy()
            return None

        ######### MISCELLANEOUS
        def energy_u(self):
            return (np.abs(self.u)**2).sum()

        def energy_v(self):
            return (np.abs(self.v)**2).sum()
        
        def energy_umv(self):
            return (np.abs(self.u-self.v)**2).sum()

        
