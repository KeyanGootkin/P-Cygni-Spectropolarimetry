import numpy as np
import pandas as pd
from sklearn.model_selection import GridSearchCV
from sklearn.neighbors import KernelDensity
from astropy.modeling.functional_models import Gaussian2D
from astropy.modeling import fitting
import matplotlib.pyplot as plt
import argparse
import warnings
import sys
import sympy as sp
from scipy.stats import pearsonr
sp.init_printing()
fitter = fitting.LevMarLSQFitter()


class Ellipse(object):
    def __init__(self,q=None,u=None):
        self.q = q
        self.u = u
        self.muq = self.q.mean()
        self.muu = self.u.mean()
        self.mu = np.array([self.muq,self.muu])
        max_range = np.max([self.q.max()-self.q.min(),self.u.max()-self.u.min()])
        self.X, self.Y = np.meshgrid(np.linspace(self.muq-max_range/2,self.muq+max_range/2,100), 
                                     np.linspace(self.muu-max_range/2,self.muu+max_range/2,100))
        self.cov = (1/self.q.size)*np.sum([np.outer(np.array([Q,U])-self.mu,np.array([Q,U])-self.mu)    
                                           for Q,U in zip(self.q,self.u)],axis=0)
        self.g2D = Gaussian2D(x_mean=self.muq,
                              y_mean=self.muu,
                              cov_matrix=self.cov)
        self.G = self.g2D.evaluate(self.X,
                                  self.Y,
                                  self.g2D.amplitude.value,
                                  self.g2D.x_mean.value,
                                  self.g2D.y_mean.value,
                                  self.g2D.x_stddev.value,
                                  self.g2D.y_stddev.value,
                                  self.g2D.theta.value)
        self.stdq = np.sqrt(np.var(self.q))
        self.stdu = np.sqrt(np.var(self.u))
        self.theta = np.rad2deg(self.g2D.theta.value)
        self.r,self.p = pearsonr(self.q,self.u)
        self.L = 2*self.stdq*self.stdu*np.sqrt(1-self.r**2)/(self.stdq**2+self.stdu**2)
        self.sig = self.L**(self.q.size-2)
        
    def KDE(self):
        data = np.vstack([self.q, self.u]).T
        #Grid search for best KDE bandwidth
        params = {'bandwidth': np.linspace(np.min(np.diff(self.q)),np.max(np.diff(self.u)),100)}
        grid = GridSearchCV(KernelDensity(), params)
        grid.fit(data)

        KDE_bandwidth = grid.best_estimator_.bandwidth

        kde = grid.best_estimator_

        xy = np.vstack([self.X.ravel(), self.Y.ravel()]).T
        #compute the KDE on a 100x100 grid of points
        self.Z = np.exp(kde.score_samples(xy)).reshape(self.X.shape)
        return self.Z
        
    def KDE_Gaussian(self):
        X=self.X
        Y=self.Y
        Z = self.KDE()
        kg2D_init = Gaussian2D(amplitude=np.max(Z), 
                          x_mean=X[np.unravel_index(np.argmax(Z),Z.shape)], 
                          y_mean=Y[np.unravel_index(np.argmax(Z),Z.shape)], 
                          x_stddev=np.std(X), 
                          y_stddev=np.std(Y), 
                          bounds={'theta': (-2*np.pi,2*np.pi),
                                  'x_mean': (np.min(X),np.max(X)),
                                  'y_mean': (np.min(Y),np.max(Y)),
                                  'x_stddev':(0.001,1),
                                  'y_stddev':(0.001,1)})
        kg2D = fitter(kg2D_init, X, Y, Z, weights = Z)
        kde_G = kg2D.evaluate(X,Y,kg2D.amplitude.value,kg2D.x_mean.value,kg2D.y_mean.value,kg2D.x_stddev.value,kg2D.y_stddev.value,kg2D.theta.value)
        return kg2D,kde_G
    
    def KDE_sig(self):
        kg2D, self.kde_G = self.KDE_Gaussian()
        r,p = pearsonr(self.q,self.u)
        L = 2*kg2D.x_stddev.value*kg2D.y_stddev.value*np.sqrt(1-r**2)/(kg2D.x_stddev.value**2+kg2D.y_stddev.value**2)
        sig = L**(self.q.size-2)
        return sig