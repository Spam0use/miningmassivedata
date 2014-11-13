#!/usr/bin/env python

import numpy as np
import pdb

class adjmat(object):
    """
    adjacency matrix
    """
    def __init__(self,nparray=np.ones((3,3)),deadends=False):
        self.mat=nparray/nparray.sum(0,dtype='float')  #make stochastic
        self.nitm=self.mat.shape[0]
        if deadends:
            #replace nans of dead end columns with 1/nitm
            self.mat[np.isnan(self.mat)]=1.0/self.nitm
        else:
            #replace w/0
            self.mat[np.isnan(self.mat)]=0


    def pru(self,rank=None,beta=1.0):
        """
        page rank update
        """
        nr=self.nitm
        if rank is None:
            rank=np.ones((nr,1))/float(nr)
        #rank1 = beta*(self.mat.dot(rank))
        #rank2 = (1-beta)*np.ones((nr,1))/float(nr)
        #rank = rank1+rank2
        
        #a=beta*self.mat + (1-beta)*np.ones((nr,nr))/float(nr)
        #rank=a.dot(rank)
        
        #wk1 final
        #r1=beta*self.mat.dot(rank)
        #sr1=sum(r1)
        #rank=r1+(1-sr1)/nr
        
        r1=beta*self.mat.dot(rank)
        r2=r1+(1-beta)/nr
        rank=r2/sum(r2)
        return rank

    def pr(self,beta=1.0):
        """ iterate pru to convergence """
        r=np.ones((self.nitm,1))/float(self.nitm)
        rnew=self.pru(r,beta=beta)
        d=np.sqrt(sum((r-rnew)**2))/sum(rnew)
        r=rnew.copy()
        while(d>0.001):
            rnew=self.pru(r,beta)
            d=np.sqrt(sum((r-rnew)**2))/sum(rnew)
            r=rnew.copy()
        return r

def q1():
    #pdb.set_trace()
    q1mat=np.array([0,1,0,0,1,0,0,0,1,0,0,1,0,0,1,0]).reshape((4,4))  #in edges
    q1adj=adjmat(q1mat)
    q1r=q1adj.pr(beta=0.7)
    #q1r=q1r*q1adj.nitm/sum(q1r)
    return q1r