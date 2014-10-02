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
        #a=beta*self.mat + (1-beta)*np.ones((nr,nr))/float(nr)
        #rank=a.dot(rank)
        r1=beta*self.mat.dot(rank)
        sr1=sum(r1)
        rank=r1+(1-sr1)/nr
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
    #q1mat=np.array([0,1,1,0,0,1,0,0,1]).reshape((3,3)) #out edges
    q1mat=np.array([0,0,0,1,0,0,1,1,1]).reshape((3,3))  #in edges
    q1adj=adjmat(q1mat)
    q1r=q1adj.pr(beta=0.7)
    q1r=q1r*q1adj.nitm/sum(q1r)
    return q1r

def q2():
    mat=np.array([0,0,1,1,0,0,1,1,0]).reshape((3,3))
    adj=adjmat(mat)
    r=adj.pr(beta=0.85)
    print r/r[0]
    print r/r[1]
    print r/r[2]
    return r

def q3():
    mat=np.array([0,0,1,1,0,0,1,1,0]).reshape((3,3))
    adj=adjmat(mat)
    out=[]
    r=np.ones((adj.nitm,1))/float(adj.nitm)
    out.append(r.copy()*3)
    for i in xrange(5):
        r=adj.pru(r)
        out.append(r.copy()*3)
    r=adj.pr()*3
    out.append(r.copy())
    return out

def primes(maxint):
    out=[2]
    for i in xrange(3,maxint+1):
        prm=True
        for j in out:
            if i % j==0:
                prm=False
        if prm:
            out.append(i)
    return out

def mapdiv(input):
    '''
    map input integers p to the set of prime divisors as (divisor1,p),(divisor2,p)
    '''
    mx=max(input)
    prm=primes(mx)
    out=[]
    for i in input:
        for j in prm:
            if j<=i and i % j==0:
                out.append((j,i))
    return out

def reducediv(input):
    '''
    reduce output of map to (divizor,sum-of-inputs) pairs
    '''
    out={}
    for i in input:
        out[i[0]]=out.get(i[0],0)+i[1]
    return out

def q4():
    return reducediv(mapdiv([15,21,24,30,49]))