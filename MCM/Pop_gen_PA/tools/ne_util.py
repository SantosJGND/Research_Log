
import numpy as np


def theta_constant(tnow ,Ne= 1000):
    '''
    constant Ne
    '''
    if type(tnow).__module__ == np.__name__:
        Narray= np.ones(tnow.shape, dtype=int)
        Narray= Narray * Ne
        return Narray
    else:
        return Ne

def theta_exp(tnow, Ne= 1000, rate= 0.038):
    '''
    exponential growth or decline.
    '''
    if type(tnow).__module__ == np.__name__:
        
        Narray= Ne * np.exp(tnow * rate)
        return Narray
    
    else:
        Nt= Ne * np.exp(tnow * rate)
        return Nt

def theta_function(tnow, theta_time_array= [], inverse= False):
    '''
    pre-determined Ne scale. 
    '''
    if len(theta_time_array) == 0:
        print('no theta_time_array.')
        return tnow
    
    prime= np.where(theta_time_array[:,0] > tnow)[0]
    
    if len(prime) == 0:
        prime= theta_time_array[-1,1]
        
        if inverse:
            prime= 1 / prime
        
        return prime
    
    prime= theta_time_array[prime[0],1]
    
    if inverse:
        prime= 1 / prime
    
    return prime