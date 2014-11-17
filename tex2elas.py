#!/usr/bin/env python
#
# Calculate pollycrystal elastic constants 
# from single crystal tensor and Euler angles
# derived from VPSC calculation.
#
# Copyright (C) Andrew Walker, 2010, 2011
#               <andrew.walker@bristol.ac.uk>
import numpy as np
import numpy.linalg as npl

def mat2tens(cij_mat, compl=False):
    """Convert from Voigt to full tensor notation 

       Convert from the 6*6 elastic constants matrix to 
       the 3*3*3*3 tensor representation. Recoded from 
       the Fortran implementation in DRex. Use the optional 
       argument "compl" for the elastic compliance (not 
       stiffness) tensor to deal with the multiplication 
       of elements needed to keep the Voigt and full 
       notation consistant.

    """
    cij_tens = np.zeros((3,3,3,3))
    m2t = np.array([[0,5,4],[5,1,3],[4,3,2]])
    if compl:
        cij_mat = cij_mat / np.array([[1.0, 1.0, 1.0, 2.0, 2.0, 2.0],
                                      [1.0, 1.0, 1.0, 2.0, 2.0, 2.0],
                                      [1.0, 1.0, 1.0, 2.0, 2.0, 2.0],
                                      [2.0, 2.0, 2.0, 4.0, 4.0, 4.0],
                                      [2.0, 2.0, 2.0, 4.0, 4.0, 4.0],
                                      [2.0, 2.0, 2.0, 4.0, 4.0, 4.0]])
    for i in range(3):
        for j in range(3):
            for k in range(3):
                for l in range(3):
                    cij_tens[i,j,k,l] = cij_mat[m2t[i,j],m2t[k,l]]

    return cij_tens

def tens2mat(cij_tens, compl=False):
    """Convert from full tensor to Voigt notation

       Convert from the 3*3*3*3 elastic constants tensor to 
       to 6*6 matrix representation. Recoded from the Fortran
       implementation in DRex. Use the optional 
       argument "compl" for the elastic compliance (not 
       stiffness) tensor to deal with the multiplication 
       of elements needed to keep the Voigt and full 
       notation consistant.

    """
    t2m = np.array([[0,1,2,1,2,0],[0,1,2,2,0,1]])
    cij_mat = np.zeros((6,6))
    # Convert back to matrix form
    for i in range(6):
        for j in range(6):
            cij_mat[i,j] = cij_tens[t2m[0,i],t2m[1,i],t2m[0,j],t2m[1,j]]
    
    if compl:
        cij_mat = cij_mat * np.array([[1.0, 1.0, 1.0, 2.0, 2.0, 2.0],
                                      [1.0, 1.0, 1.0, 2.0, 2.0, 2.0],
                                      [1.0, 1.0, 1.0, 2.0, 2.0, 2.0],
                                      [2.0, 2.0, 2.0, 4.0, 4.0, 4.0],
                                      [2.0, 2.0, 2.0, 4.0, 4.0, 4.0],
                                      [2.0, 2.0, 2.0, 4.0, 4.0, 4.0]])

    return cij_mat

def rotT(T, g):
    """Rotate a rank 4 tensor, T, using a rotation matrix, g
       
       Tensor rotation involves a summation over all combinations
       of products of elements of the unrotated tensor and the 
       rotation matrix. Like this for a rank 3 tensor:

           T'(ijk) -> Sum g(i,p)*g(j,q)*g(k,r)*T(pqr)
       
       with the summation over p, q and r. The obvious implementation
       involves (2*rank) length 3 loops building up the summation in the
       inner set of loops. This optimized implementation >100 times faster 
       than that obvious implementaton using 8 nested loops. Returns a 
       3*3*3*3 numpy array representing the rotated tensor, Tprime. 

    """
    gg = np.outer(g, g) # Flatterns input and returns 9*9 array
                        # of all possible products
    gggg = np.outer(gg, gg).reshape(4 * g.shape)
                        # 81*81 array of double products reshaped
                        # to 3*3*3*3*3*3*3*3 array...
    axes = ((0, 2, 4, 6), (0, 1, 2, 3)) # We only need a subset 
                                        # of gggg in tensordot...
    return np.tensordot(gggg, T, axes)

def rotT_2(T, g):
    """Rotate a rank 2 tensor, T, using a rotation matrix, g

       This uses the obvious (looping) implementation. Scope for 
       speedup (maybe factor of 20?) using scheme in rotT, above.

    """
    Tprime = np.zeros((3,3))
    for i in range(3):
        for j in range(3):
            for ii in range(3):
                for jj in range(3):
                    Tprime[i,j] = Tprime[i,j] + g[ii,i]*g[jj,j]*T[ii,jj]

    return(Tprime)


def avarage_tensor(gs, fracs, T):
    """Use an ODF stored in a list of rotation matrices to calculate
       the Voigt avarage of a second or fourth rank tensor (assuming 
       one phase). Arguments are:

       gs: a list of rotation matrices (list of 3*3 numpy arrays)
       fracs: a list of volume fractions (list float)
       T: 3*3*3*3 or 3*3 numpy array of the single X-tal tensor
      
    """
    if (T.ndim == 4):
        # 4th order tensor (e.g. elastic constants)
        T_av = np.zeros((3,3,3,3))
        for g, frac in zip(gs, fracs):    
            # Because we are going from an array to a list in the driver
            #routine, we need to strip off an empty dimension in g...
            # which is why the [0] hack can be seen!
            Tprime = rotT(T, g[0])
            T_av = T_av + (Tprime * frac) # * XOl... possibly

    return T_av


def calc_cij(cijs_mat, gs, fracs, scheme='Voigt'):
    """Use an ODF stored in a VPSC style _TEX file to calculate
       the Voigt avarage of the elastic constants tensor (assuming 
       one phase). Arguments are:
       cijs_mat: 6*6 numpy array of the single X-tal elastic constants matrix
       gs: n*3*3 numpy array of rotation matricies for the n crystals
       fracs: length n array of volume fractions
       scheme: optional string indicating if Voigt, Reuss or Hill (i.e. 
               Voigt-Reuss-Hill) avaraging should be used. Defaults to 
               Voigt, which is fastest and provides an upper bound 
               estimate.
      
       Code is based on the Fortran 'voigt' subroutine distributed as part of
       DRex"""

    if (scheme == 'Voigt'):
        method = 0
    elif (scheme == 'Reuss'):
        method = 1
    elif (scheme == 'Hill'):
        method = 2
    else:
        raise ValueError('Scheme argument must be one of Voigt, Reuss or Hill')
 
    # Chop up input and make lists
    gs = np.split(gs, np.size(fracs))
    fracs = np.split(fracs, np.size(fracs))

    # Do the avaraging (of cijs or sijs or both)
    if ((method == 0) or (method == 2)):
        if (cijs_mat.shape == (6,6)):
            cijs = mat2tens(cijs_mat)
        else:
            cijs = cijs_mat
        cij_voigt = avarage_tensor(gs, fracs, cijs)
        if (cijs_mat.shape == (6,6)):
            cij_voigt_mat = tens2mat(cij_voigt)
        else:
            cij_voigt_mat = cij_voigt
        if (method == 0): return cij_voigt_mat
    if ((method == 1) or (method == 2)):
        if (cijs_mat.shape == (6,6)):
            sijs = mat2tens(npl.inv(cijs_mat), compl=True)
        elif (cijs_mat.shape == (3,3)):
            sijs = npl.inv(cijs_mat)
        sij_reuss = avarage_tensor(gs, fracs, sijs)
        if (cijs_mat.shape == (6,6)):
            cij_reuss_mat = npl.inv(tens2mat(sij_reuss, compl=True))
        else:
            cij_reuss_mat = npl.inv(sij_reuss)
        if (method == 1): return cij_reuss_mat
    if (method == 2):
        cij_av_mat = (cij_voigt_mat + cij_reuss_mat) / 2.0
    elif (method == 1):
        cij_av_mat = cij_reuss_mat
    elif (method == 0):
        cij_av_mat = cij_voigt_mat

    return cij_av_mat



