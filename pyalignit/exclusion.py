# Copyright (C) 2021 OliverBScott
"""
pyalignit.exclusion
-------------------
Utilities for creating exclusion spheres

"""
import numpy as np

from .cpyalignit import PharmacophorePoint, FuncGroup


# Use 0.7 as default exclusion sphere alpha
DEFAULT_EXCL_ALPHA = 0.7


def create_exclusion_spheres(
    pharmacophore, receptor, ligand=None, cutoff=4.5,
    receptor_confid=-1, ligand_confid=-1
):
    """Create and add exclusion spheres to a pharmacophore.

    Exclusion spheres are added to the Pharmacophore object,
    created from each receptor atom. If a ligand is also
    provided exclusion spheres are only built from receptor 
    atoms within a distance threshold from any ligand atom
    (default = 4.5A). The provided pharmacophore is modified
    inplace. 
    
    Parameters
    ----------
    pharmacophore : Pharmacophore
        The pharmacophore to add exclusion spheres to.
    receptor : RDKit Mol
        An RDKit Mol object from which exclusion spheres
        will be built (the receptor).
    ligand : RDKit Mol, optional
        An RDKit Mol containing the ligand the pharmacophore
        represents. When provided only the receptor atoms within
        a distance threshold, `cutoff`, from any ligand atom are
        used to build exclusion spheres.
    cutoff : float, default=4.5
        Distance cutoff, used for distance based receptor atom 
        filtering when a ligand is also provided.
    receptor_confid : int, default=-1
        Conformer ID to use for the provided receptor.
    ligand_confid : int, default=-1
        Conformer ID to use for the ligand if provided.

    Notes
    -----
    Exclusion spheres have a different role to other pharmacophore
    points during alignment, indicating regions in the pharmacophore
    model where no pharmacophore points are allowed during alignment.
    The use of exclusion spheres in a pharmacophore model nicely 
    mimics the spatial constraints of an active site. Exclusion spheres
    have a default alpha of 0.7.

    """
    # get positions of receptor atoms
    rconf = receptor.GetConformer(receptor_confid)
    rpos = rconf.GetPositions()

    # filter receptor positions using distance cutoff 
    if ligand is not None:
        lconf = ligand.GetConformer(ligand_confid)
        lpos = lconf.GetPositions()
        dists = _distsq_nxm(lpos, rpos)
        uix = np.unique(np.where(dists < cutoff)[1])
        rpos = rpos[uix]
    
    # add pharmacophore at each receptor position
    for pos in rpos:
        p = PharmacophorePoint()
        p.func = FuncGroup.EXCL
        p.alpha = DEFAULT_EXCL_ALPHA
        p.point.x = pos[0]
        p.point.y = pos[1]
        p.point.z = pos[2]
        pharmacophore.append(p)


def _distsq_nxm(X, Y):
    sx = np.sum(X**2, axis=1, keepdims=True)
    sy = np.sum(Y**2, axis=1, keepdims=True)
    return np.sqrt(-2 * X.dot(Y.T) + sx +sy.T)
