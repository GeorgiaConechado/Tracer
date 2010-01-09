# -*- coding: utf-8 -*-
# An assembly modeling the parabolic dish for building-integrated installations,
# as developped in Tel Aviv University's Faculty of Engineering.
#
# References:
# [1] Kribus A., et al, A miniature concentrating photovoltaic and thermal 
#     system, Energy Conversion and Management, Volume 47, Issue 20, December
#     2006, Pages 3582-3590, DOI: 10.1016/j.enconman.2006.01.013.
# [2] Harald Ries et al., High-flux photovoltaic solar concentrators with 
#     kaleidoscope-based optical designs, 1997, Solar Energy, 
#     DOI: 10.1016/S0038-092X(96)00159-4

from .. import optics_callables as opt
from ..surface import Surface
from ..paraboloid import ParabolicDishGM
from .homogenized_local_receiver import HomogenizedLocalReceiver

class MiniDish(HomogenizedLocalReceiver):
    def __init__(self, diameter, focal_length, dish_opt_eff,\
        receiver_pos, receiver_side, homogenizer_depth, homog_opt_eff):
        """
        Arguments:
        diameter, focal_length - of the parabolic dish
        dish_opt_eff - the optical efficiency of the dish
        receiver_pos - the distance along the optical axis from the dish to the
            receiver's end surface - the PV panel (should be about the focal 
            length)
        receiver_side - the receiver is square, with this side length.
        homogenizer_depth - the homogenizer has base dimensions to fit the PV
            square, and this height.
        homog_opt_eff - the optical efficiency of each mirror in the homogenizer
        """
        dish_surf = Surface(ParabolicDishGM(diameter, focal_length), 
            opt.Reflective(1 - dish_opt_eff))
        HomogenizedLocalReceiver.__init__(self, dish_surf, receiver_pos, \
            receiver_side, homogenizer_depth, homog_opt_eff)
        
        # for later interrogation:
        self._ext_dims = (diameter, receiver_pos)
        
    def get_external_dimensions(self):
        """
        Returns the external dimensions (the ones you would put in an assembly
        drawing, the bounding cylinder dimensions) of the entire assembly.
        
        Returns:
        diameter - of the dish
        full_height - from the dish base to the receiver surface.
        """
        return self._ext_dims

# Utility functions:
def standard_minidish(diameter, concentration, reflections, 
    dish_opt_eff=0.9, homog_opt_eff=0.9):
    """
    Create a minidish assembly with dimensions based on a dish with 45 deg.
    acceptance angle in the receiver, using dimensioning rules from [2].
    
    Arguments:
    diameter - of the dish aperture.
    concentrations - ratio of dish aperture to receiver aperture.
    reflections - number of times an edge ray gets reflected in the homogenizer
    dish_opt_eff, homog_opt_eff - passed directly to the minidish constructor.
    
    Returns:
    minidish - a MiniDish instance with the correct sizing of components.
    f, W, H - the focal length, homogenizer width and receiver distance from
        focal point that were used for the dish.
    """
    f = (N.sqrt(5) - 1)*diameter/4.
    W = diameter/2. * N.sqrt(N.pi/concentration)
    n = reflections + 1
    H = n*W*f/(diameter - n*W)
    minidish = MiniDish(diameter, f, dish_opt_eff, f + H, W, H, homog_opt_eff)
    return minidish, f, W, H
