# If you are not familiar with Python, this is where we load all the modules that we will need in this code.
import numpy as N # Numpy is basically the matlab equivalent in Python (but better)
from tracer.sources import oblique_solar_rect_bundle # a source model
from tracer.assembly import Assembly # the Assembly class that we use to regroup all the things in our model and give to the ray-tracing engine
from tracer.object import AssembledObject # The AssembledObject class that is used to create objects that are composed of surfaces.
from tracer.surface import Surface # The surface is what holds the informatipon of the geometry and the optical properties of the things we want to model
from tracer.flat_surface import RectPlateGM # The rectangular plate geometry manager
from tracer.optics_callables import LambertianReceiver, LambertianReflector # the opticla properties manager that can manage diffuse reflections and keep track of the loctaion of the hist and the energy absorbed.
from tracer.spatial_geometry import * # These are useful functiosn to generate rotation matrices etc. The '*' means that we import everything
from tracer.tracer_engine import TracerEngine # The ray-tracing engine that coordinates the ray-tracing process
from tracer.CoIn_rendering.rendering import Renderer

'''
This simple script will generate rays from a rectangular source towards a rectangular object on which we will gather data.

The type of source chosen here is able to have an orientation of the source plane different to the general propagation direction of the rays; which is convenient for problmes where the surfaces are not directly facing the sun and for problems where the angular range of the incoming rays can vary widely from the normal (like with diffuse sources etc.

The basic logic is the following:
1 - We create a source
2 - We create an assembly which is a collection of assembled objects that have a geometry and optical properties
3 - We create a tracer engine that has the assembly as an argument
4 - We use the ray-tracer method form the tracer engine with the source as an argument to perfrom the ray-trace.
5 - We collect the results and analyse them

In this example, we will model how DNI coming from the sun at a position of 10 degrees zenith angle and 90 degrees azimuth angle is irradiating a 1 m by 2 m rectangular flat plate at 20 degrees zenith angle and the azimuthal same orientation.
'''
theta_zenith_source = 10./180.*N.pi #everything is done in rafdians in tracer.
theta_zenith_plate = 20./180.*N.pi #everything is done in rafdians in tracer.
theta_azimuth = N.pi/2.

####################
# 1: Source creation:

num_rays = 1000000 # Number of rays to be traced.
center = N.vstack([0, N.sin(-theta_zenith_plate), 3*N.cos(theta_zenith_plate)]) # center of the source geometry. vstack makes the vector vertical. We use a bit of trigonometry to place the center of the source on the normal vector from the surface we will model.
zenith_direction = [0,0,1] # verical direction in 3D space.
source_direction = N.dot(zenith_direction, rotz(theta_azimuth)[:3,:3]) # We rotate around z to align with the azimuth angle of the surface
source_direction = N.dot(source_direction, rotx(N.pi-theta_zenith_plate)[:3,:3]) # We rotate around x to align with the zenith angle of the plate, plus pi to come from that direction instead of going there. Now we have the direction of the source normal. This has to be a vector of unit length. We rotate the source to be parallel to the flat surface.
rays_direction = N.dot(zenith_direction, rotz(theta_azimuth)[:3,:3])
rays_direction = N.dot(rays_direction, rotx(N.pi-theta_zenith_source)[:3,:3]) # General direction of propagation of the rays, this time woith the source angle.
x = 4 # Length in x direction
y = 4 # Length in y direction
ang_range = 4.65e-3 #(rad) this is the typical angular range of the solar disk. This way of modelling the DNI is called the "Pillbox" unshape model.
flux = 1000. # (W/m2) the solar irradiance modelled by the source. Here we need to be careful about aaht we model; namely DNI or GHI. If we model DNI, it sis only the Direct Normal Irradiance fraction.

# source declaration:
source = oblique_solar_rect_bundle(num_rays, center, source_direction, rays_direction, x, y, ang_range, flux)
####################

####################
# 2: Assembly creation:
panel_absorptivity = 0.9 # absorptivity of the surface

panel_width = 1. # width of the rectangle plate in m
panel_height = 2. # height of the rectangle plate in m

panel_geometry = RectPlateGM(panel_width, panel_height)
panel_optics = LambertianReceiver(panel_absorptivity)

plate_rotation = N.dot(rotx(theta_zenith_plate)[:3,:3], rotz(theta_azimuth)[:3,:3]) # Here we directly do the composed rotation matrix to have it as an argument in the following instantiation of the Surface class.

flat_plate_surface = Surface(geometry=panel_geometry, optics=panel_optics, location=N.array([0,0,1]), rotation=plate_rotation)

flat_plate_object = AssembledObject(surfs=[flat_plate_surface])

#ground

ground_absorptivity = 0.8 # absorptivity of the surface

ground_width = 4. # width of the rectangle plate in m
ground_height = 4. # height of the rectangle plate in m

ground_geometry = RectPlateGM(ground_width, ground_height)
ground_optics = LambertianReflector(ground_absorptivity)

ground_rotation = N.dot(rotx(0)[:3,:3], rotz(0)[:3,:3])

ground_surface = Surface(geometry=ground_geometry, optics=ground_optics, location=N.array([0,0,0]), rotation = ground_rotation)

ground_object = AssembledObject(surfs=[ground_surface]) # replaced "flat_plate" by "ground"

scene = Assembly(objects=[flat_plate_object,ground_object])
####################

####################
# 3: Build tracer assemply
engine = TracerEngine(scene)
#Renderer(engine).show_geom() # This rendering function checks teh geometry only prior to ray tracing.
####################

####################
# 4: Ray-trace:
engine.ray_tracer(bundle=source) # bundle is a ray-bundle object from tracer. Source functions return ray-bundle objects.
Renderer(engine).show_rays() # If you want to skip the rendering, just comment this line by adding a '#' charater at the start.
####################

####################
# 5: Analyse results:
absorbed_energy, hit_locations = panel_optics.get_all_hits()
print hit_locations.shape
hit_locations = flat_plate_surface.global_to_local(hit_locations) # This is VERY usefu. you can use the global_to_local from a surface to transpose the cooridinates in the global referential to the ones in the local coordinate system of the surface. It makes binning much easier.


import matplotlib.pyplot as plt # The library to plot data
from sys import path # this is something to find your local directory without having to type the whole address

plt.figure()
plt.subplot(111, aspect='equal')
bins, xedges, yedges = N.histogram2d(hit_locations[0,:], y=hit_locations[1,:], bins=(6,12), weights=absorbed_energy)
X, Y = N.meshgrid(xedges, yedges) # This meshgrid is a bit annoying but necessary to plot the 2d fluxmap. Mor eon the online documentation of pyplot.
plt.pcolormesh(X, Y, bins.T) # the .T is the transpose operation. Some functions do some flipping of the axes of teh arrays due to how the code is processing the data.
plt.xlabel('x (m)')
plt.xlabel('y (m)')

plt.colorbar(label='Absorbed power (W)')
plt.savefig(path[0]+'/absorbed_power.png')
####################


