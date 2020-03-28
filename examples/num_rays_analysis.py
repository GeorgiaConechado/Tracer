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
import matplotlib.pyplot as plt

# created a function to use in the analysis
def run (num_rays):

    theta_zenith_source = 10./180.*N.pi #everything is done in rafdians in tracer.
    theta_zenith_plate = 20./180.*N.pi #everything is done in rafdians in tracer.
    theta_azimuth = N.pi/2.

    ####################
    # 1: Source creation:

    # num_rays = 1000000 # Number of rays to be traced.
    center = N.vstack([0, N.sin(-theta_zenith_plate), N.cos(theta_zenith_plate)]) # center of the source geometry. vstack makes the vector vertical. We use a bit of trigonometry to place the center of the source on the normal vector from the surface we will model.
    zenith_direction = [0,0,1] # verical direction in 3D space.
    source_direction = N.dot(zenith_direction, rotz(theta_azimuth)[:3,:3]) # We rotate around z to align with the azimuth angle of the surface
    source_direction = N.dot(source_direction, rotx(N.pi-theta_zenith_plate)[:3,:3]) # We rotate around x to align with the zenith angle of the plate, plus pi to come from that direction instead of going there. Now we have the direction of the source normal. This has to be a vector of unit length. We rotate the source to be parallel to the flat surface.
    rays_direction = N.dot(zenith_direction, rotz(theta_azimuth)[:3,:3])
    rays_direction = N.dot(rays_direction, rotx(N.pi-theta_zenith_source)[:3,:3]) # General direction of propagation of the rays, this time woith the source angle.
    x = 1.5 # Length in x direction
    y = 2.5 # Length in y direction
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

    flat_plate_surface = Surface(geometry=panel_geometry, optics=panel_optics, location=N.array([0,0,0]), rotation=plate_rotation)

    flat_plate_object = AssembledObject(surfs=[flat_plate_surface])

    scene = Assembly(objects=[flat_plate_object])
    ####################

    ####################
    # 3: Build tracer assemply
    engine = TracerEngine(scene)
    ####################

    ####################
    # 4: Ray-trace:
    engine.ray_tracer(bundle=source) # bundle is a ray-bundle object from tracer. Source functions return ray-bundle objects.
    #Renderer(engine).show_rays() # If you want to skip the rendering, just comment this line by adding a '#' charater at the start.
    ####################

    ####################
    # 5: Analyse results:
    absorbed_energy, hit_locations = panel_optics.get_all_hits()
    #print hit_locations.shape
    hit_locations = flat_plate_surface.global_to_local(hit_locations) # This is VERY usefu. you can use the global_to_local from a surface to transpose the cooridinates in the global referential to the ones in the local coordinate system of the surface. It makes binning much easier.


    
    from sys import path # this is something to find your local directory without having to type the whole address

    bins, xedges, yedges = N.histogram2d(hit_locations[0,:], y=hit_locations[1,:], bins=(6,12), weights=absorbed_energy)

    ####################

    flat_bins = [item for sublist in bins for item in sublist] #single array of values of absorbed power

    difference = max(flat_bins)-min(flat_bins) #finds the difference in values of absorbed power (W)
    return difference


num_rays_list = [1000, 10000, 100000, 1000000, 5000000, 10000000] #different values of num_rays to analyse 

ranges = map(run, num_rays_list)

#there is a lot of randomness so i want to find the average difference for each value of num_rays
def findAvg (num_rays) : 
    sum = 0
    n = 10
    for i in range(n) :
       sum = sum + run(num_rays)
    return sum/n


#avg_ranges = map(findAvg, num_rays_list)


from sys import path

print ranges
plt.plot(num_rays_list, ranges)
plt.xscale('log')
plt.grid(True)
plt.title('Range of values vs Number of rays')
plt.xlabel('Number of rays')
plt.ylabel('Range of values of absorbed power (W)')
#plt.show()

plt.savefig(path[0]+'/num_rays_analysis.png')