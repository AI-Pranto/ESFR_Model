import openmc
import openmc.mgxs as mgxs
import numpy as np

#################################
#          MATERIALS
################################
fuel = openmc.Material(name= 'core fuel', temperature = 900.0)
fuel.add_nuclide('U234', 1.6447E-06)
fuel.add_nuclide('U235', 2.2148E-05)
fuel.add_nuclide('U236', 2.7420E-06)
fuel.add_nuclide('U238', 1.5602E-02)
fuel.add_nuclide('Np237', 2.8455E-05)
fuel.add_nuclide('Pu236', 2.6787E-10)
fuel.add_nuclide('Pu238', 1.0478E-04)
fuel.add_nuclide('Pu239', 1.7966E-03)
fuel.add_nuclide('Pu240', 1.2659E-03)
fuel.add_nuclide('Pu241', 2.0302E-04)
fuel.add_nuclide('Pu242', 2.7492E-04)
fuel.add_nuclide('Am241', 1.0157E-04)
fuel.add_nuclide('Am242_m1', 8.5597E-06)
fuel.add_nuclide('Am243', 9.0780E-05)
fuel.add_nuclide('Cm242', 5.2317E-06)
fuel.add_nuclide('Cm243', 6.2354E-07)
fuel.add_nuclide('Cm244', 6.8704E-05)
fuel.add_nuclide('Cm245', 2.0198E-05)
fuel.add_nuclide('Cm246', 1.2227E-05)
fuel.add_nuclide('O16',   4.1265E-02)
fuel.add_element('Mo', 9.9952E-04)

ht = openmc.Material(name= 'HT-9 clad mat', temperature = 600.0)
ht.add_element('Fe', 6.9715E-02)
ht.add_element('Ni', 4.2984E-04)
ht.add_element('Cr', 1.0366E-02)
ht.add_nuclide('Mn55', 4.5921E-04)
ht.add_element('Mo', 4.9007E-04)

coolant = openmc.Material(name= 'coolant Na mat', temperature = 600.0)
coolant.add_element('Na', 2.2272E-02)

materials = openmc.Materials([fuel, ht, coolant])

#######################################
#        GEOMETRY
#######################################
# Dimensions
fuel_radius       = 0.3322
clad_outer_radius = 0.3928
pin_pitch         = 0.8966

fuel_inner = openmc.ZCylinder(r = fuel_radius)
clad_or    = openmc.ZCylinder(r = clad_outer_radius)
pin_universe = openmc.model.pin(
    [fuel_inner, clad_or],
    [fuel, ht, coolant]
)

na_cell = openmc.Cell(fill= coolant)
na_uni = openmc.Universe(cells=(na_cell,))

lattice = openmc.HexLattice()
lattice.center = (0., 0.)
lattice.pitch = (pin_pitch,)
lattice.orientation = 'x'
lattice.universes = [
    [pin_universe for _ in range(max(1, 6*ring_index))]
    for ring_index in reversed(range(10))
]
lattice.outer = na_uni

subassembly_duct_outer     = 15.8123
subassembly_duct_thickness = 0.3966
subassembly_duct_inner     = subassembly_duct_outer - 2*subassembly_duct_thickness
subassembly_pitch          = 16.2471

duct_inner_hex = openmc.model.hexagonal_prism(edge_length= subassembly_duct_inner / np.sqrt(3.), orientation='x')
duct_outer_hex = openmc.model.hexagonal_prism(edge_length = subassembly_duct_outer / np.sqrt(3.), orientation='x')
outer_hex = openmc.model.hexagonal_prism(edge_length=subassembly_pitch / np.sqrt(3.), orientation='x', boundary_type='periodic')

lattice_cell  = openmc.Cell(fill=lattice, region= duct_inner_hex)
duct_ht       = openmc.Cell(fill=ht, region=~duct_inner_hex & duct_outer_hex)
duct_outside  = openmc.Cell(fill=coolant, region= ~duct_outer_hex & outer_hex)
fuel_universe = openmc.Universe(cells= [lattice_cell, duct_ht, duct_outside])

geometry = openmc.Geometry(root=fuel_universe)

#############################
#        SETTINGS
############################
lower_left  = [-11.98, -11.98, -1.0]
upper_right = [11.98,  11.98,   1.0]
uniform_dist = openmc.stats.Box(lower_left, upper_right, only_fissionable=True)
src = openmc.Source(space=uniform_dist)
settings = openmc.Settings()
settings.source = src
settings.batches = 500
settings.inactive = 100
settings.particles = 24000
settings.output = {'summary': True}

#############################################
#      Group constant generation
#############################################
# 24-group
groups = openmc.mgxs.EnergyGroups(np.array([1.00001E-11, 1.48625E-04, 3.04325E-04, 4.53999E-04, 7.48518E-04, 1.23410E-03, 2.03468E-03, 3.35463E-03,
                    5.53084E-03, 9.11882E-03, 1.50344E-02, 2.47875E-02, 4.08677E-02, 6.73795E-02, 1.11090E-01,
                    1.83156E-01, 3.01974E-01, 4.97871E-01, 8.20850E-01, 1.35335E+00, 2.23130E+00, 3.67879E+00,
                    6.06531E+00, 1.00000E+01, 1.96403E+01])*1e6)

# 4 group
# groups = openmc.mgxs.EnergyGroups(np.array([1e-5, 2.67030e+1, 1.72230e+3, 1.11090e+5, 1.41910e+07])*1e6)

total        = mgxs.TotalXS(domain=fuel_universe, energy_groups=groups, by_nuclide=False)
transport    = mgxs.TransportXS(domain=fuel_universe, energy_groups=groups, by_nuclide=False)
absorption   = mgxs.AbsorptionXS(domain=fuel_universe, energy_groups=groups, by_nuclide=False)
scattering_m = mgxs.ScatterMatrixXS(domain=fuel_universe, energy_groups=groups, by_nuclide=False)
scattering   = mgxs.ScatterXS(domain=fuel_universe, energy_groups=groups, by_nuclide=False)
fission      = mgxs.FissionXS(domain=fuel_universe, energy_groups=groups, by_nuclide=False, nu=True)
chi          = mgxs.Chi(domain=fuel_universe, energy_groups=groups, by_nuclide=False)

# Instantiate an empty Tallies object
tallies  = openmc.Tallies()
tallies += total.tallies.values()
tallies += transport.tallies.values()
tallies += absorption.tallies.values()
tallies += scattering_m.tallies.values()
tallies += scattering.tallies.values()
tallies += fission.tallies.values()
tallies += chi.tallies.values()

model = openmc.model.Model(geometry, materials, settings, tallies)
model.run()
