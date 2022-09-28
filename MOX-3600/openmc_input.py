import openmc
import openmc.mgxs as mgxs
import numpy as np

#########################################
#         MATERIALS
#########################################
fuel = openmc.Material(name='core fuel', temperature=900.0)
fuel.add_element('O', 4.2825e-02)
fuel.add_nuclide('U234', 1.7602e-06)
fuel.add_nuclide('U235', 3.3412e-05)
fuel.add_nuclide('U236', 4.0736e-06)
fuel.add_nuclide('U238', 1.8692e-02)
fuel.add_nuclide('Np237', 3.7863e-06)
fuel.add_nuclide('Np239', 3.5878e-06)
fuel.add_nuclide('Pu238', 9.4366e-05)
fuel.add_nuclide('Pu239', 1.8178e-03)
fuel.add_nuclide('Pu240', 1.0177e-03)
fuel.add_nuclide('Pu241', 2.1797e-04)
fuel.add_nuclide('Pu242', 3.2651e-04)
fuel.add_nuclide('Am241', 4.0395e-05)
fuel.add_nuclide('Am242', 1.2900e-08)
fuel.add_nuclide('Am242_m1', 1.2243e-06)
fuel.add_nuclide('Am243', 2.4048e-05)
fuel.add_nuclide('Cm242', 2.2643e-06)
fuel.add_nuclide('Cm243', 1.0596e-07)
fuel.add_nuclide('Cm244', 3.2454e-06)
fuel.add_nuclide('Cm245', 1.4745e-07)
fuel.add_nuclide('Cm246', 3.4906e-09)
fuel.add_element('Mo', 2.7413e-03)

em10 = openmc.Material(name='EM10', temperature=600.0)
em10.add_element('C', 3.8254e-04)
em10.add_element('Si', 4.9089e-04)
em10.add_element('Ti', 1.9203e-05)
em10.add_element('Cr', 7.5122e-03)
em10.add_element('Fe', 7.3230e-02)
em10.add_element('Ni', 3.9162e-04)
em10.add_element('Mo', 4.7925e-04)
em10.add_element('Mn', 4.1817e-04)

ods = openmc.Material(name='ODS', temperature=600.0)
ods.add_element('C', 3.5740e-04)
ods.add_element('O', 3.9924e-04)
ods.add_element('Ti', 5.3824e-04)
ods.add_element('Cr', 1.7753e-02)
ods.add_element('Fe', 5.3872e-02)
ods.add_element('Ni', 3.6588e-04)
ods.add_element('Mn', 2.3441e-04)
ods.add_element('P', 2.7718e-05)
ods.add_element('Al', 9.1482e-03)
ods.add_element('Co', 2.1852e-04)
ods.add_element('Cu', 1.0135e-04)
ods.add_element('Y', 2.6616e-04)

coolant = openmc.Material(name= 'coolant Na mat', temperature = 600.0)
coolant.add_element('Na', 2.1924e-02)

helium = openmc.Material(name='helium mat', temperature = 600.0)
helium.add_element('He', 1.0e-6)

materials = openmc.Materials([fuel, coolant, em10, ods, helium])

##########################
#     GEOMETRY
##########################
# Dimensions
inner_hole_radius = 0.1257
fuel_radius       = 0.4742
clad_inner_radius = 0.4893
clad_outer_radius = 0.5419
pin_pitch         = 1.1897

fuel_inner = openmc.ZCylinder(r= inner_hole_radius)
fuel_outer = openmc.ZCylinder(r= fuel_radius)
clad_inner = openmc.ZCylinder(r= clad_inner_radius)
clad_outer = openmc.ZCylinder(r= clad_outer_radius)
pin_universe = openmc.model.pin(
    [fuel_inner, fuel_outer, clad_inner, clad_outer],
    [helium, fuel, helium, ods, coolant]
)

na_cell = openmc.Cell(fill= coolant)
na_uni = openmc.Universe(cells=(na_cell,))

lattice = openmc.HexLattice()
lattice.center = (0., 0.)
lattice.pitch = (pin_pitch,)
lattice.orientation = 'x'
lattice.universes = [[pin_universe]]
lattice.universes = [
    [pin_universe for _ in range(max(1, 6*ring_index))]
    for ring_index in reversed(range(10))
]
lattice.outer = na_uni

subassembly_duct_outer     = 20.7468
subassembly_duct_thickness = 0.4525
subassembly_duct_inner     = subassembly_duct_outer - 2*subassembly_duct_thickness
subassembly_pitch          = 21.2205

duct_inner_hex = openmc.model.hexagonal_prism(edge_length= subassembly_duct_inner / np.sqrt(3.), orientation='x')
duct_outer_hex = openmc.model.hexagonal_prism(edge_length = subassembly_duct_outer / np.sqrt(3.), orientation='x')
outer_hex = openmc.model.hexagonal_prism(edge_length=subassembly_pitch / np.sqrt(3.), orientation='x', boundary_type='periodic')

lattice_cell = openmc.Cell(fill=lattice, region= duct_inner_hex)
duct_em      = openmc.Cell(fill=em10, region=~duct_inner_hex & duct_outer_hex)
duct_outside = openmc.Cell(fill=coolant, region= ~duct_outer_hex & outer_hex)
fuel_universe = openmc.Universe(cells= [lattice_cell, duct_em, duct_outside])

geometry = openmc.Geometry(root=fuel_universe)

#############################
#     SETTINGS
#############################
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
