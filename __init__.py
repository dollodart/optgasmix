import os.path as osp
data_dir = osp.join(osp.dirname(osp.abspath(__file__)), 'data')

lennard_jones_data_file = osp.join(data_dir, 'bsl-lj.csv')
lennard_jones_heat_capacities_data_file = osp.join(data_dir, 'bsl-koretsky-merged.csv')
collision_cross_section_data_file = osp.join(data_dir, 'bsl-colcross.csv')

heat_capacities_constant_factors = {'B':1e3, 'C':1e6, 'D': 1e-5, 'E':1e9}
