from capture_check import plot_neut_pos
import pandas as pd
from matplotlib import pyplot
from amuse.plot import *
from amuse.lab import units


lines = plot_neut_pos()
lines.get_dataframe_from_pkl(dat_str='neut_stars_positions.pkl')
lines.plot(use_all=True, fix='y')
