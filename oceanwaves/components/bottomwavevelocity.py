"""Calculate maximum bottom wave velocity and period for wind waves.

Examples
--------

Create a file that contains a time-series of wave height, period [or wave 
spectrum] and water depth.

>>> data = np.array([
...     [0., 1., 10., 10.],
...     [2., 2., 10., 10.],
...     [4., 3., 10., 10.],
... ])
>>> np.savetxt('wavetest.txt', data)

Create an instance of BottomWaveVelocity and initialize it with the input file.

>>> ubr = BottomWaveVelocity()
>>> ubr.initialize('wavetest.txt')
>>> ubr.get_start_time()
0.0
>>> ubr.get_end_time()
4.0

Advance the component and get bottom wave velocity and period.

>>> ubr.update()
>>> ubr.get_current_time()
1.0
>>> btmorbvel = ubr.get_grid_values('sea_bottom_water_wave__max_of_orbital_speed')
>>> btmperiod = ubr.get_grid_values('sea_bottom_water_wave__period')
>>> print round(btmorbvel, 6)
0.079255
>>> print round(btmperiod, 6)
1.330892
"""
from bisect import bisect
import numpy as np

from oceanwaves import ubspecpar


class BottomWaveVelocity(object):
    _name = 'BottomWaveVelocity'

    _input_var_names = [
        'sea_surface_water_wave__height',
        'sea_surface_water_wave__period',
        'sea_water__depth',
    ]

    _output_var_names = [
        'sea_bottom_water_wave__max_of_orbital_speed',
        'sea_bottom_water_wave__period',
    ]

    def __init__(self):
        self._water_depth = 0.
        self._wave_height = 0.
        self._wave_period = 0.
        self._wave_btmorbvel = 0.
        self._wave_btmperiod = 0.
        self._time = 0.
        self._time_step = 1.
        self._time_units = 's'
        self._data = np.zeros((1, 4))

    def get_component_name(self):
        return self._name

    def get_input_var_names(self):
        return self._input_var_names

    def get_output_var_names(self):
        return self._output_var_names

    def initialize(self, filename):
        self._data = np.genfromtxt(filename)
        self._wave_btmorbvel, self._wave_btmperiod = self.calculate_ubr_vars(0.)

    def calculate_ubr_vars(self, time):
        index = bisect(self._data[:, 0], self._time)
        (self._wave_height,
         self._wave_period,
         self._water_depth) = self._data[index, 1:]

        (wave_btmorbvel,
         wave_btmperiod) = ubspecpar(wind_speed, water_depth, fetch)

        return wave_btmorbvel, wave_btmperiod

    def update(self):
        self._time = self._time + self._time_step

        (self._wave_btmorbvel,
         self._wave_btmperiod) = self.calculate_ubr_vars(self._time)

    def update_until(self, time):
        self._time = time

        (self._wave_btmorbvel,
         self._wave_btmperiod) = self.calculate_ubr_vars(self._time)

    def finalize(self):
        self.__init__()

    def get_grid_values(self, name):
        if name == 'sea_bottom_water_wave__max_of_orbital_speed':
            return self._wave_btmorbvel
        elif name == 'sea_bottom_water_wave__period':
            return self._wave_btmperiod
        else:
            raise TypeError('name not understood')

    def set_grid_values(self, name, value):
        if name == 'sea_surface_water_wave__height':
            self._wave_height = value
        elif name == 'sea_surface_water_wave__period':
            self._wave_period = value
        elif name == 'sea_water__depth':
            self._water_depth = value
        else:
            raise TypeError('name not understood')

    def get_start_time(self):
        return self._data[0, 0]

    def get_current_time(self):
        return self._time

    def get_end_time(self):
        return self._data[-1, 0]

    def get_time_step(self):
        return self._time_step

    def get_time_units(self):
        return self._time_units
