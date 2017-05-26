"""Calculate wave heights and periods for ocean storms.

Examples
--------

Create a file that contains a time-series for wind speeds.

>>> data = np.array([
...     [0., 1., 50., 10000.],
...     [2., 2., 50., 10000.],
...     [4., 3., 50., 10000.],
... ])
>>> np.savetxt('windtest.txt', data)

Create an instance of WindWaves and initialize it with the input file.

>>> waves = WindWaves()
>>> waves.initialize('windtest.txt')
>>> waves.get_start_time()
0.0
>>> waves.get_end_time()
4.0

Advance the component and get wave height and peak period.

>>> waves.update()
>>> waves.get_current_time()
1.0
>>> height = waves.get_value('sea_surface_water_wave__height')
>>> period = waves.get_value('sea_surface_water_wave__period')
>>> spectrum = waves.get_value('sea_surface_water_wave__spectrum')
>>> print round(height, 6)
0.079255
>>> print round(period, 6)
1.330892
"""
from bisect import bisect
import numpy as np
import oceanwaves as ow


class WindWaves(object):
    _input_var_names = [
        'land_surface_10m-above_air_flow__speed',
        'sea_surface_air_flow__fetch_length',
        'sea_water__depth',
    ]

    _output_var_names = [
        'sea_surface_water_wave__height',
        'sea_surface_water_wave__period',
        'sea_surface_water_wave__spectrum',
    ]

    _var_units = {
        'land_surface_10m-above_air_flow__speed': 'm / s',
        'sea_surface_air_flow__fetch_length': 'm',
        'sea_water__depth': 'm',
        'sea_surface_water_wave__height': 'm',
        'sea_surface_water_wave__period': 's',
        'sea_surface_water_wave__spectrum': '-',
    }

    _values = {
        'sea_surface_water_wave__height': '_wave_height',
        'sea_surface_water_wave__period': '_wave_period',
        'sea_surface_water_wave__spectrum': '_wave_spectrum',
    }

    def __init__(self):
        self._wind_speed = 0.
        self._water_depth = 0.
        self._fetch = 0.
        self._wave_height = 0.
        self._wave_period = 0.
        self._wave_spectrum = []
        self._wave_specfrequencies = []
        self._time = 0.
        self._data = np.zeros((1, 4))

    def get_input_var_names(self):
        return self._input_var_names

    def get_output_var_names(self):
        return self._output_var_names

    def initialize(self, filename):
        self._data = np.genfromtxt(filename)
        self._wave_height, self._wave_period = self.calculate_wave_vars(0.)

    def calculate_wave_vars(self, time):
        index = bisect(self._data[:, 0], self._time)
        (wind_speed, water_depth, fetch) = self._data[index, 1:]

        wave_height, wave_period = ow.jonswap_hs(wind_speed, fetch)

        return wave_height, wave_period

    def update(self):
        self._time = self._time + 1

        (self._wave_height,
         self._wave_period) = self.calculate_wave_vars(self._time)

    def update_until(self, time):
        self._time = time

        (self._wave_height,
         self._wave_period) = self.calculate_wave_vars(self._time)

    def finalize(self):
        self.__init__()

    def get_var_units(self, name):
        return self._var_units[name]

    def get_var_grid(self, name):
        return 0

    def get_grid_type(self, gid):
        return 'point'

    def get_grid_x(self, gid):
        return 0.

    def get_grid_y(self, gid):
        return 0.

    def get_value(self, name):
        return getattr(self, self._values[name])

    def get_value_ref(self, name):
        raise NotImplementedError('get_value_ref')

    def set_value(self, name, value):
        if name in self._input_var_names:
            setattr(self, self._values[name], value)
        else:
            raise ValueError('{name}: not an input item'.format(name=name))

    def get_start_time(self):
        return self._data[0, 0]

    def get_current_time(self):
        return self._time

    def get_end_time(self):
        return self._data[-1, 0]
