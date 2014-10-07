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
>>> height = waves.get_grid_values('sea_surface_water_wave__height')
>>> period = waves.get_grid_values('sea_surface_water_wave__period')
>>> spectrum = waves.get_grid_values('sea_surface_water_wave__spectrum')
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

        (wave_height, wave_period, wave_spectrum, wave_specfrequencies) = ow.jonswap_hs(wind_speed, water_depth, fetch)

        return wave_height, wave_period, wave_spectrum, wave_specfrequencies

    def update(self):
        self._time = self._time + 1

        (self._wave_height,
         self._wave_period,
         self._wave_spectrum) = self.calculate_wave_vars(self._time)

    def update_until(self, time):
        self._time = time

        (self._wave_height,
         self._wave_period,
         self._wave_spectrum) = self.calculate_wave_vars(self._time)

    def finalize(self):
        self.__init__()

    def get_grid_values(self, name):
        if name == 'sea_surface_water_wave__height':
            return self._wave_height
        elif name == 'sea_surface_water_wave__period':
            return self._wave_period
        elif name == 'sea_surface_water_wave__spectrum':
            return self._wave_spectrum
        elif name == 'sea_surface_water_wave__specfrequencies':
            return self._wave_specfrequencies
        else:
            raise TypeError('name not understood')

    def set_grid_values(self, name, value):
        if name == 'land_surface_10m-above_air_flow__speed':
            self._wind_speed = value
        elif name == 'sea_surface_air_flow__fetch_length':
            self._fetch = value
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
