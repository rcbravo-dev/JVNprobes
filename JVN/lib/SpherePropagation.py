import numpy as np
from typing import Any
import matplotlib.pyplot as plt
class Galaxy:
    '''One probe can explore one star.'''
    c = 299792.458 # km/s
    sec_per_year = 60 * 60 * 24 * 365.25 # seconds in a year
    km_in_ly = sec_per_year * c # kilometers in a light year

    def __init__(self, galaxy_radius: int = 50000, star_count: int = 100, build_time=100.0, verbose: bool = False):
        self.galaxy_radius = galaxy_radius  # light years
        self.galaxy_volume = self._volume(galaxy_radius)    # cubic light years
        self.star_count = star_count *  1e9     # billions of stars, range is 100-400
        self.build_time = build_time # time to build one probe in years
        self.verbose = verbose

        self.step = 0
        self.step_time = 10.0  # years
        self.time_elapsed = 0.0
        self.stars_to_explore = star_count * 1e9     # billions of stars, range is 100-400

        self.sphere_radius = 1.0
        self.sphere_volume = self._volume(self.sphere_radius)
        self.sphere_stars = 1
        self.sphere_probes = 1

        self.shell_volume = self.sphere_volume
        self.shell_stars = 1
        self.shell_probes = 1

    def __setattr__(self, name: str, value: Any) -> None:
        self.__dict__[name] = value

    def __call__(self, speed = 0.2):
        # Plan for the next expansion
        shell_thickness, shell_stars, volume = self.calculate_next_shell_params()
        build_time = self.calculate_time_to_build_probes(shell_stars, self.shell_probes)
        travel_time = shell_thickness / speed

        # Build the probes
        self.time_elapsed += build_time

        # Launch the probes
        self.time_elapsed += travel_time

        # Update the sphere
        self.sphere_radius += shell_thickness
        self.sphere_volume += volume
        self.sphere_stars += shell_stars
        self.sphere_probes += shell_stars

        # Update the shell
        self.shell_stars = shell_stars
        self.shell_probes = shell_stars
        self.shell_volume = volume
        self.step_time = build_time + travel_time

        # Misc updates
        self.stars_to_explore -= shell_stars
        self.step += 1

    def calculate_next_shell_params(self, dt = 0.1):
        radius = self.sphere_radius
        shell_stars = 0

        # I want the shell to contain at least one star per probe, so i gradully increase the shell 
        # size (radius) until the shell contains more (or same) stars than probes.
        while self.shell_probes >= shell_stars:
            radius += dt
            shell_volume = self._volume(radius) - self.sphere_volume
            shell_stars = int(self._stars_per_volume(shell_volume))

            if shell_stars >= self.stars_to_explore:
                shell_stars = self.stars_to_explore
                break

        shell_thickness = radius - self.sphere_radius

        return shell_thickness, shell_stars, shell_volume
    
    def calculate_time_to_build_probes(self, shell_stars: int, probes: int) -> float:
        '''Each probe takes build_time years to build a daughter probe. If more probes are available
        than stars in the shell, the extra probes can assist.'''
        steps = 0.0

        while True:
            # Build probes, double each step and check if we have enough probes or too many.
            probes *= 2
            
            if probes >= shell_stars:
                probes //= 2
                # If there are more probes than probes required, the additional probes 
                # will help the build, reducing the time required to build the probes.
                probes_required = shell_stars - probes
                steps += (probes_required / probes)
                break
            else:
                steps += 1.

        return steps * self.build_time

    def _volume(self, radius = 0.0) -> float:
        return 4/3 * np.pi * radius**3
    
    def _area(self, radius = 0.0) -> float:
        return 4 * np.pi * radius**2
    
    def _stars_per_volume(self, volume) -> float:
        return self.star_count * (volume / self.galaxy_volume)

