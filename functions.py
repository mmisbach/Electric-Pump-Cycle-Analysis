import numpy as np
from numpy.core.einsumfunc import _find_contraction
import pandas as pd
import math
from scipy.constants import g


def calc_mass_chamber(characteristic_length, radius_throat, pressure_chamber, density_material, yield_strength, theta_chamber):
    # contraction_ratio = 8 * diam_throat ** -0.6 + 1.25

    # radius_throat = diam_throat/2
    radius_chamber = radius_throat * np.sqrt(12.126 * radius_throat.magnitude**-0.6 + 1.25)
 

    length_chamber = characteristic_length * (radius_throat / radius_chamber) **2

    thickness_chamber = (pressure_chamber * radius_chamber) / yield_strength

    mass_chamber = math.pi * density_material * thickness_chamber * ( 2* radius_chamber * length_chamber +
                             ((radius_chamber**2 - radius_throat**2)/math.tan(np.deg2rad(theta_chamber))))

    return mass_chamber, thickness_chamber


def calc_mass_nozzle(radius_throat, expansion_ratio, theta_nozzle, safety_factor_nozzle, thickness_nozzle, density_nozzle):
    radius_exhaust = math.sqrt(expansion_ratio) * radius_throat

    length_nozzle = (radius_exhaust - radius_throat) / math.tan(np.deg2rad(theta_nozzle))

    area_nozzle = math.pi * length_nozzle * (radius_exhaust + radius_throat) * safety_factor_nozzle

    mass_nozzle = thickness_nozzle * area_nozzle * density_nozzle

    return mass_nozzle


def calc_feed_mass(mass_flow_rate, density_fluid, vel_fluid, length_pipe, 
                yield_strength, density_material, viscosity, pressure_pipe, 
                units, safety_factor = 3.0, alpha=90):
    # Pipe
    diam_pipe =np.sqrt((4 * mass_flow_rate) / (math.pi * density_fluid * vel_fluid))
    # reynolds_num = (density_fluid * vel_fluid * diam_pipe) / viscosity
    # print('reynolds', reynolds_num)

    # if reynolds_num.magnitude <= 2000: friction_factor = 64/reynolds_num
    # elif reynolds_num.magnitude <= 10e5: friction_factor = 0.3164 * reynolds_num ** -0.25
    # else: friction_factor = 0.0032 + 0.221 * reynolds_num ** -0.237

    # pressure_change_pipe = friction_factor * (length_pipe / diam_pipe) * 0.5 * density_fluid * vel_fluid ** 2
    # pressure_change_pipe.ito('pascal')

    thickness_pipe = (pressure_pipe * diam_pipe) / (2 * yield_strength * safety_factor)
    thickness_pipe.ito('meter')
    
    radius_pipe = diam_pipe/2
    mass_pipe = density_material * length_pipe * math.pi * (radius_pipe**2 - (radius_pipe - thickness_pipe) **2) # not in paper TODO check units
    # print('PIPE THICKNESS', thickness_pipe)
    # print('DIAMETER PIPE', diam_pipe)
    mass_pipe.ito('kilogram')


    # Elbow
    # assumptions
    diam_elbow = diam_pipe
    thickness_elbow = thickness_pipe
    radius_interior = diam_pipe # Assumption
    density_elbow = density_material

    # if alpha < 90:
    #     loss_elbow = 30 * friction_factor * alpha * (0.0142 - 3.703e-5 * alpha)
    # else: loss_elbow = 30 * friction_factor

    # pressure_change_elbow = loss_elbow * 0.5 * density_fluid * vel_fluid**2

    radius_exterior = radius_interior + 2 * diam_elbow # Assumption same as pipe

    vol_elbow = (4 * math.pi * thickness_elbow * (radius_exterior + radius_interior) *
                (radius_exterior - radius_interior - thickness_elbow) * (alpha/360))

    mass_elbow = vol_elbow * density_elbow
    mass_elbow.ito('kilogram')

    # print('PRESSURE PIPE ', pressure_pipe)
    mass_valve = 0.05 * diam_pipe.to('millimeter').magnitude * pressure_pipe.to('megapascal').magnitude * units('kilogram')

    return mass_pipe * 2 + mass_elbow * 6 + mass_valve * 2


def calc_mass_pump(head_rise, vol_flow_rate, rotating_speed, units):
    specific_speed = rotating_speed * np.sqrt(vol_flow_rate) / (g * units('meter/second^2') * head_rise)**(3/4) # eg 21
    diameter_pump = (3.72 / specific_speed.magnitude)**(1/1.1429) * units.meter # (22)
    mass_pump = 0.4703 * np.exp(0.01072 * diameter_pump.magnitude) * units('kilogram')
    return mass_pump


def calc_mass_motor(power_pump_out, power_density):
    return power_pump_out / power_density # power_pump = power_motor_out


def calc_mass_inverter(power_motor, power_density, efficiency_motor):
    return power_motor /(efficiency_motor * power_density) # power_inverter_out = p


def calc_mass_battery(structural_margin, energy_efficiency_battery, efficiency_inverter, 
                    power_density, energy_density, power_inverter_out, burn_time):
    mass_power = (1+structural_margin) * power_inverter_out / (
                        power_density * efficiency_inverter)
    mass_energy = (1+structural_margin) * power_inverter_out * burn_time / (
                    energy_density * energy_efficiency_battery * efficiency_inverter)
    # mass_battery = max(mass_power.to('kilogram').magnitude, mass_energy.to('kilogram').magnitude) TODO fix units 
    mass_battery = max(mass_power, mass_energy)
    return mass_battery


def calc_mass_prop(mass_flow_rate, burn_time, prop_margin):
    mass_prop = (1 + prop_margin) * burn_time * mass_flow_rate
    return mass_prop


def calc_mass_tank(ullage, mass_prop, density_material, pressure_tank, yield_strength, density_prop):
    vol_tank = (1 + ullage) * (mass_prop / density_prop)
    length_tank = 2 * (3 * vol_tank / (4 * math.pi))**(1/3)

    mass_tank = (density_material * 4 * math.pi * 
                ((3 * vol_tank) / (4* math.pi))**(2/3)
                * pressure_tank / (2 * yield_strength)
                * ((3 * vol_tank) / (4 * math.pi))**(1/3))

    return mass_tank, length_tank

