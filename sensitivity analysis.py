# %% imports
import numpy as np
import pandas as pd
import math
import pint
import os
import sys

from functions import *
from rocketcea.cea_obj import CEA_Obj
from SALib.analyze import sobol
from SALib.sample import saltelli
from SALib.sample import latin

DISPLAY = False

class HidePrints:
    def __enter__(self):
        self._original_stdout = sys.stdout
        sys.stdout = open(os.devnull, 'w')

    def __exit__(self, et, ev, etb):
        sys.stdout.close()
        sys.stdout = self._original_stdout

def weight_analysis(power_density_motor, power_density_inverter, efficiency_motor, 
                    structural_margin, energy_efficiency_battery, efficiency_inverter,
                    power_density_battery, energy_density_battery, efficiency_pump, rotating_speed):

    # %% '''------------------Assumptions-------------------- '''
    global g
    global units

    characteristic_length = 1.145 * units.meter # from Lee


    # %% '''------------------Design Variables-------------------- '''
    thrust = 500 * units.newton
    pressure_chamber = 20.0 * units.bar # a propellant supply pressure of 4–30 bar
    expansion_ratio = 200.0 # 200 # 150-400 typically for high altitude nozzle 
    mixture_ratio = 2.6 # 2.6
    burn_time = 600 * units.second
    prop_margin = 0.1 # from Lee
    ullage = 0.2 # from Lee

    theta_chamber = 15 * units.degree # this is the angle of the chamber frustrum going into the throat. Assumed same as nozzle angle

    vel_prop = 1.2 * units('meters/second') # In single-phase, liquid lines, the flowspeed of 0.9–4.5m/s is recommended empirically. Set to max to minimize pipe diameter


    # %% '''------------------Material Properties-------------------- '''
    density_chamber_material = (8.02 + 7.83) / 2 * units('megagram/meter**3') # 301 stainless steel per Lee paper https://www.azom.com/properties.aspx?ArticleID=960
    yield_strength_chamber = 276 * units('megapascal') # from Lee

    yield_strength_pipe = 55 * units('megapascal') # https://www.fergusonperf.com/the-perforating-process/material-information/specialized-aluminum/6061-aluminium-alloy/#:~:text=elongation%20of%2016%25.-,6061%2DT6,it%20has%20elongation%20of%2010%25.
    safety_factor_pipe = 3 # assumption from Lee
    density_pipe = 2.7 * units('gram/centimeter**3') # from Lee
    # length_pipe = 5 * units.meter # calculated from tank properties
    alpha_elbow =  45 * units.degree

    pressure_tank = 3.3 * units('bar') # technically different pressures for oxygen and fuel


    # %% '''------------------Technologies-------------------- '''
    power_density_motor = 0.875 * units('kilowatt/kilogram') # tech
    power_density_inverter = 0.650 * units('kilowatt/kilogram') # tech
    efficiency_motor = 0.87 # tech
    structural_margin = 0.2 # tech
    energy_efficiency_battery = 0.925 # tech
    efficiency_inverter = 0.85 # tech
    power_density_battery = 0.650 * units('kilowatt/kilogram') # tech
    energy_density_battery = 325 * units('watt * hour / kilogram') # tech
    efficiency_pump = 0.66 # tech
    rotating_speed = 57940 # tech

    # %% '''------------------NASA Chemical Equilibrium with Application (CEA)-------------------- '''
    with HidePrints():
        c = CEA_Obj(oxName='LOX', fuelName='RP_1')
        isp = c.get_Isp(Pc=pressure_chamber.magnitude, eps=expansion_ratio, MR=mixture_ratio) * units.seconds
        characteristic_velocity = c.get_Cstar(Pc=pressure_chamber.magnitude, MR=mixture_ratio) * units.feet * units.seconds**-1

    characteristic_velocity.ito('meter/second')

    if DISPLAY: print('isp ', isp)
    if DISPLAY: print('c* ', characteristic_velocity)
    

    # %% '''------------------Main Chamber, Nozzle-------------------- '''
    mass_flow_rate = thrust / (g * isp) # alternatively mdot = (pressure chamber * area throat) / c*
    mass_flow_rate.ito('kilogram/second')
    if DISPLAY: print('Mass flow rate: ', mass_flow_rate)

    radius_throat = np.sqrt((mass_flow_rate * characteristic_velocity) / (pressure_chamber * math.pi)).to('meters')

    mass_chamber, thickness_chamber = calc_mass_chamber(characteristic_length, radius_throat, pressure_chamber, 
                                density_chamber_material, yield_strength_chamber, theta_chamber)

    mass_chamber.ito('kilogram')
    if DISPLAY: print('Mass Chamber: ', mass_chamber)

    # %% Nozzle
    theta_nozzle = 15 * units.degree
    safety_factor_nozzle = 3 # from Lee
    thickness_nozzle = 0.8 * thickness_chamber # 0.50 * units('centimeter') # TODO research. this is a random guess
    density_nozzle_material = density_chamber_material # assumed 301 stainless like chamber
    mass_nozzle = calc_mass_nozzle(radius_throat, expansion_ratio, theta_nozzle,  
                                    safety_factor_nozzle, thickness_nozzle, density_nozzle_material)

    mass_nozzle.ito('kilogram')
    if DISPLAY: print('Mass Nozzle: ', mass_nozzle)

    if DISPLAY: print('Combined nozzle combustor mass: ', mass_chamber + mass_nozzle)


    # %% '''------------------Pump-------------------- '''
    with HidePrints():
        density_prop = c.get_Chamber_Density(Pc=pressure_chamber.magnitude, eps=expansion_ratio, MR=mixture_ratio) * units.pounds * units.feet**-3

    # density_prop.ito('kilogram/meter^3')
    # if DISPLAY: print('DENSITY PROP CEA: ', density_prop)
    density_prop = 915 * units('kilogram/meter^3') # from https://en.wikipedia.org/wiki/RP-1

    density_prop.ito('kilogram/meter^3')

    with HidePrints():
        viscosity_chamber = c.get_Chamber_Transport(Pc=pressure_chamber.magnitude, eps=expansion_ratio, MR=mixture_ratio)[1] * units('millipoise')

    pressure_change_pump = 2 * (pressure_chamber.to('pascal') - pressure_tank.to('pascal')) # TODO include real pressure losses bottom pg 6. Losses from injector?
    vol_flow_rate = mass_flow_rate / density_prop

    head_rise = pressure_change_pump / (density_prop * g)
    head_rise.ito('meter')
    mass_pump = calc_mass_pump(head_rise, vol_flow_rate, rotating_speed, units)

    mass_pump.ito('kilogram')
    if DISPLAY: print('Mass pump: ', mass_pump.to('kilogram'))

    power_pump_req = (pressure_change_pump * mass_flow_rate) / (density_prop * efficiency_pump) # head rise or pressure rise? seeem to be used synonymously in paper
    power_pump_req.ito('watts')
    if DISPLAY: 
        print('Power required pump: ', power_pump_req)
        print('    Pressure change pump ',pressure_change_pump)
        print('    Mass flow rate ', mass_flow_rate)
        print('    Propellant density ', density_prop)
        print('    Pump efficiency ', efficiency_pump)


    # %% '''------------------Motor-------------------- '''
    mass_motor = calc_mass_motor(power_pump_req, power_density_motor)

    mass_motor.ito('kilogram')
    if DISPLAY: print('Mass motor: ', mass_motor)
    power_motor_out = power_pump_req


    # %% '''------------------Inverter-------------------- '''
    mass_inverter = calc_mass_inverter(power_motor_out, power_density_inverter, efficiency_motor)
    if DISPLAY: print('Mass inverter: ', mass_inverter)


    # %% '''------------------Battery-------------------- '''
    power_inverter_out = power_pump_req/(efficiency_inverter * efficiency_motor) # TODO this should include the motor efficiency

    mass_battery = calc_mass_battery(structural_margin, energy_efficiency_battery, efficiency_inverter, 
                                    power_density_battery, energy_density_battery, power_inverter_out, burn_time)
    if DISPLAY: print('Mass battery: ', mass_battery)


    # %% '''------------------Propellant, Tank-------------------- '''
    # flow_rate_ox = mass_flow_rate / (1 + 1/mixture_ratio)
    mass_prop = calc_mass_prop(mass_flow_rate, burn_time, prop_margin)
    mass_prop.ito('kilogram')
    if DISPLAY: print('Mass propellant: ', mass_prop)

    density_tank = density_pipe
    yield_strength_tank = yield_strength_pipe
    mass_tank, length_pipe = calc_mass_tank(ullage, mass_prop, density_tank, 
                                pressure_tank, yield_strength_tank, density_prop)
    mass_tank.ito('kilogram')
    if DISPLAY: print('Mass tank: ', mass_tank)


    # %% '''------------------Feed Line-------------------- '''
    # assume length_pipe = height of propellant tank
    pressure_pipe = 2 * pressure_chamber + pressure_tank
    mass_feed_line = calc_feed_mass(mass_flow_rate, density_prop, vel_prop,
                                    length_pipe, yield_strength_pipe, density_pipe, viscosity_chamber, 
                                    pressure_pipe, units, safety_factor_pipe, alpha_elbow, )
    mass_feed_line.ito('kilogram')

    mass_power_system = mass_pump + mass_motor + mass_inverter + mass_battery

    if DISPLAY: 
        print('Mass feed system: ', mass_feed_line)
        print('Combined power mass (pump, motor, inverter, battery): ', mass_power_system)
        total_mass = (mass_chamber + mass_nozzle + mass_feed_line + mass_pump + mass_motor +
                    mass_inverter + mass_battery + mass_prop + mass_tank)
        print('Total mass: ', total_mass )
        print('Dry mass', total_mass - mass_prop )
        print('=======================================================')

    return(mass_power_system.magnitude)

# %%
if __name__ == '__main__':
    units = pint.UnitRegistry()
    g = g * units('meter/second**2')

    problem = {
        'num_vars': 10,
        'names': [
            'power_density_motor',
            'power_density_inverter',
            'efficiency_motor',
            'structural_margin',
            'energy_efficiency_battery',
            'efficiency_inverter',
            'power_density_battery',
            'energy_density_battery',
            'efficiency_pump',
            'rotating_speed'],
        'bounds':[
            [0.875, 2* 0.875],
            [0.650, 2* 0.650],
            [0.87, 1],
            [0.0, 0.2],
            [0.925, 1],
            [0.85, 1],
            [0.650, 2* 0.650],
            [0.325  , 2* 0.325 ],
            [0.66, 1],
            [57940, 2* 57940 ],
        ]
    }

    # %% SALib calculations
    param_values = saltelli.sample(problem, 1024)
    # print(param_values.shape)
    # print(param_values)

    Y = np.zeros([param_values.shape[0]])

    for i, X in enumerate(param_values):
        Y[i] = weight_analysis(*X)

    si = sobol.analyze(problem, Y)

    total_Si, first_Si, second_Si = si.to_df()
    total_Si.to_csv('Total.csv')
    first_Si.to_csv('First.csv')
    second_Si.to_csv('Second.csv')

    print('S1: ', si['S1'])
    print('S2: ', si['S2'])
    print('ST: ', si['ST'])

    # %% Data for JMP
    jmp_vals = latin.sample(problem, 5096)
    # print(jmp_vals.shape)
    # print(jmp_vals)

    Yj = np.zeros([jmp_vals.shape[0]])
    for i, X in enumerate(jmp_vals):
        Yj[i] = weight_analysis(*X)

    vals = pd.DataFrame(jmp_vals)
    responses = pd.DataFrame(Yj)
    out = pd.concat([vals, responses], axis='columns')
    out.columns = problem['names'] + ['Response']
    out.to_csv('to_model_data.csv')

    