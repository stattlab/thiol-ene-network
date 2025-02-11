#The pressure tensor will NOT be output to the deform.log file. In hoomd 4, it must be 
# accessed via a gsd file. It will be output as pressures.gsd. Use these lines to read:
#   import gsd.hoomd
#   logs = gsd.hoomd.read_log(***Path to pressures.gsd***)
#   pressure_tensors = logs["log/md/compute/ThermodynamicQuantities/pressure_tensor"]

import numpy as np
import sys
import gsd,gsd.hoomd
import sys, os 
sys.path.insert(0,'/Users/statt/Programs/hoomd_fork_4.0.1/')
import hoomd
# from timeit import default_timer as timer
# from scripts.system import System
from pathlib import Path

class BoxWriter(hoomd.custom.Action):
    # def __init__(self):
    #     super().attach(simulation)

    def attach(self,simulation):
        super().attach(simulation)

    @hoomd.logging.log(category='scalar')
    def Lx(self):
        snap = self._state.get_snapshot()
        if snap.communicator.rank == 0:
            return snap.configuration.box[0]
        
    @hoomd.logging.log(category='scalar')
    def Ly(self):
        snap = self._state.get_snapshot()
        if snap.communicator.rank == 0:
            return snap.configuration.box[1]
        
    @hoomd.logging.log(category='scalar')
    def Lz(self):
        snap = self._state.get_snapshot()
        if snap.communicator.rank == 0:
            return snap.configuration.box[2]

    def act(self, timestep):
        pass



def main(gsd_file_path,frame_number):
    parent_path = str(Path(gsd_file_path).parent.absolute())

    try:
        device = hoomd.device.GPU()
    except:
        device = hoomd.device.CPU()

    sim = hoomd.Simulation(device=device, seed=1)
    sim.create_state_from_gsd(filename=gsd_file_path, frame=int(frame_number))

    # ------------------- Box resizer -------------------
    # Parameters taken from Brandon's previous parameter scans for NPT deformation of 
    # triblock system, except step size scaled down by a fact or 10. 1000 -> 100
    step_size = 100
    delta_lam = 0.001
    # NVT deformation
    lam = 3.0
    deform_time = (lam - 1) / delta_lam * step_size
    box_variant = hoomd.variant.Ramp(A=0,B=1,t_start=sim.timestep,t_ramp=int(deform_time))
    initial_box = sim.state.box
    final_box = [initial_box.Lx*lam,
                 initial_box.Ly*1/np.sqrt(lam),
                 initial_box.Lz*1/np.sqrt(lam)]
    deform_box_interpolate = hoomd.variant.box.Interpolate(initial_box=initial_box,
                                                           final_box=final_box,
                                                           variant=box_variant)

    box_resize = hoomd.update.BoxResize(trigger=hoomd.trigger.Periodic(step_size),
                                        box=deform_box_interpolate)
    sim.operations.updaters.append(box_resize)

    #-----------------------Force fields------------------------
    particle_types = sim.state.particle_types
    cell = hoomd.md.nlist.Cell(buffer=0.4)

    lj = hoomd.md.pair.LJ(nlist=cell)
    lj.params[(particle_types, particle_types)] = dict(epsilon=1.0,sigma=1.0)
    lj.r_cut[particle_types, particle_types] = 2.5

    lj.params[(particle_types, 'Dummy')] = dict(epsilon=0.0,sigma=0.0)
    lj.r_cut[particle_types, 'Dummy'] = 0
    lj.mode = 'shift'

    fene = hoomd.md.bond.FENEWCA()
    fene.params[sim.state.bond_types] = dict(k=30,r0=1.5,epsilon=1.0, sigma=1.0, delta=0.0)
    fene.params['Dummy'] = dict(k=0,r0=1.5,epsilon=0.0, sigma=1.0, delta=0.0)

    if len(sim.state.angle_types)==0:
        FJ_system=True
    else:
        FJ_system=False

    if FJ_system==False:
        angle_constant = 100
        print("Could read angle_constant, but assuming 100")
        cosinesq = hoomd.md.angle.CosineSquared()
        cosinesq.params[sim.state.angle_types] = dict(k=angle_constant, t0=np.pi*110/180)# https://www.sciencedirect.com/science/article/pii/S0032386110003642?ref=cra_js_challenge&fr=RR-1
        cosinesq.params['Dummy'] = dict(k=0.0001, t0=np.pi)  # k>0 to make warning go away (should not do anything)


    #----------------------Integration parameters-------------------------

    kT = 1 # Could read, but just assuming job.sp["temperature"]
    print(f"Could read, but assuming temperature = {kT}")
    dt = 0.005 # Also taken from Brandon's deformation parameter scan
    types_to_integrate =  hoomd.filter.Type(sim.state.particle_types[:-1]) # everything but "Dummy" particles

    integrator = hoomd.md.Integrator(dt=dt)
    
    if FJ_system==False:
        integrator.forces = [lj,fene,cosinesq]
    else:
        integrator.forces = [lj,fene]

    mttk = hoomd.md.methods.thermostats.MTTK(kT=kT,tau=100*dt)
    nvt = hoomd.md.methods.ConstantVolume(filter=types_to_integrate,thermostat=mttk)
    integrator.methods = [nvt]

    sim.operations.integrator = integrator

    #------------------------Configure output---------------------------

    logger = hoomd.logging.Logger(categories=['scalar','string'])
    logger.add(sim, quantities=['timestep','tps'])

    thermo = hoomd.md.compute.ThermodynamicQuantities(types_to_integrate)
    sim.operations.computes.append(thermo)
    logger.add(thermo, ['kinetic_temperature','kinetic_energy', 'potential_energy','pressure','volume'])
    
    box_writer_action = BoxWriter()
    logger.add(box_writer_action, ['Lx','Ly','Lz'])
    box_writer = hoomd.write.CustomWriter(action=box_writer_action,
                                             trigger=hoomd.trigger.Periodic(10))
    sim.operations.writers.append(box_writer)

    log_writer = hoomd.write.Table(output=open(parent_path + '/deform.log','w'),
                                    trigger=hoomd.trigger.Periodic(10),
                                    logger=logger,)
    sim.operations.writers.append(log_writer)


    logger = hoomd.logging.Logger(categories=['sequence'])
    logger.add(thermo, ['pressure_tensor'])
    gsd_writer = hoomd.write.GSD(filename=parent_path + '/deform.gsd',
                                    trigger=hoomd.trigger.Periodic(int(deform_time/50)),
                                    dynamic=['property','momentum','topology','attribute'],
                                    mode='wb',
                                    logger=logger)
    sim.operations.writers.append(gsd_writer)

    #add pressure to a new gsd writer
    logger_pressure = hoomd.logging.Logger(categories = ['scalar','sequence'])
    logger_pressure.add(thermo, ['pressure_tensor'])
    gsd_writer_pressure = hoomd.write.GSD(filename=parent_path+"/pressures.gsd",
                                          filter = hoomd.filter.Null(),
                                          mode = "wb", 
                                          trigger=hoomd.trigger.Periodic(10),
                                          logger=logger_pressure)
    sim.operations.writers.append(gsd_writer_pressure)

    #--------------------------Run the simulation----------------------

    sim.run(int(deform_time/4))

    # sim.run(deform_time+1)


    print("Deformation done! Written gsd file to: ",parent_path + "/deform.gsd")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python3 deform.py gsd_file_path frame_number")
        sys.exit(1)
    main(str(sys.argv[1]), str(sys.argv[2]))