# Diffusion Analysis Script
# This script calculates the diffusion lengths of monomer species during the diffusion
# simulation.

import gsd.hoomd
import numpy as np
from datetime import datetime
import signac

class DiffusionAnalyzer():
    def __init__(self,job,dt):
        self.dt = dt
        self.job = job
        self.kT = job.sp['temperature']
        self.rho = job.sp['density']
        self.crosslinker_percent = job.sp['crosslinker_percent']
        self.N_monomers = job.sp['N_monomers']
        self.monomer_size = job.sp['monomer_size']
        self.extender_size = job.sp['extender_size']
        self.angle_constant = job.sp['angle_constant']

        self.radical_number_percent = job.sp['radical_percent']
        self.r_cut_reaction = job.sp['r_cut_reaction']
        
        self.chain_transfer_probability= job.sp['chain_transfer_probability']
        self.chain_side_reaction_probability= job.sp['chain_side_reaction_probability']

        self.polymerize_period = job.sp['polymerize_period'] 

        # file names 
        self.equi_gsd_file = job.fn('equi.gsd')
        self.polymerize_gsd_file = job.fn('polymerize.gsd')
        self.output_txt =job.fn('polymerize.txt')
        self.output_run_txt =job.fn('run.txt')
        self.run_gsd_file = job.fn('run.gsd')
        self.contract_bonds_gsd_file = job.fn('contract_bonds.gsd')
        self.diffuse_gsd_file = job.fn('diffuse.gsd')
        self.diffuse_output_txt = job.fn('diffuse.txt')

        if "custom_action" in self.job.sp["polymerization_method"]:
            import scripts.customAction.polymerize as polymerize
    '''
    def calculate_MSD(self):
        print("starting MSD analysis on job: ",self.job,flush=True)
        traj = gsd.hoomd.open(self.job.fn('diffuse.gsd'),mode='r')
        timestep = traj[0].configuration.step
        print(timestep)

        N_frames = len(traj)
        N_frames_analyzed = int(N_frames/2)
        thiol_id = traj[0].particles.types.index('Thiol')
        ene_id = traj[0].particles.types.index('Ene_C')
        tetra_id = traj[0].particles.types.index('C')

        #account for dummy particles
        dummy_id = traj[0].particles.types.index('Dummy')
        frame_ids = traj[0].particles.typeid[traj[0].particles.typeid != dummy_id]
        n_reactant_particles = np.sum(traj[0].particles.typeid == thiol_id) + np.sum(traj[0].particles.typeid == ene_id) + np.sum(traj[0].particles.typeid == tetra_id)

        print("Number of reactant particles:", n_reactant_particles, flush=True)
        dithiol_mask = np.zeros((n_reactant_particles,),dtype=bool)
        diene_mask = np.zeros((n_reactant_particles,),dtype=bool)
        tetrathiol_mask = np.zeros((n_reactant_particles,),dtype=bool)
        print("length of masks:", len(dithiol_mask), len(diene_mask), len(tetrathiol_mask), flush=True)

        #diene's mask is easy, just look for Ene_C
        diene_mask[frame_ids == ene_id] = True
        print("length of diene_mask:", len(diene_mask), flush=True)

        #tetrathiol's mask is easy, just look for C, and then include the four thiols before it
        tetrathiol_mask[frame_ids == tetra_id] = True
        C_indexes = np.where(frame_ids == tetra_id)[0]
        for i in C_indexes:
            tetrathiol_mask[i-4:i] = True
        print("length of tetrathiol_mask:", len(tetrathiol_mask), flush=True)

        # dithiol's mask is the thiol index mask, minus the tetrathiol mask
        thiol_mask = np.zeros((n_reactant_particles,),dtype=bool)
        thiol_mask[frame_ids == thiol_id] = True
        print("length of thiol_mask:", len(thiol_mask), flush=True)
        dithiol_mask = thiol_mask & ~tetrathiol_mask #(bitwise AND NOT operation)

        # check the masks. The xor of all masks should all be true
        # if np.any(not(np.logical_xor(dithiol_mask,np.logical_xor(diene_mask,tetrathiol_mask)))):
        #     print("Inconsistent masks detected", flush=True)
        #     print("dithiol_mask:", dithiol_mask, flush=True)
        #     print("diene_mask:", diene_mask, flush=True)
        #     print("tetrathiol_mask:", tetrathiol_mask, flush=True)
        #     exit(1)

        MSDArr_dithiol = np.zeros((N_frames_analyzed,))
        MSDArr_diene = np.zeros((N_frames_analyzed,))
        MSDArr_tetrathiol = np.zeros((N_frames_analyzed,))

        dithiol_pos_arr = []
        diene_pos_arr = []
        tetrathiol_pos_arr = []

        for j,frame in enumerate(traj):
            # Quantities in the frame are numpy arrays
            Box = frame.configuration.box[0:3]
            pos = frame.particles.position[frame.particles.typeid != dummy_id]  # remove dummy particles
            # vel = frame.particles.velocity
            ids = frame.particles.typeid[frame.particles.typeid != dummy_id]  # remove dummy particles
            images = frame.particles.image[frame.particles.typeid != dummy_id]  # remove dummy particles

            # account for periodic boundary conditions
            pos = pos + images * Box

            #down_select positions
            tmp_dithiol_pos = pos[dithiol_mask]
            tmp_dithiol_pos = tmp_dithiol_pos.reshape((-1,2,3)) # reshape to pairs of thiols
            tmp_diene_pos = pos[diene_mask]
            tmp_diene_pos = tmp_diene_pos.reshape((-1,4,3)) # reshape to 4 members of dienes
            tmp_tetrathiol_pos = pos[tetrathiol_mask]
            tmp_tetrathiol_pos = tmp_tetrathiol_pos.reshape((-1,5,3)) # reshape to 5 members of tetrathiols

            # average the positions of the components over pbc (using the images considers PBC)
            dithiol_pos = np.mean(tmp_dithiol_pos, axis=1)
            diene_pos = np.mean(tmp_diene_pos, axis=1)
            tetrathiol_pos = np.mean(tmp_tetrathiol_pos, axis=1)

            # save the positions
            dithiol_pos_arr.append(dithiol_pos)
            diene_pos_arr.append(diene_pos)
            tetrathiol_pos_arr.append(tetrathiol_pos)

        start_time = datetime.now()
        # Calculate the MSD
        del_frame_min = 1000
        for delT in range(0,N_frames_analyzed):
            if delT % 10 == 0:
                print("Calculating MSD for delT = ",delT,"/",N_frames_analyzed,flush=True)
                print("Time elapsed: ",datetime.now()-start_time,flush=True)
            if delT < 10:
                MSDArr_dithiol[delT] = 0
                MSDArr_diene[delT] = 0
                MSDArr_tetrathiol[delT] = 0
                continue
            
            tmp_dithiol = []
            tmp_diene = []
            tmp_tetrathiol = []
            for t0 in range(0,N_frames_analyzed-delT):
                tmp_dithiol.append(calc_MSD(dithiol_pos_arr[t0], dithiol_pos_arr[t0+delT]))
                tmp_diene.append(calc_MSD(diene_pos_arr[t0], diene_pos_arr[t0+delT]))
                tmp_tetrathiol.append(calc_MSD(tetrathiol_pos_arr[t0], tetrathiol_pos_arr[t0+delT]))

            MSDArr_dithiol[delT] = np.average(tmp_dithiol)
            MSDArr_diene[delT] = np.average(tmp_diene)
            MSDArr_tetrathiol[delT] = np.average(tmp_tetrathiol)

        # Average, save, ... results
        gsd_period = traj[1].configuration.step - traj[0].configuration.step
        deltat_arr = np.arange(0,N_frames_analyzed*self.dt*gsd_period,self.dt*gsd_period)

        np.savetxt(self.job.fn('MSD_dithiol.txt'), np.c_[deltat_arr,MSDArr_dithiol], header="time, MSD", fmt='%f')
        np.savetxt(self.job.fn('MSD_diene.txt'), np.c_[deltat_arr,MSDArr_diene], header="time, MSD", fmt='%f')
        np.savetxt(self.job.fn('MSD_tetrathiol.txt'), np.c_[deltat_arr,MSDArr_tetrathiol], header="time, MSD", fmt='%f')

        print("MSD analysis done in ",datetime.now()-start_time,flush=True)

    def calculate_diffusion_coeff(self):
        # read the MSD data
        MSD_dithiol = np.loadtxt(self.job.fn('MSD_dithiol.txt'), skiprows=1)
        MSD_diene = np.loadtxt(self.job.fn('MSD_diene.txt'), skiprows=1)
        MSD_tetrathiol = np.loadtxt(self.job.fn('MSD_tetrathiol.txt'), skiprows=1)

        # calculate the diffusion coefficients
        # D = MSD / (6 * t) for 3D diffusion

        #I can either make a linear fit to the MSD data or just average the derivative
        time_minimum = 0
        time_maximum = np.inf

        MSD_dithiol = MSD_dithiol[(MSD_dithiol[:,0] >= time_minimum) & (MSD_dithiol[:,0] <= time_maximum)]
        MSD_diene = MSD_diene[(MSD_diene[:,0] >= time_minimum) & (MSD_diene[:,0] <= time_maximum)]
        MSD_tetrathiol = MSD_tetrathiol[(MSD_tetrathiol[:,0] >= time_minimum) & (MSD_tetrathiol[:,0] <= time_maximum)]

        # average the derivative
        D_dithiol = np.mean(np.gradient(MSD_dithiol[:,1], MSD_dithiol[:,0])) / 6
        D_diene = np.mean(np.gradient(MSD_diene[:,1], MSD_diene[:,0])) / 6
        D_tetrathiol = np.mean(np.gradient(MSD_tetrathiol[:,1], MSD_tetrathiol[:,0])) / 6

        # save the diffusion coefficients
        np.savetxt(self.job.fn('diffusion_coeff.txt'), np.c_[D_dithiol, D_diene, D_tetrathiol], header="D_dithiol, D_diene, D_tetrathiol", fmt='%f')
    '''
    def calculate_diffusion_length_at_delta_timestep(self, timesteps):
        """
        Calculates the average diffusion length for dithiols, dienes, and tetrathiols
        over the given time interval specified by timesteps. This is done by computing
        the Mean Squared Displacement (MSD) during the diffusion simulation.

        The diffusion lengths are saved in 'diffusion_length.txt' file within the
        job directory.
        """
        print("starting MSD analysis on job: ",self.job, f"at {timesteps} timesteps",flush=True)
        traj = gsd.hoomd.open(self.job.fn('diffuse.gsd'),mode='r')
        timestep = traj[0].configuration.step
        print(timestep)

        N_frames = len(traj)
        thiol_id = traj[0].particles.types.index('Thiol')
        ene_id = traj[0].particles.types.index('Ene_C')
        tetra_id = traj[0].particles.types.index('C')

        #account for dummy particles
        dummy_id = traj[0].particles.types.index('Dummy')
        frame_ids = traj[0].particles.typeid[traj[0].particles.typeid != dummy_id]
        n_reactant_particles = np.sum(traj[0].particles.typeid == thiol_id) + np.sum(traj[0].particles.typeid == ene_id) + np.sum(traj[0].particles.typeid == tetra_id)

        print("Number of reactant particles:", n_reactant_particles, flush=True)
        dithiol_mask = np.zeros((n_reactant_particles,),dtype=bool)
        diene_mask = np.zeros((n_reactant_particles,),dtype=bool)
        tetrathiol_mask = np.zeros((n_reactant_particles,),dtype=bool)
        print("length of masks:", len(dithiol_mask), len(diene_mask), len(tetrathiol_mask), flush=True)

        #diene's mask is easy, just look for Ene_C
        diene_mask[frame_ids == ene_id] = True
        print("length of diene_mask:", len(diene_mask), flush=True)

        #tetrathiol's mask is easy, just look for C, and then include the four thiols before it
        tetrathiol_mask[frame_ids == tetra_id] = True
        C_indexes = np.where(frame_ids == tetra_id)[0]
        for i in C_indexes:
            tetrathiol_mask[i-4:i] = True
        print("length of tetrathiol_mask:", len(tetrathiol_mask), flush=True)

        # dithiol's mask is the thiol index mask, minus the tetrathiol mask
        thiol_mask = np.zeros((n_reactant_particles,),dtype=bool)
        thiol_mask[frame_ids == thiol_id] = True
        print("length of thiol_mask:", len(thiol_mask), flush=True)
        dithiol_mask = thiol_mask & ~tetrathiol_mask #(bitwise AND NOT operation)

        # check the masks. The xor of all masks should all be true
        if np.sum(np.logical_xor(dithiol_mask,np.logical_xor(diene_mask,tetrathiol_mask))) != len(dithiol_mask):
            print("Inconsistent masks detected", flush=True)
            print("dithiol_mask:", dithiol_mask, flush=True)
            print("diene_mask:", diene_mask, flush=True)
            print("tetrathiol_mask:", tetrathiol_mask, flush=True)
            exit(1)

        MSD_dithiol = np.zeros((1,))
        MSD_diene = np.zeros((1,))
        MSD_tetrathiol = np.zeros((1,))

        dithiol_pos_arr = []
        diene_pos_arr = []
        tetrathiol_pos_arr = []

        for j,frame in enumerate(traj):
            # Quantities in the frame are numpy arrays
            Box = frame.configuration.box[0:3]
            pos = frame.particles.position[frame.particles.typeid != dummy_id]  # remove dummy particles
            # vel = frame.particles.velocity
            ids = frame.particles.typeid[frame.particles.typeid != dummy_id]  # remove dummy particles
            images = frame.particles.image[frame.particles.typeid != dummy_id]  # remove dummy particles

            # account for periodic boundary conditions
            pos = pos + images * Box

            #down_select positions
            tmp_dithiol_pos = pos[dithiol_mask]
            tmp_dithiol_pos = tmp_dithiol_pos.reshape((-1,2,3)) # reshape to pairs of thiols
            tmp_diene_pos = pos[diene_mask]
            tmp_diene_pos = tmp_diene_pos.reshape((-1,4,3)) # reshape to 4 members of dienes
            tmp_tetrathiol_pos = pos[tetrathiol_mask]
            tmp_tetrathiol_pos = tmp_tetrathiol_pos.reshape((-1,5,3)) # reshape to 5 members of tetrathiols

            # average the positions of the components over pbc (using the images considers PBC)
            dithiol_pos = np.mean(tmp_dithiol_pos, axis=1)
            diene_pos = np.mean(tmp_diene_pos, axis=1)
            tetrathiol_pos = np.mean(tmp_tetrathiol_pos, axis=1)

            # save the positions
            dithiol_pos_arr.append(dithiol_pos)
            diene_pos_arr.append(diene_pos)
            tetrathiol_pos_arr.append(tetrathiol_pos)

        start_time = datetime.now()
        # Calculate the MSD
        # del_frame_min = 1000
        for delT in [timesteps]:
            if delT % 10 == 0:
                print("Calculating MSD for delT = ",delT,"/",N_frames,flush=True)
                print("Time elapsed: ",datetime.now()-start_time,flush=True)
            
            tmp_dithiol = []
            tmp_diene = []
            tmp_tetrathiol = []
            for t0 in range(0,N_frames-delT):
                tmp_dithiol.append(calc_MSD(dithiol_pos_arr[t0], dithiol_pos_arr[t0+delT]))
                tmp_diene.append(calc_MSD(diene_pos_arr[t0], diene_pos_arr[t0+delT]))
                tmp_tetrathiol.append(calc_MSD(tetrathiol_pos_arr[t0], tetrathiol_pos_arr[t0+delT]))

            MSD_dithiol[0] = np.average(tmp_dithiol)
            MSD_diene[0] = np.average(tmp_diene)
            MSD_tetrathiol[0] = np.average(tmp_tetrathiol)

        diff_length_dithiol = np.sqrt(MSD_dithiol[0])
        diff_length_diene = np.sqrt(MSD_diene[0])
        diff_length_tetrathiol = np.sqrt(MSD_tetrathiol[0])

        # Average, save, ... results
        np.savetxt(self.job.fn(f'diffusion_length_at_delta_timesteps_{timesteps}.txt'), np.c_[diff_length_dithiol, diff_length_diene, diff_length_tetrathiol], header="dithiol, diene, tetrathiol", fmt='%f')

        print(f"diffusion length at deltaTimesteps {timesteps} analysis done in ",datetime.now()-start_time,flush=True)


# from azobenzene project
def calc_MSD(positions, ref_positions):
    """
    Calculates the Mean Squared Displacement (MSD) between two sets of positions.
    Assumes that both position arrays are of the same length, and that the row indices
    correspond to the same particles.
    """
    sqDistSum = 0
    for i,particle_pos in enumerate(positions):
        vec = positions[i] - ref_positions[i]
        sqDistSum = sqDistSum + np.dot(vec,vec)
    MSD = sqDistSum/len(positions)
    return MSD