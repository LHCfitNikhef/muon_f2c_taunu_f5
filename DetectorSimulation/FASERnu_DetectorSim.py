import os
import numpy as np
import json
from skhep.math.vectors import LorentzVector, Vector3D
import pyhepmc

class FASERvDetectorSimulation():
    
    def __init__(
            self,
            debug=False,
        ):
        
        # save inputs
        self.debug=debug
        self.maxevent_print=5
       
        # initialize default detector
        self.set_detector()

    def set_detector(
            self,
            xmin=-0.125,
            xmax=0.125,
            ymin=-0.15,
            ymax=0.15,
            zmin=0.2,
            zmax=0.8,
            zfront=0.2,
            zback=0.2,
            particle_minimum_momentum=0.3,
            particle_maximum_angle=1,
            sigma_track_momentum=0.46,
            sigma_shower_energy=0.5,
            charmtau_id_probability=0.75,
            interaction_length=0.1,
        ):
    
        # dimensions of detector
        self.xmin, self.xmax = xmin, xmax
        self.ymin, self.ymax = ymin, ymax
        self.zmin, self.zmax = zmin, zmax
        self.zfront, self.zback = zfront, zback
        self.interaction_length = interaction_length
        
        # acceptance
        self.particle_minimum_momentum = particle_minimum_momentum
        self.particle_maximum_angle = particle_maximum_angle
        
        # momentum/energy resolution
        self.sigma_track_momentum=sigma_track_momentum
        self.sigma_shower_energy=sigma_shower_energy
        
        # charm ID
        self.charmtau_id_probability = charmtau_id_probability
        
        # detectable particles
        self.detectable = [22,11,13,15,211,321,2212,3222,3112,3312,411,431,4122,111]
        self.charged = [11,13,15,211,321,2212,3222,3112,3312,411,431,4122]
        self.tracks = [13,15,211,321,2212,3222,3112,3312,411,431,4122]
        self.showers = [11, 22, 111]
        self.charmtau = [15, 411,431,4122]
        
    def simulate(
            self,
            inputfilename,
            outputfilename,
            do_return = True
        ):
        
        # define input/output
        self.inputfilename = inputfilename
        self.outputfilename = outputfilename
        self.events = {}
        
        # delete outputfile if it already esists
        if os.path.exists(self.outputfilename):
            os.remove(self.outputfilename)

        # read and process HEPMC file
        print ("-- open file", self.inputfilename)
        with pyhepmc.open(self.inputfilename, "r") as reader:
            #initialize event counter
            ievent=0
            # loop over events in file
            for event in reader:
                # increase counter
                ievent+=1
                # update on how many events have been processed
                if (ievent % 1000 == 0): print ("   processed", ievent, " Events")
                # analyse events
                self.read_event(event, ievent)
        
        # write output
        with open(outputfilename, 'w') as f:
            json.dump(self.events, f, indent=2)
        print ("-- output is stored in", self.outputfilename)
        
        if do_return: return self.events

    def smear_track(self, central):
        factor=0
        while factor<=0: factor=np.random.normal(1, self.sigma_track_momentum)
        return central*factor
    
    def smear_shower(self, central):
        factor=0
        while factor<=0: factor=np.random.normal(1, self.sigma_shower_energy)
        return central*factor
    
    def charm_tau_ID(self, apid, momentum):
        if apid not in self.charmtau: return 0
        return self.charmtau_id_probability
    
    def muon_ID(self, apid, status, momentum, zposition):
        if apid not in self.tracks: return 0
        if apid == 13: return 1
        if status == 1: distance = self.zback + self.zmax - zposition
        if status == 4: distance = self.zfront + zposition - self.zmin
        prob_no_int = np.exp(-distance/self.interaction_length)
        rand = np.random.uniform(0,1)
        if rand < prob_no_int: return 1
        else: return 0
    
    def read_event(self,hepmc_event,ievent,):
        
        # initialize event particles
        self.hepmc_particles=hepmc_event.particles
        self.hepmc_vertices=hepmc_event.vertices
        self.particles = []
        
        # get incoming particle energy
        self.initial_particle = self.hepmc_particles[1]
        
        # read all visible particles at truth level
        self.events_detector=[]
        
        #randomly choose location of event
        posx = np.random.uniform(self.xmin,self.xmax)
        posy = np.random.uniform(self.ymin,self.ymax)
        posz = np.random.uniform(self.zmin,self.zmax)
        self.position = [posx,posy,posz]
        
        # output input particles for debugging
        if self.debug and ievent<=self.maxevent_print:
            print ("-- original particles for event", ievent)
            for particle in self.hepmc_particles:
                print ("  ", particle)
            print ("")

        # get truth leptons
        self.truth_leptons = []
        for iparticle,particle in enumerate(self.hepmc_particles):
            if particle.status not in [1,4]: continue
            if abs(particle.pid) not in [13]: continue
            p=particle.momentum
            self.truth_leptons.append([p[0],p[1],p[2],p[3]])
            if len(self.truth_leptons)>1: break
            
        # output input muons for debugging
        if self.debug and ievent<=self.maxevent_print:
            print ("-- original muons", ievent)
            print ("  ", self.truth_leptons[0])
            print ("  ", self.truth_leptons[1])
            print ("")
        
        #get particles after detector simulation
        self.event = []
        self.reco_leptons = []
        ntracks = 0
        for iparticle,particle in enumerate(self.hepmc_particles):
                        
            # ignore unstable particles, except tau, D+ and Ds+
            if particle.status not in [1,4]: continue
            apid = abs(particle.pid)
        
            # continue when not detectable particle
            if apid not in self.detectable: continue
                
            # extract momentum
            p=particle.momentum
            momentum=LorentzVector(p[0],p[1],p[2],p[3])
            
            # ignore soft and large-angle particles
            if momentum.p<0.5*self.particle_minimum_momentum: continue
            if momentum.pt/np.abs(momentum.pz) > self.particle_maximum_angle: continue

            # smear (only absolute momentum)
            sintheta_x  = momentum.px/momentum.p
            sintheta_y  = momentum.py/momentum.p
            if apid in self.tracks: p_measured = self.smear_track(momentum.p)
            if apid in self.showers: p_measured = self.smear_shower(momentum.p)
            px_measured = p_measured * sintheta_x
            py_measured = p_measured * sintheta_y
            pz_measured = np.sqrt(p_measured**2 - px_measured**2 - py_measured**2)
            if p_measured < self.particle_minimum_momentum: continue
                
            # count charged tracks
            if (apid in self.charged) and (particle.status==1): ntracks+=1
            
            # particle ID
            is_charm = self.charm_tau_ID(apid, momentum)
            is_muon = self.muon_ID(apid, particle.status, momentum, self.position[2])
            if   apid==11: pid_measured = 'electron'
            elif apid==111: pid_measured = 'pi0'
            elif is_charm: pid_measured = 'charm'
            elif is_muon: pid_measured = 'muon'
            else: pid_measured = 'track'
                
            # status
            outgoing = 0 if particle.status==4 else 1

            # save smeared particle
            self.event.append([outgoing, pid_measured, px_measured, py_measured, pz_measured])
            
            # save smeared muons
            if is_muon:
                energy_smeared = np.sqrt(px_measured**2+ py_measured**2+ pz_measured**2+ 0.1057**2)
                if len(self.reco_leptons)==2 and self.reco_leptons[1][3] < energy_smeared: self.reco_leptons.pop()
                self.reco_leptons.append([px_measured, py_measured, pz_measured, energy_smeared])
              
        # add None's in case of missing reco muons
        while len(self.reco_leptons)<2: self.reco_leptons.append(None)
        
        # output reco muons for debugging
        if self.debug and ievent<=self.maxevent_print:
            print ("-- smeared muons", ievent)
            for particle in self.reco_leptons: print ("  ", particle)
            print ("")
            
        # kinematics
        l0, l1 = np.array(self.truth_leptons[0]), np.array(self.truth_leptons[1])
        truth_m  = 0.938
        truth_q2 = (l0[0]-l1[0])**2+(l0[1]-l1[1])**2+(l0[2]-l1[2])**2-(l0[3]-l1[3])**2
        truth_v  = l0[3]-l1[3]
        truth_x  = truth_q2 / (2*truth_m*truth_v)
        truth_y  = truth_v / l0[3]
        truth_w2 = truth_m**2 + truth_q2*(1-truth_x)/truth_x
        
        reco_m, reco_q2, reco_v, reco_x, reco_y, reco_w2 = None, None, None, None, None, None
        if len(self.reco_leptons)==2:
            l0, l1 = np.array(self.reco_leptons[0]), np.array(self.reco_leptons[1])
            reco_m  = 0.938
            reco_q2 = (l0[0]-l1[0])**2+(l0[1]-l1[1])**2+(l0[2]-l1[2])**2-(l0[3]-l1[3])**2
            reco_v  = l0[3]-l1[3]
            reco_x  = reco_q2 / (2*reco_m*reco_v)
            reco_y  = reco_v / l0[3]
            reco_w2 = truth_m**2 + reco_q2*(1-reco_x)/reco_x
        
        # write event
        event = {}
        event['ntracks'] = ntracks
        event['truth_muons'] = {'in': self.truth_leptons[0], 'out': self.truth_leptons[1]}
        event['reco_muons'] = {'in': self.reco_leptons[0], 'out': self.reco_leptons[1]}
        event['truth_kinematics'] = {'x': truth_x, 'y': truth_y, 'q2': truth_q2}
        event['reco_kinematics'] = {'x': reco_x, 'y': reco_y, 'q2': reco_q2}
        event['reco_event'] = self.event
        event['vertex_position'] = self.position
        self.events[ievent] = event

    def write_event(self):

        # write file
        with open(self.outputfilename, "ab") as outputfile:
            outputfile.write("New Event with Enu = " + str(self.neutrino_energy) + " \n")
            outputfile.write("position: "+str(self.position.x)+" "+str(self.position.y)+" "+str(self.position.z)+" \n")
            if len(self.events_detector)>0: np.savetxt(outputfile, self.events_detector, fmt='%12.i %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f' ,delimiter=' ')
