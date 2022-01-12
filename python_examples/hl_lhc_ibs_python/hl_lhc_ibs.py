#
# --- HL-LHC IBS Study
#
import os
import sys
import json

import numpy as np
import pandas as pd

from config_ibs import configuration

from cpymad.madx import Madx

def twiss2df(mad):
	_tmp = mad.table['twiss']
	_df = pd.DataFrame(dict(_tmp))
	_df = _df.set_intex('name', drop=False)
	_df.index.name = ''
	return _df

def summ2df(mad):
	_tmp = mad.table['summ']
	_df = pd.DataFrame(dict(_tmp))
	_df.index = [_tmp._name]
	return _df

def tsummvar(mad, variable):
	return mad.table['summ'][variable][0]

def build_sequenceIBS(mad, beam):

    mad.input('''

        ! Specify machine version
        ver_lhc_run = 0;
        ver_hllhc_optics = 1.5;

        value,mylhcbeam;

        option, -echo, -warn, -info;
        ! Load the main sequence
        call,file="optics_repository/runIII/lhc.seq";
		
		! Install HL-LHC
        call, file="optics_repository/HLLHCV1.5/hllhc_sequence.madx";
        ! exec, disable_sext(ms.10) ! disable ms10 in the sequence

        ! Get the toolkit
        call, file="optics_repository/HLLHCV1.5/toolkit/macro.madx";

        !Cycling w.r.t. to IP3 (mandatory to find closed orbit in collision in the presence of errors)
        seqedit, sequence=lhcb1; flatten; cycle, start=IP3; flatten; endedit;
        seqedit, sequence=lhcb2; flatten; cycle, start=IP3; flatten; endedit;
        ''')

def attach_beams_to_sequencesIBS(mad, beamparams):
    mad.globals.nrj = beamparams['beam_energy_tot']
    particle_type = 'proton'

    particle_mass = mad.globals.pmass # proton mass
    particle_charge = 1.

    gamma_rel = (particle_charge*beamparams['beam_energy_tot'])/particle_mass
    for ss in mad.sequence.keys():
        # bv and bv_aux flags
        if ss == 'lhcb1':
            ss_beam_bv, ss_bv_aux = 1, 1
        elif ss == 'lhcb2':
            ss_beam_bv, ss_bv_aux = -1, 1

        mad.globals['bv_aux'] = ss_bv_aux
        print(f'>>>> adding beam for sequence {ss}')
        print(f'''
        beam, particle={particle_type},sequence={ss},
            energy={beamparams['beam_energy_tot']*particle_charge},
            sigt={beamparams['beam_sigt']},
            bv={ss_beam_bv},
            npart={beamparams['beam_npart']},
            sige={beamparams['beam_sige']},
            ex={beamparams['beam_norm_emit_x'] * 1e-6 / gamma_rel},
            ey={beamparams['beam_norm_emit_y'] * 1e-6 / gamma_rel},
            mass={particle_mass},
            charge={particle_charge};
        ''')
        mad.input(f'''
        beam, particle={particle_type},sequence={ss},
            energy={beamparams['beam_energy_tot']*particle_charge},
            sigt={beamparams['beam_sigt']},
            bv={ss_beam_bv},
            npart={beamparams['beam_npart']},
            sige={beamparams['beam_sige']},
            ex={beamparams['beam_norm_emit_x'] * 1e-6 / gamma_rel},
            ey={beamparams['beam_norm_emit_y'] * 1e-6 / gamma_rel},
            mass={particle_mass},
            charge={particle_charge};
        ''')

def apply_optics(mad, optics_file):
    mad.call(optics_file)
    if optics_file.find('thin')>=0:
    	print (f'>>>> thin detected in file name, applying slicing!')
    	mad.input('exec, myslice;')

def apply_RF(mad, seq, vtot):
    # Switch on/off RF cavities
    mad.globals['vrf400'] = vtot
    if seq == 'lhcb1':
        mad.globals['lagrf400.b1'] = 0.5
    elif seq == 'lhcb2':
        mad.globals['lagrf400.b2'] = 0.


# Make links
links = configuration['links']
for kk in links.keys():
    if os.path.exists(kk):
        os.remove(kk)
    os.symlink(os.path.abspath(links[kk]), kk)

# Create empty temp folder
os.system('rm -r temp')
os.system('mkdir temp')

#
# --- Start MAD-X here 
#
mad = Madx(command_log='hl_lhc_ibs.cmdlog')


# Select beam
mad.input('mylhcbeam =1;')

incpymad = False
if incpymad :
	build_sequenceIBS(mad, beam=1)
else:
	mad.input(f'''
	option -echo,-warn,-info;
	call,file='/afs/cern.ch/eng/lhc/optics/runIII/lhc.seq';
	call,file='/afs/cern.ch/eng/lhc/optics/HLLHCV1.5/hllhc_sequence.madx';
	call,file='/afs/cern.ch/eng/lhc/optics/HLLHCV1.5/toolkit/macro.madx';
	''')

attach_beams_to_sequencesIBS(mad, configuration['beam_params'])
apply_optics(mad, configuration['optics_file'])


if configuration['check_betas_at_ips']:
	mad.input('exec,check_ip(b1); exec,check_ip(b2);')

ring_length = tsummvar(mad, 'length')
clight = mad.globals['clight']
ring_period = ring_length/clight

apply_RF(mad, 'lhcb1', configuration['vrf_total'])

# --- initial evaluation to get teh RSynch parameters

mad.input(f'''
	use,sequence=lhcb1;
	beam,sequence=lhcb1,radiate;
	show,beam%lhcb1;
	emit,deltap=0;
	twiss,chrom;

	print,text=">>>> twiss completed Radiation Integrals: ";
	I1=table(summ,synch_1);
	I2=table(summ,synch_2);
	I3=table(summ,synch_3);
	I4=table(summ,synch_4);
	I5=table(summ,synch_5);

	value,I1,I2,I3,I4,I5;
	''')

# mad.twiss(sequence='lhcb1', chrom=True)


