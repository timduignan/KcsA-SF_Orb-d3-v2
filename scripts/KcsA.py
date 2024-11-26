import torch
torch.set_float32_matmul_precision('high')  # Set high precision for float32 matrix multiplication
from ase.io import read
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.io.trajectory import Trajectory
from ase.md import MDLogger
from ase import units
from orb_models.forcefield import atomic_system
from orb_models.forcefield import pretrained
from ase.calculators.calculator import Calculator, all_properties
import numpy as np
import os
from ase.constraints import FixAtoms
from ase.md.langevin import Langevin  # Base Langevin thermostat for NVT



# Define input file name and cell size at the top
input_file = "kcsa-mod-filtered-cylinder-4K.xyz"
cell_size = 72.950196  # in Angstroms

# Generate output file name
output_file = os.path.splitext(input_file)[0] + ".traj"


            
# Define a custom Orb-d3-v1 calculator with CUDA support
class OrbD3Calculator(Calculator):
    implemented_properties = ['energy', 'forces']

    def __init__(self, model, additional_forces=None, **kwargs):
        super().__init__(**kwargs)
        self.model = model
        self.additional_forces = additional_forces

        # Check if CUDA is available and move the model to GPU if possible
        if torch.cuda.is_available():
            self.model.cuda()
        else:
            print("CUDA is not available, running on CPU")

    def calculate(self, atoms=None, properties=['energy'], system_changes=all_properties):
        # Ensure the atoms object is set correctly
        super().calculate(atoms, properties, system_changes)

        # Convert the ASE atoms to a graph representation
        graph = atomic_system.ase_atoms_to_atom_graphs(atoms)

        # Move the graph to GPU if CUDA is available
        if torch.cuda.is_available():
            graph = graph.to('cuda')

        # Perform the prediction using orb-d3-v1
        result = self.model.predict(graph)

        # Store results in self.results dictionary for ASE compatibility
        self.results['energy'] = result["graph_pred"].detach().cpu().numpy().sum()
        forces = result["node_pred"].detach().cpu().numpy()

        # Add additional forces if provided
        if self.additional_forces is not None:
            forces += self.additional_forces

        self.results['forces'] = forces

# Load the checkpoint manually
local_model_path = "/scratch/tn51/ttd110/orb/orb-d3-v2-20241011.ckpt"

# Load the model using the local checkpoint path
orbff = pretrained.orb_d3_v2(weights_path=local_model_path)

# Load XYZ file into ASE Atoms object
atoms = read(input_file)

# Update cell size
atoms.set_cell([cell_size, cell_size, cell_size])
atoms.set_pbc([True, True, True])  # Apply periodic boundary conditions


# --- Add the following code to identify potassium ions and define the force ---

# Identify the indices of potassium ions
potassium_indices = [atom.index for atom in atoms if atom.symbol == 'K']



# Define the external force magnitude in eV/Angstrom (ASE's default unit for force)
F_magnitude = 0.05  # Replace with the calculated force magnitude
force_vector = np.array([0.0, 0.0, F_magnitude])

# Create per-atom forces array
num_atoms = len(atoms)
per_atom_forces = np.zeros((num_atoms, 3))
for i in potassium_indices:
    per_atom_forces[i] += force_vector  # Add force to potassium ions

# Attach the custom Orb-d3-v1 calculator to the atoms object with additional forces
atoms.calc = OrbD3Calculator(model=orbff, additional_forces=per_atom_forces)


# Open and read the file containing indices of atoms to fix
with open('outerind.dat', 'r') as f:
    outind = [int(x) for x in f.read().split()]

# Apply constraints to fix specified atoms
constraint = FixAtoms(indices=outind)
atoms.set_constraint(constraint)


# Set up the initial velocities corresponding to a temperature of 300 K
temperature_K = 300  # Temperature in Kelvin
#commenting out cause its a restart
MaxwellBoltzmannDistribution(atoms, temperature_K=temperature_K)

# Zero velocities of the fixed atoms
velocities = atoms.get_velocities()
velocities[outind] = 0.0
atoms.set_velocities(velocities)



# Determine the last step from the restart file
# Assuming the restart file is a trajectory file, you can use the Trajectory object to read it
#restart_traj = Trajectory(input_file)
#last_step = len(restart_traj)  # Get the number of frames in the trajectory


# Define the Langevin thermostat for NVT
friction = 0.01 / units.fs  # Friction coefficient in inverse time units (1/fs)
timestep = 0.5 * units.fs

# Initialize the dynamics with the correct starting step
dyn = Langevin(atoms, timestep, temperature_K=temperature_K, friction=friction)

# Function to wrap atoms back into the simulation box (if needed)
def wrap_atoms(atoms=atoms):
    atoms.wrap()

# Function to write XYZ files with custom comment line
def write_xyz(atoms, step, F_magnitude, potassium_count, cell_size, timestep):
    xyz_file = f"KcsA_{potassium_count}K_{F_magnitude:.1f}F.xyz"
    
    time_fs = step * timestep / units.fs  # Convert time to femtoseconds
    comment = f"KcsA with {potassium_count} potassiums, force: {F_magnitude} eV/Å, cell size: {cell_size} Å, time: {time_fs:.2f} fs"
    with open(xyz_file, 'a') as f:  # Open the file in append mode
        f.write(f"{len(atoms)}\n")
        f.write(f"{comment}\n")
        for atom in atoms:
            pos = atom.position
            f.write(f"{atom.symbol} {pos[0]} {pos[1]} {pos[2]}\n")

# Attach the XYZ writer to the dynamics
potassium_count = len(potassium_indices)
xyz_interval = 20  # Save XYZ every 20 steps
dyn.attach(lambda: write_xyz(atoms, dyn.get_number_of_steps(), F_magnitude, potassium_count,cell_size,timestep), interval=xyz_interval)

# Save the trajectory to a .traj file for restart purposes
restart_trajectory = Trajectory("restart.traj", "w", atoms)

# Function to write restart file
def write_restart(atoms):
    restart_trajectory.write(atoms)

# Attach the restart writer to the dynamics
restart_interval = 10000  # Save restart file every 10,000 steps
dyn.attach(lambda: write_restart(atoms), interval=restart_interval)

# Optional: Add an MDLogger to print energy and forces
dyn.attach(MDLogger(dyn, atoms, "md_nvt.log", header=True, stress=False, peratom=True), interval=20)

# Run the simulation for 1,000,000 time steps
dyn.run(250000)
