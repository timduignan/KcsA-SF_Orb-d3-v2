import math

# Given times in nanoseconds
times_ns = [0.19, 0.43, 0.72, 0.25, 0.05, 0.87]
n = len(times_ns)

# Calculate mean time and standard error
mean_time_ns = sum(times_ns) / n
variance_time = sum((x - mean_time_ns) ** 2 for x in times_ns) / (n - 1)
std_dev_time = math.sqrt(variance_time)
std_error_time = std_dev_time / math.sqrt(n)

# Convert mean time and standard error to seconds
mean_time_s = mean_time_ns * 1e-9
std_error_time_s = std_error_time * 1e-9

# Constants
q = 1.602176634e-19  # Elementary charge in Coulombs

# Compute current using mean time
I = q / mean_time_s
# Compute error in current
delta_I = I * (std_error_time_s / mean_time_s)

# Given force and membrane width
F_eV_per_A = 0.1  # Force in eV/Å
d_A = 50  # Membrane width in Ångströms

# Compute voltage
V = F_eV_per_A * d_A  # Voltage in Volts (since 1 eV = 1 V for electrons)

# Compute conductance and its error
G = I / V
delta_G = delta_I / V

# Convert conductance to picoSiemens for readability
G_pS = G * 1e12  # Convert S to pS
delta_G_pS = delta_G * 1e12  # Convert S to pS

print(f"Mean time (ns): {mean_time_ns:.4f} ± {std_error_time:.4f}")
print(f"Mean current (A): {I:.4e} ± {delta_I:.4e}")
print(f"Voltage across the membrane (V): {V}")
print(f"Mean Conductance (pS): {G_pS:.2f} ± {delta_G_pS:.2f}")
