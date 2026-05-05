using Korg, Test

# Test 1: Simple tuple creation
println("Test 1: Creating line with tuple parameters...")
line1 = Korg.Line(5000.0, 0.5, Korg.Species("H"), 1.0; gamma_mol_lorentz=(1e4, 2e4), n_exp=(0.5, 0.6))
println("✓ Line created: gamma_mol_lorentz = ", line1.gamma_mol_lorentz, ", n_exp = ", line1.n_exp)

# Test 2: Copy constructor with mode switching
println("\nTest 2: Copy constructor with mode switching...")
line2 = Korg.Line(5000.0, 0.0, Korg.Species("Fe"), 1.0)  # Mode 1
println("  Original line: gamma_mol_lorentz = ", line2.gamma_mol_lorentz)
line3 = Korg.Line(line2; gamma_mol_lorentz=(1e-8, 2e-8), n_exp=(0.5, 0.6))  # Switch to Mode 2
println("✓ After copy with mode switching: gamma_mol_lorentz = ", line3.gamma_mol_lorentz, ", n_exp = ", line3.n_exp)

# Test 3: HDF5 roundtrip
println("\nTest 3: HDF5 roundtrip...")
filename = tempname() * ".h5"
Korg.save_linelist(filename, [line1])
loaded_lines = Korg.read_linelist(filename)
println("✓ Loaded line: gamma_mol_lorentz = ", loaded_lines[1].gamma_mol_lorentz, ", n_exp = ", loaded_lines[1].n_exp)
@test loaded_lines[1].gamma_mol_lorentz == (1e4, 2e4)
@test loaded_lines[1].n_exp == (0.5, 0.6)
rm(filename)

println("\nAll tests passed!")
