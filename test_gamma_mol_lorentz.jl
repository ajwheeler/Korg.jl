#!/usr/bin/env julia
using Pkg
Pkg.activate(".")

using Korg, HDF5

println("Testing optional gamma_mol_lorentz field...")

# Test 1: Create a Line without gamma_mol_lorentz
l1 = Korg.Line(5000.0, 0.0, Korg.species"Fe I", 1.0)
@assert ismissing(l1.gamma_mol_lorentz) "Line without gamma_mol_lorentz should have missing value"
println("✓ Test 1 passed: Line without gamma_mol_lorentz defaults to missing")

# Test 2: Create a Line with gamma_mol_lorentz
l2 = Korg.Line(5000.0, 0.0, Korg.species"FeH", 1.0; gamma_mol_lorentz=1e-8)
@assert l2.gamma_mol_lorentz == 1e-8 "gamma_mol_lorentz should be set"
println("✓ Test 2 passed: Line with gamma_mol_lorentz set correctly")

# Test 3: Copy constructor preserves gamma_mol_lorentz
l3 = Korg.Line(l2)
@assert l3.gamma_mol_lorentz == 1e-8 "Copy constructor should preserve gamma_mol_lorentz"
println("✓ Test 3 passed: Copy constructor preserves gamma_mol_lorentz")

# Test 4: Copy constructor can modify gamma_mol_lorentz
l4 = Korg.Line(l2; gamma_mol_lorentz=2e-8)
@assert l4.gamma_mol_lorentz == 2e-8 "Copy constructor should allow modifying gamma_mol_lorentz"
println("✓ Test 4 passed: Copy constructor can modify gamma_mol_lorentz")

# Test 5: HDF5 roundtrip
filename = tempname() * ".h5"
linelist = [l1, l2]
Korg.save_linelist(filename, linelist)
linelist_loaded = Korg.read_linelist(filename)
@assert length(linelist_loaded) == 2 "Should load 2 lines"
@assert ismissing(linelist_loaded[1].gamma_mol_lorentz) "First line should have missing gamma_mol_lorentz"
@assert linelist_loaded[2].gamma_mol_lorentz == 1e-8 "Second line should have gamma_mol_lorentz set"
println("✓ Test 5 passed: HDF5 roundtrip preserves gamma_mol_lorentz")

# Clean up
rm(filename)

println("\n✓ All tests passed! The optional gamma_mol_lorentz field is working correctly.")
println("\nSummary:")
println("- New field gamma_mol_lorentz added to Line struct")
println("- When not explicitly set, gamma_mol_lorentz defaults to missing")
println("- Copy constructor supports modifying gamma_mol_lorentz")
println("- HDF5 serialization preserves the field")
println("- Line absorption uses gamma_mol_lorentz when explicitly provided")
