#!/usr/bin/env julia
using Pkg
Pkg.activate(".")

using Korg, HDF5

println("Testing dual-mode Line broadening implementation...")
println()

# Test 1: Create a Line in Mode 1 (default, backward compatible)
println("Test 1: Mode 1 - Default Stark + vdW broadening")
l1 = Korg.Line(5000.0, 0.0, Korg.species"Fe I", 1.0)
@assert !ismissing(l1.gamma_stark) "Mode 1: gamma_stark should be approximated"
@assert !ismissing(l1.vdW) "Mode 1: vdW should be approximated"
@assert ismissing(l1.gamma_mol_lorentz) "Mode 1: gamma_mol_lorentz should be missing"
@assert ismissing(l1.n_exp) "Mode 1: n_exp should be missing"
println("✓ Line created in Mode 1 (Stark + vdW)")
println("  gamma_stark = $(l1.gamma_stark)")
println("  vdW = $(l1.vdW)")
println()

# Test 2: Create a Line in Mode 2 (molecular Lorentz with temperature dependence)
println("Test 2: Mode 2 - Molecular Lorentz with temperature exponent")
l2 = Korg.Line(5000.0, 0.0, Korg.species"H2O", 1.0; gamma_mol_lorentz=1e-8, n_exp=0.5)
@assert ismissing(l2.gamma_stark) "Mode 2: gamma_stark should be missing"
@assert ismissing(l2.vdW) "Mode 2: vdW should be missing"
@assert l2.gamma_mol_lorentz == 1e-8 "Mode 2: gamma_mol_lorentz should be set"
@assert l2.n_exp == 0.5 "Mode 2: n_exp should be set to 0.5"
println("✓ Line created in Mode 2 (gamma_mol_lorentz + n_exp)")
println("  gamma_mol_lorentz = $(l2.gamma_mol_lorentz)")
println("  n_exp = $(l2.n_exp)")
println()

# Test 3: Error when mixing modes (providing both Mode 1 and Mode 2 parameters)
println("Test 3: Error handling when mixing modes")
try
    l_bad = Korg.Line(5000.0, 0.0, Korg.species"H2O", 1.0, missing, 1e-7, missing; 
                      gamma_mol_lorentz=1e-8, n_exp=0.5)
    error("Should have raised an error when mixing modes!")
catch e
    if contains(string(e), "Cannot provide both")
        println("✓ Correctly rejected mixed-mode parameters")
        println("  Error message: $(sprint(showerror, e))")
    else
        rethrow()
    end
end
println()

# Test 4: Error when providing n_exp without gamma_mol_lorentz
println("Test 4: Error when providing n_exp without gamma_mol_lorentz")
try
    l_bad = Korg.Line(5000.0, 0.0, Korg.species"H2O", 1.0; n_exp=0.5)
    error("Should have raised an error!")
catch e
    if contains(string(e), "When providing n_exp, gamma_mol_lorentz must also be provided")
        println("✓ Correctly rejected n_exp without gamma_mol_lorentz")
        println("  Error message: $(sprint(showerror, e))")
    else
        rethrow()
    end
end
println()

# Test 5: Copy constructor in Mode 1
println("Test 5: Copy constructor - Mode 1")
l1_copy = Korg.Line(l1)
@assert !ismissing(l1_copy.gamma_stark) "Copy should preserve gamma_stark"
@assert !ismissing(l1_copy.vdW) "Copy should preserve vdW"
@assert ismissing(l1_copy.gamma_mol_lorentz) "Copy should preserve missing gamma_mol_lorentz"
println("✓ Mode 1 copy constructor works correctly")
println()

# Test 6: Copy constructor in Mode 2
println("Test 6: Copy constructor - Mode 2")
l2_copy = Korg.Line(l2)
@assert ismissing(l2_copy.gamma_stark) "Copy should preserve missing gamma_stark"
@assert ismissing(l2_copy.vdW) "Copy should preserve missing vdW"
@assert l2_copy.gamma_mol_lorentz == 1e-8 "Copy should preserve gamma_mol_lorentz"
@assert l2_copy.n_exp == 0.5 "Copy should preserve n_exp"
println("✓ Mode 2 copy constructor works correctly")
println()

# Test 7: Modify Mode 1 line
println("Test 7: Modifying Mode 1 line via copy constructor")
l1_mod = Korg.Line(l1; gamma_stark=2e-7)
@assert l1_mod.gamma_stark == 2e-7 "Should be able to modify gamma_stark"
@assert ismissing(l1_mod.gamma_mol_lorentz) "gamma_mol_lorentz should remain missing"
println("✓ Can modify gamma_stark in Mode 1")
println()

# Test 8: Modify Mode 2 line
println("Test 8: Modifying Mode 2 line via copy constructor")
l2_mod = Korg.Line(l2; n_exp=0.7)
@assert ismissing(l2_mod.gamma_stark) "gamma_stark should remain missing"
@assert l2_mod.n_exp == 0.7 "n_exp should be modified"
println("✓ Can modify n_exp in Mode 2")
println()

# Test 9: HDF5 roundtrip with both modes
println("Test 9: HDF5 serialization and deserialization")
filename = tempname() * ".h5"
linelist = [l1, l2, l1_mod, l2_mod]
Korg.save_linelist(filename, linelist)
println("  Saved $(length(linelist)) lines to HDF5")

linelist_loaded = Korg.read_linelist(filename)
@assert length(linelist_loaded) == 4 "Should load 4 lines"
println("  Loaded $(length(linelist_loaded)) lines from HDF5")

# Check Mode 1 line
l1_loaded = linelist_loaded[1]
@assert !ismissing(l1_loaded.gamma_stark) "Loaded Mode 1 line should have gamma_stark"
@assert !ismissing(l1_loaded.vdW) "Loaded Mode 1 line should have vdW"
@assert ismissing(l1_loaded.gamma_mol_lorentz) "Loaded Mode 1 line should have missing gamma_mol_lorentz"
println("  ✓ Mode 1 line correctly loaded")

# Check Mode 2 line
l2_loaded = linelist_loaded[2]
@assert ismissing(l2_loaded.gamma_stark) "Loaded Mode 2 line should have missing gamma_stark"
@assert ismissing(l2_loaded.vdW) "Loaded Mode 2 line should have missing vdW"
@assert l2_loaded.gamma_mol_lorentz == 1e-8 "Loaded Mode 2 line should have gamma_mol_lorentz"
@assert l2_loaded.n_exp == 0.5 "Loaded Mode 2 line should have n_exp"
println("  ✓ Mode 2 line correctly loaded")

rm(filename)
println("✓ HDF5 roundtrip successful")
println()

# Test 10: Check that the implementation is backward compatible
println("Test 10: Backward compatibility check")
# Create a line the old way and verify it still works
l_old = Korg.Line(5000.0, 0.0, Korg.species"Fe I", 1.0)
@assert !ismissing(l_old.gamma_stark)
@assert !ismissing(l_old.vdW)
println("✓ Old-style line creation still works")
println()

println(repeat("=", 60))
println("✓ ALL TESTS PASSED!")
println(repeat("=", 60))
println()
println("Summary:")
println("========")
println("• Mode 1 (default): Uses gamma_stark and vdW broadening")
println("  - Backward compatible")
println("  - Parameters auto-approximated if not provided")
println()
println("• Mode 2 (new): Uses gamma_mol_lorentz with temperature exponent n_exp")
println("  - Enables custom molecular broadening with temperature dependence")
println("  - Scaling: γ = gamma_mol_lorentz * (T/T_ref)^n_exp")
println("  - Both parameters required when using Mode 2")
println()
println("• Mode mixing is prevented: Cannot provide both Mode 1 and Mode 2 parameters")
println("• HDF5 serialization preserves both modes")
println("• Copy constructor supports both modes")
