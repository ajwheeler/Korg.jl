# Implementation Summary: Optional `gamma_mol_lorentz` Parameter for `Line` Type

## Overview
Successfully added an optional `gamma_mol_lorentz` parameter to the `Line` type that allows replacing `gamma_stark` and `vdW` broadening with a unified molecular Lorentz broadening width when explicitly set.

## Changes Made

### 1. **Line Type Definition** ([src/linelist.jl](src/linelist.jl))
   - Added new field: `gamma_mol_lorentz::Union{F6,Missing}`
   - Updated struct signature from 6 to 7 type parameters
   - Field defaults to `missing` when not explicitly provided

### 2. **Line Constructors** ([src/linelist.jl](src/linelist.jl))
   - **Primary constructor**: Added keyword-only parameter `gamma_mol_lorentz::Union{F6,Missing}=missing`
   - **Copy constructor**: Updated to support modifying `gamma_mol_lorentz` via keyword arguments
   - Both preserve backward compatibility—existing code continues to work unchanged

### 3. **Line Absorption Logic** ([src/line_absorption.jl](src/line_absorption.jl))
   - **Conditional branching**: When `!ismissing(line.gamma_mol_lorentz)`:
     - Always starts with `gamma_rad`
     - Adds `gamma_mol_lorentz` instead of `gamma_stark` + `vdW` broadening
   - **Default behavior** (when `gamma_mol_lorentz` is `missing`):
     - Uses the original behavior: `gamma_rad` + `gamma_stark` + `vdW` broadening
   - vdW contribution is handled conditionally for non-molecular species

### 4. **HDF5 Serialization** ([src/linelist.jl](src/linelist.jl))
   - **save_linelist**: Stores `gamma_mol_lorentz` as NaN for missing values
   - **read_korg_linelist**: Reconstructs missing values from NaN, uses `map` to properly handle keyword arguments
   - Backward compatible: reads files without the field by filling with `missing`

### 5. **Tests** ([test/linelist.jl](test/linelist.jl))
   - Test: Line created without `gamma_mol_lorentz` has `missing` value
   - Test: Line created with `gamma_mol_lorentz` preserves the value
   - Test: Copy constructor preserves and allows modifying `gamma_mol_lorentz`
   - Test: HDF5 roundtrip preserves the field correctly

## Usage Examples

### Without `gamma_mol_lorentz` (default behavior unchanged)
```julia
line = Korg.Line(5000e-8, -1.0, species"Fe I", 0.3)
# Internally: gamma_mol_lorentz = missing
# In absorption: uses gamma_rad + gamma_stark + vdW
```

### With `gamma_mol_lorentz` (new feature)
```julia
line = Korg.Line(5000e-8, -1.0, species"H2O", 0.3; gamma_mol_lorentz=1e-8)
# Internally: gamma_mol_lorentz = 1e-8
# In absorption: uses gamma_rad + gamma_mol_lorentz (replaces gamma_stark + vdW)
```

### Modifying via copy constructor
```julia
line2 = Korg.Line(line; gamma_mol_lorentz=2e-8)
```

## Backward Compatibility

✓ **Fully backward compatible**
- Existing code that doesn't use `gamma_mol_lorentz` continues to work unchanged
- All Line constructors have `gamma_mol_lorentz` default to `missing`
- Line absorption defaults to original behavior when `gamma_mol_lorentz` is `missing`
- HDF5 files without the field are still readable (missing values are filled with `missing`)

## Design Rationale

The implementation uses a keyword-only parameter with `missing` as the default sentinel value because:
1. **Optional activation**: Only affects behavior when explicitly set
2. **Clear intent**: `missing` unambiguously indicates "use default behavior"
3. **Type-safe**: Works with Julia's type system and broadcasting
4. **Serializable**: `missing` can be converted to NaN in HDF5 and back
5. **Non-intrusive**: Doesn't require changes to any calling code

## Testing

Run validation with:
```bash
julia --project=. test_gamma_mol_lorentz.jl
```

All tests pass, confirming:
- Field initialization and access
- Copy constructor preservation and modification
- HDF5 roundtrip serialization
- Correct integration with line absorption logic
