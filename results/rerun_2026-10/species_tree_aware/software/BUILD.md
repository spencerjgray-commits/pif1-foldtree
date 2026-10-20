# Software builds used for results/rerun_2026-10/species_tree_aware (2026-10-16)

The build trees lived in a temporary session directory; the recipes below rebuild them. The binaries actually used are
kept, outside git, in `tools/generax-2.0.4-arm64/`, `tools/ale-b59dfaf-x86_64-patched/` and `tools/ale-b59dfaf-x86_64-unpatched/`
(they link the `strecon` / `strecon_x86` conda envs' libmpi, libomp and libc++). Build and validation logs are in
`build_logs/`; raw validation runs are in `../validation/generax204_raw/` and `../validation/ale_raw/`. Two new conda envs were created
for the toolchains (no existing env or global config was changed):
`strecon` (osx-arm64: generax=2.1.3 [bioconda, only for comparison], cmake, openmpi, cxx-compiler, make, llvm-openmp, git, bison, flex)
and `strecon_x86` (osx-64: cmake, openmpi, cxx-compiler, c-compiler, make, llvm-openmp, bison, flex).
SDK: /Library/Developer/CommandLineTools/SDKs/MacOSX15.sdk (exported as SDKROOT and CONDA_BUILD_SYSROOT).

## GeneRax 2.0.4 (arm64, native)
    git clone https://github.com/BenoitMorel/GeneRax.git generax204 && cd generax204 && git checkout 2.0.4   # e4fab40, 2021-08-03
    git submodule update --init --recursive
    (cd ext/pll-modules/libs/libpll && git apply generax-2.0.4_libpll_arm64.patch)   # no x86 cpuid / x86intrin on arm64
    mkdir build_upstream && cd build_upstream
    cmake .. -DCMAKE_POLICY_VERSION_MINIMUM=3.5 -DENABLE_SSE=False -DENABLE_AVX=False -DENABLE_AVX2=False \
      -DCMAKE_C_COMPILER=clang -DCMAKE_CXX_COMPILER=clang++ -DCMAKE_OSX_SYSROOT=$SDKROOT \
      -DBISON_EXECUTABLE=$(which bison) -DFLEX_EXECUTABLE=$(which flex)
    make -j6 generax
**Do not set CMAKE_BUILD_TYPE=Release.** GeneRax's own install.sh builds without a build type (asserts on). With
-DNDEBUG, `src/core/search/Moves.cpp:103` (`assert(PLL_SUCCESS == pllmod_utree_spr(...))`) and
`src/core/search/Rollbacks.cpp:10` (`assert(PLL_SUCCESS == pllmod_tree_rollback(...))`) are compiled out, so SPR moves
and rollbacks are never applied and the SPR worker segfaults at the first move. A Release build of 2.0.4 crashed
this way on arm64 and on x86_64 (Rosetta); the bioconda GeneRax 2.1.3 (arm64) crashed the same way at its first SPR
move (bus error) and, in EVAL mode, in `Scenario::savePerSpeciesEventsCounts`.

Validation of the build used: `--strategy EVAL` on the archived `results/reconciliation/aa3di` input reproduces the
archived 2.0.4 `pif1_eventCounts.txt` and `pif1_speciesEventCounts.txt` byte for byte; on the new 957-gene / v2fix
input, its event counts equal workflow/35 on all 1,439 species nodes (`validation/compare_generax204_vs_35_aa_deep_v2fix.txt`).
GeneRax 2.0.4 (like 2.1.3) supports one substitution model per family (`PLLTreeInfo::buildModel` asserts one
submodel; only `partitions[0]` is used), so the AA+3Di partition model cannot be run in GeneRax.

## ALE (ssolo/ALE master b59dfaf3, 2025-05-14), x86_64 under Rosetta
    git clone https://github.com/ssolo/ALE.git ale_src && cd ale_src && git apply ale-b59dfaf_undated_noT_sampler.patch
    see build_ale_x86.sh (cmake + make ALEobserve ALEml_undated in the strecon_x86 env, -DCMAKE_OSX_ARCHITECTURES=x86_64)
- arm64 builds do not work for this family: on arm64 macOS `long double` is a 64-bit double, so ALEobserve aborts
  (boost::math tgamma overflow for the (2n-5)!! normaliser, n = 957), and ALE's unscaled probabilities underflow.
- The patch is one line in `sample_undated` (src/undated.cpp:995): `if (r*uq_sum<uq_resum or no_T)` ->
  `if (r*uq_sum<uq_resum)`. Unpatched, with tau=0 (DL only, `no_T`) the sampler returns a duplication for the
  first partition whenever the speciation draw fails, so sampled reconciliations are biased toward duplications.
  The likelihood is not affected (logL with tau=0 is -1432.53 patched and unpatched on the reduced AA sample).
- tau must be 0 for a DL-only model: logL rises with tau (-1432.53 at 0 or 1e-12, -1428.22 at 1e-9, -1413.16 at 1e-6),
  i.e. even tiny transfer rates absorb gene-tree/species-tree discordance.

`build_logs/*.sh` are the verbatim original build scripts. They point at the temporary build directory used at the time, and
`build_gr204_x86.sh` is the failed Release-type x86_64 build described above. They are kept as records, not for re-use.
