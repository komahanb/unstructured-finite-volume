#!/bin/bash
set -euo pipefail
here="$(cd "$(dirname "$0")" && pwd)"
root="$(cd "$here/../.." && pwd)"
if [ "${UFVM_SKIP_LIBRARY_BUILD:-0}" != 1 ]; then
    (cd "$root" && ./build.sh && ./application/build.sh)
fi
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT
module_dir=${UFVM_GTI_MODULE_DIR:-$work}
compiler=${F90:-gfortran-15}
if [ ! -f "$module_dir/gti_demos.mod" ]; then
    # verify.sh supplies a fresh private directory shared by the GTI suites.
    sed '/^program graph_time_integrator/,$d' "$root/application/module_graph_time_integrator.f90" > "$module_dir/modules.f90"
    "$compiler" -std=f2023 -fcoarray=single -cpp -fbounds-check -O2 -I"$root/lib" -J"$module_dir" \
        -c "$module_dir/modules.f90" -o "$module_dir/modules.o"
fi
"$compiler" -std=f2023 -fcoarray=single -fbounds-check -O2 -I"$root/lib" -I"$module_dir" -J"$work" \
    "$here/test.f90" "$module_dir/modules.o" "$root/lib/libufvm.a" -o "$work/run"
"$work/run" accuracy
"$work/run" status
"$work/run" grid_stationary
"$work/run" functional_error
for mode in adaptation_unmet adaptive_failure minimum_step forward reverse linear_forward linear_reverse; do
    case "$mode" in
        adaptation_unmet) expected='ADAPTATION_UNMET' ;;
        adaptive_failure) expected='step solve did not converge' ;;
        minimum_step) expected='minimum step cannot meet the tolerance' ;;
        forward|reverse) expected='primal march must converge' ;;
        linear_forward|linear_reverse) expected='derivative solve did not converge' ;;
    esac
    if "$work/run" "$mode" > "$work/$mode.log" 2>&1; then
        echo "FAIL: $mode was accepted"
        exit 1
    fi
    if ! grep -q "$expected" "$work/$mode.log"; then
        cat "$work/$mode.log"
        exit 1
    fi
    echo "PASS: $mode reports its failure"
done
# The application discovers an adaptive grid for the one configured family at
# max_discretization_order, and refuses the superseded word, two families, and step
# doubling of a multistep family, which has no step from one state.
application="$root/application"
common='--grid=adaptive --time_duration=2 --design=0 --max_derivative_degree=0 --max_discretization_order=2'
for case in dirk_stationary bdf_stationary dirk_doubling superseded_word two_families multistep_doubling \
            dirk_functional_error budget_unmet; do
    families='--families=dirk'
    case "$case" in
        dirk_stationary)    arguments='--adaptive_check=grid_stationarity'
                            expected='8 steps of dirk2 at grid-stationarity tolerance  1.00E-12 (0 rejected)'; accepted=1 ;;
        # BDF-2's F_h is not stationary on a uniform grid: e/|F_h| = 1.82e-3 on the seed and
        # 1.35e-3 after one halving round, so 1.5e-3 requires exactly one rejection: the halving,
        # the re-march through the startup block and the re-differentiation all run once
        bdf_stationary)     arguments='--adaptive_check=grid_stationarity --grid_stationarity_tolerance=1.5e-3'
                            families='--families=bdf'
                            expected='12 steps of bdf2 at grid-stationarity tolerance  1.50E-03 (1 rejected)'; accepted=1 ;;
        dirk_doubling)      arguments='--adaptive_check=step_doubling --tolerance=1e-6'
                            expected='steps of dirk2 at estimated local-error tolerance  1.00E-06'; accepted=1 ;;
        superseded_word)    arguments='--adaptive_check=goal_oriented'
                            expected='adaptive_check = grid_stationarity'; accepted=0 ;;
        two_families)       arguments='--adaptive_check=grid_stationarity'; families='--families=bdf dirk'
                            expected='discovered for one family'; accepted=0 ;;
        multistep_doubling) arguments='--adaptive_check=step_doubling'; families='--families=bdf'
                            expected='self-starting scheme'; accepted=0 ;;
        # the estimator-driven grid from the configured 21-instant seed: implicit midpoint's
        # E - E_h = 2.5e-3 at h = 0.1 is above 1e-3, every step is halved once, and the
        # 40-step grid has E - E_h = 6.2e-4 within the tolerance
        dirk_functional_error) arguments='--adaptive_check=functional_error --functional_error_tolerance=1e-3 --instants=21'
                            expected='40 steps of dirk2 at functional-error tolerance  1.00E-03 (1 rejected)'; accepted=1 ;;
        # the same loop with adaptation_instants below the 41 instants it needs: no
        # acceptance, the outcome is reported and the program stops
        budget_unmet)       arguments='--adaptive_check=functional_error --functional_error_tolerance=1e-3 --instants=21 --adaptation_instants=30'
                            expected='ADAPTATION_UNMET'; accepted=0 ;;
    esac
    if (cd "$application" && ./graph_time_integrator $common $arguments "$families") > "$work/$case.log" 2>&1; then
        status=1
    else
        status=0
    fi
    if [ "$status" != "$accepted" ]; then
        cat "$work/$case.log"
        echo "FAIL: $case exit status"
        exit 1
    fi
    if ! grep -qF "$expected" "$work/$case.log"; then
        cat "$work/$case.log"
        echo "FAIL: $case does not report: $expected"
        exit 1
    fi
    echo "PASS: $case"
done
# The names a law brings with it are refused where the law does not admit them:
# the radial oscillator is of second order in time, van der Pol admits no square
# integral, and the radial oscillator admits no dissipation.
radial='--physics=radial_oscillator --design=2.0 --time_duration=2 --instants=11 --max_derivative_degree=0'
for case in radial_degree unadmitted_functional unadmitted_dissipation mode_without_eigenfunction; do
    case "$case" in
        radial_degree)         arguments="$radial --state_degree=1 --functionals=energy"
                               expected='the radial oscillator must be of second order in time' ;;
        unadmitted_functional) arguments='--physics=vanderpol --time_duration=2 --instants=11 --functionals=square_integral'
                               expected='vanderpol admits no functional named square_integral' ;;
        unadmitted_dissipation) arguments="$radial --functionals=dissipation"
                               expected='radial_oscillator admits no functional named dissipation' ;;
        mode_without_eigenfunction)
                               arguments='--config=disc --spatial_geometry=elliptical --initial_field=mode --check=none'
                               expected='requires the box or the disc (the ellipse has no eigenfunction of the laplacian)' ;;
    esac
    if (cd "$application" && ./graph_time_integrator $arguments) > "$work/$case.log" 2>&1; then
        cat "$work/$case.log"
        echo "FAIL: $case was accepted"
        exit 1
    fi
    if ! grep -qF "$expected" "$work/$case.log"; then
        cat "$work/$case.log"
        echo "FAIL: $case does not report: $expected"
        exit 1
    fi
    echo "PASS: $case"
done
