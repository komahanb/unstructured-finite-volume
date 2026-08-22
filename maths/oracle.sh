#!/bin/bash
# Restore the verification suites the report cites from revision
# ddcf592 (the last commit before test/ was deleted), apply the three
# mechanical API changes of August 2026, build them against the
# working tree's library, and run them. Nothing is written outside
# maths/oracle/.
#
#   usage: bash maths/oracle.sh            (from the repository root)
#
# The suites' own run.sh expect lib/, src/ and build.sh two levels
# up, so maths/oracle/ links those to the repository.
set -e

root="$(cd "$(dirname "$0")/.." && pwd)"
base=ddcf592
suites="graph-binary graph-field graph-ordinary graph-dense-direct graph-minimization \
        graph-differentiation graph-marching graph-multigrid graph-partition"

oracle="$root/maths/oracle"
rm -rf "$oracle"
mkdir -p "$oracle"
cd "$root"
git archive "$base" $(for s in $suites; do echo test/$s; done) | tar -x -C "$oracle"
ln -sfn "$root/lib" "$oracle/lib"
ln -sfn "$root/src" "$oracle/src"
ln -sfn "$root/Makefile.in" "$oracle/Makefile.in"
printf '#!/bin/bash\ncd %s && ./build.sh\n' "$root" > "$oracle/build.sh"
chmod +x "$oracle/build.sh"
cp "$root/check_naming.sh" "$oracle/"

#---------------------------------------------------------------------
# The API changes since $base, applied to the test sources:
#   x % pattern                     ->  call x % dependencies(p)
#   call x % weights % real_vector  ->  w = x % weights   (constants likewise)
#   dual_by_basis                   ->  the compiled transpose
#---------------------------------------------------------------------
python3 - "$oracle/test" <<'EOF'
import re, sys, glob
root = sys.argv[1]

# graph-dense-direct: derived pattern, array weights/constants
for p in glob.glob(f'{root}/graph-dense-direct/*.f90'):
    s = open(p).read(); o = s
    names = set(re.findall(r'(\w+) % pattern', s))
    for n in names:
        s = s.replace(f'{n} % pattern', f'{n}_pattern')
        s = re.sub(r'(\n(\s*)call \w+ % attach\(' + n + r',)',
                   r'\n\2call ' + n + r' % dependencies(' + n + r'_pattern)\1', s, count=1)
        s = re.sub(r'(\n(\s*)type\(stencil\)\s*::\s*[^\n]*\b' + n + r'\b[^\n]*)',
                   r'\1\n\2class(directed_graph), allocatable :: ' + n + r'_pattern', s, count=1)
    if names and 'use view_directed' not in s.split('implicit none')[0]:
        s = s.replace('  implicit none', '  use view_directed, only : directed_graph\n  implicit none', 1)
    s = re.sub(r'call (\w+) % weights % real_vector\((\w+)\)', r'\2 = \1 % weights', s)
    s = re.sub(r'call (\w+) % constants % real_vector\((\w+)\)', r'\2 = \1 % constants', s)
    # the affine stencil is compiled without an attach: derive its pattern first
    s = s.replace("  compiled = stencil(affine, affine_pattern, 2, 'compiled affine')\n",
                  "  call affine % dependencies(affine_pattern)\n  compiled = stencil(affine, affine_pattern, 2, 'compiled affine')\n")
    # the contained column_of derives its own
    # the contained column_of derives its own pattern, and builds the
    # basis on the stencil's declared column domain when it has one:
    # a compiled stencil refuses inputs on any other set, by identity
    s = s.replace("""    real(dp) :: e(2)

    e    = 0.0_dp
    e(j) = 1.0_dp
    basis = stored_field('e', statement_pattern % vertex_set(), 2)
""", """    real(dp) :: e(2)
    class(directed_graph), allocatable :: statement_pattern

    call statement % dependencies(statement_pattern)
    e    = 0.0_dp
    e(j) = 1.0_dp
    if (statement % column_entries > 0) then
       basis = stored_field('e', statement % column_domain, 2)
    else
       basis = stored_field('e', statement_pattern % vertex_set(), 2)
    end if
""")
    if s != o: open(p, 'w').write(s)

# graph-differentiation: the dual by the compiled transpose
p = f'{root}/graph-differentiation/test.f90'
s = open(p).read()
s = s.replace("use operation_linearization, only : linearization, tangent_of, dual_by_basis",
              "use operation_linearization, only : linearization, tangent_of")
s = s.replace("""    call dual_by_basis(tangent_q, three, lambda, g)
    rhs = dot_product([1.0_dp, -2.0_dp, 0.5_dp], g)
    compiled = stencil(tangent_q, three, 3)
    adjoint  = compiled % transpose()
    call adjoint % apply(three, [lf], output)
    call output % real_vector(gt)
    call report(near(lhs, rhs, 1.0e-12_dp) .and. maxval(abs(g - 2.0_dp * q3 * lambda)) < 1.0e-12_dp &
         & .and. maxval(abs(g - gt)) < 1.0e-12_dp, &
         & "state block: <J v, lambda> = <v, J^T lambda>, J^T lambda = 2 q lambda, &
         &and the dual by basis equals the compiled transpose", nfail)
""", """    compiled = stencil(tangent_q, three, 3)
    adjoint  = compiled % transpose()
    call adjoint % apply(three, [lf], output)
    call output % real_vector(g)
    rhs = dot_product([1.0_dp, -2.0_dp, 0.5_dp], g)
    call report(near(lhs, rhs, 1.0e-12_dp) .and. maxval(abs(g - 2.0_dp * q3 * lambda)) < 1.0e-12_dp, &
         & "state block: <J v, lambda> = <v, J^T lambda>, J^T lambda = 2 q lambda, &
         &through the compiled transpose", nfail)
""")
s = s.replace("""    call dual_by_basis(tangent_xi, three, lambda, g)
    rhs = 0.75_dp * g(1)
""", """    compiled = stencil(tangent_xi, three, 1, column_domain=xif % domain(), column_entries=1)
    adjoint  = compiled % transpose()
    call adjoint % apply(three, [lf], output)
    call output % real_vector(g)
    rhs = 0.75_dp * g(1)
""")
open(p, 'w').write(s)
EOF

#---------------------------------------------------------------------
# Build the library once, then every suite; run each and count.
#---------------------------------------------------------------------
( cd "$root" && ./build.sh > /dev/null )

total_pass=0; total_fail=0; status=0
for s in $suites; do
    d="$oracle/test/$s"
    if ! make -C "$d" > "$d/build.log" 2>&1; then
        echo " FAIL : $s does not build (see $d/build.log)"; status=1; continue
    fi
    out=$( cd "$d" && ./run 2>&1 ) || true
    p=$(echo "$out" | grep -c ' PASS' || true); f=$(echo "$out" | grep -c ' FAIL' || true)
    total_pass=$((total_pass + p)); total_fail=$((total_fail + f))
    printf ' %-24s PASS %3d  FAIL %2d\n' "$s" "$p" "$f"
    [ "$f" -eq 0 ] || status=1
done

# the refusal cases of graph-differentiation: each must stop, loudly
d="$oracle/test/graph-differentiation"
if [ -x "$d/refusal" ]; then
    r=$( cd "$d" && sed -n '/declare -A/,$p' run.sh > refusals.sh && bash refusals.sh 2>&1 | grep -c ' PASS' || true )
    printf ' %-24s PASS %3d\n' "refusals" "$r"
fi

echo " TOTAL : PASS $total_pass  FAIL $total_fail"
exit $status
