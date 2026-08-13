#!/bin/bash
# The adjoint tower's RESULT CONTRACT, in one place: what an
# ADJOINT_RESULT marker must look like, and a self-test proving the
# check accepts genuine real values and refuses malformed ones.
#
# The contract is shape and syntax only:
#
#     exactly one marker line
#     exactly one token - this tower's answer is a scalar sensitivity
#     that token a finite real numeric literal
#
# and NOTHING about what the number should be. df/dp is a real
# quantity and is serialized as it stands: the computed value prints
# as 6.9999999999999991E+00, which is the honest image of the
# arithmetic. Whether that IS the mathematical 7 is the Level-9
# test's business, not this checker's - a checker that demanded the
# integer 7 would be demanding a rounded answer.

# marker_ok <marker-count> <expected-width> <tokens>
marker_ok() {
    local marks="$1" width="$2" values="$3"
    local real='^[+-]?([0-9]+(\.[0-9]*)?|\.[0-9]+)([eEdDqQ][+-]?[0-9]+)?$'
    local count tok

    [ "$marks" -eq 1 ] || return 1
    count=$(echo $values | wc -w)
    [ "$count" -eq "$width" ] || return 1
    for tok in $values; do
        [[ $tok =~ $real ]] || return 1
    done
    return 0
}

if [ "$1" = "--selftest" ]; then
    fail=0

    accept() {
        if marker_ok "${2:-1}" 1 "$1"; then :; else
            echo " FAIL : the result contract refused '$1'"
            fail=1
        fi
    }
    refuse() {
        if marker_ok "${2:-1}" 1 "$1"; then
            echo " FAIL : the result contract accepted '$1'"
            fail=1
        fi
    }

    accept "6.9999999999999991E+00"   # the value this tower actually prints
    accept "7.0000000000000000E+00"
    accept "-3.5"
    accept "7"
    accept ".5"
    accept "2."
    accept "7.0D+00"
    accept "-2.5e-03"

    refuse "7 7"                      # too wide
    refuse "seven"                    # not a number
    refuse "1.2.3"                    # not a number
    refuse ""                         # nothing at all
    refuse "7" 0                      # no marker
    refuse "7" 2                      # two markers

    if [ "$fail" -ne 0 ]; then
        echo "RESULT CONTRACT: the marker check is wrong"
        exit 1
    fi
    echo "result contract: one real sensitivity accepted, malformed markers refused"
fi
