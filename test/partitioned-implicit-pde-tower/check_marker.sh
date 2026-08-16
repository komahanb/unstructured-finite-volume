#!/bin/bash
# The partitioned PDE tower's RESULT CONTRACT, in one place: what a
# PARTITIONED_PDE_RESULT marker must look like, and a self-test
# proving the check accepts genuine real values and refuses
# malformed ones.
#
# The contract is shape and syntax only:
#
#     exactly one marker line
#     exactly six tokens - one per global vertex, in global order
#     every token a real numeric literal
#
# The token check is LEXICAL, not evaluative: the grammar admits a
# signed decimal with an optional exponent, which excludes NaN and
# Inf spellings because they carry letters it does not accept - but
# it does not evaluate, so an overflowing literal such as 1e999999
# would pass the syntax and parse to infinity. Calling this a
# FINITE check would be an overclaim.
#
# And it says NOTHING about what the numbers should be. This
# tower's answer is a real-valued field; the marker carries the
# field as it stands - 1.0000000000000002E+00 and its neighbours -
# which is the honest image of the arithmetic. Whether that IS q*
# is the Gate-C test's business, never this script's.

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
        if marker_ok "${2:-1}" 6 "$1"; then :; else
            echo " FAIL : the result contract refused '$1'"
            fail=1
        fi
    }
    refuse() {
        if marker_ok "${2:-1}" 6 "$1"; then
            echo " FAIL : the result contract accepted '$1'"
            fail=1
        fi
    }

    accept "1.0000000000000002E+00 2.0 4.0 7.0 1.1E+01 1.6E+01"
    accept "1 2 4 7 11 16"
    accept "-1.5 .5 2. 3e-4 -2.5E+03 7.0D+00"

    refuse "1 2 4 7 11"          # too narrow
    refuse "1 2 4 7 11 16 21"    # too wide
    refuse "1 2 4 7 11 six"      # not a number
    refuse "1 2 4 7 11 1.2.3"    # not a number
    refuse ""                    # nothing at all
    refuse "1 2 4 7 11 16" 0     # no marker
    refuse "1 2 4 7 11 16" 2     # two markers

    if [ "$fail" -ne 0 ]; then
        echo "RESULT CONTRACT: the marker check is wrong"
        exit 1
    fi
    echo "result contract: one six-member real field accepted, malformed markers refused"
fi
