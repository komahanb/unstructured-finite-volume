#!/bin/bash
# The derivative tower's RESULT CONTRACT, in one place: what a
# DERIVATIVE_RESULT marker must look like, and a self-test proving
# the check accepts genuine real values and refuses malformed ones.
#
# The contract is shape and syntax only:
#
#     exactly one marker line
#     exactly one token per member of the independent domain
#     every token a real numeric literal
#
# The token check is LEXICAL, not evaluative: the grammar admits a
# signed decimal with an optional exponent, which excludes NaN and
# Inf spellings because they carry letters it does not accept - but
# it does not evaluate, so an overflowing literal such as 1e999999
# would pass the syntax and parse to infinity. Calling this a
# FINITE check would be an overclaim, and the value is not this
# script's business in any case:
#
# and NOTHING about what the numbers should be. This tower's answer
# is a real-valued field; the marker carries the field as it stands,
# so the check must admit decimals and exponents - never integers
# alone, which would silently demand a rounded answer.

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
        if marker_ok "${2:-1}" 2 "$1"; then :; else
            echo " FAIL : the result contract refused '$1'"
            fail=1
        fi
    }
    refuse() {
        if marker_ok "${2:-1}" 2 "$1"; then
            echo " FAIL : the result contract accepted '$1'"
            fail=1
        fi
    }

    # A real-valued field must survive serialization intact.
    accept "3.0000000000000000E+00  3.0000000000000000E+00"
    accept "2.5 -1.25"
    accept ".5 2."
    accept "3e-4 -2.5E+03"
    accept "3.0D+00 -1.0D-02"
    accept "3 3"                 # integral values, honestly written
    accept "-10 -4"

    refuse "3"                   # too narrow
    refuse "3 3 3"               # too wide
    refuse "3 oops"              # not a number
    refuse "1.2.3 4"             # not a number
    refuse "e5 3"                # not a number
    refuse "- 3"                 # not a number
    refuse ""                    # nothing at all
    refuse "3 3" 0               # no marker
    refuse "3 3" 2               # two markers

    if [ "$fail" -ne 0 ]; then
        echo "RESULT CONTRACT: the marker check is wrong"
        exit 1
    fi
    echo "result contract: real vectors accepted, malformed markers refused"
fi
