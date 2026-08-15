#!/bin/bash
# The time integration tower's RESULT CONTRACT, in one place: what a
# TIME_INTEGRATION_RESULT marker must look like, and a self-test
# proving the check accepts genuine real values and refuses
# malformed ones.
#
# The contract is shape and syntax only:
#
#     exactly one marker line
#     exactly two tokens - one per state coordinate, x then y, in
#     Q's declaration order
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
# tower's answer is q(t4) on the two-member state set Q; the
# marker carries what the program COMPUTED, unrounded. Whether
# those numbers are 7/24 and 83/144 is the Level-9 test's business,
# never this script's.
#
# TWO tokens, not five. The result lives on Q, and the tower spent
# nine levels establishing that Q is not the time axis, not the
# steps and not any graph's vertex set. A five-token marker would
# quietly undo that at the last line of output.

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

    accept "2.9166666666666674E-01 5.7638888888888895E-01"
    accept "0.2916666666666667 0.5763888888888889"
    accept "-1.5 .5"
    accept "2. 7.0D+00"

    refuse "2.9166666666666674E-01"          # too narrow
    refuse "0.29 0.58 0.11"                  # too wide - Q has two members
    refuse "0.29 five"                       # not a number
    refuse "0.29 1.2.3"                      # not a number
    refuse "0.29 NaN"                        # letters the grammar refuses
    refuse ""                                # nothing at all
    refuse "0.29 0.58" 0                     # no marker
    refuse "0.29 0.58" 2                     # two markers

    if [ "$fail" -ne 0 ]; then
        echo "RESULT CONTRACT: the marker check is wrong"
        exit 1
    fi
    echo "result contract: one two-member real state accepted, malformed markers refused"
fi
