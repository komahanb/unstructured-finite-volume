#=====================================================================#
# Contains a derived type 'string' and implemented procedures
#
# Author: Komahan Boopathy (komahan@gatech.edu)
#=====================================================================#


class String:
    """Derived type for string."""

    #===================================================================#
    # Construct a string object from the supplied literal, find its
    # length, initialize its hashcode as zero.
    #
    # Fortran: interface string => create(str)
    #===================================================================#
    def __init__(self, str):
        self.str = str            # character array
        self.count = len(str)     # length

    #===================================================================#
    # Tokenize the string object and return an array of tokens
    #
    # Fortran: subroutine tokenize(this, delimiter, num_tokens, tokens)
    # The Fortran out-arguments num_tokens and tokens are returned here.
    # `tokens` is built (as a list of String) only when want_tokens is set,
    # mirroring the optional `tokens` argument.
    #===================================================================#
    def tokenize(self, delimiter, want_tokens=True):
        num_tokens = 0
        if len(delimiter) == 0:
            return num_tokens, None              # doesnt match
        if self.str.find(delimiter) == -1:       # index(...) .eq. 0, doesnt match
            return num_tokens, None

        # Lower and upper index of tokens (1-based inclusive, as in Fortran)
        tidx = []

        # Initialize
        sidx = 1
        eidx = len(self.str)
        token_ctr = 0
        while True:                              # do while (len(str(sidx:eidx)) >= 0)

            # Get the -th index of delimiter (1-based; 0 if not found)
            token_idx = self.str[sidx - 1:eidx].find(delimiter) + 1

            if token_idx != 0:

                token_ctr = token_ctr + 1
                tidx.append((sidx, token_idx - 1 + sidx))

                # We found the match record the index
                sidx = sidx + token_idx

            else:

                # Check if its the last substring
                if token_ctr >= 1:

                    # Yes, this is a token
                    token_ctr = token_ctr + 1

                    token_idx = 1

                    tidx.append((sidx, eidx))

                break

        # Set the return arguments
        num_tokens = token_ctr
        tokens = None
        if want_tokens:
            tokens = [String(self.str[tidx[i][0] - 1:tidx[i][1]])
                      for i in range(num_tokens)]

        return num_tokens, tokens

    #===================================================================#
    # Overridden string equality logic. Based on comparison of entries
    #===================================================================#
    def equals(self, element):
        # string objects are equal if their values are equal
        return element.str == self.str

    #===================================================================#
    # Returns the string representation of the object
    #===================================================================#
    def print(self):
        if self.str is not None:
            print(" string : ", self.str)
        else:
            print(" string : ", "NULL")

    #===================================================================#
    # Return the integer evaluation of string
    #===================================================================#
    def asinteger(self):
        return int(self.str)

    #===================================================================#
    # Get the real number from string
    #===================================================================#
    def asreal(self):
        return float(self.str)
