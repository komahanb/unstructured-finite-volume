# Check if compiler exists
# command -v gfortran-7 >/dev/null 2>&1 || { echo "I require gfortran-7 but it's not installed.  Aborting. Try sudo apt-get install gfortran-7 and retry." >&2; exit 1; }

# Execute the sequence of commands to build the library
mkdir -p lib
cd src
make clean
make
make install
cd ../
echo "Now proceed to tests to 'make' and run examples"
