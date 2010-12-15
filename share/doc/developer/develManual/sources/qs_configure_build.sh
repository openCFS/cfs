# Define variable for devel directory and go there
DEVDIR=$HOME/Documents/dev
cd $DEVDIR

# Create two build directories
mkdir -p build/release build/debug

# Set environment variables to desired compilers
export CC=gcc; export CXX=g++; export FC=gfortran

# Go to release directory
cd $DEVDIR/build/release

# Configure release build tree
cmake -DDEBUG:BOOL=OFF \
      -DCFS_DEPS_ROOT:PATH=$DEVDIR/CFSDEPS \
      -DCFS_DEPS_CACHE_DIR:PATH=$DEVDIR/CFSDEPSCACHE \
      -G "Unix Makefiles"
      $DEVDIR/CFS_TRUNK

# Build executables (available in bin/...) using
# two processors
make -j2

# Go to debug build tree, configure and build.
cd $DEVDIR/build/debug

cmake -DDEBUG:BOOL=ON \
      -DCFS_DEPS_ROOT:PATH=$DEVDIR/CFSDEPS \
      -DCFS_DEPS_CACHE_DIR:PATH=$DEVDIR/CFSDEPSCACHE \
      -G "Unix Makefiles"
      $DEVDIR/CFS_TRUNK

make -j2
