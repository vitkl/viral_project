# install R from source
wget http://cran.rstudio.com/src/base/R-latest.tar.gz
tar xvf R-latest.tar.gz
mv R-3.4.1 /hps/nobackup/research/petsalaki/users/vitalii/R
cd /hps/nobackup/research/petsalaki/users/vitalii/R
./configure --prefix=/hps/nobackup/research/petsalaki/users/vitalii/R  --enable-R-shlib --with-x=yes# to make sure you compile R shared library. Otherwise, RStudio will complain.

# next, edit or create .bash_profile

# Get the aliases and functions
if [ -f ~/.bashrc ]; then
. ~/.bashrc
fi

# User specific environment and startup programs
# set PATH so it includes user's private bin if it exists

if [ -d "$HOME/bin" ];then
PATH="$HOME/bin:$PATH"
fi
PATH=/hps/nobackup/research/petsalaki/users/vitalii/R/bin:$PATH
export PATH
## end of edit or create .bash_profile

# compile R
make && make install


./configure --prefix=/nfs/research2/thornton/veronika/R/ --enable-R-shlib

# https://blogs.msdn.microsoft.com/gpalem/2013/02/06/building-zeromq-package-for-r-rzmq-from-source/
# configure zeromq
./configure --prefix=/ebi/research/thornton/veronika/zeromq-4.2.2/ --exec-prefix=/nfs/research2/thornton/veronika/zeromq-4.2.2/

# install rzmq
export PKG_CONFIG_PATH=/ebi/research/thornton/veronika/zeromq-4.2.2/src/
LD_LIBRARY_PATH=/ebi/research/thornton/veronika/zeromq-4.2.2/lib/:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH

install.packages("rzmq", configure.vars = "INCLUDE_DIR=/ebi/research/thornton/veronika/zeromq-4.2.2/include LIB_DIR=/ebi/research/thornton/veronika/zeromq-4.2.2/lib/")

# load rzmq
LD_LIBRARY_PATH=/ebi/research/thornton/veronika/zeromq-4.2.2/lib/:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH
library(rzmq)
