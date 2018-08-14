###############
# install R from source
wget http://cran.rstudio.com/src/base/R-latest.tar.gz
tar xvf R-latest.tar.gz
mv R-3.4.4 /nfs/research1/petsalaki/users/vitalii/R
cd /nfs/research1/petsalaki/users/vitalii/R
./configure --prefix=/nfs/research1/petsalaki/users/vitalii/R  --enable-R-shlib --with-x=yes# to make sure you compile R shared library. Otherwise, RStudio will complain.

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
############# EBI cluster
## installing clustermq
## install zeromq
wget https://github.com/zeromq/libzmq/releases/download/v4.2.3/zeromq-4.2.3.tar.gz
tar xzvf zeromq-4.2.3.tar.gz
mv zeromq-4.2.3 /nfs/research1/petsalaki/users/vitalii/zeromq
cd /nfs/research1/petsalaki/users/vitalii/zeromq
./configure --prefix=/homes/vitalii/zeromq
make
make check
make install
## Edit ~/.bashrc NOT ~/.bash_profile (doesn't start the right version of R on EBI cluster) file:
PATH=/homes/vitalii/zeromq/bin:$PATH
export PKG_CONFIG_PATH="$PKG_CONFIG_PATH:/homes/vitalii/zeromq/lib/pkgconfig"
## One need to put the directory that contains file libzmq.pc to the PKG_CONFIG_PATH
export LD_LIBRARY_PATH="/nfs/research1/petsalaki/users/vitalii/R/lib64:/homes/vitalii/zeromq/lib:$LD_LIBRARY_PATH"
# install zrmq and clustermq
install.packages("rzmq")
install.packages('clustermq')
#C. If everything is going OK so far …and you managed to successfully install clustermq in your R. Then try to give it a test example:
library("clustermq")
fx = function(x) x * 2
Q(fx, x=1:3, n_jobs=1)
##########

############# Sanger cluster
## installing clustermq
## install zeromq
wget https://github.com/zeromq/libzmq/releases/download/v4.2.5/zeromq-4.2.5.tar.gz
tar xzvf zeromq-4.2.5.tar.gz
#mv zeromq-4.2.5 $vk7/software/zeromq
cd zeromq-4.2.5
./configure --prefix=$vk7/software/zeromq
make
make check
make install
## Edit ~/.bashrc NOT ~/.bash_profile (doesn't start the right version of R on EBI cluster) file:
PATH=$vk7/software/zeromq/bin:$PATH
export PKG_CONFIG_PATH="$PKG_CONFIG_PATH:$vk7/software/zeromq/lib/pkgconfig"
## One need to put the directory that contains file libzmq.pc to the PKG_CONFIG_PATH
export LD_LIBRARY_PATH="$vk7/software/R/lib64:$vk7/software/zeromq/lib:$LD_LIBRARY_PATH"
# install zrmq and clustermq
install.packages("rzmq")
install.packages('clustermq')
#C. If everything is going OK so far …and you managed to successfully install clustermq in your R. Then try to give it a test example:
library("clustermq")
fx = function(x) x * 2
Q(fx, x=1:3, n_jobs=1)
##########

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
