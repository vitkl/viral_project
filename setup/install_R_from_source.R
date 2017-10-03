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
