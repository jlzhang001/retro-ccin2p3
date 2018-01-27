export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export OMP_NUM_THREADS=1

# Clean the path from MATLAB dynamic libraries
p=
for item in $(echo $LD_LIBRARY_PATH | tr : '\n')
do
    if [[ $item != *"matlab"*  ]]
    then
        p=$item:$p
    fi
done
export LD_LIBRARY_PATH=$p

# Disable certificate check
git config --global http.sslverify false

# Use iRODS
. irods_env.sh
