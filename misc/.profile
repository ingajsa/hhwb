#
# Standard Path
#
PATH=.:$PATH:/usr/lpp/X11/bin:/usr/sbin

PATH=$PATH:/usr/lpp/bin
export PATH="${HOME}/hhwb/misc:$PATH"
export PYTHONPATH="${HOME}/hhwb:$PYTHONPATH"

# Ergänzung um ~/.kshrc für jede Kornshell zu aktivieren #

ENV=$HOME/.kshrc
export ENV

source /etc/profile.d/modules.sh
module load anaconda/5.0.0
module load anaconda/5.0.0_py3
module load nco 2>/dev/null
module load cdo/1.7.1 2>/dev/null
module load netcdf-c/4.6.1/intel/serial
module load intel 2>/dev/null