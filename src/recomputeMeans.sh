while getopts r:s: option
do
case "${option}"
in
r) RISK=${OPTARG};;
s) SEA=${OPTARG};;
esac
done


for dir in /DATA/These/Projects/WaterFlowSimulationModel/output/simu/*/     # list directories in the form "/tmp/dirname/"
do
    dir=${dir%*/}      # remove the trailing "/"
    echo ${dir##*/}    # print everything after the final "/"
    python MeanExtraction.py -sea $SEA -risk $RISK -m ${dir##*/}
done

#python MeanExtraction.py -sea 0.63 -risk 1 -m initPeriod10000Timestep1Sealvl30
# Syntax of array in bash
# ARRAY=()
# ARRAY+=('foo')
# ARRAY+=('bar')

# Ref options : https://www.lifewire.com/pass-arguments-to-bash-script-2200571