export CurrDir=$PWD

mkdir -p /scratch/$USER/${SLURM_JOBID}
export ScrDir=/scratch/$USER/${SLURM_JOBID}

cp $CurrDir/$1 $ScrDir/
cd $ScrDir

# optional PAL modification
check_par=` grep -c -i pal $1 `
if [ $check_par -eq 0 ]; then
          if [ ! -z $SLURM_NTASKS ]; then
                   if [ $SLURM_NTASKS -ne 1 ]; then
                      nprocs=$SLURM_NTASKS
                      check_par=1
                   fi
           fi
           if [ $check_par -eq 1 ]; then
           echo "%pal nprocs $nprocs end" >> $1
    fi
fi

# copy zyx files  
for line in `cat $1`
do
        if [[ $line == *'xyz"'* ]]
        then
                echo copying $line to workdir
                filename=$line|tr '"' ''
                cp ./$filename $ScrDir
        fi
done

cp $CurrDir/*.xyz $ScrDir/

# copy tmp gbws to the ScrDir
inpfile=$1
input_gbw_file=` grep -c -i 'tmp.gbw' $inpfile `

if [ $input_gbw_file -eq 1 ]
then
 gwb_file=` grep -i moinp $inpfile | cut -f 2 -d '"' `
 cp $CurrDir/tmp.gbw $ScrDir
fi

cp $CurrDir/tmp.gbw $ScrDir
ls $ScrDir


# execute ORCA
path/to/orca $ScrDir/$1 > $CurrDir/$1.out

# Copy and delete ScrDir
cd $CurrDir
cp $ScrDir/* .

cd $CurrDir
rm -rf $ScrDir
