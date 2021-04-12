#!/bin/bash
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
echo "Hi from bash"
#declare -a dyn_noise
#declare -a obs_noise

#Dataset="long"

Dataset='long'
attribute='_'
tau_min=0
tau_max=5
mask_type=('y')
extent=("inter") #inter or intra # irrelevant for MovingWindow
minlength=5 #minimal timeseries length in years
method=('parcor')
step=(91)
processing=("_anomalysmooth") # "_MovingWindow_normalized")   # "_anomaly" "_meanseason" "_log" "_combi" "_MovingWindow" "_normalized")  at least a _ must be there _ , _log or _combi; options: _ ; _log ; _combi ; _anomaly ; _meanseason ; _MovingWindow
event_mask_use=('None') # 'Jul-Sep-2003')) # 'Jun-Oct-2003' 'Jun-Oct-2004' 'Jun-Oct-2002') # 'Jul-Sep-2003') #('heatwave-2003') #'Jun-Oct-2003'
quality_value=0.9
attribute2='_'
Network2="interMeanNetwork"
#pc_alpha=('None') #Has to be changed in script, set to None

  for p in "${processing[@]}"
      do
        for m in "${method[@]}"
        do
          for s in "${step[@]}"
          do
            for mt in "${mask_type[@]}"
            do
              for e in "${extent[@]}"
              do
                  for ev in "${event_mask_use[@]}"
                    do
                echo "Start Network construction"
            #   bsub -M 8000000 -J "$Network2" python 01_NetworkConstruction.py $t short $tau_min $tau_max $minlength $mt $e $a $m $s $p
                 bsub -M 8000000 -R "rusage[mem=8000]" -J "$Network2" python 01_NetworkConstruction.py long $attribute $tau_min $tau_max $mt $e $minlength $m $s $p $ev $quality_value $attribute2
        #        python 01_NetworkConstruction.py long $attribute $tau_min $tau_max $mt $e $minlength $m $s $p $ev $quality_value $attribute2
                 echo "Start Network Analysis"
              #   bsub -M 8000000 -w "ended($Network2)" julia 02_NetworkAnalysis.jl $t both $tau_min $tau_max $minlength $mt $a $m $s $p
              #   bsub -M 8000000 julia 02_NetworkAnalysis.jl $t both $tau_min $tau_max $minlength $mt $a $m $s $p
                  # -m "node31"
            #    bsub julia 02_NetworkAnalysis.jl  500 $ModelName $tau_min 10 $h [$d,$d,$d] [$o,$o,$o] $minlength $a $m 0 $p
            #    bsub julia 02_NetworkAnalysis.jl  1000 $ModelName $tau_min 10 $h [$d,$d,$d] [$o,$o,$o] $minlength $a $m 0 $p

                  done
              done
            done
          done
        done
      done
