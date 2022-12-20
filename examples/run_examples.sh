#! /bin/bash
# run this file to make sure all examples work
FILES="Grover QAOA HHL QCBM QuGAN Shor VQE GateLearning GroundStateSolvers BernsteinVazirani MerminMagic"
for folder in $FILES
do
    echo "executing example: $folder"
    julia --project=$folder -e "using Pkg; Pkg.instantiate()"
    #julia --project=$folder -e "using Pkg; Pkg.update()"
    julia --project=$folder "$folder/main.jl"
done
