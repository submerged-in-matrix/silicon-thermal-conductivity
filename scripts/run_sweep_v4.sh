#!/bin/bash
# Parametric sweep for v4 thermal conductivity
# dt=0.5fs, swap/100 steps, 11 slabs
echo "Starting v4 thermal conductivity sweep..."
echo "dt=0.0005 ps, swap every 100 steps, 20 slabs"
echo ""

for LEN in 20 30 40 60; do
    echo "Running LEN = ${LEN}..."
    cd ~/silicon-thermal-conductivity/inputs && lmp -in in.si_thermal_v4 -var LEN $LEN > ../outputs/log_v4_L${LEN}.txt 2>&1

    if grep -q "RESULT_v4" ../outputs/log_v4_L${LEN}.txt; then
        echo "  -> Done: $(grep RESULT_v4 ../outputs/log_v4_L${LEN}.txt)"
    else
        echo "  -> FAILED! Check ../outputs/log_v4_L${LEN}.txt"
    fi
done
echo "All v4 runs complete."
