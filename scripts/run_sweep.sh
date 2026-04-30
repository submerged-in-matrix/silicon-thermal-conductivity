#!/bin/bash
echo "Starting thermal conductivity sweep..."
for LEN in 20 30 40 60; do
    echo "Running LEN = ${LEN}..."
    lmp -in in.si_thermal_v3 -var LEN $LEN > ../outputs/log_L${LEN}.txt 2>&1
    
    # Check if it completed
    if grep -q "RESULT" ../outputs/log_L${LEN}.txt; then
        echo "  -> Done: $(grep RESULT ../outputs/log_L${LEN}.txt)"
    else
        echo "  -> FAILED! Check ../outputs/log_L${LEN}.txt"
    fi
done
echo "Sweep complete."
