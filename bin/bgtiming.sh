#!/bin/bash
if [ ! -d testo ]; then
    mkdir testo
fi
# A=tcp://127.0.0.1:9467
time python bin/chloe.py annotate --workers=4 -o testo $(ls testfa/*.fa)
rm -rf testo