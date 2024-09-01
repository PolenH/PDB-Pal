#!/bin/bash

# Loop through directories 0 to 19
for i in {0..19}
do
    cd Out/$i
    
    python final2.py
    
    cd ../../
done