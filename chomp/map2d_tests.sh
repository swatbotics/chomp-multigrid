SETTINGS_V1="-n 127 -g 0.5 -a 0.02 -e 1e-6 -m 200 -o vel"
SETTINGS_V0="-n  63 -g 0.8 -a 0.03 -e 1e-6 -m 200 -o vel"

SETTINGS_A1="-n 127 -g 0.03 -a 0.001 -e 1e-12 -m 400 -o accel"
SETTINGS_A0="-n  63 -g 0.08 -a 0.002 -e 1e-12 -m 400 -o accel"


MAP3A="-c 2.7,-2.7,-2.7,2.7 ../maps/map3.txt"
MAP3B="-c -2.7,-2.7,2.7,2.7 ../maps/map3.txt"

../build/map2d_demo -p10 $SETTINGS_V1 $MAP3A
../build/map2d_demo -p10 $SETTINGS_V1 $MAP3B
../build/map2d_demo -p10 $SETTINGS_V0 $MAP3A
../build/map2d_demo -p10 $SETTINGS_V0 $MAP3B

../build/map2d_demo -p20 $SETTINGS_A1 $MAP3A
../build/map2d_demo -p20 $SETTINGS_A1 $MAP3B
../build/map2d_demo -p20 $SETTINGS_A0 $MAP3A
../build/map2d_demo -p20 $SETTINGS_A0 $MAP3B

