#!/bin./bash

PROGRAM="mono./bin/cRho.exe"
DIR1="./JSON"
DIR2="./local"

#run
$PROGRAM &
PID=$!

echo ":::: cRho running : PID: $PID"

# wait unitl DIR1, DIR2 are created
while [ ! -d "$DIR1" ] || [ ! -d "$DIR2" ]; do
	echo "::::: waiting folder creation"
	sleep 1
done

chmod a+x ./elris/triangle_*
chmod a+x update.sh

kill $PID

echo "::::: Basic setting done. Reboot OR run by mono ./bin/cRho.exe"