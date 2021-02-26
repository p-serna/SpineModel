#! /bin/bash
echo "We configure stuff"

echo "We build iv"
cd /neuron/iv
./build.sh
./configure --prefix=/neuron/iv/
make
make install
#! /bin/bash
echo "We compile methods nrniv"

cd /SpineModel/mod/
nrnivmodl

mkdir /SpineModel/fconditions
mkdir /SpineModel/fconditions/spatial
echo "FINISHED!"

# docker exec -it b91c015c9706 /bin/bash
