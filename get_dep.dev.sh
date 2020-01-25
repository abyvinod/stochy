#!/bin/sh
if ! [ -x "$(command -v sudo)" ]; then
  echo 'Warning: sudo is not installed.' >&2
  alias sudo=''
fi

sudo apt-get -y install cmake libboost-all-dev libopenblas-dev liblapack-dev 

echo "--------------------------------"
echo "working on armadillo ..."
mkdir -p armadillo
cd armadillo
git checkout 9.200.x
./configure
make
sudo make install
echo "--------------------------------"
cd ../

sudo apt-get -y install libmatio-dev libginac-dev
