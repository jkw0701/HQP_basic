# HQP_basic
* This project is the whole-body controller for non-holonomic mobile manipulator.
* This is based on Hierarchical Quadratic Programming(HQP).
* This is developed for Ubuntu 16.04.

## Related Works
* A. Escande et al, "Hierarchical quadratic programming: Fast online humanoid-robot motion generation", IJRR, 2014.
* S. Kim et al, "Continuous Task Transition Approach for Robot Controller Based on Hierarchical Quadratic Programming", IEEE Robotics and Automation Letters, 2019.

## Dependencies
* cmake(version >= 2.8)
* eigen3
* qpOASES
* RBDL(Rigid Body Dynamics Library)
* V-REP

### RBDL Setup 

#### Installing
```sh
wget https://bitbucket.org/rbdl/rbdl/get/default.zip
unzip default.zip
cd rbdl-rbdl-0879ee8c548a
mkdir build
cd build
cmake ..
make all
sudo make install
```

### qpOASES setup
Download qpOASES [Link](http://www.qpoases.org/go/release) 
```sh
cd qpOASES-3.2.1
mkdir build
cd build
cmake ..
make all
sudo make install
```

#### qpOASES error handling
if error occures, add following line to qpOASES-3.2.1/CMakeLists.txt, below PROJECT(qpOASES CXX), which is line 34

```
add_compile_options(-fPIC)
```

### V-REP setup
#### How to start with V-REP ###
Download V-REP [Link](http://www.coppeliarobotics.com/downloads.html) 

```sh
cd V-REP_PRO_EDU_.. 
./vrep.sh -gREMOTEAPISERVERSERVICE_-1000_FALSE_TRUE
```

To increase the speed of the simulation, we modify the initial setup when the connection between 
the server and the client is made. Please refer to https://github.com/BenjaminNavarro/vrep_remote_api_shm


## Compile 
```sh
cd HQP_basic
mkdir build
cd build
cmake ..
make 
```
## Run the code
1. Start the V-REP with scene Husky_franka_basic.ttt

2. Run executable
```sh
cd build
./test_hqp
```

