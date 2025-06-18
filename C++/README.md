# S3CA_CPlusVersion
S3CA C++ Version

## Required Library (Single)
Please install [FFTW](http://www.fftw.org/) and add the FFTWLibrary folder when compile.   
If you do not know how to install it you can follow [this](https://github.com/Jingyi-li/SFFT_Guide).  

## GNUmakefile
Please modify and link the FFTW library \
make local library:`make` \
make local library of autossca:`make ssca`

## Run
compile and run run_autossca.cpp
## Results
![image](https://github.com/Jingyi-li/S3CA/blob/main/Figures/speedup_sum_80_fix_2.png)
