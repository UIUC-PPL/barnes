./charmrun +p2\
        ./barnes \
        -killat=4 \
        -dtime=0.00001 \
        -ppc=5 -p=8 \
        -b=3 \
        -yield=4 \
        -in=20.dat \
        -balancePeriod=4 \
        +balancer Orb3dLB_notopo \
        +noAnytimeMigration \
        +traceroot $PWD/logs \
        +logsize 2000000 \
        +LBDebug 2 \
        +LBPeriod 0.01 +LBCommOff ++debug-no-pause

