./charmrun +p8 ./barnes -killat=10 -dtime=0.00001 -ppc=2000 -p=140 -in=100k.dat -balancePeriod=10 +balancer Orb3dLB_notopo +noAnytimeMigration +LBPeriod 0.01 +LBCommOff out.nolb
./charmrun +p8 ./barnes -killat=10 -dtime=0.00001 -ppc=2000 -p=140 -in=100k.dat -balancePeriod=3 +balancer Orb3dLB_notopo +noAnytimeMigration +LBPeriod 0.01 +LBCommOff out.lb
